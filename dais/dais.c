// Copyright 2016 Robert W. Fuller <hydrologiccycle@gmail.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// dais.c
// written by Robert W. Fuller on 160621
//

#include "r.h"


static RVector modelParms;

static RMatrix forcings;
#define getForcing(c, r) getMatrixElem(&forcings, (c), (r))
#define Ta(r)   getForcing(1, (r))
#define SL(r)   getForcing(2, (r))
#define GSL(r)  getForcing(3, (r))
#define Toc(r)  getForcing(4, (r))

static RVector  output[3];
#define getOut(v, r)    output[ (v) ].comn.dbl_arr[ (r) - 1 ]
#define Rad(r)          getOut(0, (r))  // Radius of ice sheet
#define Vais(r)         getOut(1, (r))  // Ice volume
#define SLE(r)          getOut(2, (r))  // Sea-level equivalent [m]


// extra parameters:  constants or parameters that are temporarily fixed
static RVector extraParms;

// logicals to enable/disable parts of models
static RIntVector switches;

static RParm parms[] = {
    { "mp",     &modelParms,    RTYPE_VECTOR,       0 },
    { "frc",    &forcings,      RTYPE_MATRIX,       0 },
    { "out",    &output,        RTYPE_VECTORS,      NELEMS(output) },
    { "ep",     &extraParms,    RTYPE_VECTOR,       0 }, 
    { "sw",     &switches,      RTYPE_INT_VECTOR,   0 }
};


static double b0, slope, mu, h0, c, P0, kappa, nu, f0, Gamma, alpha, Tf, rho_w, rho_i, rho_m, Toc_0, Rad0;

static RNamedReal realParms[] = {
    { "b0",     &b0 },
    { "slope",  &slope },
    { "mu",     &mu },
    { "h0",     &h0 },
    { "c",      &c },
    { "P0",     &P0 },
    { "kappa",  &kappa },
    { "nu",     &nu },
    { "f0",     &f0 },
    { "Gamma",  &Gamma },
    { "alpha",  &alpha },
    { "Tf",     &Tf },
    { "rho_w",  &rho_w },
    { "rho_i",  &rho_i },
    { "rho_m",  &rho_m },
    { "Toc_0",  &Toc_0 },
    { "Rad0",   &Rad0 }
};


static int sw_simple_vol;

static RNamedInt swParms[] = {
    { "simple_vol",        &sw_simple_vol }
};


static DL_FUNC get_deSolve_gparms;

void R_init_dais(DllInfo *dll)
{
    get_deSolve_gparms = R_GetCCallable("deSolve", "get_deSolve_gparms");

    sortNamedStructs(realParms);
    sortNamedStructs(parms);
    sortNamedStructs(swParms);
}


SEXP daisInit(SEXP gparms)
{
    modelParms.comn.s_ptr = NULL;
    extraParms.comn.s_ptr = NULL;
    switches.comn.s_ptr   = NULL;

    initParms(parms, NELEMS(parms), gparms);

    if (modelParms.comn.s_ptr != NULL) {
        initNamedReals(realParms, NELEMS(realParms), &modelParms);
    }
    if (extraParms.comn.s_ptr != NULL) {
        initNamedReals(realParms, NELEMS(realParms), &extraParms);
    }
    if (switches.comn.s_ptr != NULL) {
        initNamedInts(swParms, NELEMS(swParms), &switches);
        //Rprintf("simple_vol is %d\n", sw_simple_vol);
    }

    return R_NilValue;
}


void daisOdeInit(void (*odeparms)(int *, double *))
{
    SEXP gparms;

    gparms = (SEXP) get_deSolve_gparms();

    daisInit(gparms);
}


#define Volo 2.4789e16

SEXP daisOdeC()
{
    double del, eps1, eps2, R, rc, hr, P, beta, rR, Btot, mit, F, ISO, Hw, Speed, fac;
    int i, np;

    // Initialize intermediate parameters
    del  = rho_w/rho_i;           // Ratio sea water and ice density [-]
    eps1 = rho_i/(rho_m - rho_i);  // Ratio ice density and density difference between rock and ice [-]
    eps2 = rho_w/(rho_m - rho_i);  // Ratio sea water density and density difference between rock and ice [-]

    np = forcings.rows;

    // Initial conditions
    R  = Rad0;                    // gets updated at end of loop
    rc = (b0 - SL(1))/slope;      // application of equation 1 (paragraph after eq3)
    mit = (R <= rc) ? 0.0 : 1.0;  // marine ice sheet or not?
    Rad(1)  = R;
    Vais(1) = M_PI * (1+eps1) * ( (8/15) * sqrt(mu) * pow(R, 2.5) - 1/3*slope*pow(R, 3))
            - mit * M_PI*eps2 * ( (2/3) * slope*(pow(R, 3)-pow(rc, 3))-b0*(R*R-rc*rc) );
    SLE(1)  = 57*(1-Vais(1)/Volo);  // Takes steady state present day volume to correspond to 57m SLE

    // Run model
    for (i = 2;  i <= np;  ++i) {
    
        hr = h0 + c * Ta(i-1);        // equation 5
        rc = (b0 - SL(i-1))/slope;    // application of equation 1 (paragraph after eq3)
        P  = P0 * exp(kappa*Ta(i-1)); // equation 6
        beta = nu * sqrt(P);          // equation 7 (corrected with respect to text)

        // Total mass accumulation on ice sheet (equation 8)
        if (hr > 0) {
            rR = R - (pow(hr - b0 + slope*R, 2) / mu);

            Btot = P * M_PI * R*R
                 - M_PI * beta * (hr - b0 + slope*R) * (R*R - rR*rR)
                 - (4 * M_PI * beta * sqrt(mu) *   pow(R-rR, 2.5)) / 5
                 + (4 * M_PI * beta * sqrt(mu) * R*pow(R-rR, 1.5)) / 3;
        } else {
            Btot = P * M_PI*R*R;
        }

        // Ice flux at grounding line (F), ISO and fac
        if (R <= rc) {

            // no grounding line / ice shelves
            mit = 0;  // marine ice term (if mit=0 marine ice sheet terms are ignored)
            F   = 0;  // no ice flux
            ISO = 0;  // (third term equation 14) NAME?
        } else {
            // marine ice sheet / grounding line
            mit = 1;

            Hw = slope*R - b0 + SL(i-1);  // (equation 10) 

            // Ice speed at grounding line (equation 11)
            // KELSEY: last term is different than in manuscript! (slope*Rad0 - b0) rather than (b0 - slope*Rad0)
            // KELSEY: this seems to corrected for in equation 9 that misses a minus sign (wrt ms)
            Speed = f0
                  * ((1-alpha) + alpha * pow((Toc(i-1) - Tf)/(Toc_0 - Tf), 2))
                  * pow(Hw, Gamma) / pow(slope*Rad0 - b0, Gamma - 1);

            F     = 2*M_PI*R * del * Hw * Speed;  // (equation 9)

            // Kelsey uses GSL instead of calculating rate herself
            //ISO   = 2*M_PI*eps2* (slope*rc*rc - (b0/slope)*rc) * (SL(i) - SL(i-1));  // third term equation 14 !! NAME?
            ISO   = 2*M_PI*eps2* (slope*rc*rc - (b0/slope)*rc) * GSL(i);  // third term equation 14 !! NAME?
        }

        // dV/dR (equation 14)
        fac = M_PI * (1+eps1) * (4/3 * sqrt(mu) * pow(R, 1.5) - slope*R*R)
            - mit * ( 2*M_PI*eps2 * (slope*R*R - b0*R) );  // in case of marine ice sheet (mit=1)

        // Ice sheet volume (equation 13)
        // KELSEY: it seems some parantheses were missing in your code
        R       = R + (Btot-F+ISO)/fac;
        Rad(i)  = R;

        if (sw_simple_vol) {
            Vais(i) = Vais(i-1) + (Btot-F+ISO);
        } else {
            Vais(i) = M_PI * (1+eps1) * ( (8/15) *  sqrt(mu) * pow(R, 2.5) - 1/3*slope*pow(R, 3))
                    - mit * M_PI*eps2 * ( (2/3) * slope*(pow(R, 3)-pow(rc, 3))-b0*(R*R-rc*rc) );
        }
        SLE(i)  = 57*(1-Vais(i)/Volo);  // Takes steady state present day volume to correspond to 57m SLE
    }

    return R_NilValue;
}
