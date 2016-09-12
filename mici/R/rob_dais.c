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
#define Toc(r)  getForcing(2, (r))
#define GSL(r)  getForcing(3, (r))
#define SL(r)   getForcing(4, (r))

static RVector output[5];
#define getOut(v, r)    output[ (v) ].comn.dbl_arr[ (r) - 1 ]
#define SLE(r)          getOut(0, (r))  // Sea-level equivalent [m]
#define Vais(r)         getOut(1, (r))  // Ice volume
#define Rad(r)          getOut(2, (r))  // Radius of ice sheet
#define Flow(r)         getOut(3, (r))  // Ice flow
#define Depth(r)        getOut(4, (r))  // Water depth


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


static double b0, slope, mu, h0, c, P0, kappa, nu, f0, Gamma, alpha, Tf, rho_w, rho_i, rho_m, Toc_0, Rad0, Tcrit, lambda;

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
    { "gamma",  &Gamma },
    { "alpha",  &alpha },
    { "Tf",     &Tf },
    { "rho_w",  &rho_w },
    { "rho_i",  &rho_i },
    { "rho_m",  &rho_m },
    { "Toc_0",  &Toc_0 },
    { "Rad0",   &Rad0 },
    { "Tcrit",  &Tcrit },
    { "lambda", &lambda }
};


static int sw_fast_dyn;

static RNamedInt swParms[] = {
    { "fast_dyn",       &sw_fast_dyn }
};


void R_init_rob_dais(DllInfo *dll)
{
    sortNamedStructs(realParms);
    sortNamedStructs(parms);
    sortNamedStructs(swParms);
}


static SEXP daisInit(SEXP gparms)
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
        //Rprintf("fast_dyn is %d\n", sw_fast_dyn);
    }

    return R_NilValue;
}


#define Volo 2.4789e16

static SEXP daisOdeSolve()
{
    double del, eps1, eps2, R, rc, hr, P, beta, rR, Btot, F, ISO, Hw, Speed, fac, disint_rate;
    int i, np;

    // Initialize intermediate parameters
    del  = rho_w/rho_i;            // Ratio sea water and ice density [-]
    eps1 = rho_i/(rho_m - rho_i);  // Ratio ice density and density difference between rock and ice [-]
    eps2 = rho_w/(rho_m - rho_i);  // Ratio sea water density and density difference between rock and ice [-]

    np = forcings.rows;

    // diff -u dais.R.ORIG dais.R 

    // Initial conditions
    R  = Rad0;                    // gets updated at end of loop
    rc = (b0 - SL(1))/slope;      // application of equation 1 (paragraph after eq3)
    Rad(1)  = R;
    Vais(1) = M_PI * (1.0+eps1) * ( (8.0/15.0) * sqrt(mu) * pow(R, 2.5) - 1.0/3.0*slope*pow(R, 3.0));
    if (R > rc) {
        Vais(1) -= M_PI*eps2 * ( (2.0/3.0) * slope*(pow(R, 3.0)-pow(rc, 3.0))-b0*(R*R-rc*rc) );
    }
    SLE(1)  = 57.0*(1.0-Vais(1)/Volo);  // Takes steady state present day volume to correspond to 57m SLE
    disint_rate = 0;

    // Run model
    for (i = 2;  i <= np;  ++i) {
    
        hr = h0 + c * Ta(i-1);        // equation 5
        rc = (b0 - SL(i-1))/slope;    // application of equation 1 (paragraph after eq3)
        P  = P0 * exp(kappa*Ta(i-1)); // equation 6
        beta = nu * sqrt(P);          // equation 7 (corrected with respect to text)

        // Total mass accumulation on ice sheet (equation 8)
        if (hr > 0) {
            rR = R - (pow(hr - b0 + slope*R, 2.0) / mu);

            Btot = P * M_PI * R*R
                 - M_PI * beta * (hr - b0 + slope*R) * (R*R - rR*rR)
                 - (4.0 * M_PI * beta * sqrt(mu) *   pow(R-rR, 2.5)) / 5.0
                 + (4.0 * M_PI * beta * sqrt(mu) * R*pow(R-rR, 1.5)) / 3.0;
        } else {
            Btot = P * M_PI*R*R;
        }

        fac = M_PI * (1.0+eps1) * (4.0/3.0 * sqrt(mu) * pow(R, 1.5) - slope*R*R);  // ratio dV/dR

        // this would go in the first case below, but want to save the water depth
        Hw = slope*R - b0 + SL(i-1);  // (equation 10)
        Depth(i-1) = Hw;

        // In case there is a marine ice sheet / grounding line
        if (R > rc) {

            fac -= ( 2.0*M_PI*eps2 * (slope*R*R - b0*R) );  // correction fac (=ratio dV/dR)

            //Hw = slope*R - b0 + SL(i-1);  // (equation 10)

            // Ice speed at grounding line (equation 11)
            // KELSEY: last term is different than in manuscript! (slope*Rad0 - b0) rather than (b0 - slope*Rad0)
            // KELSEY: this seems to corrected for in equation 9 that misses a minus sign (wrt ms)
            Speed = f0
                  * ((1.0-alpha) + alpha * pow((Toc(i-1) - Tf)/(Toc_0 - Tf), 2.0))
                  * pow(Hw, Gamma) / pow(slope*Rad0 - b0, Gamma - 1.0);

            F     = 2.0*M_PI*R * del * Hw * Speed;  // Ice flux (equation 9)

            // Alex now uses GSL instead of calculating rate;  Kelsey always used GSL
          //ISO   = 2.0*M_PI*eps2* (slope*rc*rc - (b0/slope)*rc) * (SL(i) - SL(i-1));  // third term equation 14 !! NAME?
            ISO   = 2.0*M_PI*eps2* (slope*rc*rc - (b0/slope)*rc) * GSL(i-1);  // third term equation 14 !! NAME?

        } else {

            // In case there is no marine ice sheet / grounding line
            F   = 0;    // no ice flux
            ISO = 0;    // (third term equation 14) NAME?
        }

        Flow(i-1) = F;

        if (sw_fast_dyn) {
            if (Ta(i-1) > Tcrit) {
                // Takes steady state present day volume to correspond to 57m SLE
                disint_rate = -lambda * Volo / 57.0;
            } else {
                disint_rate = 0;
            }
        }

        // Ice sheet volume (equation 13)
        R      += (Btot-F+ISO+disint_rate)/fac;
        Rad(i)  = R;

        Vais(i) = Vais(i-1) + (Btot-F+ISO+disint_rate);
        SLE(i)  = 57.0*(1.0-Vais(i)/Volo);  // Takes steady state present day volume to correspond to 57m SLE
    }

    // calculate final depth and ice flux
    Hw = slope*R - b0 + SL(np);  // (equation 10)
    Depth(np) = Hw;
    rc = (b0 - SL(np))/slope;
    if (R > rc) {
        Speed = f0
              * ((1.0-alpha) + alpha * pow((Toc(np) - Tf)/(Toc_0 - Tf), 2.0))
              * pow(Hw, Gamma) / pow(slope*Rad0 - b0, Gamma - 1.0);
        F     = 2.0*M_PI*R * del * Hw * Speed;  // Ice flux (equation 9)
    } else {
        F     = 0;
    }
    Flow(np)  = F;

    return R_NilValue;
}


SEXP daisRobOdeC(SEXP gparms)
{
    SEXP rc;

    daisInit(gparms);
    rc = daisOdeSolve();

    return rc;
}
