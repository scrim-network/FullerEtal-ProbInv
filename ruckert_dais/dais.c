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

static RVector  output[3];
#define getOut(v, r)    output[ (v) ].comn.dbl_arr[ (r) - 1 ]
#define SLE(r)          getOut(0, (r))  // Sea-level equivalent [m]
#define Vais(r)         getOut(1, (r))  // Ice volume
#define Rad(r)          getOut(2, (r))  // Radius of ice sheet


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
    { "gamma",  &Gamma },
    { "alpha",  &alpha },
    { "Tf",     &Tf },
    { "rho_w",  &rho_w },
    { "rho_i",  &rho_i },
    { "rho_m",  &rho_m },
    { "Toc_0",  &Toc_0 },
    { "Rad0",   &Rad0 }
};


static int sw_complex_vol;

static RNamedInt swParms[] = {
    { "complex_vol",        &sw_complex_vol }
};


static DL_FUNC get_deSolve_gparms;

void R_init_dais(DllInfo *dll)
{
    //get_deSolve_gparms = R_GetCCallable("deSolve", "get_deSolve_gparms");

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
        //Rprintf("complex_vol is %d\n", sw_complex_vol);
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

SEXP daisOdeCdeSolve()
{
    double del, eps1, eps2, R, rc, hr, P, beta, rR, Btot, F, ISO, Hw, f, fac, dR, V, Vsea, Volt;
    int i, np;

    // Initialize intermediate parameters
    del  = rho_w/rho_i;            // Ratio sea water and ice density [-]
    eps1 = rho_i/(rho_m - rho_i);  // Ratio ice density and density difference between rock and ice [-]
    eps2 = rho_w/(rho_m - rho_i);  // Ratio sea water density and density difference between rock and ice [-]
    R  = Rad0;                    // gets updated at end of loop

    np = forcings.rows;

    // diff -u dais.R.ORIG dais.R
    
    // Run model
    for (i = 1;  i <= np;  ++i) {
        // function used in calculating Ice speed at grounding line and Ice Flux, F, (part of equation 11)
        f = f0 * ((1.0-alpha) + alpha * pow((Toc(i) - Tf)/(Toc_0 - Tf), 2.0)) / pow((slope*Rad0 - b0), (Gamma - 1.0));
        hr = h0 + c * Ta(i);        // equation 5
        rc = (b0-SL(i)) / slope;     // application of equation 1 (paragraph after eq3)
        P = P0 * exp(kappa*Ta(i));   // equation 6
        beta = nu * pow(P, 0.5);         // equation 7 (corrected with respect to text)
        rR = R - (fabs(hr - b0 + slope * R) * (hr - b0 + slope * R))/mu;
        
        if (R <= rR && R <= rc) {
            // Total mass accumulation on ice sheet (equation 8)
            Btot = M_PI * P * R * R;
            // In case there is no marine ice sheet / grounding line
            F = 0;       // no ice flux
            ISO = 0;     // (third term equation 14) NAME?
            fac = M_PI * (1.0 + eps1) * (4.0/3.0 * pow(mu, 0.5) * pow(R, 1.5) - slope*R*R);   // ratio dV/dR
            
        } else if (R > rR && R <= rc) {
            // Total mass accumulation on ice sheet minus runoff
            Btot = M_PI * P * R * R
            - M_PI * beta * (hr - b0 + slope * R) * (R*R - rR*rR)
            - (4.0 * M_PI * beta * pow(mu, 0.5))/5.0 * pow(R - rR, 2.5)
            + (4.0 * M_PI * beta * pow(mu, 0.5))/3.0 * R * pow(R-rR, 1.5);
            // In case there is no marine ice sheet / grounding line
            F = 0;        // no ice flux
            ISO = 0;      // (third term equation 14) NAME?
            fac = M_PI * (1.0 + eps1) * (4.0/3.0 * pow(mu, 0.5) * pow(R, 1.5) - slope*R*R);   // ratio dV/dR
            
        } else if (R <= rR && R >= rc) {
            // Total mass accumulation with marine ice sheet / grounding line
            Btot = M_PI * P * R * R;
            Hw = slope * R - b0 + SL(i);  // (equation 10)
            F = 2.0 * M_PI * R * f * del * pow(Hw, Gamma + 1.0);   // Ice flux (equation 9)
            ISO = 2.0 * M_PI * eps2 * (slope * rc * rc - (b0/slope) * rc) * GSL(i);  // third term equation 14 !! NAME?
            fac = M_PI * (1.0 + eps1) * (4.0/3.0 * pow(mu, 0.5) * pow(R, 1.5) - slope*R*R)
            - 2.0 * M_PI * eps2 * (slope * R * R - b0 * R);
            
        } else {
            // Total mass accumulation minus runoff with marine ice sheet / grounding line
            Btot = M_PI * P * R * R
            - M_PI * beta * (hr - b0 + slope * R) * (R*R - rR*rR)
            - (4.0 * M_PI * beta * pow(mu, 0.5))/5.0 * pow(R - rR, 2.5) + (4.0 * M_PI * beta * pow(mu, 0.5))/3.0 * R * pow(R-rR, 1.5);
            Hw = slope * R - b0 + SL(i);  // (equation 10)
            F = 2.0 * M_PI * R * f * del * pow(Hw, Gamma + 1.0);   // Ice flux (equation 9)
            ISO = 2.0 * M_PI * eps2 * (slope * rc * rc - (b0/slope) * rc) * GSL(i);  // third term equation 14 !! NAME?
            fac = M_PI * (1.0 + eps1) * (4.0/3.0 * pow(mu, 0.5) * pow(R, 1.5) - slope*R*R)
            - 2.0 * M_PI * eps2 * (slope * R * R - b0 * R);
        }
        dR = (Btot - F + ISO)/fac;
        R = R + dR;                // Esitmate new radius
        V = 8.0/15.0 * M_PI * pow(mu, 0.5) * pow(R, 2.5) - 1.0/3.0 * M_PI * slope * pow(R, 3.0);
        //V =  M_PI * (8.0/15.0 * pow(mu, 0.5) * pow(R, 2.5) - 1.0/3.0 * slope * pow(R, 3.0));
        
        // Calculate sea volume
        Vsea = M_PI * (2.0/3.0 * slope * (pow(R, 3.0) - pow(rc, 3.0))
                       - b0 * (pow(R, 2.0) - pow(rc, 2.0)));

        // Calulate the volume change over time
        if(R <= rc){
            Volt = (1.0 + eps1) * V;
        }
        else{
            Volt = (1.0 + eps1) * V - eps2 * Vsea;
        }
        
        // Ice sheet volume (equation 13)
        Rad(i)  = R;
        
        Vais(i) = Volt;
        SLE(i)  = 57.0 * (1.0 - Vais(i)/Volo);  // Takes steady state present day volume to correspond to 57m SLE
    }
    
    return R_NilValue;
}


SEXP daisOdeC(SEXP gparms)
{
    SEXP rc;

    daisInit(gparms);
    rc = daisOdeCdeSolve();

    return rc;
}
