// Copyright 2009, 2010 Robert W. Fuller <hydrologiccycle@gmail.com>
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
// slrmodel.c
// written by Robert W. Fuller on 090809
//
// check for warnings:  gcc -Wall -Wextra -Wc++-compat -c -I/usr/lib/R/include slrmodel.c && rm slrmodel.o
// build from R:  system("rm slrmodel.o r.o; R CMD SHLIB slrmodel.c r.c"); dynReload("slrmodel", "grinstedOdeInit")
//

#include "r.h"


static RVector modelParms;

static RIntVector samples;
#define gmstCol (getIntVectorElem(&samples, 1))

static RMatrix forcings[5];
#define gmstMat     (forcings[0])
#define massGmstMat (forcings[1])  // TODO:  move this to extraParms?

// extra parameters:  constants or parameters that are temporarily fixed
static RVector extraParms;

// logicals to enable/disable parts of models
static RIntVector switches;

static RParm parms[] = {
    { "mp",     &modelParms,    RTYPE_VECTOR,       0 },
    { "frc",    forcings,       RTYPE_MATRICES,     NELEMS(forcings) },
    { "spl",    &samples,       RTYPE_INT_VECTOR,   0 },
    { "ep",     &extraParms,    RTYPE_VECTOR,       0 },
    { "sw",     &switches,      RTYPE_INT_VECTOR,   0 }
};


static double mp_s0, mp_a, mp_b, mp_tau, mp_a2, mp_b2, mp_c, mp_d, mp_tau2, max_sle, gis_scl, gis_s0, gis_temp;

static RNamedReal realParms[] = {
    { "s0",         &mp_s0 },
    { "a",          &mp_a },
    { "b",          &mp_b },
    { "tau",        &mp_tau },
    { "c",          &mp_c },
    { "d",          &mp_d },
    { "tau2",       &mp_tau2 },
    { "max_sle",    &max_sle },
    { "gis_scl",    &gis_scl },
    { "gis_s0",     &gis_s0 },
    { "gis_temp",   &gis_temp },
    { "a2",         &mp_a2 },
    { "b2",         &mp_b2 }
};


static int sw_tide_ts, sw_rignot_ts, sw_tau2_prior, sw_alley_prior, sw_log_tau, sw_log_tau2, sw_old_ref;

static RNamedInt swParms[] = {
    { "tide_ts",        &sw_tide_ts },
    { "rignot_ts",      &sw_rignot_ts },
    { "tau2_prior",     &sw_tau2_prior },
    { "alley_prior",    &sw_alley_prior },
    { "log_tau",        &sw_log_tau },
    { "log_tau2",       &sw_log_tau2 },
    { "old_ref",        &sw_old_ref }
};


static DL_FUNC get_deSolve_gparms;

void R_init_slrmodel(DllInfo *dll)
{
    get_deSolve_gparms = R_GetCCallable("deSolve", "get_deSolve_gparms");

    sortNamedStructs(realParms);
    sortNamedStructs(parms);
    sortNamedStructs(swParms);
}


void grinstedOdeInit(void (*odeparms)(int *, double *))
{
    SEXP gparms;

    gparms = (SEXP) get_deSolve_gparms();

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
        //Rprintf("log_tau is %d\n", sw_log_tau);
    }
}


void grinstedOdeC(
    int *neq, double *t, double *y, double *ydot,
    double *yout, int *ip)
{
    double s, temp, ds;

    s = y[0];
    temp = tsFindByDate(&gmstMat, *t, gmstCol);
    ds = (mp_a * temp + mp_b - s) / mp_tau;

    ydot[0] = ds;
}


void gringisFitOdeC(
    int *neq, double *t, double *y, double *ydot,
    double *yout, int *ip)
{
    double s, temp, ds;

    s = y[0];
    temp = tsFindByDate(&gmstMat, *t, gmstCol);
    ds = (mp_a2 * temp + mp_b2 - s) / mp_tau2;

    ydot[0] = ds;

    if (1 != ip[0]) {
        error("slrmodel.c expects ode(nout=1) for gringisFitOdeC()");
    }
    yout[0] = ds;
}


void gringisOdeC(
    int *neq, double *t, double *y, double *ydot,
    double *yout, int *ip)
{
    double s_other, s_gis, temp, ds_other, ds_gis, ds_total;

    //s_total = y[0];
    s_other   = y[1];
    s_gis     = y[2];

    temp = tsFindByDate(&gmstMat, *t, gmstCol);

    ds_other = (mp_a  * temp + mp_b  - s_other) / mp_tau;
    ds_gis   = (mp_a2 * temp + mp_b2 - s_gis)   / mp_tau2;
    ds_total = ds_other + ds_gis;

    ydot[0] = ds_total;
    ydot[1] = ds_other;
    ydot[2] = ds_gis;

    if (1 != ip[0]) {
        error("slrmodel.c expects ode(nout=1) for gringisOdeC()");
    }
    yout[0] = ds_gis;
}


#if 0

// TODO:  note that it is now broken to assume that gmstCol will give
// corresponding columns from gmstMat and massGmstMat;  this is due to
// the order of magnitude speed up achieved in assim.R::runPredict()
// by dicing up these matrices and always setting the column to 2;
// using an offset such as gis_temp is a better alternative anyway
//

void quadgisOdeC(
    int *neq, double *t, double *y, double *ydot,
    double *yout, int *ip)
{
    double s, temp, masstemp, dm, ds;

    s = y[0];
    temp     = tsFindByDate(    &gmstMat, *t, gmstCol);
    masstemp = tsFindByDate(&massGmstMat, *t, gmstCol);

    //dm = (mp_c * masstemp + mp_d * masstemp * masstemp) / -360 / 1000;
    dm = masstemp * (mp_c + mp_d * masstemp) / -360 / 1000;
    ds = (mp_a * temp + mp_b - s) / mp_tau + dm;

    ydot[0] = ds;
}


void qdgrgisOdeC(
    int *neq, double *t, double *y, double *ydot,
    double *yout, int *ip)
{
    double s_other, s_gis, temp, masstemp, ds_other, ds_gis, ds_total, max_gis;

    //s_total = y[0];
    s_other   = y[1];
    s_gis     = y[2];
    max_gis   = max_sle;
    if (s_gis > max_gis) {
        s_gis = max_gis;
    }

    temp     = tsFindByDate(    &gmstMat, *t, gmstCol);
    masstemp = tsFindByDate(&massGmstMat, *t, gmstCol);

    ds_other = (mp_a * temp + mp_b                  - s_other) / mp_tau;
    ds_gis   = (masstemp * (mp_c + mp_d * masstemp) - s_gis)   / mp_tau2;
    ds_total = ds_other + ds_gis;

    ydot[0] = ds_total;
    ydot[1] = ds_other;
    ydot[2] = ds_gis;

    if (1 != ip[0]) {
        error("slrmodel.c expects ode(nout=1) for gringisOdeC()");
    }
    yout[0] = ds_gis;
}


void qdgrgisFitOdeC(
    int *neq, double *t, double *y, double *ydot,
    double *yout, int *ip)
{
    double s, temp, ds, max_gis;

    s = y[0];
    max_gis = max_sle;
    if (s > max_gis) {
        s = max_gis;
    }

    temp = tsFindByDate(&gmstMat, *t, gmstCol);
    ds = (temp * (mp_c + mp_d * temp) - s) / mp_tau2;

    ydot[0] = ds;

    if (1 != ip[0]) {
        error("slrmodel.c expects ode(nout=1) for gringisFitOdeC()");
    }
    yout[0] = ds;
}

#endif


void allgrgisOdeC(
    int *neq, double *t, double *y, double *ydot,
    double *yout, int *ip)
{
    double s_other, s_gis, temp, masstemp, ds_other, ds_gis, ds_total, max_gis;

    //s_total = y[0];
    s_other   = y[1];
    s_gis     = y[2];
    max_gis   = max_sle;
    if (s_gis > max_gis) {
        s_gis = max_gis;
    }

    temp         = tsFindByDate(    &gmstMat, *t, gmstCol);
    if (sw_old_ref) {
        masstemp = tsFindByDate(&massGmstMat, *t, gmstCol) * gis_scl;
    } else {
        masstemp = (temp + gis_temp) * gis_scl;
    }

    ds_other = (mp_a * temp + mp_b                  - s_other) / mp_tau;

    if (masstemp <= 0) {
        ds_gis   = (masstemp *  mp_c                    - s_gis)   / mp_tau2;
    } else {
        ds_gis   = (masstemp * (mp_c + mp_d * masstemp) - s_gis)   / mp_tau2;
    }
    ds_total = ds_other + ds_gis;

    ydot[0] = ds_total;
    ydot[1] = ds_other;
    ydot[2] = ds_gis;

    if (2 != ip[0]) {
        error("slrmodel.c expects ode(nout=2) for gringisOdeC()");
    }
    yout[0] = ds_gis;
    yout[1] = ds_total;
}


void allgrgisFitOdeC(
    int *neq, double *t, double *y, double *ydot,
    double *yout, int *ip)
{
    double s, temp, ds, max_gis;

    s = y[0];
    max_gis = max_sle;
    if (s > max_gis) {
        s = max_gis;
    }

    if (sw_old_ref) {
        // could get the same effect be setting gis_temp to zero
        temp =  tsFindByDate(&gmstMat, *t, gmstCol) * gis_scl;
    } else {
        temp = (tsFindByDate(&gmstMat, *t, gmstCol) + gis_temp) * gis_scl;
    }
    if (temp <= 0) {
        ds = (temp *  mp_c                - s) / mp_tau2;
    } else {
        ds = (temp * (mp_c + mp_d * temp) - s) / mp_tau2;
    }

    ydot[0] = ds;

    if (1 != ip[0]) {
        error("slrmodel.c expects ode(nout=1) for gringisFitOdeC()");
    }
    yout[0] = ds;
}
