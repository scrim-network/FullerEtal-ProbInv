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
// r.h
// written by Robert W. Fuller on 090809
// updated 091010
// check for warnings:  gcc -Wall -Wextra -Wc++-compat -c -I/usr/lib/R/include r.h && rm r.h.gch
//

#ifndef R_H_INCLUDED
#define R_H_INCLUDED

// need this for ssize_t on hammer
#include <sys/types.h>
#include <stdlib.h> // qsort

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

struct _RIntVector;
struct _RMatrix;
struct _RParm;
struct _RType;
struct _RNamedReal;
struct _RVector;
struct _RNamedInt;


void initIntVector(struct _RIntVector *v, SEXP s_vec);
void initMatrix(struct _RMatrix *m, SEXP s_mat);
int initParms(struct _RParm *parms, int count, SEXP s_list);
int initNamedReals(struct _RNamedReal *parms, int count, struct _RVector *v);
int initNamedInts(struct _RNamedInt *parms, int count, struct _RIntVector *v);
int cmpNamedStruct(const void *a, const void *b);

#define  NELEMS(vector) (sizeof(vector) / sizeof(vector[0]))
#define SNELEMS(vector) ((ssize_t) NELEMS(vector))

#define sortNamedStructs(x) (qsort(x, NELEMS(x), sizeof(x[0]), cmpNamedStruct))


typedef enum _RParmType {

    RTYPE_MATRIX,
    RTYPE_VECTOR,
    RTYPE_INT_VECTOR,
    RTYPE_MATRICES,
    RTYPE_VECTORS,
    RTYPE_INT_VECTORS,
    RTYPE_FUNCTION

} RParmType;


typedef void (*RParmFun)(struct _RParm *, SEXP);


typedef struct _NamedStruct {

    const char *name;

} NamedStruct;


typedef struct _RParm {

    const char *name;
    void *elem;
    RParmType elemType;
    int count;

} RParm;


typedef struct _RNamedReal {

    const char *name;
    double *elem;

} RNamedReal;


typedef struct _RNamedInt {

    const char *name;
    int *elem;

} RNamedInt;


typedef struct _RCommon {

    SEXP s_ptr;
    union {
        double *dbl_arr;
        int    *int_arr;
    };
    int length;

    // on 64-bit, there are 4 unused bytes here

} RCommon;


typedef struct _RMatrix {

    RCommon comn;
    int rows, cols;
    int findHint;

} RMatrix;


typedef struct _RIntVector {

    RCommon comn;
    int preAlloc[5];

} RIntVector;


typedef struct _RVector {

    RCommon comn;

} RVector;


static inline void initCommon(RCommon *comn, SEXP s_ptr)
{
    comn->s_ptr  = s_ptr;
    comn->length = LENGTH(s_ptr);
}


static inline double *getMatrixColumn(RMatrix *m, int col)
{
    int index = (col - 1) * m->rows;

    if (col < 1 || col > m->cols) {
        error("illegal index in getMatrixColumn()");
    }

    return &m->comn.dbl_arr[index];
}


static inline double getMatrixElem(RMatrix *m, int col, int row)
{
    if (row < 1 || row > m->rows) {
        error("illegal index in getMatrixElem()");
    }

    return getMatrixColumn(m, col)[row - 1];
}


static inline void resetMatrixFindHint(RMatrix *m)
{
    m->findHint = 0;
}


static inline double tsFindByDate(RMatrix *m, double time, int col)
{
    int mflag;

    m->findHint = findInterval(m->comn.dbl_arr, m->rows, time, FALSE, FALSE, m->findHint, &mflag);
    if (time != m->comn.dbl_arr[m->findHint - 1]) {
        error("did not find time %g in tsFindByDate()", time);
    }

    // should check mflag here

    return getMatrixElem(m, col, m->findHint);
}


static inline void initVector(RVector *v, SEXP s_vec)
{
    initCommon(&v->comn, s_vec);
    v->comn.dbl_arr = REAL(s_vec);
}


static inline double getVectorElem(RVector *v, int index)
{
    if (index < 1 || index > v->comn.length) {
        error("illegal index in getVectorElem()");
    }

    return v->comn.dbl_arr[index - 1];
}


static inline int getIntVectorElem(RIntVector *v, int index)
{
    if (index < 1 || index > v->comn.length) {
        error("illegal index in getIntVectorElem()");
    }

    return v->comn.int_arr[index - 1];
}


#endif // R_H_INCLUDED
