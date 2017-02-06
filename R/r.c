// Copyright (C) 2009, 2010 Robert W. Fuller <hydrologiccycle@gmail.com>
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
// r.c
// written by Robert W. Fuller on 091014
// check for warnings:  gcc -Wall -Wextra -Wc++-compat -c -I/usr/lib/R/include r.c && rm r.o
//

#include "r.h"


void initIntVector(RIntVector *v, SEXP s_vec)
{
    double *d_vec;
    int i;

    initCommon(&v->comn, s_vec);
    if (isReal(s_vec)) {
        d_vec = REAL(s_vec);

        if (SNELEMS(v->preAlloc) < v->comn.length) {
            v->comn.int_arr = (int *) R_alloc(v->comn.length, sizeof(int));
            if (!v->comn.int_arr) {
                error("out of memory in initIntVector()");
            }
        } else {
            v->comn.int_arr = v->preAlloc;
        }

        for (i = 0;  i < v->comn.length;  ++i) {
            v->comn.int_arr[i] = (int) d_vec[i];
        }
    } else {

        // this works for LOGICAL as well
        v->comn.int_arr = INTEGER(s_vec);
    }
}


void initMatrix(RMatrix *m, SEXP s_mat)
{
    SEXP rdim;
    int *dims;

    initCommon(&m->comn, s_mat);
    m->comn.dbl_arr = REAL(s_mat);

    rdim = getAttrib(s_mat, R_DimSymbol);
    dims = INTEGER(rdim);
    if (2 != LENGTH(rdim)) {
        error("not a matrix in initMatrix()");
    }
    if (dims[0] * dims[1] != m->comn.length) {
        error("length of matrix does not match dimensions in initMatrix()");
    }
    m->rows = dims[0];
    m->cols = dims[1];

    resetMatrixFindHint(m);
}


typedef void (*initFun)(void *, SEXP);


int initFromList(void *elem, initFun fn, int elemSize, int elemCount, SEXP s_list)
{
    int len, i;

    len = LENGTH(s_list);
    if (len > elemCount) {
        error("too many list elements in initFromList()");
    }

    for (i = 0;  i < len;  ++i) {
        fn((char *) elem + i * elemSize, VECTOR_ELT(s_list, i));
    }

    return len;
}


int cmpNamedStruct(const void *a, const void *b)
{
    return strcmp(((NamedStruct *) a)->name, ((NamedStruct *) b)->name);
}


int initNamedReals(RNamedReal *parms, int count, RVector *v)
{
    NamedStruct key;
    RNamedReal *match;
    SEXP s_names;
    int i, namesLen;

    s_names = getAttrib(v->comn.s_ptr, R_NamesSymbol);
    namesLen = LENGTH(s_names);
    for (i = 0;  i < namesLen;  ++i) {

        key.name = CHAR(STRING_ELT(s_names, i));
        if (!key.name[0]) {
            error("all arguments must be named in initNamedReals()");
        }

        match = (RNamedReal *) bsearch(&key, parms, count, sizeof(RNamedReal), cmpNamedStruct);
        if (match) {
            *(match->elem) = getVectorElem(v, i + 1);
            continue;
        }

        error("no match found for argument %s in initNamedReals()", key.name);
    }

    return namesLen;
}


int initNamedInts(RNamedInt *parms, int count, RIntVector *v)
{
    NamedStruct key;
    RNamedInt *match;
    SEXP s_names;
    int i, namesLen;

    s_names = getAttrib(v->comn.s_ptr, R_NamesSymbol);
    namesLen = LENGTH(s_names);
    for (i = 0;  i < namesLen;  ++i) {

        key.name = CHAR(STRING_ELT(s_names, i));
        if (!key.name[0]) {
            error("all arguments must be named in initNamedInts()");
        }

        match = (RNamedInt *) bsearch(&key, parms, count, sizeof(RNamedInt), cmpNamedStruct);
        if (match) {
            *(match->elem) = getIntVectorElem(v, i + 1);
            continue;
        }

        error("no match found for argument %s in initNamedInts()", key.name);
    }

    return namesLen;
}


int initParms(RParm *parms, int count, SEXP s_list)
{
    NamedStruct key;
    RParm *match;
    void *elem;
    SEXP s_names, s_elem;
    int i, namesLen;

    s_names = getAttrib(s_list, R_NamesSymbol);
    namesLen = LENGTH(s_names);
    for (i = 0;  i < namesLen;  ++i) {

        key.name = CHAR(STRING_ELT(s_names, i));
        if (!key.name[0]) {
            error("all arguments must be named in initParms()");
        }

        match = (RParm *) bsearch(&key, parms, count, sizeof(RParm), cmpNamedStruct);
        if (match) {

            s_elem = VECTOR_ELT(s_list, i);
            elem = match->elem; 
            switch (match->elemType) {

                case RTYPE_MATRIX:
                    initMatrix((RMatrix *) elem, s_elem);
                    break;

                case RTYPE_VECTOR:
                    initVector((RVector *) elem, s_elem);
                    break;

                case RTYPE_INT_VECTOR:
                    initIntVector((RIntVector *) elem, s_elem);
                    break;

                case RTYPE_MATRICES:
                    initFromList(elem, (initFun) initMatrix, sizeof(RMatrix), match->count, s_elem);
                    break;

                case RTYPE_VECTORS:
                    initFromList(elem, (initFun) initVector, sizeof(RVector), match->count, s_elem);
                    break;

                case RTYPE_INT_VECTORS:
                    initFromList(elem, (initFun) initIntVector, sizeof(RIntVector), match->count, s_elem);
                    break;

                case RTYPE_FUNCTION:
                    ((RParmFun) elem)(match, s_elem);
                    break;

                default:
                    error("invalid type %d in initParms()", match->elemType);
                    break;
            }

            continue;
        }

        error("no match found for argument %s in initParms()", key.name);
    }

    return namesLen;
}
