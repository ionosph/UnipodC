/*
 *  UpcGlass.h
 *  UnipodC
 *
 *  Created by ionosph on 11/07/07.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _UPCGLASS_H
#define _UPCGLASS_H

#include "UpcBase.h"

struct _UpcDictDD;
struct _UpcSpline;

typedef struct _UpcGlass UpcGlass;
struct _UpcGlass {
    int _refCount;
    void (*_dealloc)(UpcGlass *);
    char *name;
    struct _UpcDictDD *cash;
    double coef[6];
    struct _UpcSpline *splineObj;
    double (*n)(UpcGlass *, double);
    int nwav;
    double *indexF;
    double *indexB;
    void (*toXml)(const UpcGlass *, int, FILE *);
    UpcErrStatus (*lset)(UpcGlass *);
};

UpcGlass *UpcGlassSellmeier_init(const char *name, double A1, double A2, double A3, double B1, double B2, double B3);
UpcGlass *UpcGlassNdVd_init(const char *name, double nd, double vd);
UpcGlass *UpcGlassConstant_init(const char *name, double ind);
UpcGlass *UpcGlassSpline_init(const char *name, int ndata, const double *wls, const double *inds);
Bool UpcGlass_setIndex(UpcGlass *self, int nwav, const double *wls);

#endif
