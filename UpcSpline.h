/*
 *  UpcSpline.h
 *  UnipodC
 *
 *  Created by ionosph on 11/07/06.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _UPCSPLINE_H
#define _UPCSPLINE_H


#include "UpcBase.h"

typedef struct _UpcSpline UpcSpline;
struct _UpcSpline {
    int _refCount;
    void (*_dealloc)(UpcSpline *);
    int ndata;
    double *xa;
    double *ya;
    double *d2ya;
    int iprev;
};

UpcSpline *UpcSpline_init(int n, const double *xa, const double *ya, double dy1, double dyN);
double UpcSpline_y(UpcSpline *self, double x);

#endif

