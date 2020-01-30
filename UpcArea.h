/*
 *  UpcArea.h
 *  UnipodC
 *
 *  Created by ionosph on 2010/10/22.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _UPCAREA_H
#define _UPCAREA_H

#include "UpcBase.h"

typedef struct _UpcArea UpcArea;
struct _UpcArea {
    int _refCount;
    void (*_dealloc)(UpcArea *);
    double center[2];
    double semiDia[2];
    double upper[2];
    double lower[2];
    UpcArea *(*copy)(UpcArea *);
    void (*print)(const UpcArea *, int, FILE *);
    void (*toXml)(const UpcArea *, int, FILE *);
    UpcErrStatus (*lset)(UpcArea *);
};

UpcArea *UpcArea_init(void);
Bool UpcArea_isActive(const UpcArea *self, double x, double y);
double UpcArea_semiDia(const UpcArea *self);
void UpcArea_setSemiDia(UpcArea *self, double sd);
double UpcArea_dia(const UpcArea *self);
void UpcArea_setDia(UpcArea *self, double dia);
Bool UpcArea_xSymmetric(const UpcArea *self);
Bool UpcArea_ySymmetric(const UpcArea *self);
Bool UpcArea_rSymmetric(const UpcArea *self);

#define UPC_INFSD 1.0e30

#endif
