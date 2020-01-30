/*
 *  UpcSurf.h
 *  UnipodC
 *
 *  Created by ionosph on 2011/04/01.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _UPCSURF_H
#define _UPCSURF_H

#include "UpcBase.h"
#include "UpcShapeSph.h"
#include "UpcShapeAsph.h"
#include "UpcShapeXYP.h"
#include "UpcArea.h"
#include "UpcDecenter.h"
#include "UpcGlass.h"
#include "UpcCoord.h"


typedef enum _UpcSurfRmode {
    RmodeRefractive,
    RmodeReflective
} UpcSurfRmode;

typedef struct _UpcSurf UpcSurf;
struct _UpcSurf {
    int _refCount;
    void (*_dealloc)(UpcSurf *);
    UpcShape *shape;
    UpcGlass *medium; // 弱参照
    UpcArea *area;
    UpcDecenter *dec;
    UpcCoord *coord;
    UpcSurfRmode mode;
    double d;
    double *_n;
};

UpcErrStatus UpcSurf_lset(UpcSurf *self);
char *UpcSurfRmode_toString(UpcSurfRmode i);
UpcSurfRmode UpcSurfRmode_fromString(const char *s);
UpcSurf *UpcSurf_init(void);
UpcSurf *UpcSurf_copy(UpcSurf *self);
void UpcSurf_print(const UpcSurf *self, int isur, int indent, FILE *ostream);
void UpcSurf_toXml(const UpcSurf *self, int isur, int indent, FILE *ostream);
double UpcSurf_e(UpcSurf *self, int iwav);
Bool UpcSurf_xSymmetric(const UpcSurf *self);
Bool UpcSurf_ySymmetric(const UpcSurf *self);
Bool UpcSurf_rSymmetric(const UpcSurf *self);

UpcErrStatus UpcSurf_setShapeTypeSph(UpcSurf *self);
UpcErrStatus UpcSurf_setShapeTypeAsph(UpcSurf *self, int max_order);
UpcErrStatus UpcSurf_setShapeTypeXYP(UpcSurf *self, int max_order, int max_Xorder, int max_Yorder);
UpcErrStatus UpcSurf_setShapeTypeZernike(UpcSurf *self, int max_orderR, int zcn_max);


#endif
