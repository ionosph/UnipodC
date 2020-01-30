/*
 *  UpcShapeAsph.h
 *  UnipodC
 *
 *  Created by ionosph on 2010/10/20.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _UPCSHAPEASPH_H
#define _UPCSHAPEASPH_H

#include "UpcBase.h"
#include "UpcShapeSph.h"

Bool UpcShapeAsph_rpolyZ(int max_order, const double *ac, double x, double y, double nrmR2inv, int mode, double *z, double *zx, double *zy, double *zxx, double *zxy, double *zyy);

#define UpcShapeAsph_NRMR(self) (1.0/sqrt((self)->nrmR2inv))

UpcShape *UpcShapeAsph_init(int max_order);
UpcShape *UpcShapeAsph_initWithShape(int max_order, UpcShape *other);

/* 以下の関数は子クラスからも使用可 */
double UpcShapeAsph_nrmR(const UpcShape *self);
Bool UpcShapeAsph_setNrmR(UpcShape *self, double nrmR);

#endif

