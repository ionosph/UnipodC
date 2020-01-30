/*
 *  UpcShapeXYP.h
 *  UnipodC
 *
 *  Created by ionosph on 2011/04/28.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _UPCSHAPEXYP_H
#define _UPCSHAPEXYP_H

#include "UpcBase.h"
#include "UpcShapeSph.h"
#include "UpcShapeAsph.h"

UpcShape *UpcShapeXYP_init(int max_order, int max_Xorder, int max_Yorder);
UpcShape *UpcShapeXYP_initWithShape(int max_orderR, int max_orderX, int max_orderY, UpcShape *other);

double UpcShapeXYP_nrmX(const UpcShape *self);
double UpcShapeXYP_nrmY(const UpcShape *self);
Bool UpcShapeXYP_setNrmX(UpcShape *self, double nrmX);
Bool UpcShapeXYP_setNrmY(UpcShape *self, double nrmY);
void UpcShapeXYP_clearXYPC(UpcShape *self);
Bool UpcShapeXYP_addZcoef(UpcShape *self, int cno, double coef);

Bool UpcShapeXYP_xypolyZ(int max_orderX, int max_orderY, double **xypc, double x, double y, double nrmXinv, double nrmYinv, int mode, double *z, double *zx, double *zy, double *zxx, double *zxy, double *zyy);

UpcShape *UpcShapeZernike_init(int max_Rorder, int zcn_max);
UpcShape *UpcShapeZernike_initWithShape(int max_orderR, int zcn_max, UpcShape *other);

double UpcShapeZernike_nrmZR(const UpcShape *self);
Bool UpcShapeZernike_setNrmZR(UpcShape *self, double nrmZR);




#endif

