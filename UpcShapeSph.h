/*
 *  UpcShapeSph.h
 *  UnipodC
 *
 *  Created by ionosph on 2010/10/16.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef _UPCSHAPESPH_H
#define _UPCSHAPESPH_H

#include "UpcBase.h"
//#include "UpcArea.h"

struct _UpcArea;

Bool UpcShapeSph_conicZ(double curv, double koni, double x, double y, int mode, double *z, double *zx, double *zy, double *zxx, double *zxy, double *zyy);
int UpcShapeSph_solveQuadEq(double a, double b, double c, double *k1, double *k2);
int UpcShapeSph_conicIntersection(double curv, double koni, const double *p, const double *d, double *p1, double *p2);

/* 面タイプ
 *     継承関係にある面タイプは割り切れるようになっている.
 *     すなわち, 子タイプ % 親タイプ == 0 が成り立つ.
 */
enum UpcShapeType {
    UpcShapeType_sph     = 2,
    UpcShapeType_asph    = 2 * 3,
    UpcShapeType_xyp     = 2 * 3 * 5,
    UpcShapeType_zernike = 2 * 3 * 7
};

typedef struct _UpcShape UpcShape;
struct _UpcShape {
    int _refCount;
    void (*_dealloc)(UpcShape *self);
    // attribute(sph)
    enum UpcShapeType shapeType;
    double curv;
    // attribute(asph)
    double coni;
    int romax;
    double *rpc;     // rpc[0..romax]
    double nrmR2inv; // = 1 / (nrmR)^2
    int itermax;
    double tol;
    // attribute(xyp)
    int xomax;
    int yomax;
    double **xypc;   // xypc[0..yomax][0..xomax]
    double nrmXinv;  // = 1 / nrmX
    double nrmYinv;  // = 1 / nrmY
    // attribute(zernike)
    int zcnmax;
    double *zc;
    double nrmZRinv;
    // method
    double (*axialPower)(const UpcShape *self);
    Bool (*getZ)(const UpcShape *self, double x, double y, double *z);
    Bool (*getCurv)(const UpcShape *self, double y, double *curvMeri, double *curvSagi);
    void (*getRBC)(const UpcShape *self, double *rv, double *bv, double *cv);
    Bool (*getNormalVector)(const UpcShape *self, double x, double y, double *nv);
    Bool (*getIntersection)(const UpcShape *self, const double *pos, const double *dir, struct _UpcArea *area, double *p, double *n, Bool *isActive);
    Bool (*xSymmetric)(const UpcShape *self);
    Bool (*ySymmetric)(const UpcShape *self);
    Bool (*rSymmetric)(const UpcShape *self);
    UpcShape *(*copy)(const UpcShape *self);
    void (*print)(const UpcShape *, int, FILE *);
    void (*toXml)(const UpcShape *, int, FILE *);
    UpcErrStatus (*lset)(UpcShape *);
};

#define UpcShapeSph_RADI(self) (((self)->curv)?(1.0/((self)->curv)):0.0)

char *UpcShapeType_string(enum UpcShapeType i);
UpcShape *UpcShapeSph_init(void);
UpcShape *UpcShapeSph_initWithShape(UpcShape *other);

/* 以下の関数は, 子クラスからも使用可 */
double UpcShapeSph_r(const UpcShape *self);
void UpcShapeSph_setR(UpcShape *self, double r);

#endif
