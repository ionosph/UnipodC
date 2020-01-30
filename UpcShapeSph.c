/*
 *  UpcShapeSph.c
 *  UnipodC
 *
 *  Created by ionosph on 2010/10/16.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcShapeSph.h"
#include "UpcVector.h"
#include "UpcArea.h"

static const double curvDefault = 0.0;

static double axialPower(const UpcShape *self);
static Bool getZ(const UpcShape *self, double x, double y, double *z);
static Bool getCurv(const UpcShape *self, double y, double *curvMeri, double *curvSagi);
static void getRBC(const UpcShape *self, double *rv, double *bv, double *cv);
static Bool getNormalVector(const UpcShape *self, double x, double y, double *nv);
static Bool getIntersection(const UpcShape *self, const double *pos, const double *dir, UpcArea *area, double *p, double *n, Bool *isActive);
static Bool xSymmetric(const UpcShape *self);
static Bool ySymmetric(const UpcShape *self);
static Bool rSymmetric(const UpcShape *self);


char *UpcShapeType_string(enum UpcShapeType i)
{
    switch (i) {
        case UpcShapeType_sph:
            return "sph";
        case UpcShapeType_asph:
            return "asph";
        case UpcShapeType_xyp:
            return "xyp";
        case UpcShapeType_zernike:
            return "zernike";
        default:
            return "unknown shape type";
    }
}

static void UpcShapeSph_dealloc(UpcShape *self)
{
    free(self);
}

static void UpcShapeSph_print(const UpcShape *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);

    fprintf(ostream, "%stype = %s\n", space, UpcShapeType_string(self->shapeType));
    fprintf(ostream, "%sr    = %13.6e\n", space, UpcShapeSph_RADI(self));
}

static void UpcShapeSph_toXml(const UpcShape *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);

    fprintf(ostream, "%s<shape type=\"%s\">\n", space, UpcShapeType_string(self->shapeType));
    if (self->curv != curvDefault)
        fprintf(ostream, "%s  <r>%23.15e</r>\n", space, UpcShapeSph_RADI(self));
    fprintf(ostream, "%s</shape>\n", space);
}

static UpcErrStatus UpcShapeSph_lset(UpcShape *self)
{
    return UpcE_NoIssues;
}

static UpcShape *UpcShapeSph_copy(const UpcShape *self)
{
    UpcShape *t = UpcShapeSph_init();

    if (t) {
        t->curv = self->curv;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeSph_copy");
    }
    return t;
}

UpcShape *UpcShapeSph_init(void)
{
    UpcShape *self = (UpcShape *)calloc(1, sizeof(UpcShape));

    if (self) {
        self->_refCount = 1;
        self->_dealloc = UpcShapeSph_dealloc;
        self->shapeType = UpcShapeType_sph;
        self->curv = curvDefault;
        self->axialPower = axialPower;
        self->getZ = getZ;
        self->getCurv = getCurv;
        self->getRBC = getRBC;
        self->getNormalVector = getNormalVector;
        self->getIntersection = getIntersection;
        self->xSymmetric = xSymmetric;
        self->ySymmetric = ySymmetric;
        self->rSymmetric = rSymmetric;
        self->copy = UpcShapeSph_copy;
        self->print = UpcShapeSph_print;
        self->toXml = UpcShapeSph_toXml;
        self->lset = UpcShapeSph_lset;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeSph_init");
    }
    return self;
}

UpcShape *UpcShapeSph_initWithShape(UpcShape *other)
{
    UpcShape *self = UpcShapeSph_init();
    enum UpcShapeType otype = other->shapeType;

    if (self) {
        if (otype % UpcShapeType_sph == 0) { // UpcShapeSphクラスまたは子クラスであれば
            self->curv = other->curv;
        }
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeSph_initWithShape");
    }
    return self;
}

double UpcShapeSph_r(const UpcShape *self)
{
    return self->curv ? 1.0 / self->curv : 0.0;
}

void UpcShapeSph_setR(UpcShape *self, double r)
{
    self->curv = r ? 1.0 / r : 0.0;
}

static double axialPower(const UpcShape *self)
{
    return self->curv;
}

static Bool getZ(const UpcShape *self, double x, double y, double *z)
{
    double dmy;

    return UpcShapeSph_conicZ(self->curv, 0.0, x, y, 0, z, &dmy, &dmy, &dmy, &dmy, &dmy);
}

static Bool getCurv(const UpcShape *self, double y, double *curvMeri, double *curvSagi)
{
    double curv = self->curv;

    if (curv) {
        if (fabs(y) <= fabs(1.0 / curv)) {
            *curvMeri = curv;
            *curvSagi = curv;
            return TRUE;
        }
        return FALSE;
    }
    else {
        *curvMeri = curv;
        *curvSagi = curv;
        return TRUE;
    }
}

static void getRBC(const UpcShape *self, double *rv, double *bv, double *cv)
{
    *rv = self->curv;
    *bv = 0.0;
    *cv = 0.0;
}

static Bool getNormalVector(const UpcShape *self, double x, double y, double *nv)
{
    double zx, zy, dmy, ilen;
    
    if (UpcShapeSph_conicZ(self->curv, 0.0, x, y, 1, &dmy, &zx, &zy, &dmy, &dmy, &dmy)) {
        ilen = 1.0 / sqrt(zx * zx + zy * zy + 1.0);
        UpcVector_SET(nv, -ilen * zx, -ilen * zy, ilen);
        return TRUE;
    }
    return FALSE;
}

/* 直線とselfとの交点を求める.
 * 入力
 *     self : UpcShapeオブジェクト
 *     pos : 直線が通る点の座標
 *     dir : 直線の方向余弦
 *     area : UpcAreaオブジェクト
 * 出力
 *     p : 交点座標
 *     n : その点での単位法線ベクトル
 *     isActive : 交点が有効領域にあるか否かを示すBool値
 *     戻り値 : 交点が求まればTRUE, 求まらなければFALSE
 */
static Bool getIntersection(const UpcShape *self, const double *pos, const double *dir, UpcArea *area, double *p, double *n, Bool *isActive)
{
    double p2[3];
    double curv = self->curv;
    int nPoint = UpcShapeSph_conicIntersection(curv, 0.0, pos, dir, p, p2);
    double zx, zy, ilen, dmy;

    if (nPoint == 2) {
        Bool a1 = UpcArea_isActive(area, p[X], p[Y]);
        Bool a2 = UpcArea_isActive(area, p2[X], p2[Y]);

        if (a1) {
            *isActive = TRUE;
        }
        else {
            if (a2) {
                UpcVector_COPY(p, p2);
                *isActive = TRUE;
            }
            else {
                *isActive = FALSE;
            }
        }
        goto pReturn;
    }
    else if (nPoint) {
        *isActive = UpcArea_isActive(area, p[X], p[Y]);
        goto pReturn;
    }
    return FALSE;

pReturn:
    UpcShapeSph_conicZ(curv, 0.0, p[X], p[Y], 1, &dmy, &zx, &zy, &dmy, &dmy, &dmy);
    ilen = 1.0 / sqrt(zx * zx + zy * zy + 1.0);
    UpcVector_SET(n, -ilen * zx, -ilen * zy, ilen);
    return TRUE;    
}

static Bool xSymmetric(const UpcShape *self)
{
    return TRUE;
}

static Bool ySymmetric(const UpcShape *self)
{
    return TRUE;
}

static Bool rSymmetric(const UpcShape *self)
{
    return TRUE;
}

/* 円錐曲線のz座標と微分値を求める.
 * 入力
 *     curv : 曲率
 *     koni : コーニック係数
 *     x : x座標
 *     y : y座標
 *     mode : = 0 : z を求める.
 *            = 1 : 1次微分まで(z, zx, zy)を求める.
 *            = 2 : 2次微分まで(z, zx, zy, zxx, zxy, zyy)を求める.
 * 出力
 *     z   : z
 *     zx  : ∂z/∂x
 *     zy  : ∂z/∂y
 *     zxx : ∂2z/∂x2
 *     zxy : ∂2z/∂x∂y
 *     zyy : ∂2z/∂y2
 *     戻り値 : 正常に計算できればTRUE, 定義域外ならFALSE
 */
Bool UpcShapeSph_conicZ(double curv, double koni, double x, double y, int mode, double *z, double *zx, double *zy, double *zxx, double *zxy, double *zyy)
{
    double r2 = x * x + y * y;
    double a = (1.0 + koni) * curv * curv;
    double b, b2 = 1.0 - a * r2;
    
    if (b2 <= 0.0) {
        // conic surface domain error
        *z = curv * r2; // FALSEを返すがzには境界値(b2 == 0 の値)を入れておく. (描画のため)
        return FALSE;
    }
    b = sqrt(b2);
    *z = curv * r2 / (1.0 + b);
    if (mode >= 1) {
        double zr_r = curv / b;

        *zx = zr_r * x;
        *zy = zr_r * y;

        if (mode == 2) {
            double e = curv / (b * b2);

            *zxx = e * (1.0 - a * y * y);
            *zxy = e * a * x * y;
            *zyy = e * (1.0 - a * x * x);
        }
    }
    return TRUE;
}

/* 2次方程式 a x^2 + b x + c = 0 の実数解を求める.
 * 入力
 *     方程式の係数 a, b, c
 * 出力
 *     戻り値 : 解の個数n (n = 0, 1, 2)
 *     k1 : 解1 (n >= 1 のとき有効な値が入る)
 *     k2 : 解2 (n == 2 のとき有効な値が入る)
 */
int UpcShapeSph_solveQuadEq(double a, double b, double c, double *k1, double *k2)
{
    if (a) {
        b /= a;
        c /= a;
        if (c) {
            double d;
            b *= 0.5;
            d = b * b - c;
            if (d == 0.0) {      // 重解
                *k1 = -b;
                *k2 = -b;
                return 2;
            }
            else if (d > 0.0) {  // x^2 + 2 b x + c = 0
                double alpha, beta;
                if (b >= 0.0)
                    alpha = - b - sqrt(d);
                else
                    alpha = - b + sqrt(d);
                beta = c / alpha;
                *k1 = beta;
                *k2 = alpha;
                return 2;
            }
            else {               // 虚数解
                return 0;
            }
        }
        else {                   // x^2 + b x = 0
            *k1 = -b;
            *k2 = 0.0;
            return 2;
        }
    }
    else {
        if (b) {                 // 1次方程式
            *k1 = -c / b;
            return 1;
        }
        else {                   // 不定または解なし
            return 0;
        }
    }
}

/* 円錐曲面と直線の交点座標を求める.
 * 入力
 *     curv : 曲率
 *     koni : コーニック係数
 *     p[3] : 直線が通る点
 *     d[3] : 直線の方向余弦
 * 出力
 *     戻り値 : 交点の個数n (n = 0, 1, 2)
 *     p1[3]  : 交点1の座標 (n >= 1 のとき有効な値が入る)
 *     p2[3]  : 交点2の座標 (n == 2 のとき有効な値が入る)
 *     |p1[Z]| ≦ |p2[Z]| となる.
 */
int UpcShapeSph_conicIntersection(double curv, double koni, const double *p, const double *d, double *p1, double *p2)
{
    double px = p[X];
    double py = p[Y];
    double pz = p[Z];
    double dx = d[X];
    double dy = d[Y];
    double dz = d[Z];
    
    if (curv == 0.0) {
        if (dz == 0.0){
            return 0;
        }
        else {
            double t = - pz / dz;

            UpcVector_SET(p1, px + dx * t, py + dy * t, 0.0);
            return 1;
        }
    }
    else {
        double kp1 = 1.0 + koni;
        double a, b, c, k1, k2;
        int nk;

        a = dx * dx + dy * dy + kp1 * dz * dz;
        b = dx * px + dy * py + kp1 * dz * pz - dz / curv;
        c = px * px + py * py + kp1 * pz * pz - 2.0 * pz / curv;
        nk = UpcShapeSph_solveQuadEq(a, 2.0 * b, c, &k1, &k2);
        if (nk == 2) {
            double z1 = pz + dz * k1;
            double z2 = pz + dz * k2;

            if (kp1 * curv * z1 <= 1.0) { // 元の円錐曲線に乗っているかの確認
                if (kp1 * curv * z2 <= 1.0) {
                    if (fabs(z1) < fabs(z2)) {
                        UpcVector_SET(p1, px + dx * k1, py + dy * k1, z1);
                        UpcVector_SET(p2, px + dx * k2, py + dy * k2, z2);
                    }
                    else {
                        UpcVector_SET(p1, px + dx * k2, py + dy * k2, z2);
                        UpcVector_SET(p2, px + dx * k1, py + dy * k1, z1);
                    }
                    return 2;
                }
                else {
                    UpcVector_SET(p1, px + dx * k1, py + dy * k1, z1);
                    return 1;
                }
            }
            else {
                if (kp1 * curv * z2 <= 1.0) {
                    UpcVector_SET(p1, px + dx * k2, py + dy * k2, z2);
                    return 1;
                }
                else
                    return 0;
            }
        }
        else if (nk == 1) {
            double z = pz + dz * k1;

            if (kp1 * curv * z <= 1.0) {
                UpcVector_SET(p1, px + dx * k1, py + dy * k1, z);
                return 1;
            }
            else {
                return 0;
            }
        }
        return 0;
    }
}


