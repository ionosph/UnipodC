/*
 *  UpcShapeAsph.c
 *  UnipodC
 *
 *  Created by ionosph on 2010/10/20.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcShapeAsph.h"
#include "UpcVector.h"
#include "UpcArea.h"

static const double curvDefault = 0.0;
static const double coniDefault = 0.0;
static const double coefDefault = 0.0;
static const double nrmR2invDefault = 1.0;
static const int itermaxDefault = 100;
static const double tolDefault = 1.0e-12;

static double axialPower(const UpcShape *self);
static Bool getZ(const UpcShape *self, double x, double y, double *z);
static Bool getCurv(const UpcShape *self, double y, double *curvMeri, double *curvSagi);
static void getRBC(const UpcShape *self, double *rv, double *bv, double *cv);
static Bool getNormalVector(const UpcShape *self, double x, double y, double *nv);
static Bool getIntersection(const UpcShape *self, const double *pos, const double *dir, UpcArea *area, double *p, double *n, Bool *isActive);
static Bool xSymmetric(const UpcShape *self);
static Bool ySymmetric(const UpcShape *self);
static Bool rSymmetric(const UpcShape *self);

static void UpcShapeAsph_dealloc(UpcShape *self)
{
    free(self->rpc);
    free(self);
}

static void UpcShapeAsph_print(const UpcShape *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);
    int i;

    fprintf(ostream, "%stype    = %s\n", space, UpcShapeType_string(self->shapeType));
    fprintf(ostream, "%sr       = %13.6e\n", space, UpcShapeSph_RADI(self));
    fprintf(ostream, "%sconi    = %13.6e\n", space, self->coni);
    for (i = 2; i <= self->romax; i++) {
        if (self->rpc[i])
            fprintf(ostream, "%sac[%2d]  = %13.6e\n", space, i, self->rpc[i]);
    }
    fprintf(ostream, "%snrmR    = %13.6e\n", space, UpcShapeAsph_NRMR(self));
    fprintf(ostream, "%sitermax = %d\n", space, self->itermax);
    fprintf(ostream, "%stol     = %13.6e\n", space, self->tol);
}

static void UpcShapeAsph_toXml(const UpcShape *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);
    int i;

    fprintf(ostream, "%s<shape type=\"%s\" r_order_max=\"%d\">\n", space, UpcShapeType_string(self->shapeType), self->romax);
    if (self->curv != curvDefault)
        fprintf(ostream, "%s  <r>%23.15e</r>\n", space, UpcShapeSph_RADI(self));
    if (self->coni != coniDefault)
        fprintf(ostream, "%s  <coni>%23.15e</coni>\n", space, self->coni);
    for (i = 2; i <= self->romax; i++) {
        if (self->rpc[i] != coefDefault)
            fprintf(ostream, "%s  <rcoef order=\"%d\">%23.15e</rcoef>\n", space, i, self->rpc[i]);
    }
    if (self->nrmR2inv != nrmR2invDefault)
        fprintf(ostream, "%s  <nrm_r>%23.15e</nrm_r>\n", space, UpcShapeAsph_NRMR(self));
    if (self->itermax != itermaxDefault)
        fprintf(ostream, "%s  <iter_max>%d</iter_max>\n", space, self->itermax);
    if (self->tol != tolDefault)
        fprintf(ostream, "%s  <tol>%23.15e</tol>\n", space, self->tol);
    fprintf(ostream, "%s</shape>\n", space);
}

static UpcErrStatus UpcShapeAsph_lset(UpcShape *self)
{
    return UpcE_NoIssues;
}

static UpcShape *UpcShapeAsph_copy(const UpcShape *self)
{
    UpcShape *t = UpcShapeAsph_init(self->romax);
    int i;
    
    if (t) {
        t->curv = self->curv;
        t->coni = self->coni;
        // t->romax = self->romax;
        t->nrmR2inv = self->nrmR2inv;
        t->itermax = self->itermax;
        t->tol = self->tol;
        for (i = 2; i <= t->romax; i++)
            t->rpc[i] = self->rpc[i];
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeAsph_copy");
    }
    return t;
}

UpcShape *UpcShapeAsph_init(int max_order)
{
    UpcShape *self;
    int i;

    if (max_order < 2) {
        UpcERRHANDLER(UpcE_ValueError, "invalid max_order (in UpcShapeAsph_init)");
        return NULL;
    }

    self = (UpcShape *)calloc(1, sizeof(UpcShape));
    if (self) {
        self->rpc = (double *)calloc(max_order + 1, sizeof(double));
        if (!(self->rpc)) {
            free(self);
            UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeAsph_init");
            return NULL;
        }
        self->_refCount = 1;
        self->_dealloc = UpcShapeAsph_dealloc;
        self->shapeType = UpcShapeType_asph;
        self->curv = curvDefault;
        self->coni = coniDefault;
        self->romax = max_order;
        if (coefDefault) {
            for (i = 2; i <= max_order; i++)
                self->rpc[i] = coefDefault;
        }
        self->nrmR2inv = nrmR2invDefault;
        self->itermax = itermaxDefault;
        self->tol = tolDefault;
        // method
        self->axialPower = axialPower;
        self->getZ = getZ;
        self->getCurv = getCurv;
        self->getRBC = getRBC;
        self->getNormalVector = getNormalVector;
        self->getIntersection = getIntersection;
        self->xSymmetric = xSymmetric;
        self->ySymmetric = ySymmetric;
        self->rSymmetric = rSymmetric;
        self->copy = UpcShapeAsph_copy;
        self->print = UpcShapeAsph_print;
        self->toXml = UpcShapeAsph_toXml;
        self->lset = UpcShapeAsph_lset;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeAsph_init");
    }
    return self;
}

UpcShape *UpcShapeAsph_initWithShape(int max_order, UpcShape *other)
{
    UpcShape *self = UpcShapeAsph_init(max_order);
    enum UpcShapeType otype = other->shapeType;
    
    if (self) {
        if (otype % UpcShapeType_sph == 0) { // UpcShapeSphクラスまたは子クラスであれば
            self->curv = other->curv;
        }
        if (otype % UpcShapeType_asph == 0) { // UpcShapeAsphクラスまたは子クラスであれば
            int i, n = (max_order < other->romax) ? max_order : other->romax;
            
            self->coni = other->coni;
            for (i = 2; i <= n; i++) {
                self->rpc[i] = other->rpc[i];
            }
            self->nrmR2inv = other->nrmR2inv;
            self->itermax = other->itermax;
            self->tol = other->tol;
        }
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeAsph_initWithShape");
    }
    return self;
}

double UpcShapeAsph_nrmR(const UpcShape *self)
{
    return 1.0 / sqrt(self->nrmR2inv);
}

Bool UpcShapeAsph_setNrmR(UpcShape *self, double nrmR)
{
    if (nrmR) {
        self->nrmR2inv = 1.0 / (nrmR * nrmR);
        return TRUE;
    }
    return FALSE;
}

static double axialPower(const UpcShape *self)
{
    return self->curv + 2.0 * self->rpc[2] * self->nrmR2inv;
}

static Bool getZ(const UpcShape *self, double x, double y, double *z)
{
    double zconi, zpoly, dmy;

    if (UpcShapeSph_conicZ(self->curv, self->coni, x, y, 0, &zconi, &dmy, &dmy, &dmy, &dmy, &dmy)) {
        UpcShapeAsph_rpolyZ(self->romax, self->rpc, x, y, self->nrmR2inv, 0, &zpoly, &dmy, &dmy, &dmy, &dmy, &dmy);
        *z = zconi + zpoly;
        return TRUE;
    }
    *z = zconi; // FALSEを返すが, zには値を入れておく. (描画のため)
    return FALSE;
}

static Bool getCurv(const UpcShape *self, double y, double *curvMeri, double *curvSagi)
{
    if (y) {
        double zyconi, zypoly, zyyconi, zyypoly, zy, zyy, dmy;
        
        y = fabs(y);
        if (UpcShapeSph_conicZ(self->curv, self->coni, 0.0, y, 2, &dmy, &dmy, &zyconi, &dmy, &dmy, &zyyconi)) {
            UpcShapeAsph_rpolyZ(self->romax, self->rpc, 0.0, y, self->nrmR2inv, 2, &dmy, &dmy, &zypoly, &dmy, &dmy, &zyypoly);
            zy = zyconi + zypoly;
            zyy = zyyconi + zyypoly;
            *curvMeri = zyy / pow(1.0 + zy * zy, 1.5);
            *curvSagi = (zy / y) / sqrt(1.0 + zy * zy);
            return TRUE;
        }
        else {
            // no intersection
            return FALSE;
        }
    }
    else {
        double cuy = self->curv + 2.0 * self->rpc[2] * self->nrmR2inv;
        
        *curvMeri = cuy;
        *curvSagi = cuy;
        return TRUE;
    }
}

static void getRBC(const UpcShape *self, double *rv, double *bv, double *cv)
{
    double k = self->coni;
    double kp1 = k + 1.0;
    double Av = self->rpc[2] * self->nrmR2inv;
    double Av2 = Av * Av;
    double Av4 = Av2 * Av2;
    double c = self->curv + 2.0 * Av;
    double cvcv = c * self->curv;
    double ac4 = (self->romax >= 4) ? (self->rpc[4] * self->nrmR2inv * self->nrmR2inv) : 0.0;
    double ac6 = (self->romax >= 6) ? (self->rpc[6] * self->nrmR2inv * self->nrmR2inv * self->nrmR2inv) : 0.0;

    *rv = c;
    *bv = k * pow(c, 3) + 8.0 * ac4 - 2.0 * Av * kp1 * (4.0 * Av2 + 3.0 * cvcv);
    *cv = k * (k + 2.0) * pow(c, 5)
        + 16.0 * ac6
        - 2.0 * Av * kp1 * kp1 * (16.0 * Av4 + 20.0 * Av2 * cvcv + 5.0 * cvcv * cvcv);
}

static Bool getNormalVector(const UpcShape *self, double x, double y, double *nv)
{
    double zxconi, zyconi, dmy;
    double zxpoly, zypoly, tx, ty, ilen;

    if (UpcShapeSph_conicZ(self->curv, self->coni, x, y, 1, &dmy, &zxconi, &zyconi, &dmy, &dmy, &dmy)) {
        UpcShapeAsph_rpolyZ(self->romax, self->rpc, x, y, self->nrmR2inv, 1, &dmy, &zxpoly, &zypoly, &dmy, &dmy, &dmy);
        tx = -(zxconi + zxpoly);
        ty = -(zyconi + zypoly);
        ilen = 1.0 / sqrt(tx * tx + ty * ty + 1.0);
        UpcVector_SET(nv, ilen * tx, ilen * ty, ilen);
        return TRUE;
    }
    return FALSE;
}

/* 逐次近似によって非球面の交点を求めるルーチン(内部関数)
 * 入力
 *     pos : 交点座標初期値(ベースの2次曲面との交点座標を与える)
 *     dir : 方向余弦
 * 出力
 *     p[3] : 交点座標
 *     n[3] : 交点における単位法線ベクトル
 */
static Bool iterAsph(const UpcShape *self, const double *pos, const double *dir, double *p, double *n)
{
    const double curv = self->curv;
    const double coni = self->coni;
    const int max_order = self->romax;
    const double *ac = self->rpc;
    const double nrmR2inv = self->nrmR2inv;
    const double iter_max = self->itermax;
    const double tol = self->tol;
    int i;
    double asphz, dpz;
    double zconi, zxconi, zyconi, zpoly, zxpoly, zypoly, dmy;
    double nx, ny, stp, ilen;
    double px = pos[X];
    double py = pos[Y];
    double pz = pos[Z];
    double dx = dir[X];
    double dy = dir[Y];
    double dz = dir[Z];
    
    for (i = 0; i < iter_max; i++) {
        if (UpcShapeSph_conicZ(curv, coni, px, py, 1, &zconi, &zxconi, &zyconi, &dmy, &dmy, &dmy)) {
            UpcShapeAsph_rpolyZ(max_order, ac, px, py, nrmR2inv, 1, &zpoly, &zxpoly, &zypoly, &dmy, &dmy, &dmy);
            asphz = zconi + zpoly;
            nx = -(zxconi + zxpoly);
            ny = -(zyconi + zypoly);
            // nz = 1.0; (正規化されていない法線ベクトル.)
        }
        else
            return FALSE;
        // stp = (asphz - pz) * nz / (nx * dx + ny * dy + nz * dz);
        stp = (asphz - pz) / (nx * dx + ny * dy + dz); // nz == 1 の場合
        dpz = dz * stp;
        if (fabs(dpz) < tol) {
            UpcVector_SET(p, px, py, pz);
            ilen = 1.0 / sqrt(nx * nx + ny * ny + 1.0);
            UpcVector_SET(n, nx * ilen, ny * ilen, ilen);
            return TRUE;
        }
        px += dx * stp;
        py += dy * stp;
        pz += dpz;
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
    double p1[3], p2[3];
    int np = UpcShapeSph_conicIntersection(self->curv, self->coni, pos, dir, p1, p2);

    if (np == 2) {
        double q1[3], q2[3], q3[3], n1[3], n2[3];
        Bool flg1 = iterAsph(self, p1, dir, q1, n1);
        Bool flg2 = iterAsph(self, p2, dir, q2, n2);
        Bool a1 = UpcArea_isActive(area, q1[X], q1[Y]);
        Bool a2 = UpcArea_isActive(area, q2[X], q2[Y]);
        Bool a3 = a1;

        if (flg1 && flg2) {
            if (fabs(q1[Z]) > fabs(q2[Z])) {
                UpcVector_COPY(q3, q1);
                UpcVector_COPY(q1, q2);
                UpcVector_COPY(q2, q3);
                a1 = a2;
                a2 = a3;
            }
            if (a1) {
                UpcVector_COPY(p, q1);
                UpcVector_COPY(n, n1);
                *isActive = TRUE;
            }
            else {
                if (a2) {
                    UpcVector_COPY(p, q2);
                    UpcVector_COPY(n, n2);
                    *isActive = TRUE;
                }
                else {
                    UpcVector_COPY(p, q1);
                    UpcVector_COPY(n, n1);
                    *isActive = FALSE;
                }
            }
            return TRUE;
        }
        else {
            if (flg1) {
                UpcVector_COPY(p, q1);
                UpcVector_COPY(n, n1);
                *isActive = a1;
                return TRUE;
            }
            if (flg2) {
                UpcVector_COPY(p, q2);
                UpcVector_COPY(n, n2);
                *isActive = a2;
                return TRUE;
            }
        }
        return FALSE;
    }
    else if (np) {
        Bool flg = iterAsph(self, p1, dir, p, n);

        *isActive = UpcArea_isActive(area, p[X], p[Y]);
        return flg;
    }
    return FALSE;
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

/* r多項式のz座標と微分を求める. (z = a2 r^2 + a3 r^3 + ... + aN r^N (N = max_order))
 * 入力
 *     max_order : R多項式の最大次数
 *     ac = ac[0..max_order] : 多項式の係数(ただし, ac[0], ac[1]は使われない)
 *     x : x座標
 *     y : y座標
 *     nrmR2inv : 正規化半径の自乗の逆数
 *     mode : = 0 : z を求める.
 *            = 1 : 1次微分まで(z, zx, zy)を求める.
 *            = 2 : 2次微分まで(z, zx, zy, zxx, zxy, zyy)を求める.
 * 出力
 *     z   : z
 *     zx  : ∂z/∂x
 *     zy  : ∂z/∂y
 *     zxx : ∂2z/∂x2
 *     zxy : ∂2z/∂x∂y
 *     zyy : ∂2z/∂y
 *     戻り値 : 常にTRUE
 */
Bool UpcShapeAsph_rpolyZ(int max_order, const double *ac, double x, double y, double nrmR2inv, int mode, double *z, double *zx, double *zy, double *zxx, double *zxy, double *zyy)
{
    double rdsh2 = (x * x + y * y) * nrmR2inv;
    double rdsh = sqrt(rdsh2);
    int i;
    
    if (!rdsh2) {
        *z   = 0.0;
        *zx  = 0.0;
        *zy  = 0.0;
        *zxx = 2.0 * ac[2] * nrmR2inv;
        *zxy = 0.0;
        *zyy = 2.0 * ac[2] * nrmR2inv;
        return TRUE;
    }

    if (mode == 0) {
        double f = ac[max_order];

        for (i = max_order - 1; i >= 2; i--)
            f = f * rdsh + ac[i];
        *z = f * rdsh2;
    }
    else if (mode == 1) {
        double f = ac[max_order];
        double fr_r = max_order * ac[max_order];

        for (i = max_order - 1; i >= 2; i--) {
            f = f * rdsh + ac[i];
            fr_r = fr_r * rdsh + i * ac[i];
        }
        *z = f * rdsh2;
        *zx = fr_r * x * nrmR2inv;
        *zy = fr_r * y * nrmR2inv;
    }
    else if (mode == 2) {
        double x2_r2 = x * x / rdsh2 * nrmR2inv;
        double xy_r2 = x * y / rdsh2 * nrmR2inv;
        double y2_r2 = y * y / rdsh2 * nrmR2inv;
        double f = ac[max_order];
        double fr_r = max_order * ac[max_order];
        double frr = max_order * (max_order - 1) * ac[max_order];

        for (i = max_order - 1; i >= 2; i--) {
            f = f * rdsh + ac[i];
            fr_r = fr_r * rdsh + i * ac[i];
            frr = frr * rdsh + i * (i - 1) * ac[i];
        }
        *z = f * rdsh2;
        *zx = fr_r * x * nrmR2inv;
        *zy = fr_r * y * nrmR2inv;
        *zxx = nrmR2inv * (x2_r2 * frr + y2_r2 * fr_r);
        *zxy = nrmR2inv * xy_r2 * (frr - fr_r);
        *zyy = nrmR2inv * (y2_r2 * frr + x2_r2 * fr_r);
    }
    return TRUE;
}

