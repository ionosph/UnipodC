/*
 *  UpcShapeXYP.c
 *  UnipodC
 *
 *  Created by ionosph on 2011/04/28.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcShapeXYP.h"
#include "UpcVector.h"
#include "UpcArea.h"
#include "UpcZernikeToXYP.h"

static const double curvDefault = 0.0;
static const double coniDefault = 0.0;
static const double coefDefault = 0.0;
static const double nrmR2invDefault = 1.0;
static const int itermaxDefault = 100;
static const double tolDefault = 1.0e-12;
static const double xyCoefDefault = 0.0;
static const double nrmXinvDefault = 1.0;
static const double nrmYinvDefault = 1.0;
static const double nrmZRinvDefault = 1.0;

static double axialPower(const UpcShape *self);
static Bool getZ(const UpcShape *self, double x, double y, double *z);
static Bool getCurv(const UpcShape *self, double y, double *curvMeri, double *curvSagi);
static void getRBC(const UpcShape *self, double *rv, double *bv, double *cv);
static Bool getNormalVector(const UpcShape *self, double x, double y, double *nv);
static Bool getIntersection(const UpcShape *self, const double *pos, const double *dir, UpcArea *area, double *p, double *n, Bool *isActive);
static Bool xSymmetric(const UpcShape *self);
static Bool ySymmetric(const UpcShape *self);
static Bool rSymmetric(const UpcShape *self);

static void UpcShapeXYP_dealloc(UpcShape *self)
{
    upc_freeD2(self->xypc);
    free(self->rpc);
    free(self);
}

static void UpcShapeXYP_print(const UpcShape *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);
    int i, ix, iy;

    fprintf(ostream, "%stype    = %s\n", space, UpcShapeType_string(self->shapeType));
    fprintf(ostream, "%sr    = %13.6e\n", space, UpcShapeSph_RADI(self));
    fprintf(ostream, "%sconi = %13.6e\n", space, self->coni);
    for (i = 2; i <= self->romax; i++) {
        if (self->rpc[i])
            fprintf(ostream, "%srpc[%2d] = %13.6e\n", space, i, self->rpc[i]);
    }
    fprintf(ostream, "%snrmR    = %13.6e\n", space, UpcShapeAsph_nrmR(self));
    fprintf(ostream, "%sitermax = %d\n", space, self->itermax);
    fprintf(ostream, "%stol     = %13.6e\n", space, self->tol);
    for (iy = 0; iy <= self->yomax; iy++) {
        for (ix = 0; ix <= self->xomax; ix++) {
            if (self->xypc[iy][ix])
                fprintf(ostream, "%sxypc[%2d][%2d] = %13.6e\n", space, iy, ix, self->xypc[iy][ix]);
        }
    }
    fprintf(ostream, "%snrmX    = %13.6e\n", space, UpcShapeXYP_nrmX(self));
    fprintf(ostream, "%snrmY    = %13.6e\n", space, UpcShapeXYP_nrmY(self));
}


static void UpcShapeXYP_toXml(const UpcShape *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);
    int i, ix, iy;

    fprintf(ostream, "%s<shape type=\"%s\" r_order_max=\"%d\" x_order_max=\"%d\" y_order_max=\"%d\">\n", space, UpcShapeType_string(self->shapeType), self->romax, self->xomax, self->yomax);
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
    for (iy = 0; iy <= self->yomax; iy++) {
        for (ix = 0; ix <= self->xomax; ix++) {
            if (self->xypc[iy][ix] != xyCoefDefault)
                fprintf(ostream, "%s  <xycoef xorder=\"%d\" yorder=\"%d\">%23.15e</xycoef>\n", space, ix, iy, self->xypc[iy][ix]);
        }
    }
    if (self->nrmXinv != nrmXinvDefault)
        fprintf(ostream, "%s  <nrm_x>%23.15e</nrm_x>\n", space, UpcShapeXYP_nrmX(self));
    if (self->nrmYinv != nrmYinvDefault)
        fprintf(ostream, "%s  <nrm_y>%23.15e</nrm_y>\n", space, UpcShapeXYP_nrmY(self));
    fprintf(ostream, "%s</shape>\n", space);
}

static UpcErrStatus UpcShapeXYP_lset(UpcShape *self)
{
    return UpcE_NoIssues;
}

static UpcShape *UpcShapeXYP_copy(const UpcShape *self)
{
    UpcShape *t = UpcShapeXYP_init(self->romax, self->xomax, self->yomax);
    int i, iy, ix;
    
    if (t) {
        t->curv = self->curv;
        t->coni = self->coni;
        // t->romax = self->romax;
        t->nrmR2inv = self->nrmR2inv;
        t->itermax = self->itermax;
        t->tol = self->tol;
        for (i = 2; i <= t->romax; i++)
            t->rpc[i] = self->rpc[i];
        t->xomax = self->xomax;
        t->yomax = self->yomax;
        t->nrmXinv = self->nrmXinv;
        t->nrmYinv = self->nrmYinv;
        for (iy = 0; iy <= t->yomax; iy++) {
            for (ix = 0; ix <= t->xomax; ix++) {
                t->xypc[iy][ix] = self->xypc[iy][ix];
            }
        }
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeXYP_copy");
    }
    return t;
}

UpcShape *UpcShapeXYP_init(int max_orderR, int max_orderX, int max_orderY)
{
    UpcShape *self;
    
    if (max_orderR < 2 || max_orderX < 0 || max_orderY < 0) {
        UpcERRHANDLER(UpcE_ValueError, "invalid parameter (in UpcShapeXYP_init)");
        return NULL;
    }
    
    self = (UpcShape *)calloc(1, sizeof(UpcShape));
    if (self) {
        self->rpc = (double *)calloc(max_orderR + 1, sizeof(double));
        if (!(self->rpc)) {
            free(self);
            UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeXYP_init");
            return NULL;
        }
        self->xypc = upc_callocD2(max_orderY + 1, max_orderX + 1);
        if (!(self->xypc)) {
            free(self->rpc);
            free(self);
            UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeXYP_init");
            return NULL;
        }
        self->_refCount = 1;
        self->_dealloc = UpcShapeXYP_dealloc;
        self->shapeType = UpcShapeType_xyp;
        self->curv = curvDefault;
        self->coni = coniDefault;
        self->romax = max_orderR;
        if (coefDefault) {
            int i;

            for (i = 2; i <= max_orderR; i++)
                self->rpc[i] = coefDefault;
        }
        self->nrmR2inv = nrmR2invDefault;
        self->itermax = itermaxDefault;
        self->tol = tolDefault;
        self->xomax = max_orderX;
        self->yomax = max_orderY;
        if (xyCoefDefault) {
            int ix, iy;

            for (iy = 0; iy <= max_orderY; iy++) {
                for (ix = 0; ix <= max_orderX; ix++) {
                    self->xypc[iy][ix] = xyCoefDefault;
                }
            }
        }
        self->nrmXinv = nrmXinvDefault;
        self->nrmYinv = nrmYinvDefault;
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
        self->copy = UpcShapeXYP_copy;
        self->print = UpcShapeXYP_print;
        self->toXml = UpcShapeXYP_toXml;
        self->lset = UpcShapeXYP_lset;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeXYP_init");
    }
    return self;
}

UpcShape *UpcShapeXYP_initWithShape(int max_orderR, int max_orderX, int max_orderY, UpcShape *other)
{
    UpcShape *self = UpcShapeXYP_init(max_orderR, max_orderX, max_orderY);
    enum UpcShapeType otype = other->shapeType;
    
    if (self) {
        if (otype % UpcShapeType_sph == 0) { // UpcShapeSphクラスまたは子クラスであれば
            self->curv = other->curv;
        }
        if (otype % UpcShapeType_asph == 0) { // UpcShapeAsphクラスまたは子クラスであれば
            int i, n = (max_orderR < other->romax) ? max_orderR : other->romax;
            
            self->coni = other->coni;
            for (i = 2; i <= n; i++) {
                self->rpc[i] = other->rpc[i];
            }
            self->nrmR2inv = other->nrmR2inv;
            self->itermax = other->itermax;
            self->tol = other->tol;
        }
        if (otype % UpcShapeType_xyp == 0) { // UpcShapeXYPクラスまたは子クラスであれば
            int ix, iy;
            int xomax = (max_orderX < other->xomax) ? max_orderX : other->xomax;
            int yomax = (max_orderY < other->yomax) ? max_orderY : other->yomax;

            for (iy = 0; iy <= yomax; iy++) {
                for (ix = 0; ix <= xomax; ix++) {
                    self->xypc[iy][ix] = other->xypc[iy][ix];
                }
            }
            self->nrmXinv = other->nrmXinv;
            self->nrmYinv = other->nrmYinv;
        }
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeXYP_initWithShape");
    }
    return self;
}

double UpcShapeXYP_nrmX(const UpcShape *self)
{
    return 1.0 / self->nrmXinv;
}

Bool UpcShapeXYP_setNrmX(UpcShape *self, double nrmX)
{
    if (nrmX) {
        self->nrmXinv = 1.0 / nrmX;
        return TRUE;
    }
    return FALSE;
}

double UpcShapeXYP_nrmY(const UpcShape *self)
{
    return 1.0 / self->nrmYinv;
}

Bool UpcShapeXYP_setNrmY(UpcShape *self, double nrmY)
{
    if (nrmY) {
        self->nrmYinv = 1.0 / nrmY;
        return TRUE;
    }
    return FALSE;
}

void UpcShapeXYP_clearXYPC(UpcShape *self)
{
    const int xomax = self->xomax;
    const int yomax = self->yomax;
    int ix, iy;

    for (iy = 0; iy <= yomax; iy++) {
        for (ix = 0; ix <= xomax; ix++) {
            self->xypc[iy][ix] = 0.0;
        }
    }
}

Bool UpcShapeXYP_addZcoef(UpcShape *self, int cno, double coef)
{
    UpcQMatrix *zMatrix;
    int ix, iy, xIsOdd, yIsOdd, ixMax, iyMax;
    double z;

    if (cno <= 0 || cno > 400) {
        UpcERRHANDLER(UpcE_ValueError, "cno is out of range");
        return FALSE;
    }

    zMatrix = UpcZernikeToXYP_zcoefxy(cno);
    ixMax = zMatrix->ixMax;
    iyMax = zMatrix->iyMax;

    if (ixMax > self->xomax) {
        UpcERRHANDLER(UpcE_ValueError, "There is not sufficient size for X.");
        return FALSE;
    }
    if (iyMax > self->yomax) {
        UpcERRHANDLER(UpcE_ValueError, "There is not sufficient size for Y.");
        return FALSE;
    }

    xIsOdd = zMatrix->xIsOdd;
    yIsOdd = zMatrix->yIsOdd;
    for (ix = xIsOdd; ix <= ixMax; ix += 2) {
        for (iy = yIsOdd; iy <= iyMax; iy += 2) {
            z = UpcQMatrix_data(zMatrix, iy, ix);
            if (z) {
                self->xypc[iy][ix] += z * coef;
            }
        }
    }
    return TRUE;
}

static double axialPower(const UpcShape *self)
{
    return self->curv + 2.0 * self->rpc[2] * self->nrmR2inv; // XYベキ非球面は近軸計算に反映されない.
}

static Bool getZ(const UpcShape *self, double x, double y, double *z)
{
    double zconi, zrp, zxyp, dmy;

    if (UpcShapeSph_conicZ(self->curv, self->coni, x, y, 0, &zconi, &dmy, &dmy, &dmy, &dmy, &dmy)) {
        UpcShapeAsph_rpolyZ(self->romax, self->rpc, x, y, self->nrmR2inv, 0, &zrp, &dmy, &dmy, &dmy, &dmy, &dmy);
        UpcShapeXYP_xypolyZ(self->xomax, self->yomax, self->xypc, x, y, self->nrmXinv, self->nrmYinv, 0, &zxyp, &dmy, &dmy, &dmy, &dmy, &dmy);
        *z = zconi + zrp + zxyp;
        return TRUE;
    }
    *z = zconi; // FALSEを返すが, zには値を入れておく. (描画のため)
    return FALSE;
}

static Bool getCurv(const UpcShape *self, double y, double *curvMeri, double *curvSagi)
{
    if (y) {
        double zyconi, zyrp, zyxyp, zyyconi, zyyrp, zyyxyp, zy, zyy, dmy;

        y = fabs(y);
        if (UpcShapeSph_conicZ(self->curv, self->coni, 0.0, y, 2, &dmy, &dmy, &zyconi, &dmy, &dmy, &zyyconi)) {
            UpcShapeAsph_rpolyZ(self->romax, self->rpc, 0.0, y, self->nrmR2inv, 2, &dmy, &dmy, &zyrp, &dmy, &dmy, &zyyrp);
            UpcShapeXYP_xypolyZ(self->xomax, self->yomax, self->xypc, 0.0, y, self->nrmXinv, self->nrmYinv, 2, &dmy, &dmy, &zyxyp, &dmy, &dmy, &zyyxyp);
            zy = zyconi + zyrp + zyxyp;
            zyy = zyyconi + zyyrp + zyyxyp;
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

        *curvMeri = cuy + 2.0 * self->xypc[2][0] * self->nrmYinv * self->nrmYinv;
        *curvSagi = cuy + 2.0 * self->xypc[0][2] * self->nrmXinv * self->nrmXinv;
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
    // XYベキ非球面は3次収差計算に反映されない.
}

static Bool getNormalVector(const UpcShape *self, double x, double y, double *nv)
{
    double zxconi, zyconi, dmy;
    double zx_rp, zy_rp, zx_xyp, zy_xyp, tx, ty, ilen;

    if (UpcShapeSph_conicZ(self->curv, self->coni, x, y, 1, &dmy, &zxconi, &zyconi, &dmy, &dmy, &dmy)) {
        UpcShapeAsph_rpolyZ(self->romax, self->rpc, x, y, self->nrmR2inv, 1, &dmy, &zx_rp, &zy_rp, &dmy, &dmy, &dmy);
        UpcShapeXYP_xypolyZ(self->xomax, self->yomax, self->xypc, x, y, self->nrmXinv, self->nrmYinv, 1, &dmy, &zx_xyp, &zy_xyp, &dmy, &dmy, &dmy);
        tx = -(zxconi + zx_rp + zx_xyp);
        ty = -(zyconi + zy_rp + zy_xyp);
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
    const int romax = self->romax;
    const double *rpc = self->rpc;
    const double nrmR2inv = self->nrmR2inv;
    const int xomax = self->xomax;
    const int yomax = self->yomax;
    double **xypc = self->xypc;
    const double nrmXinv = self->nrmXinv;
    const double nrmYinv = self->nrmYinv;
    const double iter_max = self->itermax;
    const double tol = self->tol;
    int i;
    double asphz, dpz;
    double z_coni, zx_coni, zy_coni, z_rp, zx_rp, zy_rp, z_xyp, zx_xyp, zy_xyp, dmy;
    double nx, ny, stp, ilen;
    double px = pos[X];
    double py = pos[Y];
    double pz = pos[Z];
    double dx = dir[X];
    double dy = dir[Y];
    double dz = dir[Z];

    for (i = 0; i < iter_max; i++) {
        if (UpcShapeSph_conicZ(curv, coni, px, py, 1, &z_coni, &zx_coni, &zy_coni, &dmy, &dmy, &dmy)) {
            UpcShapeAsph_rpolyZ(romax, rpc, px, py, nrmR2inv, 1, &z_rp, &zx_rp, &zy_rp, &dmy, &dmy, &dmy);
            UpcShapeXYP_xypolyZ(xomax, yomax, xypc, px, py, nrmXinv, nrmYinv, 1, &z_xyp, &zx_xyp, &zy_xyp, &dmy, &dmy, &dmy);
            asphz = z_coni + z_rp + z_xyp;
            nx = -(zx_coni + zx_rp + zx_xyp);
            ny = -(zy_coni + zy_rp + zy_xyp);
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
    int ix, iy;

    for (iy = 1; iy <= self->yomax; iy += 2) {
        for (ix = 0; ix <= self->xomax; ix++) {
            if (self->xypc[iy][ix])
                return FALSE;
        }
    }
    return TRUE;
}

static Bool ySymmetric(const UpcShape *self)
{
    int ix, iy;

    for (iy = 0; iy <= self->yomax; iy++) {
        for (ix = 1; ix <= self->xomax; ix += 2) {
            if (self->xypc[iy][ix])
                return FALSE;
        }
    }
    return TRUE;
}

static Bool rSymmetric(const UpcShape *self)
{
    return FALSE;
}

/* xy多項式のz座標と微分を求める.
 * 入力
 *     max_orderX : 多項式のX最大次数
 *     max_orderY : 多項式のY最大次数
 *     xypc[0..max_orderY][0..max_orderX] : 多項式の係数
 *     x : x座標
 *     y : y座標
 *     nrmXinv : X正規化係数の逆数
 *     nrmYinv : Y正規化係数の逆数
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
Bool UpcShapeXYP_xypolyZ(int max_orderX, int max_orderY, double **xypc, double x, double y, double nrmXinv, double nrmYinv, int mode, double *z, double *zx, double *zy, double *zxx, double *zxy, double *zyy)
{
    double xdsh = x * nrmXinv;
    double ydsh = y * nrmYinv;
    int ix, iy;
    const double *xpc;

    if (mode == 0) {
        double yc, zt;

        zt = 0.0;
        for (iy = max_orderY; iy >= 0; iy--) {
            xpc = xypc[iy];
            yc = xpc[max_orderX];
            for (ix = max_orderX - 1; ix >= 0; ix--)
                yc = yc * xdsh + xpc[ix];
            zt = zt * ydsh + yc;
        }
        *z = zt;
    }
    else if (mode == 1) {
        double yc, yxc, zt, zxt, zyt;

        zt = zxt = zyt = 0.0;
        for (iy = max_orderY; iy >= 0; iy--) {
            xpc = xypc[iy];
            yc = yxc = 0.0;
            for (ix = max_orderX; ix > 0; ix--) {
                yc = yc * xdsh + xpc[ix];
                yxc = yxc * xdsh + ix * xpc[ix];
            }
            yc = yc * xdsh + xpc[0]; // 最後の一回
            zt = zt * ydsh + yc;
            zxt = zxt * ydsh + yxc;
            if (iy)
                zyt = zyt * ydsh + iy * yc;
        }
        *z = zt;
        *zx = zxt * nrmXinv;
        *zy = zyt * nrmYinv;
    }
    else if (mode == 2) {
        double yc, yxc, yxxc;
        double zt, zxt, zyt, zxxt, zxyt, zyyt;
        
        zt = zxt = zyt = zxxt = zxyt = zyyt = 0.0;
        for (iy = max_orderY; iy >= 0; iy--) {
            xpc = xypc[iy];
            yc = yxc = yxxc = 0.0;
            for (ix = max_orderX; ix > 1; ix--) {
                yc = yc * xdsh + xpc[ix];
                yxc = yxc * xdsh + ix * xpc[ix];
                yxxc = yxxc * xdsh + ix * (ix - 1) * xpc[ix];
            }
            if (max_orderX) {
                yc = yc * xdsh + xpc[1];
                yxc = yxc * xdsh + ix * xpc[1];
            }
            yc = yc * xdsh + xpc[0]; // 最後の一回
            zt = zt * ydsh + yc;
            zxt = zxt * ydsh + yxc;
            zxxt = zxxt * ydsh + yxxc;
            if (iy) {
                zyt = zyt * ydsh + iy * yc;
                zxyt = zxyt * ydsh + iy * yxc;
                if (iy - 1)
                    zyyt = zyyt * ydsh + iy * (iy - 1) * yc;
            }
        }
        *z = zt;
        *zx = zxt * nrmXinv;
        *zy = zyt * nrmYinv;
        *zxx = zxxt * nrmXinv * nrmXinv;
        *zxy = zxyt * nrmXinv * nrmYinv;
        *zyy = zyyt * nrmYinv * nrmYinv;
    }
    return TRUE;
}

/*--------------------------------------------------------------------------------*/
/* UpcShapeZernikeオブジェクト */

static void UpcShapeZernike_dealloc(UpcShape *self)
{
    free(self->zc);
    upc_freeD2(self->xypc);
    free(self->rpc);
    free(self);
}

static void UpcShapeZernike_print(const UpcShape *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);
    int i;
    
    fprintf(ostream, "%stype    = %s\n", space, UpcShapeType_string(self->shapeType));
    fprintf(ostream, "%sr    = %13.6e\n", space, UpcShapeSph_RADI(self));
    fprintf(ostream, "%sconi = %13.6e\n", space, self->coni);
    for (i = 2; i <= self->romax; i++) {
        if (self->rpc[i])
            fprintf(ostream, "%srpc[%2d] = %13.6e\n", space, i, self->rpc[i]);
    }
    fprintf(ostream, "%snrmR    = %13.6e\n", space, UpcShapeAsph_nrmR(self));
    fprintf(ostream, "%sitermax = %d\n", space, self->itermax);
    fprintf(ostream, "%stol     = %13.6e\n", space, self->tol);

    for (i = 1; i <= self->zcnmax; i++) {
        if (self->zc[i])
            fprintf(ostream, "%szc[%2d] = %13.6e\n", space, i, self->zc[i]);
    }
    fprintf(ostream, "%snrmZR   = %13.6e\n", space, UpcShapeZernike_nrmZR(self));
}

static void UpcShapeZernike_toXml(const UpcShape *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);
    int i;
    
    fprintf(ostream, "%s<shape type=\"%s\" r_order_max=\"%d\" zcn_max=\"%d\">\n", space, UpcShapeType_string(self->shapeType), self->romax, self->zcnmax);
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

    for (i = 1; i <= self->zcnmax; i++) {
        if (self->rpc[i] != coefDefault)
            fprintf(ostream, "%s  <zcoef cno=\"%d\">%23.15e</rcoef>\n", space, i, self->zc[i]);
    }
    if (self->nrmZRinv != nrmZRinvDefault)
        fprintf(ostream, "%s  <nrm_zr>%23.15e</nrm_zr>\n", space, UpcShapeZernike_nrmZR(self));

    fprintf(ostream, "%s</shape>\n", space);
}

static UpcErrStatus UpcShapeZernike_lset(UpcShape *self)
{
    int i;

    UpcShapeXYP_clearXYPC(self);
    for (i = 1; i <= self->zcnmax; i++) {
        if (self->zc[i]) {
            UpcShapeXYP_addZcoef(self, i, self->zc[i]);
        }
    }
    return UpcE_NoIssues;
}

static UpcShape *UpcShapeZernike_copy(const UpcShape *self)
{
    UpcShape *t = UpcShapeZernike_init(self->romax, self->zcnmax);
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

        // t->zcnmax = self->zcnmax;
        t->nrmZRinv = self->nrmZRinv;
        for (i = 0; i <= t->zcnmax; i++)
            t->zc[i] = self->zc[i];
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeZernike_copy");
    }
    return t;
}

UpcShape *UpcShapeZernike_init(int max_orderR, int zcn_max)
{
    UpcShape *self;
    int xomax, yomax;
    
    if (max_orderR < 2 || zcn_max < 1) {
        UpcERRHANDLER(UpcE_ValueError, "invalid parameter (in UpcShapeZernike_init)");
        return NULL;
    }

    xomax = yomax = UpcZernikeToXYP_maxorder(zcn_max);

    self = UpcShapeXYP_init(max_orderR, xomax, yomax);
    if (self) {
        self->zc = (double *)calloc(zcn_max + 1, sizeof(double));
        if (!self->zc) {
            upc_freeD2(self->xypc);
            free(self->rpc);
            free(self);
            UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeZernike_init");
            return NULL;
        }
        if (coefDefault) {
            int i;
            
            for (i = 1; i <= zcn_max; i++)
                self->zc[i] = coefDefault;
        }
        self->_dealloc = UpcShapeZernike_dealloc;

        self->shapeType = UpcShapeType_zernike;
        self->zcnmax = zcn_max;
        self->copy = UpcShapeZernike_copy;
        self->print = UpcShapeZernike_print;
        self->toXml = UpcShapeZernike_toXml;
        self->lset = UpcShapeZernike_lset;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeZernike_init");
    }
    return self;
}

UpcShape *UpcShapeZernike_initWithShape(int max_orderR, int zcn_max, UpcShape *other)
{
    UpcShape *self = UpcShapeZernike_init(max_orderR, zcn_max);
    enum UpcShapeType otype = other->shapeType;
    
    if (self) {
        if (otype % UpcShapeType_sph == 0) { // UpcShapeSphクラスまたは子クラスであれば
            self->curv = other->curv;
        }
        if (otype % UpcShapeType_asph == 0) { // UpcShapeAsphクラスまたは子クラスであれば
            int i, n = (max_orderR < other->romax) ? max_orderR : other->romax;
            
            self->coni = other->coni;
            for (i = 2; i <= n; i++) {
                self->rpc[i] = other->rpc[i];
            }
            self->nrmR2inv = other->nrmR2inv;
            self->itermax = other->itermax;
            self->tol = other->tol;
        }
        if (otype % UpcShapeType_zernike == 0) { // UpcShapeZernikeクラスまたは子クラスであれば
            int i;
            int zcnmax = (zcn_max < other->zcnmax) ? zcn_max : other->zcnmax;
            
            for (i = 1; i <= zcnmax; i++) {
                self->zc[i] = other->zc[i];
            }
            self->nrmZRinv = other->nrmZRinv;
        }
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcShapeZernike_initWithShape");
    }
    return self;
}


double UpcShapeZernike_nrmZR(const UpcShape *self)
{
    return 1.0 / self->nrmZRinv;
}

Bool UpcShapeZernike_setNrmZR(UpcShape *self, double nrmZR)
{
    if (nrmZR) {
        self->nrmZRinv = 1.0 / nrmZR;
        return TRUE;
    }
    return FALSE;
}





