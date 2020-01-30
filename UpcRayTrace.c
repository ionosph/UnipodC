/*
 *  UpcRayTrace.c
 *  UnipodC
 *
 *  Created by ionosph on 2011/04/05.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcRayTrace.h"

static const int ITER_MAX = 100; // 光線サーチ計算の最大繰り返し数
static const double HugeVal = 1.0e30;

/* Snellの法則にしたがって光線の屈折を計算する
 * [入力]
 *     vi[3] : 入射光線の方向余弦
 *     vn[3] : 面の単位法線ベクトル
 *     ni : 面の前の屈折率
 *     no : 面の後の屈折率
 * [出力]
 *     vo[3] : 射出光線の方向余弦
 *     cosj : 入射角のcosθ
 *     cosk : 射出角のcosθ
 *     戻り値 : 正常に計算できたらTRUE, 全反射した場合はFALSE.
 */
static Bool UpcRayTrace_snell_refraction(const double *vi, const double *vn, double ni, double no, double *vo, double *cosj, double *cosk)
{
    double nd, t, gtld;
    
    *cosj = vi[X] * vn[X] + vi[Y] * vn[Y] + vi[Z] * vn[Z];
    nd = fabs(ni / no);
    t = 1.0 - nd * nd * (1.0 - *cosj * *cosj);
    if (t >= 0.0) {
        *cosk = (ni * no * *cosj > 0.0) ? sqrt(t) : -sqrt(t);
        gtld = *cosk - nd * *cosj;
        vo[X] = nd * vi[X] + gtld * vn[X];
        vo[Y] = nd * vi[Y] + gtld * vn[Y];
        vo[Z] = nd * vi[Z] + gtld * vn[Z];
        return TRUE;
    }
    return FALSE; // Total Reflection
}

/* 第j面の光線から次の第k面の光線を計算する.(順トレース)
 * [入力]
 *     sj : 第j面を表すUpcSurfオブジェクト
 *     sk : 第k面を表すUpcSurfオブジェクト
 *     rj : 第j面の光線を表すUpcRayComponentオブジェクト
 *     iwav : 波長番号
 * [出力]
 *     rk->pos  : 第k面光線の位置
 *     rk->dir  : 第k面光線の射出方向余弦
 *     rk->cosi : 第k面光線の入射角の余弦
 *     rk->coso : 第k面光線の射出角の余弦
 *     rj->opl  : 第j面から第k面までの光路に沿った距離
 *     戻り値 : 光線追跡ステータス
 */
static UpcRayStatus UpcRayTrace_cal_next(const UpcSurf *sj, const UpcSurf *sk, UpcRayComponent *rj, UpcRayComponent *rk, int iwav)
{
    double rjpos[3], rjdir[3], nv[3], t[3];
    Bool flg, isActive;

    // rj->pos, rj->dirをsk座標で表現
    t[X] = rj->pos[X];
    t[Y] = rj->pos[Y];
    t[Z] = rj->pos[Z] - sj->d;
    if (sk->dec->decenterFlg) {
        UpcDecenter_convForward(sk->dec, t, rjpos);
        UpcDecenter_rotForward(sk->dec, rj->dir, rjdir);
    }
    else {
        UpcVector_COPY(rjpos, t);
        UpcVector_COPY(rjdir, rj->dir);
    }
    // 面と光線の交点を求める
    flg = sk->shape->getIntersection(sk->shape, rjpos, rjdir, sk->area, rk->pos, nv, &isActive);
    if (!flg)
        return rk->status = UpcRayE_NoIntersectionError;
    // 射出光線の方向余弦を求める
    flg = UpcRayTrace_snell_refraction(rjdir, nv, sj->_n[iwav], sk->_n[iwav], rk->dir, &(rk->cosI), &(rk->cosO));
    if (!flg)
        return rk->status = UpcRayE_TotalReflectionError;
    // 光路に沿った距離を求める
    rj->opl = (rk->pos[Z] - rjpos[Z]) / rjdir[Z];
    return rk->status = isActive ? UpcRayE_Succeeded : UpcRayE_OutOfArea;
}

/* 物体面上座標とターゲット座標を与えて光線追跡する.
 * [入力]
 *     lobj : UpcLensオブジェクト
 *     s0pos[3] : 物体面上座標
 *     tpos[3] : ターゲット座標 (物体面座標系)
 *     iwav : 波長番号
 *     isur_end : この面まで光線追跡する.(0にすると像面まで)
 * [出力]
 *     ray : UpcRayオブジェクト (r->nsur == lobj->nsur でなければならない)
 *     戻り値 : isur_end面まで光線が通ればTRUE, 途中で失敗したらFALSE.
 */
Bool UpcRayTrace_trace0(const UpcLens *lobj, const double *s0pos, const double *tpos, int iwav, int isur_end, UpcRay *ray)
{
    int i;
    double v[3], len2, ilen;
    UpcSurf *sk, *sj = lobj->s[0];
    UpcRayComponent *rk, *rj = ray->s;
    UpcRayStatus status;

    if (!isur_end)
        isur_end = lobj->nsur + 1;
    UpcRay_resetStatus(ray);
    UpcVector_COPY(ray->s[0].pos, s0pos);
    UpcVector_SUB(v, tpos, s0pos);
    len2 = v[X] * v[X] + v[Y] * v[Y] + v[Z] * v[Z];
    if (!len2) {
        ray->status = UpcRayE_SettingError;
        return FALSE;
    }
    ilen = 1.0 / sqrt(len2);
    ray->s[0].dir[X] = ilen * v[X];
    ray->s[0].dir[Y] = ilen * v[Y];
    ray->s[0].dir[Z] = ilen * v[Z];
    ray->status = ray->s[0].status = UpcArea_isActive(sj->area, s0pos[X], s0pos[Y]) ? UpcRayE_Succeeded : UpcRayE_OutOfArea;
    for (i = 1; i <= isur_end; i++) {
        sk = lobj->s[i];
        rk = ray->s + i;
        status = UpcRayTrace_cal_next(sj, sk, rj, rk, iwav);
        if (status == UpcRayE_NoIntersectionError || status == UpcRayE_TotalReflectionError) {
            ray->status = status;
            return FALSE;
        }
        else if (status == UpcRayE_OutOfArea) {
            ray->status = status;
        }
        sj = sk;
        rj = rk;
    }
    return TRUE;
}

/* 物体面を出射する光線の座標と方向余弦を与えて光線追跡する.
 * [入力]
 *     lobj : UpcLensオブジェクト
 *     s0pos[3] : 物体面上座標
 *     s0dir[3] : 出射方向余弦
 *     iwav : 波長番号
 *     isur_end : この面まで光線追跡する.(0にすると像面まで)
 * [出力]
 *     ray : UpcRayオブジェクト (ray->nsur == lobj.nsur でなければならない)
 *     戻り値 : isur_end面まで光線が通ればTRUE, 途中で失敗したらFALSE.
 */
Bool UpcRayTrace_trace0d(const UpcLens *lobj, const double *s0pos, const double *s0dir, int iwav, int isur_end, UpcRay *ray)
{
    int i;
    UpcSurf *sk, *sj = lobj->s[0];
    UpcRayComponent *rk, *rj = ray->s;
    UpcRayStatus status;

    if (!isur_end)
        isur_end = lobj->nsur + 1;
    UpcRay_resetStatus(ray);
    UpcVector_COPY(ray->s[0].pos, s0pos);
    UpcVector_COPY(ray->s[0].dir, s0dir);
    ray->status = ray->s[0].status = UpcArea_isActive(sj->area, s0pos[X], s0pos[Y]) ? UpcRayE_Succeeded : UpcRayE_OutOfArea;
    for (i = 1; i <= isur_end; i++) {
        sk = lobj->s[i];
        rk = ray->s + i;
        status = UpcRayTrace_cal_next(sj, sk, rj, rk, iwav);
        if (status == UpcRayE_NoIntersectionError || status == UpcRayE_TotalReflectionError) {
            ray->status = status;
            return FALSE;
        }
        else if (status == UpcRayE_OutOfArea) {
            ray->status = status;
        }
        sj = sk;
        rj = rk;
    }
    return TRUE;
}

/* 主光線の物体面上座標と出射方向余弦, および瞳座標を与えて光線追跡する.
 * [入力]
 *     lobj : UpcLensオブジェクト
 *     s0pos[3] : 物体面上座標
 *     s0dir[3] : 基準光線(主光線, 瞳中心光線)の出射方向余弦
 *     sinTx, sinTy : 瞳座標(sinθx, sinθy)
 *     iwav : 波長番号
 * [出力]
 *     ray : UpcRayオブジェクト (ray->nsur == lobj.nsur でなければならない)
 *     戻り値 : 成功すればTRUE, 途中で失敗したらFALSE.
 */
Bool UpcRayTrace_trace0p(const UpcLens *lobj, const double *s0pos, const double *s0dir, double sinTx, double sinTy, int iwav, UpcRay *ray)
{
    double tdir[3], sq;
    
    tdir[X] = s0dir[X] + sinTx;
    tdir[Y] = s0dir[Y] + sinTy;
    sq = tdir[X] * tdir[X] + tdir[Y] * tdir[Y];
    if (sq >= 1.0) {
        UpcRay_resetStatus(ray);
        ray->status = UpcRayE_SettingError;
        return FALSE;
    }
    tdir[Z] = sqrt(1.0 - sq);
    return UpcRayTrace_trace0d(lobj, s0pos, tdir, iwav, 0, ray);
}

/* 方向余弦とターゲット座標を与えて光線追跡する.
 * [入力]
 *     lobj : UpcLensオブジェクト
 *     s0dir[3] : 方向余弦 (物体面座標系)
 *     tpos[3] : ターゲット座標 (第1面座標系. ただし第1面が偏心している場合は偏心前の座標系)
 *     iwav : 波長番号
 *     isur_end : この面まで光線追跡する.(0にすると像面まで)
 * [出力]
 *     ray : UpcRayオブジェクト (ray->nsur == lobj.nsur でなければならない)
 *     戻り値 : isur_end面まで光線が通ればTRUE, 途中で失敗したらFALSE.
 */
Bool UpcRayTrace_trace1(const UpcLens *lobj, const double *s0dir, const double *tpos, int iwav, int isur_end, UpcRay *ray)
{
    int i, flg;
    double p[3], d[3], dmmy[3], nv[3];
    UpcSurf **sur = lobj->s;
    UpcSurf *sk, *sj = sur[1];
    UpcRayComponent *rk, *rj = ray->s + 1;
    Bool isActive;
    UpcRayStatus status;

    if (!isur_end)
        isur_end = lobj->nsur + 1;
    UpcRay_resetStatus(ray);
    UpcVector_COPY(ray->s[0].dir, s0dir);
    p[X] = tpos[X];
    p[Y] = tpos[Y];
    p[Z] = tpos[Z] + sur[0]->d;
    flg = sur[0]->shape->getIntersection(sur[0]->shape, p, ray->s[0].dir, sur[0]->area, ray->s[0].pos, dmmy, &isActive);
    if (!flg) {
        ray->status = ray->s[0].status = UpcRayE_NoIntersectionError;
        return FALSE;
    }
    ray->status = ray->s[0].status = isActive ? UpcRayE_Succeeded : UpcRayE_OutOfArea;
    if (sj->dec->decenterFlg) {
        UpcDecenter_convForward(sj->dec, tpos, p);
        UpcDecenter_rotForward(sj->dec, s0dir, d);
    }
    else {
        UpcVector_COPY(p, tpos);
        UpcVector_COPY(d, s0dir);
    }
    flg = sj->shape->getIntersection(sj->shape, p, d, sj->area, ray->s[1].pos, nv, &isActive);
    if (!flg) {
        ray->status = ray->s[1].status = UpcRayE_NoIntersectionError;
        return FALSE;
    }
    flg = UpcRayTrace_snell_refraction(d, nv, sur[0]->_n[iwav], sj->_n[iwav], ray->s[1].dir, &(ray->s[1].cosI), &(ray->s[1].cosO));
    if (!flg) {
        ray->status = ray->s[1].status = UpcRayE_TotalReflectionError;
        return FALSE;
    }
    if (isActive)
        ray->s[1].status = UpcRayE_Succeeded;
    else
        ray->status = ray->s[1].status = UpcRayE_OutOfArea;

    for (i = 2; i <= isur_end; i++) {
        sk = sur[i];
        rk = ray->s + i;
        status = UpcRayTrace_cal_next(sj, sk, rj, rk, iwav);
        if (status == UpcRayE_NoIntersectionError || status == UpcRayE_TotalReflectionError) {
            ray->status = status;
            return FALSE;
        }
        else if (status == UpcRayE_OutOfArea) {
            ray->status = status;
        }
        sj = sk;
        rj = rk;
    }
    return TRUE;
}

/* 像面に入射する光線の座標と方向余弦を与えて逆トレースする.
 * [入力]
 *     lobj : UpcLensオブジェクト
 *     sipos[3] : 像面上座標
 *     skdir[3] : 像面に入射する光線の方向余弦
 *     iwav : 波長番号
 * [出力]
 *     ray : UpcRayオブジェクト (ray->nsur == lobj.nsur でなければならない)
 *     戻り値 : 光線が通ればTRUE, 失敗したらFALSE.
 */
Bool UpcRayTrace_trace_rev(const UpcLens *lobj, const double *sipos, const double *skdir, int iwav, UpcRay *ray)
{
    const int isur_img = lobj->nsur + 1;
    UpcSurf **sur = lobj->s;
    UpcSurf *sj, *sk = sur[isur_img];
    UpcRayComponent *rj, *rk = ray->s + isur_img;
    double nv[3], rkpos[3], t[3];
    int j, flg;
    Bool isActive;

    UpcRay_resetStatus(ray);
    UpcVector_COPY(rk->pos, sipos);
    UpcVector_COPY(rk->dir, skdir);
    if (!(sk->shape->getNormalVector(sk->shape, rk->pos[X], rk->pos[Y], nv))) {
        ray->status = rk->status = UpcRayE_NoIntersectionError;
        return FALSE;
    }
    ray->status = rk->status = UpcArea_isActive(sk->area, rk->pos[X], rk->pos[Y]) ? UpcRayE_Succeeded : UpcRayE_OutOfArea;

    for (j = lobj->nsur; j >= 0; j--) {
        sj = sur[j];
        rj = ray->s + j;
        // 入射光線の方向余弦を求める
        flg = UpcRayTrace_snell_refraction(rk->dir, nv, sk->_n[iwav], sj->_n[iwav], rj->dir, &(rk->cosO), &(rk->cosI)); // (注) 前のステップで求めたnvを利用
        if (!flg) {
            ray->status = rj->status = UpcRayE_TotalReflectionError;
            return FALSE;
        }
        // rk->posをsj座標で表現
        if (sk->dec->decenterFlg) {
            UpcDecenter_convReverse(sk->dec, rk->pos, rkpos);
            rkpos[Z] += sj->d;
            UpcDecenter_rotReverse(sk->dec, rj->dir, t);
            UpcVector_COPY(rj->dir, t);
        }
        else {
            UpcVector_SET(rkpos, rk->pos[X], rk->pos[Y], rk->pos[Z] + sj->d);
        }
        // 面と光線の交点を求める
        flg = sj->shape->getIntersection(sj->shape, rkpos, rj->dir, sj->area, rj->pos, nv, &isActive); // (注) nvは次のステップで再利用する.
        if (!flg) {
            ray->status = rj->status = UpcRayE_NoIntersectionError;
            return FALSE;
        }
        // 光路に沿った距離を求める
        rj->opl = (rkpos[Z] - rj->pos[Z]) / rj->dir[Z];
        if (isActive)
            rj->status = UpcRayE_Succeeded;
        else
            ray->status = rj->status = UpcRayE_OutOfArea;
        sk = sj;
        rk = rj;
    }
    return TRUE;
}

/* 物体面上(ox, oy)から射出し,第isur面上(x, y)を通る光線を探す.
 * [入力]
 *     lobj : UpcLensオブジェクト
 *     ox, oy, isur, x, y : 光線通過条件
 *     iwav : 波長番号
 *     tol : 収束精度(同一点とみなす距離)
 * [出力]
 *     ray : UpcRayオブジェクト (ray->nsur == lobj.nsur でなければならない)
 *     戻り値 : 光線が見つかればTRUE, 失敗したらFALSE.
 */
Bool UpcRayTrace_psearch0(const UpcLens *lobj, double ox, double oy, int isur, double x, double y, int iwav, double tol, UpcRay *ray)
{
    UpcSurf **sur = lobj->s;
    double s0pos[3], tpos[3], t[3];
    double c[4], ekdsh, hk, akdsh, h1, a1, enpz;
    int i;
    double dx, dy;
    double A, B, C, D, E, F, det, idet, x0, y0, sumsq_t;
    
    if (isur < 1 || isur > lobj->nsur) {
        ray->status = UpcRayE_SettingError;
        return FALSE;
    }
    s0pos[X] = ox;
    s0pos[Y] = oy;
    if (!(sur[0]->shape->getZ(sur[0]->shape, ox, oy, s0pos + Z))) {
        ray->status = ray->s[0].status = UpcRayE_NoIntersectionError;
        return FALSE;
    }

    if (isur == 1) {
        t[X] = x;
        t[Y] = y;
        if (!(sur[1]->shape->getZ(sur[1]->shape, x, y, t + Z))) {
            ray->status = ray->s[1].status = UpcRayE_NoIntersectionError;
            return FALSE;
        }
        UpcDecenter_convReverse(sur[1]->dec, t, tpos);
        tpos[Z] += sur[0]->d;
        return UpcRayTrace_trace0(lobj, s0pos, tpos, iwav, 0, ray);
    }
    // 物体から見た指定面の像の位置を求める.
    if (isur == lobj->isur_sto) {
        enpz = sur[0]->d + lobj->parax[iwav].enp;
    }
    else {
        UpcLens_paraxMatrix(lobj, 1, isur - 1, iwav, c);
        UpcMatrix2_inverse(c);
        ekdsh = UpcSurf_e(sur[isur - 1], iwav);
        if (ekdsh) {
            hk = 1.0;
            akdsh = hk / ekdsh;
        }
        else {
            hk = 0.0;
            akdsh = 1.0;
        }
        h1 = c[0] * hk + c[1] * akdsh;
        a1 = c[2] * hk + c[3] * akdsh;
        if (fabs(a1) < 1.0e-6)
            enpz = 1.0e6;
        else
            enpz = sur[0]->d + sur[0]->_n[iwav] * h1 / a1;
    }
    // A, B, C, D, E, F の初期値を求める.
    UpcVector_SET(tpos, 0.0, 0.0, enpz);
    if (!UpcRayTrace_trace0(lobj, s0pos, tpos, iwav, isur, ray))
        return FALSE;
    E = ray->s[isur].pos[X];
    F = ray->s[isur].pos[Y];
    tpos[X] = 1.0;
    tpos[Y] = 0.0;
    if (!UpcRayTrace_trace0(lobj, s0pos, tpos, iwav, isur, ray))
        return FALSE;
    A = ray->s[isur].pos[X] - E;
    C = ray->s[isur].pos[Y] - F;
    tpos[X] = 0.0;
    tpos[Y] = 1.0;
    if (!UpcRayTrace_trace0(lobj, s0pos, tpos, iwav, isur, ray))
        return FALSE;
    B = ray->s[isur].pos[X] - E;
    D = ray->s[isur].pos[Y] - F;
    // 繰り返し計算により光線をサーチする.
    x0 = x - E;
    y0 = y - F;
    for (i = 0; i < ITER_MAX; i++) {
        det = A * D - B * C;
        if (!det) {
            ray->status = UpcRayE_IterationError;
            return FALSE;
        }
        idet = 1.0 / det;
        tpos[X] = idet * ( D * x0 - B * y0);
        tpos[Y] = idet * (-C * x0 + A * y0);
        if (!UpcRayTrace_trace0(lobj, s0pos, tpos, iwav, isur, ray))
            return FALSE;
        dx = ray->s[isur].pos[X] - x;
        dy = ray->s[isur].pos[Y] - y;
        if (fabs(dx) < tol && fabs(dy) < tol) {
            if (tpos[Z] - s0pos[Z] < 0.0) {
                tpos[X] = 2.0 * s0pos[X] - tpos[X];
                tpos[Y] = 2.0 * s0pos[Y] - tpos[Y];
                tpos[Z] = 2.0 * s0pos[Z] - tpos[Z];
            }
            return UpcRayTrace_trace0(lobj, s0pos, tpos, iwav, 0, ray);
        }
        sumsq_t = tpos[X] * tpos[X] + tpos[Y] * tpos[Y];
        if (!sumsq_t) {
            ray->status = UpcRayE_IterationError;
            return FALSE;
        }
        sumsq_t = 1.0 / sumsq_t;
        A += tpos[X] * dx * sumsq_t;
        B += tpos[Y] * dx * sumsq_t;
        C += tpos[X] * dy * sumsq_t;
        D += tpos[Y] * dy * sumsq_t;
    }
    ray->status = UpcRayE_IterationError;
    return FALSE;
}

/* 第1面に方向余弦(d1x, d1y)で入射し,第isur面上(x, y)を通る光線を探す.
 * [入力]
 *     lobj : UpcLensオブジェクト
 *     d1x, d1y, isur, x, y : 光線通過条件
 *     iwav : 波長番号
 *     tol : 収束精度(同一点とみなす距離)
 * [出力]
 *     ray : UpcRayオブジェクト (ray->nsur == lobj.nsur でなければならない)
 *     戻り値 : 光線が見つかればTRUE, 失敗したらFALSE.
 */
Bool UpcRayTrace_psearch1(const UpcLens *lobj, double d1x, double d1y, int isur, double x, double y, int iwav, double tol, UpcRay *ray)
{
    UpcSurf **sur = lobj->s;
    double t = 1.0 - d1x * d1x - d1y * d1y;
    double s0dir[3], tpos[3], p[3];
    double c[4], ekdsh, hk, akdsh, h1, a1, enpz;
    int i;
    double dx, dy;
    double A, B, C, D, E, F, det, idet, x0, y0, sumsq_t;
    
    if (isur < 1 || isur > lobj->nsur) {
        ray->status = UpcRayE_SettingError;
        return FALSE;
    }
    if (t <= 0.0){
        ray->status = UpcRayE_SettingError;
        return FALSE;  // 方向余弦初期値が正しくない.
    }
    UpcVector_SET(s0dir, d1x, d1y, sqrt(t));
    if (isur == 1) {
        p[X] = x;
        p[Y] = y;
        if (!(sur[1]->shape->getZ(sur[1]->shape, x, y, p + Z))) {
            ray->status = ray->s[1].status = UpcRayE_NoIntersectionError;
            return FALSE;
        }
        UpcDecenter_convReverse(sur[1]->dec, p, tpos);
        return UpcRayTrace_trace1(lobj, s0dir, tpos, iwav, 0, ray);
    }
    // 第1面から見た指定面の像の位置を求める.
    if (isur == lobj->isur_sto) {
        enpz = lobj->parax[iwav].enp;
    }
    else {
        UpcLens_paraxMatrix(lobj, 1, isur - 1, iwav, c);
        UpcMatrix2_inverse(c);
        ekdsh = UpcSurf_e(sur[isur - 1], iwav);
        if (ekdsh) {
            hk = 1.0;
            akdsh = hk / ekdsh;
        }
        else {
            hk = 0.0;
            akdsh = 1.0;
        }
        h1 = c[0] * hk + c[1] * akdsh;
        a1 = c[2] * hk + c[3] * akdsh;
        if (fabs(a1) < 1.0e-6)
            enpz = 1.0e6;
        else
            enpz = sur[0]->_n[iwav] * h1 / a1;
    }
    // A, B, C, D, E, F の初期値を求める.
    UpcVector_SET(tpos, 0.0, 0.0, enpz);
    if (!(UpcRayTrace_trace1(lobj, s0dir, tpos, iwav, isur, ray)))
        return FALSE;
    E = ray->s[isur].pos[X];
    F = ray->s[isur].pos[Y];
    tpos[X] = 1.0;
    tpos[Y] = 0.0;
    if (!(UpcRayTrace_trace1(lobj, s0dir, tpos, iwav, isur, ray)))
        return FALSE;
    A = ray->s[isur].pos[X] - E;
    C = ray->s[isur].pos[Y] - F;
    tpos[X] = 0.0;
    tpos[Y] = 1.0;
    if (!(UpcRayTrace_trace1(lobj, s0dir, tpos, iwav, isur, ray)))
        return FALSE;
    B = ray->s[isur].pos[X] - E;
    D = ray->s[isur].pos[Y] - F;
    // 繰り返し計算により光線をサーチする.
    x0 = x - E;
    y0 = y - F;
    for (i = 0; i < ITER_MAX; i++) {
        det = A * D - B * C;
        if (!det) {
            ray->status = UpcRayE_IterationError;
            return FALSE;
        }
        idet = 1.0 / det;
        tpos[X] = idet * ( D * x0 - B * y0);
        tpos[Y] = idet * (-C * x0 + A * y0);
        if (!(UpcRayTrace_trace1(lobj, s0dir, tpos, iwav, isur, ray)))
            return FALSE;
        dx = ray->s[isur].pos[X] - x;
        dy = ray->s[isur].pos[Y] - y;
        if (fabs(dx) < tol && fabs(dy) < tol) {
            return UpcRayTrace_trace1(lobj, s0dir, tpos, iwav, 0, ray);
        }
        sumsq_t = tpos[X] * tpos[X] + tpos[Y] * tpos[Y];
        if (!sumsq_t) {
            ray->status = UpcRayE_IterationError;
            return FALSE;
        }
        sumsq_t = 1.0 / sumsq_t;
        A += tpos[X] * dx * sumsq_t;
        B += tpos[Y] * dx * sumsq_t;
        C += tpos[X] * dy * sumsq_t;
        D += tpos[Y] * dy * sumsq_t;
    }
    ray->status = UpcRayE_IterationError;
    return FALSE;
}

/* 物体面上(ox, oy)から射出し,射出方向余弦が(cosx, cosy)になる光線を探す.
 * [入力]
 *     lobj : UpcLensオブジェクト
 *     ox, oy, cosx, cosy : 光線通過条件
 *     iwav : 波長番号
 *     tol : 収束精度(同一点とみなす距離)
 * [出力]
 *     ray : UpcRayオブジェクト (ray->nsur == lobj.nsur でなければならない)
 *     戻り値 : 光線が見つかればTRUE, 失敗したらFALSE.
 */
Bool UpcRayTrace_psearch0_xt(const UpcLens *lobj, double ox, double oy, double cosx, double cosy, int iwav, double tol, UpcRay *ray)
{
    int isur_img = lobj->nsur + 1;
    UpcSurf *si = lobj->s[isur_img];
    double skdir[3], sipos[3];
    int i;
    double dx, dy;
    double A, B, C, D, E, F, det, idet, x0, y0, sumsq_t;
    double t = 1.0 - cosx * cosx - cosy * cosy;

    if (t <= 0.0) {
        ray->status = UpcRayE_SettingError;
        return FALSE;
    }
    UpcVector_SET(skdir, cosx, cosy, sqrt(t));
    // A, B, C, D, E, F の初期値を求める.
    sipos[X] = 0.0;
    sipos[Y] = 0.0;
    if (!(si->shape->getZ(si->shape, sipos[X], sipos[Y], sipos + Z))) {
        ray->status = ray->s[isur_img].status = UpcRayE_NoIntersectionError;
        return FALSE;
    }
    if (!UpcRayTrace_trace_rev(lobj, sipos, skdir, iwav, ray))
        return FALSE;
    E = ray->s[0].pos[X];
    F = ray->s[0].pos[Y];
    sipos[X] = 1.0;
    sipos[Y] = 0.0;
    if (!(si->shape->getZ(si->shape, sipos[X], sipos[Y], sipos + Z))) {
        ray->status = ray->s[isur_img].status = UpcRayE_NoIntersectionError;
        return FALSE;
    }
    if (!UpcRayTrace_trace_rev(lobj, sipos, skdir, iwav, ray))
        return FALSE;
    A = ray->s[0].pos[X] - E;
    C = ray->s[0].pos[Y] - F;
    sipos[X] = 0.0;
    sipos[Y] = 1.0;
    if (!(si->shape->getZ(si->shape, sipos[X], sipos[Y], sipos + Z))) {
        ray->status = ray->s[isur_img].status = UpcRayE_NoIntersectionError;
        return FALSE;
    }
    if (!UpcRayTrace_trace_rev(lobj, sipos, skdir, iwav, ray))
        return FALSE;
    B = ray->s[0].pos[X] - E;
    D = ray->s[0].pos[Y] - F;
    // 繰り返し計算により光線をサーチする.
    x0 = ox - E;
    y0 = oy - F;
    for (i = 0; i < ITER_MAX; i++) {
        det = A * D - B * C;
        if (!det) {
            ray->status = UpcRayE_IterationError;
            return FALSE;
        }
        idet = 1.0 / det;
        sipos[X] = idet * ( D * x0 - B * y0);
        sipos[Y] = idet * (-C * x0 + A * y0);
        if (!(si->shape->getZ(si->shape, sipos[X], sipos[Y], sipos + Z))) {
            ray->status = ray->s[isur_img].status = UpcRayE_NoIntersectionError;
            return FALSE;
        }
        if (!UpcRayTrace_trace_rev(lobj, sipos, skdir, iwav, ray))
            return FALSE;
        dx = ray->s[0].pos[X] - ox;
        dy = ray->s[0].pos[Y] - oy;
        if (fabs(dx) < tol && fabs(dy) < tol)
            return TRUE;
        sumsq_t = sipos[X] * sipos[X] + sipos[Y] * sipos[Y];
        if (!sumsq_t) {
            ray->status = UpcRayE_IterationError;
            return FALSE;
        }
        sumsq_t = 1.0 / sumsq_t;
        A += sipos[X] * dx * sumsq_t;
        B += sipos[Y] * dx * sumsq_t;
        C += sipos[X] * dy * sumsq_t;
        D += sipos[Y] * dy * sumsq_t;
    }
    ray->status = UpcRayE_IterationError;
    return FALSE;
}

/* 出射光線の射出瞳球面上での正弦座標を求める.
 * [入力]
 *     ppos[3] : 主光線の像面上座標
 *     pdir[3] : 主光線の最終面出射方向余弦
 *     rpos[3] : 光線の像面上座標
 *     rdir[3] : 光線の最終面出射方向余弦
 *     expz : 射出瞳z位置(像面基準)  ※通常は負の数
 * [出力]
 *     sinTx, sinTy : 瞳座標(sinθx, sinθy)
 *     戻り値 : 計算できればTRUE, 失敗したらFALSE.
 */
static Bool UpcRayTrace_expCoord(const double *ppos, const double *pdir, const double *rpos, const double *rdir, double expz, double *sinTx, double *sinTy)
{
    double p0[3], p1[3], p2[3];

    if (!expz)
        return FALSE;
    p0[X] = rpos[X] - ppos[X];
    p0[Y] = rpos[Y] - ppos[Y];
    p0[Z] = rpos[Z] - ppos[Z] - expz;
    if (!UpcShapeSph_conicIntersection(-1.0 / expz, 0.0, p0, rdir, p1, p2))
        return FALSE;
    *sinTx = (p1[X] / -expz) + pdir[X];
    *sinTy = (p1[Y] / -expz) + pdir[Y];
    return TRUE;
}

/* 有効径端を通るマージナル光線を探す.(有限物体距離版)
 * [入力]
 *     lobj : UpcLensオブジェクト
 *     pray : 主光線(UpcRayオブジェクト)(正常にトレースできたものを与えること)
 *     dx, dy : = {(0.0, 1.0): YUpper, (0.0, -1.0): YLower, (1.0, 0.0): XUpper, (-1.0, 0.0): XLower}
 *     iwav : 波長番号
 *     tol : 収束精度(同一点とみなす距離 @ 絞り面)
 * [出力]
 *     ray : UpcRayオブジェクト (ray->nsur == lobj.nsur でなければならない)
 *     戻り値 : 光線が見つかればTRUE, 失敗したらFALSE.
 */
Bool UpcRayTrace_msearch0_ea(const UpcLens *lobj, const UpcRay *pray, double dx, double dy, int iwav, double tol, UpcRay *ray)
{
    const double *s0pos = pray->s[0].pos;
    const double *s0dir = pray->s[0].dir;
    int i, sto;
    double ex, ey;
    double p_upp, p_low, p_now, stoUX, stoUY, stoLX, stoLY;

    if (pray->status != UpcRayE_Succeeded) {
        ray->status = UpcRayE_SettingError; // 初期光線が正常でない.
        return FALSE;
    }
    if (lobj->pray == UpcPrayMode_itel) {
        sto = lobj->nsur;
    }
    else if (lobj->pray == UpcPrayMode_otel) {
        sto = 1;
    }
    else {
        sto = lobj->isur_sto;
    }
    p_low = 0.0;
    stoLX = pray->s[sto].pos[X];
    stoLY = pray->s[sto].pos[Y];
    p_now = 0.0001;
    // 囲い込み
    for (i = 0; i < ITER_MAX; i++) {
        ex = p_now * dx;
        ey = p_now * dy;
        if (UpcRayTrace_trace0p(lobj, s0pos, s0dir, ex, ey, iwav, ray)) {
            if (ray->status == UpcRayE_Succeeded) {
                p_low = p_now;
                stoLX = ray->s[sto].pos[X];
                stoLY = ray->s[sto].pos[Y];
                p_now *= 2.0;
            }
            else {
                p_upp = p_now;
                stoUX = ray->s[sto].pos[X];
                stoUY = ray->s[sto].pos[Y];
                break;
            }
        }
        else {
            p_upp = p_now;
            stoUX = HugeVal;
            stoUY = HugeVal;
            break;
        }
    }
    if (i == ITER_MAX) {
        ray->status = UpcRayE_IterationError;
        return FALSE;
    }
    // 収束
    for (i = 0; i < ITER_MAX; i++) {
        p_now = 0.5 * (p_low + p_upp);
        ex = p_now * dx;
        ey = p_now * dy;
        if (UpcRayTrace_trace0p(lobj, s0pos, s0dir, ex, ey, iwav, ray)) {
            if (ray->status == UpcRayE_Succeeded) {
                p_low = p_now;
                stoLX = ray->s[sto].pos[X];
                stoLY = ray->s[sto].pos[Y];
            }
            else {
                p_upp = p_now;
                stoUX = ray->s[sto].pos[X];
                stoUY = ray->s[sto].pos[Y];
            }
            if (fabs(stoUX - stoLX) < tol && fabs(stoUY - stoLY) < tol) {
                if (ray->status != UpcRayE_Succeeded)
                    UpcRayTrace_trace0p(lobj, s0pos, s0dir, p_low * dx, p_low * dy, iwav, ray);
                return TRUE;
            }
        }
        else {
            p_upp = p_now;
            stoUX = HugeVal;
            stoUY = HugeVal;
        }
    }
    ray->status = UpcRayE_IterationError;
    return FALSE;
}

/* 有効径端を通るマージナル光線を探す.(無限遠物体版)
 * [入力]
 *     lobj : Lensオブジェクト
 *     pray : 主光線(UpcRayオブジェクト)(正常にトレースできたものを与えること)
 *     dx, dy : = {(0.0, 1.0): YUpper, (0.0, -1.0): YLower, (1.0, 0.0): XUpper, (-1.0, 0.0): XLower}
 *     iwav : 波長番号
 *     tol : 収束精度(同一点とみなす距離 @ 絞り面)
 * [出力]
 *     ray : UpcRayオブジェクト (ray->nsur == lobj.nsur でなければならない)
 *     戻り値 : 光線が見つかればTRUE, 失敗したらFALSE.
 */
Bool UpcRayTrace_msearch1_ea(const UpcLens *lobj, const UpcRay *pray, double dx, double dy, int iwav, double tol, UpcRay *ray)
{
    const double *s0dir = pray->s[0].dir;
    int i, sto;
    double s1x_ini, s1y_ini, tpos[3];
    double p_upp, p_low, p_now, stoUX, stoUY, stoLX, stoLY;

    if (pray->status != UpcRayE_Succeeded) {
        ray->status = UpcRayE_SettingError; // 初期光線が正常でない.
        return FALSE;
    }
    if (!s0dir[Z]) {
        ray->status = UpcRayE_SettingError; // 初期光線が正常でない.
        return FALSE;
    }
    if (lobj->pray == UpcPrayMode_itel) {
        sto = lobj->nsur;
    }
    else if (lobj->pray == UpcPrayMode_otel) {
        sto = 1;
    }
    else {
        sto = lobj->isur_sto;
    }
    s1x_ini = pray->s[1].pos[X] - pray->s[1].pos[Z] * s0dir[X] / s0dir[Z];
    s1y_ini = pray->s[1].pos[Y] - pray->s[1].pos[Z] * s0dir[Y] / s0dir[Z];
    p_low = 0.0;
    stoLX = pray->s[sto].pos[X];
    stoLY = pray->s[sto].pos[Y];
    p_now = 0.01;
    tpos[Z] = 0.0;
    // 囲い込み
    for (i = 0; i < ITER_MAX; i++) {
        tpos[X] = s1x_ini + p_now * dx;
        tpos[Y] = s1y_ini + p_now * dy;
        if (UpcRayTrace_trace1(lobj, s0dir, tpos, iwav, 0, ray)) {
            if (ray->status == UpcRayE_Succeeded) {
                p_low = p_now;
                stoLX = ray->s[sto].pos[X];
                stoLY = ray->s[sto].pos[Y];
                p_now *= 2.0;
            }
            else {
                p_upp = p_now;
                stoUX = ray->s[sto].pos[X];
                stoUY = ray->s[sto].pos[Y];
                break;
            }
        }
        else {
            p_upp = p_now;
            stoUX = HugeVal;
            stoUY = HugeVal;
            break;
        }
    }
    if (i == ITER_MAX) {
        ray->status = UpcRayE_IterationError;
        return FALSE;
    }
    // 収束
    for (i = 0; i < ITER_MAX; i++) {
        p_now = 0.5 * (p_low + p_upp);
        tpos[X] = s1x_ini + p_now * dx;
        tpos[Y] = s1y_ini + p_now * dy;
        if (UpcRayTrace_trace1(lobj, s0dir, tpos, iwav, 0, ray)) {
            if (ray->status == UpcRayE_Succeeded) {
                p_low = p_now;
                stoLX = ray->s[sto].pos[X];
                stoLY = ray->s[sto].pos[Y];
            }
            else {
                p_upp = p_now;
                stoUX = ray->s[sto].pos[X];
                stoUY = ray->s[sto].pos[Y];
            }
            if (fabs(stoUX - stoLX) < tol && fabs(stoUY - stoLY) < tol) {
                if (ray->status != UpcRayE_Succeeded) {
                    tpos[X] = s1x_ini + p_low * dx;
                    tpos[Y] = s1y_ini + p_low * dy;
                    UpcRayTrace_trace1(lobj, s0dir, tpos, iwav, 0, ray);
                }
                return TRUE;
            }
        }
        else {
            p_upp = p_now;
            stoUX = HugeVal;
            stoUY = HugeVal;
        }
    }
    ray->status = UpcRayE_IterationError;
    return FALSE;
}

/* 像面において主光線との角度が指定sinθとなるマージナル光線を探す.(有限物体距離版)
 * [入力]
 *     lobj : UpcLensオブジェクト
 *     pray : 主光線(UpcRayオブジェクト)(正常にトレースできたものを与えること)
 *     sin_theta : sinθ
 *     dx, dy : = {(0.0, 1.0): YUpper, (0.0, -1.0): YLower, (1.0, 0.0): XUpper, (-1.0, 0.0): XLower}
 *     expz : 射出瞳z位置(像面基準)  ※通常は負の数
 *     iwav : 波長番号
 *     tol : 収束精度(同一点とみなす距離 @ 絞り面)
 * [出力]
 *     ray : UpcRayオブジェクト (ray->nsur == lobj.nsur でなければならない)
 *     戻り値 : 光線が見つかればTRUE, 失敗したらFALSE.
 */
Bool UpcRayTrace_msearch0_na(const UpcLens *lobj, const UpcRay *pray, double sin_theta, double dx, double dy, double expz, int iwav, double tol, UpcRay *ray)
{
    const int nsur = lobj->nsur;
    const double *s0pos = pray->s[0].pos;
    const double *s0dir = pray->s[0].dir;
    const double *ppos = pray->s[nsur + 1].pos;
    const double *pdir = pray->s[nsur].dir;
    int i, sto;
    double ex, ey;
    double p_low, p_upp, p_now, stoUX, stoUY, stoLX, stoLY;
    double sinTx, sinTy, sinT;
    
    if (sin_theta >= 1.0) {
        ray->status = UpcRayE_SettingError; // 初期光線が正常でない.
        return FALSE;
    }
    if (lobj->pray == UpcPrayMode_itel) {
        sto = lobj->nsur;
    }
    else if (lobj->pray == UpcPrayMode_otel) {
        sto = 1;
    }
    else {
        sto = lobj->isur_sto;
    }
    p_low = 0.00;
    stoLX = pray->s[sto].pos[X];
    stoLY = pray->s[sto].pos[Y];
    p_now = 0.0001;
    // 囲い込み
    for (i = 0; i < ITER_MAX; i++) {
        ex = p_now * dx;
        ey = p_now * dy;
        if (UpcRayTrace_trace0p(lobj, s0pos, s0dir, ex, ey, iwav, ray)) {
            if (!UpcRayTrace_expCoord(ppos, pdir, ray->s[nsur + 1].pos, ray->s[nsur].dir, expz, &sinTx, &sinTy)) {
                p_upp = p_now;
                stoUX = HugeVal;
                stoUY = HugeVal;
                break;
            }
            // sinT = sqrt(sinTx * sinTx + sinTy * sinTy);
            sinT = upc_hypot(sinTx, sinTy);
            if (sinT > sin_theta) {
                p_upp = p_now;
                stoUX = ray->s[sto].pos[X];
                stoUY = ray->s[sto].pos[Y];
                break;
            }
            else {
                p_low = p_now;
                stoLX = ray->s[sto].pos[X];
                stoLY = ray->s[sto].pos[Y];
                p_now *= 2.0;
            }
        }
        else {
            p_upp = p_now;
            stoUX = HugeVal;
            stoUY = HugeVal;
            break;
        }
    }
    if (i == ITER_MAX) {
        ray->status = UpcRayE_IterationError;
        return FALSE;
    }
    // 収束
    for (i = 0; i < ITER_MAX; i++) {
        p_now = 0.5 * (p_low + p_upp);
        ex = p_now * dx;
        ey = p_now * dy;
        if (UpcRayTrace_trace0p(lobj, s0pos, s0dir, ex, ey, iwav, ray)) {
            if (!UpcRayTrace_expCoord(ppos, pdir, ray->s[nsur + 1].pos, ray->s[nsur].dir, expz, &sinTx, &sinTy)) {
                p_upp = p_now;
                stoUX = HugeVal;
                stoUY = HugeVal;
                continue;
            }
            // sinT = sqrt(sinTx * sinTx + sinTy * sinTy);
            sinT = upc_hypot(sinTx, sinTy);
            if (sinT <= sin_theta) {
                p_low = p_now;
                stoLX = ray->s[sto].pos[X];
                stoLY = ray->s[sto].pos[Y];
            }
            else {
                p_upp = p_now;
                stoUX = ray->s[sto].pos[X];
                stoUY = ray->s[sto].pos[Y];
            }
            if (fabs(stoUX - stoLX) < tol && fabs(stoUY - stoLY) < tol) {
                return TRUE;
            }
        }
        else {
            p_upp = p_now;
            stoUX = HugeVal;
            stoUY = HugeVal;
        }
    }
    ray->status = UpcRayE_IterationError;
    return FALSE;
}

/* 像面において主光線との角度が指定sinθとなるマージナル光線を探す.(無限遠物体版)
 * [入力]
 *     lobj : UpcLensオブジェクト
 *     pray : 主光線(UpcRayオブジェクト)(正常にトレースできたものを与えること)
 *     sin_theta : sinθ
 *     dx, dy : = {(0.0, 1.0): YUpper, (0.0, -1.0): YLower, (1.0, 0.0): XUpper, (-1.0, 0.0): XLower}
 *     expz : 射出瞳z位置(像面基準)  ※通常は負の数
 *     iwav : 波長番号
 *     tol : 収束精度(同一点とみなす距離 @ 絞り面)
 * [出力]
 *     ray : UpcRayオブジェクト (ray->nsur == lobj.nsur でなければならない)
 *     戻り値 : 光線が見つかればTRUE, 失敗したらFALSE.
 */
Bool UpcRayTrace_msearch1_na(const UpcLens *lobj, const UpcRay *pray, double sin_theta, double dx, double dy, double expz, int iwav, double tol, UpcRay *ray)
{
    const int nsur = lobj->nsur;
    const double *s0dir = pray->s[0].dir;
    const double *ppos = pray->s[nsur + 1].pos;
    const double *pdir = pray->s[nsur].dir;
    int i, sto;
    double tpos[3], s1x_ini, s1y_ini;
    double p_low, p_upp, p_now, stoUX, stoUY, stoLX, stoLY;
    double sinTx, sinTy, sinT;
    
    if (sin_theta >= 1.0) {
        ray->status = UpcRayE_SettingError; // 初期光線が正常でない.
        return FALSE;
    }
    if (!s0dir[Z]) {
        ray->status = UpcRayE_SettingError; // 初期光線が正常でない.
        return FALSE;
    }
    if (lobj->pray == UpcPrayMode_itel) {
        sto = lobj->nsur;
    }
    else if (lobj->pray == UpcPrayMode_otel) {
        sto = 1;
    }
    else {
        sto = lobj->isur_sto;
    }
    s1x_ini = pray->s[1].pos[X] - pray->s[1].pos[Z] * s0dir[X] / s0dir[Z];
    s1y_ini = pray->s[1].pos[Y] - pray->s[1].pos[Z] * s0dir[Y] / s0dir[Z];
    p_low = 0.00;
    stoLX = pray->s[sto].pos[X];
    stoLY = pray->s[sto].pos[Y];
    p_now = 0.01;
    tpos[Z] = 0.0;
    // 囲い込み
    for (i = 0; i < ITER_MAX; i++) {
        tpos[X] = s1x_ini + p_now * dx;
        tpos[Y] = s1y_ini + p_now * dy;
        if (UpcRayTrace_trace1(lobj, s0dir, tpos, iwav, 0, ray)) {
            if (!UpcRayTrace_expCoord(ppos, pdir, ray->s[nsur + 1].pos, ray->s[nsur].dir, expz, &sinTx, &sinTy)) {
                p_upp = p_now;
                stoUX = HugeVal;
                stoUY = HugeVal;
                break;
            }
            // sinT = sqrt(sinTx * sinTx + sinTy * sinTy);
            sinT = upc_hypot(sinTx, sinTy);
            if (sinT > sin_theta) {
                p_upp = p_now;
                stoUX = ray->s[sto].pos[X];
                stoUY = ray->s[sto].pos[Y];
                break;
            }
            else {
                p_low = p_now;
                stoLX = ray->s[sto].pos[X];
                stoLY = ray->s[sto].pos[Y];
                p_now *= 2.0;
            }
        }
        else {
            p_upp = p_now;
            stoUX = HugeVal;
            stoUY = HugeVal;
            break;
        }
    }
    if (i == ITER_MAX) {
        ray->status = UpcRayE_IterationError;
        return FALSE;
    }
    // 収束
    for (i = 0; i < ITER_MAX; i++) {
        p_now = 0.5 * (p_low + p_upp);
        tpos[X] = s1x_ini + p_now * dx;
        tpos[Y] = s1y_ini + p_now * dy;
        if (UpcRayTrace_trace1(lobj, s0dir, tpos, iwav, 0, ray)) {
            if (!UpcRayTrace_expCoord(ppos, pdir, ray->s[nsur + 1].pos, ray->s[nsur].dir, expz, &sinTx, &sinTy)) {
                p_upp = p_now;
                stoUX = HugeVal;
                stoUY = HugeVal;
                continue;
            }
            // sinT = sqrt(sinTx * sinTx + sinTy * sinTy);
            sinT = upc_hypot(sinTx, sinTy);
            if (sinT <= sin_theta) {
                p_low = p_now;
                stoLX = ray->s[sto].pos[X];
                stoLY = ray->s[sto].pos[Y];
            }
            else {
                p_upp = p_now;
                stoUX = ray->s[sto].pos[X];
                stoUY = ray->s[sto].pos[Y];
            }
            if (fabs(stoUX - stoLX) < tol && fabs(stoUY - stoLY) < tol) {
                return TRUE;
            }
        }
        else {
            p_upp = p_now;
            stoUX = HugeVal;
            stoUY = HugeVal;
        }
    }
    ray->status = UpcRayE_IterationError;
    return FALSE;
}

/* 物体面において主光線との角度が指定sinθとなるマージナル光線を探す.(有限物体距離版)
 * [入力]
 *     lobj : UpcLensオブジェクト
 *     pray : 主光線(UpcRayオブジェクト)(正常にトレースできたものを与えること)
 *     sin_theta : sinθ
 *     dx, dy : = {(0.0, 1.0): YUpper, (0.0, -1.0): YLower, (1.0, 0.0): XUpper, (-1.0, 0.0): XLower}
 *     iwav : 波長番号
 * [出力]
 *     ray : UpcRayオブジェクト (ray->nsur == lobj.nsur でなければならない)
 *     戻り値 : 光線が見つかればTRUE, 失敗したらFALSE.
 */
Bool UpcRayTrace_msearch0_nao(const UpcLens *lobj, const UpcRay *pray, double sin_theta, double dx, double dy, int iwav, UpcRay *ray)
{
    double ex, ey;
    double len2 = dx * dx + dy * dy;
    
    if (len2) {
        double ilen = 1.0 / sqrt(len2);

        dx *= ilen;
        dy *= ilen;
        ex = sin_theta * dx;
        ey = sin_theta * dy;
        return UpcRayTrace_trace0p(lobj, pray->s[0].pos, pray->s[0].dir, ex, ey, iwav, ray);
    }
    ray->status = UpcRayE_SettingError;
    return FALSE;
}

/* 計算済みの光路に沿って波面収差を計算する.
 * [入力]
 *     lobj     : UpcLensオブジェクト
 *     ray      : UpcRayオブジェクト
 *     ensPos[3]: 入射側参照球面の表面上にある任意点の座標 (第1面座標系)
 *     ensNV[3] : ensPosにおける入射側参照球面の単位法線ベクトル (ensNV[Z] > 0)
 *     ensCurv  : 入射側参照球面の曲率 (符号も考慮. =0で平面)
 *     exsPos[3]: 出射側参照球面の表面上にある任意点の座標 (最終面座標系)
 *     exsNV[3] : exsPosにおける出射側参照球面の単位法線ベクトル (exsNV[Z] > 0)
 *     exsCurv  : 出射側参照球面の曲率 (符号も考慮. =0で平面)
 *     iwav     : 波長番号
 * [出力]
 *     wab : 波面収差(波長単位)
 *     戻り値 : 正常に計算できればTRUE, 失敗したらFALSE.
 */
Bool UpcRayTrace_calWab(const UpcLens *lobj, const UpcRay *ray, const double *ensPos, const double *ensNV, double ensCurv, const double *exsPos, const double *exsNV, double exsCurv, int iwav, double *wab)
{
    const int nsur = lobj->nsur;
    int i;
    double tpos[3], tdir[3], t;
    double w1, w, wk;
    UpcSurf **sur = lobj->s;
    
    // 入射瞳から第1面までの計算
    UpcVector_SUB(tpos, ray->s[1].pos, ensPos);
    if (sur[1]->dec->decenterFlg)
        UpcDecenter_rotForward(sur[1]->dec, ray->s[0].dir, tdir);
    else
        UpcVector_COPY(tdir, ray->s[0].dir);
    if (ensCurv) {
        double b = UpcVector_DOT(tdir, tpos) - UpcVector_DOT(tdir, ensNV) / ensCurv;
        double c = UpcVector_DOT(tpos, tpos) - 2.0 * UpcVector_DOT(tpos, ensNV) / ensCurv;
        double k1, k2, z1, z2;

        if (UpcShapeSph_solveQuadEq(1.0, 2.0 * b, c, &k1, &k2) == 2) {
            z1 = tpos[Z] + k1 * tdir[Z];
            z2 = tpos[Z] + k2 * tdir[Z];
            t = (fabs(z1) < fabs(z2)) ? k1 : k2;
        }
        else
            return FALSE; // no intersection
    }
    else {
        t = - UpcVector_DOT(ensNV, tpos) / UpcVector_DOT(ensNV, tdir);
    }
    w1 = - fabs(sur[0]->_n[iwav]) * t;
    // 第1面から最終面までの計算
    w = 0.0;
    for (i = 1; i < nsur; i++)
        w += ray->s[i].opl * fabs(sur[i]->_n[iwav]);
    // 最終面から射出瞳までの計算
    UpcVector_SUB(tpos, ray->s[nsur].pos, exsPos);
    UpcVector_COPY(tdir, ray->s[nsur].dir);
    if (exsCurv) {
        double b = UpcVector_DOT(tdir, tpos) - UpcVector_DOT(tdir, exsNV) / exsCurv;
        double c = UpcVector_DOT(tpos, tpos) - 2.0 * UpcVector_DOT(tpos, exsNV) / exsCurv;
        double k1, k2, z1, z2;

        if (UpcShapeSph_solveQuadEq(1.0, 2.0 * b, c, &k1, &k2) == 2) {
            z1 = tpos[Z] + k1 * tdir[Z];
            z2 = tpos[Z] + k2 * tdir[Z];
            t = (fabs(z1) < fabs(z2)) ? k1 : k2;
        }
        else
            return FALSE; // no intersection
    }
    else {
        t = - UpcVector_DOT(exsNV, tpos) / UpcVector_DOT(exsNV, tdir);
    }
    wk = fabs(sur[nsur]->_n[iwav]) * t;
    // 算出
    *wab = (w1 + w + wk) * 1.0e6 / UpcLens_WL(lobj, iwav);
    return TRUE;
}

/* 計算済みの光路に沿って非点収差を追跡する.
 * [入力]
 *     lobj  : UpcLensオブジェクト
 *     ray   : UpcRayオブジェクト
 *     h1, u0: 追跡初期値
 *     iwav  : 波長番号
 * [出力]
 *     hm, um, hs, us : 追跡値
 *     戻り値 : 正常に計算できればTRUE, 失敗したらFALSE.
 */
Bool UpcRayTrace_calAstigmatism(const UpcLens *lobj, const UpcRay *ray, double h1, double u0, int iwav, double *h_m, double *u_m, double *h_s, double *u_s)
{
    const int nsur = lobj->nsur;
    UpcSurf **sur = lobj->s;
    double um = u0;
    double us = u0;
    double hm = h1 / ray->s[1].cosI;
    double hs = h1;
    double cm, cs, nd, gt;
    UpcSurf *sj, *sk;
    const UpcRayComponent *rj, *rk;
    int i;
    
    if (!(sur[1]->shape->getCurv(sur[1]->shape, ray->s[1].pos[Y], &cm, &cs)))
        return FALSE;
    nd = fabs(sur[0]->_n[iwav] / sur[1]->_n[iwav]);
    gt = ray->s[1].cosO - nd * ray->s[1].cosI;
    um = (nd * ray->s[1].cosI * um + hm * gt * cm) / ray->s[1].cosO;
    us = nd * us + hs * gt * cs;
    sj = sur[1];
    rj = ray->s + 1;
    for (i = 2; i <= nsur; i++) {
        sk = sur[i];
        rk = ray->s + i;
        if (!(sk->shape->getCurv(sk->shape, rk->pos[Y], &cm, &cs)))
            return FALSE;
        nd = fabs(sj->_n[iwav] / sk->_n[iwav]);
        gt = rk->cosO - nd * rk->cosI;
        hm = (hm * rj->cosO - rj->opl * um) / rk->cosI;
        um = (nd * rk->cosI * um + hm * gt * cm) / rk->cosO;
        hs = hs - rj->opl * us;
        us = nd * us + hs * gt * cs;
        sj = sk;
        rj = rk;
    }
    *h_m = hm;
    *u_m = um;
    *h_s = hs;
    *u_s = us;
    return TRUE;
}

/*------------------------------------------------------------------------*/
// UpcTraceFactoryオブジェクト

static void UpcTraceFactory_psearch(UpcTraceFactory *self, const UpcLens *lobj);
static void UpcTraceFactory_msearch(UpcTraceFactory *self, const UpcLens *lobj);
static void UpcTraceFactory_setPInfo_infObjD(UpcTraceFactory *self, const UpcLens *lobj);
static Bool UpcTraceFactory_traceP_infObjD(const UpcTraceFactory *self, const UpcLens *lobj, int iwav, int ihgt, double px, double py, UpcRay *ray);
static void UpcTraceFactory_setPInfo_finiteObjD(UpcTraceFactory *self);
static Bool UpcTraceFactory_traceP_finiteObjD(const UpcTraceFactory *self, const UpcLens *lobj, int iwav, int ihgt, double px, double py, UpcRay *ray);
static Bool UpcTraceFactory_setWInfo(UpcTraceFactory *self, const UpcLens *lobj);

static void UpcTraceFactory_dealloc(UpcTraceFactory *self)
{
    int i, n = self->nwav * self->nhgt * self->nray;
    UpcRay **t = self->refRay[0][0];
    
    for (i = 0; i < n; i++)
        UpcObj_RELEASE(t[i]);
    free(self->refRay[0][0]);
    free(self->refRay[0]);
    free(self->refRay);
    upc_freeD3(self->pInfo);
    upc_freeD3(self->wInfo);
    upc_freeB2(self->refRayStatus);
    upc_freeB2(self->readyWabc);
    free(self);
}

/* (注) lset()が通ったレンズオブジェクトを与えること.
 */
UpcTraceFactory *UpcTraceFactory_init(const UpcLens *lobj)
{
    UpcTraceFactory *self;
    UpcRay ***refRay0, **refRay00;
    
    self = (UpcTraceFactory *)malloc(sizeof(UpcTraceFactory));
    if (self) {
        const int nsur = lobj->nsur;
        const int nwav = lobj->nwav;
        const int nhgt = lobj->nhgt;
        const int nray = 5;
        const int npInfo = 15;
        const int nwInfo = 15;
        const int n = nwav * nhgt * nray;
        int iwav, ihgt, i, j;
        
        self->_refCount = 1;
        self->_dealloc = UpcTraceFactory_dealloc;
        self->refRay = (UpcRay ****)malloc(sizeof(UpcRay ***) * nwav);
        refRay0 = (UpcRay ***)malloc(sizeof(UpcRay **) * nwav * nhgt);
        refRay00 = (UpcRay **)malloc(sizeof(UpcRay *) * n);
        self->pInfo = upc_mallocD3(nwav, nhgt, npInfo);
        self->wInfo = upc_mallocD3(nwav, nhgt, nwInfo);
        self->refRayStatus = upc_mallocB2(nwav, nhgt);
        self->readyWabc = upc_mallocB2(nwav, nhgt);
        if (!self->refRay || !refRay0 || !refRay00 || !self->pInfo || !self->wInfo || !self->refRayStatus || !self->readyWabc) {
            goto fail_return;
        }
        for (i = 0; i < n; i++) {
            refRay00[i] = UpcRay_init(nsur);
            if (!refRay00[i]) {
                for (j = 0; j < i; j++)
                    refRay00[j]->_dealloc(refRay00[j]);
                goto fail_return;
            }
        }
        for (iwav = 0; iwav < nwav; iwav++) {
            self->refRay[iwav] = refRay0 + iwav * nhgt;
            for (ihgt = 0; ihgt < nhgt; ihgt++) {
                self->refRay[iwav][ihgt] = refRay00 + iwav * (nhgt * nray) + ihgt * nray;
            }
        }
        self->nwav = nwav;
        self->nhgt = nhgt;
        self->nray = nray;
        
        UpcTraceFactory_psearch(self, lobj);
        UpcTraceFactory_msearch(self, lobj);
        if (lobj->infObjD) {
            UpcTraceFactory_setPInfo_infObjD(self, lobj);
            self->traceP = UpcTraceFactory_traceP_infObjD;
        }
        else {
            UpcTraceFactory_setPInfo_finiteObjD(self);
            self->traceP = UpcTraceFactory_traceP_finiteObjD;
        }
        UpcTraceFactory_setWInfo(self, lobj);
    }
    return self;
    
fail_return:
    free(refRay00);
    free(refRay0);
    free(self->refRay);
    upc_freeD3(self->pInfo);
    upc_freeD3(self->wInfo);
    upc_freeB2(self->refRayStatus);
    upc_freeB2(self->readyWabc);
    free(self);
    return NULL;
}

/* 与えられたレンズの主光線をサーチする.
 */
static void UpcTraceFactory_psearch(UpcTraceFactory *self, const UpcLens *lobj)
{
    int i, iwav, ihgt;
    const int n = self->nwav * self->nhgt;
    const int nhgt = self->nhgt;
    const double tol = lobj->tol_pray;
    UpcRay *pray;

    for (i = 0; i < n; i++) {
        ihgt = i % nhgt;
        iwav = i / nhgt;
        pray = self->refRay[iwav][ihgt][PRAY];
        if (lobj->infObjD) {
            if (lobj->pray == UpcPrayMode_stop) {
                double x = tan(RADIANS(UpcLens_HGTX(lobj, ihgt)));
                double y = tan(RADIANS(UpcLens_HGTY(lobj, ihgt)));
                double ilen = 1.0 / sqrt(x * x + y * y + 1.0);
                
                x *= ilen;
                y *= ilen;
                UpcRayTrace_psearch1(lobj, x, y, lobj->isur_sto, 0.0, 0.0, iwav, tol, pray);
            }
            else {
                pray->status = UpcRayE_SettingError; // 無限遠物体のときはテレセントリック指定できない
            }
        }
        else {
            double x = UpcLens_HGTX(lobj, ihgt);
            double y = UpcLens_HGTY(lobj, ihgt);

            if (lobj->pray == UpcPrayMode_stop) {
                UpcRayTrace_psearch0(lobj, x, y, lobj->isur_sto, 0.0, 0.0, iwav, tol, pray);
            }
            else if (lobj->pray == UpcPrayMode_itel) {
                UpcRayTrace_psearch0_xt(lobj, x, y, 0.0, 0.0, iwav, tol, pray);
            }
            else if (lobj->pray == UpcPrayMode_otel) {
                double z, s0pos[3], s0dir[] = {0.0, 0.0, 1.0};
                
                lobj->s[0]->shape->getZ(lobj->s[0]->shape, x, y, &z);
                UpcVector_SET(s0pos, x, y, z);
                UpcRayTrace_trace0d(lobj, s0pos, s0dir, iwav, 0, pray);
            }
        }
    }
}

/* 与えられたレンズの参照光線（マージナル光線）をサーチする.
 * (注) 主光線をサーチしてから実行すること.
 */
static void UpcTraceFactory_msearch(UpcTraceFactory *self, const UpcLens *lobj)
{
    int i, iwav, ihgt;
    const int n = self->nwav * self->nhgt;
    const int nhgt = self->nhgt;
    const double tol = lobj->tol_mray;
    double expz = lobj->parax[lobj->iwav_pri].exp - lobj->s[lobj->nsur]->d;
    double ape = lobj->aperture;
    Bool xSymmetric = lobj->xSymmetric;
    Bool ySymmetric = lobj->ySymmetric;
    UpcRay **refray;
    double hgt[2], vig[5];

    for (i = 0; i < n; i++) {
        ihgt = i % nhgt;
        iwav = i / nhgt;
        refray = self->refRay[iwav][ihgt];
        hgt[X] = UpcLens_HGTX(lobj, ihgt);
        hgt[Y] = UpcLens_HGTY(lobj, ihgt);
        vig[YUPP] = UpcLens_VIG(lobj, ihgt, YUPP);
        vig[YLOW] = UpcLens_VIG(lobj, ihgt, YLOW);
        vig[XUPP] = UpcLens_VIG(lobj, ihgt, XUPP);
        vig[XLOW] = UpcLens_VIG(lobj, ihgt, XLOW);

        if (lobj->mray == UpcMrayMode_eape) {
            if (refray[PRAY]->status != UpcRayE_Succeeded) {
                refray[YUPP]->status = UpcRayE_UnCalculated;
                refray[YLOW]->status = UpcRayE_UnCalculated;
                refray[XUPP]->status = UpcRayE_UnCalculated;
                refray[XLOW]->status = UpcRayE_UnCalculated;
                self->refRayStatus[iwav][ihgt] = FALSE;
                continue;
            }
            if (lobj->infObjD) {
                UpcRayTrace_msearch1_ea(lobj, refray[PRAY], 0.0, 1.0, iwav, tol, refray[YUPP]);
                if (xSymmetric && !hgt[Y])
                    UpcRay_copy_XZmirror(refray[YLOW], refray[YUPP]);
                else
                    UpcRayTrace_msearch1_ea(lobj, refray[PRAY], 0.0, -1.0, iwav, tol, refray[YLOW]);

                UpcRayTrace_msearch1_ea(lobj, refray[PRAY], 1.0, 0.0, iwav, tol, refray[XUPP]);
                if (ySymmetric && !hgt[X])
                    UpcRay_copy_YZmirror(refray[XLOW], refray[XUPP]);
                else
                    UpcRayTrace_msearch1_ea(lobj, refray[PRAY], -1.0, 0.0, iwav, tol, refray[XLOW]);
            }
            else {
                UpcRayTrace_msearch0_ea(lobj, refray[PRAY], 0.0, 1.0, iwav, tol, refray[YUPP]);
                if (xSymmetric && !hgt[Y])
                    UpcRay_copy_XZmirror(refray[YLOW], refray[YUPP]);
                else
                    UpcRayTrace_msearch0_ea(lobj, refray[PRAY], 0.0, -1.0, iwav, tol, refray[YLOW]);

                UpcRayTrace_msearch0_ea(lobj, refray[PRAY], 1.0, 0.0, iwav, tol, refray[XUPP]);
                if (ySymmetric && !hgt[X])
                    UpcRay_copy_YZmirror(refray[XLOW], refray[XUPP]);
                else
                    UpcRayTrace_msearch0_ea(lobj, refray[PRAY], -1.0, 0.0, iwav, tol, refray[XLOW]);
            }
        }
        else if (lobj->mray == UpcMrayMode_isin) {
            if (refray[PRAY]->status < UpcRayE_OutOfArea) {
                refray[YUPP]->status = UpcRayE_UnCalculated;
                refray[YLOW]->status = UpcRayE_UnCalculated;
                refray[XUPP]->status = UpcRayE_UnCalculated;
                refray[XLOW]->status = UpcRayE_UnCalculated;
                self->refRayStatus[iwav][ihgt] = FALSE;
                continue;
            }
            if (lobj->infObjD) {
                UpcRayTrace_msearch1_na(lobj, refray[PRAY], vig[YUPP] * ape, 0.0, 1.0, expz, iwav, tol, refray[YUPP]);
                if (xSymmetric && !hgt[Y] && vig[YUPP] == vig[YLOW])
                    UpcRay_copy_XZmirror(refray[YLOW], refray[YUPP]);
                else
                    UpcRayTrace_msearch1_na(lobj, refray[PRAY], vig[YLOW] * ape, 0.0, -1.0, expz, iwav, tol, refray[YLOW]);

                UpcRayTrace_msearch1_na(lobj, refray[PRAY], vig[XUPP] * ape, 1.0, 0.0, expz, iwav, tol, refray[XUPP]);
                if (ySymmetric && !hgt[X] && vig[XUPP] == vig[XLOW])
                    UpcRay_copy_YZmirror(refray[XLOW], refray[XUPP]);
                else
                    UpcRayTrace_msearch1_na(lobj, refray[PRAY], vig[XLOW] * ape, -1.0, 0.0, expz, iwav, tol, refray[XLOW]);
            }
            else {
                UpcRayTrace_msearch0_na(lobj, refray[PRAY], vig[YUPP] * ape, 0.0, 1.0, expz, iwav, tol, refray[YUPP]);
                if (xSymmetric && !hgt[Y] && vig[YUPP] == vig[YLOW])
                    UpcRay_copy_XZmirror(refray[YLOW], refray[YUPP]);
                else
                    UpcRayTrace_msearch0_na(lobj, refray[PRAY], vig[YLOW] * ape, 0.0, -1.0, expz, iwav, tol, refray[YLOW]);

                UpcRayTrace_msearch0_na(lobj, refray[PRAY], vig[XUPP] * ape, 1.0, 0.0, expz, iwav, tol, refray[XUPP]);
                if (ySymmetric && !hgt[X] && vig[XUPP] == vig[XLOW])
                    UpcRay_copy_YZmirror(refray[XLOW], refray[XUPP]);
                else
                    UpcRayTrace_msearch0_na(lobj, refray[PRAY], vig[XLOW] * ape, -1.0, 0.0, expz, iwav, tol, refray[XLOW]);
            }
        }
        else if (lobj->mray == UpcMrayMode_osin) {
            if (refray[PRAY]->status < UpcRayE_OutOfArea) {
                refray[YUPP]->status = UpcRayE_UnCalculated;
                refray[YLOW]->status = UpcRayE_UnCalculated;
                refray[XUPP]->status = UpcRayE_UnCalculated;
                refray[XLOW]->status = UpcRayE_UnCalculated;
                self->refRayStatus[iwav][ihgt] = FALSE;
                continue;
            }
            if (lobj->infObjD) {
                // 無限遠物体のときはNAO指定できない
                refray[YUPP]->status = UpcRayE_UnCalculated;
                refray[YLOW]->status = UpcRayE_UnCalculated;
                refray[XUPP]->status = UpcRayE_UnCalculated;
                refray[XLOW]->status = UpcRayE_UnCalculated;
                self->refRayStatus[iwav][ihgt] = FALSE;
                continue;
            }
            else {
                UpcRayTrace_msearch0_nao(lobj, refray[PRAY], vig[YUPP] * ape, 0.0, 1.0, iwav, refray[YUPP]);
                if (xSymmetric && !hgt[Y] && vig[YUPP] == vig[YLOW])
                    UpcRay_copy_XZmirror(refray[YLOW], refray[YUPP]);
                else
                    UpcRayTrace_msearch0_nao(lobj, refray[PRAY], vig[YLOW] * ape, 0.0, -1.0, iwav, refray[YLOW]);

                UpcRayTrace_msearch0_nao(lobj, refray[PRAY], vig[XUPP] * ape, 1.0, 0.0, iwav, refray[XUPP]);
                if (ySymmetric && !hgt[X] && vig[XUPP] == vig[XLOW])
                    UpcRay_copy_YZmirror(refray[XLOW], refray[XUPP]);
                else
                    UpcRayTrace_msearch0_nao(lobj, refray[PRAY], vig[XLOW] * ape, -1.0, 0.0, iwav, refray[XLOW]);
            }
        }
        if (refray[YUPP]->status > UpcRayE_UnCalculated &&
            refray[YLOW]->status > UpcRayE_UnCalculated && 
            refray[XUPP]->status > UpcRayE_UnCalculated && 
            refray[XLOW]->status > UpcRayE_UnCalculated)
            self->refRayStatus[iwav][ihgt] = TRUE;
        else
            self->refRayStatus[iwav][ihgt] = FALSE;
    }
}

/* 入射基準面の設定(無限遠物体用)
 * (注) 主光線, マージナル光線をサーチしてから実行すること.
 */
static void UpcTraceFactory_setPInfo_infObjD(UpcTraceFactory *self, const UpcLens *lobj)
{
    Bool s1dec = lobj->s[1]->dec->decenterFlg;
    int i, iwav, ihgt, iray;
    const int n = lobj->nwav * lobj->nhgt;
    const int nhgt = lobj->nhgt;
    const int nray = 5;
    UpcRay **refray;
    double tpos[3], tdir[3], t[3];
    double x, y, sd, sdmax;
    double s1x[5], s1y[5];
    double *pInfo;

    for (i = 0; i < n; i++) {
        ihgt = i % nhgt;
        iwav = i / nhgt;
        refray = self->refRay[iwav][ihgt];
        sdmax = 0.0;
        pInfo = self->pInfo[iwav][ihgt];
            
        if (!(self->refRayStatus[iwav][ihgt]))
            continue;
        for (iray = 0; iray < nray; iray++) {
            UpcVector_COPY(tpos, refray[iray]->s[1].pos);
            UpcVector_COPY(tdir, refray[iray]->s[0].dir);
            if (s1dec) {
                UpcDecenter_convReverse(lobj->s[1]->dec, tpos, t);
                UpcVector_COPY(tpos, t);
            }
            s1x[iray] = tpos[X] - tpos[Z] * tdir[X] / tdir[Z];
            s1y[iray] = tpos[Y] - tpos[Z] * tdir[Y] / tdir[Z];
            //printf("s1x[%d], s1y[%d] = %13.6e, %13.6e\n", iray, iray, s1x[iray], s1y[iray]);
            if (iray > PRAY) {
                x = s1x[iray] - s1x[PRAY];
                y = s1y[iray] - s1y[PRAY];
                // sd = sqrt(x * x + y * y);
                sd = upc_hypot(x, y);
                sdmax = (sd > sdmax) ? sd : sdmax;
                pInfo[iray * 2 + 5] = x; // PX
                pInfo[iray * 2 + 6] = y; // PY
            }
        }
        pInfo[0]  = refray[PRAY]->s[0].dir[X]; // s0dir[X]
        pInfo[1]  = refray[PRAY]->s[0].dir[Y]; // s0dir[Y]
        pInfo[2]  = refray[PRAY]->s[0].dir[Z]; // s0dir[Z]
        pInfo[3]  = s1x[PRAY]; // ox
        pInfo[4]  = s1y[PRAY]; // oy
        pInfo[5]  = sdmax;     // sd
        pInfo[6]  = 0.0;       // dmmy
        if (sdmax) {
            pInfo[ 7] /= sdmax; // YUPP PX
            pInfo[ 8] /= sdmax; //      PY
            pInfo[ 9] /= sdmax; // YLOW PX
            pInfo[10] /= sdmax; //      PY
            pInfo[11] /= sdmax; // XUPP PX
            pInfo[12] /= sdmax; //      PY
            pInfo[13] /= sdmax; // XLOW PX
            pInfo[14] /= sdmax; //      PY
        }
    }
}

/* 瞳座標(px, py)を与えて光線追跡する.(無限遠物体用)
 * (注) self->refRayStatus[iwav][ihgt]がFALSEでないことを確かめてから実行すること.
 */
static Bool UpcTraceFactory_traceP_infObjD(const UpcTraceFactory *self, const UpcLens *lobj, int iwav, int ihgt, double px, double py, UpcRay *ray)
{
    double *s0dir = self->pInfo[iwav][ihgt];
    double ox = self->pInfo[iwav][ihgt][3];
    double oy = self->pInfo[iwav][ihgt][4];
    double sd = self->pInfo[iwav][ihgt][5];
    double tpos[3];

    tpos[X] = sd * px + ox;
    tpos[Y] = sd * py + oy;
    tpos[Z] = 0.0;
    return UpcRayTrace_trace1(lobj, s0dir, tpos, iwav, 0, ray);
}

/* 入射基準面の設定(有限物体距離用)
 * (注) 主光線, マージナル光線をサーチしてから実行すること.
 */
static void UpcTraceFactory_setPInfo_finiteObjD(UpcTraceFactory *self)
{
    int i, iwav, ihgt, iray;
    const int n = self->nwav * self->nhgt;
    const int nhgt = self->nhgt;
    const int nray = 5;
    UpcRay **refray;
    UpcRayComponent *pray0; // 主光線の第0面の情報
    double x, y, sd, sdmax;
    double *pInfo;

    for (i = 0; i < n; i++) {
        ihgt = i % nhgt;
        iwav = i / nhgt;
        refray = self->refRay[iwav][ihgt];
        pray0 = refray[PRAY]->s;
        sdmax = 0.0;
        pInfo = self->pInfo[iwav][ihgt];
            
        if (!(self->refRayStatus[iwav][ihgt]))
            continue;
        for (iray = 1; iray < nray; iray++) {
            x = refray[iray]->s[0].dir[X] - pray0->dir[X];
            y = refray[iray]->s[0].dir[Y] - pray0->dir[Y];
            // sd = sqrt(x * x + y * y);
            sd = upc_hypot(x, y);
            sdmax = (sd > sdmax) ? sd : sdmax;
            pInfo[iray * 2 + 5] = x; // PX
            pInfo[iray * 2 + 6] = y; // PY
        }
        pInfo[0]  = pray0->pos[X]; // s0pos[X]
        pInfo[1]  = pray0->pos[Y]; // s0pos[Y]
        pInfo[2]  = pray0->pos[Z]; // s0pos[Z]
        pInfo[3]  = pray0->dir[X]; // s0dir[X]
        pInfo[4]  = pray0->dir[Y]; // s0dir[Y]
        pInfo[5]  = pray0->dir[Z]; // s0dir[Z]
        pInfo[6]  = sdmax; // sd
        if (sdmax) {
            pInfo[ 7] /= sdmax; // YUPP PX
            pInfo[ 8] /= sdmax; //      PY
            pInfo[ 9] /= sdmax; // YLOW PX
            pInfo[10] /= sdmax; //      PY
            pInfo[11] /= sdmax; // XUPP PX
            pInfo[12] /= sdmax; //      PY
            pInfo[13] /= sdmax; // XLOW PX
            pInfo[14] /= sdmax; //      PY
        }
    }
}

/* 瞳座標(px, py)を与えて光線追跡する.(有限物体距離用)
 * (注) self->refRayStatus[iwav][ihgt]がFALSEでないことを確かめてから実行すること.
 */
static Bool UpcTraceFactory_traceP_finiteObjD(const UpcTraceFactory *self, const UpcLens *lobj, int iwav, int ihgt, double px, double py, UpcRay *ray)
{
    double *s0pos = self->pInfo[iwav][ihgt];
    double *s0dir = self->pInfo[iwav][ihgt] + 3;
    double sd = self->pInfo[iwav][ihgt][6];

    return UpcRayTrace_trace0p(lobj, s0pos, s0dir, sd * px, sd * py, iwav, ray);
}    

/* 波面収差計算用パラメータの設定
 * (注) 近軸量と主光線サーチが終わっていること.
 */
static Bool UpcTraceFactory_setWInfo(UpcTraceFactory *self, const UpcLens *lobj)
{
    int i, iwav, ihgt;
    const int n = lobj->nwav * lobj->nhgt;
    const int nhgt = lobj->nhgt;
    double ensPos[3], ensNV[3], ensCurv, exsPos[3], exsNV[3], exsCurv, len, len_inv;
    double *prayIpos;
    double *wInfo;
    UpcRay *pray;

    for (i = 0; i < n; i++) {
        ihgt = i % nhgt;
        iwav = i / nhgt;
        UpcVector_SET(ensPos, 0.0, 0.0, lobj->parax[iwav].enp);
        UpcVector_SET(exsPos, 0.0, 0.0, lobj->parax[iwav].exp);
        pray = self->refRay[iwav][ihgt][PRAY];
        prayIpos = pray->s[lobj->nsur + 1].pos;
        wInfo = self->wInfo[iwav][ihgt];

        if (!(self->refRayStatus[iwav][ihgt])) {
            self->readyWabc[iwav][ihgt] = FALSE;
            continue;
        }
        if (lobj->infObjD) {
            double *prayOdir = pray->s[0].dir;
                
            UpcVector_COPY(ensNV, prayOdir);
            ensCurv = 0.0;
        }
        else {
            double *prayOpos = pray->s[0].pos;

            ensNV[X] = ensPos[X] - prayOpos[X];
            ensNV[Y] = ensPos[Y] - prayOpos[Y];
            ensNV[Z] = ensPos[Z] - prayOpos[Z] + lobj->s[0]->d;
            len = sqrt(UpcVector_DOT(ensNV, ensNV));
            if (len)
                len_inv = ((ensNV[Z] > 0.0) ? 1.0 : -1.0) / len;
            else {
                self->readyWabc[iwav][ihgt] = FALSE;
                continue;
            }
            ensNV[X] *= len_inv;
            ensNV[Y] *= len_inv;
            ensNV[Z] *= len_inv;
            ensCurv = -len_inv;
        }
        exsNV[X] = exsPos[X] - prayIpos[X];
        exsNV[Y] = exsPos[Y] - prayIpos[Y];
        exsNV[Z] = exsPos[Z] - prayIpos[Z] - lobj->s[lobj->nsur]->d;
        len = sqrt(UpcVector_DOT(exsNV, exsNV));
        if (len)
            len_inv = ((exsNV[Z] > 0.0) ? 1.0 : -1.0) / len;
        else {
            self->readyWabc[iwav][ihgt] = FALSE;
            continue;
        }
        exsNV[X] *= len_inv;
        exsNV[Y] *= len_inv;
        exsNV[Z] *= len_inv;
        exsCurv = -len_inv;
        wInfo[ 0] = ensPos[X];
        wInfo[ 1] = ensPos[Y];
        wInfo[ 2] = ensPos[Z];
        wInfo[ 3] = ensNV[X];
        wInfo[ 4] = ensNV[Y];
        wInfo[ 5] = ensNV[Z];
        wInfo[ 6] = ensCurv;
        wInfo[ 7] = exsPos[X];
        wInfo[ 8] = exsPos[Y];
        wInfo[ 9] = exsPos[Z];
        wInfo[10] = exsNV[X];
        wInfo[11] = exsNV[Y];
        wInfo[12] = exsNV[Z];
        wInfo[13] = exsCurv;
        if (UpcRayTrace_calWab(lobj, pray, ensPos, ensNV, ensCurv, exsPos, exsNV, exsCurv, iwav, wInfo + 14))
            self->readyWabc[iwav][ihgt] = TRUE;
        else
            self->readyWabc[iwav][ihgt] = FALSE;
    }
    return TRUE;
}

/* 波長番号, 画角番号, 瞳座標を与えて波面収差を計算する.
 * (注) lobj->raytracef->readyWabc[iwav][ihgt]がFALSEでないことを確かめてから実行すること.
 */
Bool UpcTraceFactory_wabc(const UpcLens *lobj, int iwav, int ihgt, double px, double py, double *wab)
{
    UpcTraceFactory *self = lobj->raytracef;
    UpcRay *ray = UpcRay_init(lobj->nsur);
    Bool flg;
    double *wInfo = self->wInfo[iwav][ihgt];
    
    if (!ray)
        return FALSE;
    
    if (!(self->traceP(self, lobj, iwav, ihgt, px, py, ray))) {
        ray->_dealloc(ray);
        return FALSE;
    }
    flg = UpcRayTrace_calWab(lobj, ray, wInfo, wInfo + 3, wInfo[6], wInfo + 7, wInfo + 10, wInfo[13], iwav, wab);
    *wab -= wInfo[14];
    ray->_dealloc(ray);
    return flg;
}

/* 与えられた画角Y座標の主光線をサーチし像面湾曲とディストーションを計算する(メリディオナル断面限定).
 * [入力]
 *     self : UpcTraceFactoryオブジェクト
 *     lobj : UpcLensオブジェクト
 *     y : 有限物体距離の場合 y物高[mm]. 無限遠物体の場合 tan(U1)
 *     iwav : 波長番号
 * [出力]
 *     f_meri : メリディオナルフォーカス[mm]
 *     f_sagi : サジタルフォーカス[mm]
 *     dist : ディストーション[%]
 * [注意]
 *     self->refRayStatus[iwav][ihgt]がTRUEであることを確かめてから実行すること.
 *     サーチした主光線がメリディオナル面内にない場合のこのルーチンの挙動は保証外.
 */
Bool UpcTraceFactory_asTrace(const UpcLens *lobj, double y, int iwav, double *f_meri, double *f_sagi, double *dist)
{
    const double dk = lobj->s[lobj->nsur]->d;
    const double sk = lobj->parax[lobj->iwav_pri].sk;
    const double exp = lobj->parax[lobj->iwav_pri].exp;
    const double tol = lobj->tol_pray * 10.0;
    UpcRay *ray = UpcRay_init(lobj->nsur);
    double h1, u0, hm, um, hs, us;
    Bool rtnFlg = FALSE;
    double himg, yideal, yimg;
    UpcRayComponent *rt;
    
    if (!ray)
        return FALSE;
    
    if (lobj->infObjD) {
        double dy = y / sqrt(y * y + 1.0);
        
        if (lobj->pray == UpcPrayMode_stop) {
            if (!UpcRayTrace_psearch1(lobj, 0.0, dy, lobj->isur_sto, 0.0, 0.0, iwav, tol, ray))
                goto finally_return;
        }
        else {
            // 無限遠物体のときはテレセントリック指定できない
            goto finally_return;
        }
        himg = lobj->parax[lobj->iwav_pri].fl * y;
        h1 = 1.0;
        u0 = 0.0;
    }
    else {
        if (lobj->pray == UpcPrayMode_stop) {
            if (!UpcRayTrace_psearch0(lobj, 0.0, y, lobj->isur_sto, 0.0, 0.0, iwav, tol, ray))
                goto finally_return;
        }
        else if (lobj->pray == UpcPrayMode_itel) {
            if (!UpcRayTrace_psearch0_xt(lobj, 0.0, y, 0.0, 0.0, iwav, tol, ray))
                goto finally_return;
        }
        else if (lobj->pray == UpcPrayMode_otel) {
            double s0pos[] = {0.0, 0.0, 0.0};
            double s0dir[] = {0.0, 0.0, 1.0};
            
            s0pos[Y] = y;
            lobj->s[0]->shape->getZ(lobj->s[0]->shape, s0pos[X], s0pos[Y], s0pos + Z);
            if (!UpcRayTrace_trace0d(lobj, s0pos, s0dir, iwav, 0, ray))
                goto finally_return;
        }
        himg = lobj->parax[lobj->iwav_pri].mag * y;
        h1 = ray->s[0].opl;
        u0 = -1.0;
    }
    
    rt = ray->s + ray->nsur; // 最終面の光線情報
    yideal = himg * (dk - exp) / (sk - exp); // S定義
    // yideal = himg * (rt->dir[Y] / rt->dir[Z]) * (dk - sk); // C定義
    yimg = rt->pos[Y] + (rt->dir[Y] / rt->dir[Z]) * (dk - rt->pos[Z]);
    *dist = (yideal) ? 100.0 * (yimg - yideal) / yideal : 0.0;
    if (!UpcRayTrace_calAstigmatism(lobj, ray, h1, u0, iwav, &hm, &um, &hs, &us))
        goto finally_return;
    if (um && us) {
        *f_meri = (hm * rt->cosO / um) * rt->dir[Z] + rt->pos[Z] - dk;
        *f_sagi = (hs / us) * rt->dir[Z] + rt->pos[Z] - dk;
        rtnFlg = TRUE;
    }
    
finally_return:
    ray->_dealloc(ray);
    return rtnFlg;
}

