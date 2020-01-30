/*
 *  UpcVector.c
 *  UnipodC
 *
 *  Created by ionosph on 2010/10/20.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcVector.h"



void UpcVector_set(double *v, double x, double y, double z)
{
    v[X] = x;
    v[Y] = y;
    v[Z] = z;
}


/*  ベクトルの線形結合を求める.
 */
void UpcVector_linComb(double *v, double c1, const double *v1, double c2, const double *v2)
{
    v[X] = c1 * v1[X] + c2 * v2[X];
    v[Y] = c1 * v1[Y] + c2 * v2[Y];
    v[Z] = c1 * v1[Z] + c2 * v2[Z];
}

/*  ベクトルの長さを求める.
 */
double UpcVector_norm(const double *v)
{
    return sqrt(v[X] * v[X] + v[Y] * v[Y] + v[Z] * v[Z]);
}

/*  ベクトルを規格化する.
 */
Bool UpcVector_normalize(double *v)
{
    double len2 = v[X] * v[X] + v[Y] * v[Y] + v[Z] * v[Z];
    
    if (len2) {
        double ilen = 1.0 / sqrt(len2);
        v[X] *= ilen;
        v[Y] *= ilen;
        v[Z] *= ilen;
        return TRUE;
    }
    return FALSE;
}

/*  ベクトルのドット積を求める.
 */
double UpcVector_dot(const double *v1, const double *v2)
{
    return v1[X] * v2[X] + v1[Y] * v2[Y] + v1[Z] * v2[Z];
}

/*  ベクトルのクロス積(v1 × v2)を求める.
 */
void UpcVector_cross(double *v, const double *v1, const double *v2)
{
    v[X] = v1[Y] * v2[Z] - v1[Z] * v2[Y];
    v[Y] = v1[Z] * v2[X] - v1[X] * v2[Z];
    v[Z] = v1[X] * v2[Y] - v1[Y] * v2[X];
}

/* ベクトルrを, ベクトルnの周りにφだけ回転させる.
 * ロドリグの公式
 *   r' = r cosφ + n (n・r) (1 - cosφ) - (r×n) sinφ
 * 入力
 *     r[3] : 回転させるベクトル
 *     n[3] : 回転軸を表すベクトル
 *     phi  : 回転角(rad)
 * 出力
 *     戻り値 : 計算できたらTRUE, できなかったらFALSE
 *     r[3] : 回転後のベクトルで上書きされる.
 *     (注) nも正規化され上書きされる.
 */
Bool UpcVector_rot(double *r, double *n, double phi)
{
    double len2 = n[X] * n[X] + n[Y] * n[Y] + n[Z] * n[Z];
    double cos_p, sin_p, t, rxn[3];
    
    if (len2) {
        if (len2 != 1.0) {
            double ilen = 1.0 / sqrt(len2);
            
            n[X] *= ilen;
            n[Y] *= ilen;
            n[Z] *= ilen;
        }
        cos_p = cos(phi);
        sin_p = sin(phi);
        t = (n[X] * r[X] + n[Y] * r[Y] + n[Z] * r[Z]) * (1.0 - cos_p);
        UpcVector_cross(rxn, r, n);
        r[X] = r[X] * cos_p + n[X] * t - rxn[X] * sin_p;
        r[Y] = r[Y] * cos_p + n[Y] * t - rxn[Y] * sin_p;
        r[Z] = r[Z] * cos_p + n[Z] * t - rxn[Z] * sin_p;
        return TRUE;
    }
    else {
        return FALSE;
    }
}

/* 2点間の距離を求める.
 */
double UpcVector_distance(double *v1, double *v2)
{
    double v[3];

    UpcVector_SUB(v, v1, v2);
    return UpcVector_norm(v);
}


#define C11 0
#define C12 1
#define C21 2
#define C22 3


/* 2x2行列の逆行列を求める.
 * m[4] = {a11, a12, a21, a22}
 */
Bool UpcMatrix2_inverse(double *m)
{
    double det = m[C11] * m[C22] - m[C12] * m[C21];
    
    if (det) {
        double tmp, idet = 1.0 / det;

        tmp    = m[C11];
        m[C11] = m[C22];
        m[C22] = tmp ;

        m[C11] *=  idet;
        m[C22] *=  idet;
        m[C12] *= -idet;
        m[C21] *= -idet;
        return TRUE;
    }
    return FALSE;
}



