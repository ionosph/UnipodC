//
//  UpcZernike.h
//  UnipodC
//
//  Created by ionosph on 2012/12/25.
//  Copyright (c) 2012年 ionosph. All rights reserved.
//

#ifndef _UPCZERNIKE_H
#define _UPCZERNIKE_H

#include "UpcBase.h"

Bool UpcZernike_getZf(int istt, int iend, double x, double y, double *zfarr);
/* 第istt項から第iend項までのZernike関数値をzfarr[istt..iend]に詰めて返す.
 * [入力]
 *     istt : この項から計算開始
 *     iend : この項まで計算する
 *     x : X座標(0 ≦ x ≦ 1)
 *     y : Y座標(0 ≦ y ≦ 1)
 * [出力]
 *     zfarr[0..iend] : zfarr[0]はダミー. zfarr[1]はZ1, zfarr[2]はZ2...
 *                      zfarr[istt]からzfarr[iend]に計算結果が格納される.
 *                      少なくともzfarr[iend]までは配列が確保されている必要がある.
 *     戻り値 : 成功すればTRUE, 失敗したらFALSE
 */

double UpcZernike_zval(int istt, int iend, double x, double y, const double *zcoef);
/* Zernike係数のセットzcoefの, (x, y)における値を返す (第istt項から第iend項までを考慮).
 * [入力]
 *     istt : この項から計算開始
 *     iend : この項まで計算する
 *     x : X座標(0 ≦ x ≦ 1)
 *     y : Y座標(0 ≦ y ≦ 1)
 *     zcoef[0..iend] : zcoef[0]はダミー. zcoef[1]はC1, zcoef[2]はC2...
 *                      少なくともzcoef[iend]までは確保されている必要がある.
 * [出力]
 *     戻り値 : sum(zcoef[i] * zf[i]) (istt <= i <= iend)
 */

#endif
