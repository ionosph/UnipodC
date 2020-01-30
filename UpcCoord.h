//
//  UpcCoord.h
//  UnipodC
//
//  Created by ionosph on 2013/01/28.
//  Copyright (c) 2013年 ionosph. All rights reserved.
//

#ifndef _UPCCOORD_H
#define _UPCCOORD_H

#include "UpcBase.h"


typedef struct _UpcCoord UpcCoord;
struct _UpcCoord {
    int _refCount;
    void (*_dealloc)(UpcCoord *);
    double eo[3]; // グローバル座標で表したローカル座標系の原点
    double ex[3]; // グローバル座標で表したローカル座標系のx軸
    double ey[3]; // グローバル座標で表したローカル座標系のy軸
    double ez[3]; // グローバル座標で表したローカル座標系のz軸
    double go[3]; // ローカル座標系で表したグローバル座標系の原点
    double gx[3]; // ローカル座標系で表したグローバル座標系のx軸
    double gy[3]; // ローカル座標系で表したグローバル座標系のy軸
    double gz[3]; // ローカル座標系で表したグローバル座標系のz軸
};

UpcCoord *UpcCoord_init(void);
/*
 * UpcCoordオブジェクトを生成する.
 */

void UpcCoord_setOrigin(UpcCoord *self, const double *p);
/*
 * 点pをローカル座標の原点にし, go, gx, gy, gzを再計算する.
 */

void UpcCoord_setGlobal(UpcCoord *self);
/*
 * 現在のローカル座標系の情報を元にgo, gx, gy, gzを算出する.
 */

void UpcCoord_toLocalPoint(UpcCoord *self, const double *xyzGlobal, double *xyzLocal);
/*
 * グローバル座標の点xyzGlobalをローカル座標に変換してxyzLocalに格納する.
 */

void UpcCoord_toLocalVector(UpcCoord *self, const double *xyzGlobal, double *xyzLocal);
/*
 * グローバル座標のベクトルxyzGlobalをローカル座標に変換してxyzLocalに格納する.
 */

void UpcCoord_toGlobalPoint(UpcCoord *self, const double *xyzLocal, double *xyzGlobal);
void UpcCoord_toGlobalVector(UpcCoord *self, const double *xyzLocal, double *xyzGlobal);
void UpcCoord_rot(UpcCoord *self, double *rotAxis, double angle_radian);
void UpcCoord_rotByName(UpcCoord *self, const char *rotAxisName, double angle_radian);
Bool UpcCoord_setZAxisTo(UpcCoord *self, const double *p);








#endif
