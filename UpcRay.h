/*
 *  UpcRay.h
 *  UnipodC
 *
 *  Created by ionosph on 2011/04/05.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _UPCRAY_H
#define _UPCRAY_H

#include "UpcBase.h"



typedef enum _UpcRayStatus {
    UpcRayE_Succeeded            =  2,
    UpcRayE_OutOfArea            =  1,
    UpcRayE_UnCalculated         =  0,
    UpcRayE_NoIntersectionError  = -1,
    UpcRayE_TotalReflectionError = -2,
    UpcRayE_IterationError       = -3,
    UpcRayE_SettingError         = -4
} UpcRayStatus;

typedef struct {
    double pos[3]; // 光線位置
    double dir[3]; // 出射光線の方向余弦
    double opl;    // 光路に沿った距離(次の面まで)
    double cosI;   // 入射角の余弦
    double cosO;   // 出射角の余弦
    UpcRayStatus status; // 光線追跡ステート
} UpcRayComponent;

typedef struct _UpcRay UpcRay;
struct _UpcRay {
    int _refCount;
    void (*_dealloc)(UpcRay *);
    int nsur;
    UpcRayStatus status;
    UpcRayComponent *s;
};

UpcRay *UpcRay_init(int nsur);
//void UpcRay_dealloc(UpcRay *self);
void UpcRay_print(const UpcRay *self, char *comment, FILE *ostream);
void UpcRay_copy(UpcRay *self, const UpcRay *other);
void UpcRay_copy_YZmirror(UpcRay *self, const UpcRay *other);
void UpcRay_copy_XZmirror(UpcRay *self, const UpcRay *other);
void UpcRay_resetStatus(UpcRay *self);

const char *UpcRayStatus_string(UpcRayStatus status);

#endif
