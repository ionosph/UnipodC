//
//  UpcCoord.c
//  UnipodC
//
//  Created by ionosph on 2013/01/28.
//  Copyright (c) 2013年 ionosph. All rights reserved.
//

#include "UpcCoord.h"
#include "UpcVector.h"


static void UpcCoord_dealloc(UpcCoord *self)
{
    free(self);
}

UpcCoord *UpcCoord_init(void)
{
    UpcCoord *self = (UpcCoord *)calloc(sizeof(UpcCoord), 1);
    
    if (self) {
        self->_refCount = 1;
        self->_dealloc = UpcCoord_dealloc;
        self->ex[X] = 1.0;
        self->ey[Y] = 1.0;
        self->ez[Z] = 1.0;
        UpcCoord_setGlobal(self);
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcCoord_init");
    }
    return self;
}

void UpcCoord_setOrigin(UpcCoord *self, const double *p)
{
    UpcVector_COPY(self->eo, p);
    UpcCoord_setGlobal(self);
}

void UpcCoord_setGlobal(UpcCoord *self)
{
    double oo[] = {0.0, 0.0, 0.0};
    double ox[] = {1.0, 0.0, 0.0};
    double oy[] = {0.0, 1.0, 0.0};
    double oz[] = {0.0, 0.0, 1.0};
    
    UpcCoord_toLocalPoint(self, oo, self->go);
    UpcCoord_toLocalVector(self, ox, self->gx);
    UpcCoord_toLocalVector(self, oy, self->gy);
    UpcCoord_toLocalVector(self, oz, self->gz);
}

void UpcCoord_toLocalPoint(UpcCoord *self, const double *xyzGlobal, double *xyzLocal)
{
    double t[3];
    const double *eo = self->eo;
    const double *ex = self->ex;
    const double *ey = self->ey;
    const double *ez = self->ez;
    
    UpcVector_SUB(t, xyzGlobal, eo);
    xyzLocal[X] = UpcVector_DOT(t, ex);
    xyzLocal[Y] = UpcVector_DOT(t, ey);
    xyzLocal[Z] = UpcVector_DOT(t, ez);
}

void UpcCoord_toLocalVector(UpcCoord *self, const double *xyzGlobal, double *xyzLocal)
{
    const double *ex = self->ex;
    const double *ey = self->ey;
    const double *ez = self->ez;
    double t[3];
    
    UpcVector_COPY(t, xyzGlobal); // xyzLocalとxyzGlobalが同じポインタだったときのための処置
    xyzLocal[X] = UpcVector_DOT(t, ex);
    xyzLocal[Y] = UpcVector_DOT(t, ey);
    xyzLocal[Z] = UpcVector_DOT(t, ez);
}

void UpcCoord_toGlobalPoint(UpcCoord *self, const double *xyzLocal, double *xyzGlobal)
{
    const double *go = self->go;
    const double *gx = self->gx;
    const double *gy = self->gy;
    const double *gz = self->gz;
    double t[3];
    
    UpcVector_SUB(t, xyzLocal, go);
    xyzGlobal[X] = UpcVector_DOT(t, gx);
    xyzGlobal[Y] = UpcVector_DOT(t, gy);
    xyzGlobal[Z] = UpcVector_DOT(t, gz);
}

void UpcCoord_toGlobalVector(UpcCoord *self, const double *xyzLocal, double *xyzGlobal)
{
    const double *gx = self->gx;
    const double *gy = self->gy;
    const double *gz = self->gz;
    double t[3];

    UpcVector_COPY(t, xyzLocal); // xyzLocalとxyzGlobalが同じポインタだったときのための処置
    xyzGlobal[X] = UpcVector_DOT(t, gx);
    xyzGlobal[Y] = UpcVector_DOT(t, gy);
    xyzGlobal[Z] = UpcVector_DOT(t, gz);
}

void UpcCoord_rot(UpcCoord *self, double *rotAxis, double angle_radian)
{
    UpcVector_rot(self->ex, rotAxis, angle_radian);
    UpcVector_rot(self->ey, rotAxis, angle_radian);
    UpcVector_rot(self->ez, rotAxis, angle_radian);
    UpcCoord_setGlobal(self);
}

void UpcCoord_rotByName(UpcCoord *self, const char *rotAxisName, double angle_radian)
{
    double axis[3];
    
    if (rotAxisName[0] == 'X' || rotAxisName[0] == 'x') {
        axis[X] = 1.0;
        axis[Y] = 0.0;
        axis[Z] = 0.0;
    }
    else if (rotAxisName[0] == 'Y' || rotAxisName[0] == 'y') {
        axis[X] = 0.0;
        axis[Y] = 1.0;
        axis[Z] = 0.0;
    }
    else if (rotAxisName[0] == 'Z' || rotAxisName[0] == 'z') {
        axis[X] = 0.0;
        axis[Y] = 0.0;
        axis[Z] = 1.0;
    }
    else if (!strcmp(rotAxisName, "localX") || !strcmp(rotAxisName, "ex")) {
        UpcVector_COPY(axis, self->ex);
    }
    else if (!strcmp(rotAxisName, "localY") || !strcmp(rotAxisName, "ey")) {
        UpcVector_COPY(axis, self->ey);
    }
    else if (!strcmp(rotAxisName, "localZ") || !strcmp(rotAxisName, "ez")) {
        UpcVector_COPY(axis, self->ez);
    }
    else {
        UpcERRHANDLER(UpcE_ValueError, "in UpcCoord_rotByName");
    }
    UpcCoord_rot(self, axis, angle_radian);
}

/* ローカル座標系のZ軸を, グローバル座標系で表された点pを向くように, ローカル座標系を回転する.
 */
Bool UpcCoord_setZAxisTo(UpcCoord *self, const double *p)
{
    double direc[3];
    double phi, theta;
    
    UpcVector_SUB(direc, p, self->eo);
    
    phi = atan2(direc[Y], direc[X]) + PI / 2.0;
    UpcCoord_rotByName(self, "Z", phi);
    theta = - atan(direc[Z] / upc_hypot(direc[X], direc[Y])) + PI / 2.0;
    UpcCoord_rotByName(self, "localX", theta);
    
    //Ux_vector_normalize(direc);
    //printf("(%14.6f, %14.6f, %14.6f)\n", direc[X], direc[Y], direc[Z]);
    //printf("(%14.6f, %14.6f, %14.6f)\n", self->ez[X], self->ez[Y], self->ez[Z]);
    //Ux_error("stop");
    return TRUE;
}




