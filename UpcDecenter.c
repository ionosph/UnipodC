/*
 *  UpcDecenter.c
 *  UnipodC
 *
 *  Created by ionosph on 2010/11/04.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcDecenter.h"
#include "UpcVector.h"

static const double shiftDefault[] = {0.0, 0.0, 0.0};
static const double omegaDefault[] = {0.0, 0.0, 0.0};

static void UpcDecenter_dealloc(UpcDecenter *self)
{
    free(self);
}

static void UpcDecenter_print(const UpcDecenter *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);

    fprintf(ostream, "%sflg   = %s\n", space, BOOLTOSTRING(self->decenterFlg));
    fprintf(ostream, "%sshift = (%13.6e, %13.6e, %13.6e)\n", space, self->shift[X], self->shift[Y], self->shift[Z]);
    fprintf(ostream, "%somega = (%13.6e, %13.6e, %13.6e)\n", space, self->omega[X], self->omega[Y], self->omega[Z]);
}

static void UpcDecenter_toXml(const UpcDecenter *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);

    fprintf(ostream, "%s<decenter>\n", space);
    if (self->shift[X] != shiftDefault[X])
        fprintf(ostream, "%s  <shift_x>%23.15e</shift_x>\n", space, self->shift[X]);
    if (self->shift[Y] != shiftDefault[Y])
        fprintf(ostream, "%s  <shift_y>%23.15e</shift_y>\n", space, self->shift[Y]);
    if (self->shift[Z] != shiftDefault[Z])
        fprintf(ostream, "%s  <shift_z>%23.15e</shift_z>\n", space, self->shift[Z]);
    if (self->omega[X] != omegaDefault[X])
        fprintf(ostream, "%s  <omega_x>%23.15e</omega_x>\n", space, self->omega[X]);
    if (self->omega[Y] != omegaDefault[Y])
        fprintf(ostream, "%s  <omega_y>%23.15e</omega_y>\n", space, self->omega[Y]);
    if (self->omega[Z] != omegaDefault[Z])
        fprintf(ostream, "%s  <omega_z>%23.15e</omega_z>\n", space, self->omega[Z]);
    fprintf(ostream, "%s</decenter>\n", space);
}

static UpcErrStatus UpcDecenter_lset(UpcDecenter *self)
{
    double t;

    if (self->omega[X]) {
        t = RADIANS(self->omega[X]);
        self->cosX = cos(t);
        self->sinX = sin(t);
    }
    else {
        self->cosX = 1.0;
        self->sinX = 0.0;
    }

    if (self->omega[Y]) {
        t = RADIANS(self->omega[Y]);
        self->cosY = cos(t);
        self->sinY = sin(t);
    }
    else {
        self->cosY = 1.0;
        self->sinY = 0.0;
    }

    if (self->omega[Z]) {
        t = RADIANS(self->omega[Z]);
        self->cosZ = cos(t);
        self->sinZ = sin(t);
    }
    else {
        self->cosZ = 1.0;
        self->sinZ = 0.0;
    }

    if (self->shift[X] || self->shift[Y] || self->shift[Z] || 
        self->omega[X] || self->omega[Y] || self->omega[Z]) {
        self->decenterFlg = TRUE;
    }
    else {
        self->decenterFlg = FALSE;
    }

    return UpcE_NoIssues;
}

static UpcDecenter *UpcDecenter_copy(UpcDecenter *self)
{
    UpcDecenter *t = (UpcDecenter *)malloc(sizeof(UpcDecenter));
    
    if (t) {
        memcpy(t, self, sizeof(UpcDecenter));
        t->_refCount = 1;
        /*
        t->_refCount = 1;
        t->_dealloc = UpcDecenter_dealloc;
        t->decenterFlg = self->decenterFlg;
        t->shift[X] = self->shift[X];
        t->shift[Y] = self->shift[Y];
        t->shift[Z] = self->shift[Z];
        t->omega[X] = self->omega[X];
        t->omega[Y] = self->omega[Y];
        t->omega[Z] = self->omega[Z];
        t->cosX = self->cosX;
        t->sinX = self->sinX;
        t->cosY = self->cosY;
        t->sinY = self->sinY;
        t->cosZ = self->cosZ;
        t->sinZ = self->sinZ;
        t->copy = UpcDecenter_copy;
        t->print = UpcDecenter_print;
        t->toXml = UpcDecenter_toXml;
        t->lset = UpcDecenter_lset;
         */
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcDecenter_copy");
    }
    return t;
}

UpcDecenter *UpcDecenter_init(void)
{
    UpcDecenter *self = (UpcDecenter *)calloc(1, sizeof(UpcDecenter));

    if (self) {
        self->_refCount = 1;
        //UpcDecenter_setDecenter(self, shiftDefault[X], shiftDefault[Y], shiftDefault[Z], omegaDefault[X], omegaDefault[Y], omegaDefault[Z]);
        self->_dealloc = UpcDecenter_dealloc;
        self->copy = UpcDecenter_copy;
        self->print = UpcDecenter_print;
        self->toXml = UpcDecenter_toXml;
        self->lset = UpcDecenter_lset;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcDecenter_init");
    }
    return self;
}

/*
void UpcDecenter_setDecenter(UpcDecenter *self, double shiftX, double shiftY, double shiftZ, double omegaX, double omegaY, double omegaZ)
{
    double t;

    self->shift[X] = shiftX;
    self->shift[Y] = shiftY;
    self->shift[Z] = shiftZ;
    self->omega[X] = omegaX;
    self->omega[Y] = omegaY;
    self->omega[Z] = omegaZ;

    if (omegaX) {
        t = RADIANS(omegaX);
        self->cosX = cos(t);
        self->sinX = sin(t);
    }
    else {
        self->cosX = 1.0;
        self->sinX = 0.0;
    }
    if (omegaY) {
        t = RADIANS(omegaY);
        self->cosY = cos(t);
        self->sinY = sin(t);
    }
    else {
        self->cosY = 1.0;
        self->sinY = 0.0;
    }
    if (omegaZ) {
        t = RADIANS(omegaZ);
        self->cosZ = cos(t);
        self->sinZ = sin(t);
    }
    else {
        self->cosZ = 1.0;
        self->sinZ = 0.0;
    }
    if (self->shift[X] || self->shift[Y] || self->shift[Z] || 
        self->omega[X] || self->omega[Y] || self->omega[Z]) {
        self->decenterFlg = TRUE;
    }
    else {
        self->decenterFlg = FALSE;
    }
}
*/

void UpcDecenter_setAllZero(UpcDecenter *self)
{
    self->shift[X] = 0.0;
    self->shift[Y] = 0.0;
    self->shift[Z] = 0.0;
    self->omega[X] = 0.0;
    self->omega[Y] = 0.0;
    self->omega[Z] = 0.0;
}

Bool UpcDecenter_xSymmetric(const UpcDecenter *self)
{
    if (self->shift[Y] || self->omega[X] || self->omega[Z])
        return FALSE;
    return TRUE;
}

Bool UpcDecenter_ySymmetric(const UpcDecenter *self)
{
    if (self->shift[X] || self->omega[Y] || self->omega[Z])
        return FALSE;
    return TRUE;
}

Bool UpcDecenter_rSymmetric(const UpcDecenter *self)
{
    if (self->decenterFlg)
        return FALSE;
    return TRUE;
}

/*  vecを(ωx)→(ωy)→(ωz)の順に回転する.
    すなわち, 偏心前の座標系で表された方位ベクトルを偏心後の座標系で表す.
*/
void UpcDecenter_rotForward(const UpcDecenter *self, const double *vec, double *nVec)
{
    double vx = vec[X];
    double vy = vec[Y];
    double vz = vec[Z];

    if (self->omega[X]) {
        double vny = self->cosX * vy + self->sinX * vz;
        double vnz = self->cosX * vz - self->sinX * vy;
        vy = vny;
        vz = vnz;
    }
    if (self->omega[Y]) {
        double vnz = self->cosY * vz + self->sinY * vx;
        double vnx = self->cosY * vx - self->sinY * vz;
        vz = vnz;
        vx = vnx;
    }
    if (self->omega[Z]) {
        double vnx = self->cosZ * vx + self->sinZ * vy;
        double vny = self->cosZ * vy - self->sinZ * vx;
        vx = vnx;
        vy = vny;
    }
    UpcVector_SET(nVec, vx, vy, vz);
}

/*  vecを(-ωz)→(-ωy)→(-ωx)の順に逆回転する.
    すなわち, 偏心後の座標系で表された方位ベクトルを偏心前の座標系で表す.
*/
void UpcDecenter_rotReverse(const UpcDecenter *self, const double *vec, double *nVec)
{
    double vx = vec[X];
    double vy = vec[Y];
    double vz = vec[Z];

    if (self->omega[Z]) {
        double vnx = self->cosZ * vx - self->sinZ * vy;
        double vny = self->cosZ * vy + self->sinZ * vx;
        vx = vnx;
        vy = vny;
    }
    if (self->omega[Y]) {
        double vnz = self->cosY * vz - self->sinY * vx;
        double vnx = self->cosY * vx + self->sinY * vz;
        vz = vnz;
        vx = vnx;
    }
    if (self->omega[X]) {
        double vny = self->cosX * vy - self->sinX * vz;
        double vnz = self->cosX * vz + self->sinX * vy;
        vy = vny;
        vz = vnz;
    }
    UpcVector_SET(nVec, vx, vy, vz);
}

/*  偏心前の座標系で表された点を偏心後の座標系で表す.
*/
void UpcDecenter_convForward(const UpcDecenter *self, const double *pos, double *nPos)
{
    double t[3];

    t[X] = pos[X] - self->shift[X];
    t[Y] = pos[Y] - self->shift[Y];
    t[Z] = pos[Z] - self->shift[Z];
    UpcDecenter_rotForward(self, t, nPos);
}

/*  偏心後の座標系で表された点を偏心前の座標系で表す.
*/
void UpcDecenter_convReverse(const UpcDecenter *self, const double *pos, double *nPos)
{
    double t[3];

    UpcDecenter_rotReverse(self, pos, t);
    nPos[X] = t[X] + self->shift[X];
    nPos[Y] = t[Y] + self->shift[Y];
    nPos[Z] = t[Z] + self->shift[Z];
}
