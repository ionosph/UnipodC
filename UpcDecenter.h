/*
 *  UpcDecenter.h
 *  UnipodC
 *
 *  Created by ionosph on 2010/11/04.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _UPCDECENTER_H
#define _UPCDECENTER_H

#include "UpcBase.h"

typedef struct _UpcDecenter UpcDecenter;
struct _UpcDecenter {
    int _refCount;
    void (*_dealloc)(UpcDecenter *);
    // attribute
    Bool decenterFlg;
    double shift[3];
    double omega[3];
    double cosX;
    double sinX;
    double cosY;
    double sinY;
    double cosZ;
    double sinZ;
    // methods
    UpcDecenter *(*copy)(UpcDecenter *);
    void (*print)(const UpcDecenter *, int, FILE *);
    void (*toXml)(const UpcDecenter *, int, FILE *);
    UpcErrStatus (*lset)(UpcDecenter *);
};

UpcDecenter *UpcDecenter_init(void);
void UpcDecenter_setAllZero(UpcDecenter *self);
//void UpcDecenter_setDecenter(UpcDecenter *self, double shiftX, double shiftY, double shiftZ, double omegaX, double omegaY, double omegaZ);
Bool UpcDecenter_xSymmetric(const UpcDecenter *self);
Bool UpcDecenter_ySymmetric(const UpcDecenter *self);
Bool UpcDecenter_rSymmetric(const UpcDecenter *self);
void UpcDecenter_rotForward(const UpcDecenter *self, const double *vec, double *nVec);
void UpcDecenter_rotReverse(const UpcDecenter *self, const double *vec, double *nVec);
void UpcDecenter_convForward(const UpcDecenter *self, const double *pos, double *nPos);
void UpcDecenter_convReverse(const UpcDecenter *self, const double *pos, double *nPos);



#endif
