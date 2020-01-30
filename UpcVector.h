/*
 *  UpcVector.h
 *  UnipodC
 *
 *  Created by ionosph on 2010/10/20.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 *  3要素のdouble型配列をベクトルとして扱い演算するための関数群
 *
 */

#ifndef _UPCVECTOR_H
#define _UPCVECTOR_H

#include "UpcBase.h"

#define UpcVector_SET(v,x,y,z) do{(v)[X]=(x);(v)[Y]=(y);(v)[Z]=(z);}while(0)
//#define UpcVector_COPY(v,w) do{(v)[X]=(w)[X];(v)[Y]=(w)[Y];(v)[Z]=(w)[Z];}while(0)

// wの内容をvにコピーする
#define UpcVector_COPY(v,w) memcpy((v),(w),sizeof(double)*3)

// v = v1 + v2
#define UpcVector_ADD(v,v1,v2) do{(v)[X]=(v1)[X]+(v2)[X];(v)[Y]=(v1)[Y]+(v2)[Y];(v)[Z]=(v1)[Z]+(v2)[Z];}while(0)

// v = v1 - v2
#define UpcVector_SUB(v,v1,v2) do{(v)[X]=(v1)[X]-(v2)[X];(v)[Y]=(v1)[Y]-(v2)[Y];(v)[Z]=(v1)[Z]-(v2)[Z];}while(0)

// v = c1 * v1 + c2 * v2
#define UpcVector_LINCOMB(v,c1,v1,c2,v2) do{(v)[X]=(c1)*(v1)[X]+(c2)*(v2)[X];(v)[Y]=(c1)*(v1)[Y]+(c2)*(v2)[Y];(v)[Z]=(c1)*(v1)[Z]+(c2)*(v2)[Z];}while(0)

// v = v1とv2の内積
#define UpcVector_DOT(v1,v2) (((v1)[X]*(v2)[X])+((v1)[Y]*(v2)[Y])+((v1)[Z]*(v2)[Z]))

void UpcVector_set(double *v, double x, double y, double z);
void UpcVector_linComb(double *v, double c1, const double *v1, double c2, const double *v2);
double UpcVector_norm(const double *v);
Bool UpcVector_normalize(double *v);
double UpcVector_dot(const double *v1, const double *v2);
void UpcVector_cross(double *v, const double *v1, const double *v2);
double UpcVector_distance(double *v1, double *v2);
Bool UpcVector_rot(double *r, double *n, double phi);

Bool UpcMatrix2_inverse(double *m);
    
#endif
