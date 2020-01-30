/*
 *  UpcZernikeToXYP.h
 *  UnipodC
 *
 *  Created by ionosph on 12/02/09.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _UPCZERNIKETOXYP_H
#define _UPCZERNIKETOXYP_H

#include "UpcBase.h"


/* 偶数のみまたは奇数のみを添字に持つ2次元配列
 */
typedef struct _UpcQMatrix UpcQMatrix;
struct _UpcQMatrix {
    int _refCount;
    void (*_dealloc)(UpcQMatrix *);
    double **data;
    int ixMax;
    int iyMax;
    int xIsOdd; // 添字が奇数のとき1, 偶数のとき0
    int yIsOdd;
};

//UpcQMatrix *UpcQMatrix_init(int ixMax, int iyMax, int xIsOdd, int yIsOdd);
double UpcQMatrix_data(UpcQMatrix *self, int iy, int ix);
//UpcErrStatus UpcQMatrix_setData(UpcQMatrix *self, int iy, int ix, double value);
//UpcErrStatus UpcQMatrix_addData(UpcQMatrix *self, int iy, int ix, double value);
//UpcErrStatus UpcQMatrix_mulData(UpcQMatrix *self, int iy, int ix, double value);
//void UpcQMatrix_print(UpcQMatrix *self, FILE *ostream);

UpcQMatrix *UpcZernikeToXYP_zcoefxy(int cno);
int UpcZernikeToXYP_maxorder(int cno);

Bool testUpcZernikeToXYP(void);


#endif
