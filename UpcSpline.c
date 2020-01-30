/*
 *  UpcSpline.c
 *  UnipodC
 *
 *  Created by ionosph on 11/07/06.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcSpline.h"

/* 各点の2階導関数を計算する.
 * [入力]
 *     n : データ数
 *     xa[0..n - 1] : x座標の配列
 *     ya[0..n - 1] : y座標の配列
 *     dy1 : 左端点の1次導関数. 1.0e30にすると自然スプライン条件
 *     dyN : 右端点の1次導関数. 1.0e30にすると自然スプライン条件
 * [出力]
 *     d2y[0..n - 1] : 各点の2階導関数
 */
static Bool cal_d2y(int n, const double *xa, const double *ya, double dy1, double dyN, double *d2y)
{
    double *c = (double *)malloc(sizeof(double) * (n - 1));
    double h, hi, hj, bj;
    int j;
    
    if (!c) {
        UpcERRHANDLER(UpcE_MemoryError, "in cal_d2y");
        return FALSE;
    }
    // 往路
    // (左端の条件)
    if (dy1 != 1.0e30) {
        h = xa[1] - xa[0];
        c[0] = 0.5;
        d2y[0] = 3.0 / h * ((ya[1] - ya[0]) / h - dy1);
    }
    else {
        c[0] = 0.0;
        d2y[0] = 0.0;
    }
    // (各接続点の条件 (j = 1, 2, ..., n - 2))
    for (j = 1; j < n - 1; j++) {
        hi = xa[j] - xa[j - 1];
        hj = xa[j + 1] - xa[j];
        bj = 2.0 * (hi + hj) - hi * c[j - 1];
        d2y[j] = 6.0 * ((ya[j + 1] - ya[j]) / hj - (ya[j] - ya[j - 1]) / hi) - hi * d2y[j - 1];
        c[j] = hj / bj;
        d2y[j] /= bj;
    }
    // (右端の条件)
    if (dyN != 1.0e30) {
        h = xa[n - 1] - xa[n - 2];
        d2y[n - 1] = 6.0 * (dyN - (ya[n - 1] - ya[n - 2]) / h) - h * d2y[n - 2];
        d2y[n - 1] /= h * (2.0 - c[n - 2]);
    }
    else
        d2y[n - 1] = 0.0;
    // 復路
    for (j = n - 2; j >= 0; j--)
        d2y[j] -= c[j] * d2y[j + 1];
    free(c);
    return TRUE;
}

static void UpcSpline_dealloc(UpcSpline *self)
{
    free(self->xa); // self->ya, self->d2ya もこれで解放される.
    free(self);
}

UpcSpline *UpcSpline_init(int n, const double *xa, const double *ya, double dy1, double dyN)
{
    UpcSpline *self;
    int i;
    
    if (n < 2)
        return NULL;

    self = (UpcSpline *)malloc(sizeof(UpcSpline));

    if (self) {
        self->_refCount = 1;
        self->_dealloc = UpcSpline_dealloc;
        self->ndata = n;
        self->xa = (double *)malloc(sizeof(double) * (3 * n)); // self->ya, self->d2ya の分も同時に確保.
        if (!self->xa) {
            free(self);
            UpcERRHANDLER(UpcE_MemoryError, "in UpcSpline_init");
            return NULL;
        }
        self->ya = self->xa + n; // self->ya は self->xa の後ろ部分を利用.
        self->d2ya = self->xa + 2 * n; // self->d2ya は self->xa の後ろ部分を利用.
        for (i = 0; i < n; i++) {
            self->xa[i] = xa[i];
            self->ya[i] = ya[i];
        }

        upc_heapSort2(n, self->xa - 1, self->ya - 1);
        for (i = 0; i < n - 1; i++) {
            if(self->xa[i] == self->xa[i + 1]) { // 重複するxがある.
                free(self->xa);
                free(self);
                UpcERRHANDLER(UpcE_RuntimeError, "redundant x value");
                return NULL;
            }
        }
        if (!cal_d2y(n, self->xa, self->ya, dy1, dyN, self->d2ya)) {
            free(self->xa);
            free(self);
            UpcERRHANDLER(UpcE_RuntimeError, "cal_d2y failed");
            return NULL;
        }
        self->iprev = 0;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcSpline_init");
        return NULL;
    }
    return self;
}

double UpcSpline_y(UpcSpline *self, double x)
{
    int ilo, ihi, i;
    double xlo, xhi, h, a, b;
    
    if (x >= self->xa[self->iprev] && x < self->xa[self->iprev + 1]) {
        ilo = self->iprev;
        ihi = ilo + 1;
    }
    else {
        ilo = 0;
        ihi = self->ndata - 1;
        while (ihi - ilo > 1) {
            i = (ihi + ilo) / 2;
            if (x < self->xa[i])
                ihi = i;
            else
                ilo = i;
        }
        self->iprev = ilo;
    }
    xlo = self->xa[ilo];
    xhi = self->xa[ihi];
    h = xhi - xlo;
    a = (xhi - x) / h;
    b = (x - xlo) / h;
    return a * self->ya[ilo] + b * self->ya[ihi] + ((a * a * a - a) * self->d2ya[ilo] + (b * b * b - b) * self->d2ya[ihi]) * (h * h) / 6.0;
}

