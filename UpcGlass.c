/*
 *  UpcGlass.c
 *  UnipodC
 *
 *  Created by ionosph on 11/07/07.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcGlass.h"
#include "UpcDictDD.h"
#include "UpcSpline.h"


UpcErrStatus UpcGlass_lsetDefault(UpcGlass *self)
{
    return UpcE_NoIssues;
}

//------------------------------------------------------------
// UpcGlassSellmeier

static void UpcGlassSellmeier_dealloc(UpcGlass *self)
{
    free(self->name);
    free(self->indexF); // indexBも同時に解放される.
    UpcObj_RELEASE(self->cash);
    free(self);
}

static void UpcGlassSellmeier_toXml(const UpcGlass *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);

    fprintf(ostream, "%s<glass type=\"sellmeier\">\n", space);
    fprintf(ostream, "%s  <name>%s</name>\n", space, self->name);
    fprintf(ostream, "%s  <a1>%23.15e</a1>\n", space, self->coef[0]);
    fprintf(ostream, "%s  <a2>%23.15e</a2>\n", space, self->coef[1]);
    fprintf(ostream, "%s  <a3>%23.15e</a3>\n", space, self->coef[2]);
    fprintf(ostream, "%s  <b1>%23.15e</b1>\n", space, self->coef[3]);
    fprintf(ostream, "%s  <b2>%23.15e</b2>\n", space, self->coef[4]);
    fprintf(ostream, "%s  <b3>%23.15e</b3>\n", space, self->coef[5]);
    fprintf(ostream, "%s</glass>\n", space);
}

/* 与えた波長における屈折率を返す.
 * (注) 波長はnm単位で与える.
 */
static double UpcGlassSellmeier_n(UpcGlass *self, double wl)
{
    double ind;

    if (UpcDictDD_search(self->cash, wl, &ind))
        return ind;
    else {
        const double wl2 = 1.0e-6 * wl * wl; // [nm] -> [um]
        const double A1 = self->coef[0];
        const double A2 = self->coef[1];
        const double A3 = self->coef[2];
        const double wl2_B1 = wl2 - self->coef[3];
        const double wl2_B2 = wl2 - self->coef[4];
        const double wl2_B3 = wl2 - self->coef[5];

        if (wl2_B1 && wl2_B2 && wl2_B3) {
            ind = sqrt(1.0 + A1 * wl2 / (wl2_B1) + A2 * wl2 / (wl2_B2) + A3 * wl2 / (wl2_B3));
            UpcDictDD_add(self->cash, wl, ind);
            return ind;
        }
    }

    return 1.0e30;
}

UpcGlass *UpcGlassSellmeier_init(const char *name, double A1, double A2, double A3, double B1, double B2, double B3)
{
    UpcGlass *self = (UpcGlass *)malloc(sizeof(UpcGlass));

    if (self) {
        self->_refCount = 1;
        self->_dealloc = UpcGlassSellmeier_dealloc;
        self->name = upc_copy_string(name);
        if (!(self->name)) {
            free(self);
            return NULL;
        }
        self->cash = UpcDictDD_init();
        if (!(self->cash)) {
            free(self->name);
            free(self);
            return NULL;
        }
        self->coef[0] = A1;
        self->coef[1] = A2;
        self->coef[2] = A3;
        self->coef[3] = B1;
        self->coef[4] = B2;
        self->coef[5] = B3;
        self->n = UpcGlassSellmeier_n;
        self->nwav = 0;
        self->indexF = NULL;
        self->indexB = NULL;
        self->toXml = UpcGlassSellmeier_toXml;
        self->lset = UpcGlass_lsetDefault;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcGlassSellmeier_init");
    }
    return self;
}

//------------------------------------------------------------
// UpcGlassNdVd

static const double WLd = 587.56e-3;
static const double WLF = 486.13e-3;
static const double WLC = 656.27e-3;
static const double WLg = 435.83e-3;

/* 以下の方程式をニュートン法で解く.
 * f(A1, A2) = sqrt(N(WLF)) - sqrt(N(WLC)) - (nd - 1) / vd = 0
 * N(WL) = nd^2 + A1 * (WL^2 - WLd^2) + A2 * (WL^-2 - WLd^-2)
 */
static Bool solveG(double a1, double nd2, double alpha, double theta_gF, int iter_max, double tol, double *a2, double *g)
{
    const double WLd2 = WLd * WLd;
    const double WLd2i = 1.0 / WLd2;
    const double WLF2_d2 = WLF * WLF - WLd2;
    const double WLC2_d2 = WLC * WLC - WLd2;
    const double WLg2_d2 = WLg * WLg - WLd2;
    const double WLF2i_d2i = 1.0 / (WLF * WLF) - WLd2i;
    const double WLC2i_d2i = 1.0 / (WLC * WLC) - WLd2i;
    const double WLg2i_d2i = 1.0 / (WLg * WLg) - WLd2i;
    double nF, nC, ng, f, df;
    int i;
    
    *a2 = 0.0;
    for (i = 0; i < iter_max; i++) {
        nF = sqrt(nd2 + a1 * WLF2_d2 + *a2 * WLF2i_d2i);
        nC = sqrt(nd2 + a1 * WLC2_d2 + *a2 * WLC2i_d2i);
        f = nF - nC - alpha;
        if (fabs(f) < tol)
            break;
        df = 0.5 * (WLF2i_d2i / nF - WLC2i_d2i / nC);
        *a2 -= f / df;
    }
    if (i == iter_max)
        return FALSE; // Newton法が収束しなかった.
    ng = sqrt(nd2 + a1 * WLg2_d2 + *a2 * WLg2i_d2i);
    nF = sqrt(nd2 + a1 * WLF2_d2 + *a2 * WLF2i_d2i);
    *g = ng - nF - alpha * theta_gF;
    return TRUE;
}

/* 与えられたNd, νdにFitするA0, A1, A2を探す.
 * n(λ) = sqrt(A0 + A1 * λ^2 + A2 * λ^-2)
 */
static Bool searchA0A1A2(double nd, double vd, double *A0, double *A1, double *A2)
{
    const double WLd2 = WLd * WLd;
    const double WLd2i = 1.0 / WLd2;
    const double theta_gF = -0.0392 / 24.23 * vd + 0.641462484;
    const double alpha = (nd - 1.0) / vd;
    const double nd2 = nd * nd;
    const double tol = 1.0e-14;
    const int iter_max = 20;
    double step = 0.001;
    double a1L = 0.0;
    double a1R = a1L + step;
    double a1, a2, g, gL, gR;
    int i;
    
    if (!solveG(a1L, nd2, alpha, theta_gF, iter_max, tol, &a2, &gL))
        return FALSE;
    if (!solveG(a1R, nd2, alpha, theta_gF, iter_max, tol, &a2, &gR))
        return FALSE;
    if (gL * gR > 0.0) {
        if (fabs(gL) > fabs(gR)) {
            for (i = 0; i < iter_max; i++) {
                a1L = a1R;
                gL = gR;
                step *= 2.0;
                a1R += step;
                if (!solveG(a1R, nd2, alpha, theta_gF, iter_max, tol, &a2, &gR))
                    return FALSE;
                if (gL * gR <= 0.0)
                    break;
            }
            if (i == iter_max)
                return FALSE; // 囲い込みに失敗
        }
        else {
            for (i = 0; i < iter_max; i++) {
                a1R = a1L;
                gR = gL;
                step *= 2.0;
                a1L -= step;
                if (!solveG(a1L, nd2, alpha, theta_gF, iter_max, tol, &a2, &gL))
                    return FALSE;
                if (gL * gR <= 0.0)
                    break;
            }
            if (i == iter_max)
                return FALSE; // 囲い込みに失敗
        }
    }
    if (gL >= 0.0) {
        while (fabs(gL - gR) >= tol) {
            a1 = 0.5 * (a1L + a1R);
            if (!solveG(a1, nd2, alpha, theta_gF, iter_max, tol, &a2, &g))
                return FALSE;
            if (g >= 0.0) {
                a1L = a1;
                gL = g;
            }
            else {
                a1R = a1;
                gR = g;
            }
        }
    }
    else {
        while (fabs(gL - gR) >= tol) {
            a1 = 0.5 * (a1L + a1R);
            if (!solveG(a1, nd2, alpha, theta_gF, iter_max, tol, &a2, &g))
                return FALSE;
            if (g >= 0.0) {
                a1R = a1;
                gR = g;
            }
            else {
                a1L = a1;
                gL = g;
            }
        }
    }
    *A1 = 0.5 * (a1L + a1R);
    if (!solveG(*A1, nd2, alpha, theta_gF, iter_max, tol, A2, &g))
        return FALSE;
    *A0 = nd2 - *A1 * WLd2 - *A2 * WLd2i;
    return TRUE;
}

static void UpcGlassNdVd_dealloc(UpcGlass *self)
{
    free(self->name);
    free(self->indexF); // indexBも同時に解放される.
    UpcObj_RELEASE(self->cash);
    free(self);
}

static void UpcGlassNdVd_toXml(const UpcGlass *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);

    fprintf(ostream, "%s<glass type=\"ndvd\">\n", space);
    fprintf(ostream, "%s  <name>%s</name>\n", space, self->name);
    fprintf(ostream, "%s  <nd>%23.15e</nd>\n", space, self->coef[0]);
    fprintf(ostream, "%s  <vd>%23.15e</vd>\n", space, self->coef[1]);
    fprintf(ostream, "%s</glass>\n", space);
}

/* 与えた波長における屈折率を返す.
 * (注) 波長はnm単位で与える.
 */
static double UpcGlassNdVd_n(UpcGlass *self, double wl)
{
    double ind;
    
    if (UpcDictDD_search(self->cash, wl, &ind))
        return ind;
    else {
        const double wl2 = 1.0e-6 * wl * wl; // [nm] -> [um]
        const double A0 = self->coef[2];
        const double A1 = self->coef[3];
        const double A2 = self->coef[4];
        
        if (wl2) {
            ind = sqrt(A0 + A1 * wl2 + A2 / wl2);
            UpcDictDD_add(self->cash, wl, ind);
            return ind;
        }
    }
    return 1.0e30;
}

UpcGlass *UpcGlassNdVd_init(const char *name, double nd, double vd)
{
    UpcGlass *self = (UpcGlass *)malloc(sizeof(UpcGlass));
    
    if (self) {
        self->_refCount = 1;
        self->_dealloc = UpcGlassNdVd_dealloc;
        self->name = upc_copy_string(name);
        if (!(self->name)) {
            free(self);
            return NULL;
        }
        self->cash = UpcDictDD_init();
        if (!(self->cash)) {
            free(self->name);
            free(self);
            return NULL;
        }
        self->coef[0] = nd;
        self->coef[1] = vd;
        if (!searchA0A1A2(nd, vd, self->coef + 2, self->coef + 3, self->coef + 4)) {
            self->cash->_dealloc(self->cash);
            free(self);
            return NULL;
        }
        self->n = UpcGlassNdVd_n;
        self->nwav = 0;
        self->indexF = NULL;
        self->indexB = NULL;
        self->toXml = UpcGlassNdVd_toXml;
        self->lset = UpcGlass_lsetDefault;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcGlassNdVd_init");
    }
    return self;
}

//------------------------------------------------------------
// UpcGlassConstant

static void UpcGlassConstant_dealloc(UpcGlass *self)
{
    free(self->name);
    free(self->indexF); // indexBも同時に解放される.
    free(self);
}

static void UpcGlassConstant_toXml(const UpcGlass *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);
    
    fprintf(ostream, "%s<glass type=\"constant\">\n", space);
    fprintf(ostream, "%s  <name>%s</name>\n", space, self->name);
    fprintf(ostream, "%s  <const>%23.15e</const>\n", space, self->coef[0]);
    fprintf(ostream, "%s</glass>\n", space);
}

/* 与えた波長における屈折率を返す.
 * (注) 波長はnm単位で与える.
 */
static double UpcGlassConstant_n(UpcGlass *self, double wl)
{
    return self->coef[0];
}

UpcGlass *UpcGlassConstant_init(const char *name, double ind)
{
    UpcGlass *self = (UpcGlass *)malloc(sizeof(UpcGlass));

    if (self) {
        self->_refCount = 1;
        self->_dealloc = UpcGlassConstant_dealloc;
        self->name = upc_copy_string(name);
        if (!(self->name)) {
            free(self);
            return NULL;
        }
        self->coef[0] = ind;
        self->n = UpcGlassConstant_n;
        self->nwav = 0;
        self->indexF = NULL;
        self->indexB = NULL;
        self->toXml = UpcGlassConstant_toXml;
        self->lset = UpcGlass_lsetDefault;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcGlassConstant_init");
    }
    return self;
}

//------------------------------------------------------------
// UpcGlassSpline

static void UpcGlassSpline_dealloc(UpcGlass *self)
{
    free(self->name);
    free(self->indexF); // indexBも同時に解放される.
    UpcObj_RELEASE(self->cash);
    UpcObj_RELEASE(self->splineObj);
    free(self);
}

static void UpcGlassSpline_toXml(const UpcGlass *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);
    int i, n = self->splineObj->ndata;

    fprintf(ostream, "%s<glass type=\"spline\" ndata=\"%d\">\n", space, n);
    fprintf(ostream, "%s  <name>%s</name>\n", space, self->name);
    for (i = 0; i < n; i++)
        fprintf(ostream, "%s  <wl i=\"%d\">%23.15e</wl>\n", space, i, self->splineObj->xa[i]);
    for (i = 0; i < n; i++)
        fprintf(ostream, "%s  <ind i=\"%d\">%23.15e</ind>\n", space, i, self->splineObj->ya[i]);
    fprintf(ostream, "%s</glass>\n", space);
}

/* 与えた波長における屈折率を返す.
 * (注) 波長はnm単位で与える.
 */
static double UpcGlassSpline_n(UpcGlass *self, double wl)
{
    double ind;
    
    if (UpcDictDD_search(self->cash, wl, &ind))
        return ind;
    else {
        ind = UpcSpline_y(self->splineObj, wl);
        UpcDictDD_add(self->cash, wl, ind);
        return ind;
    }
}

UpcGlass *UpcGlassSpline_init(const char *name, int ndata, const double *wls, const double *inds)
{
    if (ndata <= 0) {
        return NULL;
    }
    else if (ndata == 1) {
        return UpcGlassConstant_init(name, inds[0]);
    }
    else {
        UpcGlass *self = (UpcGlass *)malloc(sizeof(UpcGlass));
        int i;
        
        if (self) {
            self->_refCount = 1;
            self->_dealloc = UpcGlassSpline_dealloc;
            self->name = upc_copy_string(name);
            if (!(self->name)) {
                free(self);
                return NULL;
            }
            self->cash = UpcDictDD_init();
            if (!(self->cash)) {
                free(self->name);
                free(self);
                return NULL;
            }
            self->splineObj = UpcSpline_init(ndata, wls, inds, 1.0e30, 1.0e30);
            if (!(self->splineObj)) {
                UpcObj_RELEASE(self->cash);
                free(self);
                return NULL;
            }
            for (i = 0; i < ndata; i++)
                UpcDictDD_add(self->cash, wls[i], inds[i]);
            self->n = UpcGlassSpline_n;
            self->nwav = 0;
            self->indexF = NULL;
            self->indexB = NULL;
            self->toXml = UpcGlassSpline_toXml;
            self->lset = UpcGlass_lsetDefault;
        }
        else {
            UpcERRHANDLER(UpcE_MemoryError, "in UpcGlassSpline_init");
        }
        return self;
    }
}

//------------------------------------------------------------
// 共通メソッド

Bool UpcGlass_setIndex(UpcGlass *self, int nwav, const double *wls)
{
    int i;
    double n;

    if (self->nwav != nwav) {
        self->nwav = nwav;
        free(self->indexF);
        self->indexF = (double *)malloc(sizeof(double) * (nwav * 2)); // indexBの分も一緒に確保.
        if (!(self->indexF)) {
            free(self->indexF);
            UpcERRHANDLER(UpcE_MemoryError, "in UpcGlass_setIndex");
            return FALSE;
        }
        self->indexB = self->indexF + nwav; // indexBはindexFの後半部分を利用.
    }
    for (i = 0; i < nwav; i++) {
        n = self->n(self, wls[i]);
        self->indexF[i] = n;
        self->indexB[i] = -n;
    }
    return TRUE;
}



