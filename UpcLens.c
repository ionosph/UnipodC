/*
 *  UpcLens.c
 *  UnipodC
 *
 *  Created by ionosph on 2011/04/02.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcLens.h"
#include "UpcVector.h"
#include "UpcRayTrace.h"

static const double INF_OBJ = 1.0e10; // 無限遠とみなす物体距離
static const int isur_sto_default = 1;
static const int iwav_pri_default = 0;
static const Bool pim_default = 0;
static const enum UpcPrayMode pray_default = UpcPrayMode_stop;
static const enum UpcMrayMode mray_default = UpcMrayMode_eape;
static const double aperture_default = 0.0;
static const double tol_pray_default = 1.0e-12;
static const double tol_mray_default = 1.0e-10;
static const double wl_default = 587.56;
static const double wtw_default = 1.0;
static const double hgt_default = 0.0;
static const double wth_default = 1.0;
static const double vig_default = 1.0;


static void UpcLens_dealloc(UpcLens *self)
{
    int i, n = self->nsur + 2;
    
    UpcObj_RELEASE(self->gcatalog);
    for (i = 0; i < n; i++)
        UpcObj_RELEASE(self->s[i]);
    free(self->title);
    free(self->comment);
    free(self->wa);
    free(self->ha);
    free(self->s);
    free(self->parax);
    if (self->raytracef)
        UpcObj_RELEASE(self->raytracef);
    free(self);
}

static void deleteUpcGlass(void *obj)
{
    UpcGlass *t = (UpcGlass *)obj;

    UpcObj_RELEASE(t);
}

UpcLens *UpcLens_init(int nsur, int nwav, int nhgt)
{
    UpcLens *self;
    int i, j, n = nsur + 2;
    UpcGlass *air;

    if (nsur < 1 || nwav < 1 || nhgt < 1) {
        UpcERRHANDLER(UpcE_ValueError, "nsur < 1 or nwav < 1 or nhgt < 1 (in UpcLens_lset)");
        return NULL;
    }

    self = (UpcLens *)calloc(1, sizeof(UpcLens));
    if (self) {
        self->_refCount = 1;
        self->_dealloc = UpcLens_dealloc;
        self->nsur = nsur;
        self->nwav = nwav;
        self->nhgt = nhgt;
        self->wa = (double *)malloc(sizeof(double) * (nwav * 2));
        self->ha = (double *)malloc(sizeof(double) * (nhgt * 7));
        self->s = (UpcSurf **)calloc(n, sizeof(UpcSurf *));
        self->parax = (UpcParaxData *)calloc(nwav, sizeof(UpcParaxData));
        self->gcatalog = UpcDict_init(upc_comp_string, free, deleteUpcGlass);
        air = UpcGlassConstant_init("air", 1.0);
        if (!(self->wa) || !(self->ha) || !(self->s) || !(self->parax) || !(self->gcatalog) || !air) {
            goto err_return;
        }
        for (i = 0; i < n; i++) {
            self->s[i] = UpcSurf_init();
            if (!(self->s[i])) {
                for (j = 0; j < i; j++)
                    self->s[i]->_dealloc(self->s[i]);
                goto err_return;
            }
            self->s[i]->medium = air; // 弱参照 (RETAINしない)
        }
        UpcDict_add(self->gcatalog, upc_copy_string(air->name), air); // 所有権はgcatalogに移る
        for (i = 0; i < nwav; i++) {
            UpcLens_WL(self, i) = wl_default;
            UpcLens_WTW(self, i) = wtw_default;
        }
        for (i = 0; i < nhgt; i++) {
            UpcLens_HGTX(self, i) = hgt_default;
            UpcLens_HGTY(self, i) = hgt_default;
            UpcLens_WTH(self, i) = wth_default;
            UpcLens_VIG(self, i, YUPP) = vig_default;
            UpcLens_VIG(self, i, YLOW) = vig_default;
            UpcLens_VIG(self, i, XUPP) = vig_default;
            UpcLens_VIG(self, i, XLOW) = vig_default;
        }
        self->title = NULL;
        self->comment = NULL;
        self->isur_sto = isur_sto_default;
        self->iwav_pri = iwav_pri_default;
        self->pim = pim_default;
        self->pray = pray_default;
        self->mray = mray_default;
        self->aperture = aperture_default;
        self->tol_pray = tol_pray_default;
        self->tol_mray = tol_mray_default;
        self->raytracef = NULL;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcLens_init");
    }
    return self;

err_return:
    free(self->wa);
    free(self->ha);
    free(self->s);
    free(self->parax);
    if (self->gcatalog)
        self->gcatalog->_dealloc(self->gcatalog);
    if (air)
        air->_dealloc(air);
    free(self);
    UpcERRHANDLER(UpcE_MemoryError, "in UpcLens_init");
    return NULL;
}

UpcLens *UpcLens_copy(UpcLens *self)
{
    UpcLens *t = (UpcLens *)calloc(1, sizeof(UpcLens));

    if (t) {
        const int n = self->nsur + 2;
        const int nwav = self->nwav;
        const int nhgt = self->nhgt;
        int i, j;

        t->_refCount = 1;
        t->_dealloc = UpcLens_dealloc;
        t->wa = (double *)malloc(sizeof(double) * (nwav * 2));
        t->ha = (double *)malloc(sizeof(double) * (nhgt * 7));
        t->s = (UpcSurf **)calloc(n, sizeof(UpcSurf *));
        t->parax = (UpcParaxData *)calloc(nwav, sizeof(UpcParaxData));
        if (!(t->wa) || !(t->ha) || !(t->s) || !(t->parax)) {
            goto err_return;
        }
        for (i = 0; i < n; i++) {
            t->s[i] = UpcSurf_copy(self->s[i]);
            if (!(t->s[i])) {
                for (j = 0; j < i; j++)
                    t->s[i]->_dealloc(t->s[i]);
                goto err_return;
            }
        }
        t->gcatalog = self->gcatalog;
        UpcObj_RETAIN(t->gcatalog);
        t->nsur = self->nsur;
        t->nwav = nwav;
        t->nhgt = nhgt;
        if (self->title)
            t->title = upc_copy_string(self->title);
        if (self->comment)
            t->comment = upc_copy_string(self->comment);
        t->isur_sto = self->isur_sto;
        t->iwav_pri = self->iwav_pri;
        t->pim = self->pim;
        t->pray = self->pray;
        t->mray = self->mray;
        t->aperture = self->aperture;
        t->tol_pray = self->tol_pray;
        t->tol_mray = self->tol_mray;
        t->xSymmetric = self->xSymmetric;
        t->ySymmetric = self->ySymmetric;
        t->rSymmetric = self->rSymmetric;
        t->infObjD = self->infObjD;
        for (i = 0; i < nwav; i++) {
            UpcLens_WL(t, i) = UpcLens_WL(self, i);
            UpcLens_WTW(t, i) = UpcLens_WTW(self, i);
            t->parax[i] = self->parax[i];
        }
        for (i = 0; i < nhgt; i++) {
            UpcLens_HGTX(t, i) = UpcLens_HGTX(self, i);
            UpcLens_HGTY(t, i) = UpcLens_HGTY(self, i);
            UpcLens_WTH(t, i) = UpcLens_WTH(self, i);
            UpcLens_VIG(t, i, YUPP) = UpcLens_VIG(self, i, YUPP);
            UpcLens_VIG(t, i, YLOW) = UpcLens_VIG(self, i, YLOW);
            UpcLens_VIG(t, i, XUPP) = UpcLens_VIG(self, i, XUPP);
            UpcLens_VIG(t, i, XLOW) = UpcLens_VIG(self, i, XLOW);
        }
        t->raytracef = self->raytracef;
        if (t->raytracef)
            UpcObj_RETAIN(t->raytracef);
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcLens_copy");
    }
    return t;

err_return:
    free(t->wa);
    free(t->ha);
    free(t->s);
    free(t->parax);
    free(t);
    UpcERRHANDLER(UpcE_MemoryError, "in UpcLens_copy");
    return NULL;
}

char *UpcPrayMode_toString(enum UpcPrayMode i)
{
    switch (i) {
        case UpcPrayMode_stop:
            return "stop";
        case UpcPrayMode_itel:
            return "itel";
        case UpcPrayMode_otel:
            return "otel";
    }
    return "unknown mode";
}

enum UpcPrayMode UpcPrayMode_fromString(const char *s)
{
    if (!strcmp(s, "stop"))
        return UpcPrayMode_stop;
    else if (!strcmp(s, "itel"))
        return UpcPrayMode_itel;
    else if (!strcmp(s, "otel"))
        return UpcPrayMode_otel;
    else {
        UpcERRHANDLER(UpcE_ValueError, "unknown pray mode (in UpcPrayMode_fromString)");
        return pray_default;
    }
}

char *UpcMrayMode_toString(enum UpcMrayMode i)
{
    switch (i) {
        case UpcMrayMode_eape:
            return "eape";
        case UpcMrayMode_isin:
            return "isin";
        case UpcMrayMode_osin:
            return "osin";
    }
    return "unknown mode";
}

enum UpcMrayMode UpcMrayMode_fromString(const char *s)
{
    if (!strcmp(s, "eape"))
        return UpcMrayMode_eape;
    else if (!strcmp(s, "isin"))
        return UpcMrayMode_isin;
    else if (!strcmp(s, "osin"))
        return UpcMrayMode_osin;
    else {
        UpcERRHANDLER(UpcE_ValueError, "unknown mray mode (in UpcMrayMode_fromString)");
        return mray_default;
    }
}

void UpcLens_print(const UpcLens *self, FILE *ostream)
{
    int isur, iwav, ihgt, isur_img = self->nsur + 1;

    fprintf(ostream, "[UpcLens]\n");
    if (self->title)
        fprintf(ostream, "  title     = \"%s\"\n", self->title);
    else
        fprintf(ostream, "  title     = NULL\n");
    if (self->comment)
        fprintf(ostream, "  comment   = \"%s\"\n", self->comment);
    else
        fprintf(ostream, "  comment   = NULL\n");
    fprintf(ostream, "  nsur      = %2d\n", self->nsur);
    fprintf(ostream, "  nwav      = %2d\n", self->nwav);
    fprintf(ostream, "  nhgt      = %2d\n", self->nhgt);
    fprintf(ostream, "  iwav_pri  = %2d\n", self->iwav_pri);
    fprintf(ostream, "  isur_sto  = %2d\n", self->isur_sto);
    fprintf(ostream, "  pim       = %s\n", (self->pim) ? "True" : "False");
    for (iwav = 0; iwav < self->nwav; iwav++)
        fprintf(ostream, "  wl [%2d]   = %13.6e\n", iwav, UpcLens_WL(self, iwav));
    for (iwav = 0; iwav < self->nwav; iwav++)
        fprintf(ostream, "  wtw[%2d]   = %13.6e\n", iwav, UpcLens_WTW(self, iwav));
    for (ihgt = 0; ihgt < self->nhgt; ihgt++)
        fprintf(ostream, "  hgt[%2d]   = %13.6e, %13.6e\n", ihgt, UpcLens_HGTX(self, ihgt), UpcLens_HGTY(self, ihgt));
    for (ihgt = 0; ihgt < self->nhgt; ihgt++)
        fprintf(ostream, "  wth[%2d]   = %13.6e\n", ihgt, UpcLens_WTH(self, ihgt));
    for (ihgt = 0; ihgt < self->nhgt; ihgt++)
        fprintf(ostream, "  vig[%2d]   = %13.6e, %13.6e, %13.6e, %13.6e\n", ihgt, UpcLens_VIG(self, ihgt, YUPP), UpcLens_VIG(self, ihgt, YLOW), UpcLens_VIG(self, ihgt, XUPP), UpcLens_VIG(self, ihgt, XLOW));
    fprintf(ostream, "  pray      = %s\n", UpcPrayMode_toString(self->pray));
    fprintf(ostream, "  mray      = %s\n", UpcMrayMode_toString(self->mray));
    fprintf(ostream, "  aperture  = %13.6e\n", self->aperture);
    fprintf(ostream, "  tol_pray  = %13.6e\n", self->tol_pray);
    fprintf(ostream, "  tol_mray  = %13.6e\n", self->tol_mray);
    for (isur = 0; isur <= isur_img; isur++) {
        fprintf(ostream, "  surf[%d]\n", isur);
        UpcSurf_print(self->s[isur], isur, 4, ostream);
    }
    fprintf(ostream, "[calculated values]\n");
    fprintf(ostream, "  xSymmetric  = %s\n", self->xSymmetric ? "True" : "False");
    fprintf(ostream, "  ySymmetric  = %s\n", self->ySymmetric ? "True" : "False");
    fprintf(ostream, "  rSymmetric  = %s\n", self->rSymmetric ? "True" : "False");
    fprintf(ostream, "  infObjD     = %s\n", self->infObjD ? "True" : "False");
    for (iwav = 0; iwav < self->nwav; iwav++) {
        fprintf(ostream, "[parax(iwav=%d)]\n", iwav);
        fprintf(ostream, "    fl        = %13.6e\n", self->parax[iwav].fl);
        fprintf(ostream, "    mag       = %13.6e\n", self->parax[iwav].mag);
        fprintf(ostream, "    sk        = %13.6e\n", self->parax[iwav].sk);
        fprintf(ostream, "    enp       = %13.6e\n", self->parax[iwav].enp);
        fprintf(ostream, "    exp       = %13.6e\n", self->parax[iwav].exp);
        fprintf(ostream, "    o1        = %13.6e\n", self->parax[iwav].o1);
        fprintf(ostream, "    ok        = %13.6e\n", self->parax[iwav].ok);
    }
}

void UpcLens_toXml(const UpcLens *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);
    UpcGlass **glasses = (UpcGlass **)UpcDict_vals(self->gcatalog);
    int i, isur_img = self->nsur + 1;

    if (!(glasses)) {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcLens_toXml");
        return;
    }
    fprintf(ostream, "%s<lens nsur=\"%d\" nwav=\"%d\" nhgt=\"%d\">\n", space, self->nsur, self->nwav, self->nhgt);
    if (self->title)
        fprintf(ostream, "%s  <title>%s</title>\n", space, self->title);
    if (self->comment)
        fprintf(ostream, "%s  <comment>%s</comment>\n", space, self->comment);
    for (i = 0; i < self->nwav; i++)
        fprintf(ostream, "%s  <wl iwav=\"%d\">%23.15e</wl>\n", space, i, UpcLens_WL(self, i));
    for (i = 0; i < self->nwav; i++) {
        if (UpcLens_WTW(self, i) != wtw_default)
            fprintf(ostream, "%s  <wtw iwav=\"%d\">%23.15e</wtw>\n", space, i, UpcLens_WTW(self, i));
    }
    fprintf(ostream, "%s  <iwav_pri>%d</iwav_pri>\n", space, self->iwav_pri);
    fprintf(ostream, "%s  <pray>%s</pray>\n", space, UpcPrayMode_toString(self->pray));
    fprintf(ostream, "%s  <mray>%s</mray>\n", space, UpcMrayMode_toString(self->mray));
    if (self->pray == UpcPrayMode_stop)
        fprintf(ostream, "%s  <isur_sto>%d</isur_sto>\n", space, self->isur_sto);
    if (self->mray != UpcMrayMode_eape)
        fprintf(ostream, "%s  <aperture>%23.15e</aperture>\n", space, self->aperture);
    if (self->pim != pim_default)
        fprintf(ostream, "%s  <pim>%s</pim>\n", space, BOOLTOSTRING(self->pim));
    if (self->tol_pray != tol_pray_default)
        fprintf(ostream, "%s  <tol_pray>%23.15e</tol_pray>\n", space, self->tol_pray);
    if (self->tol_mray != tol_mray_default)
        fprintf(ostream, "%s  <tol_mray>%23.15e</tol_mray>\n", space, self->tol_mray);
    for (i = 0; i < self->nhgt; i++) {
        if (UpcLens_HGTX(self, i) != hgt_default)
            fprintf(ostream, "%s  <hgtx ihgt=\"%d\">%23.15e</hgtx>\n", space, i, UpcLens_HGTX(self, i));
        if (UpcLens_HGTY(self, i) != hgt_default)
            fprintf(ostream, "%s  <hgty ihgt=\"%d\">%23.15e</hgty>\n", space, i, UpcLens_HGTY(self, i));
    }
    for (i = 0; i < self->nhgt; i++) {
        if (UpcLens_WTH(self, i) != wth_default)
            fprintf(ostream, "%s  <wth ihgt=\"%d\">%23.15e</wth>\n", space, i, UpcLens_WTH(self, i));
    }
    if (self->mray != UpcMrayMode_eape) {
        for (i = 0; i < self->nhgt; i++) {
            if (UpcLens_VIG(self, i, YUPP) != vig_default)
                fprintf(ostream, "%s  <vig ihgt=\"%d\" ray=\"yupp\">%23.15e</vig>\n", space, i, UpcLens_VIG(self, i, YUPP));
            if (UpcLens_VIG(self, i, YLOW) != vig_default)
                fprintf(ostream, "%s  <vig ihgt=\"%d\" ray=\"ylow\">%23.15e</vig>\n", space, i, UpcLens_VIG(self, i, YLOW));
            if (UpcLens_VIG(self, i, XUPP) != vig_default)
                fprintf(ostream, "%s  <vig ihgt=\"%d\" ray=\"xupp\">%23.15e</vig>\n", space, i, UpcLens_VIG(self, i, XUPP));
            if (UpcLens_VIG(self, i, XLOW) != vig_default)
                fprintf(ostream, "%s  <vig ihgt=\"%d\" ray=\"xlow\">%23.15e</vig>\n", space, i, UpcLens_VIG(self, i, XLOW));
        }
    }
    fprintf(ostream, "%s  <gcatalog>\n", space);
    for (i = 0; glasses[i]; i++)
        glasses[i]->toXml(glasses[i], indent + 4, ostream);
    fprintf(ostream, "%s  </gcatalog>\n", space);
    for (i = 0; i <= isur_img; i++) {
        UpcSurf_toXml(self->s[i], i, indent + 2, ostream);
    }
    fprintf(ostream, "%s</lens>\n", space);
}

UpcErrStatus UpcLens_saveXml(UpcLens *self, char *filename)
{
    FILE *ostream = fopen(filename, "w");

    if (!ostream) {
        UpcERRHANDLER(UpcE_IOError, "can't file open (in UpcLens_save)");
        return UpcE_IOError;
    }
    fprintf(ostream, "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n");
    fprintf(ostream, "<!-- lens data for unipodc(ver 2.0) -->\n");
    UpcLens_toXml(self, 0, ostream);
    fclose(ostream);
    return UpcE_NoIssues;
}

void UpcLens_setTitle(UpcLens *self, const char *title)
{
    free(self->title);
    self->title = upc_copy_string(title);
}

void UpcLens_setComment(UpcLens *self, const char *comment)
{
    free(self->comment);
    self->comment = upc_copy_string(comment);
}

/* 波長数を設定する
 */
Bool UpcLens_setNwav(UpcLens *self, int nwav_new)
{
    int nwav_old = self->nwav;

    if (nwav_new < 1) {
        UpcERRHANDLER(UpcE_ValueError, "in UpcLens_setNwav");
        return FALSE;
    }

    if (nwav_new <= nwav_old) {
        self->nwav = nwav_new;
        return TRUE;
    }
    else {
        int i;
        double *wa_old = self->wa;
        double *wa_new = (double *)malloc(sizeof(double) * (nwav_new * 2));
        UpcParaxData *parax_new = (UpcParaxData *)malloc(sizeof(UpcParaxData) * nwav_new);
        
        if (!wa_new || !parax_new) {
            free(wa_new);
            free(parax_new);
            UpcERRHANDLER(UpcE_MemoryError, "in UpcLens_setNwav");
            return FALSE;
        }
        memcpy(wa_new, wa_old, sizeof(double) * (nwav_old * 2));
        /*
        for (i = 0; i < nwav_old; i++) {
            wa_new[2 * i] = wa_old[2 * i];
            wa_new[2 * i + 1] = wa_old[2 * i + 1];
        }
         */
        for (i = nwav_old; i < nwav_new; i++) {
            wa_new[2 * i] = wl_default;
            wa_new[2 * i + 1] = wtw_default;
        }
        free(wa_old);
        self->wa = wa_new;
        self->nwav = nwav_new;
        free(self->parax);
        self->parax = parax_new;
        return TRUE;
    }
}

/* 画角数を設定する
 */
Bool UpcLens_setNhgt(UpcLens *self, int nhgt_new)
{
    int nhgt_old = self->nhgt;
    
    if (nhgt_new < 1) {
        UpcERRHANDLER(UpcE_ValueError, "in UpcLens_setNhgt");
        return FALSE;
    }
    
    if (nhgt_new <= nhgt_old) {
        self->nhgt = nhgt_new;
        return TRUE;
    }
    else {
        int i;
        double *ha_old = self->ha;
        double *ha_new = (double *)malloc(sizeof(double) * (nhgt_new * 7));
        
        if (!ha_new) {
            UpcERRHANDLER(UpcE_MemoryError, "in UpcLens_setNhgt");
            return FALSE;
        }
        memcpy(ha_new, ha_old, sizeof(double) * (nhgt_old * 7));
        for (i = nhgt_old; i < nhgt_new; i++) {
            ha_new[7 * i + 0] = hgt_default;
            ha_new[7 * i + 1] = hgt_default;
            ha_new[7 * i + 2] = wth_default;
            ha_new[7 * i + 3] = vig_default;
            ha_new[7 * i + 4] = vig_default;
            ha_new[7 * i + 5] = vig_default;
            ha_new[7 * i + 6] = vig_default;
        }
        free(ha_old);
        self->ha = ha_new;
        self->nhgt = nhgt_new;
        return TRUE;
    }
}

/* 第isur_insert面の前に新しい面を挿入する
 */
Bool UpcLens_insertSurf(UpcLens *self, int isur_insert)
{
    int i, nsur = self->nsur;
    UpcSurf **newSa;
    UpcSurf **oldSa = self->s;
    UpcGlass *air;
    UpcSurf *newS;
    
    if (isur_insert < 1 || isur_insert > nsur + 1) {
        UpcERRHANDLER(UpcE_ValueError, "in UpcLens_insertSurf");
        return FALSE;
    }
    if (!UpcDict_search(self->gcatalog, "air", (void **)&air)) {
        UpcERRHANDLER(UpcE_RuntimeError, "'air' not found(in UpcLens_insertSurf)");
        return FALSE;
    }

    newSa = (UpcSurf **)malloc(sizeof(UpcSurf *) * (nsur + 2 + 1));
    if (!newSa) {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcLens_insertSurf");
        return FALSE;
    }
    newS = UpcSurf_init();
    if (!newS) {
        free(newSa);
        UpcERRHANDLER(UpcE_MemoryError, "in UpcLens_insertSurf");
        return FALSE;
    }
    newS->medium = air;

    for (i = 0; i < isur_insert; i++) {
        newSa[i] = oldSa[i];
    }
    newSa[isur_insert] = newS;
    for (i = isur_insert; i <= nsur + 1; i++) {
        newSa[i + 1] = oldSa[i];
    }
    free(oldSa);
    self->s = newSa;
    self->nsur += 1;
    return TRUE;
}

/* 第isur_delete面を削除する
 */
Bool UpcLens_deleteSurf(UpcLens *self, int isur_delete)
{
    int i, nsur = self->nsur;
    UpcSurf **sa = self->s;
    UpcSurf *sur;

    if (isur_delete < 1 || isur_delete > nsur) {
        UpcERRHANDLER(UpcE_ValueError, "in UxLens_deleteSur");
        return FALSE;
    }

    sur = sa[isur_delete];
    UpcObj_RELEASE(sur);
    for (i = isur_delete; i <= nsur; i++) {
        sa[i] = sa[i + 1];
    }
    self->nsur -= 1;
    return TRUE;
}

/* 面のパワーを返す.
*/
static double UpcLens_phi(const UpcLens *self, int isur, int iwav)
{
    UpcSurf *si = self->s[isur - 1];
    UpcSurf *sj = self->s[isur];

    return (sj->_n[iwav] - si->_n[iwav]) * sj->shape->axialPower(sj->shape);
}


static void UpcLens_setGlobalCoord(const UpcLens *self);


/* データの整合性をチェックし, 近軸量を算出する.
 */
UpcErrStatus UpcLens_lset(UpcLens *self)
{
    int isur, iwav;
    const int nwav = self->nwav;
    const int nsur = self->nsur;
    const int isur_img = nsur + 1;
    int direc;
    Bool xSymmetric = TRUE;
    Bool ySymmetric = TRUE;
    Bool rSymmetric = TRUE;
    UpcGlass **glss;
    double *wls;
    UpcSurf **sur = self->s;
    char err_msg[100];
    UpcErrStatus stts;
//  void UpcLens_setGlobalCoord(const UpcLens *self);

    // 物体面
    if (!(sur[0])) {
        UpcERRHANDLER(UpcE_ValueError, "s[0] is NULL");
        return UpcE_ValueError;
    }
    if (sur[0]->d >= INF_OBJ) {
        sur[0]->d = INF_OBJ;
        self->infObjD = TRUE;
    }
    else {
        self->infObjD = FALSE;
    }
    sur[0]->mode = RmodeRefractive;
    if (!(sur[0]->dec)) {
        UpcERRHANDLER(UpcE_ValueError, "s[0]->dec is NULL");
        return UpcE_ValueError;
    }
    UpcDecenter_setAllZero(sur[0]->dec);
    if ((stts = UpcSurf_lset(sur[0])) != UpcE_NoIssues) {
        UpcERRHANDLER(stts, "in s[0]");
        return stts;
    }
    // 各属性間の整合性チェック
    if (self->iwav_pri < 0 || self->iwav_pri >= nwav) {
        UpcERRHANDLER(UpcE_ValueError, "'iwav_pri' is out of range.");
        return UpcE_ValueError;
    }
    if (self->pray == UpcPrayMode_stop && (self->isur_sto <= 0 || self->isur_sto > nsur)) {
        UpcERRHANDLER(UpcE_ValueError, "'isur_sto' is out of range.");
        return UpcE_ValueError;
    }
    if (self->mray != UpcMrayMode_eape) {
        if (self->aperture >= 1.0 || !(self->aperture)) {
            UpcERRHANDLER(UpcE_ValueError, "'aperture' value is wrong.");
            return UpcE_ValueError;
        }
        if (self->mray == UpcMrayMode_osin && self->infObjD) {
            UpcERRHANDLER(UpcE_ValueError, "'mray = osin' conflicts with infinite object distance.");
            return UpcE_ValueError;
        }
    }
    if (self->tol_pray <= 0.0) {
        UpcERRHANDLER(UpcE_ValueError, "'tol_pray' must be positive value.");
        return UpcE_ValueError;
    }
    if (self->tol_mray <= 0.0) {
        UpcERRHANDLER(UpcE_ValueError, "'tol_mray' must be positive value.");
        return UpcE_ValueError;
    }
    for (iwav = 0; iwav < nwav; iwav++) {
        if (!UpcLens_WL(self, iwav)) {
            sprintf(err_msg, "'wl[%d]' is zero.", iwav);
            UpcERRHANDLER(UpcE_ValueError, err_msg);
            return UpcE_ValueError;
        }
    }
    // 各面
    for (isur = 1; isur <= nsur; isur++) {
        if (!(sur[isur])) {
            sprintf(err_msg, "s[%d] is NULL", isur);
            UpcERRHANDLER(UpcE_ValueError, err_msg);
            return UpcE_ValueError;
        }
        if ((stts = UpcSurf_lset(sur[isur])) != UpcE_NoIssues) {
            sprintf(err_msg, "in s[%d]", isur);
            UpcERRHANDLER(stts, err_msg);
            return stts;
        }
    }
    for (isur = 1; isur <= nsur; isur++) {
        if (sur[isur]->mode == RmodeReflective) {
            sur[isur]->medium = sur[isur - 1]->medium;
        }
    }
    // 像面
    if (!(sur[isur_img])) {
        UpcERRHANDLER(UpcE_ValueError, "s[isur_img] is NULL");
        return UpcE_ValueError;
    }
    sur[isur_img]->medium = sur[nsur]->medium;
    sur[isur_img]->mode = RmodeRefractive;
    sur[isur_img]->d = 0.0;
    if ((stts = UpcSurf_lset(sur[isur_img])) != UpcE_NoIssues) {
        UpcERRHANDLER(stts, "in s[isur_img]");
        return stts;
    }
    // 屈折率算出
    glss = (UpcGlass **)UpcDict_vals(self->gcatalog);
    wls = (double *)malloc(sizeof(double) * self->nwav);
    if (glss && wls) {
        int i;

        for (iwav = 0; iwav < nwav; iwav++)
            wls[iwav] = UpcLens_WL(self, iwav);
        for (i = 0; glss[i]; i++) {
            if (!UpcGlass_setIndex(glss[i], nwav, wls)) {
                free(glss);
                free(wls);
                UpcERRHANDLER(UpcE_MemoryError, "setIndex memory error (in UpcLens_lset).");
                return UpcE_MemoryError;
            }
        }
        free(glss);
        free(wls);
    }
    else {
        free(glss);
        free(wls);
        UpcERRHANDLER(UpcE_MemoryError, "in UpcLens_lset");
        return UpcE_MemoryError;
    }
    direc = 1;
    for (isur = 0; isur <= isur_img; isur++) {
        if (sur[isur]->mode == RmodeReflective)
            direc *= -1;
        sur[isur]->_n = (direc == 1) ? sur[isur]->medium->indexF : sur[isur]->medium->indexB;
    }
    for (iwav = 0; iwav < nwav; iwav++) {
        for (isur = 0; isur <= nsur; isur++) {
            if (!(sur[isur]->_n[iwav])) {
                sprintf(err_msg, "'s[%d]->n[%d]' is zero.", isur, iwav);
                UpcERRHANDLER(UpcE_ValueError, err_msg);
                return UpcE_ValueError;
            }
        }
        UpcLens_getParax(self, iwav);
    }
    // 対称性の確認
    for (isur = 0; isur <= nsur; isur++) {
        if (xSymmetric && !UpcSurf_xSymmetric(sur[isur]))
            xSymmetric = FALSE;
        if (ySymmetric && !UpcSurf_ySymmetric(sur[isur]))
            ySymmetric = FALSE;
        if (rSymmetric && !UpcSurf_rSymmetric(sur[isur]))
            rSymmetric = FALSE;
    }
    self->xSymmetric = xSymmetric;
    self->ySymmetric = ySymmetric;
    self->rSymmetric = rSymmetric;
    // 像面位置調整
    if (self->pim)
        sur[nsur]->d = self->parax[self->iwav_pri].sk;
    // グローバル座標の構築
    UpcLens_setGlobalCoord(self);
    // 光線追跡オブジェクトの構築
    if (self->raytracef)
        UpcObj_RELEASE(self->raytracef);
    self->raytracef = NULL;
    self->raytracef = UpcTraceFactory_init(self);
    if (!(self->raytracef)) {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcLens_lset");
        return UpcE_MemoryError;
    }
    return UpcE_NoIssues;
}

/* isur_stt面からisur_end面までの近軸行列を計算する.
 [入力]
 self : Lensオブジェクト
 isur_stt : 開始面 (0にすると第1面から計算)
 isur_end : 終了面 (0にすると最終面まで計算)
 iwav : 波長番号
 [出力]
 m[4] : c11, c12, c21, c22 
 [参考] 近軸行列 m の使い方
 パワーφ = m.c21
 前側主点位置Δ  = (1.0 - c22) / φ
 後側主点位置Δ' = (c11 - 1.0) / φ
 ｈk  = c11 * ｈ1 + c12 * α1
 αk' = c21 * ｈ1 + c22 * α1
 */
void UpcLens_paraxMatrix(const UpcLens *self, int isur_stt, int isur_end, int iwav, double *m)
{
    int i;
    double c12, c21, a[4], b[4];

    if (!isur_stt)
        isur_stt = 1;
    if (!isur_end)
        isur_end = self->nsur;

    m[0] = 1.0;
    m[1] = 0.0;
    m[2] = UpcLens_phi(self, isur_stt, iwav);
    m[3] = 1.0;
    for (i = isur_stt + 1; i <= isur_end; i++) {
        c12 = - UpcSurf_e(self->s[i - 1], iwav);
        c21 = UpcLens_phi(self, i, iwav);
        a[0] = 1.0;
        a[1] = c12;
        a[2] = c21;
        a[3] = c12 * c21 + 1.0;
        b[0] = a[0] * m[0] + a[1] * m[2];
        b[1] = a[0] * m[1] + a[1] * m[3];
        b[2] = a[2] * m[0] + a[3] * m[2];
        b[3] = a[2] * m[1] + a[3] * m[3];
        m[0] = b[0];
        m[1] = b[1];
        m[2] = b[2];
        m[3] = b[3];
    }
}

/* 近軸諸量の計算
 * [入力]
 *     self : Lensオブジェクト
 *     iwav : 波長番号
 * [出力]
 *     self->parax属性に計算結果が格納される.
 *     戻り値 : 計算できればTRUE, 失敗したらFALSE
 */
Bool UpcLens_getParax(UpcLens *self, int iwav)
{
    double m[4], m2[4];
    double h1, a1, hk, ak, e1, ek;
    const int istop = self->isur_sto;
    const int nsur = self->nsur;
    
    if (iwav < 0 || iwav >= self->nwav)
        return FALSE;
    UpcLens_paraxMatrix(self, 1, nsur, iwav, m);
    if (m[2]) {
        self->parax[iwav].fl = 1.0 / m[2];
        self->parax[iwav].o1 = self->s[0]->_n[iwav] * ((1.0 - m[3]) / m[2]);
        self->parax[iwav].ok = self->s[nsur]->_n[iwav] * ((m[0] - 1.0) / m[2]);
    }
    else {
        self->parax[iwav].fl = 1.0e30;
        self->parax[iwav].o1 = 1.0e30;
        self->parax[iwav].ok = 1.0e30;
    }
    if (self->infObjD) {
        h1 = 1.0;
        a1 = 0.0;
    }
    else {
        h1 = UpcSurf_e(self->s[0], iwav);
        a1 = -1.0;
    }
    hk = m[0] * h1 + m[1] * a1;
    ak = m[2] * h1 + m[3] * a1;
    if (ak) {
        self->parax[iwav].mag = a1 / ak;
        self->parax[iwav].sk  = hk / ak * self->s[nsur]->_n[iwav];
    }
    else {
        self->parax[iwav].mag = 1.0e30;
        self->parax[iwav].sk  = 1.0e30;
    }
    // 入射瞳位置の計算
    if (self->pray == UpcPrayMode_otel) {
        self->parax[iwav].enp = 1.0e6;
    }
    else {
        if (self->pray == UpcPrayMode_itel) {
            m2[0] = m[0];
            m2[1] = m[1];
            m2[2] = m[2];
            m2[3] = m[3];
            UpcMatrix2_inverse(m2);
            hk = 1.0;
            ak = 0.0;
        }
        else {
            if (istop == 1) {
                self->parax[iwav].enp = 0.0;
                goto cal_exp;
            }
            UpcLens_paraxMatrix(self, 1, istop - 1, iwav, m2);
            UpcMatrix2_inverse(m2);
            ek = UpcSurf_e(self->s[istop - 1], iwav);
            ak = 1.0;
            hk = ek * ak;
        }
        h1 = m2[0] * hk + m2[1] * ak;
        a1 = m2[2] * hk + m2[3] * ak;
        if (fabs(a1) < 1.0e-6)
            self->parax[iwav].enp = 1.0e6;
        else
            self->parax[iwav].enp = self->s[0]->_n[iwav] * h1 / a1;
    }
    // 射出瞳位置の計算
cal_exp:
    if (self->pray == UpcPrayMode_itel) {
        self->parax[iwav].exp = -1.0e6;
    }
    else {
        if (self->pray == UpcPrayMode_otel) {
            m2[0] = m[0];
            m2[1] = m[1];
            m2[2] = m[2];
            m2[3] = m[3];
            h1 = 1.0;
            a1 = 0.0;
        }
        else {
            if (istop == nsur) {
                self->parax[iwav].exp = 0.0;
                return TRUE;
            }
            UpcLens_paraxMatrix(self, istop + 1, nsur, iwav, m2);
            e1 = UpcSurf_e(self->s[istop], iwav);
            a1 = 1.0;
            h1 = - a1 * e1;
        }
        hk = m2[0] * h1 + m2[1] * a1;
        ak = m2[2] * h1 + m2[3] * a1;
        if (fabs(ak) < 1.0e-6)
            self->parax[iwav].exp = -1.0e6;
        else
            self->parax[iwav].exp = self->s[0]->_n[iwav] * hk / ak;
    }
    return TRUE;
}

/* 近軸光線追跡を行う.(物体面から)
 * [入力]
 *     self : UpcLensオブジェクト
 *     h0   : 物体面における光線高さ
 *     a0d  : 物体面を射出する光線の換算傾角
 *     iwav ; 波長番号
 * [出力]
 *     hv[0..nsur + 1]  : 第v面の光線高さ
 *     avd[0..nsur + 1] : 第v面を射出する光線の換算傾角
 */
Bool UpcLens_paraxTrace0(const UpcLens *self, double h0, double a0d, int iwav, double *hv, double *avd)
{
    double a1, h1, ai, hi;
    int i;

    if (iwav < 0 || iwav >= self->nwav)
        return FALSE;
    hv[0] = h0;
    avd[0] = a0d;
    a1 = a0d;
    h1 = h0 - UpcSurf_e(self->s[0], iwav) * a0d;
    ai = a1;
    hi = h1;
    for (i = 1; i <= self->nsur; i++) {
        avd[i] = ai + hi * UpcLens_phi(self, i, iwav);
        hv[i] = hi;
        hi = hi - UpcSurf_e(self->s[i], iwav) * avd[i];
        ai = avd[i];
    }
    hv[self->nsur + 1] = hi;
    avd[self->nsur + 1] = 0.0;
    return TRUE;
}

/* 近軸光線追跡を行う.(第1面から)
 * [入力]
 *     self : UpcLensオブジェクト
 *     h1   : 第1面に入射する光線高さ
 *     a1   : 第1面に入射する換算傾角
 *     iwav ; 波長番号
 * [出力]
 *     hv[0..nsur + 1]  : 第v面の光線高さ
 *     avd[0..nsur + 1] : 第v面を射出する光線の換算傾角
 */
Bool UpcLens_paraxTrace1(const UpcLens *self, double h1, double a1, int iwav, double *hv, double *avd)
{
    double ai, hi;
    int i;

    if (iwav < 0 || iwav >= self->nwav)
        return FALSE;
    hv[0] = 0.0;
    avd[0] = a1;
    hi = h1;
    ai = a1;
    for (i = 1; i <= self->nsur + 1; i++) {
        avd[i] = ai + hi * UpcLens_phi(self, i, iwav);
        hv[i] = hi;
        hi = hi - UpcSurf_e(self->s[i], iwav) * avd[i];
        ai = avd[i];
    }
    hv[self->nsur + 1] = hi;
    avd[self->nsur + 1] = 0.0;
    return TRUE;
}

/* 近軸2光線を追跡する.
 * [入力]
 *     self : UpcLensオブジェクト
 *     iwav ; 波長番号
 * [出力]
 *     pray_hv[0..nsur+1]  : 
 *     pray_avd[0..nsur+1] : 
 *     mray_hv[0..nsur+1]  : 
 *     mray_avd[0..nsur+1] : 
 */
Bool UpcLens_paraxTrace2(UpcLens *self, int iwav, double *pray_hv, double *pray_avd, double *mray_hv, double *mray_avd)
{
    double N1, s1, m[4], delta, g1h, t1, g1;
    double hm1, am1, ap1, hp1;

    if (iwav < 0 || iwav >= self->nwav)
        return FALSE;

    N1 = self->s[0]->_n[iwav];                // 第1面の前の屈折率
    s1 = - self->s[0]->d;                     // 第1面から物体面までの距離
    UpcLens_paraxMatrix(self, 1, self->nsur, iwav, m);
    delta = m[2] ? (1.0 - m[3]) / m[2] : 0.0; // 第1面から前側主点までの距離
    g1h = s1 - delta;                         // 前側主点から物体面までの距離
    UpcLens_getParax(self, iwav);
    t1 = self->parax[iwav].enp;               // 第1面から入射瞳までの距離
    g1 = s1 - t1;                             // 入射瞳から物体面までの距離

    if (self->infObjD) {
        hm1 =  1.0;
        am1 =  0.0;
        ap1 = -1.0;
        hp1 = -t1 / N1;
    }
    else if (fabs(t1) >= 1.0e6) {
        hm1 = s1 / g1h;
        am1 = N1 / g1h;
        ap1 = 0.0;
        hp1 = g1h / N1;
    }
    else {
        hm1 = s1 / g1h;
        am1 = N1 / g1h;
        ap1 = -g1h / g1;
        hp1 = ap1 * (t1 / N1);
    }
    UpcLens_paraxTrace1(self, hp1, ap1, iwav, pray_hv, pray_avd);
    UpcLens_paraxTrace1(self, hm1, am1, iwav, mray_hv, mray_avd);
    return TRUE;
}

/* グローバル座標で表現された各面の頂点位置とx, y, z軸のリストを返す.
 * グローバル座標系は第1面の座標系と一致する.
 * [入力]
 *     UpcLensオブジェクト
 * [出力]
 *     pos[0..nsur + 1][X, Y, Z] : 面頂点位置
 *     xax[0..nsur + 1][X, Y, Z] : X軸方向余弦
 *     yax[0..nsur + 1][X, Y, Z] : Y軸方向余弦
 *     zax[0..nsur + 1][X, Y, Z] : Z軸方向余弦
 */
void UpcLens_globalCoord_old(const UpcLens *self, double **pos, double **xax, double **yax, double **zax)
{
    int isur, isur_img = self->nsur + 1;
    UpcDecenter *sdec;
    
    UpcVector_SET(pos[1], 0.0, 0.0, 0.0);
    UpcVector_SET(xax[1], 1.0, 0.0, 0.0);
    UpcVector_SET(yax[1], 0.0, 1.0, 0.0);
    UpcVector_SET(zax[1], 0.0, 0.0, 1.0);
    for (isur = 2; isur <= isur_img; isur++) {
        sdec = self->s[isur]->dec;
        UpcVector_LINCOMB(pos[isur], 1.0, pos[isur - 1], self->s[isur - 1]->d, zax[isur - 1]);
        if (sdec->decenterFlg) {
            UpcVector_LINCOMB(pos[isur], 1.0, pos[isur], sdec->shift[X], xax[isur]);
            UpcVector_LINCOMB(pos[isur], 1.0, pos[isur], sdec->shift[Y], yax[isur]);
            UpcVector_LINCOMB(pos[isur], 1.0, pos[isur], sdec->shift[Z], zax[isur]);
            UpcDecenter_rotReverse(sdec, xax[isur - 1], xax[isur]);
            UpcDecenter_rotReverse(sdec, yax[isur - 1], yax[isur]);
            UpcDecenter_rotReverse(sdec, zax[isur - 1], zax[isur]);
        }
        else {
            UpcVector_COPY(xax[isur], xax[isur - 1]);
            UpcVector_COPY(yax[isur], yax[isur - 1]);
            UpcVector_COPY(zax[isur], zax[isur - 1]);
        }
    }
    sdec = self->s[1]->dec;
    if (sdec->decenterFlg) {
        UpcDecenter_rotForward(sdec, xax[1], xax[0]);
        UpcDecenter_rotForward(sdec, yax[1], yax[0]);
        UpcDecenter_rotForward(sdec, zax[1], zax[0]);
        UpcVector_LINCOMB(pos[0], 1.0, pos[1], -(sdec->shift[X]), xax[0]);
        UpcVector_LINCOMB(pos[0], 1.0, pos[0], -(sdec->shift[Y]), yax[0]);
        UpcVector_LINCOMB(pos[0], 1.0, pos[0], -(sdec->shift[Z]), zax[0]);
    }
    else {
        UpcVector_COPY(pos[0], pos[1]);
        UpcVector_COPY(xax[0], xax[1]);
        UpcVector_COPY(yax[0], yax[1]);
        UpcVector_COPY(zax[0], zax[1]);
    }
    UpcVector_LINCOMB(pos[0], 1.0, pos[0], -(self->s[0]->d), zax[0]);
}

/* 各面のcoordオブジェクトにグローバル座標系をセットする.
 * グローバル座標系は第1面のローカル座標系と一致する.
 * [入力]
 *     UpcLensオブジェクト
 */
static void UpcLens_setGlobalCoord(const UpcLens *self)
{
    int isur, isur_img = self->nsur + 1;
    UpcSurf **sur = self->s;
    UpcDecenter *sdec;
    UpcCoord *current;
    UpcCoord *prev;
    
    current = sur[1]->coord;
    UpcVector_set(current->eo, 0.0, 0.0, 0.0);
    UpcVector_set(current->ex, 1.0, 0.0, 0.0);
    UpcVector_set(current->ey, 0.0, 1.0, 0.0);
    UpcVector_set(current->ez, 0.0, 0.0, 1.0);
    UpcCoord_setGlobal(current);

    for (isur = 2; isur <= isur_img; isur++) {
        current = sur[isur]->coord;
        prev = sur[isur - 1]->coord;
        sdec = sur[isur]->dec;

        UpcVector_linComb(current->eo, 1.0, prev->eo, sur[isur - 1]->d, prev->ez);
        if (sdec->decenterFlg) {
            if (sdec->shift[X]) {
                //UpcVector_linComb(current->eo, 1.0, current->eo, sdec->shift[X], current->ex);
                UpcVector_linComb(current->eo, 1.0, current->eo, sdec->shift[X], prev->ex);
            }
            if (sdec->shift[Y]) {
                //UpcVector_linComb(current->eo, 1.0, current->eo, sdec->shift[Y], current->ey);
                UpcVector_linComb(current->eo, 1.0, current->eo, sdec->shift[Y], prev->ey);
            }
            if (sdec->shift[Z]) {
                //UpcVector_linComb(current->eo, 1.0, current->eo, sdec->shift[Z], current->ez);
                UpcVector_linComb(current->eo, 1.0, current->eo, sdec->shift[Z], prev->ez);
            }
            UpcDecenter_rotReverse(sdec, prev->ex, current->ex);
            UpcDecenter_rotReverse(sdec, prev->ey, current->ey);
            UpcDecenter_rotReverse(sdec, prev->ez, current->ez);
        }
        else {
            UpcVector_COPY(current->ex, prev->ex);
            UpcVector_COPY(current->ey, prev->ey);
            UpcVector_COPY(current->ez, prev->ez);
        }
        UpcCoord_setGlobal(current);
    }
    sdec = sur[1]->dec;
    prev = sur[1]->coord;
    current = sur[0]->coord;
    if (sdec->decenterFlg) {
        UpcDecenter_rotForward(sdec, prev->ex, current->ex);
        UpcDecenter_rotForward(sdec, prev->ey, current->ey);
        UpcDecenter_rotForward(sdec, prev->ez, current->ez);
        UpcVector_linComb(current->eo, 1.0, prev->eo   , -(sdec->shift[X]), current->ex);
        UpcVector_linComb(current->eo, 1.0, current->eo, -(sdec->shift[Y]), current->ey);
        UpcVector_linComb(current->eo, 1.0, current->eo, -(sdec->shift[Z]), current->ez);
    }
    else {
        UpcVector_COPY(current->eo, prev->eo);
        UpcVector_COPY(current->ex, prev->ex);
        UpcVector_COPY(current->ey, prev->ey);
        UpcVector_COPY(current->ez, prev->ez);
    }
    UpcVector_linComb(current->eo, 1.0, current->eo, -(sur[0]->d), current->ez);
    UpcCoord_setGlobal(sur[0]->coord);
}

/* 3次収差係数を返す.
 * [入力]
 *     UpcLensオブジェクト
 * [出力]
 *     I[0..nsur + 1]   : 球面収差
 *     II[0..nsur + 1]  : コマ
 *     III[0..nsur + 1] : 非点収差
 *     IV[0..nsur + 1]  : 球欠像面湾曲
 *     V[0..nsur + 1]   : 歪曲
 *     P[0..nsur + 1]   : ペッツバール項
 *     Is[0..nsur + 1]  : 瞳の球面収差
 *     IIs[0..nsur + 1] : 瞳のコマ
 *     IIIs[0..nsur + 1]: 瞳の非点収差
 *     IVs[0..nsur + 1] : 瞳の球欠像面湾曲
 *     Vs[0..nsur + 1]  : 瞳の歪曲
 *     各々, 第(nsur + 1)要素には全系の和が格納される.
 */
Bool UpcLensSeidelCoef(UpcLens *self, double *I, double *II, double *III, double *IV, double *V, double *P, double *Is, double *IIs, double *IIIs, double *IVs, double *Vs)
{
    const int nsur = self->nsur;
    const int iwav = self->iwav_pri;
    UpcSurf **sur = self->s;
    double *axpHv, *axpAvd, *axmHv, *axmAvd;
    int isur;
    double Ni, No, hm, ami, amo, hp, api, apo, curv, bv, cv;
    double hQm, hQp, dltm, dltp, psi;
    double hm2, hp2, hm_psi, hp_psi;
    double sumI    = 0.0;
    double sumII   = 0.0;
    double sumIII  = 0.0;
    double sumIV   = 0.0;
    double sumV    = 0.0;
    double sumP    = 0.0;
    double sumIs   = 0.0;
    double sumIIs  = 0.0;
    double sumIIIs = 0.0;
    double sumIVs  = 0.0;
    double sumVs   = 0.0;
    
    axpHv  = (double *)malloc(sizeof(double) * (nsur + 2) * 4);
    if (!axpHv) {
        free(axpHv);
        return FALSE;
    }
    axpAvd = axpHv + (nsur + 2) * 1;
    axmHv  = axpHv + (nsur + 2) * 2;
    axmAvd = axpHv + (nsur + 2) * 3;
    
    UpcLens_paraxTrace2(self, iwav, axpHv, axpAvd, axmHv, axmAvd);
    
    for (isur = 1; isur <= nsur; isur++) {
        Ni  = sur[isur - 1]->_n[iwav];
        No  = sur[isur]->_n[iwav];
        hp  = axpHv[isur];
        apo = axpAvd[isur];
        api = axpAvd[isur - 1];
        hm  = axmHv[isur];
        amo = axmAvd[isur];
        ami = axmAvd[isur - 1];
        sur[isur]->shape->getRBC(sur[isur]->shape, &curv, &bv, &cv);
        // 補助量の計算
        hQm    = hm * Ni * curv - ami;
        hQp    = hp * Ni * curv - api;
        dltm   = amo / (No * No) - ami / (Ni * Ni);
        dltp   = apo / (No * No) - api / (Ni * Ni);
        psi    = bv * (No - Ni);
        hm2    = hm * hm;
        hp2    = hp * hp;
        hm_psi = hm * psi;
        hp_psi = hp * psi;
        // 物体の収差係数
        sumI   += (I[isur]   = (hQm * hQm * dltm + hm2 * hm_psi) * hm);
        sumII  += (II[isur]  = (hQm * hQp * dltm + hm2 * hp_psi) * hm);
        sumIII += (III[isur] = (hQp * hQp * dltm + hp2 * hm_psi) * hm);
        sumP   += (P[isur]   = - curv * (1.0 / No - 1.0 / Ni));
        sumIV  += (IV[isur]  = III[isur] + P[isur]);
        sumV   += (V[isur]   = (hQp * hQp * dltm + hp2 * hm_psi) * hp + hQp * dltp);
        // 瞳の収差係数
        sumIs   += (Is[isur]   = (hQp * hQp * dltp + hp2 * hp_psi) * hp);
        sumIIs  += (IIs[isur]  = (hQm * hQp * dltp + hp2 * hm_psi) * hp);
        sumIIIs += (IIIs[isur] = (hQm * hQm * dltp + hm2 * hp_psi) * hp);
        sumIVs  += (IVs[isur]  = IIIs[isur] + P[isur]);
        sumVs   += (Vs[isur]   = (hQm * hQm * dltp + hm2 * hp_psi) * hm - hQm * dltm);
    }

    I  [nsur + 1] = sumI;
    II [nsur + 1] = sumII;
    III[nsur + 1] = sumIII;
    IV [nsur + 1] = sumIV;
    V  [nsur + 1] = sumV;
    P  [nsur + 1] = sumP;
    Is  [nsur + 1] = sumIs;
    IIs [nsur + 1] = sumIIs;
    IIIs[nsur + 1] = sumIIIs;
    IVs [nsur + 1] = sumIVs;
    Vs  [nsur + 1] = sumVs;

    free(axpHv);
    return TRUE;
}

