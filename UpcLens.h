/*
 *  UpcLens.h
 *  UnipodC
 *
 *  Created by ionosph on 2011/04/02.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _UPCLENS_H
#define _UPCLENS_H

#include "UpcBase.h"
#include "UpcSurf.h"
#include "UpcGlass.h"
#include "UpcDict.h"
//#include "UpcRayTrace.h"

enum UpcPrayMode {
    UpcPrayMode_stop,
    UpcPrayMode_itel,
    UpcPrayMode_otel
};

enum UpcMrayMode {
    UpcMrayMode_eape,
    UpcMrayMode_isin,
    UpcMrayMode_osin
};

typedef struct UpcParaxData {
    double fl;  // 焦点距離
    double mag; // 結像倍率
    double sk;  // 最終面から近軸像面までの距離
    double enp; // 第1面から入射瞳までの距離
    double exp; // 最終面から射出瞳までの距離
    double o1;  // 第1面から前側主点までの距離
    double ok;  // 最終面から後側主点までの距離
} UpcParaxData;

struct UpcTraceFactory;

typedef struct _UpcLens UpcLens;
struct _UpcLens {
    int _refCount;
    void (*_dealloc)(UpcLens *);
    // attribute
    int nsur;
    int nwav;
    int nhgt;
    char *title;
    char *comment;
    int isur_sto;
    int iwav_pri;
    Bool pim;
    double *wa;
    double *ha;
    enum UpcPrayMode pray;
    enum UpcMrayMode mray;
    double aperture;
    double tol_pray;
    double tol_mray;
    UpcDict *gcatalog;
    UpcSurf **s;
    Bool xSymmetric;
    Bool ySymmetric;
    Bool rSymmetric;
    Bool infObjD;
    UpcParaxData *parax;
    struct _UpcTraceFactory *raytracef;
};

UpcLens *UpcLens_init(int nsur, int nwav, int nhgt);
UpcLens *UpcLens_copy(UpcLens *self);
void UpcLens_print(const UpcLens *self, FILE *ostream);
void UpcLens_toXml(const UpcLens *self, int indent, FILE *ostream);
UpcErrStatus UpcLens_lset(UpcLens *self);
void UpcLens_setTitle(UpcLens *self, const char *title);
void UpcLens_setComment(UpcLens *self, const char *comment);

Bool UpcLens_setNwav(UpcLens *self, int nwav_new);
Bool UpcLens_setNhgt(UpcLens *self, int nhgt_new);
Bool UpcLens_insertSurf(UpcLens *self, int isur_insert);
Bool UpcLens_deleteSurf(UpcLens *self, int isur_delete);

void UpcLens_paraxMatrix(const UpcLens *self, int isur_stt, int isur_end, int iwav, double *m);
Bool UpcLens_getParax(UpcLens *self, int iwav);
char *UpcPrayMode_toString(enum UpcPrayMode i);
char *UpcMrayMode_toString(enum UpcMrayMode i);
enum UpcPrayMode UpcPrayMode_fromString(const char *s);
enum UpcMrayMode UpcMrayMode_fromString(const char *s);
UpcErrStatus UpcLens_saveXml(UpcLens *self, char *filename);

#define UpcLens_WL(self,iwav) ((self)->wa[2*(iwav)])
#define UpcLens_WTW(self,iwav) ((self)->wa[2*(iwav)+1])
#define UpcLens_HGTX(self,ihgt) ((self)->ha[7*(ihgt)])
#define UpcLens_HGTY(self,ihgt) ((self)->ha[7*(ihgt)+1])
#define UpcLens_WTH(self,ihgt) ((self)->ha[7*(ihgt)+2])
#define UpcLens_VIG(self,ihgt,ixyul) ((self)->ha[7*(ihgt)+((ixyul)+2)])

#endif
