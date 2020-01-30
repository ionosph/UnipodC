/*
 *  UpcRayTrace.h
 *  UnipodC
 *
 *  Created by ionosph on 2011/04/05.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _UPCRAYTRACE_H
#define _UPCRAYTRACE_H

#include "UpcBase.h"
#include "UpcVector.h"
#include "UpcDecenter.h"
#include "UpcSurf.h"
#include "UpcLens.h"
#include "UpcRay.h"


typedef struct _UpcTraceFactory UpcTraceFactory;
struct _UpcTraceFactory {
    int _refCount;
    void (*_dealloc)(UpcTraceFactory *);
    // attribute
    int nwav;
    int nhgt;
    int nray;
    // double tol_pray;
    // double tol_mray;
    UpcRay ****refRay; // refray[nwav][nhgt][nray]
    double ***pInfo;
    double ***wInfo;
    Bool **refRayStatus;
    Bool **readyWabc;
    // method
    Bool (*traceP)(const UpcTraceFactory *self, const UpcLens *lobj, int iwav, int ihgt, double px, double py, UpcRay *ray);
};


Bool UpcRayTrace_trace0(const UpcLens *lobj, const double *s0pos, const double *tpos, int iwav, int isur_end, UpcRay *ray);
Bool UpcRayTrace_trace0d(const UpcLens *lobj, const double *s0pos, const double *s0dir, int iwav, int isur_end, UpcRay *ray);
Bool UpcRayTrace_trace0p(const UpcLens *lobj, const double *s0pos, const double *s0dir, double sinTx, double sinTy, int iwav, UpcRay *ray);
Bool UpcRayTrace_trace1(const UpcLens *lobj, const double *s0dir, const double *tpos, int iwav, int isur_end, UpcRay *ray);
Bool UpcRayTrace_trace_rev(const UpcLens *lobj, const double *sipos, const double *skdir, int iwav, UpcRay *ray);
Bool UpcRayTrace_psearch0(const UpcLens *lobj, double ox, double oy, int isur, double x, double y, int iwav, double tol, UpcRay *ray);
Bool UpcRayTrace_psearch1(const UpcLens *lobj, double d1x, double d1y, int isur, double x, double y, int iwav, double tol, UpcRay *ray);
Bool UpcRayTrace_psearch0_xt(const UpcLens *lobj, double ox, double oy, double cosx, double cosy, int iwav, double tol, UpcRay *ray);
Bool UpcRayTrace_msearch0_ea(const UpcLens *lobj, const UpcRay *pray, double dx, double dy, int iwav, double tol, UpcRay *ray);
Bool UpcRayTrace_msearch1_ea(const UpcLens *lobj, const UpcRay *pray, double dx, double dy, int iwav, double tol, UpcRay *ray);
Bool UpcRayTrace_msearch0_na(const UpcLens *lobj, const UpcRay *pray, double sin_theta, double dx, double dy, double expz, int iwav, double tol, UpcRay *ray);
Bool UpcRayTrace_msearch1_na(const UpcLens *lobj, const UpcRay *pray, double sin_theta, double dx, double dy, double expz, int iwav, double tol, UpcRay *ray);
Bool UpcRayTrace_msearch0_nao(const UpcLens *lobj, const UpcRay *pray, double sin_theta, double dx, double dy, int iwav, UpcRay *ray);
Bool UpcRayTrace_calWab(const UpcLens *lobj, const UpcRay *ray, const double *si, const double *ni, double ci, const double *so, const double *no, double co, int iwav, double *wab);
Bool UpcRayTrace_calAstigmatism(const UpcLens *lobj, const UpcRay *ray, double h1, double u0, int iwav, double *h_m, double *u_m, double *h_s, double *u_s);

UpcTraceFactory *UpcTraceFactory_init(const UpcLens *lobj);
Bool UpcTraceFactory_asTrace(const UpcLens *lobj, double y, int iwav, double *f_meri, double *f_sagi, double *dist);
Bool UpcTraceFactory_wabc(const UpcLens *lobj, int iwav, int ihgt, double px, double py, double *wab);


#endif

