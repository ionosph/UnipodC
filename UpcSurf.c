/*
 *  UpcSurf.c
 *  UnipodC
 *
 *  Created by ionosph on 2011/04/01.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcSurf.h"

static const UpcSurfRmode surfModeDefault = RmodeRefractive;

char *UpcSurfRmode_toString(UpcSurfRmode i)
{
    switch (i) {
        case RmodeRefractive:
            return "refractive";
        case RmodeReflective:
            return "reflective";
    }
    return "unknown mode";
}

UpcSurfRmode UpcSurfRmode_fromString(const char *s)
{
    if (!strcmp(s, "reflective"))
        return RmodeReflective;
    return RmodeRefractive;
}

static void UpcSurf_dealloc(UpcSurf *self)
{
    UpcObj_RELEASE(self->shape);
    UpcObj_RELEASE(self->area);
    UpcObj_RELEASE(self->dec);
    UpcObj_RELEASE(self->coord);
    free(self);
}

UpcSurf *UpcSurf_init(void)
{
    UpcSurf *self = (UpcSurf *)malloc(sizeof(UpcSurf));

    if (self) {
        self->_refCount = 1;
        self->_dealloc = UpcSurf_dealloc;
        self->d = 0.0;
        self->mode = surfModeDefault;
        self->shape = UpcShapeSph_init();
        self->medium = NULL;
        self->_n = NULL;
        self->area = UpcArea_init();
        self->dec = UpcDecenter_init();
        self->coord = UpcCoord_init();
        if (!(self->shape) || !(self->area) || !(self->dec) || !(self->coord)) {
            if (self->shape) {
                self->shape->_dealloc(self->shape);
            }
            if (self->area) {
                self->area->_dealloc(self->area);
            }
            if (self->dec) {
                self->dec->_dealloc(self->dec);
            }
            if (self->coord) {
                self->coord->_dealloc(self->coord);
            }
            free(self);
            UpcERRHANDLER(UpcE_MemoryError, "in UpcSurf_init");
            return NULL;
        }
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcSurf_init");
    }
    return self;
}

UpcSurf *UpcSurf_copy(UpcSurf *self)
{
    UpcSurf *t = (UpcSurf *)malloc(sizeof(UpcSurf));
    
    if (t) {
        t->_refCount = 1;
        t->_dealloc = UpcSurf_dealloc;
        t->d = self->d;
        t->mode = self->mode;
        t->medium = self->medium; // 弱参照
        t->_n = self->_n; // 弱参照
        t->shape = self->shape->copy(self->shape);
        t->area = self->area->copy(self->area);
        t->dec = self->dec->copy(self->dec);
        if (!(t->shape) || !(t->area) || !(t->dec)) {
            if (t->shape) {
                t->shape->_dealloc(t->shape);
            }
            if (t->area) {
                t->area->_dealloc(t->area);
            }
            if (t->dec) {
                t->dec->_dealloc(t->dec);
            }
            if (self->coord) {
                t->coord->_dealloc(t->coord);
            }
            free(t);
            UpcERRHANDLER(UpcE_MemoryError, "in UpcSurf_copy");
            return NULL;
        }
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcSurf_copy");
    }
    return t;
}

UpcErrStatus UpcSurf_lset(UpcSurf *self)
{
    UpcErrStatus stts;

    if (!(self->shape)) {
        UpcERRHANDLER(UpcE_ValueError, "shape object is NULL");
        return UpcE_ValueError;
    }
    if (!(self->medium)) {
        UpcERRHANDLER(UpcE_ValueError, "medium object is NULL");
        return UpcE_ValueError;
    }
    if (!(self->area)) {
        UpcERRHANDLER(UpcE_ValueError, "area object is NULL");
        return UpcE_ValueError;
    }
    if (!(self->dec)) {
        UpcERRHANDLER(UpcE_ValueError, "dec object is NULL");
        return UpcE_ValueError;
    }
    if (!(self->coord)) {
        UpcERRHANDLER(UpcE_ValueError, "coord object is NULL");
        return UpcE_ValueError;
    }

    if ((stts = (self->shape->lset(self->shape))) != UpcE_NoIssues) {
        UpcERRHANDLER(stts, "in UpcSurf_lset");
        return stts;
    }
    if ((stts = (self->medium->lset(self->medium))) != UpcE_NoIssues) {
        UpcERRHANDLER(stts, "in UpcSurf_lset");
        return stts;
    }
    if ((stts = (self->area->lset(self->area))) != UpcE_NoIssues) {
        UpcERRHANDLER(stts, "in UpcSurf_lset");
        return stts;
    }
    if ((stts = (self->dec->lset(self->dec))) != UpcE_NoIssues) {
        UpcERRHANDLER(stts, "in UpcSurf_lset");
        return stts;
    }
    return UpcE_NoIssues;
}

void UpcSurf_print(const UpcSurf *self, int isur, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);

    fprintf(ostream, "%sd      = %13.6e\n", space, self->d);
    if (self->medium)
        fprintf(ostream, "%smedium = %s\n", space, self->medium->name);
    else
        fprintf(ostream, "%smedium = NULL\n", space);
    fprintf(ostream, "%sshape\n", space);
    if (self->shape)
        self->shape->print(self->shape, indent + 2, ostream);
    else
        fprintf(ostream, "%s  NULL\n", space);
    fprintf(ostream, "%sarea\n", space);
    if (self->area)
        self->area->print(self->area, indent + 2, ostream);
    else
        fprintf(ostream, "%s  NULL\n", space);
    fprintf(ostream, "%sdec\n", space);
    if (self->dec)
        self->dec->print(self->dec, indent + 2, ostream);
    else
        fprintf(ostream, "%s  NULL\n", space);
}

void UpcSurf_toXml(const UpcSurf *self, int isur, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);

    fprintf(ostream, "%s<surf isur=\"%d\">\n", space, isur);

    fprintf(ostream, "%s  <d>%23.15e</d>\n", space, self->d);
    if (self->mode != surfModeDefault)
        fprintf(ostream, "%s  <mode>%s</mode>\n", space, UpcSurfRmode_toString(self->mode));
    if (self->shape) {
        self->shape->toXml(self->shape, indent + 2, ostream);
    }
    else {
        fprintf(ostream, "%s  <!-- warning: shape object is NULL -->\n", space);
    }
    if (self->medium) {
        fprintf(ostream, "%s  <medium>%s</medium>\n", space, self->medium->name);
    }
    else {
        fprintf(ostream, "%s  <!-- warning: medium object is NULL -->\n", space);
    }
    if (self->area) {
        self->area->toXml(self->area, indent + 2, ostream);
    }
    else {
        fprintf(ostream, "%s  <!-- warning: area object is NULL -->\n", space);
    }
    if (self->dec) {
        self->dec->toXml(self->dec, indent + 2, ostream);
    }
    else {
        fprintf(ostream, "%s  <!-- warning: decenter object is NULL -->\n", space);
    }
    fprintf(ostream, "%s</surf>\n", space);
}

double UpcSurf_e(UpcSurf *self, int iwav)
{
    if (self->_n[iwav]) {
        return self->d / self->_n[iwav];
    }
    return HUGE_VAL;
}

Bool UpcSurf_xSymmetric(const UpcSurf *self)
{
    if (self->shape->xSymmetric(self->shape) && UpcArea_xSymmetric(self->area) && UpcDecenter_xSymmetric(self->dec))
        return TRUE;
    return FALSE;
}

Bool UpcSurf_ySymmetric(const UpcSurf *self)
{
    if (self->shape->ySymmetric(self->shape) && UpcArea_ySymmetric(self->area) && UpcDecenter_ySymmetric(self->dec))
        return TRUE;
    return FALSE;
}

Bool UpcSurf_rSymmetric(const UpcSurf *self)
{
    if (self->shape->rSymmetric(self->shape) && UpcArea_rSymmetric(self->area) && UpcDecenter_rSymmetric(self->dec))
        return TRUE;
    return FALSE;
}

UpcErrStatus UpcSurf_setShapeTypeSph(UpcSurf *self)
{
    UpcShape *t = UpcShapeSph_initWithShape(self->shape);
    
    if (t) {
        UpcObj_RELEASE(self->shape);
        self->shape = t;
        return UpcE_NoIssues;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcSurf_setShapeTypeSph");
        return UpcE_MemoryError;
    }
}

UpcErrStatus UpcSurf_setShapeTypeAsph(UpcSurf *self, int max_orderR)
{
    UpcShape *t = UpcShapeAsph_initWithShape(max_orderR, self->shape);
    
    if (t) {
        UpcObj_RELEASE(self->shape);
        self->shape = t;
        return UpcE_NoIssues;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcSurf_setShapeTypeAsph");
        return UpcE_MemoryError;
    }
}

UpcErrStatus UpcSurf_setShapeTypeXYP(UpcSurf *self, int max_orderR, int max_orderX, int max_orderY)
{
    UpcShape *t = UpcShapeXYP_initWithShape(max_orderR, max_orderX, max_orderY, self->shape);
    
    if (t) {
        UpcObj_RELEASE(self->shape);
        self->shape = t;
        return UpcE_NoIssues;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcSurf_setShapeTypeXYP");
        return UpcE_MemoryError;
    }
}

UpcErrStatus UpcSurf_setShapeTypeZernike(UpcSurf *self, int max_orderR, int zcn_max)
{
    UpcShape *t = UpcShapeZernike_initWithShape(max_orderR, zcn_max, self->shape);
    
    if (t) {
        UpcObj_RELEASE(self->shape);
        self->shape = t;
        return UpcE_NoIssues;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcSurf_setShapeTypeZernike");
        return UpcE_MemoryError;
    }
}





