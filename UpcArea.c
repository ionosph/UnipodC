/*
 *  UpcArea.c
 *  UnipodC
 *
 *  Created by ionosph on 2010/10/22.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcArea.h"

static const double centerDefault  =  0.0;
static const double semiDiaDefault =  UPC_INFSD;
static const double upperDefault   =  UPC_INFSD;
static const double lowerDefault   = -UPC_INFSD;

static void UpcArea_dealloc(UpcArea *self)
{
    free(self);
}

static void UpcArea_print(const UpcArea *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);

    fprintf(ostream, "%scenter = (%13.6e, %13.6e)\n", space, self->center[X] , self->center[Y] );
    fprintf(ostream, "%ssd     = (%13.6e, %13.6e)\n", space, self->semiDia[X], self->semiDia[Y]);
    fprintf(ostream, "%supper  = (%13.6e, %13.6e)\n", space, self->upper[X]  , self->upper[Y]  );
    fprintf(ostream, "%slower  = (%13.6e, %13.6e)\n", space, self->lower[X]  , self->lower[Y]  );
}

static void UpcArea_toXml(const UpcArea *self, int indent, FILE *ostream)
{
    const char *space = upc_space_string(indent);

    fprintf(ostream, "%s<area>\n", space);
    if (self->center[X] != centerDefault)
        fprintf(ostream, "%s  <center_x>%23.15e</center_x>\n", space, self->center[X]);
    if (self->center[Y] != centerDefault)
        fprintf(ostream, "%s  <center_y>%23.15e</center_y>\n", space, self->center[Y]);
    if (self->semiDia[X] != semiDiaDefault)
        fprintf(ostream, "%s  <sd_x>%23.15e</sd_x>\n", space, self->semiDia[X]);
    if (self->semiDia[Y] != semiDiaDefault)
        fprintf(ostream, "%s  <sd_y>%23.15e</sd_y>\n", space, self->semiDia[Y]);
    if (self->upper[X] != upperDefault)
        fprintf(ostream, "%s  <upper_x>%23.15e</upper_x>\n", space, self->upper[X]);
    if (self->upper[Y] != upperDefault)
        fprintf(ostream, "%s  <upper_y>%23.15e</upper_y>\n", space, self->upper[Y]);
    if (self->lower[X] != lowerDefault)
        fprintf(ostream, "%s  <lower_x>%23.15e</lower_x>\n", space, self->lower[X]);
    if (self->lower[Y] != lowerDefault)
        fprintf(ostream, "%s  <lower_y>%23.15e</lower_y>\n", space, self->lower[Y]);
    fprintf(ostream, "%s</area>\n", space);
}

static UpcErrStatus UpcArea_lset(UpcArea *self)
{
    return UpcE_NoIssues;
}

static UpcArea *UpcArea_copy(UpcArea *self)
{
    UpcArea *t = (UpcArea *)malloc(sizeof(UpcArea));
    
    if (t) {
        memcpy(t, self, sizeof(UpcArea));
        t->_refCount = 1;
        /*
        t->_refCount = 1;
        t->_dealloc = UpcArea_dealloc;
        t->center[X] = self->center[X];
        t->center[Y] = self->center[Y];
        t->semiDia[X] = self->semiDia[X];
        t->semiDia[Y] = self->semiDia[Y];
        t->upper[X] = self->upper[X];
        t->upper[Y] = self->upper[Y];
        t->lower[X] = self->lower[X];
        t->lower[Y] = self->lower[Y];
        t->copy = UpcArea_copy;
        t->print = UpcArea_print;
        t->toXml = UpcArea_toXml;
        t->lset = UpcArea_lset;
         */
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcArea_copy");
    }
    return t;
}

UpcArea *UpcArea_init(void)
{
    UpcArea *self = (UpcArea *)malloc(sizeof(UpcArea));

    if (self) {
        self->_refCount = 1;
        self->_dealloc = UpcArea_dealloc;
        self->center[X] = centerDefault;
        self->center[Y] = centerDefault;
        self->semiDia[X] = semiDiaDefault;
        self->semiDia[Y] = semiDiaDefault;
        self->upper[X] = upperDefault;
        self->upper[Y] = upperDefault;
        self->lower[X] = lowerDefault;
        self->lower[Y] = lowerDefault;
        self->copy = UpcArea_copy;
        self->print = UpcArea_print;
        self->toXml = UpcArea_toXml;
        self->lset = UpcArea_lset;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcArea_init");
    }
    return self;
}

Bool UpcArea_isActive(const UpcArea *self, double x, double y)
{
    double xd = (x - self->center[X]) / self->semiDia[X];
    double yd = (y - self->center[Y]) / self->semiDia[Y];
    
    if (xd * xd + yd * yd <= 1.0 &&
        x >= self->lower[X] &&
        x <= self->upper[X] &&
        y >= self->lower[Y] &&
        y <= self->upper[Y]) {
        return TRUE;
    }
    return FALSE;
}

double UpcArea_semiDia(const UpcArea *self)
{
    if (self->semiDia[X] == self->semiDia[Y])
        return self->semiDia[X];
    else
        return 0.5 * (self->semiDia[X] + self->semiDia[Y]);
}

void UpcArea_setSemiDia(UpcArea *self, double sd)
{
    self->semiDia[X] = self->semiDia[Y] = sd;
}

double UpcArea_dia(const UpcArea *self)
{
    if (self->semiDia[X] == self->semiDia[Y])
        return 2.0 * self->semiDia[X];
    else
        return 2.0 * 0.5 * (self->semiDia[X] + self->semiDia[Y]);
}

void UpcArea_setDia(UpcArea *self, double dia)
{
    self->semiDia[X] = self->semiDia[Y] = 0.5 * dia;
}

Bool UpcArea_xSymmetric(const UpcArea *self)
{
    if (self->center[Y] || (self->upper[Y] + self->lower[Y]))
        return FALSE;
    return TRUE;
}

Bool UpcArea_ySymmetric(const UpcArea *self)
{
    if (self->center[X] || (self->upper[X] + self->lower[X]))
        return FALSE;
    return TRUE;
}

Bool UpcArea_rSymmetric(const UpcArea *self)
{
    if (self->center[X] || self->center[Y])
        return FALSE;
    if (self->semiDia[X] - self->semiDia[Y])
        return FALSE;
    if (self->upper[X] < 1.0e30 || self->upper[Y] < 1.0e30)
        return FALSE;
    if (self->lower[X] > -1.0e30 || self->lower[Y] > -1.0e30)
        return FALSE;
    return TRUE;
}

