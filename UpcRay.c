/*
 *  UpcRay.c
 *  UnipodC
 *
 *  Created by ionosph on 2011/04/05.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcRay.h"

/* (注) NULLを渡さないこと.
 */
static void UpcRay_dealloc(UpcRay *self)
{
    free(self->s);
    free(self);
}

/* (注) nsur >= 1 であること.(呼び出し側がチェックすること)
 */
UpcRay *UpcRay_init(int nsur)
{
    UpcRay *self = (UpcRay *)malloc(sizeof(UpcRay));

    if (self) {
        self->_refCount = 1;
        self->_dealloc = UpcRay_dealloc;
        self->s = (UpcRayComponent *)malloc(sizeof(UpcRayComponent) * (nsur + 2));
        if (!self->s) {
            free(self);
            UpcERRHANDLER(UpcE_MemoryError, "in UpcRay_init");
            return NULL;
        }
        self->nsur = nsur;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcRay_init");
    }
    return self;
}

void UpcRay_print(const UpcRay *self, char *comment, FILE *ostream)
{
    int i, n = self->nsur + 2;
    char strPosX[16], strPosY[16], strPosZ[16], strOpl[16];
    
    fprintf(ostream, "[%s]\n", comment);
    fprintf(ostream, "  status, %s\n", UpcRayStatus_string(self->status));
    fprintf(ostream, "  sur,    pos[X] ,    pos[Y] ,    pos[Z] ,  dir[X]  ,  dir[Y]  ,  dir[Z]  ,    opl    , status\n");
    for (i = 0; i < n; i++) {
        UpcRayComponent s = self->s[i];

        sprintf(strPosX, (fabs(s.pos[X]) < 1000.0) ? "%11.6f" : "%11.3e", s.pos[X]);
        sprintf(strPosY, (fabs(s.pos[Y]) < 1000.0) ? "%11.6f" : "%11.3e", s.pos[Y]);
        sprintf(strPosZ, (fabs(s.pos[Z]) < 1000.0) ? "%11.6f" : "%11.3e", s.pos[Z]);
        sprintf(strOpl , (fabs(s.opl   ) < 1000.0) ? "%11.6f" : "%11.3e", s.opl);

        fprintf(ostream, "  %3d,%s,%s,%s,%10.7f,%10.7f,%10.7f,%s, %s\n"
                , i, strPosX, strPosY, strPosZ, s.dir[X], s.dir[Y], s.dir[Z], strOpl, UpcRayStatus_string(s.status));
    }
}

/* other を self にコピーする.
 */
void UpcRay_copy(UpcRay *self, const UpcRay *other)
{
    int i, n = self->nsur + 2;

    self->status = other->status;
    for (i = 0; i < n; i++)
        self->s[i] = other->s[i]; // 構造体の実体の代入
}

/* other をYZ面に関して反転して self にコピーする.
 */
void UpcRay_copy_YZmirror(UpcRay *self, const UpcRay *other)
{
    int i, n = self->nsur + 2;
    UpcRayComponent *s, *o;
    
    self->status = other->status;
    for (i = 0; i < n; i++) {
        s = self->s + i;
        o = other->s + i;
        s->pos[X] = - o->pos[X];
        s->pos[Y] =   o->pos[Y];
        s->pos[Z] =   o->pos[Z];
        s->dir[X] = - o->dir[X];
        s->dir[Y] =   o->dir[Y];
        s->dir[Z] =   o->dir[Z];
        s->opl = o->opl;
        s->cosI = o->cosI;
        s->cosO = o->cosO;
        s->status = o->status;
    }
}

/* other をXZ面に関して反転して self にコピーする.
 */
void UpcRay_copy_XZmirror(UpcRay *self, const UpcRay *other)
{
    int i, n = self->nsur + 2;
    UpcRayComponent *s, *o;

    self->status = other->status;
    for (i = 0; i < n; i++) {
        s = self->s + i;
        o = other->s + i;
        s->pos[X] =   o->pos[X];
        s->pos[Y] = - o->pos[Y];
        s->pos[Z] =   o->pos[Z];
        s->dir[X] =   o->dir[X];
        s->dir[Y] = - o->dir[Y];
        s->dir[Z] =   o->dir[Z];
        s->opl = o->opl;
        s->cosI = o->cosI;
        s->cosO = o->cosO;
        s->status = o->status;
    }
}

void UpcRay_resetStatus(UpcRay *self)
{
    int i, n = self->nsur + 2;

    self->status = UpcRayE_UnCalculated;
    for (i = 0; i < n; i++)
        self->s[i].status = UpcRayE_UnCalculated;
}

const char *UpcRayStatus_string(UpcRayStatus status)
{
    switch (status) {
        case UpcRayE_Succeeded:
            return "Succeeded";
        case UpcRayE_OutOfArea:
            return "OutOfArea";
        case UpcRayE_UnCalculated:
            return "UnCalculated";
        case UpcRayE_NoIntersectionError:
            return "NoIntersectionError";
        case UpcRayE_TotalReflectionError:
            return "TotalReflectionError";
        case UpcRayE_SettingError:
            return "SettingError";
        case UpcRayE_IterationError:
            return "IterationError";
    }
    return "unknown";
}




