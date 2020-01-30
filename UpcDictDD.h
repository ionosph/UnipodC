/*
 *  UpcDictDD.h
 *  UnipodC
 *
 *  Created by ionosph on 12/01/27.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef _UPCDICTDD_H
#define _UPCDICTDD_H

#include "UpcBase.h"


typedef double keyType;
typedef double valType;

typedef struct _UpcRBNodeDD {
    keyType key;
    valType val;
    int color;
    struct _UpcRBNodeDD *left;
    struct _UpcRBNodeDD *right;
} UpcRBNodeDD;


typedef struct _UpcDictDD UpcDictDD;
struct _UpcDictDD {
    int _refCount;
    void (*_dealloc)(UpcDictDD *);
    UpcRBNodeDD *root;
};

UpcDictDD *UpcDictDD_init(void);
Bool UpcDictDD_search(UpcDictDD *self, keyType key, valType *val);
void UpcDictDD_add(UpcDictDD *self, keyType key, valType val);
void UpcDictDD_delete(UpcDictDD *self, keyType key);
void UpcDictDD_print(UpcDictDD *self, FILE *ostream);


Bool TestUpcDictDD(void);

#endif
