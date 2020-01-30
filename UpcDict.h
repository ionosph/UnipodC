/*
 *  UpcDict.h
 *  UnipodC
 *
 *  Created by ionosph on 12/01/27.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef _UPCDICT_H
#define _UPCDICT_H

#include "UpcBase.h"


typedef struct _UpcRBNode {
    char *key;
    void *val;
    int color;
    struct _UpcRBNode *left;
    struct _UpcRBNode *right;
} UpcRBNode;


typedef struct _UpcDict UpcDict;
struct _UpcDict {
    int _refCount;
    void (*_dealloc)(UpcDict *);
    UpcRBNode *root;
    int (*compKey)(const void *, const void *);
    void (*delKey)(void *);
    void (*delVal)(void *);
};

UpcDict *UpcDict_init(int (*compKey)(const void *, const void *), void (*delKey)(void *), void (*delVal)(void *));
Bool UpcDict_search(UpcDict *self, const void *key, void **val);
void UpcDict_add(UpcDict *self, void *key, void *val);
void UpcDict_delete(UpcDict *self, void *key);
void UpcDict_print(UpcDict *self, FILE *ostream);
void **UpcDict_keys(UpcDict *self);
void **UpcDict_vals(UpcDict *self);

#endif
