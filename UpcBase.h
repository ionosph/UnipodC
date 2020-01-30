/*
 *  UpcBase.h
 *  UnipodC
 *
 *  Created by ionosph on 2011/04/28.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _UPCBASE_H
#define _UPCBASE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// 次のマクロはC99のhypotが使えない環境で定義する. 自作のhypotが使えるようになる.
//#define UPC_USE_MY_HYPOT

typedef int Bool;
static const int TRUE  = 1;
static const int FALSE = 0;

static const int X = 0;
static const int Y = 1;
static const int Z = 2;

static const int PRAY = 0;
static const int YUPP = 1;
static const int YLOW = 2;
static const int XUPP = 3;
static const int XLOW = 4;

static const double PI      =         3.1415926535897932;
static const double DEG2RAD =         3.1415926535897932 / 180.0;
static const double RAD2DEG = 180.0 / 3.1415926535897932;

#define RADIANS(x) ((x)*DEG2RAD)
#define DEGREES(x) ((x)*RAD2DEG)

#define UpcObj_RETAIN(obj) (++((obj)->_refCount))
#define UpcObj_RELEASE(obj) do{if((--((obj)->_refCount))==0){(obj)->_dealloc(obj);(obj)=NULL;}}while(0)
#define BOOLTOSTRING(i) ((i)?"true":"false")

typedef enum _UpcErrStatus {
    UpcE_NoIssues          =  0,
    UpcE_MemoryError       = -1,
    UpcE_ArithmeticError   = -2,
    UpcE_ZeroDivisionError = -3,
    UpcE_IndexError        = -4,
    UpcE_ValueError        = -5,
    UpcE_RuntimeError      = -6,
    UpcE_IOError           = -7
} UpcErrStatus;

#define UpcERRHANDLER(err,msg) UpcErrHandler((err),(msg),__FILE__,__LINE__)

extern void (*UpcErrHandler)(UpcErrStatus, const char *, const char *, int);
extern double (*upc_hypot)(double, double);

const char *UpcErrStatus_toString(UpcErrStatus err);
void upc_setErrHandler(void (*func)(UpcErrStatus, const char *, const char *, int));

Bool upc_boolFromString(const char *s);

double ***upc_mallocD3(int n3, int n2, int n1);
void upc_freeD3(double ***m);
double **upc_mallocD2(int n2, int n1);
double **upc_callocD2(int n2, int n1);
void upc_freeD2(double **m);
Bool **upc_mallocB2(int n2, int n1);
void upc_freeB2(Bool **m);
void upc_heapSort1(int n, double *ra);
void upc_heapSort2(int n, double *ra, double *rb);

char *upc_copy_string(const char *str);
int upc_comp_string(const void *key1, const void *key2);

const char *upc_space_string(int n);

double upc_roundD(double x);

void upc_rot2D(double *p, const double *q, double theta);

void Upc_random_init(void);
double Upc_random_uniform(double min, double max);

#endif
