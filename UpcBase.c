/*
 *  UpcBase.c
 *  UnipodC
 *
 *  Created by ionosph on 11/07/06.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcBase.h"
#include <time.h>

// エラーハンドラにデフォルトエラーハンドラをセット
static void UpcErrHandler_default(UpcErrStatus en, const char *msg, const char *filenm, int line);
void (*UpcErrHandler)(UpcErrStatus, const char *, const char *, int) = UpcErrHandler_default;

/* デフォルトのエラーハンドラ
 */
static void UpcErrHandler_default(UpcErrStatus en, const char *msg, const char *filenm, int line)
{
    if (en) {
        fprintf(stderr, "  Unipod Error: %s, %s (file \"%s\", line %d)\n", UpcErrStatus_toString(en), msg, filenm, line);
        exit(1);
    }
}

/* エラーハンドラをセットする.
 */
void upc_setErrHandler(void (*func)(UpcErrStatus, const char *, const char *, int))
{
    UpcErrHandler = func;
}

/* エラーステータスに対応する文字列を返す.
 */
const char *UpcErrStatus_toString(UpcErrStatus err)
{
    switch (err) {
        case UpcE_NoIssues:
            return "NoIssues";
        case UpcE_MemoryError:
            return "MemoryError";
        case UpcE_ArithmeticError:
            return "ArithmeticError";
        case UpcE_ZeroDivisionError:
            return "ZeroDivisionError";
        case UpcE_IndexError:
            return "IndexError";
        case UpcE_ValueError:
            return "ValueError";
        case UpcE_RuntimeError:
            return "RuntimeError";
        case UpcE_IOError:
            return "IOError";
        default:
            return "UnknownError";
    }
}

/* 文字列からブール値への変換.
 */
Bool upc_boolFromString(const char *s)
{
    if (!strcmp(s, "true") || !strcmp(s, "True") || !strcmp(s, "TRUE"))
        return TRUE;
    return FALSE;
}

/* double型の3次元配列m[0..(n3 - 1)][0..(n2 - 1)][0..(n1 - 1)]を確保して返す.
 */
double ***upc_mallocD3(int n3, int n2, int n1)
{
    double ***m = (double ***)malloc(sizeof(double **) * n3);
    double **m0 = (double **)malloc(sizeof(double *) * n3 * n2);
    double *m00 = (double *)malloc(sizeof(double) * n3 * n2 * n1);
    int i3, i2;
    
    if (!m || !m0 || !m00) {
        free(m);
        free(m0);
        free(m00);
        UpcERRHANDLER(UpcE_MemoryError, "in allocD3");
        return NULL;
    }
    for (i3 = 0; i3 < n3; i3++) {
        m[i3] = m0 + i3 * n2;
        for (i2 = 0; i2 < n2; i2++)
            m[i3][i2] = m00 + i3 * (n2 * n1) + i2 * n1;
    }
    return m;
}

/* 上の関数で作成した配列を解放する.
 */
void upc_freeD3(double ***m)
{
    if (m) {
        free(m[0][0]);
        free(m[0]);
        free(m);
    }
}

/* double型の2次元配列m[0..(n2 - 1)][0..(n1 - 1)]を確保して返す.
 */
double **upc_mallocD2(int n2, int n1)
{
    double **m = (double **)malloc(sizeof(double *) * n2);
    double *m0 = (double *)malloc(sizeof(double) * n2 * n1);
    int i;
    
    if (!m || !m0) {
        free(m);
        free(m0);
        UpcERRHANDLER(UpcE_MemoryError, "in allocD2");
        return NULL;
    }
    m[0] = m0;
    for (i = 1; i < n2; i++)
        m[i] = m[i - 1] + n1;
    return m;
}

/* double型の2次元配列m[0..(n2 - 1)][0..(n1 - 1)]を確保し, 0に初期化してて返す.
 */
double **upc_callocD2(int n2, int n1)
{
    double **m = (double **)malloc(sizeof(double *) * n2);
    double *m0 = (double *)calloc(n2 * n1, sizeof(double));
    int i;
    
    if (!m || !m0) {
        free(m);
        free(m0);
        UpcERRHANDLER(UpcE_MemoryError, "in callocD2");
        return NULL;
    }
    m[0] = m0;
    for (i = 1; i < n2; i++)
        m[i] = m[i - 1] + n1;
    return m;
}

/* 上の関数で作成した配列を解放する.
 */
void upc_freeD2(double **m)
{
    if (m) {
        free(m[0]);
        free(m);
    }
}

/* Bool型の配列m[0..(n2 - 1)][0..(n1 - 1)]を作成して返す.
 */
Bool **upc_mallocB2(int n2, int n1)
{
    Bool **m = (Bool **)malloc(sizeof(Bool *) * n2);
    Bool *m0 = (Bool *)malloc(sizeof(Bool) * n2 * n1);
    int i;
    
    if (!m || !m0) {
        free(m);
        free(m0);
        UpcERRHANDLER(UpcE_MemoryError, "in allocB2");
        return NULL;
    }
    m[0] = m0;
    for (i = 1; i < n2; i++)
        m[i] = m[i - 1] + n1;
    return m;
}

/* 上の関数で作成した配列を解放する.
 */
void upc_freeB2(Bool **m)
{
    if (m) {
        free(m[0]);
        free(m);
    }
}

/* ヒープソートを用いて, 配列 ra[1..n] を昇順(小さい順)に並べ換える.
 * (注) この関数は配列が添字1から始まることを想定している. a[0..n-1]をソートするには, a-1を与える必要がある.
 */
void upc_heapSort1(int n, double *ra)
{
    int i, ir, j, l;
    double rra;
    
    if (n < 2)
        return;
    l = (n >> 1) + 1;
    ir = n;
    while (1) {
        if (l > 1) {
            rra = ra[--l];
        }
        else {
            rra = ra[ir];
            ra[ir] = ra[1];
            if (--ir == 1) {
                ra[1] = rra;
                break;
            }
        }
        i = l;
        j = l + l;
        while (j <= ir) {
            if (j < ir && ra[j] < ra[j + 1])
                j++;
            if (rra < ra[j]) {
                ra[i] = ra[j];
                i = j;
                j <<= 1;
            }
            else
                break;
        }
        ra[i] = rra;
    }
}

/* ヒープソートを用いて, 配列 ra[1..n] を昇順(小さい順)に並べ換える.(添字が1から始まることに注意)
 * 配列 rb[1..n] も対応して並べ換える.
 * (注) この関数は配列が添字1から始まることを想定している. a[0..n-1], b[0..n-1]をソートするには, a-1, b-1を与える必要がある.
 */
void upc_heapSort2(int n, double *ra, double *rb)
{
    int i, ir, j, l;
    double rra, rrb;
    
    if (n < 2)
        return;
    l = (n >> 1) + 1;
    ir = n;
    while (1) {
        if (l > 1) {
            rra = ra[--l];
            rrb = rb[l];
        }
        else {
            rra = ra[ir];
            rrb = rb[ir];
            ra[ir] = ra[1];
            rb[ir] = rb[1];
            if (--ir == 1) {
                ra[1] = rra;
                rb[1] = rrb;
                break;
            }
        }
        i = l;
        j = l + l;
        while (j <= ir) {
            if (j < ir && ra[j] < ra[j + 1])
                j++;
            if (rra < ra[j]) {
                ra[i] = ra[j];
                rb[i] = rb[j];
                i = j;
                j <<= 1;
            }
            else
                break;
        }
        ra[i] = rra;
        rb[i] = rrb;
    }
}

/* 文字列をコピーして返す. 使い終わったら解放すること.
 */
char *upc_copy_string(const char *str)
{
    char *cpy = (char *)malloc(sizeof(char) * (strlen(str) + 1));
    
    if (!cpy) {
        UpcERRHANDLER(UpcE_MemoryError, "in upc_copy_string");
        return NULL;
    }
    strcpy(cpy, str);
    return cpy;
}

/* 文字列の比較. NULLでもよい版
 */
int upc_comp_string(const void *key1, const void *key2)
{
    char *s1 = (char *)key1;
    char *s2 = (char *)key2;
    
    if (s1)
        return s2 ? strcmp(s1, s2) : 1;
    else
        return s2 ? -1 : 0;
}

/* 長さnのスペース文字列を返す.
 */
const char *upc_space_string(int n)
{
    static const char *s = "                                                                                ";
    //                      01234567890123456789012345678901234567890123456789012345678901234567890123456789
    //                                1         2         3         4         5         6         7
    if (n > 80) {
        UpcERRHANDLER(UpcE_ValueError, "in upc_spaceString (n > 80)");
    }
    return s + 80 - n;
}

/* 四捨五入
 */
double upc_roundD(double x)
{
    if (x >= 0.0) {
        return floor(x + 0.5);
    }
    else {
        return -1.0 * floor(fabs(x) + 0.5);
    }
}

/* 2次元の点pを点qを中心に角度theta[radian]だけ回転する.
 * pは回転後の座標で上書きされる.
 */
void upc_rot2D(double *p, const double *q, double theta)
{
    const double cost = cos(theta);
    const double sint = sin(theta);
    const double dpX = p[X] - q[X];
    const double dpY = p[Y] - q[Y];
    
    p[X] = q[X] + dpX * cost - dpY * sint;
    p[Y] = q[Y] + dpX * sint + dpY * cost;
}

#ifdef UPC_USE_MY_HYPOT
// C99のhypotがない環境

static double my_hypot(double x, double y)
{
    double absx = fabs(x);
    double absy = fabs(y);

    if (absx > absy) {
        double r = absy / absx;

        return absx * sqrt(1.0 + r * r);
    }
    else if (absx < absy) {
        double r = absx / absy;

        return absy * sqrt(1.0 + r * r);
    }
    else { // absx == absy (== 0 の場合もあり得る)
        return sqrt(2.0) * absx;
    }
}

double (*upc_hypot)(double, double) = my_hypot;

#else
// C99のhypotがある環境

double (*upc_hypot)(double, double) = hypot;

#endif






/* 乱数の初期化
 */
void Upc_random_init(void)
{
    void init_genrand(unsigned long s);
    time_t t = time(NULL);
    
    init_genrand(t);
}

/* min <= x <= max に一様分布する乱数を返す.
 */
double Upc_random_uniform(double min, double max)
{
    double genrand_real1(void);
    double x = genrand_real1();
    
    return min + x * (max - min);
}



