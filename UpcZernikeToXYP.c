/*
 *  UpcZernikeToXYP.c
 *  UnipodC
 *
 *  Created by ionosph on 12/02/09.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcZernikeToXYP.h"

#define SWAP(x1,x2,t) do{(t)=(x1);(x1)=(x2);(x2)=(t);}while(0)
#define ZER_J(a,i) (a)[3*(i)]
#define ZER_K(a,i) (a)[3*(i)+1]
#define ZER_C(a,i) (a)[3*(i)+2]
// ZER_C = {1:sin項, 0:cos項} (軸対称項はcos項に含む)
#define POW_N1(i) (((i)%2)?-1:1)

static Bool NEED_INIT = 1;
static const int CNO_MAX = 400;
static double **BCOEF = NULL;
static int *ZERJKC = NULL;
static UpcQMatrix **ZCOEFMATRIX = NULL;


static void UpcQMatrix_dealloc(UpcQMatrix *self)
{
    upc_freeD2(self->data);
    free(self);
}

UpcQMatrix *UpcQMatrix_init(int ixMax, int iyMax, int xIsOdd, int yIsOdd)
{
    UpcQMatrix *self = (UpcQMatrix *)malloc(sizeof(UpcQMatrix));

    if (self) {
        int nx = (xIsOdd) ? (ixMax + 1) / 2 : ixMax / 2 + 1;
        int ny = (yIsOdd) ? (iyMax + 1) / 2 : iyMax / 2 + 1;
        
        self->_refCount = 1;
        self->_dealloc = UpcQMatrix_dealloc;
        self->ixMax = ixMax;
        self->iyMax = iyMax;
        self->xIsOdd = (xIsOdd) ? 1 : 0;
        self->yIsOdd = (yIsOdd) ? 1 : 0;
        self->data = upc_callocD2(ny, nx);
        if (!(self->data)) {
            free(self);
            UpcERRHANDLER(UpcE_MemoryError, "in UpcQMatrix_init");
            return NULL;
        }
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcQMatrix_init");
    }
    return self;
}

/* 指定index(ix, iy)が有効かチェックする.
 */
static Bool UpcQMatrix_checkIndex(UpcQMatrix *self, int iy, int ix)
{
    if (ix >= 0 &&
        ix <= self->ixMax &&
        iy >= 0 &&
        iy <= self->iyMax &&
        ix % 2 == self->xIsOdd &&
        iy % 2 == self->yIsOdd)
        return TRUE;
    return FALSE;
}

/* 指定index(ix, iy)のデータを返す. 指定indexが無効なら0.0を返す.
 */
double UpcQMatrix_data(UpcQMatrix *self, int iy, int ix)
{
    return (UpcQMatrix_checkIndex(self, iy, ix)) ? self->data[iy / 2][ix / 2] : 0.0;
}

UpcErrStatus UpcQMatrix_setData(UpcQMatrix *self, int iy, int ix, double value)
{
    if (UpcQMatrix_checkIndex(self, iy, ix)) {
        self->data[iy / 2][ix / 2] = value;
        return UpcE_NoIssues;
    }
    UpcERRHANDLER(UpcE_IndexError, "in UpcQMatrix_setData");
    return UpcE_IndexError;
}

UpcErrStatus UpcQMatrix_addData(UpcQMatrix *self, int iy, int ix, double value)
{
    if (UpcQMatrix_checkIndex(self, iy, ix)) {
        self->data[iy / 2][ix / 2] += value;
        return UpcE_NoIssues;
    }
    UpcERRHANDLER(UpcE_IndexError, "in UpcQMatrix_addData");
    return UpcE_IndexError;
}

UpcErrStatus UpcQMatrix_mulData(UpcQMatrix *self, int iy, int ix, double value)
{
    if (UpcQMatrix_checkIndex(self, iy, ix)) {
        self->data[iy / 2][ix / 2] *= value;
        return UpcE_NoIssues;
    }
    UpcERRHANDLER(UpcE_IndexError, "in UpcQMatrix_mulData");
    return UpcE_IndexError;
}

void UpcQMatrix_print(UpcQMatrix *self, FILE *ostream)
{
    int iy, ix;
    double z;

    fprintf(ostream, "   ix,  iy,    data\n");
    for (ix = self->xIsOdd; ix <= self->ixMax; ix += 2) {
        for (iy = self->yIsOdd; iy <= self->iyMax; iy += 2) {
            z = UpcQMatrix_data(self, iy, ix);
            if (z)
                fprintf(ostream, "  %3d, %3d, %23.15e\n", ix, iy, z);
        }
    }
}

//----------------------------------------------------------------------

/* Zernike番号 → m, n, c
 * m, n : Zernike関数のパラメータ
 *        (mはθ依存性, nは動径多項式の最大次数)
 * c : cos項(軸対称項含), sin項 を表す項
 */
/*
static void cno_to_mnc(int izer, int *m, int *n, int *c)
{
    int j = (int)ceil(sqrt(izer));
    int l = j * j - izer;

    *m = (l + 1) / 2;
    *n = 2 * (j - 1) - *m;
    *c = (l % 2); // 1:sin項, 0:cos項 (軸対称項はcos項に含む)
}
*/

/* n次までの2項係数の表を作る.
 * a[n][k] = nCk (n >= k >= 0)
 */
static double **makeBinomialCoefMatrix(int n)
{
    int i, j;
    double **a = (double **)malloc(sizeof(double *) * (n + 1));
    double *ai = (double *)malloc(sizeof(double) * (((n + 2) * (n + 1)) / 2));
    double *t;

    if (!a || !ai) {
        free(a);
        free(ai);
        UpcERRHANDLER(UpcE_MemoryError, "in makeBinomialCoefMatrix");
        return NULL;
    }
    a[0] = ai;
    a[0][0] = 1;
    for (i = 1; i <= n; i++) {
        t = a[i] = a[i - 1] + i;
        t[0] = t[i] = 1;
        for (j = 1; j < i; j++)
            t[j] = a[i - 1][j - 1] + a[i - 1][j];
    }
    return a;
}

/* Zernike第cno項までのj, k, cのリストを作る.
 */
static int *makeJKCArray(int cno)
{
    int i, A, B;
    int *a = (int *)malloc(sizeof(int) * (3 * (cno + 1)));

    if (!a) {
        UpcERRHANDLER(UpcE_MemoryError, "in makeJKCArray");
        return NULL;
    }
    ZER_J(a, 0) = ZER_K(a, 0) = 0;
    for (i = 1; i <= cno; i++) {
        A = (int)(ceil(sqrt(i)));
        B = A * A - i;
        ZER_J(a, i) = A - 1;
        ZER_K(a, i) = ZER_J(a, i) - (B + 1) / 2;
        ZER_C(a, i) = B % 2;
    }
    return a;
}

static UpcErrStatus initParams(void)
{
    NEED_INIT = FALSE;
    ZERJKC = makeJKCArray(CNO_MAX);
    if (!ZERJKC) {
        UpcERRHANDLER(UpcE_MemoryError, "in initParams");
        return UpcE_MemoryError;
    }
    BCOEF = makeBinomialCoefMatrix(ZER_J(ZERJKC, CNO_MAX) + ZER_K(ZERJKC, CNO_MAX));
    if (!BCOEF) {
        free(ZERJKC);
        UpcERRHANDLER(UpcE_MemoryError, "in initParams");
        return UpcE_MemoryError;
    }
    ZCOEFMATRIX = (UpcQMatrix **)calloc(CNO_MAX + 1, sizeof(UpcQMatrix *));
    if (!ZCOEFMATRIX) {
        free(ZERJKC);
        free(BCOEF[0]);
        free(BCOEF);
        UpcERRHANDLER(UpcE_MemoryError, "in initParams");
        return UpcE_MemoryError;
    }
    return UpcE_NoIssues;
}

//----------------------------------------------------------------------

/* Zernike第cno項のXYべき表現を構築して返す.
 */
static UpcQMatrix *makeZCoefXY(int cno)
{
    int j, k, c;
    int s, p, q;
    int js, ks, jsks, ip, iq, pq;
    double coefA, coefB, coefC;
    UpcQMatrix *coefMatrix;

    j = ZER_J(ZERJKC, cno);
    k = ZER_K(ZERJKC, cno);
    c = ZER_C(ZERJKC, cno);
    coefMatrix = UpcQMatrix_init(j + k, j + k, (j + k - c) % 2, c);
    if (!coefMatrix) {
        UpcERRHANDLER(UpcE_MemoryError, "in makeZCoefXY");
        return NULL;
    }

    for (s = 0; s <= k; s++) {
        js = j - s;
        ks = k - s;
        jsks = js + ks;
        for (ip = 0; ip <= js; ip++) {
            p = (ip % 2) ? js - ip / 2 : ip / 2;
            coefB = BCOEF[js][p];
            for (iq = 0; iq <= ks; iq++) {
                q = (iq % 2) ? ks - iq / 2 : iq / 2;
                pq = p + q;
                if (pq % 2 == c) {
                    coefC = POW_N1(q) * POW_N1(pq / 2) * BCOEF[ks][q];
                    UpcQMatrix_addData(coefMatrix, pq, jsks - pq, coefB * coefC);
                }
            }
        }
        coefA = POW_N1(s) * BCOEF[j + ks][s] * BCOEF[jsks][js];
        for (pq = c; pq <= jsks; pq += 2) {
            UpcQMatrix_mulData(coefMatrix, pq, jsks - pq, coefA);
        }
    }
    return coefMatrix;
}

/* Zernike第cno項のXYべき表現を返す.
 * この関数が返したオブジェクトは借りているだけ. 解放してはならない.
 */
UpcQMatrix *UpcZernikeToXYP_zcoefxy(int cno)
{
    if (cno < 1 || cno > CNO_MAX) {
        UpcERRHANDLER(UpcE_ValueError, "cno must be 1 <= cno <= 400 (in UpcQMatrix_zcoefxy)");
    }

    if (NEED_INIT)
        initParams();

    if (!ZCOEFMATRIX[cno]) {
        ZCOEFMATRIX[cno] = makeZCoefXY(cno);
        if (!ZCOEFMATRIX[cno]) {
            UpcERRHANDLER(UpcE_MemoryError, "in UpcQMatrix_zcoefxy");
            return NULL;
        }
    }
    return ZCOEFMATRIX[cno];
}

/* Zernike第cno項のXYべき表現の最大次数を返す.
 */
int UpcZernikeToXYP_maxorder(int cno)
{
    if (cno < 1 || cno > CNO_MAX) {
        UpcERRHANDLER(UpcE_ValueError, "cno must be 1 <= cno <= 400 (in UpcQMatrix_zcoefxy)");
    }

    if (NEED_INIT)
        initParams();
    
    return ZER_J(ZERJKC, cno) + ZER_K(ZERJKC, cno);
}

//----------------------------------------------------------------------

Bool testUpcZernikeToXYP(void)
{
    int i, j;
    double **bcoef;
    int *zjka;
    int k, s;
    double coef;
    UpcQMatrix *zcoef;

    if (FALSE) {
        bcoef = makeBinomialCoefMatrix(38);
        for (i = 0; i <= 38; i++) {
            for (j = 0; j <= i; j++)
                printf(",%.6f", bcoef[i][j]);
            printf("\n");
        }
    }

    if (FALSE) {
        zjka = makeJKCArray(400);
        for (i = 1; i <= 400; i++) {
            j = ZER_J(zjka, i);
            k = ZER_K(zjka, i);
            printf("%3d,%2d,%2d", i, j, k);
            for (s = 0; s <= k; s++) {
                coef = POW_N1(s) * bcoef[j + k - s][s] * bcoef[j + k - 2 * s][j - s];
                printf(",%.6f", coef);
            }
            printf("\n");
        }
    }

    for (i = 1; i <= 36; i++) {
        zcoef = UpcZernikeToXYP_zcoefxy(i);
        printf("C%03d, %d\n", i, UpcZernikeToXYP_maxorder(i));
        UpcQMatrix_print(zcoef, stdout);
    }
    return TRUE;
}







