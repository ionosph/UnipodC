/*------------------------------------------------------------------------------
    zernike(400項までのZernike関数を扱うモジュール)
------------------------------------------------------------------------------*/
#include "UpcBase.h"

// 中間変数を格納する構造体
typedef struct {
    double rf;
    double tf[19 + 1][2]; // r^n cos(nθ), r^n sin(nθ)
    int prev_term;
    long double r2;
} mparams;

// r^n cos(nθ), r^n sin(nθ) を指定次数まで計算
static void cal_tf(double x, double y, int iend, mparams *p)
{
    double c1 = x;
    double s1 = y;
    double ct, st;
    int i;

    p->tf[1][0] = s1;
    p->tf[1][1] = c1;
    for (i = 2; i <= iend; i++) {
        st = y * c1 + x * s1;
        ct = x * c1 - y * s1;
        s1 = st;
        c1 = ct;
        p->tf[i][0] = s1;
        p->tf[i][1] = c1;
    }
}

static double zer000(mparams *p)
{
    return 0.0;
}

static double zer001(mparams *p)
{
    return 1.0;
}

static double zer002(mparams *p)
{
    return p->tf[1][1];
}

static double zer003(mparams *p)
{
    return p->tf[1][0];
}

static double zer004(mparams *p)
{
    long double r2 = p->r2;
    long double rf =              2.0;
    rf = rf * r2 +               -1.0;
    return rf;
}

static double zer005(mparams *p)
{
    return p->tf[2][1];
}

static double zer006(mparams *p)
{
    return p->tf[2][0];
}

static double zer007(mparams *p)
{
    long double r2 = p->r2;
    long double rf =              3.0;
    rf = rf * r2 +               -2.0;
    p->rf = (double)rf;
    p->prev_term = 7;
    return p->rf * p->tf[1][1];
}

static double zer008(mparams *p)
{
    if (p->prev_term != 7)  zer007(p);
    return p->rf * p->tf[1][0];
}

static double zer009(mparams *p)
{
    long double r2 = p->r2;
    long double rf =              6.0;
    rf = rf * r2 +               -6.0;
    rf = rf * r2 +                1.0;
    return rf;
}

static double zer010(mparams *p)
{
    return p->tf[3][1];
}

static double zer011(mparams *p)
{
    return p->tf[3][0];
}

static double zer012(mparams *p)
{
    long double r2 = p->r2;
    long double rf =              4.0;
    rf = rf * r2 +               -3.0;
    p->rf = (double)rf;
    p->prev_term = 12;
    return p->rf * p->tf[2][1];
}

static double zer013(mparams *p)
{
    if (p->prev_term != 12)  zer012(p);
    return p->rf * p->tf[2][0];
}

static double zer014(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             10.0;
    rf = rf * r2 +              -12.0;
    rf = rf * r2 +                3.0;
    p->rf = (double)rf;
    p->prev_term = 14;
    return p->rf * p->tf[1][1];
}

static double zer015(mparams *p)
{
    if (p->prev_term != 14)  zer014(p);
    return p->rf * p->tf[1][0];
}

static double zer016(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             20.0;
    rf = rf * r2 +              -30.0;
    rf = rf * r2 +               12.0;
    rf = rf * r2 +               -1.0;
    return rf;
}

static double zer017(mparams *p)
{
    return p->tf[4][1];
}

static double zer018(mparams *p)
{
    return p->tf[4][0];
}

static double zer019(mparams *p)
{
    long double r2 = p->r2;
    long double rf =              5.0;
    rf = rf * r2 +               -4.0;
    p->rf = (double)rf;
    p->prev_term = 19;
    return p->rf * p->tf[3][1];
}

static double zer020(mparams *p)
{
    if (p->prev_term != 19)  zer019(p);
    return p->rf * p->tf[3][0];
}

static double zer021(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             15.0;
    rf = rf * r2 +              -20.0;
    rf = rf * r2 +                6.0;
    p->rf = (double)rf;
    p->prev_term = 21;
    return p->rf * p->tf[2][1];
}

static double zer022(mparams *p)
{
    if (p->prev_term != 21)  zer021(p);
    return p->rf * p->tf[2][0];
}

static double zer023(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             35.0;
    rf = rf * r2 +              -60.0;
    rf = rf * r2 +               30.0;
    rf = rf * r2 +               -4.0;
    p->rf = (double)rf;
    p->prev_term = 23;
    return p->rf * p->tf[1][1];
}

static double zer024(mparams *p)
{
    if (p->prev_term != 23)  zer023(p);
    return p->rf * p->tf[1][0];
}

static double zer025(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             70.0;
    rf = rf * r2 +             -140.0;
    rf = rf * r2 +               90.0;
    rf = rf * r2 +              -20.0;
    rf = rf * r2 +                1.0;
    return rf;
}

static double zer026(mparams *p)
{
    return p->tf[5][1];
}

static double zer027(mparams *p)
{
    return p->tf[5][0];
}

static double zer028(mparams *p)
{
    long double r2 = p->r2;
    long double rf =              6.0;
    rf = rf * r2 +               -5.0;
    p->rf = (double)rf;
    p->prev_term = 28;
    return p->rf * p->tf[4][1];
}

static double zer029(mparams *p)
{
    if (p->prev_term != 28)  zer028(p);
    return p->rf * p->tf[4][0];
}

static double zer030(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             21.0;
    rf = rf * r2 +              -30.0;
    rf = rf * r2 +               10.0;
    p->rf = (double)rf;
    p->prev_term = 30;
    return p->rf * p->tf[3][1];
}

static double zer031(mparams *p)
{
    if (p->prev_term != 30)  zer030(p);
    return p->rf * p->tf[3][0];
}

static double zer032(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             56.0;
    rf = rf * r2 +             -105.0;
    rf = rf * r2 +               60.0;
    rf = rf * r2 +              -10.0;
    p->rf = (double)rf;
    p->prev_term = 32;
    return p->rf * p->tf[2][1];
}

static double zer033(mparams *p)
{
    if (p->prev_term != 32)  zer032(p);
    return p->rf * p->tf[2][0];
}

static double zer034(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            126.0;
    rf = rf * r2 +             -280.0;
    rf = rf * r2 +              210.0;
    rf = rf * r2 +              -60.0;
    rf = rf * r2 +                5.0;
    p->rf = (double)rf;
    p->prev_term = 34;
    return p->rf * p->tf[1][1];
}

static double zer035(mparams *p)
{
    if (p->prev_term != 34)  zer034(p);
    return p->rf * p->tf[1][0];
}

static double zer036(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            252.0;
    rf = rf * r2 +             -630.0;
    rf = rf * r2 +              560.0;
    rf = rf * r2 +             -210.0;
    rf = rf * r2 +               30.0;
    rf = rf * r2 +               -1.0;
    return rf;
}

static double zer037(mparams *p)
{
    return p->tf[6][1];
}

static double zer038(mparams *p)
{
    return p->tf[6][0];
}

static double zer039(mparams *p)
{
    long double r2 = p->r2;
    long double rf =              7.0;
    rf = rf * r2 +               -6.0;
    p->rf = (double)rf;
    p->prev_term = 39;
    return p->rf * p->tf[5][1];
}

static double zer040(mparams *p)
{
    if (p->prev_term != 39)  zer039(p);
    return p->rf * p->tf[5][0];
}

static double zer041(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             28.0;
    rf = rf * r2 +              -42.0;
    rf = rf * r2 +               15.0;
    p->rf = (double)rf;
    p->prev_term = 41;
    return p->rf * p->tf[4][1];
}

static double zer042(mparams *p)
{
    if (p->prev_term != 41)  zer041(p);
    return p->rf * p->tf[4][0];
}

static double zer043(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             84.0;
    rf = rf * r2 +             -168.0;
    rf = rf * r2 +              105.0;
    rf = rf * r2 +              -20.0;
    p->rf = (double)rf;
    p->prev_term = 43;
    return p->rf * p->tf[3][1];
}

static double zer044(mparams *p)
{
    if (p->prev_term != 43)  zer043(p);
    return p->rf * p->tf[3][0];
}

static double zer045(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            210.0;
    rf = rf * r2 +             -504.0;
    rf = rf * r2 +              420.0;
    rf = rf * r2 +             -140.0;
    rf = rf * r2 +               15.0;
    p->rf = (double)rf;
    p->prev_term = 45;
    return p->rf * p->tf[2][1];
}

static double zer046(mparams *p)
{
    if (p->prev_term != 45)  zer045(p);
    return p->rf * p->tf[2][0];
}

static double zer047(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            462.0;
    rf = rf * r2 +            -1260.0;
    rf = rf * r2 +             1260.0;
    rf = rf * r2 +             -560.0;
    rf = rf * r2 +              105.0;
    rf = rf * r2 +               -6.0;
    p->rf = (double)rf;
    p->prev_term = 47;
    return p->rf * p->tf[1][1];
}

static double zer048(mparams *p)
{
    if (p->prev_term != 47)  zer047(p);
    return p->rf * p->tf[1][0];
}

static double zer049(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            924.0;
    rf = rf * r2 +            -2772.0;
    rf = rf * r2 +             3150.0;
    rf = rf * r2 +            -1680.0;
    rf = rf * r2 +              420.0;
    rf = rf * r2 +              -42.0;
    rf = rf * r2 +                1.0;
    return rf;
}

static double zer050(mparams *p)
{
    return p->tf[7][1];
}

static double zer051(mparams *p)
{
    return p->tf[7][0];
}

static double zer052(mparams *p)
{
    long double r2 = p->r2;
    long double rf =              8.0;
    rf = rf * r2 +               -7.0;
    p->rf = (double)rf;
    p->prev_term = 52;
    return p->rf * p->tf[6][1];
}

static double zer053(mparams *p)
{
    if (p->prev_term != 52)  zer052(p);
    return p->rf * p->tf[6][0];
}

static double zer054(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             36.0;
    rf = rf * r2 +              -56.0;
    rf = rf * r2 +               21.0;
    p->rf = (double)rf;
    p->prev_term = 54;
    return p->rf * p->tf[5][1];
}

static double zer055(mparams *p)
{
    if (p->prev_term != 54)  zer054(p);
    return p->rf * p->tf[5][0];
}

static double zer056(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            120.0;
    rf = rf * r2 +             -252.0;
    rf = rf * r2 +              168.0;
    rf = rf * r2 +              -35.0;
    p->rf = (double)rf;
    p->prev_term = 56;
    return p->rf * p->tf[4][1];
}

static double zer057(mparams *p)
{
    if (p->prev_term != 56)  zer056(p);
    return p->rf * p->tf[4][0];
}

static double zer058(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            330.0;
    rf = rf * r2 +             -840.0;
    rf = rf * r2 +              756.0;
    rf = rf * r2 +             -280.0;
    rf = rf * r2 +               35.0;
    p->rf = (double)rf;
    p->prev_term = 58;
    return p->rf * p->tf[3][1];
}

static double zer059(mparams *p)
{
    if (p->prev_term != 58)  zer058(p);
    return p->rf * p->tf[3][0];
}

static double zer060(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            792.0;
    rf = rf * r2 +            -2310.0;
    rf = rf * r2 +             2520.0;
    rf = rf * r2 +            -1260.0;
    rf = rf * r2 +              280.0;
    rf = rf * r2 +              -21.0;
    p->rf = (double)rf;
    p->prev_term = 60;
    return p->rf * p->tf[2][1];
}

static double zer061(mparams *p)
{
    if (p->prev_term != 60)  zer060(p);
    return p->rf * p->tf[2][0];
}

static double zer062(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           1716.0;
    rf = rf * r2 +            -5544.0;
    rf = rf * r2 +             6930.0;
    rf = rf * r2 +            -4200.0;
    rf = rf * r2 +             1260.0;
    rf = rf * r2 +             -168.0;
    rf = rf * r2 +                7.0;
    p->rf = (double)rf;
    p->prev_term = 62;
    return p->rf * p->tf[1][1];
}

static double zer063(mparams *p)
{
    if (p->prev_term != 62)  zer062(p);
    return p->rf * p->tf[1][0];
}

static double zer064(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           3432.0;
    rf = rf * r2 +           -12012.0;
    rf = rf * r2 +            16632.0;
    rf = rf * r2 +           -11550.0;
    rf = rf * r2 +             4200.0;
    rf = rf * r2 +             -756.0;
    rf = rf * r2 +               56.0;
    rf = rf * r2 +               -1.0;
    return rf;
}

static double zer065(mparams *p)
{
    return p->tf[8][1];
}

static double zer066(mparams *p)
{
    return p->tf[8][0];
}

static double zer067(mparams *p)
{
    long double r2 = p->r2;
    long double rf =              9.0;
    rf = rf * r2 +               -8.0;
    p->rf = (double)rf;
    p->prev_term = 67;
    return p->rf * p->tf[7][1];
}

static double zer068(mparams *p)
{
    if (p->prev_term != 67)  zer067(p);
    return p->rf * p->tf[7][0];
}

static double zer069(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             45.0;
    rf = rf * r2 +              -72.0;
    rf = rf * r2 +               28.0;
    p->rf = (double)rf;
    p->prev_term = 69;
    return p->rf * p->tf[6][1];
}

static double zer070(mparams *p)
{
    if (p->prev_term != 69)  zer069(p);
    return p->rf * p->tf[6][0];
}

static double zer071(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            165.0;
    rf = rf * r2 +             -360.0;
    rf = rf * r2 +              252.0;
    rf = rf * r2 +              -56.0;
    p->rf = (double)rf;
    p->prev_term = 71;
    return p->rf * p->tf[5][1];
}

static double zer072(mparams *p)
{
    if (p->prev_term != 71)  zer071(p);
    return p->rf * p->tf[5][0];
}

static double zer073(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            495.0;
    rf = rf * r2 +            -1320.0;
    rf = rf * r2 +             1260.0;
    rf = rf * r2 +             -504.0;
    rf = rf * r2 +               70.0;
    p->rf = (double)rf;
    p->prev_term = 73;
    return p->rf * p->tf[4][1];
}

static double zer074(mparams *p)
{
    if (p->prev_term != 73)  zer073(p);
    return p->rf * p->tf[4][0];
}

static double zer075(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           1287.0;
    rf = rf * r2 +            -3960.0;
    rf = rf * r2 +             4620.0;
    rf = rf * r2 +            -2520.0;
    rf = rf * r2 +              630.0;
    rf = rf * r2 +              -56.0;
    p->rf = (double)rf;
    p->prev_term = 75;
    return p->rf * p->tf[3][1];
}

static double zer076(mparams *p)
{
    if (p->prev_term != 75)  zer075(p);
    return p->rf * p->tf[3][0];
}

static double zer077(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           3003.0;
    rf = rf * r2 +           -10296.0;
    rf = rf * r2 +            13860.0;
    rf = rf * r2 +            -9240.0;
    rf = rf * r2 +             3150.0;
    rf = rf * r2 +             -504.0;
    rf = rf * r2 +               28.0;
    p->rf = (double)rf;
    p->prev_term = 77;
    return p->rf * p->tf[2][1];
}

static double zer078(mparams *p)
{
    if (p->prev_term != 77)  zer077(p);
    return p->rf * p->tf[2][0];
}

static double zer079(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           6435.0;
    rf = rf * r2 +           -24024.0;
    rf = rf * r2 +            36036.0;
    rf = rf * r2 +           -27720.0;
    rf = rf * r2 +            11550.0;
    rf = rf * r2 +            -2520.0;
    rf = rf * r2 +              252.0;
    rf = rf * r2 +               -8.0;
    p->rf = (double)rf;
    p->prev_term = 79;
    return p->rf * p->tf[1][1];
}

static double zer080(mparams *p)
{
    if (p->prev_term != 79)  zer079(p);
    return p->rf * p->tf[1][0];
}

static double zer081(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          12870.0;
    rf = rf * r2 +           -51480.0;
    rf = rf * r2 +            84084.0;
    rf = rf * r2 +           -72072.0;
    rf = rf * r2 +            34650.0;
    rf = rf * r2 +            -9240.0;
    rf = rf * r2 +             1260.0;
    rf = rf * r2 +              -72.0;
    rf = rf * r2 +                1.0;
    return rf;
}

static double zer082(mparams *p)
{
    return p->tf[9][1];
}

static double zer083(mparams *p)
{
    return p->tf[9][0];
}

static double zer084(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             10.0;
    rf = rf * r2 +               -9.0;
    p->rf = (double)rf;
    p->prev_term = 84;
    return p->rf * p->tf[8][1];
}

static double zer085(mparams *p)
{
    if (p->prev_term != 84)  zer084(p);
    return p->rf * p->tf[8][0];
}

static double zer086(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             55.0;
    rf = rf * r2 +              -90.0;
    rf = rf * r2 +               36.0;
    p->rf = (double)rf;
    p->prev_term = 86;
    return p->rf * p->tf[7][1];
}

static double zer087(mparams *p)
{
    if (p->prev_term != 86)  zer086(p);
    return p->rf * p->tf[7][0];
}

static double zer088(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            220.0;
    rf = rf * r2 +             -495.0;
    rf = rf * r2 +              360.0;
    rf = rf * r2 +              -84.0;
    p->rf = (double)rf;
    p->prev_term = 88;
    return p->rf * p->tf[6][1];
}

static double zer089(mparams *p)
{
    if (p->prev_term != 88)  zer088(p);
    return p->rf * p->tf[6][0];
}

static double zer090(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            715.0;
    rf = rf * r2 +            -1980.0;
    rf = rf * r2 +             1980.0;
    rf = rf * r2 +             -840.0;
    rf = rf * r2 +              126.0;
    p->rf = (double)rf;
    p->prev_term = 90;
    return p->rf * p->tf[5][1];
}

static double zer091(mparams *p)
{
    if (p->prev_term != 90)  zer090(p);
    return p->rf * p->tf[5][0];
}

static double zer092(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           2002.0;
    rf = rf * r2 +            -6435.0;
    rf = rf * r2 +             7920.0;
    rf = rf * r2 +            -4620.0;
    rf = rf * r2 +             1260.0;
    rf = rf * r2 +             -126.0;
    p->rf = (double)rf;
    p->prev_term = 92;
    return p->rf * p->tf[4][1];
}

static double zer093(mparams *p)
{
    if (p->prev_term != 92)  zer092(p);
    return p->rf * p->tf[4][0];
}

static double zer094(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           5005.0;
    rf = rf * r2 +           -18018.0;
    rf = rf * r2 +            25740.0;
    rf = rf * r2 +           -18480.0;
    rf = rf * r2 +             6930.0;
    rf = rf * r2 +            -1260.0;
    rf = rf * r2 +               84.0;
    p->rf = (double)rf;
    p->prev_term = 94;
    return p->rf * p->tf[3][1];
}

static double zer095(mparams *p)
{
    if (p->prev_term != 94)  zer094(p);
    return p->rf * p->tf[3][0];
}

static double zer096(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          11440.0;
    rf = rf * r2 +           -45045.0;
    rf = rf * r2 +            72072.0;
    rf = rf * r2 +           -60060.0;
    rf = rf * r2 +            27720.0;
    rf = rf * r2 +            -6930.0;
    rf = rf * r2 +              840.0;
    rf = rf * r2 +              -36.0;
    p->rf = (double)rf;
    p->prev_term = 96;
    return p->rf * p->tf[2][1];
}

static double zer097(mparams *p)
{
    if (p->prev_term != 96)  zer096(p);
    return p->rf * p->tf[2][0];
}

static double zer098(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          24310.0;
    rf = rf * r2 +          -102960.0;
    rf = rf * r2 +           180180.0;
    rf = rf * r2 +          -168168.0;
    rf = rf * r2 +            90090.0;
    rf = rf * r2 +           -27720.0;
    rf = rf * r2 +             4620.0;
    rf = rf * r2 +             -360.0;
    rf = rf * r2 +                9.0;
    p->rf = (double)rf;
    p->prev_term = 98;
    return p->rf * p->tf[1][1];
}

static double zer099(mparams *p)
{
    if (p->prev_term != 98)  zer098(p);
    return p->rf * p->tf[1][0];
}

static double zer100(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          48620.0;
    rf = rf * r2 +          -218790.0;
    rf = rf * r2 +           411840.0;
    rf = rf * r2 +          -420420.0;
    rf = rf * r2 +           252252.0;
    rf = rf * r2 +           -90090.0;
    rf = rf * r2 +            18480.0;
    rf = rf * r2 +            -1980.0;
    rf = rf * r2 +               90.0;
    rf = rf * r2 +               -1.0;
    return rf;
}

static double zer101(mparams *p)
{
    return p->tf[10][1];
}

static double zer102(mparams *p)
{
    return p->tf[10][0];
}

static double zer103(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             11.0;
    rf = rf * r2 +              -10.0;
    p->rf = (double)rf;
    p->prev_term = 103;
    return p->rf * p->tf[9][1];
}

static double zer104(mparams *p)
{
    if (p->prev_term != 103)  zer103(p);
    return p->rf * p->tf[9][0];
}

static double zer105(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             66.0;
    rf = rf * r2 +             -110.0;
    rf = rf * r2 +               45.0;
    p->rf = (double)rf;
    p->prev_term = 105;
    return p->rf * p->tf[8][1];
}

static double zer106(mparams *p)
{
    if (p->prev_term != 105)  zer105(p);
    return p->rf * p->tf[8][0];
}

static double zer107(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            286.0;
    rf = rf * r2 +             -660.0;
    rf = rf * r2 +              495.0;
    rf = rf * r2 +             -120.0;
    p->rf = (double)rf;
    p->prev_term = 107;
    return p->rf * p->tf[7][1];
}

static double zer108(mparams *p)
{
    if (p->prev_term != 107)  zer107(p);
    return p->rf * p->tf[7][0];
}

static double zer109(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           1001.0;
    rf = rf * r2 +            -2860.0;
    rf = rf * r2 +             2970.0;
    rf = rf * r2 +            -1320.0;
    rf = rf * r2 +              210.0;
    p->rf = (double)rf;
    p->prev_term = 109;
    return p->rf * p->tf[6][1];
}

static double zer110(mparams *p)
{
    if (p->prev_term != 109)  zer109(p);
    return p->rf * p->tf[6][0];
}

static double zer111(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           3003.0;
    rf = rf * r2 +           -10010.0;
    rf = rf * r2 +            12870.0;
    rf = rf * r2 +            -7920.0;
    rf = rf * r2 +             2310.0;
    rf = rf * r2 +             -252.0;
    p->rf = (double)rf;
    p->prev_term = 111;
    return p->rf * p->tf[5][1];
}

static double zer112(mparams *p)
{
    if (p->prev_term != 111)  zer111(p);
    return p->rf * p->tf[5][0];
}

static double zer113(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           8008.0;
    rf = rf * r2 +           -30030.0;
    rf = rf * r2 +            45045.0;
    rf = rf * r2 +           -34320.0;
    rf = rf * r2 +            13860.0;
    rf = rf * r2 +            -2772.0;
    rf = rf * r2 +              210.0;
    p->rf = (double)rf;
    p->prev_term = 113;
    return p->rf * p->tf[4][1];
}

static double zer114(mparams *p)
{
    if (p->prev_term != 113)  zer113(p);
    return p->rf * p->tf[4][0];
}

static double zer115(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          19448.0;
    rf = rf * r2 +           -80080.0;
    rf = rf * r2 +           135135.0;
    rf = rf * r2 +          -120120.0;
    rf = rf * r2 +            60060.0;
    rf = rf * r2 +           -16632.0;
    rf = rf * r2 +             2310.0;
    rf = rf * r2 +             -120.0;
    p->rf = (double)rf;
    p->prev_term = 115;
    return p->rf * p->tf[3][1];
}

static double zer116(mparams *p)
{
    if (p->prev_term != 115)  zer115(p);
    return p->rf * p->tf[3][0];
}

static double zer117(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          43758.0;
    rf = rf * r2 +          -194480.0;
    rf = rf * r2 +           360360.0;
    rf = rf * r2 +          -360360.0;
    rf = rf * r2 +           210210.0;
    rf = rf * r2 +           -72072.0;
    rf = rf * r2 +            13860.0;
    rf = rf * r2 +            -1320.0;
    rf = rf * r2 +               45.0;
    p->rf = (double)rf;
    p->prev_term = 117;
    return p->rf * p->tf[2][1];
}

static double zer118(mparams *p)
{
    if (p->prev_term != 117)  zer117(p);
    return p->rf * p->tf[2][0];
}

static double zer119(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          92378.0;
    rf = rf * r2 +          -437580.0;
    rf = rf * r2 +           875160.0;
    rf = rf * r2 +          -960960.0;
    rf = rf * r2 +           630630.0;
    rf = rf * r2 +          -252252.0;
    rf = rf * r2 +            60060.0;
    rf = rf * r2 +            -7920.0;
    rf = rf * r2 +              495.0;
    rf = rf * r2 +              -10.0;
    p->rf = (double)rf;
    p->prev_term = 119;
    return p->rf * p->tf[1][1];
}

static double zer120(mparams *p)
{
    if (p->prev_term != 119)  zer119(p);
    return p->rf * p->tf[1][0];
}

static double zer121(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         184756.0;
    rf = rf * r2 +          -923780.0;
    rf = rf * r2 +          1969110.0;
    rf = rf * r2 +         -2333760.0;
    rf = rf * r2 +          1681680.0;
    rf = rf * r2 +          -756756.0;
    rf = rf * r2 +           210210.0;
    rf = rf * r2 +           -34320.0;
    rf = rf * r2 +             2970.0;
    rf = rf * r2 +             -110.0;
    rf = rf * r2 +                1.0;
    return rf;
}

static double zer122(mparams *p)
{
    return p->tf[11][1];
}

static double zer123(mparams *p)
{
    return p->tf[11][0];
}

static double zer124(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             12.0;
    rf = rf * r2 +              -11.0;
    p->rf = (double)rf;
    p->prev_term = 124;
    return p->rf * p->tf[10][1];
}

static double zer125(mparams *p)
{
    if (p->prev_term != 124)  zer124(p);
    return p->rf * p->tf[10][0];
}

static double zer126(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             78.0;
    rf = rf * r2 +             -132.0;
    rf = rf * r2 +               55.0;
    p->rf = (double)rf;
    p->prev_term = 126;
    return p->rf * p->tf[9][1];
}

static double zer127(mparams *p)
{
    if (p->prev_term != 126)  zer126(p);
    return p->rf * p->tf[9][0];
}

static double zer128(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            364.0;
    rf = rf * r2 +             -858.0;
    rf = rf * r2 +              660.0;
    rf = rf * r2 +             -165.0;
    p->rf = (double)rf;
    p->prev_term = 128;
    return p->rf * p->tf[8][1];
}

static double zer129(mparams *p)
{
    if (p->prev_term != 128)  zer128(p);
    return p->rf * p->tf[8][0];
}

static double zer130(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           1365.0;
    rf = rf * r2 +            -4004.0;
    rf = rf * r2 +             4290.0;
    rf = rf * r2 +            -1980.0;
    rf = rf * r2 +              330.0;
    p->rf = (double)rf;
    p->prev_term = 130;
    return p->rf * p->tf[7][1];
}

static double zer131(mparams *p)
{
    if (p->prev_term != 130)  zer130(p);
    return p->rf * p->tf[7][0];
}

static double zer132(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           4368.0;
    rf = rf * r2 +           -15015.0;
    rf = rf * r2 +            20020.0;
    rf = rf * r2 +           -12870.0;
    rf = rf * r2 +             3960.0;
    rf = rf * r2 +             -462.0;
    p->rf = (double)rf;
    p->prev_term = 132;
    return p->rf * p->tf[6][1];
}

static double zer133(mparams *p)
{
    if (p->prev_term != 132)  zer132(p);
    return p->rf * p->tf[6][0];
}

static double zer134(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          12376.0;
    rf = rf * r2 +           -48048.0;
    rf = rf * r2 +            75075.0;
    rf = rf * r2 +           -60060.0;
    rf = rf * r2 +            25740.0;
    rf = rf * r2 +            -5544.0;
    rf = rf * r2 +              462.0;
    p->rf = (double)rf;
    p->prev_term = 134;
    return p->rf * p->tf[5][1];
}

static double zer135(mparams *p)
{
    if (p->prev_term != 134)  zer134(p);
    return p->rf * p->tf[5][0];
}

static double zer136(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          31824.0;
    rf = rf * r2 +          -136136.0;
    rf = rf * r2 +           240240.0;
    rf = rf * r2 +          -225225.0;
    rf = rf * r2 +           120120.0;
    rf = rf * r2 +           -36036.0;
    rf = rf * r2 +             5544.0;
    rf = rf * r2 +             -330.0;
    p->rf = (double)rf;
    p->prev_term = 136;
    return p->rf * p->tf[4][1];
}

static double zer137(mparams *p)
{
    if (p->prev_term != 136)  zer136(p);
    return p->rf * p->tf[4][0];
}

static double zer138(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          75582.0;
    rf = rf * r2 +          -350064.0;
    rf = rf * r2 +           680680.0;
    rf = rf * r2 +          -720720.0;
    rf = rf * r2 +           450450.0;
    rf = rf * r2 +          -168168.0;
    rf = rf * r2 +            36036.0;
    rf = rf * r2 +            -3960.0;
    rf = rf * r2 +              165.0;
    p->rf = (double)rf;
    p->prev_term = 138;
    return p->rf * p->tf[3][1];
}

static double zer139(mparams *p)
{
    if (p->prev_term != 138)  zer138(p);
    return p->rf * p->tf[3][0];
}

static double zer140(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         167960.0;
    rf = rf * r2 +          -831402.0;
    rf = rf * r2 +          1750320.0;
    rf = rf * r2 +         -2042040.0;
    rf = rf * r2 +          1441440.0;
    rf = rf * r2 +          -630630.0;
    rf = rf * r2 +           168168.0;
    rf = rf * r2 +           -25740.0;
    rf = rf * r2 +             1980.0;
    rf = rf * r2 +              -55.0;
    p->rf = (double)rf;
    p->prev_term = 140;
    return p->rf * p->tf[2][1];
}

static double zer141(mparams *p)
{
    if (p->prev_term != 140)  zer140(p);
    return p->rf * p->tf[2][0];
}

static double zer142(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         352716.0;
    rf = rf * r2 +         -1847560.0;
    rf = rf * r2 +          4157010.0;
    rf = rf * r2 +         -5250960.0;
    rf = rf * r2 +          4084080.0;
    rf = rf * r2 +         -2018016.0;
    rf = rf * r2 +           630630.0;
    rf = rf * r2 +          -120120.0;
    rf = rf * r2 +            12870.0;
    rf = rf * r2 +             -660.0;
    rf = rf * r2 +               11.0;
    p->rf = (double)rf;
    p->prev_term = 142;
    return p->rf * p->tf[1][1];
}

static double zer143(mparams *p)
{
    if (p->prev_term != 142)  zer142(p);
    return p->rf * p->tf[1][0];
}

static double zer144(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         705432.0;
    rf = rf * r2 +         -3879876.0;
    rf = rf * r2 +          9237800.0;
    rf = rf * r2 +        -12471030.0;
    rf = rf * r2 +         10501920.0;
    rf = rf * r2 +         -5717712.0;
    rf = rf * r2 +          2018016.0;
    rf = rf * r2 +          -450450.0;
    rf = rf * r2 +            60060.0;
    rf = rf * r2 +            -4290.0;
    rf = rf * r2 +              132.0;
    rf = rf * r2 +               -1.0;
    return rf;
}

static double zer145(mparams *p)
{
    return p->tf[12][1];
}

static double zer146(mparams *p)
{
    return p->tf[12][0];
}

static double zer147(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             13.0;
    rf = rf * r2 +              -12.0;
    p->rf = (double)rf;
    p->prev_term = 147;
    return p->rf * p->tf[11][1];
}

static double zer148(mparams *p)
{
    if (p->prev_term != 147)  zer147(p);
    return p->rf * p->tf[11][0];
}

static double zer149(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             91.0;
    rf = rf * r2 +             -156.0;
    rf = rf * r2 +               66.0;
    p->rf = (double)rf;
    p->prev_term = 149;
    return p->rf * p->tf[10][1];
}

static double zer150(mparams *p)
{
    if (p->prev_term != 149)  zer149(p);
    return p->rf * p->tf[10][0];
}

static double zer151(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            455.0;
    rf = rf * r2 +            -1092.0;
    rf = rf * r2 +              858.0;
    rf = rf * r2 +             -220.0;
    p->rf = (double)rf;
    p->prev_term = 151;
    return p->rf * p->tf[9][1];
}

static double zer152(mparams *p)
{
    if (p->prev_term != 151)  zer151(p);
    return p->rf * p->tf[9][0];
}

static double zer153(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           1820.0;
    rf = rf * r2 +            -5460.0;
    rf = rf * r2 +             6006.0;
    rf = rf * r2 +            -2860.0;
    rf = rf * r2 +              495.0;
    p->rf = (double)rf;
    p->prev_term = 153;
    return p->rf * p->tf[8][1];
}

static double zer154(mparams *p)
{
    if (p->prev_term != 153)  zer153(p);
    return p->rf * p->tf[8][0];
}

static double zer155(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           6188.0;
    rf = rf * r2 +           -21840.0;
    rf = rf * r2 +            30030.0;
    rf = rf * r2 +           -20020.0;
    rf = rf * r2 +             6435.0;
    rf = rf * r2 +             -792.0;
    p->rf = (double)rf;
    p->prev_term = 155;
    return p->rf * p->tf[7][1];
}

static double zer156(mparams *p)
{
    if (p->prev_term != 155)  zer155(p);
    return p->rf * p->tf[7][0];
}

static double zer157(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          18564.0;
    rf = rf * r2 +           -74256.0;
    rf = rf * r2 +           120120.0;
    rf = rf * r2 +          -100100.0;
    rf = rf * r2 +            45045.0;
    rf = rf * r2 +           -10296.0;
    rf = rf * r2 +              924.0;
    p->rf = (double)rf;
    p->prev_term = 157;
    return p->rf * p->tf[6][1];
}

static double zer158(mparams *p)
{
    if (p->prev_term != 157)  zer157(p);
    return p->rf * p->tf[6][0];
}

static double zer159(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          50388.0;
    rf = rf * r2 +          -222768.0;
    rf = rf * r2 +           408408.0;
    rf = rf * r2 +          -400400.0;
    rf = rf * r2 +           225225.0;
    rf = rf * r2 +           -72072.0;
    rf = rf * r2 +            12012.0;
    rf = rf * r2 +             -792.0;
    p->rf = (double)rf;
    p->prev_term = 159;
    return p->rf * p->tf[5][1];
}

static double zer160(mparams *p)
{
    if (p->prev_term != 159)  zer159(p);
    return p->rf * p->tf[5][0];
}

static double zer161(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         125970.0;
    rf = rf * r2 +          -604656.0;
    rf = rf * r2 +          1225224.0;
    rf = rf * r2 +         -1361360.0;
    rf = rf * r2 +           900900.0;
    rf = rf * r2 +          -360360.0;
    rf = rf * r2 +            84084.0;
    rf = rf * r2 +           -10296.0;
    rf = rf * r2 +              495.0;
    p->rf = (double)rf;
    p->prev_term = 161;
    return p->rf * p->tf[4][1];
}

static double zer162(mparams *p)
{
    if (p->prev_term != 161)  zer161(p);
    return p->rf * p->tf[4][0];
}

static double zer163(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         293930.0;
    rf = rf * r2 +         -1511640.0;
    rf = rf * r2 +          3325608.0;
    rf = rf * r2 +         -4084080.0;
    rf = rf * r2 +          3063060.0;
    rf = rf * r2 +         -1441440.0;
    rf = rf * r2 +           420420.0;
    rf = rf * r2 +           -72072.0;
    rf = rf * r2 +             6435.0;
    rf = rf * r2 +             -220.0;
    p->rf = (double)rf;
    p->prev_term = 163;
    return p->rf * p->tf[3][1];
}

static double zer164(mparams *p)
{
    if (p->prev_term != 163)  zer163(p);
    return p->rf * p->tf[3][0];
}

static double zer165(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         646646.0;
    rf = rf * r2 +         -3527160.0;
    rf = rf * r2 +          8314020.0;
    rf = rf * r2 +        -11085360.0;
    rf = rf * r2 +          9189180.0;
    rf = rf * r2 +         -4900896.0;
    rf = rf * r2 +          1681680.0;
    rf = rf * r2 +          -360360.0;
    rf = rf * r2 +            45045.0;
    rf = rf * r2 +            -2860.0;
    rf = rf * r2 +               66.0;
    p->rf = (double)rf;
    p->prev_term = 165;
    return p->rf * p->tf[2][1];
}

static double zer166(mparams *p)
{
    if (p->prev_term != 165)  zer165(p);
    return p->rf * p->tf[2][0];
}

static double zer167(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        1352078.0;
    rf = rf * r2 +         -7759752.0;
    rf = rf * r2 +         19399380.0;
    rf = rf * r2 +        -27713400.0;
    rf = rf * r2 +         24942060.0;
    rf = rf * r2 +        -14702688.0;
    rf = rf * r2 +          5717712.0;
    rf = rf * r2 +         -1441440.0;
    rf = rf * r2 +           225225.0;
    rf = rf * r2 +           -20020.0;
    rf = rf * r2 +              858.0;
    rf = rf * r2 +              -12.0;
    p->rf = (double)rf;
    p->prev_term = 167;
    return p->rf * p->tf[1][1];
}

static double zer168(mparams *p)
{
    if (p->prev_term != 167)  zer167(p);
    return p->rf * p->tf[1][0];
}

static double zer169(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        2704156.0;
    rf = rf * r2 +        -16224936.0;
    rf = rf * r2 +         42678636.0;
    rf = rf * r2 +        -64664600.0;
    rf = rf * r2 +         62355150.0;
    rf = rf * r2 +        -39907296.0;
    rf = rf * r2 +         17153136.0;
    rf = rf * r2 +         -4900896.0;
    rf = rf * r2 +           900900.0;
    rf = rf * r2 +          -100100.0;
    rf = rf * r2 +             6006.0;
    rf = rf * r2 +             -156.0;
    rf = rf * r2 +                1.0;
    return rf;
}

static double zer170(mparams *p)
{
    return p->tf[13][1];
}

static double zer171(mparams *p)
{
    return p->tf[13][0];
}

static double zer172(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             14.0;
    rf = rf * r2 +              -13.0;
    p->rf = (double)rf;
    p->prev_term = 172;
    return p->rf * p->tf[12][1];
}

static double zer173(mparams *p)
{
    if (p->prev_term != 172)  zer172(p);
    return p->rf * p->tf[12][0];
}

static double zer174(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            105.0;
    rf = rf * r2 +             -182.0;
    rf = rf * r2 +               78.0;
    p->rf = (double)rf;
    p->prev_term = 174;
    return p->rf * p->tf[11][1];
}

static double zer175(mparams *p)
{
    if (p->prev_term != 174)  zer174(p);
    return p->rf * p->tf[11][0];
}

static double zer176(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            560.0;
    rf = rf * r2 +            -1365.0;
    rf = rf * r2 +             1092.0;
    rf = rf * r2 +             -286.0;
    p->rf = (double)rf;
    p->prev_term = 176;
    return p->rf * p->tf[10][1];
}

static double zer177(mparams *p)
{
    if (p->prev_term != 176)  zer176(p);
    return p->rf * p->tf[10][0];
}

static double zer178(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           2380.0;
    rf = rf * r2 +            -7280.0;
    rf = rf * r2 +             8190.0;
    rf = rf * r2 +            -4004.0;
    rf = rf * r2 +              715.0;
    p->rf = (double)rf;
    p->prev_term = 178;
    return p->rf * p->tf[9][1];
}

static double zer179(mparams *p)
{
    if (p->prev_term != 178)  zer178(p);
    return p->rf * p->tf[9][0];
}

static double zer180(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           8568.0;
    rf = rf * r2 +           -30940.0;
    rf = rf * r2 +            43680.0;
    rf = rf * r2 +           -30030.0;
    rf = rf * r2 +            10010.0;
    rf = rf * r2 +            -1287.0;
    p->rf = (double)rf;
    p->prev_term = 180;
    return p->rf * p->tf[8][1];
}

static double zer181(mparams *p)
{
    if (p->prev_term != 180)  zer180(p);
    return p->rf * p->tf[8][0];
}

static double zer182(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          27132.0;
    rf = rf * r2 +          -111384.0;
    rf = rf * r2 +           185640.0;
    rf = rf * r2 +          -160160.0;
    rf = rf * r2 +            75075.0;
    rf = rf * r2 +           -18018.0;
    rf = rf * r2 +             1716.0;
    p->rf = (double)rf;
    p->prev_term = 182;
    return p->rf * p->tf[7][1];
}

static double zer183(mparams *p)
{
    if (p->prev_term != 182)  zer182(p);
    return p->rf * p->tf[7][0];
}

static double zer184(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          77520.0;
    rf = rf * r2 +          -352716.0;
    rf = rf * r2 +           668304.0;
    rf = rf * r2 +          -680680.0;
    rf = rf * r2 +           400400.0;
    rf = rf * r2 +          -135135.0;
    rf = rf * r2 +            24024.0;
    rf = rf * r2 +            -1716.0;
    p->rf = (double)rf;
    p->prev_term = 184;
    return p->rf * p->tf[6][1];
}

static double zer185(mparams *p)
{
    if (p->prev_term != 184)  zer184(p);
    return p->rf * p->tf[6][0];
}

static double zer186(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         203490.0;
    rf = rf * r2 +         -1007760.0;
    rf = rf * r2 +          2116296.0;
    rf = rf * r2 +         -2450448.0;
    rf = rf * r2 +          1701700.0;
    rf = rf * r2 +          -720720.0;
    rf = rf * r2 +           180180.0;
    rf = rf * r2 +           -24024.0;
    rf = rf * r2 +             1287.0;
    p->rf = (double)rf;
    p->prev_term = 186;
    return p->rf * p->tf[5][1];
}

static double zer187(mparams *p)
{
    if (p->prev_term != 186)  zer186(p);
    return p->rf * p->tf[5][0];
}

static double zer188(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         497420.0;
    rf = rf * r2 +         -2645370.0;
    rf = rf * r2 +          6046560.0;
    rf = rf * r2 +         -7759752.0;
    rf = rf * r2 +          6126120.0;
    rf = rf * r2 +         -3063060.0;
    rf = rf * r2 +           960960.0;
    rf = rf * r2 +          -180180.0;
    rf = rf * r2 +            18018.0;
    rf = rf * r2 +             -715.0;
    p->rf = (double)rf;
    p->prev_term = 188;
    return p->rf * p->tf[4][1];
}

static double zer189(mparams *p)
{
    if (p->prev_term != 188)  zer188(p);
    return p->rf * p->tf[4][0];
}

static double zer190(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        1144066.0;
    rf = rf * r2 +         -6466460.0;
    rf = rf * r2 +         15872220.0;
    rf = rf * r2 +        -22170720.0;
    rf = rf * r2 +         19399380.0;
    rf = rf * r2 +        -11027016.0;
    rf = rf * r2 +          4084080.0;
    rf = rf * r2 +          -960960.0;
    rf = rf * r2 +           135135.0;
    rf = rf * r2 +           -10010.0;
    rf = rf * r2 +              286.0;
    p->rf = (double)rf;
    p->prev_term = 190;
    return p->rf * p->tf[3][1];
}

static double zer191(mparams *p)
{
    if (p->prev_term != 190)  zer190(p);
    return p->rf * p->tf[3][0];
}

static double zer192(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        2496144.0;
    rf = rf * r2 +        -14872858.0;
    rf = rf * r2 +         38798760.0;
    rf = rf * r2 +        -58198140.0;
    rf = rf * r2 +         55426800.0;
    rf = rf * r2 +        -34918884.0;
    rf = rf * r2 +         14702688.0;
    rf = rf * r2 +         -4084080.0;
    rf = rf * r2 +           720720.0;
    rf = rf * r2 +           -75075.0;
    rf = rf * r2 +             4004.0;
    rf = rf * r2 +              -78.0;
    p->rf = (double)rf;
    p->prev_term = 192;
    return p->rf * p->tf[2][1];
}

static double zer193(mparams *p)
{
    if (p->prev_term != 192)  zer192(p);
    return p->rf * p->tf[2][0];
}

static double zer194(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        5200300.0;
    rf = rf * r2 +        -32449872.0;
    rf = rf * r2 +         89237148.0;
    rf = rf * r2 +       -142262120.0;
    rf = rf * r2 +        145495350.0;
    rf = rf * r2 +        -99768240.0;
    rf = rf * r2 +         46558512.0;
    rf = rf * r2 +        -14702688.0;
    rf = rf * r2 +          3063060.0;
    rf = rf * r2 +          -400400.0;
    rf = rf * r2 +            30030.0;
    rf = rf * r2 +            -1092.0;
    rf = rf * r2 +               13.0;
    p->rf = (double)rf;
    p->prev_term = 194;
    return p->rf * p->tf[1][1];
}

static double zer195(mparams *p)
{
    if (p->prev_term != 194)  zer194(p);
    return p->rf * p->tf[1][0];
}

static double zer196(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       10400600.0;
    rf = rf * r2 +        -67603900.0;
    rf = rf * r2 +        194699232.0;
    rf = rf * r2 +       -327202876.0;
    rf = rf * r2 +        355655300.0;
    rf = rf * r2 +       -261891630.0;
    rf = rf * r2 +        133024320.0;
    rf = rf * r2 +        -46558512.0;
    rf = rf * r2 +         11027016.0;
    rf = rf * r2 +         -1701700.0;
    rf = rf * r2 +           160160.0;
    rf = rf * r2 +            -8190.0;
    rf = rf * r2 +              182.0;
    rf = rf * r2 +               -1.0;
    return rf;
}

static double zer197(mparams *p)
{
    return p->tf[14][1];
}

static double zer198(mparams *p)
{
    return p->tf[14][0];
}

static double zer199(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             15.0;
    rf = rf * r2 +              -14.0;
    p->rf = (double)rf;
    p->prev_term = 199;
    return p->rf * p->tf[13][1];
}

static double zer200(mparams *p)
{
    if (p->prev_term != 199)  zer199(p);
    return p->rf * p->tf[13][0];
}

static double zer201(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            120.0;
    rf = rf * r2 +             -210.0;
    rf = rf * r2 +               91.0;
    p->rf = (double)rf;
    p->prev_term = 201;
    return p->rf * p->tf[12][1];
}

static double zer202(mparams *p)
{
    if (p->prev_term != 201)  zer201(p);
    return p->rf * p->tf[12][0];
}

static double zer203(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            680.0;
    rf = rf * r2 +            -1680.0;
    rf = rf * r2 +             1365.0;
    rf = rf * r2 +             -364.0;
    p->rf = (double)rf;
    p->prev_term = 203;
    return p->rf * p->tf[11][1];
}

static double zer204(mparams *p)
{
    if (p->prev_term != 203)  zer203(p);
    return p->rf * p->tf[11][0];
}

static double zer205(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           3060.0;
    rf = rf * r2 +            -9520.0;
    rf = rf * r2 +            10920.0;
    rf = rf * r2 +            -5460.0;
    rf = rf * r2 +             1001.0;
    p->rf = (double)rf;
    p->prev_term = 205;
    return p->rf * p->tf[10][1];
}

static double zer206(mparams *p)
{
    if (p->prev_term != 205)  zer205(p);
    return p->rf * p->tf[10][0];
}

static double zer207(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          11628.0;
    rf = rf * r2 +           -42840.0;
    rf = rf * r2 +            61880.0;
    rf = rf * r2 +           -43680.0;
    rf = rf * r2 +            15015.0;
    rf = rf * r2 +            -2002.0;
    p->rf = (double)rf;
    p->prev_term = 207;
    return p->rf * p->tf[9][1];
}

static double zer208(mparams *p)
{
    if (p->prev_term != 207)  zer207(p);
    return p->rf * p->tf[9][0];
}

static double zer209(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          38760.0;
    rf = rf * r2 +          -162792.0;
    rf = rf * r2 +           278460.0;
    rf = rf * r2 +          -247520.0;
    rf = rf * r2 +           120120.0;
    rf = rf * r2 +           -30030.0;
    rf = rf * r2 +             3003.0;
    p->rf = (double)rf;
    p->prev_term = 209;
    return p->rf * p->tf[8][1];
}

static double zer210(mparams *p)
{
    if (p->prev_term != 209)  zer209(p);
    return p->rf * p->tf[8][0];
}

static double zer211(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         116280.0;
    rf = rf * r2 +          -542640.0;
    rf = rf * r2 +          1058148.0;
    rf = rf * r2 +         -1113840.0;
    rf = rf * r2 +           680680.0;
    rf = rf * r2 +          -240240.0;
    rf = rf * r2 +            45045.0;
    rf = rf * r2 +            -3432.0;
    p->rf = (double)rf;
    p->prev_term = 211;
    return p->rf * p->tf[7][1];
}

static double zer212(mparams *p)
{
    if (p->prev_term != 211)  zer211(p);
    return p->rf * p->tf[7][0];
}

static double zer213(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         319770.0;
    rf = rf * r2 +         -1627920.0;
    rf = rf * r2 +          3527160.0;
    rf = rf * r2 +         -4232592.0;
    rf = rf * r2 +          3063060.0;
    rf = rf * r2 +         -1361360.0;
    rf = rf * r2 +           360360.0;
    rf = rf * r2 +           -51480.0;
    rf = rf * r2 +             3003.0;
    p->rf = (double)rf;
    p->prev_term = 213;
    return p->rf * p->tf[6][1];
}

static double zer214(mparams *p)
{
    if (p->prev_term != 213)  zer213(p);
    return p->rf * p->tf[6][0];
}

static double zer215(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         817190.0;
    rf = rf * r2 +         -4476780.0;
    rf = rf * r2 +         10581480.0;
    rf = rf * r2 +        -14108640.0;
    rf = rf * r2 +         11639628.0;
    rf = rf * r2 +         -6126120.0;
    rf = rf * r2 +          2042040.0;
    rf = rf * r2 +          -411840.0;
    rf = rf * r2 +            45045.0;
    rf = rf * r2 +            -2002.0;
    p->rf = (double)rf;
    p->prev_term = 215;
    return p->rf * p->tf[5][1];
}

static double zer216(mparams *p)
{
    if (p->prev_term != 215)  zer215(p);
    return p->rf * p->tf[5][0];
}

static double zer217(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        1961256.0;
    rf = rf * r2 +        -11440660.0;
    rf = rf * r2 +         29099070.0;
    rf = rf * r2 +        -42325920.0;
    rf = rf * r2 +         38798760.0;
    rf = rf * r2 +        -23279256.0;
    rf = rf * r2 +          9189180.0;
    rf = rf * r2 +         -2333760.0;
    rf = rf * r2 +           360360.0;
    rf = rf * r2 +           -30030.0;
    rf = rf * r2 +             1001.0;
    p->rf = (double)rf;
    p->prev_term = 217;
    return p->rf * p->tf[4][1];
}

static double zer218(mparams *p)
{
    if (p->prev_term != 217)  zer217(p);
    return p->rf * p->tf[4][0];
}

static double zer219(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        4457400.0;
    rf = rf * r2 +        -27457584.0;
    rf = rf * r2 +         74364290.0;
    rf = rf * r2 +       -116396280.0;
    rf = rf * r2 +        116396280.0;
    rf = rf * r2 +        -77597520.0;
    rf = rf * r2 +         34918884.0;
    rf = rf * r2 +        -10501920.0;
    rf = rf * r2 +          2042040.0;
    rf = rf * r2 +          -240240.0;
    rf = rf * r2 +            15015.0;
    rf = rf * r2 +             -364.0;
    p->rf = (double)rf;
    p->prev_term = 219;
    return p->rf * p->tf[3][1];
}

static double zer220(mparams *p)
{
    if (p->prev_term != 219)  zer219(p);
    return p->rf * p->tf[3][0];
}

static double zer221(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        9657700.0;
    rf = rf * r2 +        -62403600.0;
    rf = rf * r2 +        178474296.0;
    rf = rf * r2 +       -297457160.0;
    rf = rf * r2 +        320089770.0;
    rf = rf * r2 +       -232792560.0;
    rf = rf * r2 +        116396280.0;
    rf = rf * r2 +        -39907296.0;
    rf = rf * r2 +          9189180.0;
    rf = rf * r2 +         -1361360.0;
    rf = rf * r2 +           120120.0;
    rf = rf * r2 +            -5460.0;
    rf = rf * r2 +               91.0;
    p->rf = (double)rf;
    p->prev_term = 221;
    return p->rf * p->tf[2][1];
}

static double zer222(mparams *p)
{
    if (p->prev_term != 221)  zer221(p);
    return p->rf * p->tf[2][0];
}

static double zer223(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       20058300.0;
    rf = rf * r2 +       -135207800.0;
    rf = rf * r2 +        405623400.0;
    rf = rf * r2 +       -713897184.0;
    rf = rf * r2 +        818007190.0;
    rf = rf * r2 +       -640179540.0;
    rf = rf * r2 +        349188840.0;
    rf = rf * r2 +       -133024320.0;
    rf = rf * r2 +         34918884.0;
    rf = rf * r2 +         -6126120.0;
    rf = rf * r2 +           680680.0;
    rf = rf * r2 +           -43680.0;
    rf = rf * r2 +             1365.0;
    rf = rf * r2 +              -14.0;
    p->rf = (double)rf;
    p->prev_term = 223;
    return p->rf * p->tf[1][1];
}

static double zer224(mparams *p)
{
    if (p->prev_term != 223)  zer223(p);
    return p->rf * p->tf[1][0];
}

static double zer225(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       40116600.0;
    rf = rf * r2 +       -280816200.0;
    rf = rf * r2 +        878850700.0;
    rf = rf * r2 +      -1622493600.0;
    rf = rf * r2 +       1963217256.0;
    rf = rf * r2 +      -1636014380.0;
    rf = rf * r2 +        960269310.0;
    rf = rf * r2 +       -399072960.0;
    rf = rf * r2 +        116396280.0;
    rf = rf * r2 +        -23279256.0;
    rf = rf * r2 +          3063060.0;
    rf = rf * r2 +          -247520.0;
    rf = rf * r2 +            10920.0;
    rf = rf * r2 +             -210.0;
    rf = rf * r2 +                1.0;
    return rf;
}

static double zer226(mparams *p)
{
    return p->tf[15][1];
}

static double zer227(mparams *p)
{
    return p->tf[15][0];
}

static double zer228(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             16.0;
    rf = rf * r2 +              -15.0;
    p->rf = (double)rf;
    p->prev_term = 228;
    return p->rf * p->tf[14][1];
}

static double zer229(mparams *p)
{
    if (p->prev_term != 228)  zer228(p);
    return p->rf * p->tf[14][0];
}

static double zer230(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            136.0;
    rf = rf * r2 +             -240.0;
    rf = rf * r2 +              105.0;
    p->rf = (double)rf;
    p->prev_term = 230;
    return p->rf * p->tf[13][1];
}

static double zer231(mparams *p)
{
    if (p->prev_term != 230)  zer230(p);
    return p->rf * p->tf[13][0];
}

static double zer232(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            816.0;
    rf = rf * r2 +            -2040.0;
    rf = rf * r2 +             1680.0;
    rf = rf * r2 +             -455.0;
    p->rf = (double)rf;
    p->prev_term = 232;
    return p->rf * p->tf[12][1];
}

static double zer233(mparams *p)
{
    if (p->prev_term != 232)  zer232(p);
    return p->rf * p->tf[12][0];
}

static double zer234(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           3876.0;
    rf = rf * r2 +           -12240.0;
    rf = rf * r2 +            14280.0;
    rf = rf * r2 +            -7280.0;
    rf = rf * r2 +             1365.0;
    p->rf = (double)rf;
    p->prev_term = 234;
    return p->rf * p->tf[11][1];
}

static double zer235(mparams *p)
{
    if (p->prev_term != 234)  zer234(p);
    return p->rf * p->tf[11][0];
}

static double zer236(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          15504.0;
    rf = rf * r2 +           -58140.0;
    rf = rf * r2 +            85680.0;
    rf = rf * r2 +           -61880.0;
    rf = rf * r2 +            21840.0;
    rf = rf * r2 +            -3003.0;
    p->rf = (double)rf;
    p->prev_term = 236;
    return p->rf * p->tf[10][1];
}

static double zer237(mparams *p)
{
    if (p->prev_term != 236)  zer236(p);
    return p->rf * p->tf[10][0];
}

static double zer238(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          54264.0;
    rf = rf * r2 +          -232560.0;
    rf = rf * r2 +           406980.0;
    rf = rf * r2 +          -371280.0;
    rf = rf * r2 +           185640.0;
    rf = rf * r2 +           -48048.0;
    rf = rf * r2 +             5005.0;
    p->rf = (double)rf;
    p->prev_term = 238;
    return p->rf * p->tf[9][1];
}

static double zer239(mparams *p)
{
    if (p->prev_term != 238)  zer238(p);
    return p->rf * p->tf[9][0];
}

static double zer240(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         170544.0;
    rf = rf * r2 +          -813960.0;
    rf = rf * r2 +          1627920.0;
    rf = rf * r2 +         -1763580.0;
    rf = rf * r2 +          1113840.0;
    rf = rf * r2 +          -408408.0;
    rf = rf * r2 +            80080.0;
    rf = rf * r2 +            -6435.0;
    p->rf = (double)rf;
    p->prev_term = 240;
    return p->rf * p->tf[8][1];
}

static double zer241(mparams *p)
{
    if (p->prev_term != 240)  zer240(p);
    return p->rf * p->tf[8][0];
}

static double zer242(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         490314.0;
    rf = rf * r2 +         -2558160.0;
    rf = rf * r2 +          5697720.0;
    rf = rf * r2 +         -7054320.0;
    rf = rf * r2 +          5290740.0;
    rf = rf * r2 +         -2450448.0;
    rf = rf * r2 +           680680.0;
    rf = rf * r2 +          -102960.0;
    rf = rf * r2 +             6435.0;
    p->rf = (double)rf;
    p->prev_term = 242;
    return p->rf * p->tf[7][1];
}

static double zer243(mparams *p)
{
    if (p->prev_term != 242)  zer242(p);
    return p->rf * p->tf[7][0];
}

static double zer244(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        1307504.0;
    rf = rf * r2 +         -7354710.0;
    rf = rf * r2 +         17907120.0;
    rf = rf * r2 +        -24690120.0;
    rf = rf * r2 +         21162960.0;
    rf = rf * r2 +        -11639628.0;
    rf = rf * r2 +          4084080.0;
    rf = rf * r2 +          -875160.0;
    rf = rf * r2 +           102960.0;
    rf = rf * r2 +            -5005.0;
    p->rf = (double)rf;
    p->prev_term = 244;
    return p->rf * p->tf[6][1];
}

static double zer245(mparams *p)
{
    if (p->prev_term != 244)  zer244(p);
    return p->rf * p->tf[6][0];
}

static double zer246(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        3268760.0;
    rf = rf * r2 +        -19612560.0;
    rf = rf * r2 +         51482970.0;
    rf = rf * r2 +        -77597520.0;
    rf = rf * r2 +         74070360.0;
    rf = rf * r2 +        -46558512.0;
    rf = rf * r2 +         19399380.0;
    rf = rf * r2 +         -5250960.0;
    rf = rf * r2 +           875160.0;
    rf = rf * r2 +           -80080.0;
    rf = rf * r2 +             3003.0;
    p->rf = (double)rf;
    p->prev_term = 246;
    return p->rf * p->tf[5][1];
}

static double zer247(mparams *p)
{
    if (p->prev_term != 246)  zer246(p);
    return p->rf * p->tf[5][0];
}

static double zer248(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        7726160.0;
    rf = rf * r2 +        -49031400.0;
    rf = rf * r2 +        137287920.0;
    rf = rf * r2 +       -223092870.0;
    rf = rf * r2 +        232792560.0;
    rf = rf * r2 +       -162954792.0;
    rf = rf * r2 +         77597520.0;
    rf = rf * r2 +        -24942060.0;
    rf = rf * r2 +          5250960.0;
    rf = rf * r2 +          -680680.0;
    rf = rf * r2 +            48048.0;
    rf = rf * r2 +            -1365.0;
    p->rf = (double)rf;
    p->prev_term = 248;
    return p->rf * p->tf[4][1];
}

static double zer249(mparams *p)
{
    if (p->prev_term != 248)  zer248(p);
    return p->rf * p->tf[4][0];
}

static double zer250(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       17383860.0;
    rf = rf * r2 +       -115892400.0;
    rf = rf * r2 +        343219800.0;
    rf = rf * r2 +       -594914320.0;
    rf = rf * r2 +        669278610.0;
    rf = rf * r2 +       -512143632.0;
    rf = rf * r2 +        271591320.0;
    rf = rf * r2 +        -99768240.0;
    rf = rf * r2 +         24942060.0;
    rf = rf * r2 +         -4084080.0;
    rf = rf * r2 +           408408.0;
    rf = rf * r2 +           -21840.0;
    rf = rf * r2 +              455.0;
    p->rf = (double)rf;
    p->prev_term = 250;
    return p->rf * p->tf[3][1];
}

static double zer251(mparams *p)
{
    if (p->prev_term != 250)  zer250(p);
    return p->rf * p->tf[3][0];
}

static double zer252(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       37442160.0;
    rf = rf * r2 +       -260757900.0;
    rf = rf * r2 +        811246800.0;
    rf = rf * r2 +      -1487285800.0;
    rf = rf * r2 +       1784742960.0;
    rf = rf * r2 +      -1472412942.0;
    rf = rf * r2 +        853572720.0;
    rf = rf * r2 +       -349188840.0;
    rf = rf * r2 +         99768240.0;
    rf = rf * r2 +        -19399380.0;
    rf = rf * r2 +          2450448.0;
    rf = rf * r2 +          -185640.0;
    rf = rf * r2 +             7280.0;
    rf = rf * r2 +             -105.0;
    p->rf = (double)rf;
    p->prev_term = 252;
    return p->rf * p->tf[2][1];
}

static double zer253(mparams *p)
{
    if (p->prev_term != 252)  zer252(p);
    return p->rf * p->tf[2][0];
}

static double zer254(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       77558760.0;
    rf = rf * r2 +       -561632400.0;
    rf = rf * r2 +       1825305300.0;
    rf = rf * r2 +      -3515402800.0;
    rf = rf * r2 +       4461857400.0;
    rf = rf * r2 +      -3926434512.0;
    rf = rf * r2 +       2454021570.0;
    rf = rf * r2 +      -1097450640.0;
    rf = rf * r2 +        349188840.0;
    rf = rf * r2 +        -77597520.0;
    rf = rf * r2 +         11639628.0;
    rf = rf * r2 +         -1113840.0;
    rf = rf * r2 +            61880.0;
    rf = rf * r2 +            -1680.0;
    rf = rf * r2 +               15.0;
    p->rf = (double)rf;
    p->prev_term = 254;
    return p->rf * p->tf[1][1];
}

static double zer255(mparams *p)
{
    if (p->prev_term != 254)  zer254(p);
    return p->rf * p->tf[1][0];
}

static double zer256(mparams *p)
{
    long double r2 = p->r2;
    long double rf =      155117520.0;
    rf = rf * r2 +      -1163381400.0;
    rf = rf * r2 +       3931426800.0;
    rf = rf * r2 +      -7909656300.0;
    rf = rf * r2 +      10546208400.0;
    rf = rf * r2 +      -9816086280.0;
    rf = rf * r2 +       6544057520.0;
    rf = rf * r2 +      -3155170590.0;
    rf = rf * r2 +       1097450640.0;
    rf = rf * r2 +       -271591320.0;
    rf = rf * r2 +         46558512.0;
    rf = rf * r2 +         -5290740.0;
    rf = rf * r2 +           371280.0;
    rf = rf * r2 +           -14280.0;
    rf = rf * r2 +              240.0;
    rf = rf * r2 +               -1.0;
    return rf;
}

static double zer257(mparams *p)
{
    return p->tf[16][1];
}

static double zer258(mparams *p)
{
    return p->tf[16][0];
}

static double zer259(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             17.0;
    rf = rf * r2 +              -16.0;
    p->rf = (double)rf;
    p->prev_term = 259;
    return p->rf * p->tf[15][1];
}

static double zer260(mparams *p)
{
    if (p->prev_term != 259)  zer259(p);
    return p->rf * p->tf[15][0];
}

static double zer261(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            153.0;
    rf = rf * r2 +             -272.0;
    rf = rf * r2 +              120.0;
    p->rf = (double)rf;
    p->prev_term = 261;
    return p->rf * p->tf[14][1];
}

static double zer262(mparams *p)
{
    if (p->prev_term != 261)  zer261(p);
    return p->rf * p->tf[14][0];
}

static double zer263(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            969.0;
    rf = rf * r2 +            -2448.0;
    rf = rf * r2 +             2040.0;
    rf = rf * r2 +             -560.0;
    p->rf = (double)rf;
    p->prev_term = 263;
    return p->rf * p->tf[13][1];
}

static double zer264(mparams *p)
{
    if (p->prev_term != 263)  zer263(p);
    return p->rf * p->tf[13][0];
}

static double zer265(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           4845.0;
    rf = rf * r2 +           -15504.0;
    rf = rf * r2 +            18360.0;
    rf = rf * r2 +            -9520.0;
    rf = rf * r2 +             1820.0;
    p->rf = (double)rf;
    p->prev_term = 265;
    return p->rf * p->tf[12][1];
}

static double zer266(mparams *p)
{
    if (p->prev_term != 265)  zer265(p);
    return p->rf * p->tf[12][0];
}

static double zer267(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          20349.0;
    rf = rf * r2 +           -77520.0;
    rf = rf * r2 +           116280.0;
    rf = rf * r2 +           -85680.0;
    rf = rf * r2 +            30940.0;
    rf = rf * r2 +            -4368.0;
    p->rf = (double)rf;
    p->prev_term = 267;
    return p->rf * p->tf[11][1];
}

static double zer268(mparams *p)
{
    if (p->prev_term != 267)  zer267(p);
    return p->rf * p->tf[11][0];
}

static double zer269(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          74613.0;
    rf = rf * r2 +          -325584.0;
    rf = rf * r2 +           581400.0;
    rf = rf * r2 +          -542640.0;
    rf = rf * r2 +           278460.0;
    rf = rf * r2 +           -74256.0;
    rf = rf * r2 +             8008.0;
    p->rf = (double)rf;
    p->prev_term = 269;
    return p->rf * p->tf[10][1];
}

static double zer270(mparams *p)
{
    if (p->prev_term != 269)  zer269(p);
    return p->rf * p->tf[10][0];
}

static double zer271(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         245157.0;
    rf = rf * r2 +         -1193808.0;
    rf = rf * r2 +          2441880.0;
    rf = rf * r2 +         -2713200.0;
    rf = rf * r2 +          1763580.0;
    rf = rf * r2 +          -668304.0;
    rf = rf * r2 +           136136.0;
    rf = rf * r2 +           -11440.0;
    p->rf = (double)rf;
    p->prev_term = 271;
    return p->rf * p->tf[9][1];
}

static double zer272(mparams *p)
{
    if (p->prev_term != 271)  zer271(p);
    return p->rf * p->tf[9][0];
}

static double zer273(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         735471.0;
    rf = rf * r2 +         -3922512.0;
    rf = rf * r2 +          8953560.0;
    rf = rf * r2 +        -11395440.0;
    rf = rf * r2 +          8817900.0;
    rf = rf * r2 +         -4232592.0;
    rf = rf * r2 +          1225224.0;
    rf = rf * r2 +          -194480.0;
    rf = rf * r2 +            12870.0;
    p->rf = (double)rf;
    p->prev_term = 273;
    return p->rf * p->tf[8][1];
}

static double zer274(mparams *p)
{
    if (p->prev_term != 273)  zer273(p);
    return p->rf * p->tf[8][0];
}

static double zer275(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        2042975.0;
    rf = rf * r2 +        -11767536.0;
    rf = rf * r2 +         29418840.0;
    rf = rf * r2 +        -41783280.0;
    rf = rf * r2 +         37035180.0;
    rf = rf * r2 +        -21162960.0;
    rf = rf * r2 +          7759752.0;
    rf = rf * r2 +         -1750320.0;
    rf = rf * r2 +           218790.0;
    rf = rf * r2 +           -11440.0;
    p->rf = (double)rf;
    p->prev_term = 275;
    return p->rf * p->tf[7][1];
}

static double zer276(mparams *p)
{
    if (p->prev_term != 275)  zer275(p);
    return p->rf * p->tf[7][0];
}

static double zer277(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        5311735.0;
    rf = rf * r2 +        -32687600.0;
    rf = rf * r2 +         88256520.0;
    rf = rf * r2 +       -137287920.0;
    rf = rf * r2 +        135795660.0;
    rf = rf * r2 +        -88884432.0;
    rf = rf * r2 +         38798760.0;
    rf = rf * r2 +        -11085360.0;
    rf = rf * r2 +          1969110.0;
    rf = rf * r2 +          -194480.0;
    rf = rf * r2 +             8008.0;
    p->rf = (double)rf;
    p->prev_term = 277;
    return p->rf * p->tf[6][1];
}

static double zer278(mparams *p)
{
    if (p->prev_term != 277)  zer277(p);
    return p->rf * p->tf[6][0];
}

static double zer279(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       13037895.0;
    rf = rf * r2 +        -84987760.0;
    rf = rf * r2 +        245157000.0;
    rf = rf * r2 +       -411863760.0;
    rf = rf * r2 +        446185740.0;
    rf = rf * r2 +       -325909584.0;
    rf = rf * r2 +        162954792.0;
    rf = rf * r2 +        -55426800.0;
    rf = rf * r2 +         12471030.0;
    rf = rf * r2 +         -1750320.0;
    rf = rf * r2 +           136136.0;
    rf = rf * r2 +            -4368.0;
    p->rf = (double)rf;
    p->prev_term = 279;
    return p->rf * p->tf[5][1];
}

static double zer280(mparams *p)
{
    if (p->prev_term != 279)  zer279(p);
    return p->rf * p->tf[5][0];
}

static double zer281(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       30421755.0;
    rf = rf * r2 +       -208606320.0;
    rf = rf * r2 +        637408200.0;
    rf = rf * r2 +      -1144066000.0;
    rf = rf * r2 +       1338557220.0;
    rf = rf * r2 +      -1070845776.0;
    rf = rf * r2 +        597500904.0;
    rf = rf * r2 +       -232792560.0;
    rf = rf * r2 +         62355150.0;
    rf = rf * r2 +        -11085360.0;
    rf = rf * r2 +          1225224.0;
    rf = rf * r2 +           -74256.0;
    rf = rf * r2 +             1820.0;
    p->rf = (double)rf;
    p->prev_term = 281;
    return p->rf * p->tf[4][1];
}

static double zer282(mparams *p)
{
    if (p->prev_term != 281)  zer281(p);
    return p->rf * p->tf[4][0];
}

static double zer283(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       67863915.0;
    rf = rf * r2 +       -486748080.0;
    rf = rf * r2 +       1564547400.0;
    rf = rf * r2 +      -2974571600.0;
    rf = rf * r2 +       3718214500.0;
    rf = rf * r2 +      -3212537328.0;
    rf = rf * r2 +       1963217256.0;
    rf = rf * r2 +       -853572720.0;
    rf = rf * r2 +        261891630.0;
    rf = rf * r2 +        -55426800.0;
    rf = rf * r2 +          7759752.0;
    rf = rf * r2 +          -668304.0;
    rf = rf * r2 +            30940.0;
    rf = rf * r2 +             -560.0;
    p->rf = (double)rf;
    p->prev_term = 283;
    return p->rf * p->tf[3][1];
}

static double zer284(mparams *p)
{
    if (p->prev_term != 283)  zer283(p);
    return p->rf * p->tf[3][0];
}

static double zer285(mparams *p)
{
    long double r2 = p->r2;
    long double rf =      145422675.0;
    rf = rf * r2 +      -1085822640.0;
    rf = rf * r2 +       3650610600.0;
    rf = rf * r2 +      -7301221200.0;
    rf = rf * r2 +       9667357700.0;
    rf = rf * r2 +      -8923714800.0;
    rf = rf * r2 +       5889651768.0;
    rf = rf * r2 +      -2804596080.0;
    rf = rf * r2 +        960269310.0;
    rf = rf * r2 +       -232792560.0;
    rf = rf * r2 +         38798760.0;
    rf = rf * r2 +         -4232592.0;
    rf = rf * r2 +           278460.0;
    rf = rf * r2 +            -9520.0;
    rf = rf * r2 +              120.0;
    p->rf = (double)rf;
    p->prev_term = 285;
    return p->rf * p->tf[2][1];
}

static double zer286(mparams *p)
{
    if (p->prev_term != 285)  zer285(p);
    return p->rf * p->tf[2][0];
}

static double zer287(mparams *p)
{
    long double r2 = p->r2;
    long double rf =      300540195.0;
    rf = rf * r2 +      -2326762800.0;
    rf = rf * r2 +       8143669800.0;
    rf = rf * r2 +     -17036182800.0;
    rf = rf * r2 +      23728968900.0;
    rf = rf * r2 +     -23201658480.0;
    rf = rf * r2 +      16360143800.0;
    rf = rf * r2 +      -8413788240.0;
    rf = rf * r2 +       3155170590.0;
    rf = rf * r2 +       -853572720.0;
    rf = rf * r2 +        162954792.0;
    rf = rf * r2 +        -21162960.0;
    rf = rf * r2 +          1763580.0;
    rf = rf * r2 +           -85680.0;
    rf = rf * r2 +             2040.0;
    rf = rf * r2 +              -16.0;
    p->rf = (double)rf;
    p->prev_term = 287;
    return p->rf * p->tf[1][1];
}

static double zer288(mparams *p)
{
    if (p->prev_term != 287)  zer287(p);
    return p->rf * p->tf[1][0];
}

static double zer289(mparams *p)
{
    long double r2 = p->r2;
    long double rf =      601080390.0;
    rf = rf * r2 +      -4808643120.0;
    rf = rf * r2 +      17450721000.0;
    rf = rf * r2 +     -38003792400.0;
    rf = rf * r2 +      55367594100.0;
    rf = rf * r2 +     -56949525360.0;
    rf = rf * r2 +      42536373880.0;
    rf = rf * r2 +     -23371634000.0;
    rf = rf * r2 +       9465511770.0;
    rf = rf * r2 +      -2804596080.0;
    rf = rf * r2 +        597500904.0;
    rf = rf * r2 +        -88884432.0;
    rf = rf * r2 +          8817900.0;
    rf = rf * r2 +          -542640.0;
    rf = rf * r2 +            18360.0;
    rf = rf * r2 +             -272.0;
    rf = rf * r2 +                1.0;
    return rf;
}

static double zer290(mparams *p)
{
    return p->tf[17][1];
}

static double zer291(mparams *p)
{
    return p->tf[17][0];
}

static double zer292(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             18.0;
    rf = rf * r2 +              -17.0;
    p->rf = (double)rf;
    p->prev_term = 292;
    return p->rf * p->tf[16][1];
}

static double zer293(mparams *p)
{
    if (p->prev_term != 292)  zer292(p);
    return p->rf * p->tf[16][0];
}

static double zer294(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            171.0;
    rf = rf * r2 +             -306.0;
    rf = rf * r2 +              136.0;
    p->rf = (double)rf;
    p->prev_term = 294;
    return p->rf * p->tf[15][1];
}

static double zer295(mparams *p)
{
    if (p->prev_term != 294)  zer294(p);
    return p->rf * p->tf[15][0];
}

static double zer296(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           1140.0;
    rf = rf * r2 +            -2907.0;
    rf = rf * r2 +             2448.0;
    rf = rf * r2 +             -680.0;
    p->rf = (double)rf;
    p->prev_term = 296;
    return p->rf * p->tf[14][1];
}

static double zer297(mparams *p)
{
    if (p->prev_term != 296)  zer296(p);
    return p->rf * p->tf[14][0];
}

static double zer298(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           5985.0;
    rf = rf * r2 +           -19380.0;
    rf = rf * r2 +            23256.0;
    rf = rf * r2 +           -12240.0;
    rf = rf * r2 +             2380.0;
    p->rf = (double)rf;
    p->prev_term = 298;
    return p->rf * p->tf[13][1];
}

static double zer299(mparams *p)
{
    if (p->prev_term != 298)  zer298(p);
    return p->rf * p->tf[13][0];
}

static double zer300(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          26334.0;
    rf = rf * r2 +          -101745.0;
    rf = rf * r2 +           155040.0;
    rf = rf * r2 +          -116280.0;
    rf = rf * r2 +            42840.0;
    rf = rf * r2 +            -6188.0;
    p->rf = (double)rf;
    p->prev_term = 300;
    return p->rf * p->tf[12][1];
}

static double zer301(mparams *p)
{
    if (p->prev_term != 300)  zer300(p);
    return p->rf * p->tf[12][0];
}

static double zer302(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         100947.0;
    rf = rf * r2 +          -447678.0;
    rf = rf * r2 +           813960.0;
    rf = rf * r2 +          -775200.0;
    rf = rf * r2 +           406980.0;
    rf = rf * r2 +          -111384.0;
    rf = rf * r2 +            12376.0;
    p->rf = (double)rf;
    p->prev_term = 302;
    return p->rf * p->tf[11][1];
}

static double zer303(mparams *p)
{
    if (p->prev_term != 302)  zer302(p);
    return p->rf * p->tf[11][0];
}

static double zer304(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         346104.0;
    rf = rf * r2 +         -1716099.0;
    rf = rf * r2 +          3581424.0;
    rf = rf * r2 +         -4069800.0;
    rf = rf * r2 +          2713200.0;
    rf = rf * r2 +         -1058148.0;
    rf = rf * r2 +           222768.0;
    rf = rf * r2 +           -19448.0;
    p->rf = (double)rf;
    p->prev_term = 304;
    return p->rf * p->tf[10][1];
}

static double zer305(mparams *p)
{
    if (p->prev_term != 304)  zer304(p);
    return p->rf * p->tf[10][0];
}

static double zer306(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        1081575.0;
    rf = rf * r2 +         -5883768.0;
    rf = rf * r2 +         13728792.0;
    rf = rf * r2 +        -17907120.0;
    rf = rf * r2 +         14244300.0;
    rf = rf * r2 +         -7054320.0;
    rf = rf * r2 +          2116296.0;
    rf = rf * r2 +          -350064.0;
    rf = rf * r2 +            24310.0;
    p->rf = (double)rf;
    p->prev_term = 306;
    return p->rf * p->tf[9][1];
}

static double zer307(mparams *p)
{
    if (p->prev_term != 306)  zer306(p);
    return p->rf * p->tf[9][0];
}

static double zer308(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        3124550.0;
    rf = rf * r2 +        -18386775.0;
    rf = rf * r2 +         47070144.0;
    rf = rf * r2 +        -68643960.0;
    rf = rf * r2 +         62674920.0;
    rf = rf * r2 +        -37035180.0;
    rf = rf * r2 +         14108640.0;
    rf = rf * r2 +         -3325608.0;
    rf = rf * r2 +           437580.0;
    rf = rf * r2 +           -24310.0;
    p->rf = (double)rf;
    p->prev_term = 308;
    return p->rf * p->tf[8][1];
}

static double zer309(mparams *p)
{
    if (p->prev_term != 308)  zer308(p);
    return p->rf * p->tf[8][0];
}

static double zer310(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        8436285.0;
    rf = rf * r2 +        -53117350.0;
    rf = rf * r2 +        147094200.0;
    rf = rf * r2 +       -235350720.0;
    rf = rf * r2 +        240253860.0;
    rf = rf * r2 +       -162954792.0;
    rf = rf * r2 +         74070360.0;
    rf = rf * r2 +        -22170720.0;
    rf = rf * r2 +          4157010.0;
    rf = rf * r2 +          -437580.0;
    rf = rf * r2 +            19448.0;
    p->rf = (double)rf;
    p->prev_term = 310;
    return p->rf * p->tf[7][1];
}

static double zer311(mparams *p)
{
    if (p->prev_term != 310)  zer310(p);
    return p->rf * p->tf[7][0];
}

static double zer312(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       21474180.0;
    rf = rf * r2 +       -143416845.0;
    rf = rf * r2 +        424938800.0;
    rf = rf * r2 +       -735471000.0;
    rf = rf * r2 +        823727520.0;
    rf = rf * r2 +       -624660036.0;
    rf = rf * r2 +        325909584.0;
    rf = rf * r2 +       -116396280.0;
    rf = rf * r2 +         27713400.0;
    rf = rf * r2 +         -4157010.0;
    rf = rf * r2 +           350064.0;
    rf = rf * r2 +           -12376.0;
    p->rf = (double)rf;
    p->prev_term = 312;
    return p->rf * p->tf[6][1];
}

static double zer313(mparams *p)
{
    if (p->prev_term != 312)  zer312(p);
    return p->rf * p->tf[6][0];
}

static double zer314(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       51895935.0;
    rf = rf * r2 +       -365061060.0;
    rf = rf * r2 +       1147334760.0;
    rf = rf * r2 +      -2124694000.0;
    rf = rf * r2 +       2574148500.0;
    rf = rf * r2 +      -2141691552.0;
    rf = rf * r2 +       1249320072.0;
    rf = rf * r2 +       -512143632.0;
    rf = rf * r2 +        145495350.0;
    rf = rf * r2 +        -27713400.0;
    rf = rf * r2 +          3325608.0;
    rf = rf * r2 +          -222768.0;
    rf = rf * r2 +             6188.0;
    p->rf = (double)rf;
    p->prev_term = 314;
    return p->rf * p->tf[5][1];
}

static double zer315(mparams *p)
{
    if (p->prev_term != 314)  zer314(p);
    return p->rf * p->tf[5][0];
}

static double zer316(mparams *p)
{
    long double r2 = p->r2;
    long double rf =      119759850.0;
    rf = rf * r2 +       -882230895.0;
    rf = rf * r2 +       2920488480.0;
    rf = rf * r2 +      -5736673800.0;
    rf = rf * r2 +       7436429000.0;
    rf = rf * r2 +      -6692786100.0;
    rf = rf * r2 +       4283383104.0;
    rf = rf * r2 +      -1963217256.0;
    rf = rf * r2 +        640179540.0;
    rf = rf * r2 +       -145495350.0;
    rf = rf * r2 +         22170720.0;
    rf = rf * r2 +         -2116296.0;
    rf = rf * r2 +           111384.0;
    rf = rf * r2 +            -2380.0;
    p->rf = (double)rf;
    p->prev_term = 316;
    return p->rf * p->tf[4][1];
}

static double zer317(mparams *p)
{
    if (p->prev_term != 316)  zer316(p);
    return p->rf * p->tf[4][0];
}

static double zer318(mparams *p)
{
    long double r2 = p->r2;
    long double rf =      265182525.0;
    rf = rf * r2 +      -2035917450.0;
    rf = rf * r2 +       7057847160.0;
    rf = rf * r2 +     -14602442400.0;
    rf = rf * r2 +      20078358300.0;
    rf = rf * r2 +     -19334715400.0;
    rf = rf * r2 +      13385572200.0;
    rf = rf * r2 +      -6731030592.0;
    rf = rf * r2 +       2454021570.0;
    rf = rf * r2 +       -640179540.0;
    rf = rf * r2 +        116396280.0;
    rf = rf * r2 +        -14108640.0;
    rf = rf * r2 +          1058148.0;
    rf = rf * r2 +           -42840.0;
    rf = rf * r2 +              680.0;
    p->rf = (double)rf;
    p->prev_term = 318;
    return p->rf * p->tf[3][1];
}

static double zer319(mparams *p)
{
    if (p->prev_term != 318)  zer318(p);
    return p->rf * p->tf[3][0];
}

static double zer320(mparams *p)
{
    long double r2 = p->r2;
    long double rf =      565722720.0;
    rf = rf * r2 +      -4508102925.0;
    rf = rf * r2 +      16287339600.0;
    rf = rf * r2 +     -35289235800.0;
    rf = rf * r2 +      51108548400.0;
    rf = rf * r2 +     -52203731580.0;
    rf = rf * r2 +      38669430800.0;
    rf = rf * r2 +     -21034470600.0;
    rf = rf * r2 +       8413788240.0;
    rf = rf * r2 +      -2454021570.0;
    rf = rf * r2 +        512143632.0;
    rf = rf * r2 +        -74070360.0;
    rf = rf * r2 +          7054320.0;
    rf = rf * r2 +          -406980.0;
    rf = rf * r2 +            12240.0;
    rf = rf * r2 +             -136.0;
    p->rf = (double)rf;
    p->prev_term = 320;
    return p->rf * p->tf[2][1];
}

static double zer321(mparams *p)
{
    if (p->prev_term != 320)  zer320(p);
    return p->rf * p->tf[2][0];
}

static double zer322(mparams *p)
{
    long double r2 = p->r2;
    long double rf =     1166803110.0;
    rf = rf * r2 +      -9617286240.0;
    rf = rf * r2 +      36064823400.0;
    rf = rf * r2 +     -81436698000.0;
    rf = rf * r2 +     123512325300.0;
    rf = rf * r2 +    -132882225840.0;
    rf = rf * r2 +     104407463160.0;
    rf = rf * r2 +     -60766248400.0;
    rf = rf * r2 +      26293088250.0;
    rf = rf * r2 +      -8413788240.0;
    rf = rf * r2 +       1963217256.0;
    rf = rf * r2 +       -325909584.0;
    rf = rf * r2 +         37035180.0;
    rf = rf * r2 +         -2713200.0;
    rf = rf * r2 +           116280.0;
    rf = rf * r2 +            -2448.0;
    rf = rf * r2 +               17.0;
    p->rf = (double)rf;
    p->prev_term = 322;
    return p->rf * p->tf[1][1];
}

static double zer323(mparams *p)
{
    if (p->prev_term != 322)  zer322(p);
    return p->rf * p->tf[1][0];
}

static double zer324(mparams *p)
{
    long double r2 = p->r2;
    long double rf =     2333606220.0;
    rf = rf * r2 +     -19835652870.0;
    rf = rf * r2 +      76938289920.0;
    rf = rf * r2 +    -180324117000.0;
    rf = rf * r2 +     285028443000.0;
    rf = rf * r2 +    -321132045780.0;
    rf = rf * r2 +     265764451680.0;
    rf = rf * r2 +    -164068870680.0;
    rf = rf * r2 +      75957810500.0;
    rf = rf * r2 +     -26293088250.0;
    rf = rf * r2 +       6731030592.0;
    rf = rf * r2 +      -1249320072.0;
    rf = rf * r2 +        162954792.0;
    rf = rf * r2 +        -14244300.0;
    rf = rf * r2 +           775200.0;
    rf = rf * r2 +           -23256.0;
    rf = rf * r2 +              306.0;
    rf = rf * r2 +               -1.0;
    return rf;
}

static double zer325(mparams *p)
{
    return p->tf[18][1];
}

static double zer326(mparams *p)
{
    return p->tf[18][0];
}

static double zer327(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             19.0;
    rf = rf * r2 +              -18.0;
    p->rf = (double)rf;
    p->prev_term = 327;
    return p->rf * p->tf[17][1];
}

static double zer328(mparams *p)
{
    if (p->prev_term != 327)  zer327(p);
    return p->rf * p->tf[17][0];
}

static double zer329(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            190.0;
    rf = rf * r2 +             -342.0;
    rf = rf * r2 +              153.0;
    p->rf = (double)rf;
    p->prev_term = 329;
    return p->rf * p->tf[16][1];
}

static double zer330(mparams *p)
{
    if (p->prev_term != 329)  zer329(p);
    return p->rf * p->tf[16][0];
}

static double zer331(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           1330.0;
    rf = rf * r2 +            -3420.0;
    rf = rf * r2 +             2907.0;
    rf = rf * r2 +             -816.0;
    p->rf = (double)rf;
    p->prev_term = 331;
    return p->rf * p->tf[15][1];
}

static double zer332(mparams *p)
{
    if (p->prev_term != 331)  zer331(p);
    return p->rf * p->tf[15][0];
}

static double zer333(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           7315.0;
    rf = rf * r2 +           -23940.0;
    rf = rf * r2 +            29070.0;
    rf = rf * r2 +           -15504.0;
    rf = rf * r2 +             3060.0;
    p->rf = (double)rf;
    p->prev_term = 333;
    return p->rf * p->tf[14][1];
}

static double zer334(mparams *p)
{
    if (p->prev_term != 333)  zer333(p);
    return p->rf * p->tf[14][0];
}

static double zer335(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          33649.0;
    rf = rf * r2 +          -131670.0;
    rf = rf * r2 +           203490.0;
    rf = rf * r2 +          -155040.0;
    rf = rf * r2 +            58140.0;
    rf = rf * r2 +            -8568.0;
    p->rf = (double)rf;
    p->prev_term = 335;
    return p->rf * p->tf[13][1];
}

static double zer336(mparams *p)
{
    if (p->prev_term != 335)  zer335(p);
    return p->rf * p->tf[13][0];
}

static double zer337(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         134596.0;
    rf = rf * r2 +          -605682.0;
    rf = rf * r2 +          1119195.0;
    rf = rf * r2 +         -1085280.0;
    rf = rf * r2 +           581400.0;
    rf = rf * r2 +          -162792.0;
    rf = rf * r2 +            18564.0;
    p->rf = (double)rf;
    p->prev_term = 337;
    return p->rf * p->tf[12][1];
}

static double zer338(mparams *p)
{
    if (p->prev_term != 337)  zer337(p);
    return p->rf * p->tf[12][0];
}

static double zer339(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         480700.0;
    rf = rf * r2 +         -2422728.0;
    rf = rf * r2 +          5148297.0;
    rf = rf * r2 +         -5969040.0;
    rf = rf * r2 +          4069800.0;
    rf = rf * r2 +         -1627920.0;
    rf = rf * r2 +           352716.0;
    rf = rf * r2 +           -31824.0;
    p->rf = (double)rf;
    p->prev_term = 339;
    return p->rf * p->tf[11][1];
}

static double zer340(mparams *p)
{
    if (p->prev_term != 339)  zer339(p);
    return p->rf * p->tf[11][0];
}

static double zer341(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        1562275.0;
    rf = rf * r2 +         -8652600.0;
    rf = rf * r2 +         20593188.0;
    rf = rf * r2 +        -27457584.0;
    rf = rf * r2 +         22383900.0;
    rf = rf * r2 +        -11395440.0;
    rf = rf * r2 +          3527160.0;
    rf = rf * r2 +          -604656.0;
    rf = rf * r2 +            43758.0;
    p->rf = (double)rf;
    p->prev_term = 341;
    return p->rf * p->tf[10][1];
}

static double zer342(mparams *p)
{
    if (p->prev_term != 341)  zer341(p);
    return p->rf * p->tf[10][0];
}

static double zer343(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        4686825.0;
    rf = rf * r2 +        -28120950.0;
    rf = rf * r2 +         73547100.0;
    rf = rf * r2 +       -109830336.0;
    rf = rf * r2 +        102965940.0;
    rf = rf * r2 +        -62674920.0;
    rf = rf * r2 +         24690120.0;
    rf = rf * r2 +         -6046560.0;
    rf = rf * r2 +           831402.0;
    rf = rf * r2 +           -48620.0;
    p->rf = (double)rf;
    p->prev_term = 343;
    return p->rf * p->tf[9][1];
}

static double zer344(mparams *p)
{
    if (p->prev_term != 343)  zer343(p);
    return p->rf * p->tf[9][0];
}

static double zer345(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       13123110.0;
    rf = rf * r2 +        -84362850.0;
    rf = rf * r2 +        239028075.0;
    rf = rf * r2 +       -392251200.0;
    rf = rf * r2 +        411863760.0;
    rf = rf * r2 +       -288304632.0;
    rf = rf * r2 +        135795660.0;
    rf = rf * r2 +        -42325920.0;
    rf = rf * r2 +          8314020.0;
    rf = rf * r2 +          -923780.0;
    rf = rf * r2 +            43758.0;
    p->rf = (double)rf;
    p->prev_term = 345;
    return p->rf * p->tf[8][1];
}

static double zer346(mparams *p)
{
    if (p->prev_term != 345)  zer345(p);
    return p->rf * p->tf[8][0];
}

static double zer347(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       34597290.0;
    rf = rf * r2 +       -236215980.0;
    rf = rf * r2 +        717084225.0;
    rf = rf * r2 +      -1274816400.0;
    rf = rf * r2 +       1470942000.0;
    rf = rf * r2 +      -1153218528.0;
    rf = rf * r2 +        624660036.0;
    rf = rf * r2 +       -232792560.0;
    rf = rf * r2 +         58198140.0;
    rf = rf * r2 +         -9237800.0;
    rf = rf * r2 +           831402.0;
    rf = rf * r2 +           -31824.0;
    p->rf = (double)rf;
    p->prev_term = 347;
    return p->rf * p->tf[7][1];
}

static double zer348(mparams *p)
{
    if (p->prev_term != 347)  zer347(p);
    return p->rf * p->tf[7][0];
}

static double zer349(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       86493225.0;
    rf = rf * r2 +       -622751220.0;
    rf = rf * r2 +       2007835830.0;
    rf = rf * r2 +      -3824449200.0;
    rf = rf * r2 +       4780561500.0;
    rf = rf * r2 +      -4118637600.0;
    rf = rf * r2 +       2498640144.0;
    rf = rf * r2 +      -1070845776.0;
    rf = rf * r2 +        320089770.0;
    rf = rf * r2 +        -64664600.0;
    rf = rf * r2 +          8314020.0;
    rf = rf * r2 +          -604656.0;
    rf = rf * r2 +            18564.0;
    p->rf = (double)rf;
    p->prev_term = 349;
    return p->rf * p->tf[6][1];
}

static double zer350(mparams *p)
{
    if (p->prev_term != 349)  zer349(p);
    return p->rf * p->tf[6][0];
}

static double zer351(mparams *p)
{
    long double r2 = p->r2;
    long double rf =      206253075.0;
    rf = rf * r2 +      -1556878050.0;
    rf = rf * r2 +       5293385370.0;
    rf = rf * r2 +     -10708457760.0;
    rf = rf * r2 +      14341684500.0;
    rf = rf * r2 +     -13385572200.0;
    rf = rf * r2 +       8923714800.0;
    rf = rf * r2 +      -4283383104.0;
    rf = rf * r2 +       1472412942.0;
    rf = rf * r2 +       -355655300.0;
    rf = rf * r2 +         58198140.0;
    rf = rf * r2 +         -6046560.0;
    rf = rf * r2 +           352716.0;
    rf = rf * r2 +            -8568.0;
    p->rf = (double)rf;
    p->prev_term = 351;
    return p->rf * p->tf[5][1];
}

static double zer352(mparams *p)
{
    if (p->prev_term != 351)  zer351(p);
    return p->rf * p->tf[5][0];
}

static double zer353(mparams *p)
{
    long double r2 = p->r2;
    long double rf =      471435600.0;
    rf = rf * r2 +      -3712555350.0;
    rf = rf * r2 +      13233463425.0;
    rf = rf * r2 +     -28231388640.0;
    rf = rf * r2 +      40156716600.0;
    rf = rf * r2 +     -40156716600.0;
    rf = rf * r2 +      29002073100.0;
    rf = rf * r2 +     -15297796800.0;
    rf = rf * r2 +       5889651768.0;
    rf = rf * r2 +      -1636014380.0;
    rf = rf * r2 +        320089770.0;
    rf = rf * r2 +        -42325920.0;
    rf = rf * r2 +          3527160.0;
    rf = rf * r2 +          -162792.0;
    rf = rf * r2 +             3060.0;
    p->rf = (double)rf;
    p->prev_term = 353;
    return p->rf * p->tf[4][1];
}

static double zer354(mparams *p)
{
    if (p->prev_term != 353)  zer353(p);
    return p->rf * p->tf[4][0];
}

static double zer355(mparams *p)
{
    long double r2 = p->r2;
    long double rf =     1037158320.0;
    rf = rf * r2 +      -8485840800.0;
    rf = rf * r2 +      31556720475.0;
    rf = rf * r2 +     -70578471600.0;
    rf = rf * r2 +     105867707400.0;
    rf = rf * r2 +    -112438806480.0;
    rf = rf * r2 +      87006219300.0;
    rf = rf * r2 +     -49717839600.0;
    rf = rf * r2 +      21034470600.0;
    rf = rf * r2 +      -6544057520.0;
    rf = rf * r2 +       1472412942.0;
    rf = rf * r2 +       -232792560.0;
    rf = rf * r2 +         24690120.0;
    rf = rf * r2 +         -1627920.0;
    rf = rf * r2 +            58140.0;
    rf = rf * r2 +             -816.0;
    p->rf = (double)rf;
    p->prev_term = 355;
    return p->rf * p->tf[3][1];
}

static double zer356(mparams *p)
{
    if (p->prev_term != 355)  zer355(p);
    return p->rf * p->tf[3][0];
}

static double zer357(mparams *p)
{
    long double r2 = p->r2;
    long double rf =     2203961430.0;
    rf = rf * r2 +     -18668849760.0;
    rf = rf * r2 +      72129646800.0;
    rf = rf * r2 +    -168302509200.0;
    rf = rf * r2 +     264669268500.0;
    rf = rf * r2 +    -296429580720.0;
    rf = rf * r2 +     243617414040.0;
    rf = rf * r2 +    -149153518800.0;
    rf = rf * r2 +      68362029450.0;
    rf = rf * r2 +     -23371634000.0;
    rf = rf * r2 +       5889651768.0;
    rf = rf * r2 +      -1070845776.0;
    rf = rf * r2 +        135795660.0;
    rf = rf * r2 +        -11395440.0;
    rf = rf * r2 +           581400.0;
    rf = rf * r2 +           -15504.0;
    rf = rf * r2 +              153.0;
    p->rf = (double)rf;
    p->prev_term = 357;
    return p->rf * p->tf[2][1];
}

static double zer358(mparams *p)
{
    if (p->prev_term != 357)  zer357(p);
    return p->rf * p->tf[2][0];
}

static double zer359(mparams *p)
{
    long double r2 = p->r2;
    long double rf =     4537567650.0;
    rf = rf * r2 +     -39671305740.0;
    rf = rf * r2 +     158685222960.0;
    rf = rf * r2 +    -384691449600.0;
    rf = rf * r2 +     631134409500.0;
    rf = rf * r2 +    -741073951800.0;
    rf = rf * r2 +     642264091560.0;
    rf = rf * r2 +    -417629852640.0;
    rf = rf * r2 +     205086088350.0;
    rf = rf * r2 +     -75957810500.0;
    rf = rf * r2 +      21034470600.0;
    rf = rf * r2 +      -4283383104.0;
    rf = rf * r2 +        624660036.0;
    rf = rf * r2 +        -62674920.0;
    rf = rf * r2 +          4069800.0;
    rf = rf * r2 +          -155040.0;
    rf = rf * r2 +             2907.0;
    rf = rf * r2 +              -18.0;
    p->rf = (double)rf;
    p->prev_term = 359;
    return p->rf * p->tf[1][1];
}

static double zer360(mparams *p)
{
    if (p->prev_term != 359)  zer359(p);
    return p->rf * p->tf[1][0];
}

static double zer361(mparams *p)
{
    long double r2 = p->r2;
    long double rf =     9075135300.0;
    rf = rf * r2 +     -81676217700.0;
    rf = rf * r2 +     337206098790.0;
    rf = rf * r2 +    -846321189120.0;
    rf = rf * r2 +    1442592936000.0;
    rf = rf * r2 +   -1767176346600.0;
    rf = rf * r2 +    1605660228900.0;
    rf = rf * r2 +   -1101024156960.0;
    rf = rf * r2 +     574241047380.0;
    rf = rf * r2 +    -227873431500.0;
    rf = rf * r2 +      68362029450.0;
    rf = rf * r2 +     -15297796800.0;
    rf = rf * r2 +       2498640144.0;
    rf = rf * r2 +       -288304632.0;
    rf = rf * r2 +         22383900.0;
    rf = rf * r2 +         -1085280.0;
    rf = rf * r2 +            29070.0;
    rf = rf * r2 +             -342.0;
    rf = rf * r2 +                1.0;
    return rf;
}

static double zer362(mparams *p)
{
    return p->tf[19][1];
}

static double zer363(mparams *p)
{
    return p->tf[19][0];
}

static double zer364(mparams *p)
{
    long double r2 = p->r2;
    long double rf =             20.0;
    rf = rf * r2 +              -19.0;
    p->rf = (double)rf;
    p->prev_term = 364;
    return p->rf * p->tf[18][1];
}

static double zer365(mparams *p)
{
    if (p->prev_term != 364)  zer364(p);
    return p->rf * p->tf[18][0];
}

static double zer366(mparams *p)
{
    long double r2 = p->r2;
    long double rf =            210.0;
    rf = rf * r2 +             -380.0;
    rf = rf * r2 +              171.0;
    p->rf = (double)rf;
    p->prev_term = 366;
    return p->rf * p->tf[17][1];
}

static double zer367(mparams *p)
{
    if (p->prev_term != 366)  zer366(p);
    return p->rf * p->tf[17][0];
}

static double zer368(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           1540.0;
    rf = rf * r2 +            -3990.0;
    rf = rf * r2 +             3420.0;
    rf = rf * r2 +             -969.0;
    p->rf = (double)rf;
    p->prev_term = 368;
    return p->rf * p->tf[16][1];
}

static double zer369(mparams *p)
{
    if (p->prev_term != 368)  zer368(p);
    return p->rf * p->tf[16][0];
}

static double zer370(mparams *p)
{
    long double r2 = p->r2;
    long double rf =           8855.0;
    rf = rf * r2 +           -29260.0;
    rf = rf * r2 +            35910.0;
    rf = rf * r2 +           -19380.0;
    rf = rf * r2 +             3876.0;
    p->rf = (double)rf;
    p->prev_term = 370;
    return p->rf * p->tf[15][1];
}

static double zer371(mparams *p)
{
    if (p->prev_term != 370)  zer370(p);
    return p->rf * p->tf[15][0];
}

static double zer372(mparams *p)
{
    long double r2 = p->r2;
    long double rf =          42504.0;
    rf = rf * r2 +          -168245.0;
    rf = rf * r2 +           263340.0;
    rf = rf * r2 +          -203490.0;
    rf = rf * r2 +            77520.0;
    rf = rf * r2 +           -11628.0;
    p->rf = (double)rf;
    p->prev_term = 372;
    return p->rf * p->tf[14][1];
}

static double zer373(mparams *p)
{
    if (p->prev_term != 372)  zer372(p);
    return p->rf * p->tf[14][0];
}

static double zer374(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         177100.0;
    rf = rf * r2 +          -807576.0;
    rf = rf * r2 +          1514205.0;
    rf = rf * r2 +         -1492260.0;
    rf = rf * r2 +           813960.0;
    rf = rf * r2 +          -232560.0;
    rf = rf * r2 +            27132.0;
    p->rf = (double)rf;
    p->prev_term = 374;
    return p->rf * p->tf[13][1];
}

static double zer375(mparams *p)
{
    if (p->prev_term != 374)  zer374(p);
    return p->rf * p->tf[13][0];
}

static double zer376(mparams *p)
{
    long double r2 = p->r2;
    long double rf =         657800.0;
    rf = rf * r2 +         -3364900.0;
    rf = rf * r2 +          7268184.0;
    rf = rf * r2 +         -8580495.0;
    rf = rf * r2 +          5969040.0;
    rf = rf * r2 +         -2441880.0;
    rf = rf * r2 +           542640.0;
    rf = rf * r2 +           -50388.0;
    p->rf = (double)rf;
    p->prev_term = 376;
    return p->rf * p->tf[12][1];
}

static double zer377(mparams *p)
{
    if (p->prev_term != 376)  zer376(p);
    return p->rf * p->tf[12][0];
}

static double zer378(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        2220075.0;
    rf = rf * r2 +        -12498200.0;
    rf = rf * r2 +         30284100.0;
    rf = rf * r2 +        -41186376.0;
    rf = rf * r2 +         34321980.0;
    rf = rf * r2 +        -17907120.0;
    rf = rf * r2 +          5697720.0;
    rf = rf * r2 +         -1007760.0;
    rf = rf * r2 +            75582.0;
    p->rf = (double)rf;
    p->prev_term = 378;
    return p->rf * p->tf[11][1];
}

static double zer379(mparams *p)
{
    if (p->prev_term != 378)  zer378(p);
    return p->rf * p->tf[11][0];
}

static double zer380(mparams *p)
{
    long double r2 = p->r2;
    long double rf =        6906900.0;
    rf = rf * r2 +        -42181425.0;
    rf = rf * r2 +        112483800.0;
    rf = rf * r2 +       -171609900.0;
    rf = rf * r2 +        164745504.0;
    rf = rf * r2 +       -102965940.0;
    rf = rf * r2 +         41783280.0;
    rf = rf * r2 +        -10581480.0;
    rf = rf * r2 +          1511640.0;
    rf = rf * r2 +           -92378.0;
    p->rf = (double)rf;
    p->prev_term = 380;
    return p->rf * p->tf[10][1];
}

static double zer381(mparams *p)
{
    if (p->prev_term != 380)  zer380(p);
    return p->rf * p->tf[10][0];
}

static double zer382(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       20030010.0;
    rf = rf * r2 +       -131231100.0;
    rf = rf * r2 +        379632825.0;
    rf = rf * r2 +       -637408200.0;
    rf = rf * r2 +        686439600.0;
    rf = rf * r2 +       -494236512.0;
    rf = rf * r2 +        240253860.0;
    rf = rf * r2 +        -77597520.0;
    rf = rf * r2 +         15872220.0;
    rf = rf * r2 +         -1847560.0;
    rf = rf * r2 +            92378.0;
    p->rf = (double)rf;
    p->prev_term = 382;
    return p->rf * p->tf[9][1];
}

static double zer383(mparams *p)
{
    if (p->prev_term != 382)  zer382(p);
    return p->rf * p->tf[9][0];
}

static double zer384(mparams *p)
{
    long double r2 = p->r2;
    long double rf =       54627300.0;
    rf = rf * r2 +       -380570190.0;
    rf = rf * r2 +       1181079900.0;
    rf = rf * r2 +      -2151252675.0;
    rf = rf * r2 +       2549632800.0;
    rf = rf * r2 +      -2059318800.0;
    rf = rf * r2 +       1153218528.0;
    rf = rf * r2 +       -446185740.0;
    rf = rf * r2 +        116396280.0;
    rf = rf * r2 +        -19399380.0;
    rf = rf * r2 +          1847560.0;
    rf = rf * r2 +           -75582.0;
    p->rf = (double)rf;
    p->prev_term = 384;
    return p->rf * p->tf[8][1];
}

static double zer385(mparams *p)
{
    if (p->prev_term != 384)  zer384(p);
    return p->rf * p->tf[8][0];
}

static double zer386(mparams *p)
{
    long double r2 = p->r2;
    long double rf =      141120525.0;
    rf = rf * r2 +      -1037918700.0;
    rf = rf * r2 +       3425131710.0;
    rf = rf * r2 +      -6692786100.0;
    rf = rf * r2 +       8605010700.0;
    rf = rf * r2 +      -7648898400.0;
    rf = rf * r2 +       4805077200.0;
    rf = rf * r2 +      -2141691552.0;
    rf = rf * r2 +        669278610.0;
    rf = rf * r2 +       -142262120.0;
    rf = rf * r2 +         19399380.0;
    rf = rf * r2 +         -1511640.0;
    rf = rf * r2 +            50388.0;
    p->rf = (double)rf;
    p->prev_term = 386;
    return p->rf * p->tf[7][1];
}

static double zer387(mparams *p)
{
    if (p->prev_term != 386)  zer386(p);
    return p->rf * p->tf[7][0];
}

static double zer388(mparams *p)
{
    long double r2 = p->r2;
    long double rf =      347373600.0;
    rf = rf * r2 +      -2681289975.0;
    rf = rf * r2 +       9341268300.0;
    rf = rf * r2 +     -19409079690.0;
    rf = rf * r2 +      26771144400.0;
    rf = rf * r2 +     -25815032100.0;
    rf = rf * r2 +      17847429600.0;
    rf = rf * r2 +      -8923714800.0;
    rf = rf * r2 +       3212537328.0;
    rf = rf * r2 +       -818007190.0;
    rf = rf * r2 +        142262120.0;
    rf = rf * r2 +        -15872220.0;
    rf = rf * r2 +          1007760.0;
    rf = rf * r2 +           -27132.0;
    p->rf = (double)rf;
    p->prev_term = 388;
    return p->rf * p->tf[6][1];
}

static double zer389(mparams *p)
{
    if (p->prev_term != 388)  zer388(p);
    return p->rf * p->tf[6][0];
}

static double zer390(mparams *p)
{
    long double r2 = p->r2;
    long double rf =      818809200.0;
    rf = rf * r2 +      -6600098400.0;
    rf = rf * r2 +      24131609775.0;
    rf = rf * r2 +     -52933853700.0;
    rf = rf * r2 +      77636318760.0;
    rf = rf * r2 +     -80313433200.0;
    rf = rf * r2 +      60235074900.0;
    rf = rf * r2 +     -33145226400.0;
    rf = rf * r2 +      13385572200.0;
    rf = rf * r2 +      -3926434512.0;
    rf = rf * r2 +        818007190.0;
    rf = rf * r2 +       -116396280.0;
    rf = rf * r2 +         10581480.0;
    rf = rf * r2 +          -542640.0;
    rf = rf * r2 +            11628.0;
    p->rf = (double)rf;
    p->prev_term = 390;
    return p->rf * p->tf[5][1];
}

static double zer391(mparams *p)
{
    if (p->prev_term != 390)  zer390(p);
    return p->rf * p->tf[5][0];
}

static double zer392(mparams *p)
{
    long double r2 = p->r2;
    long double rf =     1855967520.0;
    rf = rf * r2 +     -15557374800.0;
    rf = rf * r2 +      59400885600.0;
    rf = rf * r2 +    -136745788725.0;
    rf = rf * r2 +     211735414800.0;
    rf = rf * r2 +    -232908956280.0;
    rf = rf * r2 +     187398010800.0;
    rf = rf * r2 +    -111865139100.0;
    rf = rf * r2 +      49717839600.0;
    rf = rf * r2 +     -16360143800.0;
    rf = rf * r2 +       3926434512.0;
    rf = rf * r2 +       -669278610.0;
    rf = rf * r2 +         77597520.0;
    rf = rf * r2 +         -5697720.0;
    rf = rf * r2 +           232560.0;
    rf = rf * r2 +            -3876.0;
    p->rf = (double)rf;
    p->prev_term = 392;
    return p->rf * p->tf[4][1];
}

static double zer393(mparams *p)
{
    if (p->prev_term != 392)  zer392(p);
    return p->rf * p->tf[4][0];
}

static double zer394(mparams *p)
{
    long double r2 = p->r2;
    long double rf =     4059928950.0;
    rf = rf * r2 +     -35263382880.0;
    rf = rf * r2 +     140016373200.0;
    rf = rf * r2 +    -336605018400.0;
    rf = rf * r2 +     546983154900.0;
    rf = rf * r2 +    -635206244400.0;
    rf = rf * r2 +     543454231320.0;
    rf = rf * r2 +    -348024877200.0;
    rf = rf * r2 +     167797708650.0;
    rf = rf * r2 +     -60766248400.0;
    rf = rf * r2 +      16360143800.0;
    rf = rf * r2 +      -3212537328.0;
    rf = rf * r2 +        446185740.0;
    rf = rf * r2 +        -41783280.0;
    rf = rf * r2 +          2441880.0;
    rf = rf * r2 +           -77520.0;
    rf = rf * r2 +              969.0;
    p->rf = (double)rf;
    p->prev_term = 394;
    return p->rf * p->tf[3][1];
}

static double zer395(mparams *p)
{
    if (p->prev_term != 394)  zer394(p);
    return p->rf * p->tf[3][0];
}

static double zer396(mparams *p)
{
    long double r2 = p->r2;
    long double rf =     8597496600.0;
    rf = rf * r2 +     -77138650050.0;
    rf = rf * r2 +     317370445920.0;
    rf = rf * r2 +    -793426114800.0;
    rf = rf * r2 +    1346420073600.0;
    rf = rf * r2 +   -1640949464700.0;
    rf = rf * r2 +    1482147903600.0;
    rf = rf * r2 +   -1009272143880.0;
    rf = rf * r2 +     522037315800.0;
    rf = rf * r2 +    -205086088350.0;
    rf = rf * r2 +      60766248400.0;
    rf = rf * r2 +     -13385572200.0;
    rf = rf * r2 +       2141691552.0;
    rf = rf * r2 +       -240253860.0;
    rf = rf * r2 +         17907120.0;
    rf = rf * r2 +          -813960.0;
    rf = rf * r2 +            19380.0;
    rf = rf * r2 +             -171.0;
    p->rf = (double)rf;
    p->prev_term = 396;
    return p->rf * p->tf[2][1];
}

static double zer397(mparams *p)
{
    if (p->prev_term != 396)  zer396(p);
    return p->rf * p->tf[2][0];
}

static double zer398(mparams *p)
{
    long double r2 = p->r2;
    long double rf =    17672631900.0;
    rf = rf * r2 +    -163352435400.0;
    rf = rf * r2 +     694247850450.0;
    rf = rf * r2 +   -1798432526880.0;
    rf = rf * r2 +    3173704459200.0;
    rf = rf * r2 +   -4039260220800.0;
    rf = rf * r2 +    3828882084300.0;
    rf = rf * r2 +   -2752560392400.0;
    rf = rf * r2 +    1513908215820.0;
    rf = rf * r2 +    -638045608200.0;
    rf = rf * r2 +     205086088350.0;
    rf = rf * r2 +     -49717839600.0;
    rf = rf * r2 +       8923714800.0;
    rf = rf * r2 +      -1153218528.0;
    rf = rf * r2 +        102965940.0;
    rf = rf * r2 +         -5969040.0;
    rf = rf * r2 +           203490.0;
    rf = rf * r2 +            -3420.0;
    rf = rf * r2 +               19.0;
    p->rf = (double)rf;
    p->prev_term = 398;
    return p->rf * p->tf[1][1];
}

static double zer399(mparams *p)
{
    if (p->prev_term != 398)  zer398(p);
    return p->rf * p->tf[1][0];
}

static double zer400(mparams *p)
{
    long double r2 = p->r2;
    long double rf =    35345263800.0;
    rf = rf * r2 +    -335780006100.0;
    rf = rf * r2 +    1470171918600.0;
    rf = rf * r2 +   -3934071152550.0;
    rf = rf * r2 +    7193730107520.0;
    rf = rf * r2 +   -9521113377600.0;
    rf = rf * r2 +    9424940515200.0;
    rf = rf * r2 +   -7110781013700.0;
    rf = rf * r2 +    4128840588600.0;
    rf = rf * r2 +   -1850332263780.0;
    rf = rf * r2 +     638045608200.0;
    rf = rf * r2 +    -167797708650.0;
    rf = rf * r2 +      33145226400.0;
    rf = rf * r2 +      -4805077200.0;
    rf = rf * r2 +        494236512.0;
    rf = rf * r2 +        -34321980.0;
    rf = rf * r2 +          1492260.0;
    rf = rf * r2 +           -35910.0;
    rf = rf * r2 +              380.0;
    rf = rf * r2 +               -1.0;
    return rf;
}

static double (*zfa[401])(mparams *p) = {
    zer000, zer001, zer002, zer003, zer004, zer005, zer006, zer007, zer008, zer009,
    zer010, zer011, zer012, zer013, zer014, zer015, zer016, zer017, zer018, zer019,
    zer020, zer021, zer022, zer023, zer024, zer025, zer026, zer027, zer028, zer029,
    zer030, zer031, zer032, zer033, zer034, zer035, zer036, zer037, zer038, zer039,
    zer040, zer041, zer042, zer043, zer044, zer045, zer046, zer047, zer048, zer049,
    zer050, zer051, zer052, zer053, zer054, zer055, zer056, zer057, zer058, zer059,
    zer060, zer061, zer062, zer063, zer064, zer065, zer066, zer067, zer068, zer069,
    zer070, zer071, zer072, zer073, zer074, zer075, zer076, zer077, zer078, zer079,
    zer080, zer081, zer082, zer083, zer084, zer085, zer086, zer087, zer088, zer089,
    zer090, zer091, zer092, zer093, zer094, zer095, zer096, zer097, zer098, zer099,
    zer100, zer101, zer102, zer103, zer104, zer105, zer106, zer107, zer108, zer109,
    zer110, zer111, zer112, zer113, zer114, zer115, zer116, zer117, zer118, zer119,
    zer120, zer121, zer122, zer123, zer124, zer125, zer126, zer127, zer128, zer129,
    zer130, zer131, zer132, zer133, zer134, zer135, zer136, zer137, zer138, zer139,
    zer140, zer141, zer142, zer143, zer144, zer145, zer146, zer147, zer148, zer149,
    zer150, zer151, zer152, zer153, zer154, zer155, zer156, zer157, zer158, zer159,
    zer160, zer161, zer162, zer163, zer164, zer165, zer166, zer167, zer168, zer169,
    zer170, zer171, zer172, zer173, zer174, zer175, zer176, zer177, zer178, zer179,
    zer180, zer181, zer182, zer183, zer184, zer185, zer186, zer187, zer188, zer189,
    zer190, zer191, zer192, zer193, zer194, zer195, zer196, zer197, zer198, zer199,
    zer200, zer201, zer202, zer203, zer204, zer205, zer206, zer207, zer208, zer209,
    zer210, zer211, zer212, zer213, zer214, zer215, zer216, zer217, zer218, zer219,
    zer220, zer221, zer222, zer223, zer224, zer225, zer226, zer227, zer228, zer229,
    zer230, zer231, zer232, zer233, zer234, zer235, zer236, zer237, zer238, zer239,
    zer240, zer241, zer242, zer243, zer244, zer245, zer246, zer247, zer248, zer249,
    zer250, zer251, zer252, zer253, zer254, zer255, zer256, zer257, zer258, zer259,
    zer260, zer261, zer262, zer263, zer264, zer265, zer266, zer267, zer268, zer269,
    zer270, zer271, zer272, zer273, zer274, zer275, zer276, zer277, zer278, zer279,
    zer280, zer281, zer282, zer283, zer284, zer285, zer286, zer287, zer288, zer289,
    zer290, zer291, zer292, zer293, zer294, zer295, zer296, zer297, zer298, zer299,
    zer300, zer301, zer302, zer303, zer304, zer305, zer306, zer307, zer308, zer309,
    zer310, zer311, zer312, zer313, zer314, zer315, zer316, zer317, zer318, zer319,
    zer320, zer321, zer322, zer323, zer324, zer325, zer326, zer327, zer328, zer329,
    zer330, zer331, zer332, zer333, zer334, zer335, zer336, zer337, zer338, zer339,
    zer340, zer341, zer342, zer343, zer344, zer345, zer346, zer347, zer348, zer349,
    zer350, zer351, zer352, zer353, zer354, zer355, zer356, zer357, zer358, zer359,
    zer360, zer361, zer362, zer363, zer364, zer365, zer366, zer367, zer368, zer369,
    zer370, zer371, zer372, zer373, zer374, zer375, zer376, zer377, zer378, zer379,
    zer380, zer381, zer382, zer383, zer384, zer385, zer386, zer387, zer388, zer389,
    zer390, zer391, zer392, zer393, zer394, zer395, zer396, zer397, zer398, zer399,
    zer400,};

static int ncaltheta[] = {
      0,   0,   1,   1,   1,   2,   2,   2,   2,   2,
      3,   3,   3,   3,   3,   3,   3,   4,   4,   4,
      4,   4,   4,   4,   4,   4,   5,   5,   5,   5,
      5,   5,   5,   5,   5,   5,   5,   6,   6,   6,
      6,   6,   6,   6,   6,   6,   6,   6,   6,   6,
      7,   7,   7,   7,   7,   7,   7,   7,   7,   7,
      7,   7,   7,   7,   7,   8,   8,   8,   8,   8,
      8,   8,   8,   8,   8,   8,   8,   8,   8,   8,
      8,   8,   9,   9,   9,   9,   9,   9,   9,   9,
      9,   9,   9,   9,   9,   9,   9,   9,   9,   9,
      9,  10,  10,  10,  10,  10,  10,  10,  10,  10,
     10,  10,  10,  10,  10,  10,  10,  10,  10,  10,
     10,  10,  11,  11,  11,  11,  11,  11,  11,  11,
     11,  11,  11,  11,  11,  11,  11,  11,  11,  11,
     11,  11,  11,  11,  11,  12,  12,  12,  12,  12,
     12,  12,  12,  12,  12,  12,  12,  12,  12,  12,
     12,  12,  12,  12,  12,  12,  12,  12,  12,  12,
     13,  13,  13,  13,  13,  13,  13,  13,  13,  13,
     13,  13,  13,  13,  13,  13,  13,  13,  13,  13,
     13,  13,  13,  13,  13,  13,  13,  14,  14,  14,
     14,  14,  14,  14,  14,  14,  14,  14,  14,  14,
     14,  14,  14,  14,  14,  14,  14,  14,  14,  14,
     14,  14,  14,  14,  14,  14,  15,  15,  15,  15,
     15,  15,  15,  15,  15,  15,  15,  15,  15,  15,
     15,  15,  15,  15,  15,  15,  15,  15,  15,  15,
     15,  15,  15,  15,  15,  15,  15,  16,  16,  16,
     16,  16,  16,  16,  16,  16,  16,  16,  16,  16,
     16,  16,  16,  16,  16,  16,  16,  16,  16,  16,
     16,  16,  16,  16,  16,  16,  16,  16,  16,  16,
     17,  17,  17,  17,  17,  17,  17,  17,  17,  17,
     17,  17,  17,  17,  17,  17,  17,  17,  17,  17,
     17,  17,  17,  17,  17,  17,  17,  17,  17,  17,
     17,  17,  17,  17,  17,  18,  18,  18,  18,  18,
     18,  18,  18,  18,  18,  18,  18,  18,  18,  18,
     18,  18,  18,  18,  18,  18,  18,  18,  18,  18,
     18,  18,  18,  18,  18,  18,  18,  18,  18,  18,
     18,  18,  19,  19,  19,  19,  19,  19,  19,  19,
     19,  19,  19,  19,  19,  19,  19,  19,  19,  19,
     19,  19,  19,  19,  19,  19,  19,  19,  19,  19,
     19,  19,  19,  19,  19,  19,  19,  19,  19,  19,
     19,};

/* 第istt項から第iend項までのZernike関数値をzfarr[istt..iend]に詰めて返す.
 * [入力]
 *     istt : この項から計算開始
 *     iend : この項まで計算する
 *     x : X座標(0 ≦ x ≦ 1)
 *     y : Y座標(0 ≦ y ≦ 1)
 * [出力]
 *     zfarr[0..iend] : zfarr[0]はダミー. zfarr[1]はZ1, zfarr[2]はZ2...
 *                      zfarr[istt]からzfarr[iend]に計算結果が格納される.
 *                      少なくともzfarr[iend]までは配列が確保されている必要がある.
 *     戻り値 : 成功すればTRUE, 失敗したらFALSE
 */
Bool UpcZernike_getZf(int istt, int iend, double x, double y, double *zfarr)
{
    int i;
    mparams p;
    
    if (iend < istt) {
        i = istt;
        istt = iend;
        iend = i;
    }
    if (istt < 1 || iend > 400) {
        return FALSE;
    }
    p.prev_term = 0;
    p.r2 = x * x + y * y;
    cal_tf(x, y, ncaltheta[iend], &p);

    for (i = istt; i <= iend; i++) {
        zfarr[i] = zfa[i](&p);
    }
    return TRUE;
}

double UpcZernike_zval(int istt, int iend, double x, double y, const double *zcoef)
{
    mparams p;
    double v = 0.0;
    int i;
    
    if (iend < istt) {
        i = istt;
        istt = iend;
        iend = i;
    }
    if (istt < 1) {
        istt = 1;
    }
    if (iend > 400) {
        iend = 400;
    }

    p.prev_term = 0;
    p.r2 = x * x + y * y;
    cal_tf(x, y, ncaltheta[iend], &p);
    
    for (i = istt; i <= iend; i++) {
        v += zcoef[i] * zfa[i](&p);
    }
    return v;
}



