#include "GaussianBasisDatabase.hpp"

static const int NSEG_1 = 2;
static const int LVEC_1[2] = { 0,0 };
static const int NVEC_1[2] = { 2,1 };
static const double CMAT_1[3] = { 0.15631167,0.90466909,1.0000000 };
static const double ZMAT_1[3] = { 4.50036231,0.68128924,0.15137639 };
static const int NSEG_2 = 2;
static const int LVEC_2[2] = { 0,0 };
static const int NVEC_2[2] = { 2,1 };
static const double CMAT_2[3] = { 0.17525420,0.89346310,1.0000000 };
static const double ZMAT_2[3] = { 13.62310172,1.99894370,0.38294258 };
static const int NSEG_3 = 5;
static const int LVEC_3[5] = { 0,0,1,0,1 };
static const int NVEC_3[5] = { 3,2,2,1,1 };
static const double CMAT_3[9] = { 0.02718392,0.19946919,0.84471834,0.56680799,0.49693207,-0.00031685,1.00020338,1.0000000,1.0000000 };
static const double ZMAT_3[9] = { 261.70980712,39.38599444,8.90009112,2.40550517,0.70564043,2.40550517,0.70564043,0.04779874,0.04779874 };
static const int NSEG_4 = 5;
static const int LVEC_4[5] = { 0,0,1,0,1 };
static const int NVEC_4[5] = { 3,2,2,1,1 };
static const double CMAT_4[9] = { 0.02656239,0.19608914,0.84744314,0.55252358,0.50897851,-0.00939323,1.00611127,1.0000000,1.0000000 };
static const double ZMAT_4[9] = { 503.79744442,75.82361701,17.18781371,4.71330223,1.41880203,4.71330223,1.41880203,0.10454891,0.10454891 };
static const int NSEG_5 = 5;
static const int LVEC_5[5] = { 0,0,1,0,1 };
static const int NVEC_5[5] = { 3,2,2,1,1 };
static const double CMAT_5[9] = { 0.03374414,0.23775108,0.81212145,0.22649943,0.82508535,0.43916998,0.68095885,1.0000000,1.0000000 };
static const double ZMAT_5[9] = { 338.23459961,50.98058060,11.39060370,0.56716863,0.14236995,0.56716863,0.14236995,3.03182778,3.03182778 };
static const int NSEG_6 = 5;
static const int LVEC_6[5] = { 0,0,1,0,1 };
static const int NVEC_6[5] = { 3,2,2,1,1 };
static const double CMAT_6[9] = { 0.03330322,0.23617745,0.81336259,0.24008573,0.81603757,0.46214684,0.66529098,1.0000000,1.0000000 };
static const double ZMAT_6[9] = { 499.24042249,75.25419194,16.86538669,0.89739483,0.21746772,0.89739483,0.21746772,4.52660451,4.52660451 };
static const int NSEG_7 = 5;
static const int LVEC_7[5] = { 0,0,1,0,1 };
static const int NVEC_7[5] = { 3,2,2,1,1 };
static const double CMAT_7[9] = { 0.03294106,0.23480293,0.81448223,0.24876782,0.80997026,0.47669382,0.65455215,1.0000000,1.0000000 };
static const double ZMAT_7[9] = { 693.92431602,104.60302038,23.49023698,1.29068964,0.30742289,1.29068964,0.30742289,6.33456675,6.33456675 };
static const int NSEG_8 = 5;
static const int LVEC_8[5] = { 0,0,1,0,1 };
static const int NVEC_8[5] = { 3,2,2,1,1 };
static const double CMAT_8[9] = { 0.03280967,0.23422391,0.81490980,0.27389659,0.79112437,0.48155753,0.65447861,1.0000000,1.0000000 };
static const double ZMAT_8[9] = { 910.10034975,137.19711335,30.85279077,1.72885887,0.39954770,1.72885887,0.39954770,8.35065975,8.35065975 };
static const int NSEG_9 = 5;
static const int LVEC_9[5] = { 0,0,1,0,1 };
static const int NVEC_9[5] = { 3,2,2,1,1 };
static const double CMAT_9[9] = { 0.03263225,0.23347259,0.81551836,0.78100062,0.28752950,0.65148993,0.48768552,1.0000000,1.0000000 };
static const double ZMAT_9[9] = { 1162.26781448,175.21411023,39.44340117,0.50746974,2.23562751,0.50746974,2.23562751,10.69908170,10.69908170 };
static const int NSEG_10 = 5;
static const int LVEC_10[5] = { 0,0,1,0,1 };
static const int NVEC_10[5] = { 3,2,2,1,1 };
static const double CMAT_10[9] = { 0.03243674,0.23265762,0.81620225,0.77597368,0.29455064,0.64847256,0.49285246,1.0000000,1.0000000 };
static const double ZMAT_10[9] = { 1451.46006052,218.80907282,49.29660639,0.63139230,2.81447474,0.63139230,2.81447474,13.38895401,13.38895401 };
static const int NSEG_11 = 7;
static const int LVEC_11[7] = { 0,0,1,0,1,0,1 };
static const int NVEC_11[7] = { 3,3,3,2,2,1,1 };
static const double CMAT_11[15] = { 0.02330366,0.17562017,0.86488331,0.35402828,0.57357508,0.16735016,0.04620771,0.21347098,0.82763041,0.55189794,0.51011058,0.61044572,0.48877599,1.0000000,1.0000000 };
static const double ZMAT_11[15] = { 7320.01017255,1099.21475690,248.47102766,68.10517280,20.86612762,6.49460139,68.10517280,20.86612762,6.49460139,2.02508269,0.60635153,2.02508269,0.60635153,0.04204431,0.04204431 };
static const int NSEG_12 = 7;
static const int LVEC_12[7] = { 0,0,1,0,1,0,1 };
static const int NVEC_12[7] = { 3,3,3,2,2,1,1 };
static const double CMAT_12[15] = { 0.02306231,0.17358024,0.86656933,0.32783062,0.57073610,0.19859933,0.04567364,0.21085454,0.82810900,0.48410538,0.57355367,0.60906292,0.48389573,1.0000000,1.0000000 };
static const double ZMAT_12[15] = { 9468.12001271,1421.53310401,321.28960847,88.08064448,27.01701160,8.58718567,88.08064448,27.01701160,8.58718567,2.74393227,0.85790107,2.74393227,0.85790107,0.07206199,0.07206199 };
static const int NSEG_13 = 7;
static const int LVEC_13[7] = { 0,0,1,0,1,0,1 };
static const int NVEC_13[7] = { 3,3,3,2,2,1,1 };
static const double CMAT_13[15] = { 0.02289401,0.17213922,0.86776492,0.31053456,0.56662400,0.22119953,0.04518149,0.20928045,0.82842769,0.44069202,0.61258277,0.61177585,0.47588924,1.0000000,1.0000000 };
static const double ZMAT_13[15] = { 11793.72923854,1770.45604036,400.06339907,109.64007023,33.64784884,10.83737466,109.64007023,33.64784884,10.83737466,3.52543096,1.14261480,3.52543096,1.14261480,0.11329401,0.11329401 };
static const int NSEG_14 = 7;
static const int LVEC_14[7] = { 0,0,1,0,1,0,1 };
static const int NVEC_14[7] = { 3,3,3,2,2,1,1 };
static const double CMAT_14[15] = { 0.02280973,0.17143469,0.86836878,0.30616268,0.56633497,0.22599843,0.02280973,0.20971725,0.82791361,0.43766794,0.61251885,0.62814268,0.45369543,1.0000000,1.0000000 };
static const double ZMAT_14[15] = { 13955.62828382,2094.81304780,473.12990622,129.47403982,39.69162349,12.82658551,129.47403982,39.69162349,12.82658551,4.22201068,1.41747555,4.22201068,1.41747555,0.16388839,0.16388839 };
static const int NSEG_15 = 7;
static const int LVEC_15[7] = { 0,0,1,0,1,0,1 };
static const int NVEC_15[7] = { 3,3,3,2,2,1,1 };
static const double CMAT_15[15] = { 0.02286522,0.17211413,0.86787870,0.32700602,0.57438945,0.19529712,0.02286522,0.21457368,0.82515660,0.52884813,0.51967461,0.68437788,0.38812323,1.0000000,1.0000000 };
static const double ZMAT_15[15] = { 15151.94544206,2274.34114677,513.08292993,139.82036038,42.66229825,13.59162817,139.82036038,42.66229825,13.59162817,4.47317336,1.55838863,4.47317336,1.55838863,0.22034376,0.22034376 };
static const int NSEG_16 = 7;
static const int LVEC_16[7] = { 0,0,1,0,1,0,1 };
static const int NVEC_16[7] = { 3,3,3,2,2,1,1 };
static const double CMAT_16[15] = { 0.02351515,0.17838697,0.86294128,0.45398660,0.56134562,0.06479410,0.05086938,0.25016598,0.79886860,0.16040437,0.87484226,0.34297956,0.75275946,1.0000000,1.0000000 };
static const double ZMAT_16[15] = { 12715.98885260,1909.10567091,429.10047704,115.28566038,34.54429948,10.43210743,115.28566038,34.54429948,10.43210743,0.78345876,0.21432503,0.78345876,0.21432503,3.22980606,3.22980606 };
static const int NSEG_17 = 7;
static const int LVEC_17[7] = { 0,0,1,0,1,0,1 };
static const int NVEC_17[7] = { 3,3,3,2,2,1,1 };
static const double CMAT_17[15] = { 0.02344138,0.17792925,0.86335580,0.45303954,0.56317356,0.06363089,0.05090772,0.25001320,0.79862678,0.23130821,0.81374263,0.39700121,0.70223876,1.0000000,1.0000000 };
static const double ZMAT_17[15] = { 14500.35608106,2176.80649430,488.98549689,131.13866026,39.21942166,11.88944993,131.13866026,39.21942166,11.88944993,0.86996225,0.24627185,0.86996225,0.24627185,3.72213855,3.72213855 };
static const int NSEG_18 = 7;
static const int LVEC_18[7] = { 0,0,1,0,1,0,1 };
static const int NVEC_18[7] = { 3,3,3,2,2,1,1 };
static const double CMAT_18[15] = { 0.02334667,0.17728956,0.86391334,0.44744675,0.56627977,0.06672096,0.05057537,0.24944607,0.79882484,0.26333730,0.78554563,0.42617459,0.67513936,1.0000000,1.0000000 };
static const double ZMAT_18[15] = { 16573.92190582,2487.85495144,558.59858877,149.60980979,44.69492563,13.60706499,149.60980979,44.69492563,13.60706499,1.00805419,0.28776986,1.00805419,0.28776986,4.30954676,4.30954676 };


void GaussianBasis::Pople321GSP( int const atn, int & NSEG, int const * & LVEC, int const * & NVEC, double const * & CMAT, double const * & ZMAT )
{
  switch( atn )
    {
       case 1:
          NSEG = NSEG_1;
          LVEC = LVEC_1;
          NVEC = NVEC_1;
          CMAT = CMAT_1;
          ZMAT = ZMAT_1;
          break;
       case 2:
          NSEG = NSEG_2;
          LVEC = LVEC_2;
          NVEC = NVEC_2;
          CMAT = CMAT_2;
          ZMAT = ZMAT_2;
          break;
       case 3:
          NSEG = NSEG_3;
          LVEC = LVEC_3;
          NVEC = NVEC_3;
          CMAT = CMAT_3;
          ZMAT = ZMAT_3;
          break;
       case 4:
          NSEG = NSEG_4;
          LVEC = LVEC_4;
          NVEC = NVEC_4;
          CMAT = CMAT_4;
          ZMAT = ZMAT_4;
          break;
       case 5:
          NSEG = NSEG_5;
          LVEC = LVEC_5;
          NVEC = NVEC_5;
          CMAT = CMAT_5;
          ZMAT = ZMAT_5;
          break;
       case 6:
          NSEG = NSEG_6;
          LVEC = LVEC_6;
          NVEC = NVEC_6;
          CMAT = CMAT_6;
          ZMAT = ZMAT_6;
          break;
       case 7:
          NSEG = NSEG_7;
          LVEC = LVEC_7;
          NVEC = NVEC_7;
          CMAT = CMAT_7;
          ZMAT = ZMAT_7;
          break;
       case 8:
          NSEG = NSEG_8;
          LVEC = LVEC_8;
          NVEC = NVEC_8;
          CMAT = CMAT_8;
          ZMAT = ZMAT_8;
          break;
       case 9:
          NSEG = NSEG_9;
          LVEC = LVEC_9;
          NVEC = NVEC_9;
          CMAT = CMAT_9;
          ZMAT = ZMAT_9;
          break;
       case 10:
          NSEG = NSEG_10;
          LVEC = LVEC_10;
          NVEC = NVEC_10;
          CMAT = CMAT_10;
          ZMAT = ZMAT_10;
          break;
       case 11:
          NSEG = NSEG_11;
          LVEC = LVEC_11;
          NVEC = NVEC_11;
          CMAT = CMAT_11;
          ZMAT = ZMAT_11;
          break;
       case 12:
          NSEG = NSEG_12;
          LVEC = LVEC_12;
          NVEC = NVEC_12;
          CMAT = CMAT_12;
          ZMAT = ZMAT_12;
          break;
       case 13:
          NSEG = NSEG_13;
          LVEC = LVEC_13;
          NVEC = NVEC_13;
          CMAT = CMAT_13;
          ZMAT = ZMAT_13;
          break;
       case 14:
          NSEG = NSEG_14;
          LVEC = LVEC_14;
          NVEC = NVEC_14;
          CMAT = CMAT_14;
          ZMAT = ZMAT_14;
          break;
       case 15:
          NSEG = NSEG_15;
          LVEC = LVEC_15;
          NVEC = NVEC_15;
          CMAT = CMAT_15;
          ZMAT = ZMAT_15;
          break;
       case 16:
          NSEG = NSEG_16;
          LVEC = LVEC_16;
          NVEC = NVEC_16;
          CMAT = CMAT_16;
          ZMAT = ZMAT_16;
          break;
       case 17:
          NSEG = NSEG_17;
          LVEC = LVEC_17;
          NVEC = NVEC_17;
          CMAT = CMAT_17;
          ZMAT = ZMAT_17;
          break;
       case 18:
          NSEG = NSEG_18;
          LVEC = LVEC_18;
          NVEC = NVEC_18;
          CMAT = CMAT_18;
          ZMAT = ZMAT_18;
          break;
      default:
         NSEG = 0;
         LVEC = NULL;
         NVEC = NULL;
         CMAT = NULL;
         ZMAT = NULL;
         break;
    };
}

