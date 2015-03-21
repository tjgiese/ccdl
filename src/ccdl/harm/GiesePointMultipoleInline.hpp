#ifndef _GiesePointMultipoleInline_H_
#define _GiesePointMultipoleInline_H_

inline void ccdl::SolidHarm_Ilm_S( double const *__restrict__ /* crd */, double const r2, double *__restrict__ Y )
{
   Y[0] = 1. / std::sqrt(r2);
}

inline void ccdl::SolidHarm_Ilm_P( double const *__restrict__ crd, double const r2, double *__restrict__ Y )
{
   double const oor2 = 1. / r2;
   Y[0] = std::sqrt(oor2);
   double a = oor2*Y[0];
   Y[1] = a * crd[2];
   Y[2] = -a * crd[0];
   Y[3] = -a * crd[1];
}

inline void ccdl::SolidHarm_Ilm_D( double const *__restrict__ crd, double const r2, double *__restrict__ Y )
{
   double const oor2 = 1. / r2;
   Y[0] = std::sqrt(oor2);
   double a = oor2*Y[0];
   double c = a*3.*oor2;
   double b = c * crd[2];
   Y[1] =  a * crd[2];
   Y[2] = -a * crd[0];
   Y[3] = -a * crd[1];
   Y[4] =  b * crd[2] - a;
   Y[5] = -b * crd[0];
   Y[6] = -b * crd[1];
   Y[7] = -c * (crd[1]*crd[1]-crd[0]*crd[0]);
   Y[8] =  c * 2. * crd[1]*crd[0];
}
inline void ccdl::SolidHarm_dIlm_S( double const *__restrict__ Y, double *__restrict__ dIlm )
{
   dIlm[0] = Y[2];
   dIlm[1] = Y[3];
   dIlm[2] = -Y[1];
}

inline void ccdl::SolidHarm_dIlm_P( double const *__restrict__ Y, double *__restrict__ dIlm )
{
   dIlm[0] = Y[2];
   dIlm[1] = Y[3];
   dIlm[2] = -Y[1];
   dIlm[3] = Y[5];
   dIlm[4] = Y[6];
   dIlm[5] = -Y[4];
   dIlm[6] = 5.0000000000000000e-01 * (Y[7]-Y[4]);
   dIlm[7] = 5.0000000000000000e-01 * Y[8];
   dIlm[8] = -Y[5];
   dIlm[9] = 5.0000000000000000e-01 * Y[8];
   dIlm[10] = -5.0000000000000000e-01 * (Y[4]+Y[7]);
   dIlm[11] = -Y[6];
}

inline void ccdl::SolidHarm_dIlm_D( double const *__restrict__ Y, double *__restrict__ dIlm )
{
   dIlm[0] = Y[2];
   dIlm[1] = Y[3];
   dIlm[2] = -Y[1];
   dIlm[3] = Y[5];
   dIlm[4] = Y[6];
   dIlm[5] = -Y[4];
   dIlm[6] = 5.0000000000000000e-01 * (Y[7]-Y[4]);
   dIlm[7] = 5.0000000000000000e-01 * Y[8];
   dIlm[8] = -Y[5];
   dIlm[9] = 5.0000000000000000e-01 * Y[8];
   dIlm[10] = -5.0000000000000000e-01 * (Y[4]+Y[7]);
   dIlm[11] = -Y[6];
   dIlm[12] = Y[10];
   dIlm[13] = Y[11];
   dIlm[14] = -Y[9];
   dIlm[15] = 5.0000000000000000e-01 * (Y[12]-Y[9]);
   dIlm[16] = 5.0000000000000000e-01 * Y[13];
   dIlm[17] = -Y[10];
   dIlm[18] = 5.0000000000000000e-01 * Y[13];
   dIlm[19] = -5.0000000000000000e-01 * (Y[9]+Y[12]);
   dIlm[20] = -Y[11];
   dIlm[21] = 5.0000000000000000e-01 * (Y[14]-Y[10]);
   dIlm[22] = 5.0000000000000000e-01 * (Y[11]+Y[15]);
   dIlm[23] = -Y[12];
   dIlm[24] = 5.0000000000000000e-01 * (Y[15]-Y[11]);
   dIlm[25] = -5.0000000000000000e-01 * (Y[10]+Y[14]);
   dIlm[26] = -Y[13];
}

inline void ccdl::IlmInteraction_SS( double const *__restrict__ /* crd */, double const r2, double *__restrict__ T )
{
   T[0] = 1. / std::sqrt(r2);
;}

inline void ccdl::IlmInteraction_PS( double const *__restrict__ crd, double const r2, double *__restrict__ T )
{
   T[0] = 1. / std::sqrt(r2);
   double a = T[0]/r2;
   T[1] = -a * crd[2];
   a *= 2.;
   T[2] = a * crd[0];
   T[3] = a * crd[1];
}

inline void ccdl::IlmInteraction_PP( double const *__restrict__ crd, double const r2, double *__restrict__ T )
{
   double const oor2 = 1. / r2;
   T[0] = std::sqrt(oor2);
   double a = oor2*T[0];
   double c = a*3.*oor2;
   double b = c * crd[2];
   double Y4 = b * crd[2] - a;
   double Y7 = -c * (crd[1]*crd[1]-crd[0]*crd[0]);
   T[1] =    -a * crd[2];
   T[2] =  2.*a * crd[0];
   T[3] =  2.*a * crd[1];
   T[4] = -T[1]; // (0,1) sym
   T[5] = -Y4;
   T[6] =  2.*b * crd[0];
   T[7] =  2.*b * crd[1];
   T[8] = -T[2]; // (0,2) sym
   T[9] =  T[6]; // (1,2) sym
   T[10] =  2. * (Y4-Y7);
   T[11] = -4.*c * crd[1]*crd[0];
   T[12] = -T[3]; // (0,3) sym
   T[13] =  T[7]; // (1,3) sym
   T[14] =  T[11]; // (2,3) sym
   T[15] =  2. * (Y4+Y7);
}

inline void ccdl::IlmInteraction_DS( double const *__restrict__ crd, double const r2, double *__restrict__ T )
{
   double const oor2 = 1. / r2;
   T[0] = std::sqrt(oor2);
   double a = 2.*oor2*T[0];
   double c = a*3.*oor2;
   double b = c * crd[2];
   T[1] = -0.5*a * crd[2];
   T[2] =      a * crd[0];
   T[3] =      a * crd[1];
   T[4] =  0.5*(b * crd[2] - a);
   T[5] = -b * crd[0];
   T[6] = -b * crd[1];
   T[7] = -c * (crd[1]*crd[1]-crd[0]*crd[0]);
   T[8] =  c * 2. * crd[1]*crd[0];
}

inline void ccdl::PtExpInt_SS( double const *__restrict__ /* crd */, double const r2, double *__restrict__ T )
{
   T[0] = 1. / std::sqrt(r2);
;}

inline void ccdl::PtExpInt_PS( double const *__restrict__ crd, double const r2, double *__restrict__ T )
{
   T[0] = 1. / std::sqrt(r2);
   double a = T[0]/r2;
   T[1] = -a * crd[2];
   T[2] = -a * crd[0];
   T[3] = -a * crd[1];
}

inline void ccdl::PtExpInt_PP( double const *__restrict__ crd, double const r2, double *__restrict__ T )
{
   double const oor2 = 1. / r2;
   T[0] = std::sqrt(oor2);
   double a = oor2*T[0];
   double c = a*3.*oor2;
   double b = c * crd[2];
   double Y4 = b * crd[2] - a;
   double Y7 = -c * (crd[1]*crd[1]-crd[0]*crd[0]);
   T[1] =    -a * crd[2];
   T[2] =    -a * crd[0];
   T[3] =    -a * crd[1];
   T[4] = -T[1]; // (0,1) sym
   T[5] = -Y4;  // 1,1
   T[6] = -b * crd[0]; // 2,1
   T[7] = -b * crd[1]; // 3,1
   T[8] = -T[2]; // (0,2) sym
   T[9] =  T[6]; // (1,2) sym
   T[10] =  0.5 * (Y4-Y7); // 2,2
   T[11] = -c * crd[1]*crd[0]; // 3,2
   T[12] = -T[3]; // (0,3) sym
   T[13] =  T[7]; // (1,3) sym
   T[14] =  T[11]; // (2,3) sym
   T[15] =  0.5 * (Y4+Y7); // 3,3
}

inline void ccdl::PtExpInt_DS( double const *__restrict__ crd, double const r2, double *__restrict__ T )
{
   double const oor2 = 1. / r2;
   T[0] = std::sqrt(oor2);
   double a = 2.*oor2*T[0];
   double c = a*3.*oor2;
   double b = c * crd[2];
   T[1] = -0.5*a * crd[2];
   T[2] = -0.5*a * crd[0];
   T[3] = -0.5*a * crd[1];
   T[4] =  0.25*(b * crd[2] - a);
   b /= 3.4641016151377544;
   T[5] = b * crd[0];
   T[6] = b * crd[1];
   c /= 6.9282032302755088;
   T[7] = -c * (crd[1]*crd[1]-crd[0]*crd[0]);
   T[8] =  c * 2. * crd[1]*crd[0];
}

inline void ccdl::PtExpGrd_SS( double const *__restrict__ crd, double const r2, double *__restrict__ T, double *__restrict__ dT )
{

   T[0] = 1. / std::sqrt(r2);
   double a = -T[0]/r2;
   dT[0] = a * crd[0];
   dT[1] = a * crd[1];
   dT[2] = a * crd[2];
}

inline void ccdl::PtExpGrd_PS( double const *__restrict__ crd, double const r2, double *__restrict__ T, double *__restrict__ dT )
{
 double Y[5];
   double const oor2 = 1. / r2;
   T[0] = std::sqrt(oor2);
   double a = -oor2*T[0];
   double c = a*3.*oor2;
   double b = c * crd[2];
   T[1] = a * crd[2];
   T[2] = a * crd[0];
   T[3] = a * crd[1];
// these are already scaled by an extra -0.5
   Y[0] =  0.5*(b * crd[2] - a); 
   Y[1] = -b * crd[0];
   Y[2] = -b * crd[1];
// these have another -0.5 in it
   Y[3] = 0.5 * c * (crd[1]*crd[1]-crd[0]*crd[0]);
   Y[4] = -c * crd[1]*crd[0];
   *(dT)   = T[2];
   *(++dT) = T[3];
   *(++dT) = T[1];

   *(++dT) = Y[1];
   *(++dT) = Y[2];
   *(++dT) = -2.*Y[0];

   *(++dT) = (Y[0]+Y[3]);
   *(++dT) = Y[4];
   *(++dT) = Y[1];

   *(++dT) = Y[4];
   *(++dT) = (Y[0]-Y[3]);
   *(++dT) = Y[2];
}

inline void ccdl::PtExpPot_SS( double const *__restrict__ /* crd */, double const r2, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb )
{
   double T = 1. / std::sqrt(r2);
   *pa += T * *qb;
   *pb += T * *qa;
}
inline void ccdl::PtExpPot_PS( double const *__restrict__ crd, double const r2, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb )
{
   double T[4];
   ccdl::PtExpInt_PS(crd,r2,T);

   for ( int i=0; i<4; ++i )
      {
        pa[i] += T[i] * *qb;
        *pb   += T[i] * qa[i];
      };
}
inline void ccdl::PtExpPot_PP( double const *__restrict__ crd, double const r2, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb )
{
   double T[16];
   ccdl::PtExpInt_PP(crd,r2,T);

  for ( int j=0; j<4; ++j )
     for ( int i=0; i<4; ++i )
        {
          pa[i] += T[i+j*4]  * qb[j];
          pb[j] += T[i+j*4]  * qa[i];
        };
}
inline void ccdl::PtExpPot_DS( double const *__restrict__ crd, double const r2, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb )
{
   double T[9];
   ccdl::PtExpInt_DS(crd,r2,T);

   for ( int i=0; i<9; ++i )
      {
        pa[i] += T[i] * *qb;
        *pb   += T[i] * qa[i];
      };
}
inline void ccdl::PtExpPotGrd_SS( double const *__restrict__ crd, double const r2, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb, double *__restrict__ g )
{
  double T = 1. / std::sqrt(r2);
  double a = - (*qa) * (*qb) * T / r2;
  *pa += T * *qb;
  *pb += T * *qa;
  g[0] = a * crd[0];
  g[1] = a * crd[1];
  g[2] = a * crd[2];
}
inline void ccdl::PtExpPotGrd_PS( double const *__restrict__ crd, double const r2, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb, double *__restrict__ g )
{
   double T[4],Y[5];
   double const oor2 = 1. / r2;
   T[0] = std::sqrt(oor2);
   double a = -oor2*T[0];
   double c = a*3.*oor2;
   double b = c * crd[2];
   T[1] = a * crd[2];
   T[2] = a * crd[0];
   T[3] = a * crd[1];
   pa[0] += T[0]**qb;
   pa[1] += T[1]**qb;
   pa[2] += T[2]**qb;
   pa[3] += T[3]**qb;
   *pb   += qa[0]*T[0] + qa[1]*T[1] + qa[2]*T[2] + qa[3]*T[3];
// these are already scaled by an extra -0.5
   Y[0] =  0.5*(b * crd[2] - a); 
   Y[1] = -b * crd[0];
   Y[2] = -b * crd[1];
// these have another -0.5 in it
   Y[3] = 0.5 * c * (crd[1]*crd[1]-crd[0]*crd[0]);
   Y[4] = -c * crd[1]*crd[0];
   g[0] = *qb * ( qa[0]*T[2] + qa[1]*Y[1] + qa[2]*(Y[0]+Y[3]) + qa[3]*Y[4] );
   g[1] = *qb * ( qa[0]*T[3] + qa[1]*Y[2] + qa[2]*Y[4] + qa[3]*(Y[0]-Y[3]) );
   g[2] = *qb * ( qa[0]*T[1] + qa[1]*(-2.*Y[0]) + qa[2]*Y[1] + qa[3]*Y[2]  );
}

#endif

