#ifndef _GieseMultipoleInline_H_
#define _GieseMultipoleInline_H_

#include "GieseBoysFcn.hpp"

#define GieseMultipole_PI 3.14159265358979323846264338328
#define GieseMultipole_TWO_OVER_SQRT_PI 1.12837916709551257389615890312

inline void ccdl::SolidHarm_Rlm_S( double const *__restrict__ /* crd */, double const /* r2 */, double *__restrict__ Y )
{
   Y[0] = 1.;
}

inline void ccdl::SolidHarm_Rlm_P( double const *__restrict__ crd, double const /* r2 */, double *__restrict__ Y )
{
   Y[0] = 1.;
   Y[1] = crd[2];
   Y[2] = -0.50000000000000000000 * crd[0];
   Y[3] = -0.50000000000000000000 * crd[1];
}

inline void ccdl::SolidHarm_Rlm_D( double const *__restrict__ crd, double const r2, double *__restrict__ Y )
{
   Y[0] = 1.;
   Y[1] = crd[2];
   Y[2] = -0.50000000000000000000 * crd[0];
   Y[3] = -0.50000000000000000000 * crd[1];
   Y[4] = 0.25000000000000000000 * ( 3.0 * crd[2] * crd[2] - r2 );
   Y[5] = crd[2] * Y[2];
   Y[6] = crd[2] * Y[3];
   Y[7] = 0.25000000000000000000 * ( crd[1]*Y[3] - crd[0]*Y[2] );
   Y[8] = -0.25000000000000000000 * ( crd[1]*Y[2] + crd[0]*Y[3] );
}

inline void ccdl::PrimGauExpInt_Overlap_SS( double const zab, double const *__restrict__ /* crd */, double const r2, double *__restrict__ X )
{
   X[0] = std::pow(zab/GieseMultipole_PI,1.5)*std::exp(-zab*r2);
}

inline void ccdl::PrimGauExpInt_Ewald_SS( double const sqrt_za, double const *__restrict__ /* crd */, double const r2, double *__restrict__ X )
{

  *X=0.;
  
  double const DELTAS[] = { 0.0625, 0.125, 0.25 };
  int const IOFFS[] = { 0, 80, 160 };

  double dT = sqrt_za*sqrt_za*r2;
  if ( dT < 30. )
    {
      int const nr = static_cast<int>(dT/10.);
      double const del = DELTAS[nr];
      int const ipt = static_cast<int>( dT/del + 0.5 );
      dT -= ipt*del;
      dT = -dT;
      double const dT2 = dT*dT;
      double const *__restrict__ f = GlobalGieseBoysFcnData + (ipt+IOFFS[nr]) * 18;
      double pO=*f;
      pO += *(++f) * dT;
      pO += *(++f) * dT2 / 2.;
      pO += *(++f) * dT2*dT / 6.;
      pO += *(++f) * dT2*dT2 / 24.;
      pO += *(++f) * dT2*dT2*dT / 120.;
      pO *= GieseMultipole_TWO_OVER_SQRT_PI * sqrt_za;
      *X = 1./sqrt(r2) - pO;
    };

}

inline void ccdl::PrimGauExpInt_Coulomb_SS( double const zab, double const *__restrict__ /* crd */, double const r2, double *__restrict__ X )
{

  
  double const DELTAS[] = { 0.0625, 0.125, 0.25 };
  int const IOFFS[] = { 0, 80, 160 };
  //double zab = za*zb/(za+zb);
  double dT = zab*r2;
  if ( dT < 30. )
    {
      int const nr = static_cast<int>(dT/10.);
      double const del = DELTAS[nr];
      int const ipt = static_cast<int>( dT/del + 0.5 );
      dT -= ipt*del;
      dT = -dT;
      double const dT2 = dT*dT;
      double const *__restrict__ f = GlobalGieseBoysFcnData + (ipt+IOFFS[nr]) * 18;
      *X =  *f;
      *X += *(++f) * dT;
      *X += *(++f) * dT2 / 2.;
      *X += *(++f) * dT2*dT / 6.;
      *X += *(++f) * dT2*dT2 / 24.;
      *X += *(++f) * dT2*dT2*dT / 120.;
      *X *= GieseMultipole_TWO_OVER_SQRT_PI * std::sqrt(zab);
    }
  else
    {
      *X= 1./std::sqrt(r2);
    };
}
inline void ccdl::ExpInt_UserAux_SS( double const *__restrict__ O, double const *__restrict__ /* crd */, double const /* r2 */, double *__restrict__ X )
{
   X[0] = O[0];
}

inline void ccdl::PrimGauExpPot_Overlap_SS( double const zab, double const *__restrict__ /* crd */, double const r2, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb )
{
   double X = std::pow(zab/GieseMultipole_PI,1.5)*std::exp(-zab*r2);
   *pa += *qb*X;
   *pb += *qa*X;
}

inline void ccdl::PrimGauExpPot_Ewald_SS( double const sqrt_za,  double const *__restrict__ /* crd */, double const r2, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb )
{
  
  double const DELTAS[] = { 0.0625, 0.125, 0.25 };
  int const IOFFS[] = { 0, 80, 160 };

  double dT = sqrt_za*sqrt_za*r2;
  if ( dT < 30. )
    {
      int const nr = static_cast<int>(dT/10.);
      double const del = DELTAS[nr];
      int const ipt = static_cast<int>( dT/del + 0.5 );
      dT -= ipt*del;
      dT = -dT;
      double const dT2 = dT*dT;
      double const *__restrict__ f = GlobalGieseBoysFcnData + (ipt+IOFFS[nr]) * 18;
      double pO=*f;
      pO += *(++f) * dT;
      pO += *(++f) * dT2 / 2.;
      pO += *(++f) * dT2*dT / 6.;
      pO += *(++f) * dT2*dT2 / 24.;
      pO += *(++f) * dT2*dT2*dT / 120.;
      pO *= GieseMultipole_TWO_OVER_SQRT_PI * sqrt_za;
      pO = 1./sqrt(r2) - pO;
      *pa += *qb*pO;
      *pb += *qa*pO;
    };
}

inline void ccdl::PrimGauExpPot_Coulomb_SS( double const zab,  double const *__restrict__ /* crd */, double const r2, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb )
{
  
  double const DELTAS[] = { 0.0625, 0.125, 0.25 };
  int const IOFFS[] = { 0, 80, 160 };
  //double const zab = za*zb/(za+zb);
  double dT = zab*r2;
  if ( dT < 30. )
    {
      int const nr = static_cast<int>(dT/10.);
      double const del = DELTAS[nr];
      int const ipt = static_cast<int>( dT/del + 0.5 );
      dT -= ipt*del;
      dT = -dT;
      double const dT2 = dT*dT;
      double const *__restrict__ f = GlobalGieseBoysFcnData + (ipt+IOFFS[nr]) * 18;
      double pO=*f;
      pO += *(++f) * dT;
      pO += *(++f) * dT2 / 2.;
      pO += *(++f) * dT2*dT / 6.;
      pO += *(++f) * dT2*dT2 / 24.;
      pO += *(++f) * dT2*dT2*dT / 120.;
      pO *= GieseMultipole_TWO_OVER_SQRT_PI * std::sqrt(zab);
      *pa += *qb*pO;
      *pb += *qa*pO;
    };
}

inline void ccdl::ExpPot_UserAux_SS( double const *__restrict__ O, double const *__restrict__ /* crd */, double const /* r2 */, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb )
{
   *pa += *qb**O;
   *pb += *qa**O;
}

inline void ccdl::PrimGauExpGrd_Overlap_SS( double const zab, double const *__restrict__ crd, double const r2, double *__restrict__ X, double *__restrict__ G )
{

   X[0] = std::pow(zab/GieseMultipole_PI,1.5)*std::exp(-zab*r2);
   double X1 =  -2.*zab*X[0];
   G[0] = X1 * crd[0];
   G[1] = X1 * crd[1];
   G[2] = X1 * crd[2];
}

inline void ccdl::PrimGauExpGrd_Ewald_SS( double const sqrt_za, double const *__restrict__ crd, double const r2, double *__restrict__ X, double *__restrict__ G )
{

   double T[2];
   double const oor = 1./std::sqrt(r2);
   double const za = sqrt_za*sqrt_za;
   double t = za*r2;
   GieseBoysFcnTemplate<1>(t,T);
   t = GieseMultipole_TWO_OVER_SQRT_PI * sqrt_za;
   *X = oor - *T*t;
   T[1] = -oor/r2+2.*za*T[1]*t;
   G[0] = T[1] * crd[0];
   G[1] = T[1] * crd[1];
   G[2] = T[1] * crd[2];
}

inline void ccdl::PrimGauExpGrd_Coulomb_SS( double const zab, double const *__restrict__ crd, double const r2, double *__restrict__ X, double *__restrict__ G )
{

   double T[2];
   //double const zab = za*zb/(za+zb);
   double t = zab*r2;
   GieseBoysFcnTemplate<1>(t,T);
   t = GieseMultipole_TWO_OVER_SQRT_PI * std::sqrt(zab);
   *X = *T*t;
   T[1] = -2.*zab*T[1]*t;
   G[0] = T[1] * crd[0];
   G[1] = T[1] * crd[1];
   G[2] = T[1] * crd[2];
}

inline void ccdl::ExpGrd_UserAux_SS( double const *__restrict__ O, double const *__restrict__ crd, double const /* r2 */, double *__restrict__ X, double *__restrict__ G )
{

   X[0] = O[0];
   double X1 =  2. * O[1]; // ?
   G[0] = X1 * crd[0];
   G[1] = X1 * crd[1];
   G[2] = X1 * crd[2];
}

inline void ccdl::PrimGauExpPotGrd_Overlap_SS( double const zab, double const *__restrict__ crd, double const r2, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb, double *__restrict__ G )
{
   double X = std::pow(zab/GieseMultipole_PI,1.5)*std::exp(-zab*r2);
   *pa += *qb*X;
   *pb += *qa*X;
   X *=  -*qa**qb*2.*zab;
   G[0] = X * crd[0];
   G[1] = X * crd[1];
   G[2] = X * crd[2];
}

inline void ccdl::PrimGauExpPotGrd_Ewald_SS( double const sqrt_za, double const *__restrict__ crd, double const r2, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb, double *__restrict__ G )
{

   double T[2];
   double const oor = 1./std::sqrt(r2);
   double const za = sqrt_za*sqrt_za;
   double t = za*r2;
   GieseBoysFcnTemplate<1>(t,T);
   t = GieseMultipole_TWO_OVER_SQRT_PI * sqrt_za;
   *T = oor - *T*t;
   *pa += *qb**T;
   *pb += *qa**T;
   T[1] = (-oor/r2+2.*za*T[1]*t)**qa**qb;
   G[0] = T[1] * crd[0];
   G[1] = T[1] * crd[1];
   G[2] = T[1] * crd[2];
}

inline void ccdl::PrimGauExpPotGrd_Coulomb_SS( double const zab, double const *__restrict__ crd, double const r2, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb, double *__restrict__ G )
{

   double T[2];
   //double const zab = za*zb/(za+zb);
   double t = zab*r2;
   GieseBoysFcnTemplate<1>(t,T);
   t = GieseMultipole_TWO_OVER_SQRT_PI * std::sqrt(zab);
   *T *= t;
   *pa += *qb**T;
   *pb += *qa**T;
   T[1] = (-2.*zab*T[1]*t)**qa**qb;
   G[0] = T[1] * crd[0];
   G[1] = T[1] * crd[1];
   G[2] = T[1] * crd[2];
}

inline void ccdl::ExpPotGrd_UserAux_SS( double const *__restrict__ O, double const *__restrict__ crd, double const /* r2 */, double const *__restrict__ qa, double const *__restrict__ qb, double *__restrict__ pa, double *__restrict__ pb, double *__restrict__ G )
{
   double X = O[0];
   *pa += *qb*X;
   *pb += *qa*X;
   X =  *qa**qb*2.*O[1];
   G[0] = X * crd[0];
   G[1] = X * crd[1];
   G[2] = X * crd[2];
}


#endif

