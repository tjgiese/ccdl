#include "harm.hpp"
#include <cmath>

void ccdl::SolidHarmRlm_v1
( int const Lmax, double const * R, double * Y )
{
  Y[0] = 1.;

  int const nl1 = Lmax+1;
  double const x = R[0];
  double const y = R[1];
  double const z = R[2];
  double const r2 = x*x+y*y+z*z;

  double c;
  double co  = 1.0;
  double coo = 0.0;
  for ( int L=1; L < nl1; ++L )
    {
      int L2=L*L;
      c      = ( (2*L-1)*z*co-r2*coo ) / L2;
      coo   = co;
      co    = c;
      Y[L2] = c;
    };

  cmmo=0.0;
  cpmo=1.0;

  for ( int M=1; M < nl1; ++M )
    {
      int TM = 2*M;
      int M2 = M*M;
      double cpm =  ( y*cmmo - x*cpmo ) / TM;
      double cmm = -( y*cpmo + x*cmmo ) / TM;
      cpmo=cpm;
      cmmo=cmm;
      int i  = M2+TM;
      Y[i-1] = cpm;
      Y[i]   = cmm;

      double cpm2o = cpm;
      double cmm2o = cmm;
      double cpm2oo= 0.0;
      double cmm2oo= 0.0;
      int LPMtLMM = 0;
      for ( int L=M+1; L < nl1; ++L )
	{
	  int TL       = L+L;
	  int TLm1     = TL-1;
	  LpMtLmM += TLm1;
	  double a    = TLm1 * z / LpMtLmM;
	  double b    = r2 / LpMtLmM;
	  double cpm2 = a * cpm2o - b * cpm2oo;
	  double cmm2 = a * cmm2o - b * cmm2oo;
	  i += TLm1;
	  Y[i-1] = cpm2;
	  Y[i]   = cmm2;
	  cpm2oo   = cpm2o;
	  cmm2oo   = cmm2o;
	  cpm2o    = cpm2;
	  cmm2o    = cmm2;
	};
    };
  
}
