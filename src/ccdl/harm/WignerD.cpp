#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>

#include "WignerD.hpp"
#include "../geom/RotationMatrix.hpp"

// fcoef(n) = n!
inline double fcoef( int n )
{
  double f = 1.;
  for ( int i=2; i <= n; ++i )
    f *= i;
  return f;
}

// bcoef(n,k) = n choose k
inline double bcoef( int n, int k )
{
  double bc = 0.;
  if ( k >= 0 and n-k >= 0 )
    bc = fcoef( n ) / ( fcoef( k ) * fcoef( n - k ) );
  return bc;
}

// pcoef(n) = (-1)^n
inline double pcoef( int n )
{
  return n%2 ? -1. : 1.;
}


namespace ccdl
{
  void WignerLittleD90( int j, double * d );
}


void ccdl::WignerLittleD90( int j, double * d )
{
  //
  //  d^j_{m,m'}(pi/2) 
  //         = (-1)^{m-m'} * (1/2^j)
  //          * sqrt( (j+m)!(j-m)! / ( (j+m')(j-m') ) )
  //          \sum_k  (-1)^k ( j+m' choose k ) ( j-m' choose k+m-m' )
  ///
  // NOTE:
  //
  //  d^j_{m,+m'}(pi/2) = (-1)^{ m+j } d^j_{m,-m'}(pi/2)
  //  d^j_{m,-m'}(pi/2) = (-1)^{ m+j } d^j_{m,+m'}(pi/2)
  //
  //
  //

  int nm = 2*j+1;

  for ( int mp=-j; mp <= j; ++mp )
    {
      for ( int m=-j; m <= j; ++m )
	{
	  double s = 0.;
	  for ( int k=0; k<2*j+1; ++k )
	    s += pcoef( k ) 
	      * bcoef( j+mp, k )
	      * bcoef( j-mp, k+m-mp );

	  d[(j+m)+(j+mp)*nm] 
	    = pcoef( m-mp ) * std::pow(0.5,j)
	    * std::sqrt( fcoef(j+m)*fcoef(j-m) / ( fcoef(j+mp) * fcoef(j-mp) ) )
	    * s;
	    

	};
    };

}

void ccdl::WignerLittleD( int j, double Y, double * d )
{
  //
  // d^j_{m,m'}(X) 
  //    = (-1)^{m-m'} 
  //      \sum_{m''} d^j_{m,m''}(pi/2) 
  //                * Re[ i^{2j-m-m'} e^{-im''X} ] 
  //                * d^j_{m'',-m'}(pi/2) 
  //
  //    = (-1)^{m-m'} 
  //      \sum_{m''} d^j_{m,m''}(pi/2) 
  //                * (-1)^{ m'' + j }
  //                * Re[ i^{2j-m-m'} e^{-im''X} ] 
  //                * d^j_{m'', m'}(pi/2) 
  //
  // but
  //
  //  Re[ i^{2j-m-m'} e^{-im''X} ] 
  //       = cos[(pi/2) (2j-m-m')] cos[m'' X]
  //       + sin[(pi/2) (2j-m-m')] sin[m'' X]
  //
  // so...
  //
  // d^j_{m,m'}(Y) 
  //  = (-1)^{m-m'} \sum_{m''} (-1)^{ m'' + j } 
  //        * d^j_{m,m''}(pi/2) 
  //        * (  cos[(pi/2) (2j-m-m')] cos[m'' Y] 
  //           + sin[(pi/2) (2j-m-m')] sin[m'' Y] )
  //        * d^j_{m'', m'}(pi/2) 
  //
  // The big D matrix is
  //   D^j_{m,m'}(Z,Y,Z') = e^{-imZ} d^j_{m,m'}(Y) e^{-im'Z'}
  //

  int nm = 2*j+1;
  std::vector<double> d90(nm*nm,0.);
  ccdl::WignerLittleD90( j, d90.data() );
  double const pio2 = 0.5 * M_PI;


  for ( int i=0; i<nm*nm; ++i )
    d[i] = 0.;

  for ( int mp=-j; mp <= j; ++mp )
    for ( int m=-j; m <= j; ++m )
      {
	double s = 0.;
	for ( int mpp = -j; mpp <= j; ++mpp )
	  s += pcoef(m-mp+mpp+j)
	    *  d90[(j+m)+(j+mpp)*nm]
	    * ( std::cos(pio2*(j+j-m-mp)) * std::cos(mpp*Y) +
		std::sin(pio2*(j+j-m-mp)) * std::sin(mpp*Y) )
	    * d90[(j+mpp)+(j+mp)*nm];
	d[(j+m)+(j+mp)*nm] = s;
      };

}



void ccdl::WignerBigD( int j, double Z, double Y, double Zp, double * c, double * s )
{
  int nm = 2*j+1;
  std::vector<double> d(nm*nm,0.);
  ccdl::WignerLittleD( j, Y, d.data() );
  for ( int mp = -j; mp <= j; ++mp )
    for ( int m = -j; m <= j; ++m )
      {
	c[(j+m)+(j+mp)*nm] = d[(j+m)+(j+mp)*nm] * std::cos(m*Z+mp*Zp);
	s[(j+m)+(j+mp)*nm] = d[(j+m)+(j+mp)*nm] * std::sin(m*Z+mp*Zp);
      };
}




//
// Let R be a 3x3 rotation matrix generated from (Z,Y,Zp)
// The rotation matrix is composed of intrinsic rotations, starting with Z and ending with Zp.
// r' = R.r = R(Zp).R(Y).R(Z).r
//
// D_{lm,jk} = \int Y_{lm}( R.r ) Y_{jk}( r ) dw
// Ylm(R.r) = \sum_{jk} D_{lm,jk}(Z,Y,Zp) . Yjk(r)
//
//
void ccdl::RealWignerBigD( int j, double Z, double Y, double Zp, double * D )
{
  int nm = 2*j+1;
  std::vector<double> c(nm*nm,0.);
  std::vector<double> s(nm*nm,0.);
  ccdl::WignerBigD( j, Z, Y, Zp, c.data(), s.data() );
  // Ymc = Re ((-1)^m/sqrt(2)) ( Ym + (-1)^m Y-m )  (m>0)
  // Yms = Im ((-1)^m/sqrt(2)) ( Ym - (-1)^m Y-m )  (m<0)
  // Ymc = Re ((-1)^m/sqrt(2)) ( Ym + Ym* )  (m>0)
  // Yms = Im ((-1)^m/sqrt(2)) ( Ym - Ym* )  (m<0)

  // cc
  //   Dmm' + (-1)^m' Dm-m' + (-1)^m D-mm' + (-1)^{m+m'} D-m-m'
  // cs
  //   Dmm' - (-1)^m' Dm-m' + (-1)^m D-mm' - (-1)^{m+m'} D-m-m'
  // sc
  //   Dmm' + (-1)^m' Dm-m' - (-1)^m D-mm' - (-1)^{m+m'} D-m-m'
  // ss
  //   Dmm' - (-1)^m' Dm-m' - (-1)^m D-mm' + (-1)^{m+m'} D-m-m'

  // cc
  //   Cmm' + (-1)^m' Cm-m' + (-1)^m C-mm' + (-1)^{m+m'} C-m-m'
  // cs
  //   Smm' - (-1)^m' Sm-m' + (-1)^m S-mm' - (-1)^{m+m'} S-m-m'
  // sc
  //   Smm' + (-1)^m' Sm-m' - (-1)^m S-mm' - (-1)^{m+m'} S-m-m'
  // ss
  //   Cmm' - (-1)^m' Cm-m' - (-1)^m C-mm' + (-1)^{m+m'} C-m-m'

  int imp = 0;
  for ( int mp = 0; mp <= j; ++mp )
    for ( int csp = 0; csp < 2; ++csp )
      {
	if ( csp == 1 and mp == 0 ) continue;
	//double pp = pcoef(mp) * 0.5 * std::sqrt( 2. - (mp?0.:1.) );
	double pp = 0.5 * std::sqrt( 2. - (mp?0.:1.) );

	int im = 0;
	for ( int m = 0; m <= j; ++m )
	  for ( int cs = 0; cs < 2; ++cs )
	    {
	      if ( cs == 1 and m == 0 ) continue;
	      //double p = cs ? 0.5 : ( m ? 0.5 : 1.0 );
	      //double p = pcoef(m) * 0.5 * std::sqrt( 2. - (m?0.:1.) );
	      double p = 0.5 * std::sqrt( 2. - (m?0.:1.) );

	      double v = 0.;
	      if ( cs == 0 )
		{
		  if ( csp == 0 )
		    {
		      v = 1.          * c[(j+m)+(j+mp)*nm] 
			+ pcoef(mp)   * c[(j+m)+(j-mp)*nm]
			+ pcoef(m)    * c[(j-m)+(j+mp)*nm]
			+ pcoef(m+mp) * c[(j-m)+(j-mp)*nm];
		    }
		  else
		    {
		      v = 1.          * s[(j+m)+(j+mp)*nm] 
			- pcoef(mp)   * s[(j+m)+(j-mp)*nm]
			+ pcoef(m)    * s[(j-m)+(j+mp)*nm]
			- pcoef(m+mp) * s[(j-m)+(j-mp)*nm];
		      v=-v; /// hmmmm..... why.....?
		    };
		}
	      else
		{
		  if ( csp == 0 )
		    {
		      v = 1.          * s[(j+m)+(j+mp)*nm] 
			+ pcoef(mp)   * s[(j+m)+(j-mp)*nm]
			- pcoef(m)    * s[(j-m)+(j+mp)*nm]
			- pcoef(m+mp) * s[(j-m)+(j-mp)*nm];
		    }
		  else
		    {
		      v = 1.          * c[(j+m)+(j+mp)*nm] 
			- pcoef(mp)   * c[(j+m)+(j-mp)*nm]
			- pcoef(m)    * c[(j-m)+(j+mp)*nm]
			+ pcoef(m+mp) * c[(j-m)+(j-mp)*nm];
		    };
		};

	      D[im+imp*nm] = p * pp * v ; // hmmmm, why doesn't this have factors of (-1)**(m+mp)

	      ++im;
	    };
	++imp;
      };

}



void ccdl::RealWignerBigD( int j, double const * R, double * D )
{
  double Z,Y,Zp;
  ccdl::ZYZFromRotationMatrix(R,Z,Y,Zp);
  ccdl::RealWignerBigD(j,Z,Y,Zp,D);
}
