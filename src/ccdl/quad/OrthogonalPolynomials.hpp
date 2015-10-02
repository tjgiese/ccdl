#ifndef _OrthogonalPolynomials_hpp_
#define _OrthogonalPolynomials_hpp_

#include <cmath>

namespace ccdl
{

  void HermitePolynomial
  ( double const x, 
    int const n, 
    double * val );

  void LaguerrePolynomial
  ( double const x, 
    double const alpha, 
    int const n, 
    double * val );

  void LegendrePolynomial
  ( double const x, 
    int const n, 
    double * val );

  void LegendrePolynomial
  ( double const x, 
    int const n,
    double * val,
    double * der );

  void LegendrePolynomial
  ( double const x, 
    int const n,
    double * val,
    double * der,
    double * der2 );

  void LegendrePolynomial
  ( double const x0, 
    double const x1, 
    double const x, 
    int const n, 
    double * p );

  void LegendrePolynomial
  ( double const x0, 
    double const x1, 
    double const x, 
    int const n, 
    double * p, 
    double * dp );

  void LegendrePolynomial
  ( double const x0, 
    double const x1, 
    double const x, 
    int const n, 
    double * p, 
    double * dp,
    double * d2p );


  void OrthonormalLegendreBasis
  ( double const x0, 
    double const x1, 
    double const x, 
    int const n, 
    double * p );

  void OrthonormalLegendreBasis
  ( double const x0, 
    double const x1, 
    double const x, 
    int const n, 
    double * p, 
    double * dp );

  void OrthonormalLegendreBasis
  ( double const x0, 
    double const x1, 
    double const x, 
    int const n, 
    double * p, 
    double * dp,
    double * d2p );

  void ChebyshevPolynomialFirstKind
  ( double const x,
    int const n,
    double * val );

  void ChebyshevPolynomialSecondKind
  ( double const x,
    int const n,
    double * val );

  void JacobiPolynomial
  ( double const x, 
    double const alpha, 
    double const beta, 
    int const n, 
    double * val );

}

inline void ccdl::HermitePolynomial
( double const x, 
  int const n, // number of elements, not max order.  max_order = n-1
  double * val )
{
  val[0] = 1.0;
  if ( n > 1 )
    {
      double const tx  = 2.0 * x;
      val[1] = tx;
      for ( int i=2; i < n; ++i )
	val[i] = tx * val[i-1] - 2*(i-1) * val[i-2];
    };
}

inline void ccdl::LaguerrePolynomial
( double const x, 
  double const a, 
  int const n, 
  double * val )
{
  val[0] = 1.;
  if ( n > 1 )
    {
      val[1] = 1.-x+a;
      for ( int i=2; i<n; ++i )
	val[i] = ( (i+i+a-x-1) * val[i-1] - (i+a-1) * val[i-2] ) / i;
    };
}

inline void ccdl::LegendrePolynomial
( double const x, 
  int const n, // number of elements, not max order.  max_order = n-1
  double * val )
{
  val[0] = 1.0;
  if ( n > 1 )
    {
      val[1] = x;
      for ( int i=2; i < n; ++i )
	val[i] = ( (i+i-1) * x * val[i-1] - (i-1) * val[i-2] ) / i;
    };
}


inline
void ccdl::LegendrePolynomial
( double const x, 
  int const n, // number of elements, not max order.  max_order = n-1
  double * val,
  double * der )
{
  val[0] = 1.0;
  der[0] = 0.0;
  if ( n > 1 )
    {
      val[1] = x;
      der[1] = 1.;
      for ( int i=2; i < n; ++i )
	{
	  val[i] = ( (i+i-1) * x * val[i-1] - (i-1) * val[i-2] ) / i;
	  der[i] = ( (i+i-1) * (x * der[i-1] + val[i-1]) - (i-1) * der[i-2] ) / i;
	}
    };
}

inline
void ccdl::LegendrePolynomial
( double const x, 
  int const n, // number of elements, not max order.  max_order = n-1
  double * val,
  double * der,
  double * der2 )
{
  val[0] = 1.0;
  der[0] = 0.0;
  der2[0] = 0.0;
  if ( n > 1 )
    {
      val[1] = x;
      der[1] = 1.;
      der2[1] = 0.;
      for ( int i=2; i < n; ++i )
	{
	  val[i] = ( (i+i-1) * x * val[i-1] - (i-1) * val[i-2] ) / i;
	  der[i] = ( (i+i-1) * (x * der[i-1] + val[i-1]) - (i-1) * der[i-2] ) / i;
	  der2[i]= ( (i+i-1) * (der[i-1] + x*der2[i-1] + der[i-1]) - (i-1) * der2[i-2] ) / i;
	}
    };
}

inline
void ccdl::LegendrePolynomial( double const x0, double const x1, double const x, 
			       int const n, double * p )
{
  double u = 2. * (x-x0)/(x1-x0) - 1.;
  // du/dx = 2./(x1-x0);
  ccdl::LegendrePolynomial( u, n, p );
}


inline
void ccdl::LegendrePolynomial( double const x0, double const x1, double const x, 
			       int const n, double * p, double * dp )
{
  double u = 2. * (x-x0)/(x1-x0) - 1.;
  // du/dx = 2./(x1-x0);
  ccdl::LegendrePolynomial( u, n, p, dp );
  double const dudx = 2./(x1-x0);
  for ( int k=0; k<n; ++k )
    dp[k] *= dudx;
}


inline
void ccdl::LegendrePolynomial( double const x0, double const x1, double const x, 
			       int const n, double * p, double * dp, double * d2p )
{
  double u = 2. * (x-x0)/(x1-x0) - 1.;
  // du/dx = 2./(x1-x0);
  ccdl::LegendrePolynomial( u, n, p, dp );
  double const dudx = 2./(x1-x0);
  for ( int k=0; k<n; ++k )
    {
      dp[k] *= dudx;
      d2p[k] *= dudx*dudx;
    };
}



inline
void ccdl::OrthonormalLegendreBasis
( double const x0, double const x1, 
  double const x, int const n, double * p )
{
  //
  // Pn(x;x0,x1) = sqrt( ( 2*n+1 ) / (x1-x0) ) Pn( 2*(x-x0)/(x1-x0) )
  // \int_x0^x1 Pn(x;x0,x1) Pm(x;x0,x1) dx = \delta_n,m
  //
  double u = 2. * (x-x0)/(x1-x0) - 1.;
  // du/dx = 2./(x1-x0);
  ccdl::LegendrePolynomial( u, n, p );
  for ( int k=0; k<n; ++k )
    p[k] *= std::sqrt((2*k+1)/(x1-x0));
}


inline
void ccdl::OrthonormalLegendreBasis
( double const x0, double const x1, 
  double const x, int const n, double * p, double * dp )
{
  //
  // Pn(x;x0,x1) = sqrt( ( 2*n+1 ) / (x1-x0) ) Pn( 2*(x-x0)/(x1-x0) )
  // \int_x0^x1 Pn(x;x0,x1) Pm(x;x0,x1) dx = \delta_n,m
  //
  double const u = 2. * (x-x0)/(x1-x0) - 1.;
  double const dudx = 2./(x1-x0);
  ccdl::LegendrePolynomial( u, n, p, dp );
  for ( int k=0; k<n; ++k )
    {
      double nrm = std::sqrt((2*k+1)/(x1-x0));
      p[k] *= nrm;
      dp[k] *= nrm * dudx;
    }
}

inline
void ccdl::OrthonormalLegendreBasis
( double const x0, double const x1, 
  double const x, int const n, double * p, double * dp, double * d2p )
{
  //
  // Pn(x;x0,x1) = sqrt( ( 2*n+1 ) / (x1-x0) ) Pn( 2*(x-x0)/(x1-x0) )
  // \int_x0^x1 Pn(x;x0,x1) Pm(x;x0,x1) dx = \delta_n,m
  //
  double const u = 2. * (x-x0)/(x1-x0) - 1.;
  double const dudx = 2./(x1-x0);
  ccdl::LegendrePolynomial( u, n, p, dp, d2p );
  for ( int k=0; k<n; ++k )
    {
      double nrm = std::sqrt((2*k+1)/(x1-x0));
      p[k] *= nrm;
      dp[k] *= nrm * dudx;
      d2p[k] *= nrm * dudx * dudx;
   }
}


inline void ccdl::ChebyshevPolynomialFirstKind
( double const x, 
  int const n, // number of elements, not max order.  max_order = n-1
  double * val )
{
  val[0] = 1.0;
  if ( n > 1 )
    {
      double const tx = 2.*x;
      val[1] = x;
      for ( int i=2; i < n; ++i )
	val[i] = tx * val[i-1] - val[i-2];
    };
}


inline void ccdl::ChebyshevPolynomialSecondKind
( double const x, 
  int const n, // number of elements, not max order.  max_order = n-1
  double * val )
{
  val[0] = 1.0;
  if ( n > 1 )
    {
      double const tx = 2.*x;
      val[1] = tx;
      for ( int i=2; i < n; ++i )
	val[i] = tx * val[i-1] - val[i-2];
    };
}


inline void ccdl::JacobiPolynomial
( double const x, 
  double const alpha, 
  double const beta, 
  int const n, 
  double * val )
{
  val[0] = 1.;
  if ( n > 1 )
    {
      double const ab = alpha+beta;
      double const amb = alpha*alpha-beta*beta;
      val[1] = 0.5 * ( 2.*(alpha+1.)+(ab+2.)*(x-1.) );
      
      for ( int i=2; i<n; ++i )
	{
	  double const t = 2*i+ab;
	  double const f = 2*i*(i+ab)*(t-2.);
	  val[i] = ( (t-1.)*( t*(t-2.)*x+amb ) * val[i-1] 
		     - 2.*(i+alpha-1.)*(i+beta-1.)*t * val[i-2] )/f;
	};
    };
}



#endif

