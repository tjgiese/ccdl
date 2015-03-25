#ifndef _OrthogonalPolynomials_hpp_
#define _OrthogonalPolynomials_hpp_

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

