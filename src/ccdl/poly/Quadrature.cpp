#include "Quadrature.hpp"
#include "../constants.hpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdlib>

namespace ccdl
{

  void arth
  ( double const first, 
    double const increment,
    int const n,
    double * rVec )
  {
    int const NPAR_ARTH = 16;
    int const NPAR2_ARTH = 8;
    
    if ( n > 0 )
      rVec[0] = first;
    
    if ( n <= NPAR_ARTH )
      {
	for ( int k=1; k < n; ++k )
	  rVec[k] = rVec[k-1] + increment;
      }
    else
      {
	for ( int k=1; k < NPAR2_ARTH; ++k )
	  rVec[k] = rVec[k-1] + increment;
	double temp = increment * NPAR2_ARTH;
	int k = NPAR2_ARTH;
	
	while (1)
	  {
	    if ( k >= n-1 )
	      {
		return;
	      }
	    else
	      {
		int k2 = k + k;
		int nk = n-k;
		int m  = std::min(k,nk);
		for ( int i=0; i<m; ++i )
		  rVec[k+i] = temp + rVec[i];
		temp = temp + temp;
		k = k2;
	      };
	  };
      };
    
  }
  
  double gammln( double const xx )
  {
    double const stp = 2.5066282746310005;
    double const coef[] = { 76.18009172947146, -86.50532032941677,
			    24.01409824083091, -1.231739572450155,
			    0.1208650973866179e-2, -0.5395239384953e-5 };
    
    if ( xx <= 0.0 )
      std::cerr << "ccdl::gammln xx <= 0.0" << std::endl;
    
    double const x = xx;
    double tmp = x + 5.5;
    tmp = (x+0.5)*std::log(tmp)-tmp;
    
    double arth[6];
    
    ccdl::arth(x+1.0,1.0,6,arth);
    double sum = 0.0;
    for ( int i=0; i<6; ++i )
      sum += coef[i] / arth[i];
    
    double g = tmp + std::log( stp * ( 1.000000000190015 + sum ) / x );
    return g; 
  }

}



void ccdl::GaussLaguerreRule
( double const rAlf,
  int const n,
  double * rVecX, 
  double * rVecW )
{
  int const MAXIT=1000;
  double const EPS = 3.0e-13;
  double const C1 = 9.084064e-01;
  double const C2 = 5.214976e-02;
  double const C3 = 2.579930e-03;
  double const C4 = 3.986126e-03;
  double const C13 = 1.0/3.0;
  

  double const anu = 4.0 * n + 2.0 * rAlf + 2.0;
  double const PiOverAnu = ccdl::PI / anu;
  std::vector<double> rhs(n);
  std::vector<double> r3(n);
  std::vector<double> r2(n);
  std::vector<double> theta(n);
  std::vector<double> z(n);

  rhs[0] = 4.0*n-1.0;
  for ( int i=1; i<n; ++i )
    rhs[i] = rhs[i-1] - 4.0;
  for ( int i=0; i<n; ++i )
    rhs[i] *= PiOverAnu;
  for ( int i=0; i<n; ++i )
    r3[i] = std::pow(rhs[i],C13);
  for ( int i=0; i<n; ++i )
    r2[i] = r3[i]*r3[i];
  for ( int i=0; i<n; ++i )
    theta[i] = r3[i]*(C1+r2[i]*(C2+r2[i]*(C3+r2[i]*C4)));
  for ( int i=0; i<n; ++i )
    {
      double const costh = std::cos(theta[i]);
      z[i] = anu * costh * costh;
    };

  std::vector<double> p1(n);
  std::vector<double> p2(n);
  std::vector<double> p3(n);
  std::vector<double> pp(n);
  std::vector<double> z1(n);
  std::vector<bool> unfinished(n);

  for ( int i=0; i<n; ++i )
    unfinished[i] = true;

  for ( int its=0; its<MAXIT; ++its )
    {
      for ( int i=0; i<n; ++i )
	{
	  if ( unfinished[i] )
	    {
	      p1[i] = 1.0;
	      p2[i] = 0.0;
	    };
	};

      for ( int j=0; j<n; ++j )
	{
	  for ( int i=0; i<n; ++i )
	    {
	      if ( unfinished[i] )
		{
		  p3[i] = p2[i];
		  p2[i] = p1[i];
		  p1[i] = ( 
			   ( 2.0*(j+1.0)-1.0+rAlf-z[i] ) * p2[i]
			   -( j+rAlf )*p3[i]
			    )/(j+1);
		};
	    };
	};

      bool HasUnfinished = false;
      for ( int i=0; i<n; ++i )
	{
	  if ( unfinished[i] )
	    {
	      pp[i] = (n*p1[i]-(n+rAlf)*p2[i])/z[i];
	      z1[i] = z[i];
	      z[i] = z1[i]-p1[i]/pp[i];

	      unfinished[i] = ( std::abs(z[i]-z1[i]) > EPS*z[i] );
	      if ( unfinished[i] )
		{
		  HasUnfinished = true;
		};
	    };
	};

      if ( ! HasUnfinished ){ break; };
      if ( its == MAXIT )
	std::cerr << "Too many iterations in ccdl::GaussLaguerreRule" << std::endl;
    };
  
  for ( int i=0; i<n; ++i )
    rVecX[i] = z[i];
  for ( int i=0; i<n; ++i )
    rVecW[i] = -std::exp( ccdl::gammln(rAlf+n) - ccdl::gammln( (double)n ) )
      / ( pp[i]*p2[i]*n );
}











void ccdl::GaussLegendreRule
( double const x1,
  double const x2,
  int const n,
  double * rVecX, 
  double * rVecW )
{
  int const MAXIT=10;
  double const EPS = 3.0e-14;
  int const m = (n+1)/2;
  double const xm = 0.5 * ( x2 + x1 );
  double const xl = 0.5 * ( x2 - x1 );  
  double const C1 = ccdl::PI / (n+0.5);
  
  std::vector<double> z(m);
  std::vector<double> z1(m);
  std::vector<double> p1(m);
  std::vector<double> p2(m);
  std::vector<double> p3(m);
  std::vector<double> pp(m);
  std::vector<bool> unfinished(m);

  for ( int i=0; i<m; ++i )
    z[i] = std::cos( C1 * ( i + 0.75 ) );
  for ( int i=0; i<m; ++i )
    unfinished[i] = true;
  for ( int its=0; its < MAXIT; ++its )
    {
      for ( int i=0; i < m; ++i )
	{
	  if ( unfinished[i] )
	    {
	      p1[i] = 1.0;
	      p2[i] = 0.0;
	    };
	};
      for ( int j=0; j < n; ++j )
	{
	  for ( int i=0; i < m; ++i )
	    {
	      if ( unfinished[i] )
		{
		  p3[i]=p2[i];
		  p2[i]=p1[i];
		  p1[i]=((2.0*j+1.0)*z[i]*p2[i]-j*p3[i])/(j+1.0);
		};
	    };
	};

      bool HasUnfinished = false;
      for ( int i=0; i < m; ++i )
	{
	  if ( unfinished[i] )
	    {
	      pp[i]=n*(z[i]*p1[i]-p2[i])/(z[i]*z[i]-1.0);
	      z1[i]=z[i];
	      z[i]=z1[i]-p1[i]/pp[i];
	      unfinished[i] = ( std::abs(z[i]-z1[i]) > EPS );
	      if ( unfinished[i] )
		HasUnfinished = true;
	    };
	};
      if ( ! HasUnfinished ){ break; };
      
      if ( its == MAXIT )
	std::cerr << "Too many iterations in ccdl::GaussLegendreRule" << std::endl;
    };

  for ( int i=0; i < m; ++i )
    rVecX[i] = xm-xl*z[i];
  for ( int i=0; i < n-m; ++i )
    rVecX[n-i-1] = xm+xl*z[i];
  for ( int i=0; i < m; ++i )
    rVecW[i] = 2.0*xl/((1.-z[i]*z[i])*pp[i]*pp[i]);
  for ( int i=0; i < n-m; ++i )
    rVecW[n-i-1] = rVecW[i];
}








void ccdl::GaussHermiteRule
( int const n,
  double * rVecX, 
  double * rVecW )
{
  int const MAXIT=10;
  double const EPS = 3.0e-13;
  double const PIM4 = 0.7511255444649425;
  double const C1 = 9.084064e-01;
  double const C2 = 5.214976e-02;
  double const C3 = 2.579930e-03;
  double const C4 = 3.986126e-03;
  double const anu = 2.0 * n + 1.0;
  int const m = (n+1)/2;

  std::vector<double> rhs(m);
  std::vector<double> r2(m);
  std::vector<double> r3(m);
  std::vector<double> theta(m);
  std::vector<double> p1(m);
  std::vector<double> p2(m);
  std::vector<double> p3(m);
  std::vector<double> pp(m);
  std::vector<double> z(m);
  std::vector<double> z1(m);
  std::vector<bool> unfinished(m);

  double const PiOverAnu = ccdl::PI/anu;
  ccdl::arth(3.,4.,m,rhs);
  for ( int i=0; i < m; ++i )
    rhs[i] *= PiOverAnu;
  double const C13 = 1.0/3.0;
  for ( int i=0; i < m; ++i )
    r3[i] = std::pow(rhs[i],C13);
  for ( int i=0; i < m; ++i )
    r2[i] = r3[i]*r3[i];
  for ( int i=0; i < m; ++i )
    theta[i] = r3[i]*(C1+r2[i]*(C2+r2[i]*(C3+r2[i]*C4)));
  double const SqrtAnu = std::sqrt(anu);
  for ( int i=0; i < m; ++i )
    z[i] = SqrtAnu * std::cos(theta[i]);
  for ( int i=0; i < m; ++i )
    unfinished[i] = true;

  for ( int its=0; its < MAXIT; ++its )
    {
      for ( int i=0; i < m; ++i )
	{
	  if ( unfinished[i] )
	    {
	      p1[i] = PIM4;
	      p2[i] = 0.0;
	    };
	};
      
      for ( int j=0; j < n; ++j )
	{
	  double const F1 = std::sqrt(2.0/(j+1.0));
	  double const F2 = std::sqrt( (double)j/(j+1.0) );
	  for ( int i=0; i < m; ++i )
	    {
	      if ( unfinished[i] )
		{
		  p3[i] = p2[i];
		  p2[i] = p1[i];
		  p1[i] = z[i] * F1 * p2[i] - F2 * p3[i];
		};
	    };
	};
      
      bool HasUnfinished = false;
      for ( int i=0; i < m; ++i )
	{
	  if ( unfinished[i] )
	    {
	      pp[i] = std::sqrt(2.0*n) * p2[i];
	      z1[i] = z[i];
	      z[i] = z1[i]-p1[i]/pp[i];
	      unfinished[i] = ( std::abs(z[i]-z1[i]) > EPS );
	      if ( unfinished[i] )
		HasUnfinished = true;
	    };
	};
      
      if ( ! HasUnfinished ){ break; };
      if ( its == MAXIT-1 )
	std::cerr << "Too many iterations in GaussHermiteRule" << std::endl;
    };
  
  for ( int i=0; i < m; ++i )
    rVecX[i] = z[i];
  for ( int i=0; i < n-m; ++i )
    rVecX[n-i-1] = -z[i];
  for ( int i=0; i < m; ++i )
    rVecW[i] = 2.0 / ( pp[i] * pp[i] );
  for ( int i=0; i < n-m; ++i )
    rVecW[n-i-1] = rVecW[i];
}







void ccdl::GaussJacobiRule
( double const alf,
  double const bet,
  int const n,
  double * rVecX, 
  double * rVecW )
{
  double const EPS=3.0e-14;
  int const MAXIT=10;
  double const alfbet = alf+bet;
  
  std::vector<double> b(n);
  std::vector<double> p1(n);
  std::vector<double> p2(n);
  std::vector<double> p3(n);
  std::vector<double> pp(n);
  std::vector<double> z(n);
  std::vector<double> z1(n);
  std::vector<bool> unfinished(n);

  for ( int i=0; i < n; ++i )
    z[i] = std::cos(
		    (ccdl::PI * (i+0.75 + 0.5*alf))
		    /
		    (n+0.5*(alfbet+1.0))
		    );
  
  for ( int i=0; i < n; ++i )
    unfinished[i] = true;
  
  double temp;
  for ( int its=0; its < MAXIT; ++its )
    {
      temp = 2.0 + alfbet; 
      for ( int i=0; i < n; ++i )
	{
	  if ( unfinished[i] )
	    {
	      p1[i] = (alf-bet+temp*z[i])/2.0;
	      p2[i] = 1.0;
	    };
	};
      
      for ( int j=1; j < n; ++j )
	{
	  double a = 2. * (j+1.) * (j+1.+alfbet) * temp;
	  temp += 2.;
	  double c = 2. * (j+alf)*(j+bet)*temp;
	  
	  for ( int i=0; i < n; ++i )
	    {
	      if ( unfinished[i] )
		{
		  p3[i]=p2[i];
		  p2[i]=p1[i];
		  b[i] = (temp-1.)*(alf*alf-bet*bet+temp*(temp-2.)*z[i]);
		  p1[i] = (b[i]*p2[i]-c*p3[i])/a;
		};
	    };
	};
      
      bool HasUnfinished = false;
      for ( int i=0; i < n; ++i )
	{
	  if ( unfinished[i] )
	    {
	      pp[i] = (n*(alf-bet-temp*z[i])*p1[i]+2.*(n+alf)*(n+bet)*p2[i])/(temp*(1.-z[i]*z[i]));
	      z1[i] = z[i];
	      z[i]=z1[i]-p1[i]/pp[i];
	      unfinished[i] = (std::abs(z[i]-z1[i]) > EPS);
	      if ( unfinished[i] )
		HasUnfinished = true;
	    };
	};
      
      if ( ! HasUnfinished ){ break; };
      
      if ( its == MAXIT-1 )
	std::cerr << "Too many iterations in GaussJacobiRule" << std::endl;
    };
  
  for ( int i=0; i < n; ++i )
    rVecX[i] = z[i];
  for ( int i=0; i < n; ++i )
    rVecW[i] = std::exp( 
			ccdl::gammln(alf+n)
			+ ccdl::gammln(bet+n)
			- ccdl::gammln(n+1.0)
			- ccdl::gammln(n+alf+bet+1.)
			 ) * temp*std::pow(2.0,alfbet)/(pp[i]*p2[i]);
}




namespace ccdl
{
  namespace Backend
  {
    
    
    /** \private */
    void _Lebedev_6_pt( std::tr1::array<double,3> * rVecPt,
			double * rVecWt );
    
    /** \private */
    void _Lebedev_14_pt( std::tr1::array<double,3> * rVecPt,
			 double * rVecWt );
    
    /** \private */
    void _Lebedev_26_pt( std::tr1::array<double,3> * rVecPt,
			 double * rVecWt );
    
    /** \private */
    void _Lebedev_38_pt( std::tr1::array<double,3> * rVecPt,
			 double * rVecWt );
    
    /** \private */
    void _Lebedev_50_pt( std::tr1::array<double,3> * rVecPt,
			 double * rVecWt );
    
    /** \private */
    void _Lebedev_74_pt( std::tr1::array<double,3> * rVecPt,
			 double * rVecWt );
    
    /** \private */
    void _Lebedev_86_pt( std::tr1::array<double,3> * rVecPt,
			 double * rVecWt );
    
    /** \private */
    void _Lebedev_110_pt( std::tr1::array<double,3> * rVecPt,
			  double * rVecWt );
    
    /** \private */
    void _Lebedev_146_pt( std::tr1::array<double,3> * rVecPt,
			  double * rVecWt );

    /** \private */
    void _Lebedev_170_pt( std::tr1::array<double,3> * rVecPt,
			  double * rVecWt );

    /** \private */
    void _Lebedev_194_pt( std::tr1::array<double,3> * rVecPt,
			  double * rVecWt );

    /** \private */
    void _Lebedev_230_pt( std::tr1::array<double,3> * rVecPt,
			  double * rVecWt );

    /** \private */
    void _Lebedev_266_pt( std::tr1::array<double,3> * rVecPt,
			  double * rVecWt );

    /** \private */
    void _Lebedev_302_pt( std::tr1::array<double,3> * rVecPt,
			  double * rVecWt );

    /** \private */
    void _Lebedev_350_pt( std::tr1::array<double,3> * rVecPt,
			  double * rVecWt );

    /** \private */
    void _Lebedev_434_pt( std::tr1::array<double,3> * rVecPt,
			  double * rVecWt );

    /** \private */
    void _Lebedev_590_pt( std::tr1::array<double,3> * rVecPt,
			  double * rVecWt );

    /** \private */
    void _Lebedev_770_pt( std::tr1::array<double,3> * rVecPt,
			  double * rVecWt );

    /** \private */
    void _Lebedev_974_pt( std::tr1::array<double,3> * rVecPt,
			  double * rVecWt );

    /** \private */
    void _Lebedev_1202_pt( std::tr1::array<double,3> * rVecPt,
			   double * rVecWt );

    /** \private */
    void _Lebedev_1454_pt( std::tr1::array<double,3> * rVecPt,
			   double * rVecWt );

    /** \private */
    void _Lebedev_1730_pt( std::tr1::array<double,3> * rVecPt,
			   double * rVecWt );

    /** \private */
    void _Lebedev_2030_pt( std::tr1::array<double,3> * rVecPt,
			   double * rVecWt );

    /** \private */
    void _Lebedev_2354_pt( std::tr1::array<double,3> * rVecPt,
			   double * rVecWt );

    /** \private */
    void _Lebedev_2702_pt( std::tr1::array<double,3> * rVecPt,
			   double * rVecWt );

    /** \private */
    void _Lebedev_3074_pt( std::tr1::array<double,3> * rVecPt,
			   double * rVecWt );

    /** \private */
    void _Lebedev_3470_pt( std::tr1::array<double,3> * rVecPt,
			   double * rVecWt );

    /** \private */
    void _Lebedev_3890_pt( std::tr1::array<double,3> * rVecPt,
			   double * rVecWt );

    /** \private */
    void _Lebedev_4334_pt( std::tr1::array<double,3> * rVecPt,
			   double * rVecWt );

    /** \private */
    void _Lebedev_4802_pt( std::tr1::array<double,3> * rVecPt,
			   double * rVecWt );

    /** \private */
    void _Lebedev_5294_pt( std::tr1::array<double,3> * rVecPt,
			   double * rVecWt );

    /** \private */
    void _Lebedev_5810_pt( std::tr1::array<double,3> * rVecPt,
			   double * rVecWt );


    /** \private */
    void _Lebedev_Oh_Mesh( unsigned int const * const nk , 
			   double const * const pPa ,
			   std::tr1::array<double,3> * pPt ,
			   double * pWt );
    
    /** \private */
    void _8pt_helper( std::tr1::array<double,3> * pt, 
		      double const rr, 
		      double const ss, 
		      double const tt );


  }
}



void ccdl::LebedevRule
( int const iRule,
  std::tr1::array<double,3> * rVecPt,
  double * rVecWt )
{

  switch ( iRule )
    {
    case 6:
      ccdl::Backend::_Lebedev_6_pt( rVecPt, rVecWt );
      break;

    case 14:
      ccdl::Backend::_Lebedev_14_pt( rVecPt, rVecWt );
      break;

    case 26:
      ccdl::Backend::_Lebedev_26_pt( rVecPt, rVecWt );
      break;

    case 38:
      ccdl::Backend::_Lebedev_38_pt( rVecPt, rVecWt );
      break;

    case 50:
      ccdl::Backend::_Lebedev_50_pt( rVecPt, rVecWt );
      break;

    case 74:
      ccdl::Backend::_Lebedev_74_pt( rVecPt, rVecWt );
      break;

    case 86:
      ccdl::Backend::_Lebedev_86_pt( rVecPt, rVecWt );
      break;

    case 110:
      ccdl::Backend::_Lebedev_110_pt( rVecPt, rVecWt );
      break;

    case 146:
      ccdl::Backend::_Lebedev_146_pt( rVecPt, rVecWt );
      break;

    case 170:
      ccdl::Backend::_Lebedev_170_pt( rVecPt, rVecWt );
      break;

    case 194:
      ccdl::Backend::_Lebedev_194_pt( rVecPt, rVecWt );
      break;

    case 230:
      ccdl::Backend::_Lebedev_230_pt( rVecPt, rVecWt );
      break;

    case 266:
      ccdl::Backend::_Lebedev_266_pt( rVecPt, rVecWt );
      break;

    case 302:
      ccdl::Backend::_Lebedev_302_pt( rVecPt, rVecWt );
      break;

    case 350:
      ccdl::Backend::_Lebedev_350_pt( rVecPt, rVecWt );
      break;

    case 434:
      ccdl::Backend::_Lebedev_434_pt( rVecPt, rVecWt );
      break;

    case 590:
      ccdl::Backend::_Lebedev_590_pt( rVecPt, rVecWt );
      break;

    case 770:
      ccdl::Backend::_Lebedev_770_pt( rVecPt, rVecWt );
      break;

    case 974:
      ccdl::Backend::_Lebedev_974_pt( rVecPt, rVecWt );
      break;

    case 1202:
      ccdl::Backend::_Lebedev_1202_pt( rVecPt, rVecWt );
      break;

    case 1454:
      ccdl::Backend::_Lebedev_1454_pt( rVecPt, rVecWt );
      break;

    case 1730:
      ccdl::Backend::_Lebedev_1730_pt( rVecPt, rVecWt );
      break;

    case 2030:
      ccdl::Backend::_Lebedev_2030_pt( rVecPt, rVecWt );
      break;

    case 2354:
      ccdl::Backend::_Lebedev_2354_pt( rVecPt, rVecWt );
      break;

    case 2702:
      ccdl::Backend::_Lebedev_2702_pt( rVecPt, rVecWt );
      break;

    case 3074:
      ccdl::Backend::_Lebedev_3074_pt( rVecPt, rVecWt );
      break;

    case 3470:
      ccdl::Backend::_Lebedev_3470_pt( rVecPt, rVecWt );
      break;

    case 3890:
      ccdl::Backend::_Lebedev_3890_pt( rVecPt, rVecWt );
      break;

    case 4334:
      ccdl::Backend::_Lebedev_4334_pt( rVecPt, rVecWt );
      break;

    case 4802:
      ccdl::Backend::_Lebedev_4802_pt( rVecPt, rVecWt );
      break;

    case 5294:
      ccdl::Backend::_Lebedev_5294_pt( rVecPt, rVecWt );
      break;

    case 5810:
      ccdl::Backend::_Lebedev_5810_pt( rVecPt, rVecWt );
      break;


    default:

      std::cerr << "Invalid LebedevRule used in ccdl_AngularQuadratureMod.cpp " << iRule << std::endl;
      std::cerr << "The following rules are OK" << std::endl
		<< std::setw(6) << 6
		<< std::setw(6) << 14
		<< std::setw(6) << 26
		<< std::setw(6) << 38
		<< std::setw(6) << 50
		<< std::setw(6) << 86
		<< std::setw(6) << 110
		<< std::setw(6) << 146
		<< std::endl
		<< std::setw(6) << 170
		<< std::setw(6) << 194
		<< std::setw(6) << 302
		<< std::setw(6) << 350
		<< std::setw(6) << 434
		<< std::setw(6) << 590
		<< std::setw(6) << 770
		<< std::setw(6) << 974
		<< std::endl
		<< std::setw(6) << 1202
		<< std::setw(6) << 1454
		<< std::setw(6) << 1730
		<< std::setw(6) << 2030
		<< std::setw(6) << 2354
		<< std::setw(6) << 2702
		<< std::endl
		<< std::setw(6) << 3074
		<< std::setw(6) << 3470
		<< std::setw(6) << 3890
		<< std::setw(6) << 4334
		<< std::setw(6) << 4802
		<< std::setw(6) << 5294
		<< std::setw(6) << 5810
		<< std::endl
		<< "The following rules are recommended" << std::endl
		<< std::setw(6) << 14
		<< std::setw(6) << 38
		<< std::setw(6) << 86
		<< std::setw(6) << 110
		<< std::setw(6) << 302
		<< std::setw(6) << 350
		<< std::setw(6) << 590
		<< std::endl;

      assert( false );
      abort();
    };

  // normalization

  double wtsum = 0.;
  for ( int i=0; i<iRule; ++i )
    wtsum += rVecWt[i];
  wtsum *= ccdl::FOUR_PI;
  for ( int i=0; i<iRule; ++i )
    rVecWt[i] *= wtsum;
}


void ccdl::Backend::_Lebedev_6_pt( std::tr1::array<double,3> * rVecPt, double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 3, 1, 0, 0, 0, 0, 0, 6 };
  double pa[] = { 0.1666666666666667e+00 };
  //rVecPt.resize( 18 );
  //rVecWt.resize( 6 );
  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );
}


void ccdl::Backend::_Lebedev_14_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 5, 1, 1, 0, 0, 0, 0, 14 };
  double pa[] = { 0.6666666666666667e-01, 0.7500000000000000e-01 };
  //rVecPt.resize( 42 );
  //rVecWt.resize( 14 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_26_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 7, 1, 1, 1, 0, 0, 0, 26 };
  double pa[] = { 0.4761904761904762e-01, 0.3214285714285714e-01, 0.3809523809523810e-01 };
  //rVecPt.resize( 78 );
  //rVecWt.resize( 26 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_38_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 9, 1, 1, 0, 0, 1, 0, 38 };
  double pa[] = { 0.9523809523809524e-02, 0.3214285714285714e-01, 0.4597008433809831e+00, 0.2857142857142857e-01 };
  //rVecPt.resize( 114 );
  //rVecWt.resize( 38 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_50_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 11, 1, 1, 1, 1, 0, 0, 50 };
  double pa[] = { 0.1269841269841270e-01, 0.2109375000000000e-01, 0.2257495590828924e-01, 0.3015113445777636e+00, 0.2017333553791887e-01 };
  //rVecPt.resize( 150 );
  //rVecWt.resize( 50 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_74_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 13, 1, 1, 1, 1, 1, 0, 74 };
  double pa[] = { 0.5130671797338464e-03, -0.2958603896103896e-01, 0.1660406956574204e-01, 0.4803844614152614e+00, 0.2657620708215946e-01, 0.3207726489807764e+00, 0.1652217099371571e-01 };
  //rVecPt.resize( 222 );
  //rVecWt.resize( 74 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_86_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 15, 1, 1, 0, 2, 1, 0, 86 };
  double pa[] = { 0.1154401154401154e-01, 0.1194390908585628e-01, 0.3696028464541502e+00, 0.1111055571060340e-01, 0.6943540066026664e+00, 0.1187650129453714e-01, 0.3742430390903412e+00, 0.1181230374690448e-01 };
  //rVecPt.resize( 258 );
  //rVecWt.resize( 86 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_110_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 17, 1, 1, 0, 3, 1, 0, 110 };
  double pa[] = { 0.3828270494937162e-02, 0.9793737512487512e-02, 0.1851156353447362e+00, 0.8211737283191111e-02, 0.6904210483822922e+00, 0.9942814891178103e-02, 0.3956894730559419e+00, 0.9595471336070963e-02, 0.4783690288121502e+00, 0.9694996361663028e-02 };
  //rVecPt.resize( 330 );
  //rVecWt.resize( 110 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_146_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 19, 1, 1, 1, 3, 0, 1, 146 };
  double pa[] = { 0.5996313688621381e-03, 0.7210515360144488e-02, 0.7372999718620756e-02, 0.6764410400114264e+00, 0.7116355493117555e-02, 0.4174961227965453e+00, 0.6753829486314477e-02, 0.1574676672039082e+00, 0.7574394159054034e-02, 0.1403553811713183e+00, 0.4493328323269557e+00, 0.6991087353303262e-02 };
  //rVecPt.resize( 438 );
  //rVecWt.resize( 146 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_170_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 21, 1, 1, 1, 3, 1, 1, 170 };
  double pa[] = { 0.5544842902037365e-02, 0.6383674773515093e-02, 0.6071332770670752e-02, 0.2551252621114134e+00, 0.5183387587747790e-02, 0.6743601460362766e+00, 0.6317929009813725e-02, 0.4318910696719410e+00, 0.6201670006589077e-02, 0.2613931360335988e+00, 0.5477143385137348e-02, 0.4990453161796037e+00, 0.1446630744325115e+00, 0.5968383987681156e-02 };
  //rVecPt.resize( 510 );
  //rVecWt.resize( 170 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_194_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 23, 1, 1, 1, 4, 1, 1, 194 };
  double pa[] = { 0.1782340447244611e-02, 0.5573383178848738e-02, 0.5716905949977102e-02, 0.6712973442695226e+00, 0.5608704082587997e-02, 0.2892465627575439e+00, 0.5158237711805383e-02, 0.4446933178717437e+00, 0.5518771467273614e-02, 0.1299335447650067e+00, 0.4106777028169394e-02, 0.3457702197611283e+00, 0.5051846064614808e-02, 0.1590417105383530e+00, 0.8360360154824589e+00, 0.5530248916233094e-02 };
  //rVecPt.resize( 582 );
  //rVecWt.resize( 194 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_230_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 25, 1, 1, 0, 5, 2, 1, 230 };
  double pa[] = { -0.5522639919727325e-01, 0.4450274607445226e-02, 0.4492044687397611e+00, 0.4496841067921404e-02, 0.2520419490210201e+00, 0.5049153450478750e-02, 0.6981906658447242e+00, 0.3976408018051883e-02, 0.6587405243460960e+00, 0.4401400650381014e-02, 0.4038544050097660e-01, 0.1724544350544401e-01, 0.5823842309715585e+00, 0.4231083095357343e-02, 0.3545877390518688e+00, 0.5198069864064399e-02, 0.2272181808998187e+00, 0.4864661535886647e+00, 0.4695720972568883e-02 };
  //rVecPt.resize( 690 );
  //rVecWt.resize( 230 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_266_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 27, 1, 1, 1, 5, 1, 2, 266 };
  double pa[] = { -0.1313769127326952e-02, 0.4186853881700583e-02, -0.2522728704859336e-02, 0.7039373391585475e+00, 0.5315167977810885e-02, 0.1012526248572414e+00, 0.4047142377086219e-02, 0.4647448726420539e+00, 0.4112482394406990e-02, 0.3277420654971629e+00, 0.3595584899758782e-02, 0.6620338663699974e+00, 0.4256131351428158e-02, 0.8506508083520399e+00, 0.4229582700647240e-02, 0.3233484542692899e+00, 0.1153112011009701e+00, 0.4080914225780505e-02, 0.2314790158712601e+00, 0.5244939240922365e+00, 0.4071467593830964e-02 };
  //rVecPt.resize( 798 );
  //rVecWt.resize( 266 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_302_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 29, 1, 1, 0, 6, 2, 2, 302 };
  double pa[] = { 0.8545911725128148e-03, 0.3599119285025571e-02, 0.3515640345570105e+00, 0.3449788424305883e-02, 0.6566329410219612e+00, 0.3604822601419882e-02, 0.4729054132581005e+00, 0.3576729661743367e-02, 0.9618308522614784e-01, 0.2352101413689164e-02, 0.2219645236294178e+00, 0.3108953122413675e-02, 0.7011766416089545e+00, 0.3650045807677255e-02, 0.2644152887060663e+00, 0.2982344963171804e-02, 0.5718955891878961e+00, 0.3600820932216460e-02, 0.2510034751770465e+00, 0.8000727494073952e+00, 0.3571540554273387e-02, 0.1233548532583327e+00, 0.4127724083168531e+00, 0.3392312205006170e-02 };
  //rVecPt.resize( 906 );
  //rVecWt.resize( 302 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_350_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 31, 1, 1, 0, 6, 2, 3, 350 };
  double pa[] = { 0.3006796749453936e-02, 0.3050627745650771e-02, 0.7068965463912316e+00, 0.1621104600288991e-02, 0.4794682625712025e+00, 0.3005701484901752e-02, 0.1927533154878019e+00, 0.2990992529653774e-02, 0.6930357961327123e+00, 0.2982170644107595e-02, 0.3608302115520091e+00, 0.2721564237310992e-02, 0.6498486161496169e+00, 0.3033513795811141e-02, 0.1932945013230339e+00, 0.3007949555218533e-02, 0.3800494919899303e+00, 0.2881964603055307e-02, 0.2899558825499574e+00, 0.7934537856582316e+00, 0.2958357626535696e-02, 0.9684121455103957e-01, 0.8280801506686862e+00, 0.3036020026407088e-02, 0.1833434647041659e+00, 0.9074658265305127e+00, 0.2832187403926303e-02 };
  //rVecPt.resize( 1050 );
  //rVecWt.resize( 350 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_434_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 35, 1, 1, 1, 7, 2, 4, 434 };
  double pa[] = { 0.5265897968224436e-03, 0.2512317418927307e-02, 0.2548219972002607e-02, 0.6909346307509111e+00, 0.2530403801186355e-02, 0.1774836054609158e+00, 0.2014279020918528e-02, 0.4914342637784746e+00, 0.2501725168402936e-02, 0.6456664707424256e+00, 0.2513267174597564e-02, 0.2861289010307638e+00, 0.2302694782227416e-02, 0.7568084367178018e-01, 0.1462495621594614e-02, 0.3927259763368002e+00, 0.2445373437312980e-02, 0.8818132877794288e+00, 0.2417442375638981e-02, 0.9776428111182649e+00, 0.1910951282179532e-02, 0.2054823696403044e+00, 0.8689460322872412e+00, 0.2416930044324775e-02, 0.5905157048925271e+00, 0.7999278543857286e+00, 0.2512236854563495e-02, 0.5550152361076807e+00, 0.7717462626915901e+00, 0.2496644054553086e-02, 0.9371809858553722e+00, 0.3344363145343455e+00, 0.2236607760437849e-02 };
  //rVecPt.resize( 1302 );
  //rVecWt.resize( 434 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_590_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 41, 1, 1, 0, 9, 3, 6, 590 };
  double pa[] = { 0.3095121295306187e-03, 0.1852379698597489e-02, 0.7040954938227469e+00, 0.1871790639277744e-02, 0.6807744066455243e+00, 0.1858812585438317e-02, 0.6372546939258752e+00, 0.1852028828296213e-02, 0.5044419707800358e+00, 0.1846715956151242e-02, 0.4215761784010967e+00, 0.1818471778162769e-02, 0.3317920736472123e+00, 0.1749564657281154e-02, 0.2384736701421887e+00, 0.1617210647254411e-02, 0.1459036449157763e+00, 0.1384737234851692e-02, 0.6095034115507196e-01, 0.9764331165051050e-03, 0.6116843442009876e+00, 0.1857161196774078e-02, 0.3964755348199858e+00, 0.1705153996395864e-02, 0.1724782009907724e+00, 0.1300321685886048e-02, 0.5610263808622060e+00, 0.3518280927733519e+00, 0.1842866472905286e-02, 0.4742392842551980e+00, 0.2634716655937950e+00, 0.1802658934377451e-02, 0.5984126497885380e+00, 0.1816640840360209e+00, 0.1849830560443660e-02, 0.3791035407695563e+00, 0.1720795225656878e+00, 0.1713904507106709e-02, 0.2778673190586244e+00, 0.8213021581932511e-01, 0.1555213603396808e-02, 0.5033564271075117e+00, 0.8999205842074875e-01, 0.1802239128008525e-02 };
  //rVecPt.resize( 1770 );
  //rVecWt.resize( 590 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_770_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 47, 1, 1, 1, 10, 3, 9, 770 };
  double pa[] = { 0.2192942088181184e-03, 0.1421940344335877e-02, 0.1436433617319080e-02, 0.5087204410502360e-01, 0.6798123511050502e-03, 0.1228198790178831e+00, 0.9913184235294912e-03, 0.2026890814408786e+00, 0.1180207833238949e-02, 0.2847745156464294e+00, 0.1296599602080921e-02, 0.3656719078978026e+00, 0.1365871427428316e-02, 0.4428264886713469e+00, 0.1402988604775325e-02, 0.5140619627249735e+00, 0.1418645563595609e-02, 0.6306401219166803e+00, 0.1421376741851662e-02, 0.6716883332022612e+00, 0.1423996475490962e-02, 0.6979792685336881e+00, 0.1431554042178567e-02, 0.1446865674195309e+00, 0.9254401499865368e-03, 0.3390263475411216e+00, 0.1250239995053509e-02, 0.5335804651263506e+00, 0.1394365843329230e-02, 0.6944024393349413e-01, 0.2355187894242326e+00, 0.1127089094671749e-02, 0.2269004109529460e+00, 0.4102182474045730e+00, 0.1345753760910670e-02, 0.8025574607775339e-01, 0.6214302417481605e+00, 0.1424957283316783e-02, 0.1467999527896572e+00, 0.3245284345717394e+00, 0.1261523341237750e-02, 0.1571507769824727e+00, 0.5224482189696630e+00, 0.1392547106052696e-02, 0.2365702993157246e+00, 0.6017546634089558e+00, 0.1418761677877656e-02, 0.7714815866765732e-01, 0.4346575516141163e+00, 0.1338366684479554e-02, 0.3062936666210730e+00, 0.4908826589037616e+00, 0.1393700862676131e-02, 0.3822477379524787e+00, 0.5648768149099500e+00, 0.1415914757466932e-02 };
  //rVecPt.resize( 2310 );
  //rVecWt.resize( 770 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_974_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 53, 1, 1, 0, 12, 4, 12, 974 };
  double pa[] = { 0.1438294190527431e-03, 0.1125772288287004e-02, 0.4292963545341347e-01, 0.4948029341949241e-03, 0.1051426854086404e+00, 0.7357990109125470e-03, 0.1750024867623087e+00, 0.8889132771304384e-03, 0.2477653379650257e+00, 0.9888347838921435e-03, 0.3206567123955957e+00, 0.1053299681709471e-02, 0.3916520749849983e+00, 0.1092778807014578e-02, 0.4590825874187624e+00, 0.1114389394063227e-02, 0.5214563888415861e+00, 0.1123724788051555e-02, 0.6253170244654199e+00, 0.1125239325243814e-02, 0.6637926744523170e+00, 0.1126153271815905e-02, 0.6910410398498301e+00, 0.1130286931123841e-02, 0.7052907007457760e+00, 0.1134986534363955e-02, 0.1236686762657990e+00, 0.6823367927109931e-03, 0.2940777114468387e+00, 0.9454158160447096e-03, 0.4697753849207649e+00, 0.1074429975385679e-02, 0.6334563241139567e+00, 0.1129300086569132e-02, 0.5974048614181342e-01, 0.2029128752777523e+00, 0.8436884500901954e-03, 0.1375760408473636e+00, 0.4602621942484054e+00, 0.1075255720448885e-02, 0.3391016526336286e+00, 0.5030673999662036e+00, 0.1108577236864462e-02, 0.1271675191439820e+00, 0.2817606422442134e+00, 0.9566475323783357e-03, 0.2693120740413512e+00, 0.4331561291720157e+00, 0.1080663250717391e-02, 0.1419786452601918e+00, 0.6256167358580814e+00, 0.1126797131196295e-02, 0.6709284600738255e-01, 0.3798395216859157e+00, 0.1022568715358061e-02, 0.7057738183256172e-01, 0.5517505421423520e+00, 0.1108960267713108e-02, 0.2783888477882155e+00, 0.6029619156159187e+00, 0.1122790653435766e-02, 0.1979578938917407e+00, 0.3589606329589096e+00, 0.1032401847117460e-02, 0.2087307061103274e+00, 0.5348666438135476e+00, 0.1107249382283854e-02, 0.4055122137872836e+00, 0.5674997546074373e+00, 0.1121780048519972e-02 };
  //rVecPt.resize( 2922 );
  //rVecWt.resize( 974 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_1202_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 59, 1, 1, 1, 13, 4, 16, 1202 };
  double pa[] = { 0.1105189233267572e-03, 0.9133159786443561e-03, 0.9205232738090741e-03, 0.3712636449657089e-01, 0.3690421898017899e-03, 0.9140060412262223e-01, 0.5603990928680660e-03, 0.1531077852469906e+00, 0.6865297629282609e-03, 0.2180928891660612e+00, 0.7720338551145630e-03, 0.2839874532200175e+00, 0.8301545958894795e-03, 0.3491177600963764e+00, 0.8686692550179628e-03, 0.4121431461444309e+00, 0.8927076285846890e-03, 0.4718993627149127e+00, 0.9060820238568219e-03, 0.5273145452842337e+00, 0.9119777254940867e-03, 0.6209475332444019e+00, 0.9128720138604181e-03, 0.6569722711857291e+00, 0.9130714935691735e-03, 0.6841788309070143e+00, 0.9152873784554116e-03, 0.7012604330123631e+00, 0.9187436274321654e-03, 0.1072382215478166e+00, 0.5176977312965694e-03, 0.2582068959496968e+00, 0.7331143682101417e-03, 0.4172752955306717e+00, 0.8463232836379928e-03, 0.5700366911792503e+00, 0.9031122694253992e-03, 0.9827986018263947e+00, 0.1771774022615325e+00, 0.6485778453163257e-03, 0.9624249230326228e+00, 0.2475716463426288e+00, 0.7435030910982369e-03, 0.9402007994128811e+00, 0.3354616289066489e+00, 0.7998527891839054e-03, 0.9320822040143202e+00, 0.3173615246611977e+00, 0.8101731497468018e-03, 0.9043674199393299e+00, 0.4090268427085357e+00, 0.8483389574594331e-03, 0.8912407560074747e+00, 0.3854291150669224e+00, 0.8556299257311812e-03, 0.8676435628462708e+00, 0.4932221184851285e+00, 0.8803208679738260e-03, 0.8581979986041619e+00, 0.4785320675922435e+00, 0.8811048182425720e-03, 0.8396753624049856e+00, 0.4507422593157064e+00, 0.8850282341265444e-03, 0.8165288564022188e+00, 0.5632123020762100e+00, 0.9021342299040653e-03, 0.8015469370783529e+00, 0.5434303569693900e+00, 0.9010091677105086e-03, 0.7773563069070351e+00, 0.5123518486419871e+00, 0.9022692938426915e-03, 0.7661621213900394e+00, 0.6394279634749102e+00, 0.9158016174693465e-03, 0.7553584143533510e+00, 0.6269805509024392e+00, 0.9131578003189435e-03, 0.7344305757559503e+00, 0.6031161693096310e+00, 0.9107813579482705e-03, 0.7043837184021765e+00, 0.5693702498468441e+00, 0.9105760258970126e-03 };
  //rVecPt.resize( 3606 );
  //rVecWt.resize( 1202 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_1454_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 65, 1, 1, 0, 15, 5, 20, 1454 };
  double pa[] = { 0.7777160743261247e-04, 0.7557646413004701e-03, 0.3229290663413854e-01, 0.2841633806090617e-03, 0.8036733271462222e-01, 0.4374419127053555e-03, 0.1354289960531653e+00, 0.5417174740872172e-03, 0.1938963861114426e+00, 0.6148000891358593e-03, 0.2537343715011275e+00, 0.6664394485800705e-03, 0.3135251434752570e+00, 0.7025039356923220e-03, 0.3721558339375338e+00, 0.7268511789249627e-03, 0.4286809575195696e+00, 0.7422637534208629e-03, 0.4822510128282994e+00, 0.7509545035841214e-03, 0.5320679333566263e+00, 0.7548535057718401e-03, 0.6172998195394274e+00, 0.7554088969774001e-03, 0.6510679849127481e+00, 0.7553147174442808e-03, 0.6777315251687360e+00, 0.7564767653292297e-03, 0.6963109410648741e+00, 0.7587991808518730e-03, 0.7058935009831749e+00, 0.7608261832033027e-03, 0.9955546194091857e+00, 0.4021680447874916e-03, 0.9734115901794209e+00, 0.5804871793945964e-03, 0.9275693732388626e+00, 0.6792151955945159e-03, 0.8568022422795103e+00, 0.7336741211286294e-03, 0.7623495553719372e+00, 0.7581866300989608e-03, 0.5707522908892223e+00, 0.4387028039889501e+00, 0.7538257859800743e-03, 0.5196463388403083e+00, 0.3858908414762617e+00, 0.7483517247053123e-03, 0.4646337531215351e+00, 0.3301937372343854e+00, 0.7371763661112059e-03, 0.4063901697557691e+00, 0.2725423573563777e+00, 0.7183448895756934e-03, 0.3456329466643087e+00, 0.2139510237495250e+00, 0.6895815529822191e-03, 0.2831395121050332e+00, 0.1555922309786647e+00, 0.6480105801792886e-03, 0.2197682022925330e+00, 0.9892878979686097e-01, 0.5897558896594636e-03, 0.1564696098650355e+00, 0.4598642910675510e-01, 0.5095708849247346e-03, 0.6027356673721295e+00, 0.3376625140173426e+00, 0.7536906428909755e-03, 0.5496032320255096e+00, 0.2822301309727988e+00, 0.7472505965575118e-03, 0.4921707755234567e+00, 0.2248632342592540e+00, 0.7343017132279698e-03, 0.4309422998598483e+00, 0.1666224723456479e+00, 0.7130871582177445e-03, 0.3664108182313672e+00, 0.1086964901822169e+00, 0.6817022032112776e-03, 0.2990189057758436e+00, 0.5251989784120085e-01, 0.6380941145604121e-03, 0.6268724013144998e+00, 0.2297523657550023e+00, 0.7550381377920310e-03, 0.5707324144834607e+00, 0.1723080607093800e+00, 0.7478646640144802e-03, 0.5096360901960365e+00, 0.1140238465390513e+00, 0.7335918720601220e-03, 0.4438729938312456e+00, 0.5611522095882537e-01, 0.7110120527658118e-03, 0.6419978471082389e+00, 0.1164174423140873e+00, 0.7571363978689501e-03, 0.5817218061802611e+00, 0.5797589531445219e-01, 0.7489908329079234e-03 };
  //rVecPt.resize( 4362 );
  //rVecWt.resize( 1454 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_1730_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 71, 1, 1, 1, 16, 5, 25, 1730 };
  double pa[] = { 0.6309049437420976e-04, 0.6357185073530720e-03, 0.6398287705571748e-03, 0.2860923126194662e-01, 0.2221207162188168e-03, 0.7142556767711522e-01, 0.3475784022286848e-03, 0.1209199540995559e+00, 0.4350742443589804e-03, 0.1738673106594379e+00, 0.4978569136522127e-03, 0.2284645438467734e+00, 0.5435036221998053e-03, 0.2834807671701512e+00, 0.5765913388219542e-03, 0.3379680145467339e+00, 0.6001200359226003e-03, 0.3911355454819537e+00, 0.6162178172717512e-03, 0.4422860353001403e+00, 0.6265218152438485e-03, 0.4907781568726057e+00, 0.6323987160974212e-03, 0.5360006153211468e+00, 0.6350767851540569e-03, 0.6142105973596603e+00, 0.6354362775297107e-03, 0.6459300387977504e+00, 0.6352302462706235e-03, 0.6718056125089225e+00, 0.6358117881417972e-03, 0.6910888533186254e+00, 0.6373101590310117e-03, 0.7030467416823252e+00, 0.6390428961368665e-03, 0.8354951166354646e-01, 0.3186913449946576e-03, 0.2050143009099486e+00, 0.4678028558591711e-03, 0.3370208290706637e+00, 0.5538829697598626e-03, 0.4689051484233963e+00, 0.6044475907190476e-03, 0.5939400424557334e+00, 0.6313575103509012e-03, 0.1394983311832261e+00, 0.4097581162050343e-01, 0.4078626431855630e-03, 0.1967999180485014e+00, 0.8851987391293348e-01, 0.4759933057812725e-03, 0.2546183732548967e+00, 0.1397680182969819e+00, 0.5268151186413440e-03, 0.3121281074713875e+00, 0.1929452542226526e+00, 0.5643048560507316e-03, 0.3685981078502492e+00, 0.2467898337061562e+00, 0.5914501076613073e-03, 0.4233760321547856e+00, 0.3003104124785409e+00, 0.6104561257874195e-03, 0.4758671236059246e+00, 0.3526684328175033e+00, 0.6230252860707806e-03, 0.5255178579796463e+00, 0.4031134861145713e+00, 0.6305618761760796e-03, 0.5718025633734589e+00, 0.4509426448342351e+00, 0.6343092767597889e-03, 0.2686927772723415e+00, 0.4711322502423248e-01, 0.5176268945737826e-03, 0.3306006819904809e+00, 0.9784487303942695e-01, 0.5564840313313692e-03, 0.3904906850594983e+00, 0.1505395810025273e+00, 0.5856426671038980e-03, 0.4479957951904390e+00, 0.2039728156296050e+00, 0.6066386925777091e-03, 0.5027076848919780e+00, 0.2571529941121107e+00, 0.6208824962234458e-03, 0.5542087392260217e+00, 0.3092191375815670e+00, 0.6296314297822907e-03, 0.6020850887375187e+00, 0.3593807506130276e+00, 0.6340423756791859e-03, 0.4019851409179594e+00, 0.5063389934378671e-01, 0.5829627677107342e-03, 0.4635614567449800e+00, 0.1032422269160612e+00, 0.6048693376081110e-03, 0.5215860931591575e+00, 0.1566322094006254e+00, 0.6202362317732461e-03, 0.5758202499099271e+00, 0.2098082827491099e+00, 0.6299005328403779e-03, 0.6259893683876795e+00, 0.2618824114553391e+00, 0.6347722390609353e-03, 0.5313795124811891e+00, 0.5263245019338556e-01, 0.6203778981238834e-03, 0.5893317955931995e+00, 0.1061059730982005e+00, 0.6308414671239979e-03, 0.6426246321215801e+00, 0.1594171564034221e+00, 0.6362706466959498e-03, 0.6511904367376113e+00, 0.5354789536565540e-01, 0.6375414170333233e-03 };
  //rVecPt.resize( 5190 );
  //rVecWt.resize( 1730 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_2030_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 77, 1, 1, 0, 18, 6, 30, 2030 };
  double pa[] = { 0.4656031899197431e-04, 0.5421549195295507e-03, 0.2540835336814348e-01, 0.1778522133346553e-03, 0.6399322800504915e-01, 0.2811325405682796e-03, 0.1088269469804125e+00, 0.3548896312631459e-03, 0.1570670798818287e+00, 0.4090310897173364e-03, 0.2071163932282514e+00, 0.4493286134169965e-03, 0.2578914044450844e+00, 0.4793728447962723e-03, 0.3085687558169623e+00, 0.5015415319164265e-03, 0.3584719706267024e+00, 0.5175127372677937e-03, 0.4070135594428709e+00, 0.5285522262081019e-03, 0.4536618626222638e+00, 0.5356832703713962e-03, 0.4979195686463577e+00, 0.5397914736175170e-03, 0.5393075111126999e+00, 0.5416899441599930e-03, 0.6115617676843916e+00, 0.5419308476889938e-03, 0.6414308435160159e+00, 0.5416936902030596e-03, 0.6664099412721607e+00, 0.5419544338703164e-03, 0.6859161771214913e+00, 0.5428983656630975e-03, 0.6993625593503890e+00, 0.5442286500098193e-03, 0.7062393387719380e+00, 0.5452250345057301e-03, 0.7479028168349763e-01, 0.2568002497728530e-03, 0.1848951153969366e+00, 0.3827211700292145e-03, 0.3059529066581305e+00, 0.4579491561917824e-03, 0.4285556101021362e+00, 0.5042003969083574e-03, 0.5468758653496526e+00, 0.5312708889976025e-03, 0.6565821978343439e+00, 0.5438401790747117e-03, 0.1253901572367117e+00, 0.3681917226439641e-01, 0.3316041873197344e-03, 0.1775721510383941e+00, 0.7982487607213301e-01, 0.3899113567153771e-03, 0.2305693358216114e+00, 0.1264640966592335e+00, 0.4343343327201309e-03, 0.2836502845992063e+00, 0.1751585683418957e+00, 0.4679415262318919e-03, 0.3361794746232590e+00, 0.2247995907632670e+00, 0.4930847981631031e-03, 0.3875979172264824e+00, 0.2745299257422246e+00, 0.5115031867540091e-03, 0.4374019316999074e+00, 0.3236373482441118e+00, 0.5245217148457367e-03, 0.4851275843340022e+00, 0.3714967859436741e+00, 0.5332041499895321e-03, 0.5303391803806868e+00, 0.4175353646321745e+00, 0.5384583126021542e-03, 0.5726197380596287e+00, 0.4612084406355461e+00, 0.5411067210798852e-03, 0.2431520732564863e+00, 0.4258040133043952e-01, 0.4259797391468714e-03, 0.3002096800895869e+00, 0.8869424306722721e-01, 0.4604931368460021e-03, 0.3558554457457432e+00, 0.1368811706510655e+00, 0.4871814878255202e-03, 0.4097782537048887e+00, 0.1860739985015033e+00, 0.5072242910074885e-03, 0.4616337666067458e+00, 0.2354235077395853e+00, 0.5217069845235350e-03, 0.5110707008417874e+00, 0.2842074921347011e+00, 0.5315785966280310e-03, 0.5577415286163795e+00, 0.3317784414984102e+00, 0.5376833708758905e-03, 0.6013060431366950e+00, 0.3775299002040700e+00, 0.5408032092069521e-03, 0.3661596767261781e+00, 0.4599367887164592e-01, 0.4842744917904866e-03, 0.4237633153506581e+00, 0.9404893773654421e-01, 0.5048926076188130e-03, 0.4786328454658452e+00, 0.1431377109091971e+00, 0.5202607980478373e-03, 0.5305702076789774e+00, 0.1924186388843570e+00, 0.5309932388325743e-03, 0.5793436224231788e+00, 0.2411590944775190e+00, 0.5377419770895208e-03, 0.6247069017094747e+00, 0.2886871491583605e+00, 0.5411696331677717e-03, 0.4874315552535204e+00, 0.4804978774953206e-01, 0.5197996293282420e-03, 0.5427337322059053e+00, 0.9716857199366665e-01, 0.5311120836622945e-03, 0.5943493747246700e+00, 0.1465205839795055e+00, 0.5384309319956951e-03, 0.6421314033564943e+00, 0.1953579449803574e+00, 0.5421859504051886e-03, 0.6020628374713980e+00, 0.4916375015738108e-01, 0.5390948355046314e-03, 0.6529222529856881e+00, 0.9861621540127005e-01, 0.5433312705027845e-03 };
  //rVecPt.resize( 6090 );
  //rVecWt.resize( 2030 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_2354_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 83, 1, 1, 1, 19, 6, 36, 2354 };
  double pa[] = { 0.3922616270665292e-04, 0.4678202801282136e-03, 0.4703831750854424e-03, 0.2290024646530589e-01, 0.1437832228979900e-03, 0.5779086652271284e-01, 0.2303572493577644e-03, 0.9863103576375984e-01, 0.2933110752447454e-03, 0.1428155792982185e+00, 0.3402905998359838e-03, 0.1888978116601463e+00, 0.3759138466870372e-03, 0.2359091682970210e+00, 0.4030638447899798e-03, 0.2831228833706171e+00, 0.4236591432242211e-03, 0.3299495857966693e+00, 0.4390522656946746e-03, 0.3758840802660796e+00, 0.4502523466626247e-03, 0.4204751831009480e+00, 0.4580577727783541e-03, 0.4633068518751051e+00, 0.4631391616615899e-03, 0.5039849474507313e+00, 0.4660928953698676e-03, 0.5421265793440747e+00, 0.4674751807936953e-03, 0.6092660230557310e+00, 0.4676414903932920e-03, 0.6374654204984869e+00, 0.4674086492347870e-03, 0.6615136472609892e+00, 0.4674928539483207e-03, 0.6809487285958127e+00, 0.4680748979686447e-03, 0.6952980021665196e+00, 0.4690449806389040e-03, 0.7041245497695400e+00, 0.4699877075860818e-03, 0.6744033088306065e-01, 0.2099942281069176e-03, 0.1678684485334166e+00, 0.3172269150712804e-03, 0.2793559049539613e+00, 0.3832051358546523e-03, 0.3935264218057639e+00, 0.4252193818146985e-03, 0.5052629268232558e+00, 0.4513807963755000e-03, 0.6107905315437531e+00, 0.4657797469114178e-03, 0.1135081039843524e+00, 0.3331954884662588e-01, 0.2733362800522836e-03, 0.1612866626099378e+00, 0.7247167465436538e-01, 0.3235485368463559e-03, 0.2100786550168205e+00, 0.1151539110849745e+00, 0.3624908726013453e-03, 0.2592282009459942e+00, 0.1599491097143677e+00, 0.3925540070712828e-03, 0.3081740561320203e+00, 0.2058699956028027e+00, 0.4156129781116235e-03, 0.3564289781578164e+00, 0.2521624953502911e+00, 0.4330644984623263e-03, 0.4035587288240703e+00, 0.2982090785797674e+00, 0.4459677725921312e-03, 0.4491671196373903e+00, 0.3434762087235733e+00, 0.4551593004456795e-03, 0.4928854782917489e+00, 0.3874831357203437e+00, 0.4613341462749918e-03, 0.5343646791958988e+00, 0.4297814821746926e+00, 0.4651019618269806e-03, 0.5732683216530990e+00, 0.4699402260943537e+00, 0.4670249536100625e-03, 0.2214131583218986e+00, 0.3873602040643895e-01, 0.3549555576441708e-03, 0.2741796504750071e+00, 0.8089496256902013e-01, 0.3856108245249010e-03, 0.3259797439149485e+00, 0.1251732177620872e+00, 0.4098622845756882e-03, 0.3765441148826891e+00, 0.1706260286403185e+00, 0.4286328604268950e-03, 0.4255773574530558e+00, 0.2165115147300408e+00, 0.4427802198993945e-03, 0.4727795117058430e+00, 0.2622089812225259e+00, 0.4530473511488561e-03, 0.5178546895819012e+00, 0.3071721431296201e+00, 0.4600805475703138e-03, 0.5605141192097460e+00, 0.3508998998801138e+00, 0.4644599059958017e-03, 0.6004763319352512e+00, 0.3929160876166931e+00, 0.4667274455712508e-03, 0.3352842634946949e+00, 0.4202563457288019e-01, 0.4069360518020356e-03, 0.3891971629814670e+00, 0.8614309758870850e-01, 0.4260442819919195e-03, 0.4409875565542281e+00, 0.1314500879380001e+00, 0.4408678508029063e-03, 0.4904893058592484e+00, 0.1772189657383859e+00, 0.4518748115548597e-03, 0.5375056138769549e+00, 0.2228277110050294e+00, 0.4595564875375116e-03, 0.5818255708669969e+00, 0.2677179935014386e+00, 0.4643988774315846e-03, 0.6232334858144959e+00, 0.3113675035544165e+00, 0.4668827491646946e-03, 0.4489485354492058e+00, 0.4409162378368174e-01, 0.4400541823741973e-03, 0.5015136875933150e+00, 0.8939009917748489e-01, 0.4514512890193797e-03, 0.5511300550512623e+00, 0.1351806029383365e+00, 0.4596198627347549e-03, 0.5976720409858000e+00, 0.1808370355053196e+00, 0.4648659016801781e-03, 0.6409956378989354e+00, 0.2257852192301602e+00, 0.4675502017157673e-03, 0.5581222330827514e+00, 0.4532173421637160e-01, 0.4598494476455523e-03, 0.6074705984161695e+00, 0.9117488031840314e-01, 0.4654916955152048e-03, 0.6532272537379033e+00, 0.1369294213140155e+00, 0.4684709779505137e-03, 0.6594761494500487e+00, 0.4589901487275583e-01, 0.4691445539106986e-03 };
  //rVecPt.resize( 7062 );
  //rVecWt.resize( 2354 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_2702_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 89, 1, 1, 0, 21, 7, 42, 2702 };
  double pa[] = { 0.2998675149888161e-04, 0.4077860529495355e-03, 0.2065562538818703e-01, 0.1185349192520667e-03, 0.5250918173022379e-01, 0.1913408643425751e-03, 0.8993480082038376e-01, 0.2452886577209897e-03, 0.1306023924436019e+00, 0.2862408183288702e-03, 0.1732060388531418e+00, 0.3178032258257357e-03, 0.2168727084820249e+00, 0.3422945667633690e-03, 0.2609528309173586e+00, 0.3612790520235922e-03, 0.3049252927938952e+00, 0.3758638229818521e-03, 0.3483484138084404e+00, 0.3868711798859953e-03, 0.3908321549106406e+00, 0.3949429933189938e-03, 0.4320210071894814e+00, 0.4006068107541156e-03, 0.4715824795890053e+00, 0.4043192149672723e-03, 0.5091984794078453e+00, 0.4064947495808078e-03, 0.5445580145650803e+00, 0.4075245619813152e-03, 0.6072575796841768e+00, 0.4076423540893566e-03, 0.6339484505755803e+00, 0.4074280862251555e-03, 0.6570718257486958e+00, 0.4074163756012244e-03, 0.6762557330090709e+00, 0.4077647795071246e-03, 0.6911161696923790e+00, 0.4084517552782530e-03, 0.7012841911659961e+00, 0.4092468459224052e-03, 0.7064559272410020e+00, 0.4097872687240906e-03, 0.6123554989894765e-01, 0.1738986811745028e-03, 0.1533070348312393e+00, 0.2659616045280191e-03, 0.2563902605244206e+00, 0.3240596008171533e-03, 0.3629346991663361e+00, 0.3621195964432943e-03, 0.4683949968987538e+00, 0.3868838330760539e-03, 0.5694479240657952e+00, 0.4018911532693111e-03, 0.6634465430993955e+00, 0.4089929432983252e-03, 0.1033958573552305e+00, 0.3034544009063584e-01, 0.2279907527706409e-03, 0.1473521412414395e+00, 0.6618803044247135e-01, 0.2715205490578897e-03, 0.1924552158705967e+00, 0.1054431128987715e+00, 0.3057917896703976e-03, 0.2381094362890328e+00, 0.1468263551238858e+00, 0.3326913052452555e-03, 0.2838121707936760e+00, 0.1894486108187886e+00, 0.3537334711890037e-03, 0.3291323133373415e+00, 0.2326374238761579e+00, 0.3700567500783129e-03, 0.3736896978741460e+00, 0.2758485808485768e+00, 0.3825245372589122e-03, 0.4171406040760013e+00, 0.3186179331996921e+00, 0.3918125171518296e-03, 0.4591677985256915e+00, 0.3605329796303794e+00, 0.3984720419937579e-03, 0.4994733831718418e+00, 0.4012147253586509e+00, 0.4029746003338211e-03, 0.5377731830445096e+00, 0.4403050025570692e+00, 0.4057428632156627e-03, 0.5737917830001331e+00, 0.4774565904277483e+00, 0.4071719274114857e-03, 0.2027323586271389e+00, 0.3544122504976147e-01, 0.2990236950664119e-03, 0.2516942375187273e+00, 0.7418304388646328e-01, 0.3262951734212878e-03, 0.3000227995257181e+00, 0.1150502745727186e+00, 0.3482634608242413e-03, 0.3474806691046342e+00, 0.1571963371209364e+00, 0.3656596681700892e-03, 0.3938103180359209e+00, 0.1999631877247100e+00, 0.3791740467794218e-03, 0.4387519590455703e+00, 0.2428073457846535e+00, 0.3894034450156905e-03, 0.4820503960077787e+00, 0.2852575132906155e+00, 0.3968600245508371e-03, 0.5234573778475101e+00, 0.3268884208674639e+00, 0.4019931351420050e-03, 0.5627318647235282e+00, 0.3673033321675939e+00, 0.4052108801278599e-03, 0.5996390607156954e+00, 0.4061211551830290e+00, 0.4068978613940934e-03, 0.3084780753791947e+00, 0.3860125523100059e-01, 0.3454275351319704e-03, 0.3589988275920223e+00, 0.7928938987104867e-01, 0.3629963537007920e-03, 0.4078628415881973e+00, 0.1212614643030087e+00, 0.3770187233889873e-03, 0.4549287258889735e+00, 0.1638770827382693e+00, 0.3878608613694378e-03, 0.5000278512957279e+00, 0.2065965798260176e+00, 0.3959065270221274e-03, 0.5429785044928199e+00, 0.2489436378852235e+00, 0.4015286975463570e-03, 0.5835939850491711e+00, 0.2904811368946891e+00, 0.4050866785614717e-03, 0.6216870353444856e+00, 0.3307941957666609e+00, 0.4069320185051913e-03, 0.4151104662709091e+00, 0.4064829146052554e-01, 0.3760120964062763e-03, 0.4649804275009218e+00, 0.8258424547294755e-01, 0.3870969564418064e-03, 0.5124695757009662e+00, 0.1251841962027289e+00, 0.3955287790534055e-03, 0.5574711100606224e+00, 0.1679107505976331e+00, 0.4015361911302668e-03, 0.5998597333287227e+00, 0.2102805057358715e+00, 0.4053836986719548e-03, 0.6395007148516600e+00, 0.2518418087774107e+00, 0.4073578673299117e-03, 0.5188456224746252e+00, 0.4194321676077518e-01, 0.3954628379231406e-03, 0.5664190707942778e+00, 0.8457661551921499e-01, 0.4017645508847530e-03, 0.6110464353283153e+00, 0.1273652932519396e+00, 0.4059030348651293e-03, 0.6526430302051563e+00, 0.1698173239076354e+00, 0.4080565809484880e-03, 0.6167551880377548e+00, 0.4266398851548864e-01, 0.4063018753664651e-03, 0.6607195418355383e+00, 0.8551925814238349e-01, 0.4087191292799671e-03 };
  //rVecPt.resize( 8106 );
  //rVecWt.resize( 2702 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_3074_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 95, 1, 1, 1, 22, 7, 49, 3074 };
  double pa[] = { 0.2599095953754734e-04, 0.3586067974412447e-03, 0.3603134089687541e-03, 0.1886108518723392e-01, 0.9831528474385880e-04, 0.4800217244625303e-01, 0.1605023107954450e-03, 0.8244922058397242e-01, 0.2072200131464099e-03, 0.1200408362484023e+00, 0.2431297618814187e-03, 0.1595773530809965e+00, 0.2711819064496707e-03, 0.2002635973434064e+00, 0.2932762038321116e-03, 0.2415127590139982e+00, 0.3107032514197368e-03, 0.2828584158458477e+00, 0.3243808058921213e-03, 0.3239091015338138e+00, 0.3349899091374030e-03, 0.3643225097962194e+00, 0.3430580688505218e-03, 0.4037897083691802e+00, 0.3490124109290343e-03, 0.4420247515194127e+00, 0.3532148948561955e-03, 0.4787572538464938e+00, 0.3559862669062833e-03, 0.5137265251275234e+00, 0.3576224317551411e-03, 0.5466764056654611e+00, 0.3584050533086076e-03, 0.6054859420813535e+00, 0.3584903581373224e-03, 0.6308106701764562e+00, 0.3582991879040586e-03, 0.6530369230179584e+00, 0.3582371187963125e-03, 0.6718609524611158e+00, 0.3584353631122350e-03, 0.6869676499894013e+00, 0.3589120166517785e-03, 0.6980467077240748e+00, 0.3595445704531601e-03, 0.7048241721250522e+00, 0.3600943557111074e-03, 0.5591105222058232e-01, 0.1456447096742039e-03, 0.1407384078513916e+00, 0.2252370188283782e-03, 0.2364035438976309e+00, 0.2766135443474897e-03, 0.3360602737818170e+00, 0.3110729491500851e-03, 0.4356292630054665e+00, 0.3342506712303391e-03, 0.5321569415256174e+00, 0.3491981834026860e-03, 0.6232956305040554e+00, 0.3576003604348932e-03, 0.9469870086838469e-01, 0.2778748387309470e-01, 0.1921921305788564e-03, 0.1353170300568141e+00, 0.6076569878628364e-01, 0.2301458216495632e-03, 0.1771679481726077e+00, 0.9703072762711040e-01, 0.2604248549522893e-03, 0.2197066664231751e+00, 0.1354112458524762e+00, 0.2845275425870697e-03, 0.2624783557374927e+00, 0.1750996479744100e+00, 0.3036870897974840e-03, 0.3050969521214442e+00, 0.2154896907449802e+00, 0.3188414832298066e-03, 0.3472252637196021e+00, 0.2560954625740152e+00, 0.3307046414722089e-03, 0.3885610219026360e+00, 0.2965070050624096e+00, 0.3398330969031360e-03, 0.4288273776062765e+00, 0.3363641488734497e+00, 0.3466757899705373e-03, 0.4677662471302948e+00, 0.3753400029836788e+00, 0.3516095923230054e-03, 0.5051333589553359e+00, 0.4131297522144286e+00, 0.3549645184048486e-03, 0.5406942145810492e+00, 0.4494423776081795e+00, 0.3570415969441392e-03, 0.5742204122576457e+00, 0.4839938958841502e+00, 0.3581251798496118e-03, 0.1865407027225188e+00, 0.3259144851070796e-01, 0.2543491329913348e-03, 0.2321186453689432e+00, 0.6835679505297343e-01, 0.2786711051330776e-03, 0.2773159142523882e+00, 0.1062284864451989e+00, 0.2985552361083679e-03, 0.3219200192237254e+00, 0.1454404409323047e+00, 0.3145867929154039e-03, 0.3657032593944029e+00, 0.1854018282582510e+00, 0.3273290662067609e-03, 0.4084376778363622e+00, 0.2256297412014750e+00, 0.3372705511943501e-03, 0.4499004945751427e+00, 0.2657104425000896e+00, 0.3448274437851510e-03, 0.4898758141326335e+00, 0.3052755487631557e+00, 0.3503592783048583e-03, 0.5281547442266309e+00, 0.3439863920645423e+00, 0.3541854792663162e-03, 0.5645346989813992e+00, 0.3815229456121914e+00, 0.3565995517909428e-03, 0.5988181252159848e+00, 0.4175752420966734e+00, 0.3578802078302898e-03, 0.2850425424471603e+00, 0.3562149509862536e-01, 0.2958644592860982e-03, 0.3324619433027876e+00, 0.7330318886871096e-01, 0.3119548129116835e-03, 0.3785848333076282e+00, 0.1123226296008472e+00, 0.3250745225005984e-03, 0.4232891028562115e+00, 0.1521084193337708e+00, 0.3355153415935208e-03, 0.4664287050829722e+00, 0.1921844459223610e+00, 0.3435847568549328e-03, 0.5078458493735726e+00, 0.2321360989678303e+00, 0.3495786831622488e-03, 0.5473779816204180e+00, 0.2715886486360520e+00, 0.3537767805534621e-03, 0.5848617133811376e+00, 0.3101924707571355e+00, 0.3564459815421428e-03, 0.6201348281584888e+00, 0.3476121052890973e+00, 0.3578464061225468e-03, 0.3852191185387871e+00, 0.3763224880035108e-01, 0.3239748762836212e-03, 0.4325025061073423e+00, 0.7659581935637135e-01, 0.3345491784174287e-03, 0.4778486229734490e+00, 0.1163381306083900e+00, 0.3429126177301782e-03, 0.5211663693009000e+00, 0.1563890598752899e+00, 0.3492420343097421e-03, 0.5623469504853703e+00, 0.1963320810149200e+00, 0.3537399050235257e-03, 0.6012718188659246e+00, 0.2357847407258738e+00, 0.3566209152659172e-03, 0.6378179206390117e+00, 0.2743846121244060e+00, 0.3581084321919782e-03, 0.4836936460214534e+00, 0.3895902610739024e-01, 0.3426522117591512e-03, 0.5293792562683797e+00, 0.7871246819312640e-01, 0.3491848770121379e-03, 0.5726281253100033e+00, 0.1187963808202981e+00, 0.3539318235231476e-03, 0.6133658776169068e+00, 0.1587914708061787e+00, 0.3570231438458694e-03, 0.6515085491865307e+00, 0.1983058575227646e+00, 0.3586207335051714e-03, 0.5778692716064976e+00, 0.3977209689791542e-01, 0.3541196205164025e-03, 0.6207904288086192e+00, 0.7990157592981152e-01, 0.3574296911573953e-03, 0.6608688171046802e+00, 0.1199671308754309e+00, 0.3591993279818963e-03, 0.6656263089489130e+00, 0.4015955957805969e-01, 0.3595855034661997e-03 };
  //rVecPt.resize( 9222 );
  //rVecWt.resize( 3074 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_3470_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 101, 1, 1, 0, 24, 8, 56, 3470 };
  double pa[] = { 0.2040382730826330e-04, 0.3178149703889544e-03, 0.1721420832906233e-01, 0.8288115128076110e-04, 0.4408875374981770e-01, 0.1360883192522954e-03, 0.7594680813878681e-01, 0.1766854454542662e-03, 0.1108335359204799e+00, 0.2083153161230153e-03, 0.1476517054388567e+00, 0.2333279544657158e-03, 0.1856731870860615e+00, 0.2532809539930247e-03, 0.2243634099428821e+00, 0.2692472184211158e-03, 0.2633006881662727e+00, 0.2819949946811885e-03, 0.3021340904916283e+00, 0.2920953593973030e-03, 0.3405594048030089e+00, 0.2999889782948352e-03, 0.3783044434007372e+00, 0.3060292120496902e-03, 0.4151194767407910e+00, 0.3105109167522192e-03, 0.4507705766443257e+00, 0.3136902387550312e-03, 0.4850346056573187e+00, 0.3157984652454632e-03, 0.5176950817792470e+00, 0.3170516518425422e-03, 0.5485384240820989e+00, 0.3176568425633755e-03, 0.6039117238943308e+00, 0.3177198411207062e-03, 0.6279956655573113e+00, 0.3175519492394733e-03, 0.6493636169568952e+00, 0.3174654952634756e-03, 0.6677644117704504e+00, 0.3175676415467654e-03, 0.6829368572115624e+00, 0.3178923417835410e-03, 0.6946195818184121e+00, 0.3183788287531909e-03, 0.7025711542057026e+00, 0.3188755151918807e-03, 0.7066004767140119e+00, 0.3191916889313849e-03, 0.5132537689946062e-01, 0.1231779611744508e-03, 0.1297994661331225e+00, 0.1924661373839880e-03, 0.2188852049401307e+00, 0.2380881867403424e-03, 0.3123174824903457e+00, 0.2693100663037885e-03, 0.4064037620738195e+00, 0.2908673382834366e-03, 0.4984958396944782e+00, 0.3053914619381535e-03, 0.5864975046021365e+00, 0.3143916684147777e-03, 0.6686711634580175e+00, 0.3187042244055363e-03, 0.8715738780835950e-01, 0.2557175233367578e-01, 0.1635219535869790e-03, 0.1248383123134007e+00, 0.5604823383376681e-01, 0.1968109917696070e-03, 0.1638062693383378e+00, 0.8968568601900765e-01, 0.2236754342249974e-03, 0.2035586203373176e+00, 0.1254086651976279e+00, 0.2453186687017181e-03, 0.2436798975293774e+00, 0.1624780150162012e+00, 0.2627551791580541e-03, 0.2838207507773806e+00, 0.2003422342683208e+00, 0.2767654860152220e-03, 0.3236787502217692e+00, 0.2385628026255263e+00, 0.2879467027765895e-03, 0.3629849554840691e+00, 0.2767731148783578e+00, 0.2967639918918702e-03, 0.4014948081992087e+00, 0.3146542308245309e+00, 0.3035900684660351e-03, 0.4389818379260225e+00, 0.3519196415895088e+00, 0.3087338237298308e-03, 0.4752331143674377e+00, 0.3883050984023654e+00, 0.3124608838860167e-03, 0.5100457318374018e+00, 0.4235613423908649e+00, 0.3150084294226743e-03, 0.5432238388954868e+00, 0.4574484717196220e+00, 0.3165958398598402e-03, 0.5745758685072442e+00, 0.4897311639255524e+00, 0.3174320440957372e-03, 0.1723981437592809e+00, 0.3010630597881105e-01, 0.2182188909812599e-03, 0.2149553257844597e+00, 0.6326031554204694e-01, 0.2399727933921445e-03, 0.2573256081247422e+00, 0.9848566980258631e-01, 0.2579796133514652e-03, 0.2993163751238106e+00, 0.1350835952384266e+00, 0.2727114052623535e-03, 0.3407238005148000e+00, 0.1725184055442181e+00, 0.2846327656281355e-03, 0.3813454978483264e+00, 0.2103559279730725e+00, 0.2941491102051334e-03, 0.4209848104423343e+00, 0.2482278774554860e+00, 0.3016049492136107e-03, 0.4594519699996300e+00, 0.2858099509982883e+00, 0.3072949726175648e-03, 0.4965640166185930e+00, 0.3228075659915428e+00, 0.3114768142886460e-03, 0.5321441655571562e+00, 0.3589459907204151e+00, 0.3143823673666223e-03, 0.5660208438582166e+00, 0.3939630088864310e+00, 0.3162269764661535e-03, 0.5980264315964364e+00, 0.4276029922949089e+00, 0.3172164663759821e-03, 0.2644215852350733e+00, 0.3300939429072552e-01, 0.2554575398967435e-03, 0.3090113743443063e+00, 0.6803887650078501e-01, 0.2701704069135677e-03, 0.3525871079197808e+00, 0.1044326136206709e+00, 0.2823693413468940e-03, 0.3950418005354029e+00, 0.1416751597517679e+00, 0.2922898463214289e-03, 0.4362475663430163e+00, 0.1793408610504821e+00, 0.3001829062162428e-03, 0.4760661812145854e+00, 0.2170630750175722e+00, 0.3062890864542953e-03, 0.5143551042512103e+00, 0.2545145157815807e+00, 0.3108328279264746e-03, 0.5509709026935597e+00, 0.2913940101706601e+00, 0.3140243146201245e-03, 0.5857711030329428e+00, 0.3274169910910705e+00, 0.3160638030977130e-03, 0.6186149917404392e+00, 0.3623081329317265e+00, 0.3171462882206275e-03, 0.3586894569557064e+00, 0.3497354386450040e-01, 0.2812388416031796e-03, 0.4035266610019441e+00, 0.7129736739757095e-01, 0.2912137500288045e-03, 0.4467775312332510e+00, 0.1084758620193165e+00, 0.2993241256502206e-03, 0.4883638346608543e+00, 0.1460915689241772e+00, 0.3057101738983822e-03, 0.5281908348434601e+00, 0.1837790832369980e+00, 0.3105319326251432e-03, 0.5661542687149311e+00, 0.2212075390874021e+00, 0.3139565514428167e-03, 0.6021450102031452e+00, 0.2580682841160985e+00, 0.3161543006806366e-03, 0.6360520783610050e+00, 0.2940656362094121e+00, 0.3172985960613294e-03, 0.4521611065087196e+00, 0.3631055365867002e-01, 0.2989400336901431e-03, 0.4959365651560963e+00, 0.7348318468484350e-01, 0.3054555883947677e-03, 0.5376815804038283e+00, 0.1111087643812648e+00, 0.3104764960807702e-03, 0.5773314480243768e+00, 0.1488226085145408e+00, 0.3141015825977616e-03, 0.6148113245575056e+00, 0.1862892274135151e+00, 0.3164520621159896e-03, 0.6500407462842380e+00, 0.2231909701714456e+00, 0.3176652305912204e-03, 0.5425151448707213e+00, 0.3718201306118944e-01, 0.3105097161023939e-03, 0.5841860556907931e+00, 0.7483616335067346e-01, 0.3143014117890550e-03, 0.6234632186851500e+00, 0.1125990834266120e+00, 0.3168172866287200e-03, 0.6602934551848843e+00, 0.1501303813157619e+00, 0.3181401865570968e-03, 0.6278573968375105e+00, 0.3767559930245720e-01, 0.3170663659156037e-03, 0.6665611711264577e+00, 0.7548443301360158e-01, 0.3185447944625510e-03 };
  //rVecPt.resize( 10410 );
  //rVecWt.resize( 3470 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_3890_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 107, 1, 1, 1, 25, 8, 64, 3890 };
  double pa[] = { 0.1807395252196920e-04, 0.2836065837530581e-03, 0.2848008782238827e-03, 0.1587876419858352e-01, 0.7013149266673816e-04, 0.4069193593751206e-01, 0.1162798021956766e-03, 0.7025888115257997e-01, 0.1518728583972105e-03, 0.1027495450028704e+00, 0.1798796108216934e-03, 0.1371457730893426e+00, 0.2022593385972785e-03, 0.1727758532671953e+00, 0.2203093105575464e-03, 0.2091492038929037e+00, 0.2349294234299855e-03, 0.2458813281751915e+00, 0.2467682058747003e-03, 0.2826545859450066e+00, 0.2563092683572224e-03, 0.3191957291799622e+00, 0.2639253896763318e-03, 0.3552621469299578e+00, 0.2699137479265108e-03, 0.3906329503406230e+00, 0.2745196420166739e-03, 0.4251028614093031e+00, 0.2779529197397593e-03, 0.4584777520111870e+00, 0.2803996086684265e-03, 0.4905711358710193e+00, 0.2820302356715842e-03, 0.5212011669847385e+00, 0.2830056747491068e-03, 0.5501878488737995e+00, 0.2834808950776839e-03, 0.6025037877479342e+00, 0.2835282339078929e-03, 0.6254572689549016e+00, 0.2833819267065800e-03, 0.6460107179528248e+00, 0.2832858336906784e-03, 0.6639541138154251e+00, 0.2833268235451244e-03, 0.6790688515667495e+00, 0.2835432677029253e-03, 0.6911338580371512e+00, 0.2839091722743049e-03, 0.6999385956126490e+00, 0.2843308178875841e-03, 0.7053037748656896e+00, 0.2846703550533846e-03, 0.4732224387180115e-01, 0.1051193406971900e-03, 0.1202100529326803e+00, 0.1657871838796974e-03, 0.2034304820664855e+00, 0.2064648113714232e-03, 0.2912285643573002e+00, 0.2347942745819741e-03, 0.3802361792726768e+00, 0.2547775326597726e-03, 0.4680598511056146e+00, 0.2686876684847025e-03, 0.5528151052155599e+00, 0.2778665755515867e-03, 0.6329386307803041e+00, 0.2830996616782929e-03, 0.8056516651369069e-01, 0.2363454684003124e-01, 0.1403063340168372e-03, 0.1156476077139389e+00, 0.5191291632545936e-01, 0.1696504125939477e-03, 0.1520473382760421e+00, 0.8322715736994519e-01, 0.1935787242745390e-03, 0.1892986699745931e+00, 0.1165855667993712e+00, 0.2130614510521968e-03, 0.2270194446777792e+00, 0.1513077167409504e+00, 0.2289381265931048e-03, 0.2648908185093273e+00, 0.1868882025807859e+00, 0.2418630292816186e-03, 0.3026389259574136e+00, 0.2229277629776224e+00, 0.2523400495631193e-03, 0.3400220296151384e+00, 0.2590951840746235e+00, 0.2607623973449605e-03, 0.3768217953335510e+00, 0.2951047291750847e+00, 0.2674441032689209e-03, 0.4128372900921884e+00, 0.3307019714169930e+00, 0.2726432360343356e-03, 0.4478807131815630e+00, 0.3656544101087634e+00, 0.2765787685924545e-03, 0.4817742034089257e+00, 0.3997448951939695e+00, 0.2794428690642224e-03, 0.5143472814653344e+00, 0.4327667110812024e+00, 0.2814099002062895e-03, 0.5454346213905650e+00, 0.4645196123532293e+00, 0.2826429531578994e-03, 0.5748739313170252e+00, 0.4948063555703345e+00, 0.2832983542550884e-03, 0.1599598738286342e+00, 0.2792357590048985e-01, 0.1886695565284976e-03, 0.1998097412500951e+00, 0.5877141038139065e-01, 0.2081867882748234e-03, 0.2396228952566202e+00, 0.9164573914691377e-01, 0.2245148680600796e-03, 0.2792228341097746e+00, 0.1259049641962687e+00, 0.2380370491511872e-03, 0.3184251107546741e+00, 0.1610594823400863e+00, 0.2491398041852455e-03, 0.3570481164426244e+00, 0.1967151653460898e+00, 0.2581632405881230e-03, 0.3949164710492144e+00, 0.2325404606175168e+00, 0.2653965506227417e-03, 0.4318617293970503e+00, 0.2682461141151439e+00, 0.2710857216747087e-03, 0.4677221009931678e+00, 0.3035720116011973e+00, 0.2754434093903659e-03, 0.5023417939270955e+00, 0.3382781859197439e+00, 0.2786579932519380e-03, 0.5355701836636128e+00, 0.3721383065625942e+00, 0.2809011080679474e-03, 0.5672608451328771e+00, 0.4049346360466055e+00, 0.2823336184560987e-03, 0.5972704202540162e+00, 0.4364538098633802e+00, 0.2831101175806309e-03, 0.2461687022333596e+00, 0.3070423166833368e-01, 0.2221679970354546e-03, 0.2881774566286831e+00, 0.6338034669281885e-01, 0.2356185734270703e-03, 0.3293963604116978e+00, 0.9742862487067941e-01, 0.2469228344805590e-03, 0.3697303822241377e+00, 0.1323799532282290e+00, 0.2562726348642046e-03, 0.4090663023135127e+00, 0.1678497018129336e+00, 0.2638756726753028e-03, 0.4472819355411712e+00, 0.2035095105326114e+00, 0.2699311157390862e-03, 0.4842513377231437e+00, 0.2390692566672091e+00, 0.2746233268403837e-03, 0.5198477629962928e+00, 0.2742649818076149e+00, 0.2781225674454771e-03, 0.5539453011883145e+00, 0.3088503806580094e+00, 0.2805881254045684e-03, 0.5864196762401251e+00, 0.3425904245906614e+00, 0.2821719877004913e-03, 0.6171484466668390e+00, 0.3752562294789468e+00, 0.2830222502333124e-03, 0.3350337830565727e+00, 0.3261589934634747e-01, 0.2457995956744870e-03, 0.3775773224758284e+00, 0.6658438928081572e-01, 0.2551474407503706e-03, 0.4188155229848973e+00, 0.1014565797157954e+00, 0.2629065335195311e-03, 0.4586805892009344e+00, 0.1368573320843822e+00, 0.2691900449925075e-03, 0.4970895714224235e+00, 0.1724614851951608e+00, 0.2741275485754276e-03, 0.5339505133960747e+00, 0.2079779381416412e+00, 0.2778530970122595e-03, 0.5691665792531440e+00, 0.2431385788322288e+00, 0.2805010567646741e-03, 0.6026387682680377e+00, 0.2776901883049853e+00, 0.2822055834031040e-03, 0.6342676150163307e+00, 0.3113881356386632e+00, 0.2831016901243473e-03, 0.4237951119537067e+00, 0.3394877848664351e-01, 0.2624474901131803e-03, 0.4656918683234929e+00, 0.6880219556291447e-01, 0.2688034163039377e-03, 0.5058857069185980e+00, 0.1041946859721635e+00, 0.2738932751287636e-03, 0.5443204666713996e+00, 0.1398039738736393e+00, 0.2777944791242523e-03, 0.5809298813759742e+00, 0.1753373381196155e+00, 0.2806011661660987e-03, 0.6156416039447128e+00, 0.2105215793514010e+00, 0.2824181456597460e-03, 0.6483801351066604e+00, 0.2450953312157051e+00, 0.2833585216577828e-03, 0.5103616577251688e+00, 0.3485560643800719e-01, 0.2738165236962878e-03, 0.5506738792580681e+00, 0.7026308631512033e-01, 0.2778365208203180e-03, 0.5889573040995292e+00, 0.1059035061296403e+00, 0.2807852940418966e-03, 0.6251641589516930e+00, 0.1414823925236026e+00, 0.2827245949674705e-03, 0.6592414921570178e+00, 0.1767207908214530e+00, 0.2837342344829828e-03, 0.5930314017533384e+00, 0.3542189339561672e-01, 0.2809233907610981e-03, 0.6309812253390175e+00, 0.7109574040369549e-01, 0.2829930809742694e-03, 0.6666296011353230e+00, 0.1067259792282730e+00, 0.2841097874111479e-03, 0.6703715271049922e+00, 0.3569455268820809e-01, 0.2843455206008783e-03 };
  //rVecPt.resize( 11670 );
  //rVecWt.resize( 3890 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_4334_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 113, 1, 1, 0, 27, 9, 72, 4334 };
  double pa[] = { 0.1449063022537883e-04, 0.2546377329828424e-03, 0.1462896151831013e-01, 0.6018432961087496e-04, 0.3769840812493139e-01, 0.1002286583263673e-03, 0.6524701904096891e-01, 0.1315222931028093e-03, 0.9560543416134648e-01, 0.1564213746876724e-03, 0.1278335898929198e+00, 0.1765118841507736e-03, 0.1613096104466031e+00, 0.1928737099311080e-03, 0.1955806225745371e+00, 0.2062658534263270e-03, 0.2302935218498028e+00, 0.2172395445953787e-03, 0.2651584344113027e+00, 0.2262076188876047e-03, 0.2999276825183209e+00, 0.2334885699462397e-03, 0.3343828669718798e+00, 0.2393355273179203e-03, 0.3683265013750518e+00, 0.2439559200468863e-03, 0.4015763206518108e+00, 0.2475251866060002e-03, 0.4339612026399770e+00, 0.2501965558158773e-03, 0.4653180651114582e+00, 0.2521081407925925e-03, 0.4954893331080803e+00, 0.2533881002388081e-03, 0.5243207068924930e+00, 0.2541582900848261e-03, 0.5516590479041704e+00, 0.2545365737525860e-03, 0.6012371927804176e+00, 0.2545726993066799e-03, 0.6231574466449819e+00, 0.2544456197465555e-03, 0.6429416514181271e+00, 0.2543481596881064e-03, 0.6604124272943595e+00, 0.2543506451429194e-03, 0.6753851470408250e+00, 0.2544905675493763e-03, 0.6876717970626160e+00, 0.2547611407344429e-03, 0.6970895061319234e+00, 0.2551060375448869e-03, 0.7034746912553310e+00, 0.2554291933816039e-03, 0.7067017217542295e+00, 0.2556255710686343e-03, 0.4382223501131123e-01, 0.9041339695118195e-04, 0.1117474077400006e+00, 0.1438426330079022e-03, 0.1897153252911440e+00, 0.1802523089820518e-03, 0.2724023009910331e+00, 0.2060052290565496e-03, 0.3567163308709902e+00, 0.2245002248967466e-03, 0.4404784483028087e+00, 0.2377059847731150e-03, 0.5219833154161411e+00, 0.2468118955882525e-03, 0.5998179868977553e+00, 0.2525410872966528e-03, 0.6727803154548222e+00, 0.2553101409933397e-03, 0.7476563943166086e-01, 0.2193168509461185e-01, 0.1212879733668632e-03, 0.1075341482001416e+00, 0.4826419281533887e-01, 0.1472872881270931e-03, 0.1416344885203259e+00, 0.7751191883575742e-01, 0.1686846601010828e-03, 0.1766325315388586e+00, 0.1087558139247680e+00, 0.1862698414660208e-03, 0.2121744174481514e+00, 0.1413661374253096e+00, 0.2007430956991861e-03, 0.2479669443408145e+00, 0.1748768214258880e+00, 0.2126568125394796e-03, 0.2837600452294113e+00, 0.2089216406612073e+00, 0.2224394603372113e-03, 0.3193344933193984e+00, 0.2431987685545972e+00, 0.2304264522673135e-03, 0.3544935442438745e+00, 0.2774497054377770e+00, 0.2368854288424087e-03, 0.3890571932288154e+00, 0.3114460356156915e+00, 0.2420352089461772e-03, 0.4228581214259090e+00, 0.3449806851913012e+00, 0.2460597113081295e-03, 0.4557387211304052e+00, 0.3778618641248256e+00, 0.2491181912257687e-03, 0.4875487950541643e+00, 0.4099086391698978e+00, 0.2513528194205857e-03, 0.5181436529962997e+00, 0.4409474925853973e+00, 0.2528943096693220e-03, 0.5473824095600661e+00, 0.4708094517711291e+00, 0.2538660368488136e-03, 0.5751263398976174e+00, 0.4993275140354637e+00, 0.2543868648299022e-03, 0.1489515746840028e+00, 0.2599381993267017e-01, 0.1642595537825183e-03, 0.1863656444351767e+00, 0.5479286532462190e-01, 0.1818246659849308e-03, 0.2238602880356348e+00, 0.8556763251425254e-01, 0.1966565649492420e-03, 0.2612723375728160e+00, 0.1177257802267011e+00, 0.2090677905657991e-03, 0.2984332990206190e+00, 0.1508168456192700e+00, 0.2193820409510504e-03, 0.3351786584663333e+00, 0.1844801892177727e+00, 0.2278870827661928e-03, 0.3713505522209120e+00, 0.2184145236087598e+00, 0.2348283192282090e-03, 0.4067981098954663e+00, 0.2523590641486229e+00, 0.2404139755581477e-03, 0.4413769993687534e+00, 0.2860812976901373e+00, 0.2448227407760734e-03, 0.4749487182516394e+00, 0.3193686757808996e+00, 0.2482110455592573e-03, 0.5073798105075426e+00, 0.3520226949547602e+00, 0.2507192397774103e-03, 0.5385410448878654e+00, 0.3838544395667890e+00, 0.2524765968534880e-03, 0.5683065353670530e+00, 0.4146810037640963e+00, 0.2536052388539425e-03, 0.5965527620663510e+00, 0.4443224094681121e+00, 0.2542230588033068e-03, 0.2299227700856157e+00, 0.2865757664057584e-01, 0.1944817013047896e-03, 0.2695752998553267e+00, 0.5923421684485993e-01, 0.2067862362746635e-03, 0.3086178716611389e+00, 0.9117817776057715e-01, 0.2172440734649114e-03, 0.3469649871659077e+00, 0.1240593814082605e+00, 0.2260125991723423e-03, 0.3845153566319655e+00, 0.1575272058259175e+00, 0.2332655008689523e-03, 0.4211600033403215e+00, 0.1912845163525413e+00, 0.2391699681532458e-03, 0.4567867834329882e+00, 0.2250710177858171e+00, 0.2438801528273928e-03, 0.4912829319232061e+00, 0.2586521303440910e+00, 0.2475370504260665e-03, 0.5245364793303812e+00, 0.2918112242865407e+00, 0.2502707235640574e-03, 0.5564369788915756e+00, 0.3243439239067890e+00, 0.2522031701054241e-03, 0.5868757697775287e+00, 0.3560536787835351e+00, 0.2534511269978784e-03, 0.6157458853519617e+00, 0.3867480821242581e+00, 0.2541284914955151e-03, 0.3138461110672113e+00, 0.3051374637507278e-01, 0.2161509250688394e-03, 0.3542495872050569e+00, 0.6237111233730755e-01, 0.2248778513437852e-03, 0.3935751553120181e+00, 0.9516223952401907e-01, 0.2322388803404617e-03, 0.4317634668111147e+00, 0.1285467341508517e+00, 0.2383265471001355e-03, 0.4687413842250821e+00, 0.1622318931656033e+00, 0.2432476675019525e-03, 0.5044274237060283e+00, 0.1959581153836453e+00, 0.2471122223750674e-03, 0.5387354077925727e+00, 0.2294888081183837e+00, 0.2500291752486870e-03, 0.5715768898356105e+00, 0.2626031152713945e+00, 0.2521055942764682e-03, 0.6028627200136111e+00, 0.2950904075286713e+00, 0.2534472785575503e-03, 0.6325039812653463e+00, 0.3267458451113286e+00, 0.2541599713080121e-03, 0.3981986708423407e+00, 0.3183291458749821e-01, 0.2317380975862936e-03, 0.4382791182133300e+00, 0.6459548193880908e-01, 0.2378550733719775e-03, 0.4769233057218166e+00, 0.9795757037087952e-01, 0.2428884456739118e-03, 0.5140823911194238e+00, 0.1316307235126655e+00, 0.2469002655757292e-03, 0.5496977833862983e+00, 0.1653556486358704e+00, 0.2499657574265851e-03, 0.5837047306512727e+00, 0.1988931724126510e+00, 0.2521676168486082e-03, 0.6160349566926879e+00, 0.2320174581438950e+00, 0.2535935662645334e-03, 0.6466185353209440e+00, 0.2645106562168662e+00, 0.2543356743363214e-03, 0.4810835158795404e+00, 0.3275917807743992e-01, 0.2427353285201535e-03, 0.5199925041324341e+00, 0.6612546183967181e-01, 0.2468258039744386e-03, 0.5571717692207494e+00, 0.9981498331474143e-01, 0.2500060956440310e-03, 0.5925789250836378e+00, 0.1335687001410374e+00, 0.2523238365420979e-03, 0.6261658523859670e+00, 0.1671444402896463e+00, 0.2538399260252846e-03, 0.6578811126669331e+00, 0.2003106382156076e+00, 0.2546255927268069e-03, 0.5609624612998100e+00, 0.3337500940231335e-01, 0.2500583360048449e-03, 0.5979959659984670e+00, 0.6708750335901803e-01, 0.2524777638260203e-03, 0.6330523711054002e+00, 0.1008792126424850e+00, 0.2540951193860656e-03, 0.6660960998103972e+00, 0.1345050343171794e+00, 0.2549524085027472e-03, 0.6365384364585819e+00, 0.3372799460737052e-01, 0.2542569507009158e-03, 0.6710994302899275e+00, 0.6755249309678028e-01, 0.2552114127580376e-03 };
  //rVecPt.resize( 13002 );
  //rVecWt.resize( 4334 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_4802_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 119, 1, 1, 1, 28, 9, 81, 4802 };
  double pa[] = { 0.9687521879420705e-04, 0.2297310852498558e-03, 0.2307897895367918e-03, 0.2335728608887064e-01, 0.7386265944001919e-04, 0.4352987836550653e-01, 0.8257977698542210e-04, 0.6439200521088801e-01, 0.9706044762057630e-04, 0.9003943631993181e-01, 0.1302393847117003e-03, 0.1196706615548473e+00, 0.1541957004600968e-03, 0.1511715412838134e+00, 0.1704459770092199e-03, 0.1835982828503801e+00, 0.1827374890942906e-03, 0.2165081259155405e+00, 0.1926360817436107e-03, 0.2496208720417563e+00, 0.2008010239494833e-03, 0.2827200673567900e+00, 0.2075635983209175e-03, 0.3156190823994346e+00, 0.2131306638690909e-03, 0.3481476793749115e+00, 0.2176562329937335e-03, 0.3801466086947226e+00, 0.2212682262991018e-03, 0.4114652119634011e+00, 0.2240799515668565e-03, 0.4419598786519751e+00, 0.2261959816187525e-03, 0.4714925949329543e+00, 0.2277156368808855e-03, 0.4999293972879466e+00, 0.2287351772128336e-03, 0.5271387221431248e+00, 0.2293490814084085e-03, 0.5529896780837761e+00, 0.2296505312376273e-03, 0.6000856099481712e+00, 0.2296793832318756e-03, 0.6210562192785175e+00, 0.2295785443842974e-03, 0.6401165879934240e+00, 0.2295017931529102e-03, 0.6571144029244334e+00, 0.2295059638184868e-03, 0.6718910821718863e+00, 0.2296232343237362e-03, 0.6842845591099010e+00, 0.2298530178740771e-03, 0.6941353476269816e+00, 0.2301579790280501e-03, 0.7012965242212991e+00, 0.2304690404996513e-03, 0.7056471428242644e+00, 0.2307027995907102e-03, 0.4595557643585895e-01, 0.9312274696671092e-04, 0.1049316742435023e+00, 0.1199919385876926e-03, 0.1773548879549274e+00, 0.1598039138877690e-03, 0.2559071411236127e+00, 0.1822253763574900e-03, 0.3358156837985898e+00, 0.1988579593655040e-03, 0.4155835743763893e+00, 0.2112620102533307e-03, 0.4937894296167472e+00, 0.2201594887699007e-03, 0.5691569694793316e+00, 0.2261622590895036e-03, 0.6405840854894251e+00, 0.2296458453435705e-03, 0.7345133894143348e-01, 0.2177844081486067e-01, 0.1006006990267000e-03, 0.1009859834044931e+00, 0.4590362185775188e-01, 0.1227676689635876e-03, 0.1324289619748758e+00, 0.7255063095690877e-01, 0.1467864280270117e-03, 0.1654272109607127e+00, 0.1017825451960684e+00, 0.1644178912101232e-03, 0.1990767186776461e+00, 0.1325652320980364e+00, 0.1777664890718961e-03, 0.2330125945523278e+00, 0.1642765374496765e+00, 0.1884825664516690e-03, 0.2670080611108287e+00, 0.1965360374337889e+00, 0.1973269246453848e-03, 0.3008753376294316e+00, 0.2290726770542238e+00, 0.2046767775855328e-03, 0.3344475596167860e+00, 0.2616645495370823e+00, 0.2107600125918040e-03, 0.3675709724070786e+00, 0.2941150728843141e+00, 0.2157416362266829e-03, 0.4001000887587812e+00, 0.3262440400919066e+00, 0.2197557816920721e-03, 0.4318956350436028e+00, 0.3578835350611916e+00, 0.2229192611835437e-03, 0.4628239056795531e+00, 0.3888751854043678e+00, 0.2253385110212775e-03, 0.4927563229773636e+00, 0.4190678003222840e+00, 0.2271137107548774e-03, 0.5215687136707969e+00, 0.4483151836883852e+00, 0.2283414092917525e-03, 0.5491402346984905e+00, 0.4764740676087880e+00, 0.2291161673130077e-03, 0.5753520160126075e+00, 0.5034021310998277e+00, 0.2295313908576598e-03, 0.1388326356417754e+00, 0.2435436510372806e-01, 0.1438204721359031e-03, 0.1743686900537244e+00, 0.5118897057342652e-01, 0.1607738025495257e-03, 0.2099737037950268e+00, 0.8014695048539634e-01, 0.1741483853528379e-03, 0.2454492590908548e+00, 0.1105117874155699e+00, 0.1851918467519151e-03, 0.2807219257864278e+00, 0.1417950531570966e+00, 0.1944628638070613e-03, 0.3156842271975842e+00, 0.1736604945719597e+00, 0.2022495446275152e-03, 0.3502090945177752e+00, 0.2058466324693981e+00, 0.2087462382438514e-03, 0.3841684849519686e+00, 0.2381284261195919e+00, 0.2141074754818308e-03, 0.4174372367906016e+00, 0.2703031270422569e+00, 0.2184640913748162e-03, 0.4498926465011892e+00, 0.3021845683091309e+00, 0.2219309165220329e-03, 0.4814146229807701e+00, 0.3335993355165720e+00, 0.2246123118340624e-03, 0.5118863625734701e+00, 0.3643833735518232e+00, 0.2266062766915125e-03, 0.5411947455119144e+00, 0.3943789541958179e+00, 0.2280072952230796e-03, 0.5692301500357246e+00, 0.4234320144403542e+00, 0.2289082025202583e-03, 0.5958857204139576e+00, 0.4513897947419260e+00, 0.2294012695120025e-03, 0.2156270284785766e+00, 0.2681225755444491e-01, 0.1722434488736947e-03, 0.2532385054909710e+00, 0.5557495747805614e-01, 0.1830237421455091e-03, 0.2902564617771537e+00, 0.8569368062950249e-01, 0.1923855349997633e-03, 0.3266979823143256e+00, 0.1167367450324135e+00, 0.2004067861936271e-03, 0.3625039627493614e+00, 0.1483861994003304e+00, 0.2071817297354263e-03, 0.3975838937548699e+00, 0.1803821503011405e+00, 0.2128250834102103e-03, 0.4318396099009774e+00, 0.2124962965666424e+00, 0.2174513719440102e-03, 0.4651706555732742e+00, 0.2445221837805913e+00, 0.2211661839150214e-03, 0.4974752649620969e+00, 0.2762701224322987e+00, 0.2240665257813102e-03, 0.5286517579627517e+00, 0.3075627775211328e+00, 0.2262439516632620e-03, 0.5586001195731895e+00, 0.3382311089826877e+00, 0.2277874557231869e-03, 0.5872229902021319e+00, 0.3681108834741399e+00, 0.2287854314454994e-03, 0.6144258616235123e+00, 0.3970397446872839e+00, 0.2293268499615575e-03, 0.2951676508064861e+00, 0.2867499538750441e-01, 0.1912628201529828e-03, 0.3335085485472725e+00, 0.5867879341903510e-01, 0.1992499672238701e-03, 0.3709561760636381e+00, 0.8961099205022284e-01, 0.2061275533454027e-03, 0.4074722861667498e+00, 0.1211627927626297e+00, 0.2119318215968572e-03, 0.4429923648839117e+00, 0.1530748903554898e+00, 0.2167416581882652e-03, 0.4774428052721736e+00, 0.1851176436721877e+00, 0.2206430730516600e-03, 0.5107446539535904e+00, 0.2170829107658179e+00, 0.2237186938699523e-03, 0.5428151370542935e+00, 0.2487786689026271e+00, 0.2260480075032884e-03, 0.5735699292556964e+00, 0.2800239952795016e+00, 0.2277098884558542e-03, 0.6029253794562866e+00, 0.3106445702878119e+00, 0.2287845715109671e-03, 0.6307998987073145e+00, 0.3404689500841194e+00, 0.2293547268236294e-03, 0.3752652273692719e+00, 0.2997145098184479e-01, 0.2056073839852528e-03, 0.4135383879344028e+00, 0.6086725898678011e-01, 0.2114235865831876e-03, 0.4506113885153907e+00, 0.9238849548435643e-01, 0.2163175629770551e-03, 0.4864401554606072e+00, 0.1242786603851851e+00, 0.2203392158111650e-03, 0.5209708076611709e+00, 0.1563086731483386e+00, 0.2235473176847839e-03, 0.5541422135830122e+00, 0.1882696509388506e+00, 0.2260024141501235e-03, 0.5858880915113817e+00, 0.2199672979126059e+00, 0.2277675929329182e-03, 0.6161399390603444e+00, 0.2512165482924867e+00, 0.2289102112284834e-03, 0.6448296482255090e+00, 0.2818368701871888e+00, 0.2295027954625118e-03, 0.4544796274917948e+00, 0.3088970405060312e-01, 0.2161281589879992e-03, 0.4919389072146628e+00, 0.6240947677636835e-01, 0.2201980477395102e-03, 0.5279313026985183e+00, 0.9430706144280313e-01, 0.2234952066593166e-03, 0.5624169925571135e+00, 0.1263547818770374e+00, 0.2260540098520838e-03, 0.5953484627093287e+00, 0.1583430788822594e+00, 0.2279157981899988e-03, 0.6266730715339185e+00, 0.1900748462555988e+00, 0.2291296918565571e-03, 0.6563363204278871e+00, 0.2213599519592567e+00, 0.2297533752536649e-03, 0.5314574716585696e+00, 0.3152508811515374e-01, 0.2234927356465995e-03, 0.5674614932298185e+00, 0.6343865291465561e-01, 0.2261288012985219e-03, 0.6017706004970264e+00, 0.9551503504223951e-01, 0.2280818160923688e-03, 0.6343471270264178e+00, 0.1275440099801196e+00, 0.2293773295180159e-03, 0.6651494599127802e+00, 0.1593252037671960e+00, 0.2300528767338634e-03, 0.6050184986005704e+00, 0.3192538338496105e-01, 0.2281893855065666e-03, 0.6390163550880400e+00, 0.6402824353962306e-01, 0.2295720444840727e-03, 0.6711199107088448e+00, 0.9609805077002909e-01, 0.2303227649026753e-03, 0.6741354429572275e+00, 0.3211853196273233e-01, 0.2304831913227114e-03 };
  //rVecPt.resize( 14406 );
  //rVecWt.resize( 4802 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_5294_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 125, 1, 1, 0, 30, 10, 90, 5294 };
  double pa[] = { 0.9080510764308163e-04, 0.2084824361987793e-03, 0.2303261686261450e-01, 0.5011105657239616e-04, 0.3757208620162394e-01, 0.5942520409683854e-04, 0.5821912033821852e-01, 0.9564394826109721e-04, 0.8403127529194872e-01, 0.1185530657126338e-03, 0.1122927798060578e+00, 0.1364510114230331e-03, 0.1420125319192987e+00, 0.1505828825605415e-03, 0.1726396437341978e+00, 0.1619298749867023e-03, 0.2038170058115696e+00, 0.1712450504267789e-03, 0.2352849892876508e+00, 0.1789891098164999e-03, 0.2668363354312461e+00, 0.1854474955629795e-03, 0.2982941279900452e+00, 0.1908148636673661e-03, 0.3295002922087076e+00, 0.1952377405281833e-03, 0.3603094918363593e+00, 0.1988349254282232e-03, 0.3905857895173920e+00, 0.2017079807160050e-03, 0.4202005758160837e+00, 0.2039473082709094e-03, 0.4490310061597227e+00, 0.2056360279288953e-03, 0.4769586160311491e+00, 0.2068525823066865e-03, 0.5038679887049750e+00, 0.2076724877534488e-03, 0.5296454286519961e+00, 0.2081694278237885e-03, 0.5541776207164850e+00, 0.2084157631219326e-03, 0.5990467321921213e+00, 0.2084381531128593e-03, 0.6191467096294587e+00, 0.2083476277129307e-03, 0.6375251212901849e+00, 0.2082686194459732e-03, 0.6540514381131168e+00, 0.2082475686112415e-03, 0.6685899064391510e+00, 0.2083139860289915e-03, 0.6810013009681648e+00, 0.2084745561831237e-03, 0.6911469578730340e+00, 0.2087091313375890e-03, 0.6988956915141736e+00, 0.2089718413297697e-03, 0.7041335794868720e+00, 0.2092003303479793e-03, 0.7067754398018567e+00, 0.2093336148263241e-03, 0.3840368707853623e-01, 0.7591708117365267e-04, 0.9835485954117399e-01, 0.1083383968169186e-03, 0.1665774947612998e+00, 0.1403019395292510e-03, 0.2405702335362910e+00, 0.1615970179286436e-03, 0.3165270770189046e+00, 0.1771144187504911e-03, 0.3927386145645443e+00, 0.1887760022988168e-03, 0.4678825918374656e+00, 0.1973474670768214e-03, 0.5408022024266935e+00, 0.2033787661234659e-03, 0.6104967445752438e+00, 0.2072343626517331e-03, 0.6760910702685738e+00, 0.2091177834226918e-03, 0.6655644120217392e-01, 0.1936508874588424e-01, 0.9316684484675566e-04, 0.9446246161270182e-01, 0.4252442002115869e-01, 0.1116193688682976e-03, 0.1242651925452509e+00, 0.6806529315354374e-01, 0.1298623551559414e-03, 0.1553438064846751e+00, 0.9560957491205369e-01, 0.1450236832456426e-03, 0.1871137110542670e+00, 0.1245931657452888e+00, 0.1572719958149914e-03, 0.2192612628836257e+00, 0.1545385828778978e+00, 0.1673234785867195e-03, 0.2515682807206955e+00, 0.1851004249723368e+00, 0.1756860118725188e-03, 0.2838535866287290e+00, 0.2160182608272384e+00, 0.1826776290439367e-03, 0.3159578817528521e+00, 0.2470799012277111e+00, 0.1885116347992865e-03, 0.3477370882791392e+00, 0.2781014208986402e+00, 0.1933457860170574e-03, 0.3790576960890540e+00, 0.3089172523515731e+00, 0.1973060671902064e-03, 0.4097938317810200e+00, 0.3393750055472244e+00, 0.2004987099616311e-03, 0.4398256572859637e+00, 0.3693322470987730e+00, 0.2030170909281499e-03, 0.4690384114718480e+00, 0.3986541005609877e+00, 0.2049461460119080e-03, 0.4973216048301053e+00, 0.4272112491408562e+00, 0.2063653565200186e-03, 0.5245681526132446e+00, 0.4548781735309936e+00, 0.2073507927381027e-03, 0.5506733911803888e+00, 0.4815315355023251e+00, 0.2079764593256122e-03, 0.5755339829522475e+00, 0.5070486445801855e+00, 0.2083150534968778e-03, 0.1305472386056362e+00, 0.2284970375722366e-01, 0.1262715121590664e-03, 0.1637327908216477e+00, 0.4812254338288384e-01, 0.1414386128545972e-03, 0.1972734634149637e+00, 0.7531734457511935e-01, 0.1538740401313898e-03, 0.2308694653110130e+00, 0.1039043639882017e+00, 0.1642434942331432e-03, 0.2643899218338160e+00, 0.1334526587117626e+00, 0.1729790609237496e-03, 0.2977171599622171e+00, 0.1636414868936382e+00, 0.1803505190260828e-03, 0.3307293903032310e+00, 0.1942195406166568e+00, 0.1865475350079657e-03, 0.3633069198219073e+00, 0.2249752879943753e+00, 0.1917182669679069e-03, 0.3953346955922727e+00, 0.2557218821820032e+00, 0.1959851709034382e-03, 0.4267018394184914e+00, 0.2862897925213193e+00, 0.1994529548117882e-03, 0.4573009622571704e+00, 0.3165224536636518e+00, 0.2022138911146548e-03, 0.4870279559856109e+00, 0.3462730221636496e+00, 0.2043518024208592e-03, 0.5157819581450322e+00, 0.3754016870282835e+00, 0.2059450313018110e-03, 0.5434651666465393e+00, 0.4037733784993613e+00, 0.2070685715318472e-03, 0.5699823887764627e+00, 0.4312557784139123e+00, 0.2077955310694373e-03, 0.5952403350947741e+00, 0.4577175367122110e+00, 0.2081980387824712e-03, 0.2025152599210369e+00, 0.2520253617719557e-01, 0.1521318610377956e-03, 0.2381066653274425e+00, 0.5223254506119000e-01, 0.1622772720185755e-03, 0.2732823383651612e+00, 0.8060669688588620e-01, 0.1710498139420709e-03, 0.3080137692611118e+00, 0.1099335754081255e+00, 0.1785911149448736e-03, 0.3422405614587601e+00, 0.1399120955959857e+00, 0.1850125313687736e-03, 0.3758808773890420e+00, 0.1702977801651705e+00, 0.1904229703933298e-03, 0.4088458383438932e+00, 0.2008799256601680e+00, 0.1949259956121987e-03, 0.4410450550841152e+00, 0.2314703052180836e+00, 0.1986161545363960e-03, 0.4723879420561312e+00, 0.2618972111375892e+00, 0.2015790585641370e-03, 0.5027843561874343e+00, 0.2920013195600270e+00, 0.2038934198707418e-03, 0.5321453674452458e+00, 0.3216322555190551e+00, 0.2056334060538251e-03, 0.5603839113834030e+00, 0.3506456615934198e+00, 0.2068705959462289e-03, 0.5874150706875146e+00, 0.3789007181306267e+00, 0.2076753906106002e-03, 0.6131559381660038e+00, 0.4062580170572782e+00, 0.2081179391734803e-03, 0.2778497016394506e+00, 0.2696271276876226e-01, 0.1700345216228943e-03, 0.3143733562261912e+00, 0.5523469316960465e-01, 0.1774906779990410e-03, 0.3501485810261827e+00, 0.8445193201626464e-01, 0.1839659377002642e-03, 0.3851430322303653e+00, 0.1143263119336083e+00, 0.1894987462975169e-03, 0.4193013979470415e+00, 0.1446177898344475e+00, 0.1941548809452595e-03, 0.4525585960458567e+00, 0.1751165438438091e+00, 0.1980078427252384e-03, 0.4848447779622947e+00, 0.2056338306745660e+00, 0.2011296284744488e-03, 0.5160871208276894e+00, 0.2359965487229226e+00, 0.2035888456966776e-03, 0.5462112185696926e+00, 0.2660430223139146e+00, 0.2054516325352142e-03, 0.5751425068101757e+00, 0.2956193664498032e+00, 0.2067831033092635e-03, 0.6028073872853596e+00, 0.3245763905312779e+00, 0.2076485320284876e-03, 0.6291338275278409e+00, 0.3527670026206972e+00, 0.2081141439525255e-03, 0.3541797528439391e+00, 0.2823853479435550e-01, 0.1834383015469222e-03, 0.3908234972074657e+00, 0.5741296374713106e-01, 0.1889540591777677e-03, 0.4264408450107590e+00, 0.8724646633650199e-01, 0.1936677023597375e-03, 0.4609949666553286e+00, 0.1175034422915616e+00, 0.1976176495066504e-03, 0.4944389496536006e+00, 0.1479755652628428e+00, 0.2008536004560983e-03, 0.5267194884346086e+00, 0.1784740659484352e+00, 0.2034280351712291e-03, 0.5577787810220990e+00, 0.2088245700431244e+00, 0.2053944466027758e-03, 0.5875563763536670e+00, 0.2388628136570763e+00, 0.2068077642882360e-03, 0.6159910016391269e+00, 0.2684308928769185e+00, 0.2077250949661599e-03, 0.6430219602956268e+00, 0.2973740761960252e+00, 0.2082062440705320e-03, 0.4300647036213646e+00, 0.2916399920493977e-01, 0.1934374486546626e-03, 0.4661486308935531e+00, 0.5898803024755659e-01, 0.1974107010484300e-03, 0.5009658555287261e+00, 0.8924162698525409e-01, 0.2007129290388658e-03, 0.5344824270447704e+00, 0.1197185199637321e+00, 0.2033736947471293e-03, 0.5666575997416371e+00, 0.1502300756161382e+00, 0.2054287125902493e-03, 0.5974457471404752e+00, 0.1806004191913564e+00, 0.2069184936818894e-03, 0.6267984444116886e+00, 0.2106621764786252e+00, 0.2078883689808782e-03, 0.6546664713575417e+00, 0.2402526932671914e+00, 0.2083886366116359e-03, 0.5042711004437253e+00, 0.2982529203607657e-01, 0.2006593275470817e-03, 0.5392127456774380e+00, 0.6008728062339922e-01, 0.2033728426135397e-03, 0.5726819437668618e+00, 0.9058227674571398e-01, 0.2055008781377608e-03, 0.6046469254207278e+00, 0.1211219235803400e+00, 0.2070651783518502e-03, 0.6350716157434952e+00, 0.1515286404791580e+00, 0.2080953335094320e-03, 0.6639177679185454e+00, 0.1816314681255552e+00, 0.2086284998988521e-03, 0.5757276040972253e+00, 0.3026991752575440e-01, 0.2055549387644668e-03, 0.6090265823139755e+00, 0.6078402297870770e-01, 0.2071871850267654e-03, 0.6406735344387661e+00, 0.9135459984176636e-01, 0.2082856600431965e-03, 0.6706397927793709e+00, 0.1218024155966590e+00, 0.2088705858819358e-03, 0.6435019674426665e+00, 0.3052608357660639e-01, 0.2083995867536322e-03, 0.6747218676375681e+00, 0.6112185773983089e-01, 0.2090509712889637e-03 };
  //rVecPt.resize( 15882 );
  //rVecWt.resize( 5294 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}


void ccdl::Backend::_Lebedev_5810_pt( std::tr1::array<double,3> * rVecPt,double * rVecWt )
{
  using ccdl::Backend::_Lebedev_Oh_Mesh;

  int nk[] = { 131, 1, 1, 1, 31, 10, 100, 5810 };
  double pa[] = { 0.9735347946175486e-05, 0.1901059546737578e-03, 0.1907581241803167e-03, 0.1182361662400277e-01, 0.3926424538919212e-04, 0.3062145009138958e-01, 0.6667905467294382e-04, 0.5329794036834243e-01, 0.8868891315019135e-04, 0.7848165532862220e-01, 0.1066306000958872e-03, 0.1054038157636201e+00, 0.1214506743336128e-03, 0.1335577797766211e+00, 0.1338054681640871e-03, 0.1625769955502252e+00, 0.1441677023628504e-03, 0.1921787193412792e+00, 0.1528880200826557e-03, 0.2221340534690548e+00, 0.1602330623773609e-03, 0.2522504912791132e+00, 0.1664102653445244e-03, 0.2823610860679697e+00, 0.1715845854011323e-03, 0.3123173966267560e+00, 0.1758901000133069e-03, 0.3419847036953789e+00, 0.1794382485256736e-03, 0.3712386456999758e+00, 0.1823238106757407e-03, 0.3999627649876828e+00, 0.1846293252959976e-03, 0.4280466458648093e+00, 0.1864284079323098e-03, 0.4553844360185711e+00, 0.1877882694626914e-03, 0.4818736094437834e+00, 0.1887716321852025e-03, 0.5074138709260629e+00, 0.1894381638175673e-03, 0.5319061304570707e+00, 0.1898454899533629e-03, 0.5552514978677286e+00, 0.1900497929577815e-03, 0.5981009025246183e+00, 0.1900671501924092e-03, 0.6173990192228116e+00, 0.1899837555533510e-03, 0.6351365239411131e+00, 0.1899014113156229e-03, 0.6512010228227200e+00, 0.1898581257705106e-03, 0.6654758363948120e+00, 0.1898804756095753e-03, 0.6778410414853370e+00, 0.1899793610426402e-03, 0.6881760887484110e+00, 0.1901464554844117e-03, 0.6963645267094598e+00, 0.1903533246259542e-03, 0.7023010617153579e+00, 0.1905556158463228e-03, 0.7059004636628753e+00, 0.1907037155663528e-03, 0.3552470312472575e-01, 0.5992997844249967e-04, 0.9151176620841283e-01, 0.9749059382456978e-04, 0.1566197930068980e+00, 0.1241680804599158e-03, 0.2265467599271907e+00, 0.1437626154299360e-03, 0.2988242318581361e+00, 0.1584200054793902e-03, 0.3717482419703886e+00, 0.1694436550982744e-03, 0.4440094491758889e+00, 0.1776617014018108e-03, 0.5145337096756642e+00, 0.1836132434440077e-03, 0.5824053672860230e+00, 0.1876494727075983e-03, 0.6468283961043370e+00, 0.1899906535336482e-03, 0.6095964259104373e-01, 0.1787828275342931e-01, 0.8143252820767350e-04, 0.8811962270959388e-01, 0.3953888740792096e-01, 0.9998859890887728e-04, 0.1165936722428831e+00, 0.6378121797722990e-01, 0.1156199403068359e-03, 0.1460232857031785e+00, 0.8985890813745037e-01, 0.1287632092635513e-03, 0.1761197110181755e+00, 0.1172606510576162e+00, 0.1398378643365139e-03, 0.2066471190463718e+00, 0.1456102876970995e+00, 0.1491876468417391e-03, 0.2374076026328152e+00, 0.1746153823011775e+00, 0.1570855679175456e-03, 0.2682305474337051e+00, 0.2040383070295584e+00, 0.1637483948103775e-03, 0.2989653312142369e+00, 0.2336788634003698e+00, 0.1693500566632843e-03, 0.3294762752772209e+00, 0.2633632752654219e+00, 0.1740322769393633e-03, 0.3596390887276086e+00, 0.2929369098051601e+00, 0.1779126637278296e-03, 0.3893383046398812e+00, 0.3222592785275512e+00, 0.1810908108835412e-03, 0.4184653789358347e+00, 0.3512004791195743e+00, 0.1836529132600190e-03, 0.4469172319076166e+00, 0.3796385677684537e+00, 0.1856752841777379e-03, 0.4745950813276976e+00, 0.4074575378263879e+00, 0.1872270566606832e-03, 0.5014034601410262e+00, 0.4345456906027828e+00, 0.1883722645591307e-03, 0.5272493404551239e+00, 0.4607942515205134e+00, 0.1891714324525297e-03, 0.5520413051846366e+00, 0.4860961284181720e+00, 0.1896827480450146e-03, 0.5756887237503077e+00, 0.5103447395342790e+00, 0.1899628417059528e-03, 0.1225039430588352e+00, 0.2136455922655793e-01, 0.1123301829001669e-03, 0.1539113217321372e+00, 0.4520926166137188e-01, 0.1253698826711277e-03, 0.1856213098637712e+00, 0.7086468177864818e-01, 0.1366266117678531e-03, 0.2174998728035131e+00, 0.9785239488772918e-01, 0.1462736856106918e-03, 0.2494128336938330e+00, 0.1258106396267210e+00, 0.1545076466685412e-03, 0.2812321562143480e+00, 0.1544529125047001e+00, 0.1615096280814007e-03, 0.3128372276456111e+00, 0.1835433512202753e+00, 0.1674366639741759e-03, 0.3441145160177973e+00, 0.2128813258619585e+00, 0.1724225002437900e-03, 0.3749567714853510e+00, 0.2422913734880829e+00, 0.1765810822987288e-03, 0.4052621732015610e+00, 0.2716163748391453e+00, 0.1800104126010751e-03, 0.4349335453522385e+00, 0.3007127671240280e+00, 0.1827960437331284e-03, 0.4638776641524965e+00, 0.3294470677216479e+00, 0.1850140300716308e-03, 0.4920046410462687e+00, 0.3576932543699155e+00, 0.1867333507394938e-03, 0.5192273554861704e+00, 0.3853307059757764e+00, 0.1880178688638289e-03, 0.5454609081136522e+00, 0.4122425044452694e+00, 0.1889278925654758e-03, 0.5706220661424140e+00, 0.4383139587781027e+00, 0.1895213832507346e-03, 0.5946286755181518e+00, 0.4634312536300553e+00, 0.1898548277397420e-03, 0.1905370790924295e+00, 0.2371311537781979e-01, 0.1349105935937341e-03, 0.2242518717748009e+00, 0.4917878059254806e-01, 0.1444060068369326e-03, 0.2577190808025936e+00, 0.7595498960495142e-01, 0.1526797390930008e-03, 0.2908724534927187e+00, 0.1036991083191100e+00, 0.1598208771406474e-03, 0.3236354020056219e+00, 0.1321348584450234e+00, 0.1659354368615331e-03, 0.3559267359304543e+00, 0.1610316571314789e+00, 0.1711279910946440e-03, 0.3876637123676956e+00, 0.1901912080395707e+00, 0.1754952725601440e-03, 0.4187636705218842e+00, 0.2194384950137950e+00, 0.1791247850802529e-03, 0.4491449019883107e+00, 0.2486155334763858e+00, 0.1820954300877716e-03, 0.4787270932425445e+00, 0.2775768931812335e+00, 0.1844788524548449e-03, 0.5074315153055574e+00, 0.3061863786591120e+00, 0.1863409481706220e-03, 0.5351810507738336e+00, 0.3343144718152556e+00, 0.1877433008795068e-03, 0.5619001025975381e+00, 0.3618362729028427e+00, 0.1887444543705232e-03, 0.5875144035268046e+00, 0.3886297583620408e+00, 0.1894009829375006e-03, 0.6119507308734495e+00, 0.4145742277792031e+00, 0.1897683345035198e-03, 0.2619733870119463e+00, 0.2540047186389353e-01, 0.1517327037467653e-03, 0.2968149743237949e+00, 0.5208107018543989e-01, 0.1587740557483543e-03, 0.3310451504860488e+00, 0.7971828470885599e-01, 0.1649093382274097e-03, 0.3646215567376676e+00, 0.1080465999177927e+00, 0.1701915216193265e-03, 0.3974916785279360e+00, 0.1368413849366629e+00, 0.1746847753144065e-03, 0.4295967403772029e+00, 0.1659073184763559e+00, 0.1784555512007570e-03, 0.4608742854473447e+00, 0.1950703730454614e+00, 0.1815687562112174e-03, 0.4912598858949903e+00, 0.2241721144376724e+00, 0.1840864370663302e-03, 0.5206882758945558e+00, 0.2530655255406489e+00, 0.1860676785390006e-03, 0.5490940914019819e+00, 0.2816118409731066e+00, 0.1875690583743703e-03, 0.5764123302025542e+00, 0.3096780504593238e+00, 0.1886453236347225e-03, 0.6025786004213506e+00, 0.3371348366394987e+00, 0.1893501123329645e-03, 0.6275291964794956e+00, 0.3638547827694396e+00, 0.1897366184519868e-03, 0.3348189479861771e+00, 0.2664841935537443e-01, 0.1643908815152736e-03, 0.3699515545855295e+00, 0.5424000066843495e-01, 0.1696300350907768e-03, 0.4042003071474669e+00, 0.8251992715430854e-01, 0.1741553103844483e-03, 0.4375320100182624e+00, 0.1112695182483710e+00, 0.1780015282386092e-03, 0.4699054490335947e+00, 0.1402964116467816e+00, 0.1812116787077125e-03, 0.5012739879431952e+00, 0.1694275117584291e+00, 0.1838323158085421e-03, 0.5315874883754966e+00, 0.1985038235312689e+00, 0.1859113119837737e-03, 0.5607937109622117e+00, 0.2273765660020893e+00, 0.1874969220221698e-03, 0.5888393223495521e+00, 0.2559041492849764e+00, 0.1886375612681076e-03, 0.6156705979160163e+00, 0.2839497251976899e+00, 0.1893819575809276e-03, 0.6412338809078123e+00, 0.3113791060500690e+00, 0.1897794748256767e-03, 0.4076051259257167e+00, 0.2757792290858463e-01, 0.1738963926584846e-03, 0.4423788125791520e+00, 0.5584136834984293e-01, 0.1777442359873466e-03, 0.4760480917328258e+00, 0.8457772087727143e-01, 0.1810010815068719e-03, 0.5085838725946297e+00, 0.1135975846359248e+00, 0.1836920318248129e-03, 0.5399513637391218e+00, 0.1427286904765053e+00, 0.1858489473214328e-03, 0.5701118433636380e+00, 0.1718112740057635e+00, 0.1875079342496592e-03, 0.5990240530606021e+00, 0.2006944855985351e+00, 0.1887080239102310e-03, 0.6266452685139695e+00, 0.2292335090598907e+00, 0.1894905752176822e-03, 0.6529320971415942e+00, 0.2572871512353714e+00, 0.1898991061200695e-03, 0.4791583834610126e+00, 0.2826094197735932e-01, 0.1809065016458791e-03, 0.5130373952796940e+00, 0.5699871359683649e-01, 0.1836297121596799e-03, 0.5456252429628476e+00, 0.8602712528554394e-01, 0.1858426916241869e-03, 0.5768956329682385e+00, 0.1151748137221281e+00, 0.1875654101134641e-03, 0.6068186944699046e+00, 0.1442811654136362e+00, 0.1888240751833503e-03, 0.6353622248024907e+00, 0.1731930321657680e+00, 0.1896497383866979e-03, 0.6624927035731797e+00, 0.2017619958756061e+00, 0.1900775530219121e-03, 0.5484933508028488e+00, 0.2874219755907391e-01, 0.1858525041478814e-03, 0.5810207682142106e+00, 0.5778312123713695e-01, 0.1876248690077947e-03, 0.6120955197181352e+00, 0.8695262371439526e-01, 0.1889404439064607e-03, 0.6416944284294319e+00, 0.1160893767057166e+00, 0.1898168539265290e-03, 0.6697926391731260e+00, 0.1450378826743251e+00, 0.1902779940661772e-03, 0.6147594390585488e+00, 0.2904957622341456e-01, 0.1890125641731815e-03, 0.6455390026356783e+00, 0.5823809152617197e-01, 0.1899434637795751e-03, 0.6747258588365477e+00, 0.8740384899884715e-01, 0.1904520856831751e-03, 0.6772135750395347e+00, 0.2919946135808105e-01, 0.1905534498734563e-03 };
  //rVecPt.resize( 17430 );
  //rVecWt.resize( 5810 );

  _Lebedev_Oh_Mesh( nk, pa, rVecPt, rVecWt );

}






/*

  =========================================================================

 */







void ccdl::Backend::_Lebedev_Oh_Mesh
( int const * const nk, 
  double const * const pPa,
  std::tr1::array<double,3> * pPt,
  double * pWt )
{
  using ccdl::SQRT3;
  using ccdl::SQRT2;
  using ccdl::Backend::_8pt_helper;

  double const c08pt = 1.0 / SQRT3;
  double const c12pt = 1.0 / SQRT2;

  //double * pt;
  double * wt;
  double const * pa;

  pa = pPa-1;
  //pt = pPt-1;
  wt = pWt-1;

  if ( nk[1] > 0 )
    {
      pPt[0][0] =  1.0;
      pPt[0][1] =  0.0;
      pPt[0][2] =  0.0;

      pPt[1][0] = -1.0;
      pPt[1][1] =  0.0;
      pPt[1][2] =  0.0;

      pPt[2][0] =  0.0;
      pPt[2][1] =  1.0;
      pPt[2][2] =  0.0;

      pPt[3][0] =  0.0;
      pPt[3][1] = -1.0;
      pPt[3][2] =  0.0;

      pPt[4][0] =  0.0;
      pPt[4][1] =  0.0;
      pPt[4][2] =  1.0;

      pPt[5][0] =  0.0;
      pPt[5][1] =  0.0;
      pPt[5][2] = -1.0;
      pPt += 6;

      ++pa;
      for ( int i=0; i<6; ++i )
	*(++wt) = *pa;
    };

  if ( nk[2] > 0 )
    {
      _8pt_helper(pPt,c08pt,c08pt,c08pt);
      pPt += 8;

      ++pa;
      for ( int i=0; i<8; ++i )
	*(++wt) = *pa;
    };


  if ( nk[3] > 0 )
    {
      pPt[0][0] =  0.;
      pPt[0][1] =  c12pt;
      pPt[0][2] =  c12pt;

      pPt[1][0] =  0.;
      pPt[1][1] = -c12pt;
      pPt[1][2] =  c12pt;

      pPt[2][0] =  0.;
      pPt[2][1] =  c12pt;
      pPt[2][2] = -c12pt;

      pPt[3][0] =  0.;
      pPt[3][1] = -c12pt;
      pPt[3][2] = -c12pt;

      pPt[4][0] =  c12pt;
      pPt[4][1] =  0.;
      pPt[4][2] =  c12pt;

      pPt[5][0] = -c12pt;
      pPt[5][1] =  0.;
      pPt[5][2] =  c12pt;

      pPt[6][0] =  c12pt;
      pPt[6][1] =  0.;
      pPt[6][2] = -c12pt;

      pPt[7][0] = -c12pt;
      pPt[7][1] =  0.;
      pPt[7][2] = -c12pt;

      pPt[8][0] =  c12pt;
      pPt[8][1] =  c12pt;
      pPt[8][2] =  0.;

      pPt[9][0] = -c12pt;
      pPt[9][1] =  c12pt;
      pPt[9][2] =  0.;

      pPt[10][0] =  c12pt;
      pPt[10][1] = -c12pt;
      pPt[10][2] =  0.;

      pPt[11][0] = -c12pt;
      pPt[11][1] = -c12pt;
      pPt[11][2] =  0.;

      pPt += 12;
      
      ++pa;
      for ( int i=0; i<12; ++i )
	*(++wt) = *pa;
    };


  for ( int jj=0; jj<nk[4]; ++jj )
    {
      double const rr = *(++pa);
      double const ss = std::sqrt( 1.0 - 2.0 * rr * rr );

      _8pt_helper(pPt,rr,rr,ss);
      pPt += 8;
      _8pt_helper(pPt,rr,ss,rr);
      pPt += 8;
      _8pt_helper(pPt,ss,rr,rr);
      pPt += 8;
      
      ++pa;
      for ( int i=0; i<24; ++i )
	*(++wt) = *pa;
    };


  for ( int jj=0; jj<nk[5]; ++jj )
    {
      double const rr = *(++pa);
      double const ss = std::sqrt( 1.0 - rr * rr );

      pPt[0][0] =  rr;
      pPt[0][1] =  ss;
      pPt[0][2] =  0.;

      pPt[1][0] = -rr;
      pPt[1][1] =  ss;
      pPt[1][2] =  0.;

      pPt[2][0] =  rr;
      pPt[2][1] = -ss;
      pPt[2][2] =  0.;

      pPt[3][0] = -rr;
      pPt[3][1] = -ss;
      pPt[3][2] =  0.;

      pPt[4][0] =  ss;
      pPt[4][1] =  rr;
      pPt[4][2] =  0.;

      pPt[5][0] = -ss;
      pPt[5][1] =  rr;
      pPt[5][2] =  0.;

      pPt[6][0] =  ss;
      pPt[6][1] = -rr;
      pPt[6][2] =  0.;

      pPt[7][0] = -ss;
      pPt[7][1] = -rr;
      pPt[7][2] =  0.;

      pPt[8][0] =  rr;
      pPt[8][1] =  0.;
      pPt[8][2] =  ss;

      pPt[9][0] = -rr;
      pPt[9][1] =  0.;
      pPt[9][2] =  ss;

      pPt[10][0] =  rr;
      pPt[10][1] =  0.;
      pPt[10][2] = -ss;

      pPt[11][0] = -rr;
      pPt[11][1] =  0.;
      pPt[11][2] = -ss;

      pPt[12][0] =  ss;
      pPt[12][1] =  0.;
      pPt[12][2] =  rr;

      pPt[13][0] = -ss;
      pPt[13][1] =  0.;
      pPt[13][2] =  rr;

      pPt[14][0] =  ss;
      pPt[14][1] =  0.;
      pPt[14][2] = -rr;

      pPt[15][0] = -ss;
      pPt[15][1] =  0.;
      pPt[15][2] = -rr;

      pPt[16][0] =  0.;
      pPt[16][1] =  rr;
      pPt[16][2] =  ss;

      pPt[17][0] =  0.;
      pPt[17][1] = -rr;
      pPt[17][2] =  ss;

      pPt[18][0] =  0.;
      pPt[18][1] =  rr;
      pPt[18][2] = -ss;

      pPt[19][0] =  0.;
      pPt[19][1] = -rr;
      pPt[19][2] = -ss;

      pPt[20][0] =  0.;
      pPt[20][1] =  ss;
      pPt[20][2] =  rr;

      pPt[21][0] =  0.;
      pPt[21][1] = -ss;
      pPt[21][2] =  rr;

      pPt[22][0] =  0.;
      pPt[22][1] =  ss;
      pPt[22][2] = -rr;

      pPt[23][0] =  0.;
      pPt[23][1] = -ss;
      pPt[23][2] = -rr;

      pPt += 24;

      ++pa;
      for ( int i=0; i<24; ++i )
	*(++wt) = *pa;
    };


  for ( int jj=0; jj<nk[6]; ++jj )
    {
      double const rr = *(++pa);
      double const ss = *(++pa);
      double const tt = std::sqrt( 1.0 - rr*rr - ss*ss );

      _8pt_helper(pPt,rr,ss,tt);
      pPt += 8;
      _8pt_helper(pPt,rr,tt,ss);
      pPt += 8;
      _8pt_helper(pPt,ss,rr,tt);
      pPt += 8;
      _8pt_helper(pPt,ss,tt,rr);
      pPt += 8;
      _8pt_helper(pPt,tt,rr,ss);
      pPt += 8;
      _8pt_helper(pPt,tt,ss,rr);
      pPt += 8;

      ++pa;
      for ( int i=0; i<48; ++i )
	*(++wt) = *pa;
    };

}


void ccdl::Backend::_8pt_helper
( std::tr1::array<double,3> * pt, 
  double const rr, 
  double const ss, 
  double const tt )
{
  pPt[0][0] =  rr;
  pPt[0][1] =  ss;
  pPt[0][2] =  tt;
  
  pPt[1][0] = -rr;
  pPt[1][1] =  ss;
  pPt[1][2] =  tt;
  
  pPt[2][0] =  rr;
  pPt[2][1] = -ss;
  pPt[2][2] =  tt;
  
  pPt[3][0] = -rr;
  pPt[3][1] = -ss;
  pPt[3][2] =  tt;
  
  pPt[4][0] =  rr;
  pPt[4][1] =  ss;
  pPt[4][2] = -tt;
  
  pPt[5][0] = -rr;
  pPt[5][1] =  ss;
  pPt[5][2] = -tt;
  
  pPt[6][0] =  rr;
  pPt[6][1] = -ss;
  pPt[6][2] = -tt;
  
  pPt[7][0] = -rr;
  pPt[7][1] = -ss;
  pPt[7][2] = -tt;
}




double ccdl::LebedevCosmoGaussianZetaScaleFactor( int const iRule )
{
  double scl=4.90792902522;
  switch(iRule)
    {
    case 6:
      scl = 4.84566077868;
      break;
    case 14:
      scl = 4.86458714334;
      break;
      
    case 26:
      scl = 4.85478226219;
      break;

    case 38:
      scl = 4.90105812685;
      break;

    case 50:
      scl = 4.89250673295;
      break;

    case 86:
      scl = 4.89741372580;
      break;

    case 110:
      scl = 4.90101060987;
      break;

    case 146:
      scl = 4.89825187392;
      break;

    case 170:
      scl = 4.90685517725;
      break;

    case 194:
      scl = 4.90337644248;
      break;

    case 302:
      scl = 4.90498088169;
      break;

    case 350:
      scl = 4.86879474832;
      break;

    case 434:
      scl = 4.90567349080;
      break;

    case 590:
      scl = 4.90624071359;
      break;

    case 770:
      scl = 4.90656435779;
      break;

    case 974:
      scl = 4.90685167998;
      break;

    case 1202:
      scl = 4.90704098216;
      break;

    case 1454:
      scl = 4.90721023869;
      break;

    case 1730:
      scl = 4.90733270691;
      break;

    case 2030:
      scl = 4.90744499142;
      break;

    case 2354:
      scl = 4.90753082825;
      break;

    case 2702:
      scl = 4.90760972766;
      break;

    case 3074:
      scl = 4.90767282394;
      break;

    case 3470:
      scl = 4.90773141371;
      break;

    case 3890:
      scl = 4.90777965981;
      break;

    case 4334:
      scl = 4.90782469526;
      break;

    case 4802:
      scl = 4.90749125553;
      break;

    case 5294:
      scl = 4.90762073452;
      break;

    case 5810:
      scl = 4.90792902522;
      break;

    default:

      std::cerr << "Invalid angular quad rule in ccdl::LebedevCosmoGaussianZetaScaleFactor" << std::endl;
      abort();
      break;

    }

  return scl;
}


#endif

