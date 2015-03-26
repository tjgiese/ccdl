#include "Quadrature.hpp"
#include "../constants.hpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cassert>

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
  ccdl::arth(3.,4.,m,rhs.data());
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




std::tr1::shared_ptr< std::vector< double > >
ccdl::AngularQuadDatabase::GetWts( int rule )
{
  SharedWt w; 
  iterator it = mDB.find( rule );
  if ( it == mDB.end() )
    {
      SharedWt wts( new WtType( rule, 0. ) );
      SharedCrd crd( new CrdType( rule ) );
      ccdl::LebedevRule( rule, crd->data(), wts->data() );
      mDB.insert( std::make_pair( rule, std::make_pair( wts, crd ) ) );
      w = wts;
    }
  else
    {
      w = it->second.first;
    };
  return w;
}

std::tr1::shared_ptr< std::vector< std::tr1::array< double,3 > > >
ccdl::AngularQuadDatabase::GetPts( int rule )
{
  SharedCrd c;
  iterator it = mDB.find( rule );
  if ( it == mDB.end() )
    {
      SharedWt wts( new WtType( rule, 0. ) );
      SharedCrd crd( new CrdType( rule ) );
      ccdl::LebedevRule( rule, crd->data(), wts->data() );
      mDB.insert( std::make_pair( rule, std::make_pair( wts, crd ) ) );
      c = crd;
    }
  else
    {
      c = it->second.second;
    };
  return c;
}




