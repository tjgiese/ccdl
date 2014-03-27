#define _USE_MATH_DEFINES
#include <cmath>
#undef _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>

#include "CubicSpline.hpp"

double const PI = M_PI;
double const FOUR_PI = 4.*PI;
double const Y00 = std::sqrt( 1. / FOUR_PI );

#define FMTE(w,p)  std::scientific << std::setw(w) << std::setprecision(p)
#define FMTF(w,p)  std::fixed << std::setw(w) << std::setprecision(p)
#define FMTW(w)  std::setw(w)


int main()
{
  int const N = 64;
  std::vector<double> X(N,0.),Spl(N,0.),work(5*N);

  double zeta = 0.75;

  for ( int i=0; i<N; ++i )
    {
      double x = (i)*(6./N);
      X[i] = x;
      Spl[i] = std::pow(zeta/PI,1.5) * std::exp( -zeta * x*x );
    };
  double lder = Spl[0] * (-2.*zeta*X[0]);
  double rder = Spl.back() * (-2.*zeta*X.back());

  //ccdl::CubicSpline cspl( X.size(), X.data(), Spl.data(), lder, rder );
  ccdl::CubicSpline cspl( X.size(), X.data(), Spl.data(), false, false );

  ccdl::cspline_transform( X.data(), Spl, lder, rder, work.data() );

  

#define L 0
  ccdl::CubicSpline pspl(cspl);
  pspl.RadialDensityToRadialPotential( 1.0 );
  //double mint = ccdl::cspline_integral<1-L>(0.,X.back(),N,X.data(),Spl.data());

#define DFMT FMTE(20,10)
  for ( int i=0; i<N; ++i )
    {
      
      double mint = ccdl::cspline_integral<1-L>(X[i],X.back(),N,X.data(),Spl.data());

      double pint = ccdl::cspline_integral<2+L>(0.,X[i],N,X.data(),Spl.data());


      if ( L > 0 )
       	mint *= std::pow(X[i],L);
      if ( X[i] > 0. )
       	pint *= std::pow(X[i],-(L+1));
      

      double pot = 2. * std::sqrt(zeta/PI);
      if ( X[i] > 0. )
  	pot = erf(std::sqrt(zeta)*X[i])/X[i];

      std::cout << DFMT << X[i]
  		<< DFMT << cspl.GetY(i)
  		<< DFMT << pot
  		<< DFMT << pspl.GetY(i)
	<< DFMT << FOUR_PI * ( mint + pint )
	<< DFMT << FOUR_PI * mint
	<< DFMT << FOUR_PI * pint
	//		<< DFMT << X[i] << DFMT << X.back()
  		<< "\n";
    };
  //std::cout << DFMT << FOUR_PI * cspl.Integrate<2>(0.,cspl.GetRightX()) << "\n";
  std::cout << "\n\n\n";


  
  for ( int i=0; i<20*8; ++i )
    {
      double x = i*(10./(20*8));
      double pot = 2. * std::sqrt(zeta/PI);
      if ( x > 0. )
	pot = erf(std::sqrt(zeta)*x)/x;

      double mint = ccdl::cspline_integral<1-L>(x,X.back(),N,X.data(),Spl.data());
      double pint = ccdl::cspline_integral<2+L>(0.,x,N,X.data(),Spl.data());

      if ( L > 0 )
      	mint *= std::pow(x,L);
      if ( x > 0. )
      	pint *= std::pow(x,-(L+1));


      std::cout << DFMT << x 
      		<< DFMT << cspl.GetValue(x) 
       		<< DFMT << pot
     		<< DFMT << pspl.GetValue(x) 
	<< DFMT << FOUR_PI * (mint + pint)
      		<< "\n";
    };
  

#undef DFMT




  return 0;
}
