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
  int const N = 200;
  std::vector<double> X(N,0.),Spl(N,0.),work(5*N);

  double zeta = 0.75;

  for ( int i=0; i<N; ++i )
    {
      double x = (i)*(6./N);
      X[i] = x;
      Spl[i] = std::pow(zeta/PI,1.5) * std::exp( -zeta * x*x );
    };
  double lder = 0.;
  double rder = Spl.back() * (-2.*zeta*X.back());

  ccdl::CubicSpline cspl( X.size(), X.data(), Spl.data(), lder, rder );
  ccdl::cspline_transform( X.data(), Spl, lder, rder, work.data() );

#define L 0
  ccdl::CubicSpline pspl(cspl);
  pspl.ConvertToPotential();


#define DFMT FMTE(20,10)
  for ( int i=0; i<N; ++i )
    {
      /*
      double mint = cspline_integral<1-L>(X[i],X.back(),N,X.data(),Spl.data());
      double pint = cspline_integral<2+L>(X[0],X[i],N,X.data(),Spl.data());

      if ( L > 0 )
	mint *= std::pow(X[i],L);
      if ( X[i] > 0. )
	pint *= std::pow(X[i],-(L+1));
      */

      double pot = 2. * std::sqrt(zeta/PI);
      if ( X[i] > 0. )
	pot = erf(std::sqrt(zeta)*X[i])/X[i];

      std::cout << DFMT << X[i]
		<< DFMT << cspl.GetY(i)
	//<< DFMT << mint
		<< DFMT << pot
		<< DFMT << pspl.GetY(i)
	//<< DFMT << Sint[0][i+N]
		<< "\n";
    };
#undef DFMT

  return 0;
}
