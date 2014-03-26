#define _USE_MATH_DEFINES
#include <cmath>
#undef _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "CubicSpline.hpp"

double const PI = M_PI;
double const FOUR_PI = 4.*PI;
double const Y00 = std::sqrt( 1. / FOUR_PI );

#define EE(w,p)  std::setw(w) << std::setprecision(p) << std::scientific
#define FF(w,p)  std::setw(w) << std::setprecision(p) << std::fixed
#define WW(w)    std::setw(w)


int main()
{
  int const N = 10;
  std::vector<double> X(N,0.),Y(N,0.);

  //double zeta = 1.;

  for ( int i=0; i<N; ++i )
    {
      double x = (i+1)*(10./N);
      X[i] = x;
      Y[i] = 0.5 * x*x;
      //Y[i] = std::pow(zeta/PI,1.5) * std::exp( -zeta * x*x );
    };
  //double lder = 0.;
  //double rder = Y.back() * (-2.*zeta*X.back());

  double lder = X[0];
  double rder = X.back();

  ccdl::CubicSpline spl( X.size(), X.data(), Y.data(), lder, rder );


  for ( int i=0; i<N; ++i )
    std::cout << EE(24,15) << spl.GetX(i)
	      << EE(24,15) << spl.GetY(i)
	      << "\n";

  std::cout << "\n";
  for ( int i=0; i<(N+2)*20; ++i )
    {
      double x = -1. + i * 12./((N+2)*20-1);
      std::cout << EE(24,15) << x
		<< EE(24,15) << spl.GetValue(x)
		<< "\n";
    };

  std::cout << EE(24,15) << spl.Integrate<0>(1.,2.) << "\n";
  std::cout << EE(24,15) << std::pow(2.,3)/6. - std::pow(1.,3)/6. << "\n";

  std::cout << EE(24,15) << spl.Integrate<0>(11.,12.) << "\n";
  std::cout << EE(24,15) << std::pow(12.,3)/6. - std::pow(11.,3)/6. << "\n";

  return 0;
}
