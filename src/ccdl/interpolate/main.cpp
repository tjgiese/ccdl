#include <vector>
#include <iostream>
#include <iomanip>

#include "cspline.hpp"


double const PI = 3.14159265358979e+00;
double const FOUR_PI = 4.*PI;

int main()
{
  int const n = 150*4;
  std::vector<double> x(n,0.), y(n,0.), work(5*n,0.);
  double zeta = 1.;
  for ( int i=0; i<n; ++i )
    {
      x[i] = i * ( 5./(n-1) );
      y[i] = std::pow(zeta/PI,1.5) * std::exp( - zeta * x[i] * x[i]  );
  //     std::cout << std::fixed 
  //     		<< std::setw(15) << std::setprecision(7) << x[i]
  //     		<< std::scientific 
  //     		<< std::setw(23) << std::setprecision(15) << y[i]
  //     		<< "\n";
  //   };
  // std::cout << "\n";
    };

  double glo = -zeta * y[0] * x[0];
  double ghi = -zeta * y[n-1] * x[n-1];
  cspline_transform( x.data(), y, glo, ghi, work.data() );

  /*
  for ( int i=0; i<1000*n; ++i )
    {
      double X = -1. + i * ( 7./(1000.*n-1.) );
      double dYdX,d2YdX2;
      //double Y = cspline_value( X, n, x.data(), y.data() );
      //double Y = cspline_value_and_deriv( X, dYdX, n, x.data(), y.data() );
      double Y = cspline_value_and_secondderiv( X, dYdX, d2YdX2, n, x.data(), y.data() );

      std::cout << std::scientific
		<< std::setw(23) << std::setprecision(15) << X
		<< std::scientific 
		<< std::setw(23) << std::setprecision(15) << Y
		<< std::scientific 
		<< std::setw(23) << std::setprecision(15) << dYdX
		<< std::scientific 
		<< std::setw(23) << std::setprecision(15) << d2YdX2

		<< "\n";
    };
  */

 
  std::vector<double> a(1024);

  for ( int k=0; k<20; ++k )
    for ( int i=0; i<a.size(); ++i )
      {
	a[i] = cspline_integral<-1>(0.,6.,n,x.data(),y.data());
      };

  std::cout << std::scientific 
	    << std::setw(23) << std::setprecision(15)
	    << a[0]
	    << "\n";

  // 0.159155

  return 0;
}
