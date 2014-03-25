#include "cspline.hpp"

#include <iostream>
#include <iomanip>
#include <cstdlib>


/*
void cspline_tridag
( int const n, 
  double const * a, double const * b, 
  double const * c, double const * r, 
  double * u, double * work )
{
  // work[n]
  double bet = b[0];
  u[0]=r[0]/bet;
  work[0] = 0.;
  for ( int j=1; j<n; ++j )
    {
      work[j]=c[j-1]/bet;
      bet = b[j]-a[j-1]*work[j];
      u[j]=(r[j]-a[j-1]*u[j-1])/bet;
    };
  for ( int j=n; j-->0; )
    u[j]=u[j]-work[j+1]*u[j+1];
}
*/


void cspline_tridag
( int const n, 
  double const * a, double const * b, 
  double const * c, double const * r, 
  double * spldata, double * work )
{
  // work[n]
  double bet = b[0];
  spldata[1]=r[0]/bet;
  work[0] = 0.;
  for ( int j=1; j<n; ++j )
    {
      work[j]=c[j-1]/bet;
      bet = b[j]-a[j-1]*work[j];
      spldata[1+j*2]=(r[j]-a[j-1]*spldata[1+(j-1)*2])/bet;
    };
  for ( int j=n; j-->0; )
    spldata[1+j*2]-=work[j+1]*spldata[1+(j+1)*2];
}


void ccdl::cspline_transform
( double const * x, 
  std::vector<double> & y, // y values on input; spldata on output
  double const lder, 
  double const rder,
  double * work_5n )
{
  int const n = y.size();
  int const nm1 = n-1;
  // work_5n[5*n]
  double * c = work_5n;
  double * r = c+n;
  double * a = r+n;
  double * b = a+n;
  double * g = b+n;


  c[0] = x[1]-x[0];
  double prev = 6. * ( y[1]-y[0] ) / c[0];
  r[0] = prev;
  for ( int i=1; i<nm1; ++i )
    {
      c[i] = x[i+1]-x[i];
      double cur = 6. * ( y[i+1]-y[i] ) / c[i];
      r[i] = cur-prev;
      prev = cur;
      a[i] = c[i-1];
      b[i] = 2.*( c[i]+a[i] );
    };
  b[0]   = 1.;
  b[nm1] = 1.;
  r[0]   = (3.0/(x[1]-x[0])) * ( (y[1]-y[0])/(x[1]-x[0]) - lder );
  c[0]   = 0.5;
  r[nm1] = (-3.0/(x[nm1]-x[nm1-1]))
    *( (y[nm1]-y[nm1-1])/(x[nm1]-x[nm1-1]) - rder );
  a[nm1] = 0.5;

  //for ( int i=0; i<n; ++i )
  //  spldata[i*2] = y[i];
  y.resize( 2*n );
  for ( int i=n; i-->1; )
    y[2*i] = y[i];

  cspline_tridag(n,a+1,b,c,r,y.data(),g);

  // for ( int i=0; i<n; ++i )
  //   {
  //     std::cout << std::scientific 
  // 		<< std::setw(14) << std::setprecision(5) << x[i]
  // 		<< std::scientific 
  // 		<< std::setw(14) << std::setprecision(5) << y[i]
  // 		<< std::scientific 
  // 		<< std::setw(14) << std::setprecision(5) << y2[i]
  // 		<< std::scientific 
  // 		<< std::setw(14) << std::setprecision(5) << a[i]
  // 		<< std::scientific 
  // 		<< std::setw(14) << std::setprecision(5) << b[i]
  // 		<< std::scientific 
  // 		<< std::setw(14) << std::setprecision(5) << c[i]
  // 		<< std::scientific 
  // 		<< std::setw(14) << std::setprecision(5) << r[i]
  // 		<< "\n";
  //   };
  // exit(0);
}






double ccdl::cspline_value_and_deriv_given_khi
( double const x,
  double & dydx,
  int k, 
  double const * X, 
  double const * spldata )
{
  double const x1 = X[k];
  double const x0 = X[--k];
  k *= 2;
  double const h = x1-x0;
  double const h16 = h/6.;
  double const h26 = h*h16;
  double const ooh = 1./h;
  double const a = (x1-x)*ooh;
  double const b = 1.-a; //(x-x0)*ooh;
  double const c = (a*a*a-a)*h26;
  double const d = (b*b*b-b)*h26;
  //double const dadx = -ooh;
  //double const dbdx =  ooh;
  double const dcdx = (1.-3.*a*a)*h16;
  double const dddx = (3.*b*b-1.)*h16;

  double y0 = spldata[  k];
  double d0 = spldata[++k];
  double y1 = spldata[++k];
  double d1 = spldata[++k];
  //dydx = dadx*y0+dbdx*y1+dcdx*d0+dddx*d1;
  dydx = ooh *( y1-y0 ) +dcdx*d0+dddx*d1;
  return    a*y0+   b*y1+   c*d0+   d*d1;
}




double ccdl::cspline_value_and_secondderiv_given_khi
( double const x,
  double & dydx,
  double & d2ydx2,
  int k, 
  double const * X, 
  double const * spldata )
{
  double const x1 = X[k];
  double const x0 = X[--k];
  k *= 2;
  double const h = x1-x0;
  double const h16 = h/6.;
  double const h26 = h*h16;
  double const ooh = 1./h;
  double const a = (x1-x)*ooh;
  double const b = 1.-a; //(x-x0)*ooh;
  double const c = (a*a*a-a)*h26;
  double const d = (b*b*b-b)*h26;
  //double const dadx = -ooh;
  //double const dbdx =  ooh;
  double const dcdx = (1.-3.*a*a)*h16;
  double const dddx = (3.*b*b-1.)*h16;

  double y0 = spldata[  k];
  double d0 = spldata[++k];
  double y1 = spldata[++k];
  double d1 = spldata[++k];
  //dydx = dadx*y0+dbdx*y1+dcdx*d0+dddx*d1;
  dydx = ooh *( y1-y0 ) +dcdx*d0+dddx*d1;
  d2ydx2 = a * d0 + b * d1;
  return    a*y0+   b*y1+   c*d0+   d*d1;
}


