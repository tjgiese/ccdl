#ifndef _cspline_hpp_
#define _cspline_hpp_

#include <cmath>
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


void cspline_transform
( double const * x, 
  std::vector<double> & y, // y values on input; spldata on output
  double const lder, double const rder,
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

inline int cspline_khi( double const x, int const n, double const * X )
{
  return std::distance( X, std::lower_bound( X+1, X+n-1, x ) );
}


inline double cspline_value_given_khi
( double const x,
  int k, double const * X, double const * spldata )
{
  double const x1 = X[k];
  double const x0 = X[--k];
  k *= 2;
  double const h = x1-x0;
  double const a = (x1-x)/h;
  double const b = (x-x0)/h;

  double y0 = spldata[  k];
  double d0 = spldata[++k];
  double y1 = spldata[++k];
  double d1 = spldata[++k];
  return a*y0+b*y1+((a*a*a-a)*d0+(b*b*b-b)*d1)*h*h/6.;
}


inline double cspline_value
( double const x,
  int const n, double const * X, double const * spldata )
{
  return cspline_value_given_khi
    (x,
     cspline_khi(x,n,X),
     X,
     spldata);
}


double cspline_value_and_deriv_given_khi
( double const x,
  double & dydx,
  int k, double const * X, double const * spldata )
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



inline double cspline_value_and_deriv
( double const x,
  double & dydx,
  int const n, double const * X, double const * spldata )
{
  return cspline_value_and_deriv_given_khi
    (x,
     dydx,
     cspline_khi(x,n,X),
     X,
     spldata);
}




double cspline_value_and_secondderiv_given_khi
( double const x,
  double & dydx,
  double & d2ydx2,
  int k, double const * X, double const * spldata )
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


inline double cspline_value_and_secondderiv
( double const x,
  double & dydx,
  double & d2ydx2,
  int const n, double const * X, double const * spldata )
{
  return cspline_value_and_secondderiv_given_khi
    (x,
     dydx,
     d2ydx2,
     cspline_khi(x,n,X),
     X,
     spldata);
}


template <int XPOW>
double cspline_integral
( double const a, double const b,
  int const n, double const * X, double const * spldata )
{
  double TotInteg = 0.0;
  int const nm1 = n-1;

  double ai,bi;
  double an1,an2,an3,an4;
  double bn1,bn2,bn3,bn4;

  int jlo = nm1;
  for ( int j=0; j<nm1; ++j )
    {
      if ( X[j+1] >= a )
	{
	  jlo = j;
	  break;
	};
    };

  if ( X[jlo] < a ) bi = a;

  double an[4];
  double bn[4] = {0.,0.,0.,0.};
  double xa[3],xb[3];

  xb[0] = X[jlo];
  xb[1] = xb[0]*X[jlo];
  xb[2] = xb[1]*X[jlo];

  if ( bi > 0. )
    {

      bn[0] = std::pow(bi,XPOW+1);
      bn[1] = bn[0]*bi;
      bn[2] = bn[1]*bi;
      bn[3] = bn[2]*bi;

      if ( XPOW == -1 )
	{
	  bn[0] = std::log( bi );
	  bn[1] /= (XPOW+2.);
	  bn[2] /= (XPOW+3.);
	  bn[3] /= (XPOW+4.);
	}
      else if ( XPOW == -2 )
	{
	  bn[0] /= (XPOW+1.);
	  bn[1] = std::log( bi );
	  bn[2] /= (XPOW+3.);
	  bn[3] /= (XPOW+4.);
	}
      else if ( XPOW == -3 )
	{
	  bn[0] /= (XPOW+1.);
	  bn[1] /= (XPOW+2.);
	  bn[2] = std::log( bi );
	  bn[3] /= (XPOW+4.);
	}
      else if ( XPOW == -4 )
	{
	  bn[0] /= (XPOW+1.);
	  bn[1] /= (XPOW+2.);
	  bn[2] /= (XPOW+3.);
	  bn[3] = std::log( bi );
	}
      else
	{
	  bn[0] /= (XPOW+1.);
	  bn[1] /= (XPOW+2.);
	  bn[2] /= (XPOW+3.);
	  bn[3] /= (XPOW+4.);
	};
    };
      
  for ( int j=jlo; j<nm1; ++j )
    {
      if ( X[j] > b ) break;

      double del = X[j+1] - X[j];
      double del2 = del*del;

      xa[0] = xb[0];
      xa[1] = xb[1];
      xa[2] = xb[2];

      xb[0] = X[j+1];
      xb[1] = xb[0] * X[j+1];
      xb[2] = xb[1] * X[j+1];

      ai = bi;
      an[0] = bn[0];
      an[1] = bn[1];
      an[2] = bn[2];
      an[3] = bn[3];

      bi = X[j+1]; 
      if ( X[j+1] > b ) bi = b;

      bn[0] = std::pow(bi,XPOW+1);
      bn[1] = bn[0]*bi;
      bn[2] = bn[1]*bi;
      bn[3] = bn[2]*bi;
      // bn[0] /= (XPOW+1);
      // bn[1] /= (XPOW+2);
      // bn[2] /= (XPOW+3);
      // bn[3] /= (XPOW+4);
      
      // if ( XPOW == -1 )
      // 	{
      // 	  bn[0] = std::log( bi );
      // 	}
      // else if ( XPOW == -2 )
      // 	{
      // 	  bn[1] = std::log( bi );
      // 	}
      // else if ( XPOW == -3 )
      // 	{
      // 	  bn[2] = std::log( bi );
      // 	}
      // else if ( XPOW == -4 )
      // 	{
      // 	  bn[3] = std::log( bi );
      // 	};


      if ( XPOW == -1 )
	{
	  bn[0] = std::log( bi );
	  bn[1] /= (XPOW+2.);
	  bn[2] /= (XPOW+3.);
	  bn[3] /= (XPOW+4.);
	}
      else if ( XPOW == -2 )
	{
	  bn[0] /= (XPOW+1.);
	  bn[1] = std::log( bi );
	  bn[2] /= (XPOW+3.);
	  bn[3] /= (XPOW+4.);
	}
      else if ( XPOW == -3 )
	{
	  bn[0] /= (XPOW+1.);
	  bn[1] /= (XPOW+2.);
	  bn[2] = std::log( bi );
	  bn[3] /= (XPOW+4.);
	}
      else if ( XPOW == -4 )
	{
	  bn[0] /= (XPOW+1.);
	  bn[1] /= (XPOW+2.);
	  bn[2] /= (XPOW+3.);
	  bn[3] = std::log( bi );
	}
      else
	{
	  bn[0] /= (XPOW+1.);
	  bn[1] /= (XPOW+2.);
	  bn[2] /= (XPOW+3.);
	  bn[3] /= (XPOW+4.);
	};

      double IntA = an[1]-bn[1]+X[j+1]*(bn[0]-an[0]);
      double IntB = bn[1]-an[1]+X[j  ]*(an[0]-bn[0]);
      
      double IntA3 = 
	( ( xb[2] * bn[0] - 3. * xb[1] * bn[1]
	  + 3. * xb[0]  * bn[2] - bn[3] )   -
	( xb[2] * an[0] - 3. * xb[1] * an[1]
	  + 3. * xb[0] * an[2] - an[3] ) );
      
      double IntB3 = 
	( ( xa[2] * an[0] - 3. * xa[1] * an[1]
	  + 3. * xa[0]  * an[2] - an[3] )  -
	( xa[2] * bn[0] - 3. * xa[1] * bn[1]
	  + 3. * xa[0] * bn[2] - bn[3] ) );

      double IntC = (IntA3 - IntA * del2 ) / 6.0;
      double IntD = (IntB3 - IntB * del2 ) / 6.0;

      int tj = 2*j;
      double d[4];
      d[0] = spldata[  tj];
      d[1] = spldata[++tj];
      d[2] = spldata[++tj];
      d[3] = spldata[++tj];

      double const IntervalInteg = 
	(IntA * d[0] + IntB * d[2] + IntC * d[1] + IntD * d[3])/del;
      
      TotInteg += IntervalInteg;
    };

  return TotInteg;
}




namespace test
{

  template <int XPOW>
  double Int1( double x, double d, double a1 )
  {
    double Int1_res = 0.;
    if ( a1 > 0. )
      {
	double a1n1,a1n2;
	if ( XPOW+1 == 0 )
	  {
	    a1n1 = std::log(a1);
	  }
	else
	  {
	    a1n1 = std::pow(a1,XPOW+1)/(XPOW+1.);
	  };
	if ( XPOW+2 == 0 )
	  {
	    a1n2 = std::log(a1);
	  }
	else
	  {
	    a1n2 = std::pow(a1,XPOW+2)/(XPOW+2.);
	  };
	Int1_res = (x * a1n1 - a1n2) / d;
      };
    return Int1_res;
  }

  template <int XPOW>
  double Int2( double x, double d, double a1 )
  {
    double Int2_res = 0.;
    if ( a1 > 0. )
      {
	double x2 = x*x;
	double x3 = x2*x;
	double a1n1,a1n2,a1n3,a1n4;
	if ( XPOW+1 == 0 )
	  {
	    a1n1 = std::log(a1);
	  }
	else
	  {
	    a1n1 = std::pow(a1,XPOW+1)/(XPOW+1.);
	  };
	if ( XPOW+2 == 0 )
	  {
	    a1n2 = std::log(a1);
	  }
	else
	  {
	    a1n2 = std::pow(a1,XPOW+2)/(XPOW+2.);
	  };
	if ( XPOW+3 == 0 )
	  {
	    a1n3 = std::log(a1);
	  }
	else
	  {
	    a1n3 = std::pow(a1,XPOW+3)/(XPOW+3.);
	  };

	if ( XPOW+4 == 0 )
	  {
	    a1n4 = std::log(a1);
	  }
	else
	  {
	    a1n4 = std::pow(a1,XPOW+4)/(XPOW+4.);
	  };
	Int2_res = ( x3 * a1n1 
		     - 3.0 * x2 * a1n2 
		     + 3.0 * x  * a1n3 
		     - a1n4 ) / (d*d*d);
      };
    return Int2_res;
  }




  template <int XPOW>
  double cspline_integral
  ( double const a, double const b,
    int const n, double const * X, double const * spldata )
  {
    double TotInteg = 0.0;
    int const nm1 = n-1;

    for ( int j=0; j<nm1; ++j )
      {
      
	if ( X[j+1] < a ) continue;
	if ( X[j] > b ) break;

	double bi = X[j+1];
	double ai = X[j];
	if ( X[j+1] > b ) bi = b;
	if ( X[j] < a ) ai = a;

	double del = X[j+1] - X[j];

	double IntA = Int1<XPOW>(X[j+1],del,bi) - Int1<XPOW>(X[j+1],del,ai);
	double IntB = Int1<XPOW>(X[j],del,ai) - Int1<XPOW>(X[j],del,bi);

	double IntA3 = Int2<XPOW>(X[j+1],del,bi) - Int2<XPOW>(X[j+1],del,ai);
	double IntB3 = Int2<XPOW>(X[j],del,ai) - Int2<XPOW>(X[j],del,bi);

	double IntC = (IntA3 - IntA) * del*del / 6.0;
	double IntD = (IntB3 - IntB) * del*del / 6.0;

	double IntervalInteg = IntA * spldata[2*j]  + IntB * spldata[2*(j+1)]
	  + IntC * spldata[1+2*j] + IntD * spldata[1+2*(j+1)];
	TotInteg += IntervalInteg;
      
      };
  
    return TotInteg;
  }

}





template <int XPOW>
double cspline_integral_bak
( double const a, double const b,
  int const n, double const * X, double const * spldata )
{
  double TotInteg = 0.0;
  int const nm1 = n-1;

  double ai,bi;
  double an1,an2,an3,an4;
  double bn1,bn2,bn3,bn4;

  int jlo = nm1;
  for ( int j=0; j<nm1; ++j )
    {
      if ( X[j+1] >= a )
	{
	  jlo = j;
	  break;
	};
    };

  if ( X[jlo] < a ) bi = a;

  double an[4];
  double bn[4];
  double xa[3],xb[3];

  xb[0] = X[jlo];
  xb[1] = xb[0]*X[jlo];
  xb[2] = xb[1]*X[jlo];

  bn[0] = std::pow(bi,XPOW+1);
  bn[1] = bn[0]*bi;
  bn[2] = bn[1]*bi;
  bn[3] = bn[2]*bi;
  bn[0] /= (XPOW+1);
  bn[1] /= (XPOW+2);
  bn[2] /= (XPOW+3);
  bn[3] /= (XPOW+4);

  if ( bi > 0. )
    {
      if ( XPOW == -1 )
	{
	  bn[0] = std::log( bi );
	}
      else if ( XPOW == -2 )
	{
	  bn[1] = std::log( bi );
	}
      else if ( XPOW == -3 )
	{
	  bn[2] = std::log( bi );
	}
      else if ( XPOW == -4 )
	{
	  bn[3] = std::log( bi );
	};
    };
      
  for ( int j=jlo; j<nm1; ++j )
    {
      if ( X[j] > b ) break;

      double del = X[j+1] - X[j];
      double del2 = del*del;
      double del3 = del2*del;

      xa[0] = xb[0];
      xa[1] = xb[1];
      xa[2] = xb[2];

      xb[0] = X[j+1];
      xb[1] = xb[0] * X[j+1];
      xb[2] = xb[1] * X[j+1];

      ai = bi;
      an[0] = bn[0];
      an[1] = bn[1];
      an[2] = bn[2];
      an[3] = bn[3];

      bi = X[j+1]; 
      if ( X[j+1] > b ) bi = b;

      bn[0] = std::pow(bi,XPOW+1);
      bn[1] = bn[0]*bi;
      bn[2] = bn[1]*bi;
      bn[3] = bn[2]*bi;
      bn[0] /= (XPOW+1);
      bn[1] /= (XPOW+2);
      bn[2] /= (XPOW+3);
      bn[3] /= (XPOW+4);
      
      if ( XPOW == -1 )
	{
	  bn[0] = std::log( bi );
	}
      else if ( XPOW == -2 )
	{
	  bn[1] = std::log( bi );
	}
      else if ( XPOW == -3 )
	{
	  bn[2] = std::log( bi );
	}
      else if ( XPOW == -4 )
	{
	  bn[3] = std::log( bi );
	};

      // double IntA = 
      // 	(X[j+1] * bn[0] - bn[1])  -
      // 	(X[j+1] * an[0] - an[1]);
      
      // double IntB = 
      // 	(X[j] * an[0] - an[1])  -
      // 	(X[j] * bn[0] - bn[1]);

      double IntA = an[1]-bn[1]+X[j+1]*(bn[0]-an[0]);
      double IntB = bn[1]-an[1]+X[j  ]*(an[0]-bn[0]);
      
      double IntA3 = 
	( ( xb[2] * bn[0] - 3. * xb[1] * bn[1]
	  + 3. * xb[0]  * bn[2] - bn[3] )   -
	( xb[2] * an[0] - 3. * xb[1] * an[1]
	  + 3. * xb[0] * an[2] - an[3] ) );
      
      double IntB3 = 
	( ( xa[2] * an[0] - 3. * xa[1] * an[1]
	  + 3. * xa[0]  * an[2] - an[3] )  -
	( xa[2] * bn[0] - 3. * xa[1] * bn[1]
	  + 3. * xa[0] * bn[2] - bn[3] ) );

      double IntC = (IntA3 - IntA * del2 ) / 6.0;
      double IntD = (IntB3 - IntB * del2 ) / 6.0;

      int tj = 2*j;
      double d[4];
      d[0] = spldata[  tj];
      d[1] = spldata[++tj];
      d[2] = spldata[++tj];
      d[3] = spldata[++tj];

      double const IntervalInteg = 
	(IntA * d[0] + IntB * d[2] + IntC * d[1] + IntD * d[3])/del;
      
      TotInteg += IntervalInteg;
    };

  return TotInteg;
}


#endif
