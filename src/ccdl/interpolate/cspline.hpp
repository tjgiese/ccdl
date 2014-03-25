#ifndef _cspline_hpp_
#define _cspline_hpp_

#include <cmath>
#include <algorithm>
#include <vector>


namespace ccdl
{

  void cspline_transform
  ( double const * x, 
    std::vector<double> & y, // y values on input; spldata on output
    double const lder, 
    double const rder,
    double * work_5n );

  int cspline_khi
  ( double const x, 
    int const n, 
    double const * X );

  double cspline_value
  ( double const x,
    int const n, 
    double const * X, 
    double const * spldata );

  double cspline_value_and_deriv
  ( double const x,
    double & dydx,
    int const n, 
    double const * X, 
    double const * spldata );

  double cspline_value_and_secondderiv
  ( double const x,
    double & dydx,
    double & d2ydx2,
    int const n, 
    double const * X, 
    double const * spldata );

  double cspline_value_given_khi
  ( double const x,
    int k, 
    double const * X, 
    double const * spldata );

  double cspline_value_and_deriv_given_khi
  ( double const x,
    double & dydx,
    int k, 
    double const * X, 
    double const * spldata );

  double cspline_value_and_secondderiv_given_khi
  ( double const x,
    double & dydx,
    double & d2ydx2,
    int k, 
    double const * X, 
    double const * spldata );

  template <int XPOW>
  double cspline_integral
  ( double const a, 
    double const b,
    int const n, 
    double const * X, 
    double const * spldata );

  template <int L>
  void cspline_lmpotential
  ( int const n, 
    double const * X, 
    std::vector<double> const * spldata,
    std::vector<double> * pots );


}





inline int ccdl::cspline_khi( double const x, int const n, double const * X )
{
  return std::distance( X, std::lower_bound( X+1, X+n-1, x ) );
}


inline double ccdl::cspline_value_given_khi
( double const x,
  int k, 
  double const * X, 
  double const * spldata )
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


inline double ccdl::cspline_value
( double const x,
  int const n, 
  double const * X, 
  double const * spldata )
{
  return ccdl::cspline_value_given_khi
    (x,
     ccdl::cspline_khi(x,n,X),
     X,
     spldata);
}






inline double ccdl::cspline_value_and_deriv
( double const x,
  double & dydx,
  int const n, 
  double const * X, 
  double const * spldata )
{
  return ccdl::cspline_value_and_deriv_given_khi
    (x,
     dydx,
     ccdl::cspline_khi(x,n,X),
     X,
     spldata);
}





inline double ccdl::cspline_value_and_secondderiv
( double const x,
  double & dydx,
  double & d2ydx2,
  int const n, double const * X, double const * spldata )
{
  return ccdl::cspline_value_and_secondderiv_given_khi
    (x,
     dydx,
     d2ydx2,
     ccdl::cspline_khi(x,n,X),
     X,
     spldata);
}



template <int XPOW>
double ccdl::cspline_integral
( double const a, double const b,
  int const n, double const * X, double const * spldata )
{
  if ( std::abs(a-b) < 1.e-10 ) return 0.;

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





template <int L>
void ccdl::cspline_lmpotential
( int const n, double const * X, 
  std::vector<double> const * spldata,
  std::vector<double> * pots )
{
  for ( int m=0; m<2*L+1; ++m )
    pots[m].resize(2*n);

  int const nm1 = n-1;

  double ai,bi;

  bi = X[0];

  double am[4],ap[4],apowm;
  double bm[4] = {0.,0.,0.,0.};
  double bp[4] = {0.,0.,0.,0.};
  double xa[3],xb[3];

  xb[0] = bi;
  xb[1] = xb[0]*bi;
  xb[2] = xb[1]*bi;


  apowm = 0.;

  if ( bi > 0. )
    {
      apowm = std::pow(bi,-(L+1));

      bm[0] = std::pow(bi,1-L);
      bm[1] = bm[0]*bi;
      bm[2] = bm[1]*bi;
      bm[3] = bm[2]*bi;

      bp[0] = std::pow(bi,2+L);
      bp[1] = bp[0]*bi;
      bp[2] = bp[1]*bi;
      bp[3] = bp[2]*bi;

      if ( (1-L) == -1 )
	{
	  bm[0] = std::log( bi );
	  bm[1] /= ((1-L)+2.);
	  bm[2] /= ((1-L)+3.);
	  bm[3] /= ((1-L)+4.);
	}
      else if ( (1-L) == -2 )
	{
	  bm[0] /= ((1-L)+1.);
	  bm[1] = std::log( bi );
	  bm[2] /= ((1-L)+3.);
	  bm[3] /= ((1-L)+4.);
	}
      else if ( (1-L) == -3 )
	{
	  bm[0] /= ((1-L)+1.);
	  bm[1] /= ((1-L)+2.);
	  bm[2] = std::log( bi );
	  bm[3] /= ((1-L)+4.);
	}
      else if ( (1-L) == -4 )
	{
	  bm[0] /= ((1-L)+1.);
	  bm[1] /= ((1-L)+2.);
	  bm[2] /= ((1-L)+3.);
	  bm[3] = std::log( bi );
	}
      else
	{
	  bm[0] /= ((1-L)+1.);
	  bm[1] /= ((1-L)+2.);
	  bm[2] /= ((1-L)+3.);
	  bm[3] /= ((1-L)+4.);
	};


      if ( (2+L) == -1 )
	{
	  bp[0] = std::log( bi );
	  bp[1] /= ((2+L)+2.);
	  bp[2] /= ((2+L)+3.);
	  bp[3] /= ((2+L)+4.);
	}
      else if ( (2+L) == -2 )
	{
	  bp[0] /= ((2+L)+1.);
	  bp[1] = std::log( bi );
	  bp[2] /= ((2+L)+3.);
	  bp[3] /= ((2+L)+4.);
	}
      else if ( (2+L) == -3 )
	{
	  bp[0] /= ((2+L)+1.);
	  bp[1] /= ((2+L)+2.);
	  bp[2] = std::log( bi );
	  bp[3] /= ((2+L)+4.);
	}
      else if ( (2+L) == -4 )
	{
	  bp[0] /= ((2+L)+1.);
	  bp[1] /= ((2+L)+2.);
	  bp[2] /= ((2+L)+3.);
	  bp[3] = std::log( bi );
	}
      else
	{
	  bp[0] /= ((2+L)+1.);
	  bp[1] /= ((2+L)+2.);
	  bp[2] /= ((2+L)+3.);
	  bp[3] /= ((2+L)+4.);
	};

    };
      
  for ( int j=0; j<nm1; ++j )
    {
      ai = bi;
      bi = X[j+1]; 
      double del = bi-ai;
      double del2 = del*del;

      xa[0] = xb[0];
      xa[1] = xb[1];
      xa[2] = xb[2];

      xb[0] = bi;
      xb[1] = xb[0] * bi;
      xb[2] = xb[1] * bi;

      am[0] = bm[0];
      am[1] = bm[1];
      am[2] = bm[2];
      am[3] = bm[3];
      ap[0] = bp[0];
      ap[1] = bp[1];
      ap[2] = bp[2];
      ap[3] = bp[3];

      bm[0] = std::pow(bi,(1-L)+1);
      bm[1] = bm[0]*bi;
      bm[2] = bm[1]*bi;
      bm[3] = bm[2]*bi;

      bp[0] = std::pow(bi,(2+L)+1);
      bp[1] = bp[0]*bi;
      bp[2] = bp[1]*bi;
      bp[3] = bp[2]*bi;


      if ( (1-L) == -1 )
	{
	  bm[0] = std::log( bi );
	  bm[1] /= ((1-L)+2.);
	  bm[2] /= ((1-L)+3.);
	  bm[3] /= ((1-L)+4.);
	}
      else if ( (1-L) == -2 )
	{
	  bm[0] /= ((1-L)+1.);
	  bm[1] = std::log( bi );
	  bm[2] /= ((1-L)+3.);
	  bm[3] /= ((1-L)+4.);
	}
      else if ( (1-L) == -3 )
	{
	  bm[0] /= ((1-L)+1.);
	  bm[1] /= ((1-L)+2.);
	  bm[2] = std::log( bi );
	  bm[3] /= ((1-L)+4.);
	}
      else if ( (1-L) == -4 )
	{
	  bm[0] /= ((1-L)+1.);
	  bm[1] /= ((1-L)+2.);
	  bm[2] /= ((1-L)+3.);
	  bm[3] = std::log( bi );
	}
      else
	{
	  bm[0] /= ((1-L)+1.);
	  bm[1] /= ((1-L)+2.);
	  bm[2] /= ((1-L)+3.);
	  bm[3] /= ((1-L)+4.);
	};

      
      if ( (2+L) == -1 )
	{
	  bp[0] = std::log( bi );
	  bp[1] /= ((2+L)+2.);
	  bp[2] /= ((2+L)+3.);
	  bp[3] /= ((2+L)+4.);
	}
      else if ( (2+L) == -2 )
	{
	  bp[0] /= ((2+L)+1.);
	  bp[1] = std::log( bi );
	  bp[2] /= ((2+L)+3.);
	  bp[3] /= ((2+L)+4.);
	}
      else if ( (2+L) == -3 )
	{
	  bp[0] /= ((2+L)+1.);
	  bp[1] /= ((2+L)+2.);
	  bp[2] = std::log( bi );
	  bp[3] /= ((2+L)+4.);
	}
      else if ( (2+L) == -4 )
	{
	  bp[0] /= ((2+L)+1.);
	  bp[1] /= ((2+L)+2.);
	  bp[2] /= ((2+L)+3.);
	  bp[3] = std::log( bi );
	}
      else
	{
	  bp[0] /= ((2+L)+1.);
	  bp[1] /= ((2+L)+2.);
	  bp[2] /= ((2+L)+3.);
	  bp[3] /= ((2+L)+4.);
	};


      double IntAm = am[1]-bm[1]+bi*(bm[0]-am[0]);
      double IntBm = bm[1]-am[1]+ai*(am[0]-bm[0]);
      
      double IntA3m = 
	( ( xb[2] * bm[0] - 3. * xb[1] * bm[1]
	    + 3. * xb[0]  * bm[2] - bm[3] )   -
	  ( xb[2] * am[0] - 3. * xb[1] * am[1]
	    + 3. * xb[0] * am[2] - am[3] ) );
      
      double IntB3m = 
	( ( xa[2] * am[0] - 3. * xa[1] * am[1]
	    + 3. * xa[0]  * am[2] - am[3] )  -
	  ( xa[2] * bm[0] - 3. * xa[1] * bm[1]
	    + 3. * xa[0] * bm[2] - bm[3] ) );

      double IntCm = (IntA3m - IntAm * del2 ) / 6.0;
      double IntDm = (IntB3m - IntBm * del2 ) / 6.0;


      double IntAp = ap[1]-bp[1]+bi*(bp[0]-ap[0]);
      double IntBp = bp[1]-ap[1]+ai*(ap[0]-bp[0]);
      
      double IntA3p = 
	( ( xb[2] * bp[0] - 3. * xb[1] * bp[1]
	    + 3. * xb[0]  * bp[2] - bp[3] )   -
	  ( xb[2] * ap[0] - 3. * xb[1] * ap[1]
	    + 3. * xb[0] * ap[2] - ap[3] ) );
      
      double IntB3p = 
	( ( xa[2] * ap[0] - 3. * xa[1] * ap[1]
	    + 3. * xa[0]  * ap[2] - ap[3] )  -
	  ( xa[2] * bp[0] - 3. * xa[1] * bp[1]
	    + 3. * xa[0] * bp[2] - bp[3] ) );

      double IntCp = (IntA3p - IntAp * del2 ) / 6.0;
      double IntDp = (IntB3p - IntBp * del2 ) / 6.0;

      double d[4];

      for ( int m=0; m<2*L+1; ++m )
	{
	  int tj = 2*j;
	  d[0] = spldata[m][  tj];
	  d[1] = spldata[m][++tj];
	  d[2] = spldata[m][++tj];
	  d[3] = spldata[m][++tj];

	  pots[m][j]   = (IntAm * d[0] + IntBm * d[2] + IntCm * d[1] + IntDm * d[3])/del;

	  pots[m][j+1+n] = pots[m][j+n] 
	    + (IntAp * d[0] + IntBp * d[2] + IntCp * d[1] + IntDp * d[3])/del;

	  pots[m][j+n] *= apowm;

	};

      apowm = std::pow(bi,-(L+1));

    }
  for ( int m=0; m<2*L+1; ++m )
    pots[m][2*n-1] *= apowm;

  for ( int j=n-1; j-->1; )
    {
      bi = X[j];
      ai = std::pow(bi,L);
      for ( int m=0; m<2*L+1; ++m )
	{
	  pots[m][j-1] += pots[m][j];
	  pots[m][j] *= ai;
	};
    };

  if ( L > 0 )
    {
      ai = std::pow(X[0],L);
      for ( int m=0; m<2*L+1; ++m )
	pots[m][0] *= ai;
    };

  double const FOUR_PI = 4. * 3.141592653589793238462643383279502884197;
  double fact = FOUR_PI / ( 2*L+1 );
  for ( int m=0; m<2*L+1; ++m )
    {
      for ( int j=0; j<n; ++j )
	{
	  pots[m][j] += pots[m][j+n];
	  pots[m][j] *= fact;
	};
      pots[m].resize(n);
    };

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






#endif
