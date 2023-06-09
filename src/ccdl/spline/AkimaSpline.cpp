#include "AkimaSpline.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>


ccdl::AkimaSpline::AkimaSpline
( int const n, double const * xs, double const * ys )
  : mX( xs, xs+n ),
    mY( ys, ys+n ),
    mB( n, 0. ),
    mC( n, 0. ),
    mD( n, 0. )
{

  //aspline(n,mX.data(),mY.data(),mB.data(),mC.data(),mD.data());
  //return;
  
  //
  // https://github.com/scipy/scipy/blob/v1.10.1/scipy/interpolate/_cubic.py#L367-L465
  //

  if ( n == 1 )
    {
      std::cerr << "AkimaSpline requires at least 2 points; received "
		<< n << std::endl;
      std::exit(EXIT_FAILURE);
    }
  else if ( n == 2 )
    {
      double dx = xs[1]-xs[0];
      if ( dx <= 0 )
	{
	  std::cerr << "AkimaSpline::AkimaSpline xs must be ascending"
		    << std::endl;
	  std::exit(EXIT_FAILURE);
	}
      else
	{
	  mB[0] = (ys[1]-ys[0]) / ( xs[1]-xs[0] );
	}
    }
  else
    {
      int n3 = n+3;
      int n2 = n+2;
      std::vector<double> m(n3,0);
      for ( int i=0; i<n-1; ++i )
	{
	  double dx = xs[i+1]-xs[i];
	  double dy = ys[i+1]-ys[i];
	  
	  if ( dx <= 0 )
	    {
	      std::cerr << "AkimaSpline::AkimaSpline xs must be ascending"
			<< std::endl;
	      std::exit(EXIT_FAILURE);
	    }
	  
	  m[i+2] = dy/dx;
	}
      m[1] = 2. * m[2] - m[3];
      m[0] = 2. * m[1] - m[2];
      m[n3-2] = 2. * m[n3-3] - m[n3-4];
      m[n3-1] = 2. * m[n3-2] - m[n3-3];

      for ( int i=0; i<n; ++i )
	{
	  mB[i] = 0.5 * ( m[3+i] + m[i] );
	}

      std::vector<double> dm(n2,0);
      for ( int i=0; i<n2; ++i )
	{
	  dm[i] = std::abs( m[i+1]-m[i] );
	}
      std::vector<double> f12(n,0);
      double f12max = -1.e+60;
      for ( int i=0; i<n; ++i )
	{
	  f12[i] = dm[2+i] + dm[i];
	  f12max = std::max(f12max,f12[i]);
	}
      std::vector<int> ind;
      for ( int i=0; i<n; ++i )
	{
	  if ( f12[i] > 1.e-9 * f12max )
	    {
	      ind.push_back( i );
	    }
	}

      for ( int ii=0, ni = ind.size(); ii<ni; ++ii )
	{
	  int i = ind[ii];
	  mB[i] = (dm[2+i] * m[i+1] + dm[i] * m[i+2])/f12[i];
	}

      for ( int i=0; i<n-1; ++i )
	{
	  double dx = xs[i+1]-xs[i];
	  double T = (mB[i] + mB[i+1] - 2 * m[i+2])/dx;
	  mD[i] = T / dx;
	  mC[i] = (m[i+2] - mB[i])/dx - T;
	}
      mB[n-1] = 0;
    }

  // for ( int i=0; i<n; ++i )
  //   {
  //     std::printf("%9.3f %9.3f %9.3f %9.3f %9.3f\n",mX[i],mY[i],mB[i],mC[i],mD[i]);
  //   };
}




int ccdl::AkimaSpline::FindIdx( double const x ) const
{
  int n = mX.size();
  int ind = 0;
  
  if ( x < mX[0] )
    {
      ind = 0;
    }
  else if ( x > mX[n-1] )
    {
      ind = n-1;
    }
  else
    {
      int low = 0;
      int high = n-1;
      while ( low <= high )
	{
	  ind = (low+high)/2;
	  if ( x < mX[ind] )
	    {
	      high = ind-1;
	    }
	  else if ( x > mX[ind+1] )
	    {
	      low = ind + 1;
	    }
	  else
	    {
	      break;
	    }
	}
    }
  return ind;
}


double ccdl::AkimaSpline::GetValue( double const x ) const
{
  double v = 0;
  if ( x < mX[0] )
    {
      v = mY[0];
    }
  else if ( x > mX.back() )
    {
      v = mY.back();
    }
  else
    {
      int ifound = FindIdx(x);
      double dx = x - mX[ifound];
      //std::printf("%.3f %2i %.3f\n",x,ifound,dx);
      v = mY[ifound] + dx*(mB[ifound] + dx*(mC[ifound] + dx*mD[ifound]));
    }
  return v;
}

