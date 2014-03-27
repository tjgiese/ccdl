#include "CubicSpline.hpp"


ccdl::CubicSpline::CubicSpline
( int const n, 
  double const * x, 
  double const * y, 
  double const lder, 
  double const rder )
  : mX(x,x+n),
    mY(y,y+n)
{
  std::vector<double> work(5*n);
  ccdl::cspline_transform(mX.data(),mY,lder,rder,work.data());
}



double GetLeftSlope( int const n, double const * X, double const * Y )
{
  double m;
  if ( n > 2 )
    {
      double x1 = 0.5 * (X[0]+X[1]);
      double x2 = 0.5 * (X[1]+X[2]);
      double m1 = (Y[0]-Y[1]) / (X[0]-X[1]);
      double m2 = (Y[1]-Y[2]) / (X[1]-X[2]);
      double mm = (m1-m2)/(x1-x2);
      double b  = m1 - mm * x1;
      m = mm * X[0] + b;
    }
  else
    {
      m = (Y[1]-Y[0]) / (X[1]-X[0]);
    };
  return m;
}

double GetRightSlope( int const n, double const * X, double const * Y )
{
  double m;
  if ( n > 2 )
    {
      double x1 = 0.5 * (X[n-1]+X[n-2]);
      double x2 = 0.5 * (X[n-2]+X[n-3]);
      double m1 = (Y[n-1]-Y[n-2]) / (X[n-1]-X[n-2]);
      double m2 = (Y[n-2]-Y[n-3]) / (X[n-2]-X[n-3]);
      double mm = (m1-m2)/(x1-x2);
      double b  = m1 - mm * x1;
      m = mm * X[n-1] + b;
    }
  else
    {
      m = (Y[n-2]-Y[n-1]) / (X[n-2]-X[n-1]);
    };
  return m;
}



    
ccdl::CubicSpline::CubicSpline
( int const n, 
  double const * x, 
  double const * y,
  bool const lzero,
  bool const rzero )
  : mX(x,x+n),
    mY(y,y+n)
{
  double lder = lzero ? 0. : GetLeftSlope(n,x,y);
  double rder = rzero ? 0. : GetRightSlope(n,x,y);
  std::vector<double> work(5*n);
  ccdl::cspline_transform(mX.data(),mY,lder,rder,work.data());
}


void ccdl::CubicSpline::Reset
( double const * y,
  double const lder, 
  double const rder )
{
  int const n = GetN();
  mY.resize( n );
  mY.assign( y, y+n );
  std::vector<double> work(5*n);
  ccdl::cspline_transform(mX.data(),mY,lder,rder,work.data());
}


void ccdl::CubicSpline::Reset
( int const n, 
  double const * x, 
  double const * y, 
  double const lder, 
  double const rder )
{
  mX.resize( n );
  mX.assign( x, x+n );
  Reset(y,lder,rder);
}


void ccdl::CubicSpline::Reset
( double const * y,
  bool const lzero,
  bool const rzero )
{
  int const n = GetN();
  double lder = lzero ? 0. : GetLeftSlope(n,mX.data(),y);
  double rder = rzero ? 0. : GetRightSlope(n,mX.data(),y);
  Reset(y,lder,rder);
}


void ccdl::CubicSpline::Reset
( int const n, 
  double const * x, 
  double const * y,
  bool const lzero,
  bool const rzero )
{
  double lder = lzero ? 0. : GetLeftSlope(n,x,y);
  double rder = rzero ? 0. : GetRightSlope(n,x,y);
  Reset(n,x,y,lder,rder);
}


