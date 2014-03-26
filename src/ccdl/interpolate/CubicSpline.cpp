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

    
ccdl::CubicSpline::CubicSpline
( int const n, 
  double const * x, 
  double const * y,
  bool const lzero,
  bool const rzero )
  : mX(x,x+n),
    mY(y,y+n)
{
  double lder = 0.;
  double rder = 0.;
  if ( not lzero ) lder = (y[1]-y[0])/(x[1]-x[0]);
  if ( not rzero ) rder = (y[n-2]-y[n-1])/(x[n-2]-x[n-1]);
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
  double lder = 0.;
  double rder = 0.;
  if ( not lzero ) lder = (y[1]-y[0]) / (mX[1]-mX[0]);
  if ( not rzero ) rder = (y[n-2]-y[n-1]) / (mX[n-2]-mX[n-1]);
  Reset(y,lder,rder);
}


void ccdl::CubicSpline::Reset
( int const n, 
  double const * x, 
  double const * y,
  bool const lzero,
  bool const rzero )
{
  double lder = 0.;
  double rder = 0.;
  if ( not lzero ) lder = (y[1]-y[0])/(x[1]-x[0]);
  if ( not rzero ) rder = (y[n-2]-y[n-1])/(x[n-2]-x[n-1]);
  Reset(n,x,y,lder,rder);
}


