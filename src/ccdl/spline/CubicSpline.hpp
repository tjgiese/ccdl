#ifndef _CubicSpline_hpp_
#define _CubicSpline_hpp_

#include "cspline.hpp"
#include "../constants.hpp"

namespace ccdl
{

  class CubicSpline
  {
  public:

    CubicSpline
    ( int const n, 
      double const * x, 
      double const * y, 
      double const lder, 
      double const rder );
    
    CubicSpline
    ( int const n, 
      double const * x, 
      double const * y,
      bool const lzero,
      bool const rzero );
    
    int GetN() const;

    double GetLeftX() const;

    double GetRightX() const;

    double GetX( int const i ) const;

    double GetY( int const i ) const;

    void GetX( double * x ) const;

    void GetY( double * y ) const;

    void Reset
    ( double const * y,
      double const lder, 
      double const rder );
    
    void Reset
    ( int const n, 
      double const * x, 
      double const * y, 
      double const lder, 
      double const rder );
    
    void Reset
    ( double const * y,
      bool const lzero,
      bool const rzero );
    
    void Reset
    ( int const n, 
      double const * x, 
      double const * y,
      bool const lzero,
      bool const rzero );

    double GetValue( double const x ) const;

    double GetValueAndDeriv
    ( double const x, 
      double & dydx ) const;
    
    double GetValueAndSecondDeriv
    ( double const x, 
      double & dydx, 
      double & d2ydx2 ) const;
    
    template <int X_POWER>
    double Integrate
    ( double const xlo,
      double const xhi ) const;

    void YlmDensityToYlmPotential();

    void YlmDensityToYlmPotential( double const charge );

    void RadialDensityToRadialPotential();

    void RadialDensityToRadialPotential( double charge );


    void RadialDensityToSphericalPotential();

    void RadialDensityToSphericalPotential( double charge );

    void operator*= ( double const d );

    void operator/= ( double const d );


  private:

    std::vector<double> mX,mY;

  };

}



inline int ccdl::CubicSpline::GetN() const 
{ 
  return mX.size(); 
}


inline double ccdl::CubicSpline::GetLeftX() const 
{ 
  return mX[0]; 
}


inline double ccdl::CubicSpline::GetRightX() const 
{ 
  return mX.back(); 
}


inline double ccdl::CubicSpline::GetX
( int const i ) const 
{ 
  return mX[i]; 
}


inline double ccdl::CubicSpline::GetY
( int const i ) const 
{ 
  return mY[2*i]; 
}


inline void ccdl::CubicSpline::GetX
( double * x ) const
{
  int const n = GetN();
  for ( int i=0; i<n; ++i )
    x[i] = mX[i];
}


inline void ccdl::CubicSpline::GetY
( double * y ) const
{
  int const n = GetN();
  for ( int i=0; i<n; ++i )
    y[i] = mY[i];
}


inline double ccdl::CubicSpline::GetValue
( double const x ) const
{
  return ccdl::cspline_value
    (x,
     GetN(),mX.data(),mY.data());
}


inline double ccdl::CubicSpline::GetValueAndDeriv
( double const x,
  double & dydx ) const
{
  return ccdl::cspline_value_and_deriv
    (x,dydx,
     GetN(),mX.data(),mY.data());
}


inline double ccdl::CubicSpline::GetValueAndSecondDeriv
( double const x,
  double & dydx,
  double & d2ydx2 ) const
{
  return ccdl::cspline_value_and_secondderiv
    (x,dydx,d2ydx2,
     GetN(),mX.data(),mY.data());
}


template <int NINJA_POWER>
inline double ccdl::CubicSpline::Integrate
( double const xlo,
  double const xhi ) const
{
  return ccdl::cspline_integral<NINJA_POWER>
    (xlo,xhi,
     GetN(),mX.data(),mY.data());
}


inline void ccdl::CubicSpline::YlmDensityToYlmPotential()
{
  std::vector<double> data(GetN(),0.);
  ccdl::cspline_lmpotential<0>
    (GetN(),mX.data(),&mY,&data,NULL);
  Reset(data.data(),false,false);
}


inline void ccdl::CubicSpline::YlmDensityToYlmPotential( double const norm )
{
  std::vector<double> data(GetN(),0.);
  ccdl::cspline_lmpotential<0>
    (GetN(),mX.data(),&mY,&data,&norm);
  Reset(data.data(),false,false);
}


inline void ccdl::CubicSpline::RadialDensityToRadialPotential()
{
  /*
  double const FOUR_PI = 4. * 3.141592653589793238462643383279502884197;
  double const SQRT_4PI = std::sqrt( FOUR_PI );
  // angular integral of density times Y00 harmonic
  (*this) *= FOUR_PI / SQRT_4PI;
  YlmDensityToYlmPotential();
  // multiply potential by Y00 harmonic
  (*this) *= 1. / SQRT_4PI;
  */
  YlmDensityToYlmPotential();
}


inline void ccdl::CubicSpline::RadialDensityToRadialPotential( double norm )
{
  //double const FOUR_PI = 4. * 3.141592653589793238462643383279502884197;
  double const SQRT_4PI = std::sqrt( ccdl::FOUR_PI );
  YlmDensityToYlmPotential( norm / SQRT_4PI );
}


inline void ccdl::CubicSpline::operator*= ( double const d )
{
  int const n = mY.size();
  for ( int i=0; i<n; ++i )
    mY[i] *= d;
}


inline void ccdl::CubicSpline::operator/= ( double const d )
{
  int const n = mY.size();
  for ( int i=0; i<n; ++i )
    mY[i] /= d;
}


#endif

