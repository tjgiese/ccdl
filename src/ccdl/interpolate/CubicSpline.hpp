#ifndef _CubicSpline_hpp_
#define _CubicSpline_hpp_

#include "cspline.hpp"

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
    
    int const GetN() const { return mX.size(); }

    double GetLeftX() const { return mX[0]; }

    double GetRightX() const { return mX.back(); }

    double GetX( int const i ) const { return mX[i]; }

    double GetY( int const i ) const { return mY[2*i]; }

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

    void ConvertToPotential();

  private:
    std::vector<double> mX,mY;
  };

}

#endif

