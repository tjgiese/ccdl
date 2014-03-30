#ifndef _MSpline_hpp_
#define _MSpline_hpp_

#include "cspline.hpp"

namespace ccdl
{

  class MSpline
  {
    MSpline
    ( int const N, 
      double const * x );

    void push_back
    ( double const * y, 
      double const lder, 
      double const rder );
    
    void push_back
    ( double const * y, 
      bool const lzero,
      bool const rzero );

    void Erase();

    int GetN() const { return mX.size(); }

    int GetNspl() const { return mY.size(); }

    double GetLeftX() const { return mX[0]; }

    double GetRightX() const { return mX.back(); }

    double GetX( int const i ) const { return mX[i]; }

    void GetX( double * x ) const;

    double GetY
    ( int const ispl, 
      int const ipt ) const 
    { 
      return mY[ispl][2*ipt]; 
    }
    
    void GetY
    ( int const ispl, 
      double * y ) const;

    void GetValues
    ( double const x, 
      double * v ) const;

    void GetValuesAndDerivs
    ( double const x, 
      double * y, 
      double * dydx ) const;

    void GetValuesAndSecondDerivs
    ( double const x, 
      double * y, 
      double * dydx, 
      double * d2ydx2 ) const;

  private:
    std::vector<double> mX;
    std::vector< std::vector<double> > mY;
  };

}

#endif

