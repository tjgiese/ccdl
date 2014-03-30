#ifndef _HarmonicSpline_hpp_
#define _HarmonicSpline_hpp_

#include "cspline.hpp"

namespace ccdl
{

  class HarmonicSpline
  {

    HarmonicSpline
    ( int const Lmax,
      int const N, 
      double const * x );

    void Reset
    ( int const Lmax, 
      int const N, 
      double const * x );

    std::vector<double> r;
    std::vector< std::vector<double> > ylm;

    void Commit();

    void ConvertToPotentials();

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

  };

}


#endif
