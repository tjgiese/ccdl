#ifndef _linemin_nograd_nopoly_hpp_
#define _linemin_nograd_nopoly_hpp_

#include "OneDCurve.hpp"
#include "minifcn.hpp"

namespace ccdl
{
  namespace mini
  {

    ccdl::mini::OneDCurve
    FindBracket_NoGrad_NoPoly
    ( int const n, 
      double const f0, 
      double const *, // g0
      double const alpha0,
      double const * s,
      ccdl::minifcn * pfcn );


    double linemin_nograd_nopoly_given_bracket
    ( int const n, 
      double const f0, 
      double const *, // g0
      ccdl::mini::OneDCurve const & BracketedCurve,
      double const * s,
      ccdl::minifcn * pfcn,
      double const TOL = 1.e-8 );


    double linemin_nograd_nopoly
    ( int const n, 
      double const f0, 
      double const *, // g0
      double const alpha0,
      double const * s,
      ccdl::minifcn * pfcn,
      double const TOL = 1.e-8 );


  }
}

#endif

