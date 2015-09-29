#ifndef _linemin_grad_poly_hpp_
#define _linemin_grad_poly_hpp_

#include "OneDCurve.hpp"
#include "minifcn.hpp"

namespace ccdl
{
  namespace mini
  {

    ccdl::mini::OneDCurve
    linemin_grad_poly
    ( int const n, 
      double const f0, 
      double const * g0,
      double const alpha0,
      double const * s,
      ccdl::minifcn * pfcn );
  }
}

#endif

