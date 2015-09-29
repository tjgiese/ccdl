#ifndef _linemin_hpp_
#define _linemin_hpp_

#include "linemin_grad_poly.hpp"
#include "linemin_nograd_nopoly.hpp"
#include "linemin_nograd_poly.hpp"

#include <cstdlib>
#include <iostream>

namespace ccdl
{
  namespace mini
  {
    enum LINET { NOGRAD_NOPOLY, NOGRAD_POLY, GRAD_POLY };

    ccdl::mini::OneDPoint linemin
    ( LINET type,
      int const n, 
      double const f0, 
      double const * g0,
      double const alpha0,
      double const * s,
      ccdl::minifcn * pfcn );

  }
}



#endif
