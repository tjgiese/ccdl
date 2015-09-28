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

    template<class FCNVAL>
    ccdl::mini::OneDPoint linemin
    ( LINET type,
      int const n, 
      double const f0, 
      double const * x0, 
      double const * g0,
      double const alpha0,
      double const * s,
      FCNVAL fcn );

  }
}

template<class FCNVAL>
ccdl::mini::OneDPoint 
ccdl::mini::linemin
( LINET type,
  int const n, 
  double const f0, 
  double const * x0, 
  double const * g0,
  double const alp0,
  double const * s,
  FCNVAL fcn )
{
  ccdl::mini::OneDPoint pt;
  pt.grd.resize( n );
  switch ( type )
    {
    case ccdl::mini::NOGRAD_NOPOLY:
      {
	pt.x = ccdl::mini::linemin_nograd_nopoly(n,f0,x0,g0,alp0,s,fcn);
	pt.y = fcn( x0, pt.x, s, pt.grd.data() );
	pt = ccdl::mini::OneDPoint( pt.x, pt.y, n, s, pt.grd.data() );
	break;
      }
    case ccdl::mini::NOGRAD_POLY:
      {
	pt.x = ccdl::mini::linemin_nograd_poly(n,f0,x0,g0,alp0,s,fcn);
	pt.y = fcn( x0, pt.x, s, pt.grd.data() );
	pt = ccdl::mini::OneDPoint( pt.x, pt.y, n, s, pt.grd.data() );
	break;
      }
    case ccdl::mini::GRAD_POLY:
      {
	ccdl::mini::OneDCurve curve = 
	  ccdl::mini::linemin_grad_poly(n,f0,x0,g0,alp0,s,fcn);
	pt = curve[0];
	break;
      }
    default:
      {
	std::cerr << "Invalid line minimizer\n";
	std::abort();
      }
    };
  return pt;
}


#endif
