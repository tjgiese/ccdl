#include "linemin.hpp"

ccdl::mini::OneDPoint 
ccdl::mini::linemin
( LINET type,
  int const n, 
  double const f0, 
  double const * g0,
  double const alp0,
  double const * s,
  ccdl::minifcn * pfcn )
{
  ccdl::minifcn & fcn = *pfcn;
  ccdl::mini::OneDPoint pt;
  pt.grd.resize( n );
  switch ( type )
    {
    case ccdl::mini::NOGRAD_NOPOLY:
      {
	pt.x = ccdl::mini::linemin_nograd_nopoly(n,f0,g0,alp0,s,pfcn);
	pt.y = fcn( pt.x, s, pt.grd.data() );
	pt = ccdl::mini::OneDPoint( pt.x, pt.y, n, s, pt.grd.data() );
	break;
      }
    // case ccdl::mini::NOGRAD_POLY:
    //   {
    // 	pt.x = ccdl::mini::linemin_nograd_poly(n,f0,x0,g0,alp0,s,fcn);
    // 	pt.y = fcn( x0, pt.x, s, pt.grd.data() );
    // 	pt = ccdl::mini::OneDPoint( pt.x, pt.y, n, s, pt.grd.data() );
    // 	break;
    //   }
    case ccdl::mini::GRAD_POLY:
      {
	ccdl::mini::OneDCurve curve = 
	  ccdl::mini::linemin_grad_poly(n,f0,g0,alp0,s,pfcn);
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
