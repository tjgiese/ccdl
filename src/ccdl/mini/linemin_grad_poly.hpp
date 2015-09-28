#ifndef _linemin_grad_poly_hpp_
#define _linemin_grad_poly_hpp_

#include "OneDCurve.hpp"
#include <cmath>
#include <cstdio>
#include "../bmath.hpp"

namespace ccdl
{
  namespace mini
  {


    template <class FCNVAL>
    ccdl::mini::OneDCurve
    linemin_grad_poly
    ( int const n, 
      double const f0, 
      double const * x0, 
      double const * g0,
      double const alpha0,
      double const * s,
      FCNVAL fcn )
    {
      ccdl::mini::OneDCurve curve;

      std::vector<double> g( n, 0. );

      double const big_gap = 0.5*(std::sqrt(5.)-1.);
      double const sml_gap = 1. - big_gap;

      // Initial lower bound
      double alo = 0.;
      double flo = f0;
      double glo = ccdl::v_dot_v( n, s, g0 );
      curve.push_back( alo, flo, n, s, g0 );


      // Initial upper bound
      double ahi = alpha0;
      //if ( glo > 0. ) ahi *= -1.;
      double fhi = fcn( x0, ahi, s, g.data() );
      curve.push_back( ahi, fhi, n, s, g.data() );

      alo = curve[0].x;
      ahi = curve[1].x;
      flo = curve[0].y;
      fhi = curve[1].y;


      bool bracket = false;

      if ( curve[0].g * curve[1].g < 0. )
	bracket = true;
 
      //for ( int i=0; i<curve.size(); ++i )
      //std::printf("%4i %20.10e %20.10e %20.10e\n",i,curve[i].x,curve[i].y,curve[i].g);

      double aold = curve[ curve.size() - 1 ].x;
      double scl = 2.;
      while ( ! bracket )
	{
	  curve.SetBoltzmannWeight( 1.e-5 );
	  ccdl::polynomial poly;
	  if ( curve.size() == 2 )
	    poly = curve.fit( 5, 3, 0. );
	  else
	    poly = curve.fit( std::min(5,curve.size()+2) );
	  double acp = 0.;
	  curve.sort_by_x();
	  alo = curve[ curve.size()-2 ].x;
	  ahi = curve[ curve.size()-1 ].x;
	  curve.sort_by_y();

	  if ( ! poly.minimum_loc( acp ) )
	    {
	      poly = curve.fit( 4 );
	      if ( ! poly.minimum_loc( acp ) )
		{
		  poly = curve.fit( 3 );
		  if ( ! poly.minimum_loc( acp ) )
		    {
		      acp = (alo-ahi)*scl;
		      scl *= 2.;
		    }
		  else
		    {
		      acp += std::min( (acp-aold)*0.1, 2. );
		    }
		}
	      else
		{
		  acp += std::min( (acp-aold)*0.1, 1.5 );
		}
	    }
	  else
	    {
	      acp += std::min( (acp-aold)*0.1, 1. );
	    }
	  aold = acp;
	  double fcp = fcn( x0, acp, s, g.data() );
	  curve.push_back( acp, fcp, n, s, g.data() );

	  // curve.sort_by_x();
	  // for ( int i=0; i<curve.size(); ++i )
	  //   std::printf("%4i %20.10e %20.10e %20.10e\n",i,curve[i].x,curve[i].y,curve[i].g);

	  bracket = curve.g_is_bracketed();
	};


      curve.SetBoltzmannWeight( 1.e-6 );

      curve.sort_by_y();
      fhi = curve[0].y;
      flo = curve.back().y; 
      while ( flo > fhi )
	{
	  double acp = 0.;
	  ccdl::polynomial poly( curve.fit( 5 ) );
	  if ( ! poly.minimum_loc( acp ) )
	    {
	      poly = curve.fit( 4 );
	      if ( ! poly.minimum_loc( acp ) )
		poly = curve.fit( 3 );
	      else
		break;
	    };
	  flo = fcn( x0, acp, s, g.data() );
	  curve.push_back( acp, flo, n, s, g.data() );
	};

      return curve;
    }


  }
}

#endif

