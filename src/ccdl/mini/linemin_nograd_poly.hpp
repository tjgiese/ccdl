#ifndef _linemin_nograd_poly_hpp_
#define _linemin_nograd_poly_hpp_

#include "OneDCurve.hpp"
#include <cmath>

namespace ccdl
{
  namespace mini
  {


    template <class FCNVAL>
    ccdl::mini::OneDCurve
    FindBracket_NoGrad_Poly
    ( int const n, 
      double const f0, 
      double const * x0, 
      double const *, // g0
      double const alpha0,
      double const * s,
      FCNVAL fcn )
    {
      ccdl::mini::OneDCurve curve;

      double const big_gap = 0.5*(std::sqrt(5.)-1.);
      double const sml_gap = 1. - big_gap;

      // Initial lower bound
      double alo = 0.;
      double flo = f0;
      curve.push_back( alo, flo );


      // Initial upper bound
      double ahi = alpha0;
      double fhi = fcn( x0, ahi, s );
      curve.push_back( ahi, fhi );


      bool bracket = false;

      // Choose a golden bisection midpoint
      // based on which of the two bounds is lower in energy
      double acp = 0., fcp = 0.;
      if ( fhi < flo )
	acp = alo + (ahi-alo)*big_gap;
      else
	acp = alo + (ahi-alo)*sml_gap;
      fcp = fcn( x0, acp, s );
      curve.push_back( acp, fcp );


      if ( fcp < flo and fcp < fhi ) bracket = true;

      double scl = alpha0;
      double sclscl = 2.;
      while ( ! bracket )
	{
	  double apoly = 0.;
	  double ahes = 0.;
	  curve.SetBoltzmannWeight( 0.001 );
	  if ( curve.size() < 5 )
	    {
	      ccdl::polynomial poly( curve.fit( 3 ) );
	      if ( poly.minimum_loc( apoly ) )
		ahes  = poly( apoly, 2 );
	      else
		{
		  poly = curve.fit( 4, 2, 0. );
		  if ( poly.minimum_loc( apoly ) )
		    ahes  = poly( apoly, 2 );
		  else
		    ahes = -10000.;
		};
	    }
	  else
	    {
	      ccdl::polynomial poly( curve.fit( 5, 3, 0. ) );
	      if ( poly.minimum_loc( apoly ) )
		{
		  ahes  = poly( apoly, 2 );
		  if ( ahes < 0. )
		    {
		      poly = curve.fit( 5, 2, 0. );
		      if ( poly.minimum_loc( apoly ) )
			ahes = poly( apoly, 2 );
		      else
			ahes = -10000.;
		    }
		}
	      else
		ahes = -10000.;
	    };


	  if ( fcp > flo and curve.size() < 4 )
	    {
	      // try a golden section that is closer to alo
	  

	      ahi = acp;
	      fhi = fcp;
	      acp = alo + sml_gap * ( ahi-alo );

	  
	      if ( ahes > -1.e-8 and apoly > alo ) acp = apoly;

	      fcp = fcn( x0, acp, s );
	      curve.push_back( acp, fcp );
	    }
	  else if ( fcp > flo and curve.size() == 4 )
	    {
	      // Try a small change above alo
	      // This should yield a lower center point if
	      // the search direction is indeed pointing toward a minimum
	      ahi = acp;
	      fhi = fcp;
	      acp = alo + 0.001 * ( ahi-alo );
	      fcp = fcn( x0, acp, s );
	      curve.push_back( acp, fcp );
	    }
	  else if ( fcp > flo and curve.size() > 4 )
	    {
	      // If the case above (curve.size() == 4) failed, then
	      // we do not trust the search direction
	      // Let us instead try to search in the negative direction 
	      ahi = acp;
	      fhi = fcp;
	      acp = alo;
	      fcp = flo;
	      alo = acp - scl;
	      scl *= sclscl;
	      sclscl *= 1.2;

	      //if ( apoly < acp and apoly > (acp - 2.*scl) ) alo = apoly;
	      if ( ahes > -1.e-8 and apoly < acp ) alo = apoly - 0.05 * (ahi-acp);

	      flo = fcn( x0, alo, s );
	      curve.push_back( alo, flo );
	      // In the unlikely event that we happened to hop
	      // over the minimum and land at a point that exactly
	      // matches our current function value, then let
	      // us bisect between us
	      if ( fhi == flo )
		{
		  alo = (alo + acp)/2.;
		  flo = fcn( x0, alo, s );
		  curve.push_back( alo, flo );
		};
	    }
	  else if ( fcp > fhi )
	    {
	      // The input step length was too short to find a minimum
	      // We need to increase the step length
	      alo = acp;
	      flo = fcp;
	      acp = ahi;
	      fcp = fhi;
	      ahi = alo + scl;
	      scl *= sclscl;
	      sclscl *= 1.2;

	      if ( ahes > -1.e-8 and apoly > acp ) ahi = apoly + 0.05 * ( acp-alo );

	      fhi = fcn( x0, ahi, s );
	      curve.push_back( ahi, fhi );
	      // In the unlikely event that we happened to hop
	      // over the minimum and land at a point that exactly
	      // matches our current function value, then let
	      // us bisect between us
	      if ( fhi == flo )
		{
		  ahi = (ahi + acp)/2.;
		  fhi = fcn( x0, ahi, s );
		  curve.push_back( ahi, fhi );
		};
	    }
	  else
	    {
	      // If we were are in a region where y(x) = const.
	      // Then there's nothing we can do.
	      // Break out of the infinite loop
	      bracket = true;
	    }

	  if ( fcp < flo and fcp < fhi ) bracket = true;
	};

  
      // std::cout << "end of bracket routine\n";
      // curve.print( std::cout );
      // std::printf("the bracket is %20.10e %20.10e %20.10e\n               %20.10e %20.10e %20.10e\n",
      // 	      alo,acp,ahi,flo,fcp,fhi);
  
      return curve;
    }



    template <class FCNVAL>
    double linemin_nograd_poly
    ( int const n, 
      double const f0, 
      double const * x0, 
      double const *, // g0
      double const alpha0,
      double const * s,
      FCNVAL fcn,
      double const TOL = 1.e-8 )
    {
      double const big_gap = 0.5*(std::sqrt(5.)-1.);
      double const sml_gap = 1. - big_gap;
      
      ccdl::mini::OneDCurve curve( ccdl::mini::FindBracket_NoGrad_Poly( n, f0, x0, NULL, alpha0, s, fcn ) );
      // we now should be bracketed by (alo,acp,ahi)
      
      ccdl::mini::OneDCurve bcurve( curve.get_bracketing_curve() );
      
      //std::cout << "the bracketing curve\n";
      //bcurve.print( std::cout );
      
      double amin = curve[0].x;
      ccdl::polynomial poly( bcurve.fit(3) );
      if ( poly.minimum_loc( amin ) )
	{
	  double ahes = poly( amin, 2 );
	  if ( ahes < 0. )
	    amin = curve[0].x;
	}
      return amin;
    }


  }
}

#endif
