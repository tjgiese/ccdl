#include "linemin_nograd_nopoly.hpp"
#include <cmath>

ccdl::mini::OneDCurve
ccdl::mini::FindBracket_NoGrad_NoPoly
( int const , //n, 
  double const f0, 
  double const *, // g0
  double const alpha0,
  double const * s,
  ccdl::minifcn * pfcn )
{
  ccdl::minifcn & fcn = *pfcn;

  ccdl::mini::OneDCurve curve;
      
  double const big_gap = 0.5*(std::sqrt(5.)-1.);
  double const sml_gap = 1. - big_gap;
      
  // Initial lower bound
  double alo = 0.;
  double flo = f0;
  curve.push_back( alo, flo );
      
      
  // Initial upper bound
  double ahi = alpha0;
  double fhi = fcn( ahi, s );
  curve.push_back( ahi, fhi );
      
      
  bool bracket = false;
      
  // Choose a golden bisection midpoint
  // based on which of the two bounds is lower in energy
  double acp = 0., fcp = 0.;
  if ( fhi < flo )
    acp = alo + (ahi-alo)*big_gap;
  else
    acp = alo + (ahi-alo)*sml_gap;
  fcp = fcn( acp, s );
  curve.push_back( acp, fcp );
      
      
  if ( fcp < flo and fcp < fhi ) bracket = true;
      
  double scl = alpha0;
  double sclscl = 2.;
  while ( ! bracket )
    {
      if ( fcp > flo and curve.size() < 4 )
	{
	  // try a golden section that is closer to alo
	  ahi = acp;
	  fhi = fcp;
	  acp = alo + sml_gap * ( ahi-alo );
	  fcp = fcn( acp, s );
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
	  fcp = fcn( acp, s );
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
	  flo = fcn( alo, s );
	  curve.push_back( alo, flo );
	  // In the unlikely event that we happened to hop
	  // over the minimum and land at a point that exactly
	  // matches our current function value, then let
	  // us bisect between us
	  if ( fhi == flo )
	    {
	      alo = (alo + acp)/2.;
	      flo = fcn( alo, s );
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
	  fhi = fcn( ahi, s );
	  curve.push_back( ahi, fhi );
	  // In the unlikely event that we happened to hop
	  // over the minimum and land at a point that exactly
	  // matches our current function value, then let
	  // us bisect between us
	  if ( fhi == flo )
	    {
	      ahi = (ahi + acp)/2.;
	      fhi = fcn( ahi, s );
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



double ccdl::mini::linemin_nograd_nopoly_given_bracket
( int const , //n, 
  double const , //f0, 
  double const *, // g0
  ccdl::mini::OneDCurve const & BracketedCurve,
  double const * s,
  ccdl::minifcn * pfcn,
  double const TOL  )
{
  ccdl::minifcn & fcn = *pfcn;

  double const big_gap = 0.5*(std::sqrt(5.)-1.);
  double const sml_gap = 1. - big_gap;

  OneDCurve bcurve( BracketedCurve.get_y_bracketing_curve() );
  bcurve.sort_by_y();

  //std::cout << "the bracketing curve\n";
  //bcurve.print( std::cout );


  double flo = bcurve[0].y;
  double fcp = bcurve[1].y;
  double fhi = bcurve[2].y;

  double alo = bcurve[0].x;
  double acp = bcurve[1].x;
  double ahi = bcurve[2].x;

  // If we were given the line y(x) = const.
  // Then the bisection failed and there's nothing more to do
  if ( flo == fhi ) return acp;

  // std::printf("bracket %20.10f %20.10f %20.10f  (%20.10f %20.10f %20.10f)\n",
  // 	      alo,acp,ahi,flo,fcp,fhi);



  // we can now to the gold bisection
  double amin = alo;
  while ( true )
    {
      double atrial = 0;
      double ftrial = 0;
      if ( (acp-alo) >= (ahi-acp) )
	{
	  atrial = alo + big_gap * ( acp-alo );
	  ftrial = fcn( atrial, s );
	  bcurve.push_back( atrial, ftrial );
	  if ( ftrial < fcp )
	    {
	      fhi = fcp;
	      ahi = acp;
	      fcp = ftrial;
	      acp = atrial;
	    }
	  else if ( ftrial >= fcp )
	    {
	      flo = ftrial;
	      alo = atrial;
	    }
	}
      else
	{
	  atrial = acp + sml_gap * ( ahi-acp );
	  ftrial = fcn( atrial, s );
	  bcurve.push_back( atrial, ftrial );
	  if ( ftrial < fcp )
	    {
	      flo = fcp;
	      alo = acp;
	      fcp = ftrial;
	      acp = atrial;
	    }
	  else if ( ftrial >= fcp )
	    {
	      fhi = ftrial;
	      ahi = atrial;
	    }
	}

      if ( ahi-alo <= std::sqrt( TOL ) * ( std::abs(acp) + std::abs(atrial) ) )
	break;

      amin = bcurve[0].x;
    };

  amin = bcurve[0].x;
  return amin;
}



double ccdl::mini::linemin_nograd_nopoly
( int const n, 
  double const f0, 
  double const *, // g0
  double const alpha0,
  double const * s,
  ccdl::minifcn * pfcn,
  double const TOL )
{
  //ccdl::minifcn & fcn = *pfcn;

  //double const big_gap = 0.5*(std::sqrt(5.)-1.);
  //double const sml_gap = 1. - big_gap;

  OneDCurve curve( ccdl::mini::FindBracket_NoGrad_NoPoly( n, f0, NULL, alpha0, s, pfcn ) );
  // we now should be bracketed by (alo,acp,ahi)

  OneDCurve bcurve( curve.get_bracketing_curve() );
  return ccdl::mini::linemin_nograd_nopoly_given_bracket( n, f0, NULL, bcurve, s, pfcn, TOL );

  /*

  std::cout << "the bracketing curve\n";
  bcurve.print( std::cout );


  double flo = bcurve[0].y;
  double fcp = bcurve[1].y;
  double fhi = bcurve[2].y;

  double alo = bcurve[0].x;
  double acp = bcurve[1].x;
  double ahi = bcurve[2].x;

  // If we were given the line y(x) = const.
  // Then the bisection failed and there's nothing more to do
  if ( flo == fhi ) return acp;

  // std::printf("bracket %20.10f %20.10f %20.10f  (%20.10f %20.10f %20.10f)\n",
  // 	      alo,acp,ahi,flo,fcp,fhi);



  // we can now to the gold bisection
  double amin = curve[0].x;
  while ( true )
    {
      double atrial = 0;
      double ftrial = 0;
      if ( (acp-alo) >= (ahi-acp) )
	{
	  atrial = alo + big_gap * ( acp-alo );
	  ftrial = fcn( atrial, s );
	  curve.push_back( atrial, ftrial );
	  if ( ftrial < fcp )
	    {
	      fhi = fcp;
	      ahi = acp;
	      fcp = ftrial;
	      acp = atrial;
	    }
	  else if ( ftrial >= fcp )
	    {
	      flo = ftrial;
	      alo = atrial;
	    }
	}
      else
	{
	  atrial = acp + sml_gap * ( ahi-acp );
	  ftrial = fcn( atrial, s );
	  curve.push_back( atrial, ftrial );
	  if ( ftrial < fcp )
	    {
	      flo = fcp;
	      alo = acp;
	      fcp = ftrial;
	      acp = atrial;
	    }
	  else if ( ftrial >= fcp )
	    {
	      fhi = ftrial;
	      ahi = atrial;
	    }
	}

      if ( ahi-alo <= std::sqrt( TOL ) * ( std::abs(acp) + std::abs(atrial) ) )
	break;

      amin = curve[0].x;
    };

  amin = curve[0].x;
  return amin;
  */
}
