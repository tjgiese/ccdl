#include "linemin_grad_poly.hpp"
#include "linemin_nograd_nopoly.hpp"
#include "../bmath.hpp"

#include <cmath>
#include <cstdio>

#undef DPRT

ccdl::mini::OneDCurve
ccdl::mini::linemin_grad_poly
( int const n, 
  double const f0, 
  double const * g0,
  double const alpha0,
  double const * s,
  ccdl::minifcn * pfcn )
{
  ccdl::minifcn & fcn = *pfcn;

  ccdl::mini::OneDCurve curve;

  std::vector<double> g( n, 0. );

  //double const big_gap = 0.5*(std::sqrt(5.)-1.);
  //double const sml_gap = 1. - big_gap;

  // Initial lower bound
  double alo = 0.;
  double flo = f0;
  //double glo = ccdl::v_dot_v( n, s, g0 );
  curve.push_back( alo, flo, n, s, g0 );


  //std::printf("s.g = %20.10e\n", glo );

  double s2 = ccdl::v_dot_v(n,s,s);
  //std::printf("s2 = %20.10e alp*s2 %20.10e\n",s2,alpha0*s2);

  // Initial upper bound
  double ahi = alpha0;
  if ( alpha0 * s2 > 2.*n )
    {
      ahi = 2.*n / s2;
      //std::printf("s2 = %20.10e alp*s2 %20.10e\n",s2,ahi*s2);
    };
  //if ( glo > 0. ) ahi *= -1.;
  double fhi = fcn( ahi, s, g.data() );
  curve.push_back( ahi, fhi, n, s, g.data() );

  alo = curve[0].x;
  ahi = curve[1].x;
  flo = curve[0].y;
  fhi = curve[1].y;


  bool bracket = false;

  if ( curve[0].g * curve[1].g < 0. )
    bracket = true;
 
  if ( bracket )
    {
#ifdef DPRT
      std::printf("bracketed\n");
#endif

      ccdl::mini::OneDPoint plo = curve[0];
      while ( curve[1].g > -10 * curve[0].g and std::abs(curve[0].g) > 0.001 )
	{
#ifdef DPRT
	  std::printf("trying to bisect big difference in gradients\n");
#endif
	  //double w = std::abs(curve[1].g) / ( std::abs(curve[0].g) + std::abs(curve[1].g) );
	  double w = 0.5;
	  double ahi0 = ahi;

	  ahi = alo*w + ahi0*(1.-w);
	      

	  fhi = fcn( ahi, s, g.data() );
	      
	  double sdot = ccdl::v_dot_v( n, s, g.data() );
	  
	  if ( curve[0].g * sdot > 0. ) // we went too far
	    {
#ifdef DPRT
	      std::printf("we went too far\n");
#endif
	      ahi = ahi0;
	      break;
	    }
	  
	  curve.resize(0);
	  curve.data.push_back( plo );
	  curve.push_back( ahi, fhi, n, s, g.data() );

	  if ( fhi < flo )
	    {
#ifdef DPRT
	      std::printf("we lowered the value on our first attempt\n");
#endif
	      return curve;
	    };
	};
    }

  //for ( int i=0; i<curve.size(); ++i )
  //std::printf("%4i %20.10e %20.10e %20.10e\n",i,curve[i].x,curve[i].y,curve[i].g);

  //double aold = curve[ curve.size() - 1 ].x;
  double fold = fhi;
  double scl = 2.;
  int iter = 0;
  while ( ! bracket )
    {
      ++iter;
#ifdef DPRT
      std::printf("search for bracket attempt %i\n",iter);
#endif

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
	  
      double minf  = curve[0].y;
      double predf = curve.back().y;

      if ( ! poly.minimum_loc( acp ) )
	{
	  poly = curve.fit( 4 );
	  if ( ! poly.minimum_loc( acp ) )
	    {
	      poly = curve.fit( 3 );
	      if ( ! poly.minimum_loc( acp ) )
		{
		  acp = ahi + (ahi-alo)*scl;
		  scl *= 1.5;
		}
	      else
		{
		  acp *= 1.1;
		  //acp = 1.5*ahi;
		  //curve.print( std::cerr );
		  //acp += std::min( std::min( (acp-aold)*0.1, 2. ) * scl, acp );
		  //scl *= 1.5;
		  predf = poly( acp );
		}
	    }
	  else
	    {
	      acp *= 1.1;
	      //acp = 1.5*ahi;
	      //curve.print( std::cerr );
	      //acp += std::min( std::min( (acp-aold)*0.1, 1.5 ) * scl, acp );
	      //scl *= 1.25;
	      predf = poly( acp );
	    }
	}
      else
	{
	  acp *= 1.1;
	  //acp = 1.5*ahi;
	  //curve.print( std::cerr );
	  //acp += std::min( std::min( (acp-aold)*0.1, 1. ) * scl, acp );
	  //scl *= 1.125;
	  predf = poly( acp );
	}

      if ( acp < alo or iter > 4 )
	{
	  if ( acp < 1.5*ahi )
	    acp = 1.5*ahi;
	}
      else if ( acp < ahi ) 
	{
	  acp = 1.5*ahi;
	  //scl *= 1.5;
	};	  

      //std::printf("predf %24.12f  predf-f0 < 1.25 * (minf-f0) : %24.12f < %24.12f\n",
      //predf,predf-f0,1.25 * (minf-f0));

      // IF INEXACT LINE SEARCH
      if ( ( f0-predf < 1.25 * (f0-minf) ) and fold == fhi ) return curve;
	    
      //std::printf("acp = %20.10e %20.10e %20.10e\n",acp,alo,ahi);
      //aold = acp;
      double fcp = fcn( acp, s, g.data() );
      fold = fcp;
      curve.push_back( acp, fcp, n, s, g.data() );

      //std::printf("eee %20.10e\n",(fcp-f0) / std::abs(f0) );
      if ( (fcp-f0) / std::abs(f0) < -2./3. ) return curve; 
      // curve.sort_by_x();
      // for ( int i=0; i<curve.size(); ++i )
      //   std::printf("%4i %20.10e %20.10e %20.10e\n",i,curve[i].x,curve[i].y,curve[i].g);

      bracket = curve.g_is_bracketed();
    };


  curve.SetBoltzmannWeight( 1.e-6 );

  curve.sort_by_y();
  fhi = curve[0].y;
  flo = curve.back().y; 

#ifdef DPRT
  std::printf("bracket found between %20.10f %20.10f\n",curve[0].x,curve.back().x);
#endif

  double aold = curve.back().x;
  iter = 0;
  while ( flo > fhi )
    {
      iter++;

#ifdef DPRT
      std::printf("poly bisection attempt %i\n",iter);
#endif


      ccdl::mini::OneDCurve bc( curve.get_bracketing_curve() );
      bc.sort_by_y();

      double minf  = bc[0].y;
      double predf = bc.back().y;

      bc.sort_by_x();
      double acp = 0.;
      ccdl::polynomial poly;
      poly = curve.fit( 5 );

      //if ( curve.size() == 2 )
      //poly = curve.fit( 5, 3, 0. );
      //else
      //poly = curve.fit( std::min(5,curve.size()+2) );
      if ( ! poly.minimum_loc( acp ) )
	{
	  poly = curve.fit( 4 );
	  if ( ! poly.minimum_loc( acp ) )
	    {
	      poly = curve.fit( 3 );
	      poly.minimum_loc( acp );
	      /*
		std::vector<double> ms( poly.minima_loc() );
		for ( std::size_t i=0; i<ms.size(); ++i )
		if ( ms[i] > bc[0].x and ms[i] < bc[1].x )
		{
		acp = ms[i];
		break;
		}
	      */
	    }
	  else
	    break;
	}
      else
	{
	  /*
	    std::vector<double> ms( poly.minima_loc() );
	    for ( std::size_t i=0; i<ms.size(); ++i )
	    if ( ms[i] > bc[0].x and ms[i] < bc[1].x )
	    {
	    acp = ms[i];
	    break;
	    }
	  */
	};

      //std::cout << "checking " << flo << " " << acp << " " << bc[0].x << " " << bc[1].x << " " << curve.size() << "\n";
	  


      if ( acp > bc[1].x or acp < bc[0].x or 
				  iter > 1 or 
	(bc[1].g > -100.*bc[0].g) or
	   std::abs(bc[0].g) > 500. )
	{
	  //std::printf("out of bounds\n");
	  double w = std::abs( bc[1].g ) / ( std::abs(bc[0].g) + std::abs(bc[1].g) );
	  acp = bc[0].x * w + bc[1].x * ( 1.-w );
	  /*
	    poly = bc.fit( 3, false );
	    bool ok = poly.minimum_loc( acp );
	    if ( (! ok) or ( acp > bc[1].x or acp < bc[0].x ) )
	    acp = 0.5 * ( bc[1].x + bc[0].x );
	  */
	};


      predf = poly( acp );
      //std::printf("predf %24.12f  predf-f0 < 1.25 * (minf-f0) : %24.12f < %24.12f\n",
      //predf,predf-f0,1.25 * (minf-f0));

      // IF INEXACT LINE SEARCH
      if ( ( f0-predf < 1.25 * (f0-minf) ) and fold == fhi ) return curve;

      if ( acp == aold )
	{

#ifdef DPRT
	  std::printf("alpha is the same as the previous, do a golden search\n");
#endif
	  acp = ccdl::mini::linemin_nograd_nopoly
	    ( n, f0, g0, bc[1].x, s, pfcn, 1.e-4 );
	  flo = fcn( acp, s, g.data() );
	  curve.push_back( acp, flo, n, s, g.data() );
	  break;
	};

#ifdef DPRT
      std::printf("alpha attempt %20.10f  (%20.10f .. %20.10f)\n",acp,bc[0].x,bc[1].x);
#endif
      flo = fcn( acp, s, g.data() );
      fold = flo;
      aold = acp;

      double sval = ccdl::v_dot_v( n, s, g.data() );


      // std::printf("checking acp %22.16f  bc[0] %22.16f  bc[1] %22.16f\n",
      // 	      acp, bc[0].x, bc[1].x );
      // std::printf("checking flo %22.16f  bc[0] %22.16f  bc[1] %22.16f\n",
      // 	      flo, bc[0].y, bc[1].y );
      // std::printf("checking glo %22.16f  bc[0] %22.16f  bc[1] %22.16f\n",
      // 	      sval, bc[0].g, bc[1].g );


      if ( sval > bc[1].g or flo > bc[1].y )
	{
	  curve.sort_by_y();
	  //std::printf("sval > bc[1].g %20.10e\n",curve[0].x);
	  if ( curve[0].x > 1.e-13 and sval < 0. ) break;
	  curve.push_back( acp, flo, n, s, g.data() );
	  if ( curve.y_is_bracketed() )
	    curve = curve.get_y_bracketing_curve();
	}
      else
	curve.push_back( acp, flo, n, s, g.data() );


      if ( iter > 4 )
	{
	  if ( bc[0].x > 1.e-4 ) break;
	  else
	    {
#ifdef DPRT
	      std::printf("do a golden search\n");
#endif
	      acp = ccdl::mini::linemin_nograd_nopoly_given_bracket
		( n, f0, g0, curve, s, pfcn, 1.e-4 );
	      flo = fcn( acp, s, g.data() );
	      curve.push_back( acp, flo, n, s, g.data() );
	      break;
	    }
	}
      
      if ( iter > 20 ) 
	{
#ifdef DPRT
	  std::printf("exit failure\n");
#endif
	  break;
	};

      //if ( flo > fhi )
      //curve = curve.get_bracketing_curve();
    }
#ifdef DPRT
  std::printf("end of line search\n");
#endif
  curve.sort_by_y();
  return curve;
}
