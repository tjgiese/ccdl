#ifndef _ccdl_minimizer_hpp_
#define _ccdl_minimizer_hpp_

#include "step.hpp"
#include "../bmath.hpp"



/*
struct fcn
{
  int n;

  fcn( int nx ) : n(nx) {}

  double operator() ( double const * x0, double const alpha, double const * s )
  {   
    std::vector<double> x( x0, x0 + n );
    for ( int i=0; i<n; ++i )
      x[i] += alpha * s[i];
   // eval f(x)
    return f;
  }

  double operator() ( double const * x0, double const alpha, double const * s, double * g )
  {   
    std::vector<double> x( x0, x0 + n );
    for ( int i=0; i<n; ++i )
      x[i] += alpha * s[i];
   // eval f(x) and g=df/dx
    return f;
  }
}

 */


namespace ccdl
{

  namespace mini
  {
    struct options
    {
      options()
	: maxiter(1000),
	  nsd(1),
	  ncg(3),
	  gmax_tol( 1.e-4 ),
	  grms_tol( 5.e-5 ),
	  ostr( &std::cout )
      {}

      int maxiter,nsd,ncg;
      double gmax_tol,grms_tol;
      std::ostream * ostr;
    };
  }

  template <class T>
  int minimize( int const n, double * x, T fcn, ccdl::mini::options opts = ccdl::mini::options() );
}



template <class FCNVAL>
int ccdl::minimize( int const n, double * x, FCNVAL fcn, ccdl::mini::options opts )
{
  int const inc = 1;
  
#define FMTE(a) std::setw(11) << std::setprecision(2) << std::scientific << (a) << " "
#define FMTF(a) std::setw(22) << std::setprecision(12) << std::fixed << (a) << " "
#define PRT cout << "GEOMOPT " << std::setw(5) << iter++ << FMTF(step.f) << FMTE(step.gmax) << FMTE(step.grms) << "\n";

  std::ostream & cout = *(opts.ostr);
  int iter=0;

  ccdl::mini::BfgsStep step(n), prev(n);
  std::vector<double> snew(n,0.);
  step.set_s( x );
  step.cpt_fg( 1., fcn );
  PRT;


  prev = step;
  if ( step.gmax < 1.e-15 )
    {
      prev.s[0] += 5.e-5;
      prev.cpt_fg( 1., fcn );
      PRT;


      if ( prev.f >= step.f )
	{ // the user gave us a minimum as the starting point -- return
	  return 0;
	};
      step = prev;
    };
  
  double alp0 = 0.;
  double alp = 0.;

  
  for ( int iloop = 0; iter < opts.maxiter; ++iloop )
    {
      
      if ( iloop < opts.nsd )
	{
	  //for ( int i=0; i<n; ++i ) step.s[i] = -step.g[i];

	  step.steepdesc_direction();

	  alp0 = 1. / ( 1. + std::sqrt( ccdl::v_dot_v( n, step.s.data(), step.s.data() ) ) );
	}
      else if ( iloop < opts.nsd + opts.ncg )
	{
	  step.hestenes_stiefel_direction(&prev,alp0);
	  //step.polak_ribiere_direction(&prev,alp0);
	  //step.fletcher_reeves_direction(&prev,alp0);
	}
      else
	{
	  step.bfgs_direction();
	}

      double sdir = ccdl::v_dot_v( n, step.s.data(), step.g.data() );
      cout << FMTE(sdir) << "\n";

      prev = step;

      alp = step.linemin( alp0, fcn );
      PRT;

      //std::printf("GREP %5i %20.10f %11.2e %11.2e (alp0=%15.5f alp=%15.5f)\n\n\n\n\n\n\n",iter,step.f,step.gmax,step.grms,alp0,alp);

      if ( step.gmax < opts.gmax_tol and step.grms < opts.grms_tol ) break;

      step.update_hessian( prev );

      alp0 = alp;
      if ( alp0 < 0. ) alp0 = 1. / ( 1. + ccdl::v_dot_v(n,step.g.data(),step.g.data()) );
      alp0 = std::min( alp0, 1. );
      //alp0 = 1.;

    };

  std::copy( step.x.data(), step.x.data()+n, x );
  
}

#endif
