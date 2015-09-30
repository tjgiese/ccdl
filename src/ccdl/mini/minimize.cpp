#include "minimize.hpp"
#include "step.hpp"
#include "../bmath.hpp"
#include <iomanip>
#include <cmath>


//
// EXAMPLES
//
/*

#include <ccdl/mini.hpp>
#include <cmath>
#include <cstdio>


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

struct Rosenbrock : public ccdl::minifcn
{
  // minimum @ (1,1,1,1,1,...)
  Rosenbrock( int n )
    : ccdl::minifcn( n ), neval(0)
  {}

  Rosenbrock()
    : ccdl::minifcn( 2 ), neval(0)
  {
    refx[0] = -1.0;
    refx[1] =  0.5;
  }
  
  virtual ~Rosenbrock() {}
  
  virtual double operator() ( double const alp, double const * s );

  virtual double operator() ( double const alp, double const * s, 
			      double * dydx );

  int neval;
};

inline
double Rosenbrock::operator() ( double const alp, double const * s )
{
  std::vector<double> x( refx );
  for ( int i=0; i<nvar; ++i )
    x[i] += alp*s[i];
  double f = 0.;
  for ( int i=0; i < nvar-1; ++i )
    f += 100. * std::pow( x[i+1] - x[i]*x[i], 2 ) + std::pow( x[i]-1., 2 );
  return f;
}

inline
double Rosenbrock::operator() ( double const alp, double const * s, double * g )
{
  std::vector<double> x( refx );
  for ( int i=0; i<nvar; ++i )
    x[i] += alp*s[i];
  std::fill( g, g+nvar, 0. );
  double f = 0.;
  for ( int i=0; i < nvar-1; ++i )
    {
      double a = x[i+1] - x[i]*x[i];
      double b = x[i]-1.;
      f += 100. * a*a + b*b;
      double dfda = 200. * a;
      double dfdb = 2 * b;
      g[i]   += dfda * (-2.*x[i]) + dfdb;
      g[i+1] += dfda;
    };

  neval++;
  std::printf("Rosenbrock %6i %23.14e\n",neval,f);
  return f;
}


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


struct Sphere : public ccdl::minifcn
{
  // minimum @ (0,0,0,0....)
  Sphere( int n )
    : ccdl::minifcn( n ), neval(0)
  {}

  Sphere()
    : ccdl::minifcn( 2 ), neval(0)
  {
    refx[0] = -1.0;
    refx[1] =  0.5;
  }
  
  virtual ~Sphere() {}
  
  virtual double operator() ( double const alp, double const * s );

  virtual double operator() ( double const alp, double const * s, 
			      double * dydx );

  int neval;
};

inline
double Sphere::operator() ( double const alp, double const * s )
{
  std::vector<double> x( refx );
  for ( int i=0; i<nvar; ++i )
    x[i] += alp*s[i];
  double f = 0.;
  for ( int i=0; i < nvar-1; ++i )
    f += x[i]*x[i];
  return f;
}

inline
double Sphere::operator() ( double const alp, double const * s, double * g )
{
  std::vector<double> x( refx );
  for ( int i=0; i<nvar; ++i )
    x[i] += alp*s[i];
  std::fill( g, g+nvar, 0. );
  double f = 0.;
  for ( int i=0; i < nvar; ++i )
    {
      f += x[i]*x[i];
      g[i] += 2. * x[i];
    };

  neval++;
  std::printf("Sphere %6i %23.14e\n",neval,f);
  return f;
}


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


struct Beale : public ccdl::minifcn
{
  // minimum @ (3,-0.5)

  Beale()
    : ccdl::minifcn( 2 ), neval(0)
  {
    refx[0] =  -0.6;
    refx[1] =  0.5;
  }
  
  virtual ~Beale() {}
  
  virtual double operator() ( double const alp, double const * s );

  virtual double operator() ( double const alp, double const * s, 
			      double * dydx );

  int neval;
};

inline
double Beale::operator() ( double const alp, double const * s )
{
  double x = refx[0] + alp*s[0];
  double y = refx[1] + alp*s[1];

  double a = 1.5 - x + x*y;
  double b = 2.25 - x + x*y*y;
  double c = 2.625 - x + x*y*y*y;
  double f = a*a+b*b+c*c;


  if ( std::abs(x) > 4.5 )
    {
      double ax = std::abs(x);
      f += std::pow( ax-4.5, 2 );
    }

  if ( std::abs(y) > 4.5 )
    {
      double ay = std::abs(y);
      f += std::pow( ay-4.5, 2 );
    }
  return f;
}

inline
double Beale::operator() ( double const alp, double const * s, double * g )
{
  double x = refx[0] + alp*s[0];
  double y = refx[1] + alp*s[1];

  double a = 1.5 - x + x*y;
  double b = 2.25 - x + x*y*y;
  double c = 2.625 - x + x*y*y*y;
  double f = a*a+b*b+c*c;
  g[0] = 2*( a*(y-1) + b*(y*y-1) + c*(y*y*y-1) );
  g[1] = 2*( a*(x) + b*(2*x*y) + c*(3*x*y*y) );

  if ( std::abs(x) > 4.5 )
    {
      double ax = std::abs(x);
      f += std::pow( ax-4.5, 2 );
      g[0] += 2. * ( ax-4.5 ) * ( x/ax );
    }

  if ( std::abs(y) > 4.5 )
    {
      double ay = std::abs(y);
      f += std::pow( ay-4.5, 2 );
      g[1] += 2. * ( ay-4.5 ) * ( y/ay );
    }


  neval++;
  std::printf("Beale %6i %23.14e\n",neval,f);
  return f;
}


int main()
{
  Beale fcn;

  fcn.testgrd();

  ccdl::minimize( &fcn );

  std::printf("minimum @ ");
  for ( int i=0; i<fcn.size(); ++i )
    std::printf("%20.10e",fcn[i]);
  std::printf("\n");

  return 0;
}

 */












int ccdl::minimize( ccdl::minifcn * pfcn, ccdl::mini::options opts )
{
  ccdl::minifcn & fcn = *pfcn;
  int const n = fcn.size();
  //int const inc = 1;
  
#define FMTE(a) std::setw(11) << std::setprecision(2) << std::scientific << (a) << " "
#define FMTF(a) std::setw(22) << std::setprecision(12) << std::fixed << (a) << " "
#define PRT if ( opts.verbosity ) cout << "GEOMOPT " << std::setw(5) << iter++ << FMTF(step.f) << FMTE(step.gmax) << FMTE(step.grms) << "\n";

  std::ostream & cout = *(opts.ostr);
  int iter=0;

  ccdl::mini::BfgsStep step( pfcn ), prev( pfcn );
  std::vector<double> snew(n,0.);
  step.cpt_fg( 1. );
  PRT;


  prev = step;
  if ( step.gmax < 1.e-15 )
    {
      //return 0;
      
      prev.s[0] += 5.e-5;
      prev.cpt_fg( 1. );
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
      
      if ( iloop < opts.nsd  )
	{
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

      //double sdir = ccdl::v_dot_v( n, step.s.data(), step.g.data() );
      //cout << FMTE(sdir) << "\n";

      prev = step;

      // don't walk uphill
      double sgrd = ccdl::v_dot_v( n, step.g.data(), step.s.data() );
      if ( sgrd > 0. ) 
	{
	  step.steepdesc_direction();
	  alp0 = 1. / ( 1. + std::sqrt( ccdl::v_dot_v( n, step.s.data(), step.s.data() ) ) );
	}
      alp = step.linemin( alp0 );
      PRT;

      //std::printf("GREP %5i %20.10f %11.2e %11.2e (alp0=%15.5f alp=%15.5f)\n\n\n\n\n\n\n",iter,step.f,step.gmax,step.grms,alp0,alp);

      if ( (step.gmax < opts.gmax_tol and 
	    step.grms < opts.grms_tol) ) break;

      if ( step.f-prev.f == 0 )
	{
	  step.steepdesc_direction();
	  alp0 = 1. / ( 1. + std::sqrt( ccdl::v_dot_v( n, step.s.data(), step.s.data() ) ) );
	  alp = step.linemin( alp0 );
	  if ( step.f-prev.f == 0 ) break;
	}

      step.update_hessian( prev );

      alp0 = alp;
      if ( alp0 < 0. ) alp0 = 1. / ( 1. + ccdl::v_dot_v(n,step.g.data(),step.g.data()) );
      alp0 = std::min( alp0, 1. );
      //alp0 = 1.;

    };
  
  return 0;
}
