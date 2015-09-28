#ifndef _mini_step_hpp_
#define _mini_step_hpp_

#include "linemin.hpp"

namespace ccdl
{
  namespace mini
  {

    struct SteepDescStep
    {
      SteepDescStep( int n );

      void set_s( double const * u );

      template <class T>
      void cpt_f( double const alp, T fcn );

      template <class T>
      void cpt_fg( double const alp, T fcn );

      void set_fg( ccdl::mini::OneDPoint const & best );

      template <class FCNVAL>
      double linemin( double alp0, FCNVAL fcn );

      void steepdesc_direction();

      int n;
      double f,grms,gmax;
      std::vector<double> s,p,x,g;
    };



    struct ConjGradStep : public ccdl::mini::SteepDescStep
    {
      ConjGradStep( int n ) : SteepDescStep(n) {}
      void hestenes_stiefel_direction( ConjGradStep const * prev, double alp0 );
      void polak_ribiere_direction( ConjGradStep const * prev, double alp0 );
      void fletcher_reeves_direction( ConjGradStep const * prev, double alp0 );
    };






    // struct RfoFcn
    // {
    //   RfoFcn( int n, double const f, double const * g, double const * h, double const * s )
    //   {
    // 	a=f;
    // 	b = ccdl::v_dot_v(n,g,s);
    // 	std::vector<double> hs(n,0.);
    // 	ccdl::sy_dot_v( n,n, h, s, hs.data() );
    // 	c = 0.5 * ccdl::v_dot_v( n, s, hs.data() );
    //   }

    //   double operator() ( double const * x0, double const alp, double const * s )
    //   {
    // 	double x = x0[0] + alp * s[0];
    // 	double f = a + (b*x + c*x*x) / ( 1. + x*x );
    // 	return f;
    //   }

    //   double operator() ( double const * x0, double const alp, double const * s, double * g )
    //   {
    // 	double x = x0[0] + alp * s[0];
    // 	double num = (b*x + c*x*x);
    // 	double den = ( 1. + x*x );
    // 	double den2 = den*den;
    // 	double f = a + num / den;
    // 	g[0] = (b + 2.*c*x)/den - (num/den2) * (2.*x);
    // 	return f;
    //   }
    //   double a,b,c;
    // };
    //
    // std::pair<double,double> RfoStepLen
    // ( int n, 
    //   double const f, 
    //   double const * g, 
    //   double const * h, 
    //   double const * s )
    // {
    //   RfoFcn rfo( n, f, g, h, s );

    //   std::vector<double> x0(1,0.), g0(1,0.), s0(1,1.);
    //   double f0 = rfo( x0.data(), 0., s0.data(), g0.data() );
    //   OneDCurve curve = BracketPolyBisection_Grad( 1, f0, x0.data(), g0.data(), 1., s0.data(), rfo );
    //   return std::make_pair( curve[0].x, curve[0].y );
    // }



    struct BfgsStep : public ccdl::mini::ConjGradStep
    {
      BfgsStep( int n );

      void update_hessian( ccdl::mini::BfgsStep const & prev );

      void bfgs_direction();

      std::vector<double> h;
    };

  }
}


template <class T>
void ccdl::mini::SteepDescStep::cpt_f( double const alp, T fcn ) 
{ 
  f = fcn( x.data(), alp, s.data() );
  for ( int i=0;i<n;++i )
    { 
      p[i]  = alp*s[i];
      x[i] += p[i];
    };
}


template <class T>
void ccdl::mini::SteepDescStep::cpt_fg( double const alp, T fcn ) 
{ 
  std::fill(g.data(),g.data()+n,0.); 
  f = fcn( x.data(), alp, s.data(), g.data() );
  for ( int i=0;i<n;++i )
    { 
      p[i]  = alp*s[i];
      x[i] += p[i];
    };
  grms = 0.;
  gmax = 0.;
  for ( int i=0; i<n; ++i )
    {
      grms += g[i]*g[i];
      gmax = std::max( gmax, std::abs(g[i]) );
    }
  grms = std::sqrt( grms / n );
}


template <class FCNVAL>
double ccdl::mini::SteepDescStep::linemin( double alp0, FCNVAL fcn )
{
  ccdl::mini::LINET type = ccdl::mini::NOGRAD_NOPOLY;
  ccdl::mini::OneDPoint pt = ccdl::mini::linemin
    ( type, n, f, x.data(), g.data(), alp0, s.data(), fcn );
  set_fg( pt );
  return pt.x;
}


#endif
