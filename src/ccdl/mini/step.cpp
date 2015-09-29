#include "step.hpp"
#include "../bmath.hpp"
#include "linemin.hpp"

inline double vdot( std::vector<double> const & a, std::vector<double> const & b )
{
  return ccdl::v_dot_v( a.size(), a.data(), b.data() );
}



ccdl::mini::SteepDescStep::SteepDescStep( ccdl::minifcn * pfcn ) 
  : pfcn(pfcn), n(pfcn->size()), f(0.), 
    grms(0.), gmax(0.), 
    s(n,0.), p(n,0.), 
    x(n,0.), g(n,0.) 
{}



void ccdl::mini::SteepDescStep::cpt_f( double const alp ) 
{ 
  ccdl::minifcn & fcn = *pfcn;
  f = fcn( alp, s.data() );
  for ( int i=0;i<n;++i )
    { 
      p[i]  = alp*s[i];
      x[i] += p[i];
    };
}


void ccdl::mini::SteepDescStep::cpt_fg( double const alp ) 
{ 
  ccdl::minifcn & fcn = *pfcn;

  std::fill(g.data(),g.data()+n,0.); 
  f = fcn( alp, s.data(), g.data() );
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


double ccdl::mini::SteepDescStep::linemin( double alp0 )
{
  ccdl::mini::LINET type = ccdl::mini::GRAD_POLY;
  ccdl::mini::OneDPoint pt = ccdl::mini::linemin
    ( type, n, f, g.data(), alp0, s.data(), pfcn );
  pfcn->SetNewReferencePt( pt.x, s.data() );
  set_fg( pt );
  return pt.x;
}





void ccdl::mini::SteepDescStep::set_s( double const * u ) 
{ 
  std::copy( u,u+n,s.data() ); 
}



void ccdl::mini::SteepDescStep::set_fg( OneDPoint const & best )
{
  if ( ! best.grd.size() )
    {
      std::cerr << "set_fg( OneDPoint const & best ); best has no grd!\n";
      std::abort();
    }
  g = best.grd;
  f = best.y;
  double alp = best.x;
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



void ccdl::mini::SteepDescStep::steepdesc_direction()
{
  for ( int i=0; i<n; ++i ) s[i] = -g[i];
}

void ccdl::mini::ConjGradStep::hestenes_stiefel_direction
( ccdl::mini::ConjGradStep const * prev, double alp0 )
{
  std::vector<double> snew(n,0.);
  for ( int i=0; i<n; ++i ) snew[i] = -g[i];
  double beta = 0.;
  // Fletcher-Reeves
  //beta = vdot( g, g ) / vdot( prev->g, prev->g );
  // Polak-Ribiere
  //beta = ( vdot( g, g ) - vdot( g, prev->g ) ) / vdot( prev->g, prev->g );
  // Hestenes-Stiefel
  beta = - ( vdot(g,g)-vdot(g,prev->g) ) / ( vdot(s,prev->g)-vdot(s,g) );
  beta = std::max( 0., beta );
  if ( alp0 < 1.e-5 ) beta = 0.;
  for ( int i=0; i<n; ++i )
    s[i] = snew[i] + beta * prev->s[i];
}

void ccdl::mini::ConjGradStep::polak_ribiere_direction
( ccdl::mini::ConjGradStep const * prev, double alp0 )
{
  std::vector<double> snew(n,0.);
  for ( int i=0; i<n; ++i ) snew[i] = -g[i];
  double beta = 0.;
  // Fletcher-Reeves
  //beta = vdot( g, g ) / vdot( prev->g, prev->g );
  // Polak-Ribiere
  beta = ( vdot( g, g ) - vdot( g, prev->g ) ) / vdot( prev->g, prev->g );
  // Hestenes-Stiefel
  //beta = - ( vdot(g,g)-vdot(g,prev->g) ) / ( vdot(s,prev->g)-vdot(s,g) );
  beta = std::max( 0., beta );
  if ( alp0 < 1.e-5 ) beta = 0.;
  for ( int i=0; i<n; ++i )
    s[i] = snew[i] + beta * prev->s[i];
}

void ccdl::mini::ConjGradStep::fletcher_reeves_direction
( ccdl::mini::ConjGradStep const * prev, double alp0 )
{
  std::vector<double> snew(n,0.);
  for ( int i=0; i<n; ++i ) snew[i] = -g[i];
  double beta = 0.;
  // Fletcher-Reeves
  beta = vdot( g, g ) / vdot( prev->g, prev->g );
  // Polak-Ribiere
  //beta = ( vdot( g, g ) - vdot( g, prev->g ) ) / vdot( prev->g, prev->g );
  // Hestenes-Stiefel
  //beta = - ( vdot(g,g)-vdot(g,prev->g) ) / ( vdot(s,prev->g)-vdot(s,g) );
  beta = std::max( 0., beta );
  if ( alp0 < 1.e-5 ) beta = 0.;
  for ( int i=0; i<n; ++i )
    s[i] = snew[i] + beta * prev->s[i];
}


ccdl::mini::BfgsStep::BfgsStep( ccdl::minifcn * pfcn ) 
  : ccdl::mini::ConjGradStep( pfcn ), h( pfcn->size() * pfcn->size(), 0. )
{
  for ( int i=0; i<n; ++i )
    h[i+i*n] = 1.;
}

  
void ccdl::mini::BfgsStep::update_hessian( BfgsStep const & prev )
{
  std::vector<double> Hp(n,0.);
  ccdl::sy_dot_v( n, n, h.data(), p.data(), Hp.data() );
  double pHp = vdot( p, Hp );
  std::vector<double> y( g );
  for ( int i=0; i<n; ++i ) y[i] -= prev.g[i];
  double yp = vdot( y, p );
  if ( yp > 0. )
    {
      ccdl::v_tensor_v( 1./yp, n, y.data(), n, y.data(), h.data() );
      ccdl::v_tensor_v( -1./pHp, n, Hp.data(), n, Hp.data(), h.data() );
    };
}

void ccdl::mini::BfgsStep::bfgs_direction()
{
  for ( int i=0; i<n; ++i ) p[i] = -g[i];
  ccdl::dsysv( n, h.data(), p.data(), s.data() );
  std::fill( p.data(), p.data() + n, 0. );
}
