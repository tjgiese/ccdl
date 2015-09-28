#include "OneDCurve.hpp"
#include "../bmath.hpp"
#include <cmath>


ccdl::mini::OneDPoint::OneDPoint
( double X, double Y, int N, double const * S, double const * GRD )
  : x(X), y(Y), wy(1.), wg(1.), grd( GRD, GRD + N )
{
  g = ccdl::v_dot_v( N, S, GRD );
}

void ccdl::mini::OneDPoint::SetBoltzmannWeight( double y0, double beta )
{
  wy = std::exp( - beta * (y-y0) );
  wg = std::min( wg, wy );
}



bool ccdl::mini::OneDCurve::g_is_bracketed() 
{ 
  sort_by_g(); 
  bool ok = ( data[0].g * data.back().g < 0. ); 
  sort_by_y();
  return ok;
}


void ccdl::mini::OneDCurve::push_back( double x, double y )
{
  data.push_back( ccdl::mini::OneDPoint( x, y ) );
  std::sort( data.begin(), data.end(), lt_y );
}

void ccdl::mini::OneDCurve::push_back( double x, double y, double g )
{
  data.push_back( ccdl::mini::OneDPoint( x, y, g ) );
  std::sort( data.begin(), data.end(), lt_y );
}

void ccdl::mini::OneDCurve::push_back( double X, double Y, int N, double const * S, double const * GRD )
{
  data.push_back( ccdl::mini::OneDPoint( X, Y, N, S, GRD ) );
  std::sort( data.begin(), data.end(), lt_y );
}


void ccdl::mini::OneDCurve::SetBoltzmannWeight( double beta )
{
  if ( ! size() ) return;
  std::sort( data.begin(), data.end(), lt_y );
  double y0 = data[0].y;
  for ( std::vector< ccdl::mini::OneDPoint >::iterator
	  p = data.begin(), pend = data.end(); 
	p != pend; ++p )
    p->SetBoltzmannWeight( y0, beta ); 
}


ccdl::polynomial ccdl::mini::OneDCurve::fit( int npoly, int icon, double conval ) const
{
  if ( icon > npoly )
    {
      std::cerr << "OneDCurve::fit icon > npoly ; programming error" << std::endl;
      std::abort();
    }
  std::vector<double> x,y,g,wy,wg;
  for ( std::vector< ccdl::mini::OneDPoint >::const_iterator
	  p = data.begin(), pend = data.end(); 
	p != pend; ++p )
    {
      x.push_back( p->x );
      y.push_back( p->y );
      g.push_back( p->g );
      wy.push_back( p->wy );
      wg.push_back( p->wg );
    };
  ccdl::polynomial poly( npoly );
  if ( icon > -1 )
    poly.wfit( x.size(), x.data(), y.data(), g.data(),
	       wy.data(), wg.data(),
	       icon, conval );
  else
    poly.wfit( x.size(), x.data(), y.data(), g.data(),
	       wy.data(), wg.data() );

  //print( std::cerr );

  std::sort( x.begin(), x.end() );
  double xlo = x[0];
  double xhi = x.back();
  double mloc = 0.;
  if ( poly.minimum_loc( mloc ) )
    {
      //std::cout << "minimum   " << std::setw(20) << std::setprecision(10) << std::scientific << mloc << "\n";
      double d = xhi-xlo;
      xlo = std::min( xlo, mloc - d );
      xhi = std::max( xhi, mloc + d );
    }
  //else std::cout << "no minimum\n";
  //std::cerr << "\n";

  /*
  int npt = 101;
  for ( int i=0; i<npt; ++i )
    {
      double u = xlo + i * (xhi-xlo)/(npt-1);
      double f = poly(u);
      std::cerr << std::setw(20) << std::setprecision(10) << std::scientific << u
		<< std::setw(20) << std::setprecision(10) << std::scientific << f
		<< "\n";

    }
  */

  return poly;
}


ccdl::mini::OneDCurve ccdl::mini::OneDCurve::get_bracketing_curve() const
{
  ccdl::mini::OneDCurve bcurve( *this );
  bcurve.sort_by_y();

  ccdl::mini::OneDPoint cp = bcurve[0];
  if ( cp.grd.size() )
    {
      ccdl::mini::OneDPoint rhs;
      if ( cp.g > 0. )
	{
	  for ( int i=1; i < bcurve.size(); ++i )
	    if ( bcurve[i].g <= 0. )
	      {
		rhs = bcurve[i];
		break;
	      };
	}
      else
	{
	  for ( int i=1; i < bcurve.size(); ++i )
	    if ( bcurve[i].g > 0. )
	      {
		rhs = bcurve[i];
		break;
	      };
	}
      bcurve.resize(2);
      bcurve[0] = cp;
      bcurve[1] = rhs;
    }
  else
    {
      ccdl::mini::OneDPoint lo,hi;
      if ( bcurve[1].x > bcurve[0].x )
	{
	  hi = bcurve[1];
	  lo = cp;
	  for ( int i=2; i < bcurve.size(); ++i )
	    if ( bcurve[i].x < bcurve[0].x )
	      {
		lo = bcurve[i];
		break;
	      }
	}
      else
	{
	  lo = bcurve[1];
	  hi = cp;
	  for ( int i=2; i < bcurve.size(); ++i )
	    if ( bcurve[i].x > bcurve[0].x )
	      {
		hi = bcurve[i];
		break;
	      }
	}
      bcurve.resize(3);
      bcurve[0] = lo;
      bcurve[1] = cp;
      bcurve[2] = hi;
    };

  bcurve.sort_by_x();
  return bcurve;
}

void ccdl::mini::OneDCurve::print( std::ostream & cout ) const
{
  //cout << std::setw(4) << size() << "\n";
  cout << "\n";
#define FMTE(f) std::setw(20) << std::setprecision(10) << std::scientific << (f)
  for ( int i=0; i<size(); ++i )
    cout << FMTE(data[i].x) << FMTE(data[i].y) << "\n";
#undef FMTE
}


