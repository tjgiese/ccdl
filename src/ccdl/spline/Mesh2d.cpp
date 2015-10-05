#include <sstream>
#include <cstdio>
#include <iomanip>
#include "bspline.hpp"
#include "Mesh2d.hpp"
#include "../mini.hpp"
#include "../bmath.hpp"

namespace ccdl
{
  struct Mesh2dOpt : public ccdl::minifcn
  {
    Mesh2dOpt( ccdl::Mesh2d * m );
    
    double operator() ( double const alp, double const * s );
    double operator() ( double const alp, double const * s, double * g );
    
    ccdl::Mesh2d * mesh;
  };
}


ccdl::Mesh2dOpt::Mesh2dOpt( ccdl::Mesh2d * m )
  : ccdl::minifcn(2), mesh(m)
{}

double ccdl::Mesh2dOpt::operator() ( double const alp, double const * s )
{
  return mesh->GetValue( refx[0] + alp*s[0], refx[1] + alp*s[1] ).f;
}

double ccdl::Mesh2dOpt::operator() ( double const alp, double const * s, double * g )
{
  ccdl::Mesh2dValue v = mesh->GetValue( refx[0] + alp*s[0], refx[1] + alp*s[1] );
  g[0] = v.dfdx;
  g[1] = v.dfdy;
  return v.f;
}




namespace ccdl
{
  struct Mesh2dStationaryOpt : public ccdl::minifcn
  {
    Mesh2dStationaryOpt( ccdl::Mesh2d * m );
    
    double operator() ( double const alp, double const * s );
    double operator() ( double const alp, double const * s, double * g );
    
    ccdl::Mesh2d * mesh;
  };
}


ccdl::Mesh2dStationaryOpt::Mesh2dStationaryOpt( ccdl::Mesh2d * m )
  : ccdl::minifcn(2), mesh(m)
{}


double ccdl::Mesh2dStationaryOpt::operator() ( double const alp, double const * s )
{
  ccdl::Mesh2dValue v( mesh->GetValue( refx[0] + alp*s[0], refx[1] + alp*s[1] ) );
  return v.dfdx*v.dfdx + v.dfdy*v.dfdy;
}

double ccdl::Mesh2dStationaryOpt::operator() ( double const alp, double const * s, double * g )
{
  ccdl::Mesh2dHessian v = mesh->GetHessian( refx[0] + alp*s[0], refx[1] + alp*s[1] );
  g[0] = 2. * v.dfdx * v.h[0] + 2. * v.dfdy * v.h[1];
  g[1] = 2. * v.dfdx * v.h[2] + 2. * v.dfdy * v.h[3];
  return v.dfdx*v.dfdx + v.dfdy*v.dfdy;
}











ccdl::Mesh2d::Mesh2d
( int const nx, int const ny, 
  double const lx, double const ly,
  bool periodic )
  : nx(nx), ny(ny), 
    lx(lx), ly(ly),
    maxx(lx),maxy(ly),
    delx(lx/nx),dely(ly/nx),
    xlow(0.),
    ylow(0.),
    order(4), 
    walled(false),
    periodic(true),
    refdata(nx*ny,0.), 
    data(nx*ny,0.)
{
  SetPeriodic( periodic );
}

ccdl::Mesh2d::Mesh2d
( int const nx, int const ny, 
  double const Minx, double const Maxx,
  double const Miny, double const Maxy,
  bool periodic )
  : nx(nx), ny(ny), 
    lx(Maxx-Minx), ly(Maxy-Miny),
    maxx(Maxx-Minx),maxy(Maxy-Miny),
    delx((Maxx-Minx)/nx),dely((Maxy-Miny)/nx),
    xlow(Minx),
    ylow(Miny),
    order(4), 
    walled(false),
    periodic(true),
    refdata(nx*ny,0.), 
    data(nx*ny,0.)
{
  SetPeriodic( periodic );
}


ccdl::Mesh2d::Mesh2d
( std::istream & cin,
  bool periodic )
  : nx(0), ny(0), 
    lx(0), ly(0),
    maxx(0),maxy(0),
    delx(0),dely(0),
    xlow(0),
    ylow(0),
    order(4), 
    walled(false),
    periodic(true)
{
  std::string line;
  std::stringstream str;
  //int num_empty_lines = 0;
  bool prev_line_was_empty = true;
  std::vector<double> d;
  while( std::getline(cin,line) )
    {
      if ( line.empty() ) 
	{
	  prev_line_was_empty = true;
	  continue;
	}
      if ( prev_line_was_empty ) 
	{ 
	  nx = 0;
	  ++ny;
	};
      ++nx;
      prev_line_was_empty = false;

      str.str( line );
      str.clear();
      double x,y,z;
      str >> x >> y >> z;
      xlow = std::min(xlow,x);
      ylow = std::min(ylow,y);
      maxx = std::max(maxx,x);
      maxy = std::max(maxy,y);
      d.push_back( x );
      d.push_back( y );
      d.push_back( z );
    };

  int n = d.size()/3;
  if ( nx*ny != (int)(d.size()/3) )
    {
      std::cerr << "ccdl::Mesh2d read " << n 
		<< " lines of data but counted " << nx << " x " << ny
		<< " elements based on the number of blank spaces\n"
		<< "Do you have extra spaces in your file?\n";
      std::abort(); // this is not recoverable
    };

  lx = maxx-xlow;
  ly = maxy-ylow;
  if ( periodic )
    {
      lx *= (double)(nx) / (double)(nx-1);
      ly *= (double)(ny) / (double)(ny-1);
    };
  SetPeriodic( periodic );

  refdata.resize( n );
  data.resize( n );
  for ( int i=0; i<nx; ++i )
    for ( int j=0; j<ny; ++j )
      refdata[j+i*ny] = d[2+(i+j*nx)*3]; // gnuplot uses transpose
  


  for ( int i=0; i<nx; ++i )
    for ( int j=0; j<ny; ++j )
      {
	int k = i+j*nx; // gnuplot uses transpose
	double x = GetX(i);
	double y = GetY(j);
	if ( std::abs( x - d[0+k*3] ) > 1.e-7 )
	  {
	    std::cerr << "ccdl::Mesh2d data point " << k+1
		      << " x-value is " 
		      << std::scientific << std::setprecision(15) 
		      << std::setw(24) << d[0+k*3]
		      << " but I was expecting "
		      << std::scientific << std::setprecision(15) 
		      << std::setw(24) << x
		      << std::endl;
	    std::cerr << "Are the data points listed in ascending order?"
		      << std::endl;
	    std::abort(); // lazy -- I could potentially figure this out
	  };
	if ( std::abs( y - d[1+k*3] ) > 1.e-7 )
	  {
	    std::cerr << "ccdl::Mesh2d data point " << k+1
		      << " y-value is " 
		      << std::scientific << std::setprecision(15) 
		      << std::setw(24) << d[1+k*3]
		      << " but I was expecting "
		      << std::scientific << std::setprecision(15) 
		      << std::setw(24) << y
		      << std::endl;
	    std::cerr << "Are the data points listed in ascending order?"
		      << std::endl;
	    std::abort(); // lazy -- I could potentially figure this out
	  };

      }
  BsplineTransform();

}


void ccdl::Mesh2d::Add( double shift )
{
  for ( std::size_t i=0; i < refdata.size(); ++i )
    refdata[i] += shift;
  BsplineTransform();
}

void ccdl::Mesh2d::Write( std::ostream & cout ) const
{

#define FF(a) std::scientific << std::setprecision(15) << std::setw(24) << (a)

  for ( int j=0; j<ny; ++j )
    {
      for ( int i=0; i<nx; ++i )
	cout << FF(GetX(i)) << FF(GetY(j)) << FF(refdata[j+i*ny]) << "\n";
      cout << "\n";
    }

#undef FF
}



// Discretizes the range
// If periodic, then the last point is not explicitly given
// ...this is the fftw3 consistent way of doing it...
void ccdl::Mesh2d::SetPeriodic( bool logical )
{
  periodic = logical;
  if ( periodic )
    {
      maxx = lx;
      maxy = ly;
    }
  else
    {
      maxx = lx * (nx)/((double)(nx-1));
      maxy = ly * (ny)/((double)(ny-1));
    };
  delx = maxx/nx;
  dely = maxy/ny;
  //delx = lx/nx;
  //dely = ly/ny;
}

void ccdl::Mesh2d::SetWalled( bool logical )
{
  walled = logical;
}

int ccdl::Mesh2d::GetSizeX() const { return nx; }

int ccdl::Mesh2d::GetSizeY() const { return ny; }

double ccdl::Mesh2d::GetX( int i ) const
{
  return GetIndex( i, nx ) * delx + xlow;
}

double ccdl::Mesh2d::GetY( int i ) const
{
  return GetIndex( i, ny ) * dely + ylow;
}



double & ccdl::Mesh2d::operator() ( int ix, int iy )       
{ 
  ix = GetIndex( ix, nx );
  iy = GetIndex( iy, nx );
  return refdata[iy+ix*ny]; 
}

double   ccdl::Mesh2d::operator() ( int ix, int iy ) const 
{ 
  ix = GetIndex( ix, nx );
  iy = GetIndex( iy, nx );
  return refdata[iy+ix*ny]; 
}


void ccdl::Mesh2d::BsplineTransform( int bspline_order )
{
  order = bspline_order;
  data = refdata;
  ccdl::bspline_renormalize( maxx, maxy, nx, ny, order, data.data() );
}


ccdl::Mesh2dValue ccdl::Mesh2d::GetValue( double x, double y ) const
{


  x -= xlow;
  y -= ylow;

  double xin = x;
  double yin = y;
  if ( walled ) // elastic collisions with walls
    {
      double xold = x;
      for ( int i=0; i<100; ++i )
	{
	  if ( x > maxx ) x = maxx - ( x-maxx );
	  if ( x < -0. ) x = -x;
	  if ( x == xold ) break;
	  xold = x;
	};
      double yold = y;
      for ( int i=0; i<100; ++i )
	{
	  if ( y > maxy ) y = maxy - ( y-maxy );
	  if ( y < -0. ) y = -y;
	  if ( y == yold ) break;
	  yold = y;
	};
    }
  else if ( ! periodic )
    {
      if ( x > maxx ) x = maxx;
      if ( x < 0. ) x = 0.;

      if ( y > maxy ) y = maxy;
      if ( y < 0. ) y = 0.;
    };

  std::vector<double> wx( order, 0. ), wy( order, 0. );
  std::vector<double> dx( order, 0. ), dy( order, 0. );
  std::vector<int> gx( order, 0 ), gy( order, 0 );
  ccdl::bspline_periodic_deriv( x, maxx, nx, order, 
				wx.data(), dx.data(), gx.data() );
  ccdl::bspline_periodic_deriv( y, maxy, ny, order, 
				wy.data(), dy.data(), gy.data() );




  ccdl::Mesh2dValue v;
  v.x = x + xlow;
  v.y = y + ylow;
  for ( int ii=0; ii<order; ++ii )
    {
      int i = gx[ii];
      for ( int jj=0; jj<order; ++jj )
	{
	  int j = gy[jj];
	  double u = data[ j+i*ny ];
	  v.f    += wx[ii] * wy[jj] * u;
	  v.dfdx += dx[ii] * wy[jj] * u;
	  v.dfdy += wx[ii] * dy[jj] * u;
	};
    }
  

  if ( ! periodic and xin != x )
    v.dfdx = 0.;
  if ( ! periodic and yin != y )
    v.dfdy = 0.;

  x = xin;
  y = yin;

  if ( ! periodic )
    {
      double k = 100.;
      double DEL = -1.e-8;
      double d = 0.;
      if ( x > maxx-DEL ) 
	{
	  d = x - (maxx-DEL);
	  v.f    += k * d*d;
	  v.dfdx += 2.*k * d;
	}
      else if ( x < DEL )
	{
	  d = DEL - x;
	  v.f    += k * d*d;
	  v.dfdx += 2.*k * d * (-1);
	}
      if ( y > maxy-DEL ) 
	{
	  d = y - (maxy-DEL);
	  v.f    += k * d*d;
	  v.dfdy += 2.*k * d;
	}
      else if ( y < DEL )
	{
	  d = DEL - y;
	  v.f    += k * d*d;
	  v.dfdy += 2.*k * d * (-1);
	}
    }
  

  //std::printf("entered %24.15e %24.15e %24.15e %24.15e %24.15e\n",v.f,x,y,v.dfdx,v.dfdy);

  return v;
}


/*

std::vector< ccdl::Mesh2dValue > ccdl::Mesh2d::GetMinima( double const TOL )
{
  std::vector< ccdl::Mesh2dValue > balls;
  ccdl::Mesh2dOpt o( this );
  ccdl::mini::options options;
  options.grms_tol = TOL;
  options.verbosity = 0;


  std::vector< ccdl::Mesh2dValue > vals(nx*ny);
  for ( int i=0; i<nx; ++i )
    for ( int j=0; j<ny; ++j )
      vals[j+i*ny] = GetValue( GetX(i), GetY(j) );
	
  int const nxn = periodic ? nx : nx-1;
  int const nyn = periodic ? ny : ny-1;
  int const nx0 = periodic ? 0 : 1;
  int const ny0 = periodic ? 0 : 1;

  for ( int i=nx0; i<nxn; ++i )
    for ( int j=ny0; j<nyn; ++j )
      {
	int i1 = GetIndex( i+1, nx );
	int j1 = GetIndex( j+1, ny );

	ccdl::Mesh2dValue & ll = vals[j +i *ny];
	ccdl::Mesh2dValue & lr = vals[j +i1*ny];
	ccdl::Mesh2dValue & ul = vals[j1+i *ny];
	ccdl::Mesh2dValue & ur = vals[j1+i1*ny];

	double g1 = ll.dfdx * lr.dfdx ;
	double g2 = ul.dfdx * ur.dfdx ;
	double g3 = ll.dfdy * ul.dfdy ;
	double g4 = lr.dfdy * ur.dfdy ;

	// if ( ll.x > -2.27 and ll.x < -2.2 and 
	//      ll.y > -0.94 and ll.y < -0.85 )
	//   {
	//     std::printf("KK  %18.10f %18.10f  %11.2e %11.2e %11.2e %11.2e\n",
	// 		ll.x,
	// 		ll.y,
	// 		g1,g2,g3,g4 );
	//     std::printf("ll  %18.10f %18.10f  %11.2e %11.2e\n",
	// 		ll.x, ll.y, ll.dfdx, ll.dfdy );
	//     std::printf("lr  %18.10f %18.10f  %11.2e %11.2e\n",
	// 		lr.x, lr.y, lr.dfdx, lr.dfdy );
	//     std::printf("ul  %18.10f %18.10f  %11.2e %11.2e\n",
	// 		ul.x, ul.y, ul.dfdx, ul.dfdy );
	//     std::printf("ur  %18.10f %18.10f  %11.2e %11.2e\n",
	// 		ur.x, ur.y, ur.dfdx, ur.dfdy );
	//     std::printf("\n");
	//   }

	int nsatisfied = (g1<=0) + (g2<=0) + (g3<=0) + (g4<=0);

	if ( nsatisfied > 2 )
	  {
	    // std::printf("exr %18.10f %18.10f  %13.4e %13.4e\n",
	    // 		ll.x + 0.5 * delx,
	    // 		ll.y + 0.5 * dely,
	    // 		ll.dfdx, ll.dfdy );
	    // there is an extremum in this quadrant
	    if ( (ll.dfdx < 0 and ll.dfdy < 0) )
	      {
		// this is a minimum

		o[0] = ll.x + 0.5 * delx;
		o[1] = ll.y + 0.5 * dely;
		ccdl::minimize( &o, options );
		balls.push_back( GetValue( o[0], o[1] ) );


		ccdl::Mesh2dValue lmin = ll;
		if ( lr.f < ll.f ) lmin = lr;
		if ( ur.f < ll.f ) lmin = ur;
		if ( ul.f < ll.f ) lmin = ul;

		o[0] = lmin.x;
		o[1] = lmin.y;
		ccdl::minimize( &o, options );
		balls.push_back( GetValue( o[0], o[1] ) );

	      }
	  }
      }


  std::sort( balls.begin(), balls.end(), ccdl::sort_value_by_position );

  typedef std::vector< ccdl::Mesh2dValue >::iterator iter;
  for ( iter p = balls.begin(); p != balls.end(); ++p )
    for ( iter q = p+1; q != balls.end(); )
      {
	double dx = std::abs(p->x - q->x);
	double dy = std::abs(p->y - q->y);
	if ( dx < 0.01 * lx and
	     dy < 0.01 * ly )
	  q = balls.erase( q );
	else
	  ++q;
      }
  // for ( iter p = balls.begin(); p != balls.end(); ++p )
  //   {
  //     double h = GetLowestFreq( p->x , p->y );
  //     std::printf("extremum %20.12f %20.12f  %20.10e %13.4e\n",
  // 		  p->x, p->y, p->f, h );
  //   }
  if ( ! periodic )
    for ( iter p = balls.begin(); p != balls.end(); )
      if ( std::abs( p->x - xlow ) / lx < 0.02 or
	   std::abs( p->y - ylow ) / lx < 0.02 or
	   std::abs( p->x - (xlow+lx) ) / lx < 0.02 or
	   std::abs( p->y - (ylow+ly) ) / ly < 0.02 )
	p = balls.erase( p );
      else
	++p;

  for ( iter p = balls.begin(); p != balls.end(); )
    if ( GetLowestFreq( p->x, p->y ) < 0.001 )
      p = balls.erase( p );
    else
      ++p;


  return balls;
}



std::vector< ccdl::Mesh2dValue > ccdl::Mesh2d::GetMaxima( double const TOL )
{
  int nxy = nx*ny;
  for ( int k=0; k<nxy; ++k )
    refdata[k] *= -1;
  BsplineTransform();
  std::vector< ccdl::Mesh2dValue > balls( GetMinima(TOL) );
  for ( int k=0; k<nxy; ++k )
    refdata[k] *= -1;
  BsplineTransform();
  typedef std::vector< ccdl::Mesh2dValue >::iterator iter;
  for ( iter p = balls.begin(); p != balls.end(); ++p )
    {
      p->f *= -1;
      p->dfdx *= -1;
      p->dfdy *= -1;
    }
  return balls;
}









std::vector< ccdl::Mesh2dHessian > ccdl::Mesh2d::GetStationaryPts( double const TOL )
{
  std::vector< ccdl::Mesh2dHessian > balls;
  ccdl::Mesh2dStationaryOpt o( this );
  ccdl::mini::options options;
  options.grms_tol = TOL;
  options.verbosity = 0;


  std::vector< ccdl::Mesh2dValue > vals(nx*ny);
  for ( int i=0; i<nx; ++i )
    for ( int j=0; j<ny; ++j )
      vals[j+i*ny] = GetValue( GetX(i), GetY(j) );
	
  int const nxn = periodic ? nx : nx-1;
  int const nyn = periodic ? ny : ny-1;
  int const nx0 = periodic ? 0 : 1;
  int const ny0 = periodic ? 0 : 1;

  for ( int i=nx0; i<nxn; ++i )
    for ( int j=ny0; j<nyn; ++j )
      {
	int i1 = GetIndex( i+1, nx );
	int j1 = GetIndex( j+1, ny );

	ccdl::Mesh2dValue & ll = vals[j +i *ny];
	ccdl::Mesh2dValue & lr = vals[j +i1*ny];
	ccdl::Mesh2dValue & ul = vals[j1+i *ny];
	ccdl::Mesh2dValue & ur = vals[j1+i1*ny];

	double l2r[2] = { delx, 0. };
	double u2l[2] = { 0., dely };
	double ul2lr[2] = { delx, -dely };
	double ll2ur[2] = { delx,  dely };

	double g1 = ( l2r[0]*ll.dfdx + l2r[1]*ll.dfdy ) * ( l2r[0]*lr.dfdx + l2r[1]*lr.dfdy );
	double g2 = ( u2l[0]*ul.dfdx + u2l[1]*ul.dfdy ) * ( u2l[0]*ll.dfdx + u2l[1]*ll.dfdy );
	double g3 = ( ll2ur[0]*ll.dfdx + ll2ur[1]*ll.dfdy ) * ( ll2ur[0]*ur.dfdx + ll2ur[1]*ur.dfdy );
	double g4 = ( ul2lr[0]*ul.dfdx + ul2lr[1]*ul.dfdy ) * ( ul2lr[0]*lr.dfdx + ul2lr[1]*lr.dfdy );


	// if ( ll.x > -2.05 and ll.x < -1.9 and 
	//      ll.y > -0.6 and ll.y < -0.4 )
	//   {
	//     std::printf("KK  %18.10f %18.10f  %11.2e %11.2e %11.2e %11.2e\n",
	// 		ll.x,
	// 		ll.y,
	// 		g1,g2,g3,g4 );
	//     std::printf("ll  %18.10f %18.10f  %11.2e %11.2e\n",
	// 		ll.x, ll.y, ll.dfdx, ll.dfdy );
	//     std::printf("lr  %18.10f %18.10f  %11.2e %11.2e\n",
	// 		lr.x, lr.y, lr.dfdx, lr.dfdy );
	//     std::printf("ul  %18.10f %18.10f  %11.2e %11.2e\n",
	// 		ul.x, ul.y, ul.dfdx, ul.dfdy );
	//     std::printf("ur  %18.10f %18.10f  %11.2e %11.2e\n",
	// 		ur.x, ur.y, ur.dfdx, ur.dfdy );
	//     std::printf("\n");
	//   }

	int nsatisfied = (g1<=0) + (g2<=0) + (g3<=0) + (g4<=0);

	if ( nsatisfied > 2 )
	  {
	     // std::printf("exr %18.10f %18.10f  %13.4e %13.4e\n",
	     // 		ll.x + 0.5 * delx,
	     // 		ll.y + 0.5 * dely,
	     // 		ll.dfdx, ll.dfdy );
	    // there is an extremum in this quadrant
	    //if ( (ll.dfdx < 0 and ll.dfdy < 0) or (  )
	      {

		o[0] = ll.x + 0.5 * delx;
		o[1] = ll.y + 0.5 * dely;
		ccdl::minimize( &o, options );
		balls.push_back( GetHessian( o[0], o[1] ) );

		// o[0] = ll.x;
		// o[1] = ll.y;
		// ccdl::minimize( &o, options );
		// balls.push_back( GetHessian( o[0], o[1] ) );

	      }
	  }
      }


  std::sort( balls.begin(), balls.end(), ccdl::sort_hessian_by_position );

  typedef std::vector< ccdl::Mesh2dHessian >::iterator iter;

  for ( iter p = balls.begin(); p != balls.end(); )
    {
      double gnrm = std::sqrt( p->dfdx*p->dfdx + p->dfdy*p->dfdy );
      if ( gnrm > 1.e-4 )
	p = balls.erase( p );
      else
	++p;
    }

  for ( iter p = balls.begin(); p != balls.end(); ++p )
    for ( iter q = p+1; q != balls.end(); )
      {
	double dx = std::abs(p->x - q->x);
	double dy = std::abs(p->y - q->y);
	if ( (dx < 0.01 * lx and dy < 0.01 * ly) )
	  q = balls.erase( q );
	else
	  ++q;
      }

  
  if ( ! periodic )
    for ( iter p = balls.begin(); p != balls.end(); )
      if ( std::abs( p->x - xlow ) / lx < 0.03 or
	   std::abs( p->y - ylow ) / lx < 0.03 or
	   std::abs( p->x - (xlow+lx) ) / lx < 0.03 or
	   std::abs( p->y - (ylow+ly) ) / ly < 0.03 )
	p = balls.erase( p );
      else
	++p;
  

  // for ( iter p = balls.begin(); p != balls.end(); )
  //   if ( GetLowestFreq( p->x, p->y ) < 0.001 )
  //     p = balls.erase( p );
  //   else
  //     ++p;


  return balls;
}



*/









double ccdl::Mesh2d::GetLowestFreq( double x, double y ) const
{
  double DEL = 5.e-5;
  ccdl::Mesh2dValue xh = GetValue( x+DEL, y );
  ccdl::Mesh2dValue xl = GetValue( x-DEL, y );
  ccdl::Mesh2dValue yh = GetValue( x, y+DEL );
  ccdl::Mesh2dValue yl = GetValue( x, y-DEL );

  std::vector<double> H(4, 0.);
  H[0] = (xh.dfdx-xl.dfdx) / ( 2.*DEL );
  H[1] = (xh.dfdy-xl.dfdy + yh.dfdx-yl.dfdx) / ( 4.*DEL );
  H[2] = H[1];
  H[3] = (yh.dfdy-yl.dfdy) / ( 2.*DEL );
  std::vector<double> U(4,0.),E(2,0.);
  ccdl::eigen( 2, H.data(), E.data(), U.data() );
  return E[0];
}



ccdl::Mesh2dHessian ccdl::Mesh2d::GetHessian( double x, double y ) const
{
  double DEL = 5.e-5;

  ccdl::Mesh2dHessian v( GetValue(x,y) );
  ccdl::Mesh2dValue xh = GetValue( x+DEL, y );
  ccdl::Mesh2dValue xl = GetValue( x-DEL, y );
  ccdl::Mesh2dValue yh = GetValue( x, y+DEL );
  ccdl::Mesh2dValue yl = GetValue( x, y-DEL );

  v.h[0] = (xh.dfdx-xl.dfdx) / ( 2.*DEL );
  v.h[1] = (xh.dfdy-xl.dfdy + yh.dfdx-yl.dfdx) / ( 4.*DEL );
  v.h[2] = v.h[1];
  v.h[3] = (yh.dfdy-yl.dfdy) / ( 2.*DEL );
  ccdl::eigen( 2, v.h.data(), v.eval.data(), v.evec.data() );

  return v;
}






void ccdl::Mesh2d::FiniteDifferenceDebug( double x, double y ) const
{
  double DEL = 2.5e-5;
  ccdl::Mesh2dValue v = GetValue(x,y);
  double xhi = GetValue(x+DEL,y).f;
  double xlo = GetValue(x-DEL,y).f;
  double yhi = GetValue(x,y+DEL).f;
  double ylo = GetValue(x,y-DEL).f;
  double gx = (xhi-xlo)/(2*DEL);
  double gy = (yhi-ylo)/(2*DEL);

  std::printf("# debug @ %8.4f %8.4f : dfdx %12.4e %12.4e (%12.4e) dfdy %12.4e %12.4e (%12.4e)\n",
	      x, y, 
	      v.dfdx, gx, v.dfdx-gx,
	      v.dfdy, gy, v.dfdy-gy );

}


