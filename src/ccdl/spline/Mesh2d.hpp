#ifndef _ccdl_mesh2d_hpp_
#define _ccdl_mesh2d_hpp_

#include <vector>
#include <tr1/array>
#include <algorithm>
#include <iostream>
#include <cmath>

namespace ccdl
{

  struct Mesh2dValue
  {
    Mesh2dValue() : x(0),y(0),f(0),dfdx(0),dfdy(0) {}
    double x,y;
    double f;
    double dfdx;
    double dfdy;
  };

  struct Mesh2dHessian : public ccdl::Mesh2dValue
  {
    Mesh2dHessian() : Mesh2dValue() {}
    Mesh2dHessian( ccdl::Mesh2dValue v ) : Mesh2dValue(v) 
    {
      std::fill(h.begin(),h.end(),0.);
      std::fill(eval.begin(),eval.end(),0.);
      std::fill(evec.begin(),evec.end(),0.);
    }
    std::tr1::array<double,4> h;
    std::tr1::array<double,2> eval;
    std::tr1::array<double,4> evec;
  };

  bool sort_value_by_position( ccdl::Mesh2dValue const & a, ccdl::Mesh2dValue const & b );
  bool sort_hessian_by_position( ccdl::Mesh2dHessian const & a, ccdl::Mesh2dHessian const & b );


  class Mesh2d
  {
  public:


    // creates a mesh from
    //  x \in [0,lx) in steps of lx/nx     (periodic -- default)
    //  x \in [0,lx] in steps of lx/(nx-1) (nonperiodic)
    Mesh2d( int const nx, int const ny, 
	    double const lx, double const ly,
	    bool periodic = true );

    // creates a mesh from
    //  x \in [minx,maxx) in steps of lx/nx     (periodic -- default)
    //  x \in [minx,maxx] in steps of lx/(nx-1) (nonperiodic)
    Mesh2d( int const nx, int const ny, 
	    double const minx, double const maxx,
	    double const miny, double const maxy,
	    bool periodic = true );

    // Read a mesh in gnuplot format
    Mesh2d( std::istream & cin, 
	    bool periodic = true );

    // Default is "not walled"
    // By setting this to true, and if x is outside the
    // the bounds of the box, then GetValue(x,y) will reflect
    // the value of x such that it is equivalent to having
    // made an elastic collision with a hard-wall.
    // In other words, 
    //   If x > xhigh, GetValue(xhigh-(x-xhigh),y), else GetValue(x,y)
    // And similarly for y, or their combination, and for the lower
    // wall, etc...
    void SetWalled( bool logical );


    int GetSizeX() const;
    int GetSizeY() const;

    double GetX( int i ) const;
    double GetY( int i ) const;

    double GetLowX() const { return xlow; }
    double GetLowY() const { return ylow; }

    double GetHighX() const { return xlow+lx; }
    double GetHighY() const { return ylow+ly; }


    double & operator() ( int ix, int iy );
    double   operator() ( int ix, int iy ) const;

    void BsplineTransform( int order = 4 );

    double GetValueOnly( double x, double y ) const;
    ccdl::Mesh2dValue GetValue( double x, double y ) const;
    ccdl::Mesh2dHessian GetHessian( double x, double y ) const;

    void Add( double f );
    void Write( std::ostream & cout ) const;

    // std::vector< ccdl::Mesh2dValue > GetMinima( double const TOL=1.e-7 );
    // std::vector< ccdl::Mesh2dValue > GetMaxima( double const TOL=1.e-7 );
    // std::vector< ccdl::Mesh2dHessian > GetStationaryPts( double const TOL=1.e-7 );

    double GetLowestFreq( double x, double y ) const;

    void FiniteDifferenceDebug( double x, double y ) const;

    static int GetIndex( int i, int n );

  private:



    // By setting this to false, you are doing 2 things:
    // (1) Changing the definition of the "last point", i.e.
    //      x \in [0,lx) in steps of lx/nx     (periodic -- default)
    //  vs  x \in [0,lx] in steps of lx/(nx-1) (nonperiodic)
    // (2) A non-periodic cell contains a large, positive
    //     quadratic function outside of the box
    //
    // This would be fine as a public function if the constructors
    // only defined ranges, but if we are reading from a file
    // that already has the x,y values defined, then we can't change
    // how to discretize the range without ruining our data.
    //
    // Therefore, one must set the periodicity flag in the constructor
    // and we make this private.
    void SetPeriodic( bool logical );

  private:

    int nx,ny;

    // these are the box lengths
    double lx,ly;
    // for periodic systems, maxx and maxy are also the box lengths,
    // for for nonperiodic systems, it includes one extra delx
    // and dely layer so we can feed to to fftw3 and have it
    // interpret the points correctly
    double maxx,maxy;
    double delx,dely;

    // these define the lower corners of the box
    // the values of high and yhigh are redundant, as these are 
    // xhigh == xlow + lx
    // yhigh == ylow + ly
    double xlow,ylow; 

    // bspline order.  4-6 should be fine
    int order;

    bool walled,periodic;

    // raw, unmolested data
    std::vector<double> refdata;

    // data molested by bspline fourier muckery
    // this is the array to perform b-spline lookups on
    std::vector<double> data;

  };


}




inline int ccdl::Mesh2d::GetIndex( int i, int n )
{
  return ((i%n)+n)%n;
}

inline bool ccdl::sort_value_by_position( ccdl::Mesh2dValue const & a, ccdl::Mesh2dValue const & b )
{
  double dx = std::abs(a.x-b.x);
  bool xlt = dx < 1.e-4;
  return (!xlt) ? a.x < b.x : a.y < b.y;
}


inline bool ccdl::sort_hessian_by_position( ccdl::Mesh2dHessian const & a, ccdl::Mesh2dHessian const & b )
{
  double dx = std::abs(a.x-b.x);
  bool xlt = dx < 1.e-4;
  return (!xlt) ? a.x < b.x : a.y < b.y;
}



#endif
