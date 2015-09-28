#ifndef _OneDCurve_hpp_
#define _OneDCurve_hpp_

#include "../spline/polynomial.hpp"
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace ccdl
{

  namespace mini
  {

    struct OneDPoint
    {
      OneDPoint() : x(0),y(0),g(0),wy(0),wg(0) {}
      OneDPoint( double X, double Y ) : x(X),y(Y),g(0),wy(1.),wg(0.) {}
      OneDPoint( double X, double Y, double G ) : x(X),y(Y),g(G),wy(1.),wg(1.) {}
      OneDPoint( double X, double Y, int N, double const * S, double const * GRD );
      void SetBoltzmannWeight( double y0, double beta, bool with_g = true );
      double x,y,g,wy,wg;
      std::vector<double> grd;
    };


    struct OneDCurve
    {
      void push_back( double x, double y );
      void push_back( double x, double y, double g );
      void push_back( double X, double Y, int N, double const * S, double const * GRD );
      
      int size() const { return data.size(); }
      void clear() { data.resize(0); }
      void resize( int n ) { data.resize(n); };
      ccdl::mini::OneDPoint const & back() const { return data.back(); }
      ccdl::mini::OneDPoint & back() { return data.back(); }
      
      ccdl::mini::OneDCurve & sort_by_x();
      ccdl::mini::OneDCurve & sort_by_y();
      ccdl::mini::OneDCurve & sort_by_g();
      
      ccdl::mini::OneDPoint & operator[] ( int i ) { return data[i]; }
      ccdl::mini::OneDPoint const & operator[] ( int i ) const { return data[i]; }
      
      void SetBoltzmannWeight( double beta, bool with_g = true );
      
      ccdl::polynomial fit( int npoly, int icon = -1, double conval = 0. ) const;
      
      ccdl::mini::OneDCurve get_bracketing_curve() const;
      ccdl::mini::OneDCurve get_y_bracketing_curve() const;

      bool y_is_bracketed();
      bool g_is_bracketed();
      
      void print( std::ostream & cout ) const;
      
      std::vector< ccdl::mini::OneDPoint > data;
    };


  }

}


inline bool lt_x
( ccdl::mini::OneDPoint const & lhs, 
  ccdl::mini::OneDPoint const & rhs ) 
{ return lhs.x < rhs.x; }

inline bool lt_y
( ccdl::mini::OneDPoint const & lhs, 
  ccdl::mini::OneDPoint const & rhs ) 
{ return lhs.y < rhs.y; }

inline bool lt_g
( ccdl::mini::OneDPoint const & lhs, 
  ccdl::mini::OneDPoint const & rhs ) 
{ return lhs.g < rhs.g; }


inline
ccdl::mini::OneDCurve & ccdl::mini::OneDCurve::sort_by_x() 
{ std::sort( data.begin(), data.end(), lt_x ); return *this; }

inline
ccdl::mini::OneDCurve & ccdl::mini::OneDCurve::sort_by_y() 
{ std::sort( data.begin(), data.end(), lt_y ); return *this; }

inline
ccdl::mini::OneDCurve & ccdl::mini::OneDCurve::sort_by_g() 
{ std::sort( data.begin(), data.end(), lt_g ); return *this; }
      
#endif
