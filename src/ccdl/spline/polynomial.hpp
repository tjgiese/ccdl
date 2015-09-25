#ifndef _ccdl_polynomial_hpp_
#define _ccdl_polynomial_hpp_

#include <vector>

namespace ccdl
{

  struct polynomial
  {
  public:

    polynomial() : n(0) {}
    polynomial( int N ) : n(N), c(N,0.) {}
    polynomial( int N, double const * C ) : n(N), c(C,C+N) {}


    // polynomial is \sum_i=0^N C[i]*pow(x,i)
    void reset( int N, double const * C );

    // clears all coefficients
    void reset( int N );

    int size() const { return n; }

    // compute value @ x
    double operator() ( double const x ) const { return value(x); }

    // compute nder'th derivative @ x
    double operator() ( double const x, int const nder ) const { return gradient( x, nder ); }

    // manually adjust polynomial coefficients
    double operator[] ( int i ) const { return c[i]; }
    double & operator[] ( int i ) { return c[i]; }

    // compute value @ x
    double value( double x ) const;

    // compute nder'th derivative @ x
    double gradient( double x, int const nder = 1 ) const;

    // return the first derivative of the polynomial
    ccdl::polynomial deriv( int const nder = 1 ) const;

    // get the roots of the polynomial
    // (the values of x where y=0)
    //
    // only the real roots are returned
    std::vector<double> roots() const;

    // the locations of the minima and maxima
    std::vector<double> extrema_loc() const;

    // the locations of the inflection points
    // up2down and down2up catagorize the inflections
    // based on whether they change the concavity
    // from concave-up to concave-down or vice-versa
    void inflection_loc( std::vector<double> & up2down, 
			 std::vector<double> & down2up ) const;

    // the locations of the minima
    std::vector<double> minima_loc() const;

    // the locations of the maxima
    std::vector<double> maxima_loc() const;

    // the location of the lowest minimum
    //
    // if the return value is false, then the polynomial
    // does not have a minimum 
    // (e.g., a line or a concave-down parabola)
    bool minimum_loc( double & x ) const;

    // the location of the lowest minimum
    //
    // if the return value is false, then the polynomial
    // does not have a minimum 
    // (e.g., a line or a concave-down parabola)
    bool maximum_loc( double & x ) const;

    // locations where two polynomials intersect
    std::vector<double> intersection_loc( ccdl::polynomial rhs ) const;

    // returns a line tangent to the curve @ x
    ccdl::polynomial tangent_line( double const x ) const;

    // returns a line normal to the curve @ x
    ccdl::polynomial normal_line( double const x ) const;

    // subtracts two polynomials (this-rhs)
    ccdl::polynomial subtract( ccdl::polynomial rhs ) const;

    // multiply two polynomials (this*rhs)
    ccdl::polynomial multiply( ccdl::polynomial rhs ) const;


    void fit( int npt, double const * x, double const * y );

    void fit( int npt, double const * x, double const * y,
	      int const icon, double const conval = 0. );

    void fit( int npt, double const * x, double const * y, 
	      double const * g );

    void fit( int npt, double const * x, double const * y, 
	      double const * g,
	      int const icon, double const conval = 0. );

    void wfit( int npt, double const * x, double const * y, 
	       double const * wy );

    void wfit( int npt, double const * x, double const * y, 
	       double const * wy,
	       int const icon, double const conval = 0. );

    void wfit( int npt, double const * x, double const * y, double const * g, 
	       double const * wy, double const * wg );

    void wfit( int npt, double const * x, double const * y, double const * g, 
	       double const * wy, double const * wg, 
	       int const icon, double const conval = 0. );


  private:
    int n;
    std::vector<double> c;
  };

}

#endif
