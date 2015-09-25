#include "polynomial.hpp"
#include "../bmath/lsfit.hpp"

#include <cstdio>
#include <cmath>
#include <algorithm>

extern "C"
{
  void
  dgeev_(const char* jobvl, const char* jobvr,
	 const int* n, double* a, const int* lda,
	 double* wr, double* wi, double* vl, const int* ldvl,
	 double* vr, const int* ldvr,
	 double* work, const int* lwork, int* info);

}



inline
void EraseSimilarValues( std::vector<double> & vec, double const TOL = 1.e-12 )
{
  if ( vec.size() > 1 )
    {
      std::sort( vec.begin(), vec.end() );
      std::vector<double>::iterator prev = vec.begin();
      for ( std::vector<double>::iterator 
	      p=vec.begin()+1; p != vec.end(); )
	{
	  if ( std::abs( *p - *prev ) < TOL )
	    p = vec.erase( p );
	  else
	    {
	      prev = p;
	      ++p;
	    };
	};
    };
}



void ccdl::polynomial::reset( int N, double const * C )
{
  n=N;
  c.assign( C, C+N );
}


void ccdl::polynomial::reset( int N )
{
  n=N;
  c.assign( N, 0. );
}


double ccdl::polynomial::value( double x ) const
{
  double f = 0.;
  double xn = 1.;
  for ( int i=0; i<n; ++i )
    {
      f += c[i] * xn;
      xn *= x;
    };
  return f;
}


double ccdl::polynomial::gradient( double x, int const nder ) const
{
  double f = 0.;
  double xn = 1.;
  if ( nder == 1 )
    {
      for ( int i=1; i<n; ++i )
	{
	  f += i * c[i] * xn;
	  xn *= x;
	};
    }
  else
    {
      std::vector<double> fact(n,1.);
      for ( int i=2; i<n; ++i )
	fact[i] = fact[i-1] * i;
      for ( int i=nder; i<n; ++i )
	{
	  f += (fact[i]/fact[i-nder]) * c[i] * xn;
	  xn *= x;
	};
    }
  return f;
}


ccdl::polynomial ccdl::polynomial::deriv( int const nder ) const
{
  ccdl::polynomial p;
  if ( n > nder )
    {
      std::vector<double> C( c );
      if ( nder == 1 )
	for ( int i=0; i<n; ++i )
	  C[i] *= i;
      else
	{
	  std::vector<double> fact( n, 1. );
	  for ( int i=2; i<n; ++i )
	    fact[i] = fact[i-1] * i;
	  for ( int i=nder; i<n; ++i )
	    C[i] *= (fact[i]/fact[i-nder]);
	}
      p.reset( n-nder, C.data() + nder );
    };
  return p;
}



std::vector<double> ccdl::polynomial::roots() const
{
  int m = n-1;
  // exclude zero-values high-order monomials
  // to avoid divide-by-zero in the companion matrix
  for ( int i=1; i<n; ++i )
    if ( std::abs( c[i] ) > 1.e-13 ) m = i;

  if ( m == 0 )
    return std::vector<double>();

  // this is the "companion matrix"
  // the eigenvalues of this matrix
  // are the roots of its "characteristic ccdl::polynomial",
  // which is, by construction, the ccdl::polynomial
  // of interest here
  std::vector<double> A( m*m, 0. );
  for ( int i=0; i<m-1; ++i )
    A[(i+1) + i*m] = 1.;
  for ( int i=0; i<m; ++i )
    A[i + (m-1)*m] = - c[i] / c[m];

  std::vector<double> xc(m,0.);
  std::vector<double> xs(m,0.);
  std::vector<double> work(1);
  int lwork = -1;
  int info = 0;
  double dummy;
  dgeev_( "N", "N",
	  &m, A.data(), &m,
	  xc.data(), xs.data(), 
	  &dummy, &m,
	  &dummy, &m,
	  work.data(), &lwork, &info);
  lwork = 1 + (int)work[0];
  work.resize( lwork );
  dgeev_( "N", "N",
	  &m, A.data(), &m,
	  xc.data(), xs.data(), 
	  &dummy, &m,
	  &dummy, &m,
	  work.data(), &lwork, &info);

  // discard the ccdl::polynomial's imaginary roots
  std::vector<double> x_real_roots(0);
  for ( int i=0; i<m; ++i )
    {
      if ( std::abs( xs[i] ) < 1.e-10 ) 
	x_real_roots.push_back( xc[i] );
    };

  EraseSimilarValues( x_real_roots );
  return x_real_roots;
}


std::vector<double> ccdl::polynomial::extrema_loc() const
{
  return deriv().roots();
}


void ccdl::polynomial::inflection_loc
( std::vector<double> & up2down,
  std::vector<double> & down2up ) const
{
  up2down.resize(0);
  down2up.resize(0);
  std::vector<double> pts( deriv(2).roots() );
  for ( std::vector<double>::iterator 
	  p=pts.begin(), pend=pts.end(); 
	p != pend; ++p )
    {
      double y = gradient( *p, 3 );
      if ( y > 0. )
	down2up.push_back( *p );
      else
	up2down.push_back( *p );
    }
}



std::vector<double> ccdl::polynomial::minima_loc() const
{
  std::vector<double> extrema( extrema_loc() );
  std::vector<double> minima(0);
  for ( std::vector<double>::iterator 
	  p=extrema.begin(), pend=extrema.end(); 
	p != pend; ++p )
    {
      double f2 = gradient( *p, 2 );
      if ( f2 > 0. )
	minima.push_back( *p );
    }
  return minima;
}



std::vector<double> ccdl::polynomial::maxima_loc() const
{
  std::vector<double> extrema( extrema_loc() );
  std::vector<double> maxima(0);
  for ( std::vector<double>::iterator 
	  p=extrema.begin(), pend=extrema.end(); 
	p != pend; ++p )
    {
      double f2 = gradient( *p, 2 );
      if ( f2 < 0. )
	maxima.push_back( *p );
    }
  return maxima;
}




bool ccdl::polynomial::minimum_loc( double & x ) const
{
  x = 0.;
  std::vector<double> minima( minima_loc() );
  //bool has_a_minimum = minima.size() > 0;
  double ymin = 1.e+100;
  for ( std::vector<double>::iterator 
	  p=minima.begin(), pend=minima.end(); 
	p != pend; ++p )
    {
      double y = value( *p );
      if ( y < ymin )
	{
	  ymin = y;
	  x = *p;
	};
    }
  return minima.size() > 0;
}



bool ccdl::polynomial::maximum_loc( double & x ) const
{
  x = 0.;
  std::vector<double> maxima( maxima_loc() );
  //bool has_a_minimum = maxima.size() > 0;
  double ymax = -1.e+100;
  for ( std::vector<double>::iterator 
	  p=maxima.begin(), pend=maxima.end(); 
	p != pend; ++p )
    {
      double y = value( *p );
      if ( y > ymax )
	{
	  ymax = y;
	  x = *p;
	};
    }
  return maxima.size() > 0;
}




// returns a line tangent to the curve @ x

ccdl::polynomial ccdl::polynomial::tangent_line( double const x ) const
{
  double y = value(x);
  double m = gradient(x);
  // y = m * x + b
  double b = y-m*x;
  ccdl::polynomial line(2);
  line[0] = b;
  line[1] = m;
  return line;
}

// returns a line normal to the curve @ x

ccdl::polynomial ccdl::polynomial::normal_line( double const x ) const
{
  double y = value(x);
  double m = gradient(x);
  if ( std::abs(m) < 1.e-40 )
    m = 1.e-40;
  m = - 1. / m;
  // y = m * x + b
  double b = y-m*x;
  ccdl::polynomial line(2);
  line[0] = b;
  line[1] = m;
  return line;
}



ccdl::polynomial ccdl::polynomial::subtract( ccdl::polynomial other ) const
{
  int nmax = std::max( n, other.n );
  ccdl::polynomial diff( nmax );
  for ( int i=0; i<n; ++i )
    diff[i] += c[i];
  for ( int i=0; i<other.n; ++i )
    diff[i] -= other[i];
  return diff;
}


ccdl::polynomial ccdl::polynomial::multiply( ccdl::polynomial other ) const
{
  int nmax = n + other.n - 1;
  ccdl::polynomial mult( nmax );
  for ( int i=0; i<n; ++i )
    for ( int j=0; j<other.n; ++j )
      mult[i+j] += c[i] * c[j];
  return mult;
}


std::vector<double> ccdl::polynomial::intersection_loc( ccdl::polynomial other ) const
{
  return subtract( other ).roots();
}





void ccdl::polynomial::fit( int npts, double const * x, double const * y ) 
{
  std::vector< double > A( npts * n, 0. );
  for ( int i=0; i<npts; ++i )
    for ( int k=0; k<n; ++k )
      A[i+k*n] = std::pow(x[i],k);
  ccdl::LeastSquaresFit( npts, n, A.data(), c.data(), y );
}

void ccdl::polynomial::fit( int npts, double const * x, double const * y, double const * g ) 
{
  std::vector< double > A( 2 * npts * n, 0. );
  for ( int i=0; i<npts; ++i )
    for ( int k=0; k<n; ++k )
      {
	A[i+(0+k*2)*npts] = std::pow(x[i],k);
	if ( k )
	  A[i+(1+k*2)*npts] = k * std::pow(x[i],k-1);
      };
  std::vector<double> obs( 2*npts, 0. );
  for ( int i=0; i<npts; ++i )
    {
      obs[i] = y[i];
      obs[i+npts] = g[i];
    }
  ccdl::LeastSquaresFit( 2*npts, n, A.data(), c.data(), obs.data() );
}



void ccdl::polynomial::wfit( int npts, double const * x, double const * y, 
			     double const * wy ) 
{
  std::vector< double > A( npts * n, 0. );
  for ( int i=0; i<npts; ++i )
    for ( int k=0; k<n; ++k )
      A[i+k*n] = std::pow(x[i],k);
  ccdl::WeightedLeastSquaresFit( npts, n, A.data(), c.data(), y, wy );
}

void ccdl::polynomial::wfit( int npts, double const * x, double const * y, double const * g,
			     double const * wy, double const * wg ) 
{
  std::vector< double > A( 2 * npts * n, 0. );
  for ( int i=0; i<npts; ++i )
    for ( int k=0; k<n; ++k )
      {
	A[i+(0+k*2)*npts] = std::pow(x[i],k);
	if ( k )
	  A[i+(1+k*2)*npts] = k * std::pow(x[i],k-1);
      };
  std::vector<double> obs( 2*npts, 0. ), w( 2*npts, 0. );
  for ( int i=0; i<npts; ++i )
    {
      obs[i] = y[i];
      obs[i+npts] = g[i];
      w[i] = wy[i];
      w[i+npts] = wg[i];
    }
  ccdl::WeightedLeastSquaresFit( 2*npts, n, A.data(), c.data(), obs.data(), w.data() );
}



void ccdl::polynomial::fit( int npts, double const * x, double const * y,
			    int const ck, double const conval ) 
{
  std::vector< double > A( npts * n, 0. );
  for ( int i=0; i<npts; ++i )
    for ( int k=0; k<n; ++k )
      A[i+k*n] = std::pow(x[i],k);

  std::vector<double> cvec( n, 0. ), cval( n, 0. );
  cvec[ck] = 1.;
  cval[ck] = conval;

  ccdl::ConstrainedLeastSquaresFit( npts, n, A.data(), c.data(), y, 
				    1, cvec.data(), cval.data() );
}

void ccdl::polynomial::fit( int npts, double const * x, double const * y, double const * g,
			    int const ck, double const conval ) 
{
  std::vector< double > A( 2 * npts * n, 0. );
  for ( int i=0; i<npts; ++i )
    for ( int k=0; k<n; ++k )
      {
	A[i+(0+k*2)*npts] = std::pow(x[i],k);
	if ( k )
	  A[i+(1+k*2)*npts] = k * std::pow(x[i],k-1);
      };
  std::vector<double> obs( 2*npts, 0. );
  for ( int i=0; i<npts; ++i )
    {
      obs[i] = y[i];
      obs[i+npts] = g[i];
    }

  std::vector<double> cvec( n, 0. ), cval( n, 0. );
  cvec[ck] = 1.;
  cval[ck] = conval;

  ccdl::ConstrainedLeastSquaresFit( 2*npts, n, A.data(), c.data(), obs.data(),
				    1, cvec.data(), cval.data() );
}


void ccdl::polynomial::wfit( int npts, double const * x, double const * y, 
			     double const * wy,
			     int const ck, double const conval ) 
{
  std::vector< double > A( npts * n, 0. );
  for ( int i=0; i<npts; ++i )
    for ( int k=0; k<n; ++k )
      A[i+k*n] = std::pow(x[i],k);

  std::vector<double> cvec( n, 0. ), cval( n, 0. );
  cvec[ck] = 1.;
  cval[ck] = conval;

  ccdl::ConstrainedWeightedLeastSquaresFit
    ( npts, n, A.data(), c.data(), y, wy, 
      1, cvec.data(), cval.data() );
}

void ccdl::polynomial::wfit( int npts, double const * x, double const * y, double const * g,
			     double const * wy, double const * wg,
			     int const ck, double const conval ) 
{
  std::vector< double > A( 2 * npts * n, 0. );
  for ( int i=0; i<npts; ++i )
    for ( int k=0; k<n; ++k )
      {
	A[i+(0+k*2)*npts] = std::pow(x[i],k);
	if ( k )
	  A[i+(1+k*2)*npts] = k * std::pow(x[i],k-1);
      };
  std::vector<double> obs( 2*npts, 0. ), w( 2*npts, 0. );
  for ( int i=0; i<npts; ++i )
    {
      obs[i] = y[i];
      obs[i+npts] = g[i];
      w[i] = wy[i];
      w[i+npts] = wg[i];
    }

  std::vector<double> cvec( n, 0. ), cval( n, 0. );
  cvec[ck] = 1.;
  cval[ck] = conval;

  ccdl::ConstrainedWeightedLeastSquaresFit
    ( 2*npts, n, A.data(), c.data(), obs.data(), w.data(), 
      1, cvec.data(), cval.data() );
}



