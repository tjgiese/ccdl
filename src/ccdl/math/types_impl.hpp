#ifndef _ccdl_types_impl_hpp_
#define _ccdl_types_impl_hpp_

#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>
#include <algorithm>
#include "types.hpp"
#include "exceptions.hpp"
#include "LAPACK.hpp"
#include "SPBLAS.hpp"

std::ostream & operator<< ( std::ostream & cout, ccdl::cv const a );
std::ostream & operator<< ( std::ostream & cout, ccdl::cge const a );
std::ostream & operator<< ( std::ostream & cout, ccdl::csy const a );
std::ostream & operator<< ( std::ostream & cout, ccdl::cdi const a );

namespace ccdl
{
  template <class M>
  void svd_driver( M const & ain, ccdl::ge & U, ccdl::di & w, ccdl::ge & VT, int const nscr, double * scr )
  {
    assert( U.nfast() == ain.nfast() );
    assert( U.nslow() == ain.nfast() );
    assert( w.nfast() == ain.nfast() );
    assert( w.nslow() == ain.nslow() );
    assert( VT.nfast() == ain.nslow() );
    assert( VT.nslow() == ain.nslow() );
    int const nf = ain.nfast();
    int const ns = ain.nslow();
    std::vector<int> iwork( 8*std::min(nf,ns) );
    int info = 0;
    ccdl::ge A(nf,ns,scr);
    int const N = nf*ns;
    for ( int i=0; i<N; ++i )
      A[i] = ain[i];
    int const lwork = nscr - N;
    FORTRAN_NAME(dgesdd)("A",&nf,&ns,A.data(),
                         &nf,w.data(),U.data(),
                         &nf,VT.data(),
                         &ns,A.end(),
                         &lwork,iwork.data(),&info);
    if ( info < 0 )
      {
	throw ccdl::illegal_argument("ccdl::svd_driver; dgesdd",-info);
      }
    else if ( info > 0 )
      {
	throw ccdl::exception("ccdl::svd_driver; dgesdd; DBDSDC did not converge, updating process failed.");
      };
  }


  template<class MA,class MB>
  ccdl::v & v_eq_ge_dot_v( ccdl::v & a, double const alpha, MA const & A, MB const & x, double const beta )
  {
    assert( a.size() == A.nfast() );
    assert( A.nslow() == x.size() );
    int const inc = 1;
    int const Af = A.nfast();
    int const As = A.nslow();

    //std::cout << Af << " " << As << " " << alpha << " " << A[0] << " " << Af << " " << x[0] << " " << inc << " " << beta << " " << a[0] << " " << inc << std::endl;
    FORTRAN_NAME(dgemv)("N",&Af,&As,&alpha,
			A.data(),&Af,
			x.data(),&inc,&beta,
			a.data(),&inc);
    return a;
  }

  template<class MA,class MB>
  ccdl::v & v_eq_gt_dot_v( ccdl::v & a, double const alpha, MA const & A, MB const & x, double const beta )
  {
    assert( a.size() == A.nslow() );
    assert( A.nfast() == x.size() );
    int const Af = A.nfast();
    int const As = A.nslow();
    int const inc = 1;
    FORTRAN_NAME(dgemv)("T",&Af,&As,&alpha,
			A.data(),&Af,
			x.data(),&inc,&beta,
			a.data(),&inc);
    return a;
  }

  template<class MA,class MB>
  ccdl::v & v_eq_sy_dot_v( ccdl::v & a, double const alpha, MA const & A, MB const & x, double const beta )
  {
    assert( a.size() == A.nfast() );
    assert( A.nslow() == x.size() );
    int const Af = A.nfast();
    int const inc = 1;
    FORTRAN_NAME(dsymv)("U",&Af,&alpha,
			A.data(),&Af,
			x.data(),&inc,&beta,
			a.data(),&inc);
    return a;
  }



  template<class MA,class MB>
  ccdl::ge & ge_eq_ge_dot_ge( ccdl::ge & a, double const alpha, MA const & A, MB const & B, double const beta )
  {
    assert( a.nfast() == A.nfast() );
    assert( a.nslow() == B.nslow() );
    assert( A.nslow() == B.nfast() );
    int const nf = a.nfast();
    int const ns = a.nslow();
    int const Af = A.nfast();
    int const As = A.nslow();
    int const Bf = B.nfast();
    FORTRAN_NAME(dgemm)("N","N",&nf,&ns,&As,
			&alpha,
			A.data(),&Af,
			B.data(),&Bf,
			&beta,
			a.data(),&nf);
    return a;
  }

  template<class MA,class MB>
  ccdl::ge & ge_eq_gt_dot_ge( ccdl::ge & a, double const alpha, MA const & A, MB const & B, double const beta )
  {
    assert( a.nfast() == A.nslow() );
    assert( a.nslow() == B.nslow() );
    assert( A.nfast() == B.nfast() );
    int const nf = a.nfast();
    int const ns = a.nslow();
    int const Af = A.nfast();
    int const Bf = B.nfast();
    FORTRAN_NAME(dgemm)("T","N",&nf,&ns,&Af,
			&alpha,
			A.data(),&Af,
			B.data(),&Bf,
			&beta,
			a.data(),&nf);
    return a;
  }

  template<class MA,class MB>
  ccdl::ge & ge_eq_ge_dot_gt( ccdl::ge & a, double const alpha, MA const & A, MB const & B, double const beta )
  {
    assert( a.nfast() == A.nfast() );
    assert( a.nslow() == B.nfast() );
    assert( A.nslow() == B.nslow() );
    int const nf = a.nfast();
    int const ns = a.nslow();
    int const As = A.nslow();
    int const Af = A.nfast();
    int const Bf = B.nfast();
    FORTRAN_NAME(dgemm)("N","T",&nf,&ns,&As,
			&alpha,
			A.data(),&Af,
			B.data(),&Bf,
			&beta,
			a.data(),&nf);
    return a;
  }

  template<class MA,class MB>
  ccdl::ge & ge_eq_gt_dot_gt( ccdl::ge & a, double const alpha, MA const & A, MB const & B, double const beta )
  {
    assert( a.nfast() == A.nslow() );
    assert( a.nslow() == B.nfast() );
    assert( A.nfast() == B.nslow() );
    int const nf = a.nfast();
    int const ns = a.nslow();
    int const Af = A.nfast();
    int const Bf = B.nfast();
    FORTRAN_NAME(dgemm)("T","T",&nf,&ns,&Af,
			&alpha,
			A.data(),&Af,
			B.data(),&Bf,
			&beta,
			a.data(),&nf);
    return a;
  }

  template<class MA,class MB>
  ccdl::ge & ge_eq_sy_dot_ge( ccdl::ge & a, double const alpha, MA const & A, MB const & B, double const beta )
  {
    assert( a.nfast() == A.nfast() );
    assert( a.nslow() == B.nslow() );
    assert( A.nslow() == B.nfast() );
    int const nf = a.nfast();
    int const ns = a.nslow();
    int const Af = A.nfast();
    int const Bf = B.nfast();
    FORTRAN_NAME(dsymm)("L","U",&nf,&ns,
			&alpha,
			A.data(),&Af,
			B.data(),&Bf,
			&beta,
			a.data(),&nf);
    return a;
  }

  template<class MA,class MB>
  ccdl::ge & ge_eq_ge_dot_sy( ccdl::ge & a, double const alpha, MA const & A, MB const & B, double const beta )
  {
    assert( a.nfast() == A.nfast() );
    assert( a.nslow() == B.nslow() );
    assert( A.nslow() == B.nfast() );
    int const nf = a.nfast();
    int const ns = a.nslow();
    int const Af = A.nfast();
    int const Bf = B.nfast();
    FORTRAN_NAME(dsymm)("R","U",&nf,&ns,
			&alpha,
			B.data(),&Bf,
			A.data(),&Af,
			&beta,
			a.data(),&nf);
    return a;
  }

  template<class MA,class MB>
  ccdl::ge & ge_eq_di_dot_ge( ccdl::ge & a, MA const & A, MB const & B )
  {
    assert( a.nfast() == A.nfast() );
    assert( a.nslow() == A.nslow() );
    assert( a.nslow() == B.nslow() );
    assert( A.nslow() == B.nfast() );
    int const ns = a.nslow();
    int const nf = a.nfast();
    int const Bf = B.nfast();
    for ( int j=0; j<ns; ++j )
      for ( int i=0; i<nf; ++i )
	a[i+j*nf] = A[i] * B[i+j*Bf];
    return a;
  }

  template<class MA,class MB>
  ccdl::ge & ge_eq_ge_dot_di( ccdl::ge & a, MA const & A, MB const & B )
  {
    assert( a.nfast() == A.nfast() );
    assert( a.nslow() == A.nslow() );
    assert( a.nslow() == B.nslow() );
    assert( A.nslow() == B.nfast() );
    int const ns = a.nslow();
    int const nf = a.nfast();
    int const Bf = B.nfast();
    for ( int j=0; j<ns; ++j )
      for ( int i=0; i<nf; ++i )
	a[i+j*nf] = A[i+j*Bf] * B[j];
    return a;
  }

  template<class MA,class MB>
  ccdl::ge & ge_eq_gt_dot_di( ccdl::ge & a, MA const & A, MB const & B )
  {
    assert( a.nfast() == A.nslow() );
    assert( a.nslow() == B.nslow() );
    assert( A.nfast() == B.nfast() );
    int const ns = a.nslow();
    int const nf = a.nfast();
    int const As = A.nslow();
    if ( ns < 1000 )
      {
	for ( int j=0; j<ns; ++j )
	  for ( int i=0; i<nf; ++i )
	    a[i+j*nf] = A[j+i*As] * B[j];
      }
    else
      {
	int const BLK = 8;
	for ( int jb=0; jb<ns; jb+=BLK )
	  {
	    int const ju = std::min(ns,jb+BLK);
	    for ( int ib=0; ib<nf; ib+=BLK )
	      {
		int const iu = std::min(nf,ib+BLK);
		for ( int j=jb; j<ju; ++j )
		  for ( int i=ib; i<iu; ++i )
		    a[i+j*nf] = A[j+i*As] * B[j];
	      };
	  };
      };
    return a;
  }


}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
inline ccdl::v::v( int const n, double * d )
  : ccdl::base::vector(n,d) 
{
}
inline ccdl::v::v( std::vector<double> & a )
  : ccdl::base::vector(a)
{
}
inline ccdl::v & ccdl::v::operator+= ( ccdl::cv const a )
{
  int const n = size();
  assert( n == a.size() );
  for ( int i=0; i<n; ++i )
    mData[i] += a[i];
  return *this;
}
inline ccdl::v & ccdl::v::operator-= ( ccdl::cv const a )
{
  int const n = size();
  assert( n == a.size() );
  for ( int i=0; i<n; ++i )
    mData[i] -= a[i];
  return *this;
}
inline ccdl::v & ccdl::v::operator= ( ccdl::cv const a )
{
  int const n = size();
  assert( n == a.size() );
  for ( int i=0; i<n; ++i )
    mData[i] = a[i];
  return *this;
}
inline ccdl::v & ccdl::v::operator= ( ccdl::v const & a )
{
  int const n = size();
  assert( n == a.size() );
  for ( int i=0; i<n; ++i )
    mData[i] = a[i];
  return *this;
}
inline ccdl::v & ccdl::v::operator= ( double const a )
{
  int const n = size();
  for ( int i=0; i<n; ++i )
    mData[i] = a;
  return *this;
}
inline ccdl::v & ccdl::v::operator*= ( ccdl::cv const a )
{
  int const n = size();
  assert( n == a.size() );
  for ( int i=0; i<n; ++i )
    mData[i] *= a[i];
  return *this;
}

inline double ccdl::v::dot( ccdl::cv const x ) const
{
  assert( size() == x.size() );
  double d = 0.;
  int const n = size();
  for ( int i=0; i<n; ++i )
    d += mData[i] * x[i];
  return d;
}
inline double ccdl::v::nrm2() const
{
  double d = 0.;
  for ( int i=0; i<size(); ++i)
    d += mData[i]*mData[i];
  return d;
}
inline int ccdl::v::iabsmax() const
{
  int imax=0;
  double mx = 0.;
  for ( int i=0; i<size(); ++i )
    {
      double t = std::abs(mData[i]);
      if ( t > mx )
        {
          mx = t;
          imax = i;
        };
    };
  return imax;
}
inline ccdl::v & ccdl::v::axpy( double const alpha, ccdl::cv const x )
{
  assert( size() == x.size() );
  int const n = size();
  for ( int i=0; i<n; ++i )
    mData[i] += alpha * x[i];
  return *this;
}
inline ccdl::v & ccdl::v::axpy( double const alpha, ccdl::cv const x, ccdl::cv const y )
{
  int const n = size();
  assert( n == x.size() );
  assert( n == y.size() );
  for ( int i=0; i<n; ++i )
    mData[i] = alpha * x[i] + y[i];
  return *this;
}

inline ccdl::v & ccdl::v::dot( double const alpha, ccdl::cge const A, ccdl::cv const  x, double const beta )
{
  return ccdl::v_eq_ge_dot_v(*this,alpha,A,x,beta);
}
inline ccdl::v & ccdl::v::dot( double const alpha, ccdl::ge const  A, ccdl::cv const  x, double const beta )
{
  return ccdl::v_eq_ge_dot_v(*this,alpha,A,x,beta);
}
inline ccdl::v & ccdl::v::dot( double const alpha, ccdl::cge const  A, ccdl::v const  x, double const beta )
{
  return ccdl::v_eq_ge_dot_v(*this,alpha,A,x,beta);
}
inline ccdl::v & ccdl::v::dot( double const alpha, ccdl::ge const  A, ccdl::v const  x, double const beta )
{
  return ccdl::v_eq_ge_dot_v(*this,alpha,A,x,beta);
}

inline ccdl::v & ccdl::v::dot( double const alpha, ccdl::cgt const  A, ccdl::cv const  x, double const beta )
{
  return ccdl::v_eq_gt_dot_v(*this,alpha,A,x,beta);
}
inline ccdl::v & ccdl::v::dot( double const alpha, ccdl::cgt const  A, ccdl::v const  x, double const beta )
{
  return ccdl::v_eq_gt_dot_v(*this,alpha,A,x,beta);
}

inline ccdl::v & ccdl::v::dot( double const alpha, ccdl::csy const  A, ccdl::cv const  x, double const beta )
{
  return ccdl::v_eq_sy_dot_v(*this,alpha,A,x,beta);
}
inline ccdl::v & ccdl::v::dot( double const alpha, ccdl::sy const  A, ccdl::cv const  x, double const beta )
{
  return ccdl::v_eq_sy_dot_v(*this,alpha,A,x,beta);
}
inline ccdl::v & ccdl::v::dot( double const alpha, ccdl::csy const  A, ccdl::v const  x, double const beta )
{
  return ccdl::v_eq_sy_dot_v(*this,alpha,A,x,beta);
}
inline ccdl::v & ccdl::v::dot( double const alpha, ccdl::sy const  A, ccdl::v const  x, double const beta )
{
  return ccdl::v_eq_sy_dot_v(*this,alpha,A,x,beta);
}

inline ccdl::v & ccdl::v::dot( ccdl::cge const  A, ccdl::cv const  x )
{
  return dot(1.,A,x,0.);
}
inline ccdl::v & ccdl::v::dot( ccdl::ge const  A, ccdl::cv const  x )
{
  return dot(1.,A,x,0.);
}
inline ccdl::v & ccdl::v::dot( ccdl::cge const  A, ccdl::v const  x )
{
  return dot(1.,A,x,0.);
}
inline ccdl::v & ccdl::v::dot( ccdl::ge const  A, ccdl::v const  x )
{
  return dot(1.,A,x,0.);
}

inline ccdl::v & ccdl::v::dot( ccdl::cgt const  A, ccdl::cv const  x )
{
  return dot(1.,A,x,0.);
}
inline ccdl::v & ccdl::v::dot( ccdl::cgt const  A, ccdl::v const  x )
{
  return dot(1.,A,x,0.);
}

inline ccdl::v & ccdl::v::dot( ccdl::csy const  A, ccdl::cv const  x )
{
  return dot(1.,A,x,0.);
}
inline ccdl::v & ccdl::v::dot( ccdl::sy const  A, ccdl::cv const  x )
{
  return dot(1.,A,x,0.);
}
inline ccdl::v & ccdl::v::dot( ccdl::csy const  A, ccdl::v const  x )
{
  return dot(1.,A,x,0.);
}
inline ccdl::v & ccdl::v::dot( ccdl::sy const  A, ccdl::v const  x )
{
  return dot(1.,A,x,0.);
}

inline int ccdl::v::query_solve( int const N ) 
{ 
  return N*N; 
}

inline int ccdl::v::query_svd_solve( int const M, int const N ) 
{ 
  return M*N + ccdl::ge::query_svd_inverse(M,N); 
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
inline ccdl::cv::cv( int const n, double const * d )
  : ccdl::base::const_vector(n,d)
{
}
inline ccdl::cv::cv( std::vector<double> const & a )
  : ccdl::base::const_vector(a)
{
}
inline ccdl::cv::cv( ccdl::v const & a )
  : ccdl::base::const_vector(a.size(),a.data())
{
}
inline double ccdl::cv::dot( ccdl::cv const x ) const
{
  assert( size() == x.size() );
  double d = 0.;
  int const n = size();
  for ( int i=0; i<n; ++i )
    d += mData[i] * x[i];
  return d;
}
inline double ccdl::cv::nrm2() const
{
  double d = 0.;
  for ( int i=0; i<size(); ++i)
    d += mData[i]*mData[i];
  return d;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
inline ccdl::di::di( int const nfast, int const nslow, double * data )
  : ccdl::base::vector(std::min(nfast,nslow),data), 
    mFast(nfast), mSlow(nslow)
{
}
inline ccdl::di::di( int const nfast, int const nslow, std::vector<double> & data )
  : ccdl::base::vector(std::min(nfast,nslow),data.data()), 
    mFast(nfast), mSlow(nslow)
{
}
inline int ccdl::di::nfast() const
{ 
  return mFast; 
}
inline int ccdl::di::nslow() const
{ 
  return mSlow; 
}
inline void ccdl::di::operator+= ( ccdl::cdi const a )
{
  int const n = size();
  assert( nfast() == a.nfast() );
  assert( nslow() == a.nslow() );
  for ( int i=0; i<n; ++i )
    mData[i] += a[i];
}
inline void ccdl::di::operator-= ( ccdl::cdi const a )
{
  int const n = size();
  assert( nfast() == a.nfast() );
  assert( nslow() == a.nslow() );
  for ( int i=0; i<n; ++i )
    mData[i] -= a[i];
}
inline ccdl::di & ccdl::di::operator= ( ccdl::cdi const a )
{
  int const n = size();
  assert( nfast() == a.nfast() );
  assert( nslow() == a.nslow() );
  for ( int i=0; i<n; ++i )
    mData[i] = a[i];
  return *this;
}
inline ccdl::di & ccdl::di::operator= ( ccdl::di const & a )
{
  int const n = size();
  assert( nfast() == a.nfast() );
  assert( nslow() == a.nslow() );
  for ( int i=0; i<n; ++i )
    mData[i] = a[i];
  return *this;
}
inline ccdl::di & ccdl::di::operator= ( double const a )
{
  int const n = size();
  for ( int i=0; i<n; ++i )
    mData[i] = a;
  return *this;
}
inline ccdl::cdi ccdl::di::t() const
{
  return ccdl::cdi(nslow(),nfast(),data());
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
inline ccdl::cdi::cdi( int const nfast, int const nslow, double const * data )
  : ccdl::base::const_vector(std::min(nfast,nslow),data),
    mFast(nfast), mSlow(nslow)
{
}
inline ccdl::cdi::cdi( int const nfast, int const nslow, std::vector<double> const & data )
  : ccdl::base::const_vector(std::min(nfast,nslow),data.data()),
    mFast(nfast), mSlow(nslow)
{
}
inline ccdl::cdi::cdi( ccdl::di const & a )
  : ccdl::base::const_vector(std::min(a.nfast(),a.nslow()),a.data()),
    mFast(a.nfast()), mSlow(a.nslow())
{
}
inline int ccdl::cdi::nfast() const
{
  return mFast;
}
inline int ccdl::cdi::nslow() const
{
  return mSlow;
}
inline ccdl::cdi ccdl::cdi::t() const
{
  return ccdl::cdi( nslow(),nfast(),data() );
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
inline ccdl::ge::ge( int m, int n, double * d )
  : ccdl::base::matrix(m,n,d) 
{
}
inline ccdl::ge::ge( int m, int n, std::vector<double> & a )
  : ccdl::base::matrix(m,n,a.data())
{
}
inline void ccdl::ge::operator+= ( ccdl::cge const a )
{
  assert( nfast() == a.nfast() );
  assert( nslow() == a.nslow() );
  int const n = size();
  for ( int i=0; i<n; ++i )
    mData[i] += a[i];
}
inline void ccdl::ge::operator-= ( ccdl::cge const a )
{
  assert( nfast() == a.nfast() );
  assert( nslow() == a.nslow() );
  int const n = size();
  for ( int i=0; i<n; ++i )
    mData[i] -= a[i];
}
inline ccdl::ge & ccdl::ge::operator= ( ccdl::cge const a )
{
  assert( nfast() == a.nfast() );
  assert( nslow() == a.nslow() );
  int const n = size();
  for ( int i=0; i<n; ++i )
    mData[i] = a[i];
  return *this;
}
inline ccdl::ge & ccdl::ge::operator= ( ccdl::ge const & a )
{
  assert( nfast() == a.nfast() );
  assert( nslow() == a.nslow() );
  int const n = size();
  for ( int i=0; i<n; ++i )
    mData[i] = a[i];
  return *this;
}
inline ccdl::ge & ccdl::ge::operator= ( double const a )
{
  int const n = size();
  for ( int i=0; i<n; ++i )
    mData[i] = a;
  return *this;
}
inline ccdl::cgt ccdl::ge::t() const
{
  return ccdl::cgt( nfast(), nslow(), data() );
}
inline ccdl::v ccdl::ge::col( int const i )
{
  return ccdl::v( nfast(), mData+i*nfast() );
}
inline ccdl::cv ccdl::ge::col( int const i ) const
{
  return ccdl::cv( nfast(), mData+i*nfast() );
}


inline ccdl::ge & ccdl::ge::dot( ccdl::sparsemat const & A, ccdl::cge const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::sparsesym const & A, ccdl::cge const B )
{
  return dot(1.,A,B,0.);
}

inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::cge const A, ccdl::cge const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_ge(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::ge const A, ccdl::cge const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_ge(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::cge const A, ccdl::ge const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_ge(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::ge const A, ccdl::ge const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_ge(*this,alpha,A,B,beta);
}

inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::cgt const A, ccdl::cge const B, double const beta )
{
  return ccdl::ge_eq_gt_dot_ge(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::cgt const A, ccdl::ge const B, double const beta )
{
  return ccdl::ge_eq_gt_dot_ge(*this,alpha,A,B,beta);
}

inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::cge const A, ccdl::cgt const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_gt(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::ge const A, ccdl::cgt const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_gt(*this,alpha,A,B,beta);
}

inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::cgt const A, ccdl::cgt const B, double const beta )
{
  return ccdl::ge_eq_gt_dot_gt(*this,alpha,A,B,beta);
}

inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::csy const A, ccdl::cge const B, double const beta )
{
  return ccdl::ge_eq_sy_dot_ge(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::sy const A, ccdl::cge const B, double const beta )
{
  return ccdl::ge_eq_sy_dot_ge(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::csy const A, ccdl::ge const B, double const beta )
{
  return ccdl::ge_eq_sy_dot_ge(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::sy const A, ccdl::ge const B, double const beta )
{
  return ccdl::ge_eq_sy_dot_ge(*this,alpha,A,B,beta);
}

inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::cge const A, ccdl::csy const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_sy(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::ge const A, ccdl::csy const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_sy(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::cge const A, ccdl::sy const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_sy(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::ge const A, ccdl::sy const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_sy(*this,alpha,A,B,beta);
}

inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::sy const A, ccdl::sy const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_sy(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::csy const A, ccdl::sy const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_sy(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::sy const A, ccdl::csy const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_sy(*this,alpha,A,B,beta);
}
inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::csy const A, ccdl::csy const B, double const beta )
{
  return ccdl::ge_eq_ge_dot_sy(*this,alpha,A,B,beta);
}


inline ccdl::ge & ccdl::ge::dot( ccdl::cge const A, ccdl::cge const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::ge const A, ccdl::cge const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::cge const A, ccdl::ge const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::ge const A, ccdl::ge const B )
{
  return dot(1.,A,B,0.);
}

inline ccdl::ge & ccdl::ge::dot( ccdl::cgt const A, ccdl::cge const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::cgt const A, ccdl::ge const B )
{
  return dot(1.,A,B,0.);
}

inline ccdl::ge & ccdl::ge::dot( ccdl::cge const A, ccdl::cgt const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::ge const A, ccdl::cgt const B )
{
  return dot(1.,A,B,0.);
}


inline ccdl::ge & ccdl::ge::dot( ccdl::cgt const A, ccdl::cgt const B )
{
  return dot(1.,A,B,0.);
}

inline ccdl::ge & ccdl::ge::dot( ccdl::csy const A, ccdl::cge const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::sy const A, ccdl::cge const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::csy const A, ccdl::ge const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::sy const A, ccdl::ge const B )
{
  return dot(1.,A,B,0.);
}

inline ccdl::ge & ccdl::ge::dot( ccdl::cge const A, ccdl::csy const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::ge const A, ccdl::csy const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::cge const A, ccdl::sy const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::ge const A, ccdl::sy const B )
{
  return dot(1.,A,B,0.);
}

inline ccdl::ge & ccdl::ge::dot( ccdl::sy const A, ccdl::sy const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::csy const A, ccdl::sy const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::sy const A, ccdl::csy const B )
{
  return dot(1.,A,B,0.);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::csy const A, ccdl::csy const B )
{
  return dot(1.,A,B,0.);
}

inline ccdl::ge & ccdl::ge::dot( ccdl::cdi const A, ccdl::cge const B )
{
  return ccdl::ge_eq_di_dot_ge(*this,A,B);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::di const A, ccdl::cge const B )
{
  return ccdl::ge_eq_di_dot_ge(*this,A,B);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::cdi const A, ccdl::ge const B )
{
  return ccdl::ge_eq_di_dot_ge(*this,A,B);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::di const A, ccdl::ge const B )
{
  return ccdl::ge_eq_di_dot_ge(*this,A,B);
}

inline ccdl::ge & ccdl::ge::dot( ccdl::cgt const A, ccdl::cdi const B )
{
  return ccdl::ge_eq_gt_dot_di(*this,A,B);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::cgt const A, ccdl::di const B )
{
  return ccdl::ge_eq_gt_dot_di(*this,A,B);
}



inline ccdl::ge & ccdl::ge::dot( ccdl::cge const A, ccdl::cdi const B )
{
  return ccdl::ge_eq_ge_dot_di(*this,A,B);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::cge const A, ccdl::di const B )
{
  return ccdl::ge_eq_ge_dot_di(*this,A,B);
}

inline ccdl::ge & ccdl::ge::dot( ccdl::ge const A, ccdl::cdi const B )
{
  return ccdl::ge_eq_ge_dot_di(*this,A,B);
}
inline ccdl::ge & ccdl::ge::dot( ccdl::ge const A, ccdl::di const B )
{
  return ccdl::ge_eq_ge_dot_di(*this,A,B);
}


inline int ccdl::ge::query_svd_inverse( int const M, int const N )
{
  return M*M + std::min(M,N) + N*N + ccdl::ge::query_svd(M,N);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
inline ccdl::cge::cge( int m, int n, double const * d )
  : ccdl::base::const_matrix(m,n,d)
{
}
inline ccdl::cge::cge( ccdl::ge const & a )
  : ccdl::base::const_matrix(a.nfast(),a.nslow(),a.data())
{
}
inline ccdl::cgt ccdl::cge::t() const
{
  return ccdl::cgt( nfast(), nslow(), data() );
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
inline ccdl::cgt::cgt( int m, int n, double const * d )
  : ccdl::base::const_matrix(m,n,d)
{
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
inline ccdl::sy::sy( int const n, double * d )
  : ccdl::base::matrix(n,n,d)
{
}
inline ccdl::sy::sy( int const n, std::vector<double> & d )
  : ccdl::base::matrix(n,n,d.data())
{
}
inline void ccdl::sy::operator+= ( ccdl::csy const a )
{
  assert( nfast() == a.nfast() );
  assert( nslow() == a.nslow() );
  int const n = size();
  for ( int i=0; i<n; ++i )
    mData[i] += a[i];
}
inline void ccdl::sy::operator-= ( ccdl::csy const a )
{
  assert( nfast() == a.nfast() );
  assert( nslow() == a.nslow() );
  int const n = size();
  for ( int i=0; i<n; ++i )
    mData[i] -= a[i];
}
inline ccdl::sy & ccdl::sy::operator= ( ccdl::csy const a )
{
  assert( nfast() == a.nfast() );
  assert( nslow() == a.nslow() );
  int const n = size();
  for ( int i=0; i<n; ++i )
    mData[i] = a[i];
  return *this;
}
inline ccdl::sy & ccdl::sy::operator= ( ccdl::sy const & a )
{
  assert( nfast() == a.nfast() );
  assert( nslow() == a.nslow() );
  int const n = size();
  for ( int i=0; i<n; ++i )
    mData[i] = a[i];
  return *this;
}
inline ccdl::sy & ccdl::sy::operator= ( double const a )
{
  int const n = size();
  for ( int i=0; i<n; ++i )
    mData[i] = a;
  return *this;
}
inline ccdl::cge ccdl::sy::ge() const
{
  return ccdl::cge(nfast(),nslow(),data());
}
inline ccdl::cge ccdl::csy::ge() const
{
  return ccdl::cge(nfast(),nslow(),data());
}
inline ccdl::ge ccdl::sy::ge()
{
  return ccdl::ge(nfast(),nslow(),data());
}
inline ccdl::sy & ccdl::sy::dot( double const alpha, ccdl::cge const A, ccdl::cge const B, double const beta )
{
  ge().dot(alpha,A,B,beta);
  return *this;
}
inline ccdl::sy & ccdl::sy::dot( ccdl::cge const A, ccdl::cge const B )
{
  ge().dot(1.,A,B,0.);
  return *this;
}
inline ccdl::sy & ccdl::sy::dot( double const alpha, ccdl::cgt const A, ccdl::cge const B, double const beta )
{
  ge().dot(alpha,A,B,beta);
  return *this;
}
inline ccdl::sy & ccdl::sy::dot( ccdl::cgt const A, ccdl::cge const B )
{
  ge().dot(1.,A,B,0.);
  return *this;
}

inline int ccdl::sy::query_svd_inverse( int const N )
{
  return ccdl::ge::query_svd_inverse(N,N);
}
inline ccdl::sy & ccdl::sy::svd_inverse( double const tol, int const nscr, double * scr )
{
  ge().svd_inverse(tol,nscr,scr);
  return *this;
}
inline ccdl::sy & ccdl::sy::svd_inverse( double const tol )
{
  std::vector<double> scr( ccdl::sy::query_svd_inverse(nfast()) );
  return svd_inverse(tol,scr.size(),scr.data());
}





inline int ccdl::sy::query_svd( int const N )
{
  return ccdl::ge::query_svd(N,N);
}
inline void ccdl::sy::svd( ccdl::ge & U, ccdl::di & w, ccdl::ge & VT, int const nscr, double * scr ) const
{
  ccdl::svd_driver(*this,U,w,VT,nscr,scr);
}
inline void ccdl::sy::svd( ccdl::ge & U, ccdl::di & w, ccdl::ge & VT ) const
{
  std::vector<double> scr(ccdl::sy::query_svd(nfast()));
  ccdl::svd_driver(*this,U,w,VT,scr.size(),scr.data());
}

inline void ccdl::sy::dsyev( ccdl::di & E, ccdl::ge & U ) const
{
  std::vector<double> scr( ccdl::sy::query_dsyev(nfast()) );
  dsyev(E,U,scr.size(),scr.data());
}
inline void ccdl::sy::dsyevd( ccdl::di & E, ccdl::ge & U ) const
{
  std::vector<double> scr( ccdl::sy::query_dsyevd(nfast()) );
  dsyevd(E,U,scr.size(),scr.data());
}
inline void ccdl::sy::dsyevr( ccdl::di & E, ccdl::ge & U ) const
{
  std::vector<double> scr( ccdl::sy::query_dsyevr(nfast()) );
  dsyevr(E,U,scr.size(),scr.data());
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
inline ccdl::csy::csy( int m, int n, double const * d )
  : ccdl::base::const_matrix(m,n,d)
{
}
inline ccdl::csy::csy( int m, int n, std::vector<double> const & d )
  : ccdl::base::const_matrix(m,n,d.data())
{
}
inline ccdl::csy::csy( ccdl::sy const & a )
  : ccdl::base::const_matrix(a.nfast(),a.nslow(),a.data())
{
}



inline void ccdl::csy::svd( ccdl::ge & U, ccdl::di & w, ccdl::ge & VT, int const nscr, double * scr ) const
{
  ccdl::svd_driver(*this,U,w,VT,nscr,scr);
}
inline void ccdl::csy::svd( ccdl::ge & U, ccdl::di & w, ccdl::ge & VT ) const
{
  std::vector<double> scr(ccdl::sy::query_svd(nfast()));
  ccdl::svd_driver(*this,U,w,VT,scr.size(),scr.data());
}

inline void ccdl::csy::dsyev( ccdl::di & E, ccdl::ge & U ) const
{
  std::vector<double> scr( ccdl::sy::query_dsyev(nfast()) );
  dsyev(E,U,scr.size(),scr.data());
}
inline void ccdl::csy::dsyevd( ccdl::di & E, ccdl::ge & U ) const
{
  std::vector<double> scr( ccdl::sy::query_dsyevd(nfast()) );
  dsyevd(E,U,scr.size(),scr.data());
}
inline void ccdl::csy::dsyevr( ccdl::di & E, ccdl::ge & U ) const
{
  std::vector<double> scr( ccdl::sy::query_dsyevr(nfast()) );
  dsyevr(E,U,scr.size(),scr.data());
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


inline ccdl::v & ccdl::v::dot( double const alpha, ccdl::sparsemat const & A, ccdl::cv const x, double const beta )
{
  if ( ! A.treat_as_transpose() )
    {
      sparse_gemv( A.nfast(), A.nslow(), alpha, 
		   A.csr_values(), A.csr_colidx(), A.csr_rowoff(),
		   beta, x.data(), data() );
    }
  else
    {
      sparse_gtmv( A.nfast(), A.nslow(), alpha, 
		   A.csr_values(), A.csr_colidx(), A.csr_rowoff(),
		   beta, x.data(), data() );
    };
  return *this;
}

inline ccdl::v & ccdl::v::dot( double const alpha, ccdl::sparsesym const & A, ccdl::cv const x, double const beta )
{
  sparse_symv( A.nfast(), alpha, 
	       A.csr_values(), A.csr_colidx(), A.csr_rowoff(),
	       beta, x.data(), data() );
  return *this;
}


inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::sparsemat const & A, ccdl::cge const B, double const beta )
{
  if ( ! A.treat_as_transpose() )
    {
      sparse_gemm( A.nfast(), A.nslow(), B.nfast(), nfast(), nslow(), 
		   alpha, A.csr_values(), A.csr_colidx(), A.csr_rowoff(),
		   beta, B.data(), data() );
    }
  else
    {
      sparse_gtmm( A.nfast(), A.nslow(), B.nfast(), nfast(), nslow(), 
		   alpha, A.csr_values(), A.csr_colidx(), A.csr_rowoff(),
		   beta, B.data(), data() );
    };
  return *this;
}


inline ccdl::ge & ccdl::ge::dot( double const alpha, ccdl::sparsesym const & A, ccdl::cge const B, double const beta )
{
  sparse_symm( A.nfast(), A.nslow(), B.nfast(), nfast(), nslow(), 
	       alpha, A.csr_values(), A.csr_colidx(), A.csr_rowoff(),
	       beta, B.data(), data() );
  return *this;
}


inline ccdl::v & ccdl::v::dot( ccdl::sparsemat const & A, ccdl::cv const x )
{
  return dot(1.,A,x,0.);
}

inline ccdl::v & ccdl::v::dot( ccdl::sparsesym const & A, ccdl::cv const x )
{
  return dot(1.,A,x,0.);
}



inline ccdl::v & ccdl::v::solve
( ccdl::sparsecholeskypreconditioner const & L, ccdl::cv const b )
{
  assert( L.nfast() == b.size() );
  assert( size() == b.size() );
  int const N = b.size();
  std::vector<double> t(N,0.);
  sparse_trisolve_Ax_eq_b(  N, 1., L.csr_values(), 
			    L.csr_colidx(), L.csr_rowoff(), 
			    b.data(), t.data() );
  sparse_trisolve_Atx_eq_b( N, 1., L.csr_values(), 
			    L.csr_colidx(), L.csr_rowoff(), 
			    t.data(), data() );
  return *this;
}



#endif

