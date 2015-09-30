#include "bmath.hpp"

#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>

#include "LAPACK.hpp"


#undef PDBG


double ccdl::v_dot_v( int n, double const * A, double const * B )
{
  int inc=1;
  return FORTRAN_NAME(ddot)( &n, A, &inc, B, &inc );
}


void ccdl::v_tensor_v
( double const a, 
  int const m, double const * x, int const n, double const * y, 
  double * A )
{
  int inc=1;
  FORTRAN_NAME(dger)( &m, &n, &a, x, &inc, y, &inc, A, &m );
}



void ccdl::axpy( double const a, int const n, double const * x, double * y )
{
  int inc=1;
  FORTRAN_NAME(daxpy)( &n, &a, x, &inc, y, &inc );
}


void ccdl::ge_dot_v
( int nfA, int nsA, double const * A, 
  double const * x,
  double * Ax,
  double alpha, double beta )
{
  int const inc = 1;
  FORTRAN_NAME(dgemv)("N",&nfA,&nsA,&alpha,A,&nfA,x,&inc,&beta,Ax,&inc);
}



void ccdl::gt_dot_v( int nfA, int nsA, double const * A,
		      double const * x,
		      double * Ax,
		      double alpha,
		      double beta )
{
  int const inc = 1;
  FORTRAN_NAME(dgemv)("T",&nfA,&nsA,&alpha,
		      A,&nfA,
		      x,&inc,&beta,
		      Ax,&inc);
}



void ccdl::sy_dot_v( int nfA, int , double const * A,
		      double const * x,
		      double * Ax,
		      double alpha,
		      double beta  )
{
  int const inc = 1;
  FORTRAN_NAME(dsymv)("U",&nfA,&alpha,
		      A,&nfA,
		      x,&inc,&beta,
		      Ax,&inc);
}


void ccdl::ge_dot_ge
( int nfA, int nsA, double const * A, 
  int nfB, int nsB, double const * B, 
  double * AB,
  double alpha, double beta )
{
  int nf = nfA;
  int ns = nsB;
  FORTRAN_NAME(dgemm)("N","N",&nf,&ns,&nsA,&alpha, A,&nfA, B,&nfB,&beta, AB,&nf);
}


void ccdl::gt_dot_ge
( int nfA, int nsA, double const * A, 
  int nfB, int nsB, double const * B, 
  double * AB,
  double alpha, double beta )
{
  int nf = nsA;
  int ns = nsB;
  FORTRAN_NAME(dgemm)("T","N",&nf,&ns,&nfA,&alpha,A,&nfA,B,&nfB,&beta,AB,&nf);
}


void ccdl::ge_dot_gt
( int nfA, int nsA, double const * A, 
  int nfB, int , double const * B, 
  double * AB,
  double alpha, double beta )
{
  int nf = nfA;
  int ns = nfB;
  FORTRAN_NAME(dgemm)("N","T",&nf,&ns,&nsA,&alpha,A,&nfA,B,&nfB,&beta,AB,&nf);
}

//ok
void ccdl::di_dot_gt
( int nfA, int nsA, double const * A, 
  int nfB, int , double const * B,
  double * AB,
  double alpha, double beta )
{
  int nf = nfA;
  int ns = nfB;
  int nmin = std::min(nfA,nsA);
  for ( int i=0; i<nf*ns; ++i )
    AB[i] *= beta;
  for ( int j=0; j<ns; ++j ) 
    for ( int i=0; i<nmin; ++i )
      AB[i+j*nf] += alpha * A[i] * B[j+i*nfB];

  // int const inc = 1;
  // FORTRAN_NAME(dgemv)("N",&ns,&nmin,&alpha,
  // 		      B,&nfB,
  // 		      A,&inc,&beta,
  // 		      AB,&inc);
}
//ok
void ccdl::ge_dot_di
( int nfA, int , double const * A, 
  int nfB, int nsB, double const * B,
  double * AB,
  double alpha, double beta )
{
  int nf = nfA;
  int ns = nsB;
  int nmin = std::min(nfB,nsB);
  for ( int i=0; i<nf*ns; ++i )
    AB[i] *= beta;
  for ( int j=0; j<nmin; ++j ) 
    for ( int i=0; i<nf; ++i )
      AB[i+j*nf] += alpha * A[i+j*nf] * B[j];
}

//ok
void ccdl::di_dot_ge
( int nfA, int nsA, double const * A, 
  int nfB, int nsB, double const * B,
  double * AB,
  double alpha, double beta )
{
  int nf = nfA;
  int ns = nsB;
  int nmin = std::min(nfA,nsA);
  for ( int i=0; i<nf*ns; ++i )
    AB[i] *= beta;
  for ( int j=0; j<ns; ++j ) 
    for ( int i=0; i<nmin; ++i )
      AB[i+j*nf] += alpha * A[i] * B[i+j*nfB];
}




void ccdl::sy_dot_ge
( int nfA, int , double const * A, 
  int nfB, int nsB, double const * B, 
  double * AB,
  double alpha, double beta )
{
  int nf = nfA;
  int ns = nsB;
  FORTRAN_NAME(dsymm)("L","U",&nf,&ns,&alpha, A,&nfA, B,&nfB,&beta, AB,&nf);
}


void ccdl::ge_dot_sy
( int nfA, int , double const * A, 
  int nfB, int nsB, double const * B, 
  double * AB,
  double alpha, double beta )
{
  int nf = nfA;
  int ns = nsB;
  FORTRAN_NAME(dsymm)("R","U",&nf,&ns,&alpha, B,&nfB, A,&nfA,&beta, AB,&nf);
}




int ccdl::query_sdd
( int const M, int const N )
{
  double * p = NULL;
  double work = 0.;
  int lwork = -1, iwork = 0, info = 0;
  FORTRAN_NAME(dgesdd)("A",&M,&N,p,&M,p,p,&M,p,&N,&work,&lwork,&iwork,&info);
  if ( info < 0 )
    {
      std::cerr << "dgesdd : SVD failed, error code: " 
		<< info << std::endl;
    }
  else if ( info > 0 )
    {
      std::cerr << "dgesdd : DBDSDC did not converge, "
		<< "updating process failed" << std::endl;
    };
  return ((int)work) + M*N;
}


int ccdl::sdd_decomp
( int const nf, int const ns, double const * A,
  double * U, double * w, double * VT )
{
  std::vector<double> T(A,A+nf*ns);

  int lwork = query_sdd(nf,ns);
  std::vector<double> work( lwork, 0. );
  std::vector<int> iwork( 8*std::min(nf,ns) );

#ifdef PDBG
  std::printf("sdd_decomp\n");
#endif

  int info = 0;
  FORTRAN_NAME(dgesdd)("A",&nf,&ns,T.data(),
		       &nf,w,U,
		       &nf,VT,
		       &ns,work.data(),
		       &lwork,iwork.data(),&info);

#ifdef PDBG
  std::printf("return %i\n",info);
#endif

  return info;

  // if ( info < 0 )
  //   {
  //     std::cerr << "dgesdd : SVD failed; error code: " 
  // 		<< info << std::endl;
  //   }
  // else if ( info > 0 )
  //   {
  //     std::cerr << "dgesdd : DBDSDC did not converge, "
  // 		<< "updating process failed" << std::endl;
  //   };
}



int ccdl::sdd_power
( double const power, 
  int const nf, int const ns, double * A, 
  double TOL )
{
  int nmin = std::max(nf,ns);
  std::vector<double> U(nf*nf,0.),w(nmin,0.),VT(ns*ns,0.);

  int info = ccdl::sdd_decomp( nf,ns, A, U.data(), w.data(), VT.data() );

  if ( info == 0 )
    {
      for ( int i=0; i<nmin; ++i )
	{
	  if ( std::abs(w[i]) > TOL )
	    w[i] = std::pow(w[i],power);
	  else
	    w[i] = 0.;
	};
      
      std::vector<double> T(ns*nf,0.);
      //T = w^T . U^T
      ccdl::di_dot_gt( ns,nf, w.data(), nf,nf, U.data(), T.data() );
      //A^{-1} = VT^T . T
      ccdl::gt_dot_ge( ns,ns, VT.data(), ns,nf, T.data(), A );  
    };
  return info;
}


int ccdl::sdd_sym_power
( double const power, 
  int const nf, int const ns, double * A, 
  double TOL )
{
  int info = ccdl::sdd_power(power,nf,ns,A,TOL);
  if ( info == 0 )
    {
      for ( int i=1; i<nf; ++i )
	for ( int j=0; j<i; ++j )
	  {
	    double x = 0.5 * ( A[i+j*nf] + A[j+i*nf] );
	    A[i+j*nf] = x;
	    A[j+i*nf] = x;
	  };
    };
  return info;
}




int ccdl::sym_inverse( int const nf, int const , double * A )
{
#ifdef PDBG
  std::printf("sym_inverse\n");
#endif
  int info;
  FORTRAN_NAME(dpotrf)( "U", &nf, A, &nf, &info );

  if ( info == 0 )
    FORTRAN_NAME(dpotri)( "U", &nf, A, &nf, &info );

#ifdef PDBG
  std::printf("return %i\n",info);
#endif

  if ( info == 0 )
    for ( int j=1; j < nf; ++j )
      for ( int i=0; i < j; ++i )
	A[j+i*nf] = A[i+j*nf];
  return info;
}











int ccdl::query_dsyev( int const N )
{
  double p = 0.;
  int info = 0;
  int lwork=-1;
  double work = 0.;
  dsyev_("V","U", &N, &p,
	 &N, &p,
	 &work, &lwork, &info );
  //assert( info == 0 );
  //std::cout << "query_dsyev " << work << "\n";
  return (int)work;
}

int ccdl::query_dsyevd( int const N )
{
  int info = 0;
  int lwork=-1;
  double work = 0.;
  int liwork=-1;
  int iwork = 0;
  double p = 0.;
  dsyevd_("V","U",
	  &N,&p,
	  &N,&p,
	  &work, &lwork,
	  &iwork, &liwork,
	  &info );
  assert( info == 0 );
  return work;
}

int ccdl::query_dsyevr( int const N )
{
  int info  =  0, ISUPPZ = 0, IL = 1, IU = N, M = N;
  double work = 0., p = 0., VL = 0., VU = 0.;
  int lwork = -1, liwork=-1, iwork = 0;
  double ABSTOL = dlamch_("Safe minimum");
  dsyevr_("V","A","U",
	  &N,&p,&N,
	  &VL,&VU,
	  &IL,&IU,
	  &ABSTOL,
	  &M,&p,&p,&N,
	  &ISUPPZ,
	  &work,&lwork,
	  &iwork,&liwork,
	  &info );
  assert( info == 0 );
  return ((int)work) + N*N;
}

int ccdl::iquery_dsyevd( int const N )
{
  int info = 0;
  int lwork=-1;
  double work = 0.;
  int liwork=-1;
  int iwork = 0;
  double p = 0.;
  dsyevd_("V","U",
	  &N,&p,
	  &N,&p,
	  &work, &lwork,
	  &iwork, &liwork,
	  &info );
  assert( info == 0 );
  return iwork;
}

int ccdl::iquery_dsyevr( int const N )
{
  int info  =  0, ISUPPZ = 0, IL = 1, IU = N, M = N;
  double work = 0., p = 0., VL = 0., VU = 0.;
  int lwork = -1, liwork=-1, iwork = 0;
  double ABSTOL = 0.;
  dsyevr_("V","A","U",
	  &N,&p,&N,
	  &VL,&VU,
	  &IL,&IU,
	  &ABSTOL,
	  &M,&p,&p,&N,
	  &ISUPPZ,
	  &work,&lwork,
	  &iwork,&liwork,
	  &info );
  assert( info == 0 );
  return iwork+2*N;
}



int ccdl::dsyev
( int const N, double const * symat, 
  double * E, double * U, int const nscr, double * scr )
{
  int const nf = N;
  std::copy( symat, symat+nf*nf, U );
  int info = 0;

#ifdef PDBG
  std::printf("dsyev_\n");
#endif
  dsyev_("V","U", &nf, U,
	 &nf, E,
	 scr, &nscr, &info );

#ifdef PDBG
  std::printf("return %i\n",info);
#endif

  return info;
}

int ccdl::dsyevd
( int const N, double const * symat, 
  double * E, double * U, int const nscr, double * scr )
{
  int const nf = N;
  std::copy( symat, symat+nf*nf, U );
  int info = 0;
  int liwork= ccdl::iquery_dsyevd(nf);
  std::vector<int> iwork(liwork);

#ifdef PDBG
  std::printf("dsyevd_\n");
#endif
  dsyevd_("V","U",
	  &nf,U,
	  &nf,E,
	  scr, &nscr,
	  iwork.data(), &liwork,
	  &info );

#ifdef PDBG
  std::printf("return %i\n",info);
#endif
  return info;
}

int ccdl::dsyevr
( int const N, double const * symat, 
  double * E, double * U, int const nscr, double * scr )
{
  int const nf = N;
  int const N2 = nf*nf;
  double * T = scr;
  int const lwork = nscr-N2;
  double * work = scr + N2;
  for ( int i=0; i<N2; ++i )
    T[i] = symat[i];
  int info,M;
  int const liwork = ccdl::iquery_dsyevr(nf);
  std::vector<int> isuppz(liwork);
  int * iwork = isuppz.data() + 2*nf;
  const double VL = 0., VU = 0.;
  const int IL = 1, IU = nf;
  const double ABSTOL = dlamch_("Safe minimum");

#ifdef PDBG
  std::printf("dsyevr_\n");
#endif
  dsyevr_("V","A","U",
	  &nf,T,&nf,
	  &VL,&VU,
	  &IL,&IU,
	  &ABSTOL,
	  &M,E,U,&nf,
	  isuppz.data(),
	  work,&lwork,
	  iwork,&liwork,
	  &info );

#ifdef PDBG
  std::printf("return %i\n",info);
#endif
  return info;
}




#define EVMAX 70
#define EVRMAX 440


int ccdl::query_eigen( int const N )
{
  int n = 1;
  if ( N < EVMAX )
    { n = ccdl::query_dsyev(N); }
  else if ( N < EVRMAX )
    { n = ccdl::query_dsyevr(N); }
  else
    { n = ccdl::query_dsyevd(N); };
  return n;
}


int ccdl::eigen
( int const n, double const * symat, 
  double * E, double * U, int const nscr, double * scr )
{
  int info = 0;
  if ( n < EVMAX )
    { info = ccdl::dsyev(n,symat,E,U,nscr,scr); }
  else if ( n < EVRMAX )
    { info = ccdl::dsyevr(n,symat,E,U,nscr,scr); }
  else
    { info = ccdl::dsyevd(n,symat,E,U,nscr,scr); };
  return info;
}


int ccdl::eigen
( int const n, double const * symat, 
  double * E, double * U )
{
  int nscr = ccdl::query_eigen( n );
  std::vector<double> scr(nscr,0.);
  return eigen( n,symat,E,U,nscr,scr.data() );
}

#undef EVMAX
#undef EVRMAX


int ccdl::query_dsysv
( int const n )
{
  double a;
  int info = 0;
  int ipiv = 0;
  double work = 0.;
  int lwork = -1;
  int one = 1;
  FORTRAN_NAME(dsysv)( "U", &n, &one, 
		       &a, &n, 
		       &ipiv, 
		       &a, &n, 
		       &work, &lwork, &info );
  return ((int)work) + n*n;
}

int ccdl::dsysv
( int const n, double const * a, double const * b, double * x, 
  int const nscr, double * scr )
{
  std::fill(scr,scr+nscr,0.);
  std::copy( a, a+n*n, scr );
  std::copy( b, b+n, x );
  double * work = scr + n*n;
  int lwork = nscr - n*n;
  std::vector<int> ipiv(n,0);
  int info = 0;
  int one = 1;

#ifdef PDBG
  std::printf("ccdl::dsysv\n");
#endif
  FORTRAN_NAME(dsysv)( "U", &n, &one, 
		       scr, &n, 
		       ipiv.data(), 
		       x, &n, 
		       work, &lwork, &info );

#ifdef PDBG
  std::printf("return %i\n",info);
#endif
  return info;
}


int ccdl::dsysv
( int const n, double const * a, double const * b, double * x )
{
  std::vector<double> scr( ccdl::query_dsysv( n ), 0. );
  return dsysv( n,a,b,x, scr.size(), scr.data() );
}
