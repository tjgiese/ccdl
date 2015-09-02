#include "bmath.hpp"

#include <vector>
#include <iostream>
#include <cmath>


#include "LAPACK.hpp"


double ccdl::v_dot_v( int n, double const * A, double const * B )
{
  int inc=1;
  return FORTRAN_NAME(ddot)( &n, A, &inc, B, &inc );
}


void ccdl::axpy( double const a, int const n, double * x, double * y )
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

  int info = 0;
  FORTRAN_NAME(dgesdd)("A",&nf,&ns,T.data(),
		       &nf,w,U,
		       &nf,VT,
		       &ns,work.data(),
		       &lwork,iwork.data(),&info);
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
  int info;
  FORTRAN_NAME(dpotrf)( "U", &nf, A, &nf, &info );

  if ( info == 0 )
    FORTRAN_NAME(dpotri)( "U", &nf, A, &nf, &info );
  if ( info == 0 )
    for ( int j=1; j < nf; ++j )
      for ( int i=0; i < j; ++i )
	A[j+i*nf] = A[i+j*nf];
  return info;
}
