#include "lsfit.hpp"

#include <vector>
#include <cstdio>
#include <cmath>

#undef PDBG

extern "C"
{
  void dgelsy_(const int* m, const int* n, const int* nrhs,
	       double* a, const int* lda, double* b, const int* ldb,
	       int* jpvt, double const* rcond, int* rank,
	       double* work, const int* lwork, int* info);
  void dsysv_(const char* uplo, const int* n,
	      const int* nrhs,
	      double* a, const int* lda, int* ipiv,
	      double* b, const int* ldb,
	      double* work, const int* lwork, int* info);
  void dgglse_(const int* m, const int* n, const int* p,
	       double* a, const int* lda,
	       double* b, const int* ldb,
	       double* c, double* d, double* x,
	       double* work, const int* lwork, int* info);
  void dgemv_(const char *trans, const int *m, const int *n,
	      double const *alpha, double const *a, const int *lda,
	      double const *x, const int *incx, double const *beta,
	      double *y, const int *incy);
  void dgemm_(const char *transa, const char *transb, const int *m,
	      const int *n, const int *k, double const *alpha,
	      double const *a, const int *lda,
	      double const *b, const int *ldb,
	      double const *beta, double *c, const int *ldc);

  void dgelsd_(int const *m, int const *n, int const * nrhs,
	       double *a, int const * lda, double *b, int const * ldb, double * s,
	       double *rcond, int *rank, double *work, int *lwork,
	       int *iwork, int *info);

}



int ccdl::LeastSquaresFit
( int const nobs, int const nparam,
  double const * A_obs_by_param,
  double * x_param,
  double const * b_obs,
  double relative_accuracy_of_the_obs )
{
  // min (xt.At-bt).(A.x-b)

  std::vector<double> A( A_obs_by_param, A_obs_by_param + nparam*nobs );
  int nmax = std::max( nparam, nobs );
  std::vector<double> X( nmax, 0. );
  std::copy( b_obs, b_obs + nobs, X.data() );
  std::vector<int> jpvt( nparam, 0 );
  int LWORK = -1;
  int INFO = 0;
  double twork = 0;
  int rank = 0;
  int nrhs=1;
  dgelsy_( &nobs, &nparam, &nrhs, 
	   A.data(), &nobs, X.data(), &nmax, 
	   jpvt.data(), &relative_accuracy_of_the_obs, &rank,
	   &twork, &LWORK, &INFO );
  if ( INFO == 0 )
    {
      LWORK = twork+1;
      std::vector<double> WORK( LWORK, 0. );
#ifdef PDBG
      std::printf("dgelsy_\n");
#endif
      dgelsy_( &nobs, &nparam, &nrhs, 
	       A.data(), &nobs, X.data(), &nmax, 
	       jpvt.data(), &relative_accuracy_of_the_obs, &rank,
	       WORK.data(), &LWORK, &INFO );
#ifdef PDBG
      std::printf("return %i\n",INFO);
#endif
      std::copy( X.data(), X.data() + nparam, x_param );
    }
  else
    {
      std::fill( x_param, x_param + nparam, 0. );
    };
  return INFO;
}




int ccdl::LeastSquaresSvdFit
( int const nobs, int const nparam,
  double const * A_obs_by_param,
  double * x_param,
  double const * b_obs,
  double TOL )
{
  // min (xt.At-bt).(A.x-b)

  std::vector<double> A( A_obs_by_param, A_obs_by_param + nparam*nobs );
  int nmax = std::max( nparam, nobs );
  std::vector<double> X( nmax, 0. );
  std::copy( b_obs, b_obs + nobs, X.data() );
  std::vector<double> S( nmax, 0. );
  int LWORK = -1;
  std::vector<int> iwork(1,0);
  int INFO = 0;
  double twork = 0;
  int rank = 0;
  int nrhs=1;

  dgelsd_( &nobs, &nparam, &nrhs,
	   A.data(), &nobs, X.data(), &nmax,
	   S.data(), &TOL, &rank, &twork, &LWORK, iwork.data(), &INFO );

  
  if ( INFO == 0 )
    {
      LWORK = twork+1;
      std::vector<double> WORK( LWORK, 0. );
      int LIWORK = iwork[0];
      iwork.resize( LIWORK );
      
#ifdef PDBG
      std::printf("dgelsd_\n");
#endif
      
      dgelsd_( &nobs, &nparam, &nrhs,
	       A.data(), &nobs, X.data(), &nmax,
	       S.data(), &TOL, &rank, WORK.data(), &LWORK, iwork.data(), &INFO );

#ifdef PDBG
      std::printf("return %i\n",INFO);
#endif
      std::copy( X.data(), X.data() + nparam, x_param );
    }
  else
    {
      std::fill( x_param, x_param + nparam, 0. );
    };
  return INFO;
}



int ccdl::WeightedLeastSquaresFit
( int const nobs, int const nparam,
  double const * A_obs_by_param,
  double * x_param,
  double const * b_obs,
  double const * w_obs,
  double relative_accuracy_of_the_obs )
{
  std::vector<double> sqrtw( w_obs, w_obs + nobs );
  for ( int i=0; i<nobs; ++i )
    if ( w_obs[i] > 1.e-14 )
      sqrtw[i] = std::sqrt( w_obs[i] );
  
  std::vector<double> WA( A_obs_by_param, A_obs_by_param + nparam*nobs );
  for ( int j=0; j<nparam; ++j )
    for ( int i=0; i<nobs; ++i )
      WA[i+j*nobs] *= sqrtw[i];

  std::vector<double> WB( b_obs, b_obs + nobs );
  for ( int i=0; i<nobs; ++i )
    WB[i] *= sqrtw[i];
  
  return ccdl::LeastSquaresFit
    ( nobs, nparam, WA.data(), x_param, WB.data(),
      relative_accuracy_of_the_obs );
  
  /*
  // min (xt.At-bt).W.(A.x-b)
  // At.W.A.x = At.W.b
  // A'.x = b'
  // A' = At.W.A
  // b' = At.W.b
  std::vector<double> WA( A_obs_by_param, A_obs_by_param + nparam*nobs );
  for ( int j=0; j<nparam; ++j )
    for ( int i=0; i<nobs; ++i )
      WA[i+j*nobs] *= w_obs[i];

  std::vector<double> AtWtb( nparam, 0. );  
  double alpha = 1.;
  double beta  = 0.;
  int inc = 1;
  dgemv_( "T", &nobs, &nparam, &alpha,
	  WA.data(), &nobs, b_obs, &inc, &beta, AtWtb.data(), &inc );

  std::vector<double> AtWtA( nparam*nparam, 0. );
  dgemm_( "T","N", &nparam, &nparam, &nobs, &alpha,
	  WA.data(), &nobs, A_obs_by_param, &nobs, &beta, AtWtA.data(), &nparam );

  int INFO = 0;
  {
    std::vector<int> ipiv( nparam, 0 );
    int LWORK = -1;
    dsysv_("U",&nparam,&inc,AtWtA.data(),&nparam,ipiv.data(),AtWtb.data(),
	   &nparam,WA.data(),&LWORK,&INFO);
    LWORK = WA[0]+1;
    WA.resize( LWORK );
#ifdef PDBG
    std::printf("dsysv_\n");
#endif
    dsysv_("U",&nparam,&inc,AtWtA.data(),&nparam,ipiv.data(),AtWtb.data(),
	   &nparam,WA.data(),&LWORK,&INFO);
#ifdef PDBG
    std::printf("return %i\n",INFO);
#endif
    if ( INFO > 0 )
      throw ccdl::SingularMatrixException("ccdl::WeightedLeastSquaresFit");
    std::copy( AtWtb.data(), AtWtb.data()+nparam, x_param );
  }
  return INFO; 
  */
}




int ccdl::ConstrainedLeastSquaresFit
( int const nobs, int const nparam,
  double const * A_obs_by_param,
  double * x_param,
  double const * b_obs,
  int const ncon,
  double const * D_con_by_param,
  double const * c_con )
{
  // min (xt.At-bt).(A.x-b)   s.t. D.x = c

  std::vector<double> A( A_obs_by_param, A_obs_by_param + nparam*nobs );
  std::vector<double> b( b_obs, b_obs + nobs );
  std::vector<double> D( D_con_by_param, D_con_by_param + ncon*nparam );
  std::vector<double> c( c_con, c_con + ncon );

  int LWORK = -1;
  int INFO = 0;
  double twork = 0;

  dgglse_( &nobs, &nparam, &ncon,
	   A.data(), &nobs,
	   D.data(), &ncon,
	   b.data(), 
	   c.data(), 
	   x_param,
	   &twork, &LWORK, &INFO );

  if ( INFO == 0 )
    {
      LWORK = twork+1;
      std::vector<double> WORK( LWORK, 0. );
#ifdef PDBG
      std::printf("dgglse_\n");
#endif
      dgglse_( &nobs, &nparam, &ncon,
	       A.data(), &nobs,
	       D.data(), &ncon,
	       b.data(), 
	       c.data(), 
	       x_param,
	       WORK.data(), &LWORK, &INFO );
#ifdef PDBG
      std::printf("return %i\n",INFO);
#endif
      if ( INFO > 0 )
	throw ccdl::SingularMatrixException("ccdl::ConstrainedLeastSquaredFit");
    }
  else
    {
      std::fill( x_param, x_param + nparam, 0. );
    };
  return INFO;
}




int ccdl::ConstrainedWeightedLeastSquaresFit
( int const nobs, int const nparam,
  double const * A_obs_by_param,
  double * x_param,
  double const * b_obs,
  double const * w_obs,
  int const ncon,
  double const * D_con_by_param,
  double const * c_con )
{

  std::vector<double> sqrtw( w_obs, w_obs + nobs );
  for ( int i=0; i<nobs; ++i )
    if ( w_obs[i] > 1.e-14 )
      sqrtw[i] = std::sqrt( w_obs[i] );
  
  std::vector<double> WA( A_obs_by_param, A_obs_by_param + nparam*nobs );
  for ( int j=0; j<nparam; ++j )
    for ( int i=0; i<nobs; ++i )
      WA[i+j*nobs] *= sqrtw[i];

  std::vector<double> WB( b_obs, b_obs + nobs );
  for ( int i=0; i<nobs; ++i )
    WB[i] *= sqrtw[i];

  return ccdl::ConstrainedLeastSquaresFit
    ( nobs, nparam, WA.data(), x_param, WB.data(),
      ncon, D_con_by_param, c_con );
  
  /*
  // min (xt.At-bt).W.(A.x-b)   s.t. D.x = c

  std::vector<double> WA( A_obs_by_param, A_obs_by_param + nparam*nobs );
  for ( int j=0; j<nparam; ++j )
    for ( int i=0; i<nobs; ++i )
      WA[i+j*nobs] *= w_obs[i];

  std::vector<double> AtWtb( nparam, 0. );  
  double alpha = 1.;
  double beta  = 0.;
  int inc = 1;
  dgemv_( "T", &nobs, &nparam, &alpha,
	  WA.data(), &nobs, b_obs, &inc, &beta, AtWtb.data(), &inc );

  std::vector<double> AtWtA( nparam*nparam, 0. );
  dgemm_( "T","N", &nparam, &nparam, &nobs, &alpha,
	  WA.data(), &nobs, A_obs_by_param, &nobs, &beta, AtWtA.data(), &nparam );

  return ccdl::ConstrainedLeastSquaresFit
    ( nparam, nparam, AtWtA.data(), x_param, AtWtb.data(),
      ncon, D_con_by_param, c_con );
  */
}

