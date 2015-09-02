extern "C"
{
  double dlamch_(const char* cmach);

  double ddot_( const int *n,
	       double const *dx, const int *incx,
	       double const *dy, const int *incy );
  
  void daxpy_( const int *n, double const *alpha,
	       double const *dx, const int *incx,
	       double *dy, const int *incy );
  
  void dger_(const int *m, const int *n, double const *alpha,
	     double const *x, const int *incx,
	     double const *y, const int *incy,
	     double *a, const int *lda);

  void dgemv_(const char *trans, const int *m, const int *n,
	      double const *alpha, double const *a, const int *lda,
	      double const *x, const int *incx, double const *beta,
	      double *y, const int *incy);

  void dsymv_(const char *uplo, const int *n, double const *alpha,
	      double const *a, const int *lda,
	      double const *x, const int *incx,
	      double const *beta, double *y, const int *incy);

  void dposv_(const char* uplo, const int* n, const int* nrhs,
	      double* a, const int* lda,
	      double* b, const int* ldb, int* info);
  
  void dsysv_(const char* uplo, const int* n,
	      const int* nrhs,
	      double* a, const int* lda, int* ipiv,
	      double* b, const int* ldb,
	      double* work, const int* lwork, int* info);
  
  void dpotrf_(const char* uplo, const int* n,
	       double* a, const int* lda, int* info);
  
  void dpotri_(const char* uplo, const int* n,
	       double* a, const int* lda, int* info);
  
  void dsymv_(const char *uplo, const int *n, double const *alpha,
	      double const *a, const int *lda,
	      double const *x, const int *incx,
	      double const *beta, double *y, const int *incy);

  void dsyev_(const char* jobz, const char* uplo,
	      const int* n, double* a, const int* lda,
	      double* w, double* work, const int* lwork, int* info);

  void dsyevd_(const char* jobz, const char* uplo,
	       const int* n, double* a, const int* lda,
	       double* w, double* work, const int* lwork,
	       int* iwork, const int* liwork, int* info);

  void dsyevr_(const char *jobz, const char *range, const char *uplo,
	       const int *n, double *a, const int *lda,
	       double const *vl, double const *vu,
	       const int *il, const int *iu,
	       double const *abstol, int *m, double *w,
	       double *z, const int *ldz, int *isuppz,
	       double *work, const int *lwork,
	       int *iwork, const int *liwork,
	       int *info);

}

#define _USE_MATH_DEFINES
#include <iostream>
#include <cstdio>
#include <vector>
#include <cassert>
#include <cmath>

#include "OptPrimitives.hpp"


namespace lpck
{

  int query_dsyev( int const N )
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

  int query_dsyevd( int const N )
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

  int query_dsyevr( int const N )
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

  int iquery_dsyevd( int const N )
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

  int iquery_dsyevr( int const N )
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



  void dsyev( int const N, double const * symat, 
	      double * E, double * U, int const nscr, double * scr )
  {
    int const nf = N;
    std::copy( symat, symat+nf*nf, U );
    int info = 0;
    dsyev_("V","U", &nf, U,
	   &nf, E,
	   scr, &nscr, &info );
    assert( info == 0 );
  }

  void dsyevd( int const N, double const * symat, 
	       double * E, double * U, int const nscr, double * scr )
  {
    int const nf = N;
    std::copy( symat, symat+nf*nf, U );
    int info = 0;
    int liwork= lpck::iquery_dsyevd(nf);
    std::vector<int> iwork(liwork);
    dsyevd_("V","U",
	    &nf,U,
	    &nf,E,
	    scr, &nscr,
	    iwork.data(), &liwork,
	    &info );
    assert( info == 0 );
  }

  void dsyevr( int const N, double const * symat, 
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
    int const liwork = lpck::iquery_dsyevr(nf);
    std::vector<int> isuppz(liwork);
    int * iwork = isuppz.data() + 2*nf;
    const double VL = 0., VU = 0.;
    const int IL = 1, IU = nf;
    const double ABSTOL = dlamch_("Safe minimum");
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
    assert( info == 0 );
  }




#define EVMAX 70
#define EVRMAX 440


  int query_eigen( int const N )
  {
    int n = 1;
    if ( N < EVMAX )
      { n = lpck::query_dsyev(N); }
    else if ( N < EVRMAX )
      { n = lpck::query_dsyevr(N); }
    else
      { n = lpck::query_dsyevd(N); };
    return n;
  }


  void eigen( int const n, double const * symat, 
	      double * E, double * U, int const nscr, double * scr )
  {
    if ( n < EVMAX )
      { lpck::dsyev(n,symat,E,U,nscr,scr); }
    else if ( n < EVRMAX )
      { lpck::dsyevr(n,symat,E,U,nscr,scr); }
    else
      { lpck::dsyevd(n,symat,E,U,nscr,scr); };
  }


  void eigen( int const n, double const * symat, 
	      double * E, double * U )
  {
    int nscr = lpck::query_eigen( n );
    std::vector<double> scr(nscr,0.);
    eigen( n,symat,E,U,nscr,scr.data() );
  }

#undef EVMAX
#undef EVRMAX

}






double ccdl::gopt::CptDeltaX_TrustRadius
( int const n, double const * H, double const * g,
  double const trust_radius, double * dx )
{
  // computes Eq (22) via Eq (23)-(24)
  double rad2 = trust_radius * trust_radius;
  std::vector<double> E(n,0.),U(n*n,0.),T(n,0.);
  lpck::eigen(n,H,E.data(),U.data());



  double alpha=1.;
  double beta =0.;
  int one     =1;
  dgemv_( "T", &n, &n, &alpha, U.data(), &n, g, &one, &beta, T.data(), &one );
  // T[i]*T[i] == the numerator of Eq (24)
  
  double lambda_hi = E[0]-1.e-8;
  if ( lambda_hi > 0. ) lambda_hi = 0.;
  double tau2_hi = 0.;
  for ( int i=0; i<n; ++i ) tau2_hi += std::pow( T[i]/(E[i]-lambda_hi), 2 ); // Eq (24)

  double lambda = lambda_hi;
  double tau2 = tau2_hi;
  if ( tau2_hi > rad2 ) // step is larger than trust_radius
    {
      // get the lower bound of the bisection
      double lambda_lo = -5./trust_radius;
      double tau2_lo = 0.;
      for ( int i=0; i<n; ++i ) tau2_lo += std::pow( T[i]/(E[i]-lambda_lo), 2 ); // Eq (24)

      // initial midpoint of the bisection
      lambda = 0.5*(lambda_hi+lambda_lo);
      tau2 = 0.;
      for ( int i=0; i<n; ++i ) tau2 += std::pow( T[i]/(E[i]-lambda), 2 ); // Eq (24)

      double error = tau2-rad2;
      // bisect until close
      while ( std::abs(error) > rad2*1.e-10 )
	{
	  if ( tau2 > rad2 )
	    lambda_hi = lambda;
	  else if ( tau2 < rad2 )
	    lambda_lo = lambda;
	  lambda = 0.5*(lambda_hi+lambda_lo);
	  tau2 = 0.;
	  for ( int i=0; i<n; ++i ) tau2 += std::pow( T[i]/(E[i]-lambda), 2 ); // Eq (24)
	  error = tau2-rad2;
	};
    };

  // now that we have lambda, we can compute the step
  std::fill(dx,dx+n,0.);
  for ( int i=0; i<n; ++i )
    {
      if ( std::abs( E[i]-lambda ) > 1.e-30 )
	alpha = -T[i]/(E[i]-lambda);
      else
	alpha = 0.;
      //alpha = -T[i]/(E[i]-lambda);
      daxpy_( &n, &alpha, U.data()+i*n, &one, dx, &one ); // Eq (23)
    };
  return std::sqrt( tau2 ); // actual step size taken
}




double ccdl::gopt::CptDeltaX_EigenFollow
( int const n, double const * H, double const * g, double trust_radius, double * dx )
{
  // computes Eq (11) from Peng_IsraelJChem_1993_v33_p449
  double rad2 = trust_radius * trust_radius;
  std::vector<double> E(n,0.),U(n*n,0.),T(n,0.);
  lpck::eigen(n,H,E.data(),U.data());

  double alpha=1.;
  double beta =0.;
  int one     =1;
  dgemv_( "T", &n, &n, &alpha, U.data(), &n, g, &one, &beta, T.data(), &one );
  double vg = ddot_(&n,U.data(),&one,g,&one);
  double lambda = ( E[0] + std::sqrt( E[0]*E[0] + 4*vg*vg ) ) / 2.; 
  std::fill(dx,dx+n,0.);
  for ( int i=0; i<1; ++i )
    {
      if ( std::abs( E[i]-lambda ) > 1.e-30 )
	alpha = -T[i]/(E[i]-lambda);
      else
	alpha = 0.;
      daxpy_( &n, &alpha, U.data()+i*n, &one, dx, &one ); 
    };

  double tau2 = ddot_(&n,dx,&one,dx,&one);
  if ( tau2 < rad2 )
    {
      rad2 -= tau2;

      double lambda_hi = E[1]-1.e-8;
      if ( lambda_hi > 0. ) lambda_hi = 0.;
      double tau2_hi = 0.;
      for ( int i=1; i<n; ++i ) tau2_hi += std::pow( T[i]/(E[i]-lambda_hi), 2 ); // Eq (24)

      double lambda = lambda_hi;
      tau2 = tau2_hi;
      if ( tau2_hi > rad2 ) // step is larger than trust_radius
	{
	  // get the lower bound of the bisection
	  double lambda_lo = -5./trust_radius;
	  double tau2_lo = 0.;
	  for ( int i=1; i<n; ++i ) tau2_lo += std::pow( T[i]/(E[i]-lambda_lo), 2 ); // Eq (24)

	  // initial midpoint of the bisection
	  lambda = 0.5*(lambda_hi+lambda_lo);
	  tau2 = 0.;
	  for ( int i=1; i<n; ++i ) tau2 += std::pow( T[i]/(E[i]-lambda), 2 ); // Eq (24)

	  double error = tau2-rad2;
	  // bisect until close
	  while ( std::abs(error) > rad2*1.e-10 )
	    {
	      if ( tau2 > rad2 )
		lambda_hi = lambda;
	      else if ( tau2 < rad2 )
		lambda_lo = lambda;
	      lambda = 0.5*(lambda_hi+lambda_lo);
	      tau2 = 0.;
	      for ( int i=1; i<n; ++i ) tau2 += std::pow( T[i]/(E[i]-lambda), 2 ); // Eq (24)
	      error = tau2-rad2;
	    };
	};

      // now that we have lambda, we can compute the step
      //std::fill(dx,dx+n,0.);
      for ( int i=1; i<n; ++i )
	{
	  if ( std::abs( E[i]-lambda ) > 1.e-30 )
	    alpha = -T[i]/(E[i]-lambda);
	  else
	    alpha = 0.;
	  //alpha = -T[i]/(E[i]-lambda);
	  daxpy_( &n, &alpha, U.data()+i*n, &one, dx, &one ); // Eq (23)
	};
    };

  tau2 = ddot_(&n,dx,&one,dx,&one);
  //std::printf("%20.10e %20.10e\n",tau2,rad2);
  return std::sqrt( tau2 ); // actual step size taken
}




double ccdl::gopt::PredictEnergyChange( int const n, double const * H, double const * g, double const * dx )
{
  // computes Eq (13)
  std::vector<double> v( g, g+n );
  double alpha = 0.5;
  double beta  = 1.0;
  int one      = 1;
  dsymv_( "U", &n, &alpha, H, &n, g, &one, &beta, v.data(), &one );
  alpha = ddot_( &n, v.data(), &one, dx, &one );
  return alpha;
}

void ccdl::gopt::UpdateHessian_BFGS( int const n, double const * H, double const * g, 
			 double const * dx, double * Hnew )
{
  // computes Eq (18)
  int one = 1;
  std::copy(H,H+n*n,Hnew);
  double beta = 0, alpha = 1. / ddot_(&n,dx,&one,g,&one);
  dger_( &n, &n, &alpha, g, &one, g, &one, Hnew, &n ); // Hnew[i+j*n] += alpha*g[i]*g[j]
  alpha = 1.;
  std::vector<double> v(n,0.);
  dsymv_("U",&n,&alpha,H,&n,dx,&one,&beta,v.data(),&one );
  alpha = - 1. / ddot_(&n,dx,&one,v.data(),&one);
  dger_( &n, &n, &alpha, v.data(), &one, v.data(), &one, Hnew, &n ); // Hnew[i+j*n] += alpha*v[i]*v[j]
}


void ccdl::gopt::UpdateHessian_MS
( int const n, double const * H, double const * g, 
  double const * dx, double * Hnew )
{
  // computes Eq (20)
  int one = 1;
  std::copy(H,H+n*n,Hnew);
  std::vector<double> v(g,g+n);
  double alpha = -1.;
  double beta  =  1.;
  dsymv_("U",&n,&alpha,H,&n,dx,&one,&beta,v.data(),&one); // v=g-H.x
  alpha = ddot_(&n,dx,&one,v.data(),&one); // alpha = 1 / dx.v
  if ( std::abs(alpha) > 1.e-30 ) alpha = 1./alpha; else alpha = 0.;
  dger_( &n, &n, &alpha, v.data(), &one, v.data(), &one, Hnew, &n ); // Hnew[i+j*n] += alpha*v[i]*v[j]
}


void ccdl::gopt::UpdateHessian_Bofill_BFGS_MS
( int const n, double const * H, double const * g, 
  double const * dx, double * Hnew )
{
  // computes Eq (19)
  int one = 1;
  std::copy(H,H+n*n,Hnew);
  std::vector<double> v(g,g+n);
  double alpha = -1.;
  double beta  =  1.;
  dsymv_("U",&n,&alpha,H,&n,dx,&one,&beta,v.data(),&one); // v=g-H.x
  double xnrm = std::sqrt( ddot_(&n,dx,&one,dx,&one) );
  double vnrm = std::sqrt( ddot_(&n,v.data(),&one,v.data(),&one) );
  double nnrm = std::abs( ddot_(&n,dx,&one,v.data(),&one) );
  double phi = 1.;
  if ( xnrm*vnrm > 1.e-30 ) phi = nnrm / ( xnrm * vnrm ); // Eq (21)
  alpha = ddot_(&n,dx,&one,g,&one); // alpha = 1 / dx.v
  if ( std::abs( alpha ) > 1.e-30 ) alpha = 1./alpha; else alpha = 0.;
  alpha *= (1. - phi); // MS contribution
  dger_( &n, &n, &alpha, v.data(), &one, v.data(), &one, Hnew, &n ); // Hnew[i+j*n] += alpha*v[i]*v[j]
  // Hnew now contains the MS contribution to the update -- still need BFGS
  beta = 0;
  //alpha = 1. / ddot_(&n,dx,&one,g,&one);
  alpha = ddot_(&n,dx,&one,g,&one);
  if ( std::abs( alpha ) > 1.e-30 ) alpha = 1./alpha; else alpha = 0.;
  alpha *= phi; // BFGS contribution
  dger_( &n, &n, &alpha, g, &one, g, &one, Hnew, &n ); // Hnew[i+j*n] += alpha*g[i]*g[j]
  alpha = 1.;
  dsymv_("U",&n,&alpha,H,&n,dx,&one,&beta,v.data(),&one );
  alpha = ddot_(&n,dx,&one,v.data(),&one);
  if ( std::abs(alpha) > 1.e-30 ) alpha = -1./alpha; else alpha = 0.;
  alpha *= phi; // BFGS contribution
  dger_( &n, &n, &alpha, v.data(), &one, v.data(), &one, Hnew, &n ); // Hnew[i+j*n] += alpha*v[i]*v[j]
}

void ccdl::gopt::UpdateHessian_PSB
( int const n, double const * H, double const * g, 
  double const * dx, double * Hnew )
{
  // computes Eq (36)
  int one = 1;
  std::copy(H,H+n*n,Hnew);
  std::vector<double> v(g,g+n);
  double alpha = -1.;
  double beta  =  1.;
  dsymv_("U",&n,&alpha,H,&n,dx,&one,&beta,v.data(),&one); // v=g-H.x
  double pref1 = ddot_(&n,dx,&one,dx,&one);
  if ( pref1 > 1.e-30 ) pref1 = 1./pref1; else pref1 = 0.;
  double pref2 = -pref1*pref1 * ddot_(&n,dx,&one,v.data(),&one); // -x.v/(x.x)^2
  dger_( &n, &n, &pref1, dx, &one, v.data(), &one, Hnew, &n ); // Hnew[i+j*n] += pref2*dx[i]*v[j]
  dger_( &n, &n, &pref1, v.data(), &one, dx, &one, Hnew, &n ); // Hnew[i+j*n] += pref2*v[i]*dx[j]
  dger_( &n, &n, &pref2, dx, &one, dx, &one, Hnew, &n ); // Hnew[i+j*n] += pref2*dx[i]*dx[j]
}

void ccdl::gopt::UpdateHessian_Bofill_PSB_MS
( int const n, double const * H, double const * g, 
  double const * dx, double * Hnew )
{
  // computes Eq (37)
  int one = 1;
  std::copy(H,H+n*n,Hnew);
  std::vector<double> v(g,g+n);
  double alpha = -1.;
  double beta  =  1.;
  dsymv_("U",&n,&alpha,H,&n,dx,&one,&beta,v.data(),&one); // v=g-H.x
  double xnrm = ddot_(&n,dx,&one,dx,&one);
  double vnrm = ddot_(&n,v.data(),&one,v.data(),&one);
  double nnrm = ddot_(&n,dx,&one,v.data(),&one);
  double phi = 1.;
  if ( std::abs(xnrm*vnrm) > 1.e-30 ) phi = nnrm*nnrm / ( xnrm * vnrm ); // Eq (38)
  alpha = ddot_(&n,dx,&one,v.data(),&one); // alpha = 1 / dx.v
  if ( std::abs( alpha ) > 1.e-30 ) alpha = 1./alpha; else alpha = 0.;
  alpha *= phi; // MS contribution
  dger_( &n, &n, &alpha, v.data(), &one, v.data(), &one, Hnew, &n ); // Hnew[i+j*n] += alpha*v[i]*v[j]
  // Hnew now contains the MS contribution to the update -- still need PSB
  double pref1 = ddot_(&n,dx,&one,dx,&one);
  if ( pref1 > 1.e-30 ) pref1 = 1./pref1; else pref1 = 0.;
  double pref2 = -pref1*pref1 * ddot_(&n,dx,&one,v.data(),&one); // -x.v/(x.x)^2
  pref1 *= (1.-phi); // PSB contribution
  pref2 *= (1.-phi); // PSB contribution
  dger_( &n, &n, &pref1, dx, &one, v.data(), &one, Hnew, &n ); // Hnew[i+j*n] += pref2*dx[i]*v[j]
  dger_( &n, &n, &pref1, v.data(), &one, dx, &one, Hnew, &n ); // Hnew[i+j*n] += pref2*v[i]*dx[j]
  dger_( &n, &n, &pref2, dx, &one, dx, &one, Hnew, &n ); // Hnew[i+j*n] += pref2*dx[i]*dx[j]
}


// this is fastest
int ccdl::gopt::CptDeltaX_dsysv
( int const n, double const * H, double const * g, 
  double const lambda, double * dx )
{
  // computes Eq (22)
  // solve for x: 
  // (H-lambda*I).x = -g
  //
  // x = - (H-lambda*I)^{-1} . g
  // x = h^{-1} . (-g)
  // where h = (H-lambda*I)
  //
  int info = 0;
  int const one = 1;
  { // x = h^{-1} . (-g) using dsysv (symmetric positive indefinite)
    std::vector<double> h(H,H+n*n);
    for ( int i=0; i<n; ++i ) h[i+i*n] -= lambda;
    for ( int i=0; i<n; ++i ) dx[i] = -g[i];
    std::vector<int> ipiv(n,0);
    int lwork = -1;
    double w = 0.;
    dsysv_( "U", &n, &one, h.data(), &n, ipiv.data(), dx, &n, &w, &lwork, &info );
    if ( info != 0 ) { std::cerr << "dsysv: workspace error, info = " << info << "\n"; }
    lwork = (int)w;
    std::vector<double> work( lwork, 0. );
    dsysv_( "U", &n, &one, h.data(), &n, ipiv.data(), dx, &n, work.data(), &lwork, &info );
    if ( info != 0 ) { std::cerr << "dsysv: error, info = " << info << "\n"; }
    //for ( int i=0; i<n; ++i ) std::printf("%3i %20.10e\n",i,dx[i]);
  }

  return info;
}


// this is reasonably fast
int ccdl::gopt::CptDeltaX_dposv
( int const n, double const * H, double const * g, double const lambda, double * dx )
{
  // computes Eq (22)
  // solve for x: 
  // (H-lambda*I).x = -g
  //
  // x = - (H-lambda*I)^{-1} . g
  // x = h^{-1} . (-g)
  // where h = (H-lambda*I)
  //
  int info = 0;
  int const one = 1;
  { // x = h^{-1} . (-g) using dposv (symmetric positive definite)
    std::vector<double> h(H,H+n*n);
    for ( int i=0; i<n; ++i ) h[i+i*n] -= lambda;
    for ( int i=0; i<n; ++i ) dx[i] = -g[i];
    dposv_( "U", &n, &one, h.data(), &n, dx, &n, &info );
    if ( info != 0 ) { std::cerr << "dposv: error, info = " << info << "\n"; }
    //for ( int i=0; i<n; ++i ) std::printf("%3i %20.10e\n",i,dx[i]);
  }
  return info;
}




// this is super-slow
int ccdl::gopt::CptDeltaX_dpotrf
( int const n, double const * H, double const * g, double const lambda, double * dx )
{
  // computes Eq (22)
  // solve for x: 
  // (H-lambda*I).x = -g
  //
  // x = - (H-lambda*I)^{-1} . g
  // x = h^{-1} . (-g)
  // where h = (H-lambda*I)
  //
  int info = 0;
  {// x = h^{-1} . (-g) by inverting a symmetric positive definite matrix and matmul
    std::vector<double> h(H,H+n*n);
    for ( int i=0; i<n; ++i ) h[i+i*n] -= lambda;
    dpotrf_( "U", &n, h.data(), &n, &info );
    if ( info != 0 ) { std::cerr << "dpotrf info = " << info << "\n"; }
    dpotri_( "U", &n, h.data(), &n, &info );
    if ( info != 0 ) { std::cerr << "dpotri info = " << info << "\n"; }
    for ( int j=1; j < n; ++j )
      for ( int i=0; i < j; ++i )
	h[j+i*n] = h[i+j*n];
    // h is now h^{-1}
    double const alpha = -1.;
    double const beta  =  0.;
    int const    one   =  1;
    // x = -h^{-1}.g
    dsymv_( "U", &n, &alpha, h.data(), &n, g, &one, &beta, dx, &one );
    //for ( int i=0; i<n; ++i ) std::printf("%3i %20.10e\n",i,dx[i]);
  }

  return info;
}



// this is just debugging (it has a call to eigen in it
void ccdl::gopt::CptDeltaX_eigen
( int const n, double const * H, double const * g, double const lambda, double * dx )
{
  // computes Eq (23)
  std::vector<double> E(n,0.),U(n*n,0.),T(n,0.);
  lpck::eigen(n,H,E.data(),U.data());

  double alpha=1.;
  double beta =0.;
  int one     =1;
  dgemv_( "T", &n, &n, &alpha, U.data(), &n, g, &one, &beta, T.data(), &one );

  std::fill(dx,dx+n,0.);
  for ( int i=0; i<n; ++i )
    {
      if ( std::abs( E[i]-lambda ) > 1.e-30 )
	alpha = -T[i]/(E[i]-lambda);
      else
	alpha = 0.;
      daxpy_( &n, &alpha, U.data()+i*n, &one, dx, &one );
    };

}

// this is just debugging (it has a call to eigen in it
double ccdl::gopt::CptDeltaXNorm
( int const n, double const * H, double const * g, double const lambda  )
{
  // computes Eq (24)  well, the sqrt of Eq (24)
  std::vector<double> E(n,0.),U(n*n,0.),T(n,0.);
  lpck::eigen(n,H,E.data(),U.data());

  double alpha=1.;
  double beta =0.;
  int one     =1;
  dgemv_( "T", &n, &n, &alpha, U.data(), &n, g, &one, &beta, T.data(), &one );
  for ( int i=0; i<n; ++i ) T[i] *= T[i];

  double tau2 = 0.;
  for ( int i=0; i<n; ++i )
    tau2 += T[i]/std::pow( E[i]-lambda, 2 );

  return std::sqrt( tau2 );
}






// SIMPLE EXAMPLE


/*
struct stepinfo
{
  stepinfo( int n ) : n(n), x(n,0.), g(n,0.), h(n*n,0.), e(0.),predicted_de(0.) {};
  int n;
  std::vector<double> x,g,h; 
  double e,predicted_de;
};

template<class T>
void optimize( int const n, double * x, T fcn )
{
  int maxiter = 50;
  double trust_radius = 0.5;
  bool calcfc = true;
  stepinfo step(n), prevstep(n);
  std::vector<double> dg(n,0.);
  std::vector<double> dx(n,0.);
  double de = 0.;
  // input starting geom
  std::copy( x,x+n, step.x.data() );
  // obtain initial hessian
  {
    for ( int i=0; i<n; ++i ) step.h[i+i*n] = 1.;
    if ( calcfc )
      {
	double DEL = 1.e-4;
	for ( int i=0; i<n; ++i )
	  {
	    x[i] += DEL;
	    fcn(n,x,dg.data());
	    x[i] -= 2*DEL;
	    fcn(n,x,dx.data());
	    x[i] += DEL;
	    for ( int j=0; j<n; ++j )
	      step.h[j+i*n] = (dg[j]-dx[j])/(2.*DEL);
	  }
	for ( int i=0; i<n; ++i )
	  for ( int j=0; j<i; ++j )
	    {
	      double avg = 0.5*(step.h[i+j*n]+step.h[j+i*n]);
	      step.h[i+j*n] = avg;
	      step.h[j+i*n] = avg;
	    };
	std::fill(dx.data(),dx.data()+n,0.);
	std::fill(dg.data(),dg.data()+n,0.);
      };
  }

  for ( int iter = 0; iter < maxiter; ++iter )
    {
      // calc energy and gradient
      step.e = fcn(n,step.x.data(),step.g.data());
      // minimize between current and previous point (?)
      // todo ?
      de = step.e - prevstep.e;
      for ( int i=0; i<n; ++i ) dx[i] = step.x[i]-prevstep.x[i];
      for ( int i=0; i<n; ++i ) dg[i] = step.g[i]-prevstep.g[i];

      std::printf("e: %11.4e de: %9.1e pe %9.1e tr: %10.2e x:",
		  step.e,de,step.predicted_de,trust_radius);
      for ( int i=0; i<n; ++i ) std::printf("%12.4e",step.x[i]);
      std::printf(" g:");
      for ( int i=0; i<n; ++i ) std::printf("%10.2e",step.g[i]);
      std::printf(" h:");
      for ( int i=0; i<n*n; ++i ) std::printf("%9.2e",step.h[i]);
      std::printf("\n");

      double ade = std::abs(de);
      double dde = std::abs(de-step.predicted_de);
      if ( ade > 1.e-30 and iter )
	{
	  double ratio = dde/ade;
	  if ( ratio < 0.3 ) trust_radius *= 1.1;
	  else if ( ratio > 0.9 ) trust_radius /= 1.1;
	  //std::printf("dde %20.10e\n",dde/ade);
	};
      prevstep = step;
      // update hessian
      // MIN SEARCH
      //ccdl::gopt::UpdateHessian_BFGS( n, prevstep.h.data(), dg.data(), dx.data(), step.h.data() );
      //ccdl::gopt::UpdateHessian_MS( n, prevstep.h.data(), dg.data(), dx.data(), step.h.data() );
      //ccdl::gopt::UpdateHessian_Bofill_BFGS_MS( n, prevstep.h.data(), dg.data(), dx.data(), step.h.data() );
      //ccdl::gopt::UpdateHessian_PSB( n, prevstep.h.data(), dg.data(), dx.data(), step.h.data() );
      // TS SEARCH
      ccdl::gopt::UpdateHessian_Bofill_PSB_MS( n, prevstep.h.data(), dg.data(), dx.data(), step.h.data() );

      // MIN SEARCH
      //double steplen = ccdl::gopt::CptDeltaX_TrustRadius( n, prevstep.h.data(), prevstep.g.data(), trust_radius, dx.data() );

      // TS SEARCH
      double steplen = ccdl::gopt::CptDeltaX_EigenFollow( n, prevstep.h.data(), prevstep.g.data(), trust_radius, dx.data() );

      for ( int i=0; i<n; ++i ) step.x[i] = prevstep.x[i] + dx[i];
      std::copy( step.x.data(), step.x.data() + n, x );
      step.predicted_de = ccdl::gopt::PredictEnergyChange( n, prevstep.h.data(), prevstep.g.data(), dx.data() );
    }
}




double fcn( int const n, double const * x, double * g )
{
  double f = std::sin(x[0]*M_PI)*std::cos(x[1]*M_PI);
  g[0] = std::cos(x[0]*M_PI)*std::cos(x[1]*M_PI);
  g[1] = -std::sin(x[0]*M_PI)*std::sin(x[1]*M_PI);
  return f;
  // double e0 = std::exp( - x[0]*x[0] );
  // double e1 = std::exp( - (x[1]-1.)*(x[1]-1.) );
  // double e01 = std::exp( - (x[0]-0.)*(x[1]-1.) );
  // double f = e0 + 2.*e1 + e01;
  // g[0] = -2 * (x[0]-0.) * e0 - (x[1]-1.) * e01;
  // g[1] = -4 * (x[1]-1.) * e1 - (x[0]-0.) * e01;
  // return f;
}

int main()
{
  int n = 2;
  std::vector<double> x(n,0.5);
  x[1]=0.;
  optimize( n, x.data(), &fcn );
}

*/

