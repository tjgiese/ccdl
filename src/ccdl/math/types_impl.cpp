#include "types_impl.hpp"


ccdl::v & ccdl::v::solve( ccdl::csy const A, ccdl::cv const b, int const nscr, double * scr )
{
  if ( nscr < size()*size() )
    throw ccdl::exception("ccdl::v::solve scratch space too small");
  return dot( (  ccdl::sy(size(),scr)=A  ).inverse() , b );
}

ccdl::v & ccdl::v::solve( ccdl::csy const A, ccdl::cv const b )
{
  std::vector<double> scr(size()*size());  
  return dot( (  ccdl::sy(size(),scr)=A  ).inverse() , b );
}

ccdl::v & ccdl::v::constrained_solve( ccdl::csy const A, ccdl::cv const b, ccdl::cv const d, double const N )
{
  assert( d.nrm2() > 1.e-10 );
  std::vector<double> scr(size()*size()+size());
  ccdl::sy Ainv(size(),scr.data());
  ccdl::v t(size(),Ainv.end());
  t.dot( (Ainv=A).inverse() , d );
  double const mu = ( N - t.dot(b) ) / t.dot(d);
  return dot( Ainv, t.axpy(mu,d,b) );
}



ccdl::v & ccdl::v::svd_solve( ccdl::cge const A, ccdl::cv const b, double const tol, int const nscr, double * scr  )
{
  int const M = A.nfast();
  int const N = A.nslow();
  ccdl::ge Ainv( M,N, scr );
  ( Ainv = A ).svd_inverse( tol, nscr-M*N, Ainv.end() );
  return dot( Ainv, b );
}

ccdl::v & ccdl::v::svd_solve( ccdl::cge const A, ccdl::cv const b, double const tol )
{
  std::vector<double> scr( query_svd_solve( A.nfast(), A.nslow() ) );
  return svd_solve( A, b, tol, scr.size(), scr.data() );
}

ccdl::v & ccdl::v::svd_solve( ccdl::csy const A, ccdl::cv const b, double const tol, int const nscr, double * scr  )
{
  return svd_solve( A.ge(), b, tol, nscr, scr );
}

ccdl::v & ccdl::v::svd_solve( ccdl::csy const A, ccdl::cv const b, double const tol )
{
  return svd_solve( A.ge(), b, tol );
}


int ccdl::cv::iabsmax() const
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

ccdl::di & ccdl::di::inverse( double const tol )
{
  //int const nlow = std::min(4,size());
  int const nlow = size();
  for ( int i=0; i<nlow; ++i )
    {
      if ( std::abs(mData[i]) > tol )
        { mData[i] = 1./mData[i]; }
      else
        { mData[i] = 0.; };
    };
  //for ( int i=nlow; i<nmin; ++i ) data[i] = 1./data[i];
  std::swap(mFast,mSlow);
  return *this;
}


ccdl::ge & ccdl::ge::graham_schmidt()
{
  double const TOL = 1./1.e-10;
  for ( int i=0; i<nslow(); ++i )
    {
      ccdl::v vi( col(i) );
      double inorm = 1./std::sqrt( vi.nrm2() );
      vi *= inorm;
      if ( inorm > TOL )
        {
          vi = 0.;
          continue;
        };
      for ( int j=0; j<i; ++j )
        {
          ccdl::v vj( col(j) );
          double f = vi.dot(vj);
          vi.axpy( -f, vj );
        };
      double onorm = 1./std::sqrt( vi.nrm2() );
      vi *= onorm;
      if ( onorm > TOL )
        vi = 0.;
    };
  return *this;
}

ccdl::ge & ccdl::ge::svd_inverse( double const tol, int const nscr, double * scr )
{
  int const nf=nfast();
  int const ns=nslow();
  ccdl::ge U(nf,nf,scr);
  ccdl::di w(nf,ns,U.end());
  ccdl::ge VT(ns,ns,w.end());
  ccdl::ge t(nf,ns,VT.end());
  assert( nscr > std::distance(scr,VT.end()) );
  svd(U,w,VT,nscr-std::distance(scr,VT.end()),VT.end());
  std::swap(mFast,mSlow);
  return dot( t.dot( w.inverse(tol).t(), VT ).t(), U.t() );
}

ccdl::ge & ccdl::ge::svd_inverse( double const tol )
{
  std::vector<double> scr( ccdl::ge::query_svd_inverse(nfast(),nslow()) );
  return svd_inverse(tol,scr.size(),scr.data());
}

int ccdl::ge::query_svd( int const M, int const N )
{
  double * p = NULL;
  double work = 0.;
  int lwork = -1, iwork = 0, info = 0;
  FORTRAN_NAME(dgesdd)("A",&M,&N,p,&M,p,p,&M,p,&N,&work,&lwork,&iwork,&info);
  if ( info < 0 )
    {
      throw ccdl::illegal_argument("query_svd_inverse; dgesdd",-info);
    }
  else if ( info > 0 )
    {
      throw ccdl::exception("query_svd_inverse; dgesdd; DBDSDC did not converge, updating process failed.");
    };
  return ((int)work) + M*N;
}

void ccdl::ge::svd( ccdl::ge & U, ccdl::di & w, ccdl::ge & VT, int const nscr, double * scr ) const
{
  ccdl::svd_driver(*this,U,w,VT,nscr,scr);
}

void ccdl::ge::svd( ccdl::ge & U, ccdl::di & w, ccdl::ge & VT ) const
{
  std::vector<double> scr(ccdl::ge::query_svd(nfast(),nslow()));
  ccdl::svd_driver(*this,U,w,VT,scr.size(),scr.data());
}

ccdl::sy & ccdl::sy::inverse()
{
  int const nf = nfast();
  int info;
  FORTRAN_NAME(dpotrf)( "U", &nf, mData, &nf, &info );
  if ( info != 0 )
    {
      if ( info < 0 )
        throw ccdl::illegal_argument("ccdl::sy::inverse dpotrf",-info);
      if ( info > 0 )
        throw ccdl::exception("ccdl::sy::inverse dpotrf says matrix is not positive definite");
    };
  FORTRAN_NAME(dpotri)( "U", &nf, mData, &nf, &info );
  if ( info != 0 )
    {
      if ( info < 0 )
        throw ccdl::illegal_argument("ccdl::sy::inverse dpotri",-info);
      if ( info > 0 )
        throw ccdl::exception("ccdl::sy::inverse dpotri says the inverse could not be computed");
    };
  for ( int j=1; j < nf; ++j )
    for ( int i=0; i < j; ++i )
      mData[j+i*nf] = mData[i+j*nf];
  return *this;
}

int ccdl::sy::query_dsyev( int const N )
{
  double p = 0.;
  int info = 0;
  int lwork=-1;
  double work = 0.;
  FORTRAN_NAME(dsyev)("V","U", &N, &p,
                      &N, &p,
                      &work, &lwork, &info );
  assert( info == 0 );
  //std::cout << "query_dsyev " << work << "\n";
  return (int)work;
}

int ccdl::sy::query_dsyevd( int const N )
{
  int info = 0;
  int lwork=-1;
  double work = 0.;
  int liwork=-1;
  int iwork = 0;
  double p = 0.;
  FORTRAN_NAME(dsyevd)("V","U",
                       &N,&p,
                       &N,&p,
                       &work, &lwork,
                       &iwork, &liwork,
                       &info );
  assert( info == 0 );
  return work;
}

int ccdl::sy::query_dsyevr( int const N )
{
  int info  =  0, ISUPPZ = 0, IL = 1, IU = N, M = N;
  double work = 0., p = 0., VL = 0., VU = 0.;
  int lwork = -1, liwork=-1, iwork = 0;
  double ABSTOL = FORTRAN_NAME(dlamch)("Safe minimum");
  FORTRAN_NAME(dsyevr)("V","A","U",
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

int ccdl::sy::iquery_dsyevd( int const N )
{
  int info = 0;
  int lwork=-1;
  double work = 0.;
  int liwork=-1;
  int iwork = 0;
  double p = 0.;
  FORTRAN_NAME(dsyevd)("V","U",
                       &N,&p,
                       &N,&p,
                       &work, &lwork,
                       &iwork, &liwork,
                       &info );
  assert( info == 0 );
  return iwork;
}

int ccdl::sy::iquery_dsyevr( int const N )
{
  int info  =  0, ISUPPZ = 0, IL = 1, IU = N, M = N;
  double work = 0., p = 0., VL = 0., VU = 0.;
  int lwork = -1, liwork=-1, iwork = 0;
  double ABSTOL = 0.;
  FORTRAN_NAME(dsyevr)("V","A","U",
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

void ccdl::sy::dsyev( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const
{
  int const nf = nfast();
  for ( int i=0; i<nf*nf; ++i )
    U[i] = mData[i];
  int info = 0;
  FORTRAN_NAME(dsyev)("V","U", &nf, U.data(),
                      &nf, E.data(),
                      scr, &nscr, &info );
  assert( info == 0 );
}

void ccdl::sy::dsyevd( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const
{
  int const nf = nfast();
  for ( int i=0; i<nf*nf; ++i )
    U[i] = mData[i];
  int info = 0;
  int liwork= ccdl::sy::iquery_dsyevd(nf);
  std::vector<int> iwork(liwork);
  FORTRAN_NAME(dsyevd)("V","U",
                       &nf,U.data(),
                       &nf,E.data(),
                       scr, &nscr,
                       iwork.data(), &liwork,
                       &info );
  assert( info == 0 );
}

void ccdl::sy::dsyevr( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const
{
  int const nf = nfast();
  int const N2 = nf*nf;
  double * T = scr;
  int const lwork = nscr-N2;
  double * work = scr + N2;
  for ( int i=0; i<N2; ++i )
    T[i] = mData[i];
  int info,M;
  int const liwork = ccdl::sy::iquery_dsyevr(nf);
  std::vector<int> isuppz(liwork);
  int * iwork = isuppz.data() + 2*nf;
  const double VL = 0., VU = 0.;
  const int IL = 1, IU = nf;
  const double ABSTOL = FORTRAN_NAME(dlamch)("Safe minimum");
  FORTRAN_NAME(dsyevr)("V","A","U",
                       &nf,T,&nf,
                       &VL,&VU,
                       &IL,&IU,
                       &ABSTOL,
                       &M,E.data(),U.data(),&nf,
                       isuppz.data(),
                       work,&lwork,
                       iwork,&liwork,
                       &info );
  assert( info == 0 );
}




void ccdl::csy::dsyev( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const
{
  int const nf = nfast();
  for ( int i=0; i<nf*nf; ++i )
    U[i] = mData[i];
  int info = 0;

  // std::cout << "ccdl::csy::dsyev " << nf << " " << U.nfast() << " " << U.nslow() << " " << E.nfast() << " " << nscr << " " << std::endl;
  // for ( int i=0; i<nf*nf; ++i )
  //   std::cout << "U[" << i << "] = " << U[i] << "\n";
  // for ( int i=0; i<nf; ++i )
  //   std::cout << "E[" << i << "] = " << E[i] << "\n";
  // for ( int i=0; i<nscr; ++i )
  //   std::cout << "scr["<<i<<"] = " << scr[i] << "\n";
  //int NSCR=10000;
  //std::vector<double> SCR(NSCR,0.);
  FORTRAN_NAME(dsyev)("V","U", &nf, U.data(),
                      &nf, E.data(),
		      scr, &nscr, &info );
  //SCR, &NSCR, &info );
  assert( info == 0 );
}

void ccdl::csy::dsyevd( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const
{
  int const nf = nfast();
  for ( int i=0; i<nf*nf; ++i )
    U[i] = mData[i];
  int info = 0;
  int liwork= ccdl::sy::iquery_dsyevd(nf);
  std::vector<int> iwork(liwork);
  FORTRAN_NAME(dsyevd)("V","U",
                       &nf,U.data(),
                       &nf,E.data(),
                       scr, &nscr,
                       iwork.data(), &liwork,
                       &info );
  assert( info == 0 );
}

void ccdl::csy::dsyevr( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const
{
  int const nf = nfast();
  int const N2 = nf*nf;
  double * T = scr;
  int const lwork = nscr-N2;
  double * work = scr + N2;
  for ( int i=0; i<N2; ++i )
    T[i] = mData[i];
  int info,M;
  int const liwork = ccdl::sy::iquery_dsyevr(nf);
  std::vector<int> isuppz(liwork);
  int * iwork = isuppz.data() + 2*nf;
  const double VL = 0., VU = 0.;
  const int IL = 1, IU = nf;
  const double ABSTOL = FORTRAN_NAME(dlamch)("Safe minimum");
  FORTRAN_NAME(dsyevr)("V","A","U",
                       &nf,T,&nf,
                       &VL,&VU,
                       &IL,&IU,
                       &ABSTOL,
                       &M,E.data(),U.data(),&nf,
                       isuppz.data(),
                       work,&lwork,
                       iwork,&liwork,
                       &info );
  assert( info == 0 );
}



#define EVMAX 70
#define EVRMAX 440

int ccdl::sy::query_eigen( int const N )
{
  int n = 1;
  if ( N < EVMAX )
    { n = ccdl::sy::query_dsyev(N); }
  else if ( N < EVRMAX )
    { n = ccdl::sy::query_dsyevr(N); }
  else
    { n = ccdl::sy::query_dsyevd(N); };
  return n;
}

void ccdl::sy::eigen( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const
{
  if ( nfast() < EVMAX )
    { dsyev(E,U,nscr,scr); }
  else if ( nfast() < EVRMAX )
    { dsyevr(E,U,nscr,scr); }
  else
    { dsyevd(E,U,nscr,scr); };
}

void ccdl::sy::eigen( ccdl::di & E, ccdl::ge & U ) const
{
  if ( nfast() < EVMAX )
    { dsyev(E,U); }
  else if ( nfast() < EVRMAX )
    { dsyevr(E,U); }
  else
    { dsyevd(E,U); };
}



void ccdl::csy::eigen( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const
{
  if ( nfast() < EVMAX )
    { dsyev(E,U,nscr,scr); }
  else if ( nfast() < EVRMAX )
    { dsyevr(E,U,nscr,scr); }
  else
    { dsyevd(E,U,nscr,scr); };
}

void ccdl::csy::eigen( ccdl::di & E, ccdl::ge & U ) const
{
  if ( nfast() < EVMAX )
    { dsyev(E,U); }
  else if ( nfast() < EVRMAX )
    { dsyevr(E,U); }
  else
    { dsyevd(E,U); };
}


std::ostream & operator<< ( std::ostream & cout, ccdl::cv const a )
{
  for ( int i=0; i<a.size(); ++i )
    cout << std::setw(11) << std::setprecision(3) << std::scientific
         << a[i];
  return cout;
}

std::ostream & operator<< ( std::ostream & cout, ccdl::cge const a )
{
  for ( int i=0; i<a.nfast(); ++i )
    {
      for ( int j=0; j<a.nslow(); ++j )
	cout << std::setw(11) << std::setprecision(3) << std::scientific
	     << a[i+j*a.nfast()];
      cout << "\n";
    };
  return cout;
}

std::ostream & operator<< ( std::ostream & cout, ccdl::csy const a )
{
  for ( int i=0; i<a.nfast(); ++i )
    {
      for ( int j=0; j<a.nslow(); ++j )
	cout << std::setw(11) << std::setprecision(3) << std::scientific
	     << a[i+j*a.nfast()];
      cout << "\n";
    };
  return cout;
}

std::ostream & operator<< ( std::ostream & cout, ccdl::cdi const a )
{
  for ( int i=0; i<a.nfast(); ++i )
    {
      for ( int j=0; j<a.nslow(); ++j )
	{
	  cout << std::setw(11) << std::setprecision(3) << std::scientific;
	  if ( i != j )
	    { cout << 0.; }
	  else
	    { cout << a[i+j*a.nfast()];}
	};
      cout << "\n";
    };
  return cout;
}


