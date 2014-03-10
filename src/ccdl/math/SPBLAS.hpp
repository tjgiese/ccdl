#ifndef SPBLAS_H
#define SPBLAS_H

#include <cmath>
#include <vector>
#include <algorithm>


// # define MKL

#ifdef __cplusplus
extern "C" {
#endif

  void mkl_dcsrmv( char const * transa, 
		   int const * nrowA, int const * ncolA, 
		   double const * alpha, char const * matdescra,
		   double const * Adata, 
		   int const * colidx, 
		   int const * rowoff, int const * rowend,
		   double const * Xdata, double const * beta, 
		   double * Ydata );
  
  void mkl_dcsrmm( char const * transa, 
		   int const * nrowA, 
		   int const * ncolC, int const * ncolA, 
		   double const * alpha, char const * matdescra, 
		   double const * Adata, 
		   int const * colidx, 
		   int const * rowoff, int const * rowend,
		   double const * Bdata, int const * nrowB, 
		   double const * beta, double * Cdata, 
		   int const * nrowC );
  
  void dcsrmm( const int transa, const int m, const int n, 
	       const int k,
	       const double alpha, const int descra[], 
	       const double val[],
	       const int indx[], const int pntrb[], const int pntre[],
	       const double b[], const int ldb,
	       const double beta, double c[], const int ldc );


  void mkl_dcsrsm( char const * transa, int const * m, int const * n, 
		   double const * alpha, char const * matdescra, 
		   double const * val, int const * indx, 
		   int const * pntrb, int const * pntre, 
		   double const * b, int const * ldb, 
		   double * c, int const * ldc );

  
  void  dcsrsm( const int transa, const int m, const int n,
		const int unitd, const double dv[],
		const double alpha, const int descra[], const double val[],
		const int indx[], const int pntrb[], const int pntre[],
		const double b[], const int ldb,
		const double beta, double c[], const int ldc,
		double work[], const int lwork );
  

#ifdef __cplusplus
}
#endif


inline void sparse_symv
( int const n, 
  double const alpha, 
  double const * a,
  int const * colidx,
  int const * rowoff,
  double const beta,
  double const * b,
  double * c )
{
#ifdef MKL
  char matdescra[] = "SLNC1 ";
  mkl_dcsrmv("N", &n, &n, &alpha, matdescra, a, colidx, rowoff, rowoff+1, b, &beta, c);
#else
  int descra[] = {1,1,0,0,1};
  dcsrmm(0,n,1,n,alpha,descra,a,colidx,rowoff,rowoff+1,b,n,beta,c,n);
#endif
}

inline void sparse_gemv
( int const nrowA, int const ncolA, 
  double const alpha, 
  double const * a,
  int const * colidx,
  int const * rowoff,
  double const beta,
  double const * b,
  double * c )
{
#ifdef MKL
  char matdescra[] = "GLNC1 ";
  mkl_dcsrmv("N", &nrowA, &ncolA, &alpha, matdescra, a, colidx, rowoff, rowoff+1, b, &beta, c);
#else
  int descra[] = {0,1,0,0,1};
  dcsrmm(0,nrowA,1,ncolA,alpha,descra,a,colidx,rowoff,rowoff+1,b,ncolA,beta,c,nrowA);
#endif
}


inline void sparse_gtmv
( int const nrowA, int const ncolA, 
  double const alpha, 
  double const * a,
  int const * colidx,
  int const * rowoff,
  double const beta,
  double const * b,
  double * c )
{
#ifdef MKL
  char matdescra[] = "GLNC1 ";
  mkl_dcsrmv("T", &nrowA, &ncolA, &alpha, matdescra, a, colidx, rowoff, rowoff+1, b, &beta, c);
#else
  int descra[] = {0,1,0,0,1};
  dcsrmm(1,ncolA,1,nrowA,alpha,descra,a,colidx,rowoff,rowoff+1,b,nrowA,beta,c,ncolA);
#endif
}

inline void sparse_symm
(
#ifdef MKL
 int const nrowA, int const ncolA,
#else
 int const, int const,
#endif 
  int const nrowB, 
  int const nrowC, int const ncolC, 
  double const alpha, 
  double const * a,
  int const * colidx,
  int const * rowoff,
  double const beta,
  double const * b,
  double * c )
{
#ifdef MKL
  char matdescra[] = "SLNC1 ";
  mkl_dcsrmm("N", &nrowA, &ncolC, &ncolA, 
	     &alpha, matdescra, a, 
	     colidx, 
	     rowoff, rowoff+1, 
	     b, &nrowB, 
	     &beta, c, &nrowC );
#else
  int descra[] = {1,1,0,0,1};
  dcsrmm(0,nrowC,ncolC,nrowB,alpha,descra,a,colidx,rowoff,rowoff+1,b,nrowB,beta,c,nrowC);
#endif
}


inline void sparse_gemm
( 
#ifdef MKL
 int const nrowA, int const ncolA, 
#else
 int const, int const,
#endif
  int const nrowB, 
  int const nrowC, int const ncolC, 
  double const alpha, 
  double const * a,
  int const * colidx,
  int const * rowoff,
  double const beta,
  double const * b,
  double * c )
{
#ifdef MKL
  char matdescra[] = "GLNC1 ";
  mkl_dcsrmm("N", &nrowA, &ncolC, &ncolA, 
	     &alpha, matdescra, a, 
	     colidx, 
	     rowoff, rowoff+1, 
	     b, &nrowB, 
	     &beta, c, &nrowC );
#else
  int descra[] = {0,1,0,0,1};
  dcsrmm(0,nrowC,ncolC,nrowB,alpha,descra,a,colidx,rowoff,rowoff+1,b,nrowB,beta,c,nrowC);
#endif
}

inline void sparse_gtmm
(
#ifdef MKL 
 int const nrowA, int const ncolA, 
#else
 int const, int const,
#endif
  int const nrowB, 
  int const nrowC, int const ncolC, 
  double const alpha, 
  double const * a,
  int const * colidx,
  int const * rowoff,
  double const beta,
  double const * b,
  double * c )
{
#ifdef MKL
  char matdescra[] = "GLNC1 ";
  mkl_dcsrmm("T", &nrowA, &ncolC, &ncolA, 
	     &alpha, matdescra, a, 
	     colidx, 
	     rowoff, rowoff+1, 
	     b, &nrowB, 
	     &beta, c, &nrowC );
#else
  int descra[] = {0,1,0,0,1};
  dcsrmm(1,nrowC,ncolC,nrowB,alpha,descra,a,colidx,rowoff,rowoff+1,b,nrowB,beta,c,nrowC);
#endif
}



inline void sparse_trisolve_Ax_eq_b
(
 int const n, 
 double const alpha, 
 double const * a,
 int const * colidx,
 int const * rowoff,
 double const * b,
 double * x )
{
#ifdef MKL
  int const o = 1;
  char matdescra[] = "TLNC  ";
  mkl_dcsrsm("N", &n, &o, 
	     &alpha, matdescra, a, 
	     colidx, rowoff, rowoff+1, 
	     b, &o, x, &o );
#else
  int descra[] = {3,1,0,0,1};
  double * dv = NULL;
  std::vector<double> work(n*n,0.);
  dcsrsm( 0, n, 1, 1, dv, alpha, descra, a, colidx, rowoff, rowoff+1,
	  b, 1, 0., x, 1, work.data(), work.size() );
#endif
}


inline void sparse_trisolve_Atx_eq_b
(
 int const n, 
 double const alpha, 
 double const * a,
 int const * colidx,
 int const * rowoff,
 double const * b,
 double * x )
{
#ifdef MKL
  int const o = 1;
  char matdescra[] = "TLNC  ";
  mkl_dcsrsm("T", &n, &o, 
	     &alpha, matdescra, a, 
	     colidx, rowoff, rowoff+1, 
	     b, &o, x, &o );
#else
  int descra[] = {3,1,0,0,1};
  double * dv = NULL;
  std::vector<double> work(n*n,0.);
  dcsrsm( 1, n, 1, 1, dv, alpha, descra, a, colidx, rowoff, rowoff+1,
	  b, 1, 0., x, 1, work.data(), work.size() );
#endif
}


#endif



