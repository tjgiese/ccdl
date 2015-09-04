#ifndef _BLAS_FortranInterface_H
#define _BLAS_FortranInterface_H

#define F77_CALL(x)    x ## _
#define FORTRAN_NAME(x)    F77_CALL(x)



#ifdef  __cplusplus
extern "C" {
#endif

  /** \private */
  typedef struct 
  {
    double r;
    double i;
  } ComplexT;
  
#ifdef  __cplusplus
}
#endif




#ifdef  __cplusplus
extern "C" {
#endif


/* Double Precision Level 1 BLAS */

/** \private */
extern double /* DASUM - sum of absolute values of a one-dimensional array */
FORTRAN_NAME(dasum)(const int *n, double const *dx, const int *incx);

/** \private */
extern void   /* DAXPY - replace y by alpha*x + y */
FORTRAN_NAME(daxpy)(const int *n, double const *alpha,
		double const *dx, const int *incx,
		double *dy, const int *incy);

/** \private */
extern void   /* DCOPY - copy x to y */
FORTRAN_NAME(dcopy)(const int *n, double const *dx, const int *incx,
		double *dy, const int *incy);

/** \private */
extern double /* DDOT - inner product of x and y */
FORTRAN_NAME(ddot)(const int *n, double const *dx, const int *incx,
	       double const *dy, const int *incy);

/** \private */
extern double /* DNRM2 - 2-norm of a vector */
FORTRAN_NAME(dnrm2)(const int *n, double const *dx, const int *incx);

/** \private */
extern void   /* DROT - apply a Given's rotation */
FORTRAN_NAME(drot)(const int *n, double *dx, const int *incx,
	       double *dy, const int *incy, double const *c, double const *s);

/** \private */
extern void   /* DROTG - generate a Given's rotation */
FORTRAN_NAME(drotg)(double const *a, double const *b, double *c, double *s);

/** \private */
extern void   /* DROTM - apply a modified Given's rotation */
FORTRAN_NAME(drotm)(const int *n, double *dx, const int *incx,
		double *dy, const int *incy, double const *dparam);

/** \private */
extern void   /* DROTMG - generate a modified Given's rotation */
FORTRAN_NAME(drotmg)(double const *dd1, double const *dd2, double const *dx1,
		 double const *dy1, double *param);

/** \private */
extern void   /* DSCAL - scale a one-dimensional array */
FORTRAN_NAME(dscal)(const int *n, double const *alpha, double *dx, const int *incx);

/** \private */
extern void   /* DSWAP - interchange one-dimensional arrays */
FORTRAN_NAME(dswap)(const int *n, double *dx, const int *incx,
		double *dy, const int *incy);

/** \private */
extern int    /* IDAMAX - return the index of the element with max abs value */
FORTRAN_NAME(idamax)(const int *n, double const *dx, const int *incx);

/* Double Precision Level 2 BLAS */

/* DGBMV - perform one of the matrix-vector operations */
/* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y, */

/** \private */
extern void
FORTRAN_NAME(dgbmv)(const char *trans, const int *m, const int *n,
		const int *kl,const int *ku,
		double const *alpha, double const *a, const int *lda,
		double const *x, const int *incx,
		double const *beta, double *y, const int *incy);
/* DGEMV - perform one of the matrix-vector operations */
/* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y,  */

/** \private */
extern void
FORTRAN_NAME(dgemv)(const char *trans, const int *m, const int *n,
		double const *alpha, double const *a, const int *lda,
		double const *x, const int *incx, double const *beta,
		double *y, const int *incy);
/* DSBMV - perform the matrix-vector operation */
/* y := alpha*A*x + beta*y, */

/** \private */
extern void
FORTRAN_NAME(dsbmv)(const char *uplo, const int *n, const int *k,
		double const *alpha, double const *a, const int *lda,
		double const *x, const int *incx,
		double const *beta, double *y, const int *incy);
/* DSPMV - perform the matrix-vector operation */
/* y := alpha*A*x + beta*y, */

/** \private */
extern void
FORTRAN_NAME(dspmv)(const char *uplo, const int *n,
		double const *alpha, double const *ap,
		double const *x, const int *incx,
		double const *beta, double *y, const int *incy);

/* DSYMV - perform the matrix-vector operation */
/*  y := alpha*A*x + beta*y, */

/** \private */
extern void
FORTRAN_NAME(dsymv)(const char *uplo, const int *n, double const *alpha,
		double const *a, const int *lda,
		double const *x, const int *incx,
		double const *beta, double *y, const int *incy);
/* DTBMV - perform one of the matrix-vector operations */
/* x := A*x, or x := A'*x, */

/** \private */
extern void
FORTRAN_NAME(dtbmv)(const char *uplo, const char *trans,
		const char *diag, const int *n, const int *k,
		double const *a, const int *lda,
		double *x, const int *incx);
/* DTPMV - perform one of the matrix-vector operations */
/* x := A*x, or x := A'*x, */

/** \private */
extern void
FORTRAN_NAME(dtpmv)(const char *uplo, const char *trans, const char *diag,
		const int *n, double const *ap,
		double *x, const int *incx);
/* DTRMV - perform one of the matrix-vector operations  */
/* x := A*x, or x := A'*x, */

/** \private */
extern void
FORTRAN_NAME(dtrmv)(const char *uplo, const char *trans, const char *diag,
		const int *n, double const *a, const int *lda,
		double *x, const int *incx);
/* DTBSV - solve one of the systems of equations */
/* A*x = b, or A'*x = b, */

/** \private */
extern void
FORTRAN_NAME(dtbsv)(const char *uplo, const char *trans,
		const char *diag, const int *n, const int *k,
		double const *a, const int *lda,
		double *x, const int *incx);
/* DTPSV - solve one of the systems of equations */
/* A*x = b, or A'*x = b, */
extern void
FORTRAN_NAME(dtpsv)(const char *uplo, const char *trans,
		const char *diag, const int *n,
		double const *ap, double *x, const int *incx);
/* DTRSV - solve one of the systems of equations */
/* A*x = b, or A'*x = b, */

/** \private */
extern void
FORTRAN_NAME(dtrsv)(const char *uplo, const char *trans,
		const char *diag, const int *n,
		double const *a, const int *lda,
		double *x, const int *incx);
/* DGER - perform the rank 1 operation   A := alpha*x*y' + A */

/** \private */
extern void
FORTRAN_NAME(dger)(const int *m, const int *n, double const *alpha,
	       double const *x, const int *incx,
	       double const *y, const int *incy,
	       double *a, const int *lda);
/* DSYR - perform the symmetric rank 1 operation A := alpha*x*x' + A */

/** \private */
extern void
FORTRAN_NAME(dsyr)(const char *uplo, const int *n, double const *alpha,
	       double const *x, const int *incx,
	       double *a, const int *lda);
/* DSPR - perform the symmetric rank 1 operation A := alpha*x*x' + A */


/** \private */
extern void
FORTRAN_NAME(dspr)(const char *uplo, const int *n, double const *alpha,
	       double const *x, const int *incx, double *ap);
/* DSYR2 - perform the symmetric rank 2 operation */
/* A := alpha*x*y' + alpha*y*x' + A, */


/** \private */
extern void
FORTRAN_NAME(dsyr2)(const char *uplo, const int *n, double const *alpha,
		double const *x, const int *incx,
		double const *y, const int *incy,
		double *a, const int *lda);
/* DSPR2 - perform the symmetric rank 2 operation */
/* A := alpha*x*y' + alpha*y*x' + A,  */


/** \private */
extern void
FORTRAN_NAME(dspr2)(const char *uplo, const int *n, double const *alpha,
		double const *x, const int *incx,
		double const *y, const int *incy, double *ap);

/* Double Precision Level 3 BLAS */

/* DGEMM - perform one of the matrix-matrix operations    */
/* C := alpha*op( A )*op( B ) + beta*C */


/** \private */
extern void
FORTRAN_NAME(dgemm)(const char *transa, const char *transb, const int *m,
		const int *n, const int *k, double const *alpha,
		double const *a, const int *lda,
		double const *b, const int *ldb,
		double const *beta, double *c, const int *ldc);
/* DTRSM - solve one of the matrix equations  */
/* op(A)*X = alpha*B, or  X*op(A) = alpha*B  */


/** \private */
extern void
FORTRAN_NAME(dtrsm)(const char *side, const char *uplo,
		const char *transa, const char *diag,
		const int *m, const int *n, double const *alpha,
		double const *a, const int *lda,
		double *b, const int *ldb);
/* DTRMM - perform one of the matrix-matrix operations */
/* B := alpha*op( A )*B, or B := alpha*B*op( A ) */


/** \private */
extern void
FORTRAN_NAME(dtrmm)(const char *side, const char *uplo, const char *transa,
		const char *diag, const int *m, const int *n,
		double const *alpha, double const *a, const int *lda,
		double *b, const int *ldb);
/* DSYMM - perform one of the matrix-matrix operations   */
/*  C := alpha*A*B + beta*C, */


/** \private */
extern void
FORTRAN_NAME(dsymm)(const char *side, const char *uplo, const int *m,
		const int *n, double const *alpha,
		double const *a, const int *lda,
		double const *b, const int *ldb,
		double const *beta, double *c, const int *ldc);
/* DSYRK - perform one of the symmetric rank k operations */
/* C := alpha*A*A' + beta*C or C := alpha*A'*A + beta*C */


/** \private */
extern void
FORTRAN_NAME(dsyrk)(const char *uplo, const char *trans,
		const int *n, const int *k,
		double const *alpha, double const *a, const int *lda,
		double const *beta, double *c, const int *ldc);
/* DSYR2K - perform one of the symmetric rank 2k operations */
/* C := alpha*A*B' + alpha*B*A' + beta*C or */
/* C := alpha*A'*B + alpha*B'*A + beta*C */


/** \private */
extern void
FORTRAN_NAME(dsyr2k)(const char *uplo, const char *trans,
		 const int *n, const int *k,
		 double const *alpha, double const *a, const int *lda,
		 double const *b, const int *ldb,
		 double const *beta, double *c, const int *ldc);
/*
  LSAME is a LAPACK support routine, not part of BLAS
*/

/* Double complex BLAS routines added for 2.3.0 */
/* #ifdef HAVE_FORTRAN_DOUBLE_COMPLEX */
    

/** \private */
extern double
    FORTRAN_NAME(dcabs1)(double *z);
    

/** \private */
extern double
    FORTRAN_NAME(dzasum)(int *n, ComplexT *zx, int *incx);
    

/** \private */
extern double
    FORTRAN_NAME(dznrm2)(int *n, ComplexT *x, int *incx);
    

/** \private */
extern int
    FORTRAN_NAME(izamax)(int *n, ComplexT *zx, int *incx);
    

/** \private */
extern void
    FORTRAN_NAME(zaxpy)(int *n, ComplexT *za, ComplexT *zx,
		    int *incx, ComplexT *zy, int *incy);
    

/** \private */
extern void
    FORTRAN_NAME(zcopy)(int *n, ComplexT *zx, int *incx,
		    ComplexT *zy, int *incy);
    

/** \private */
extern ComplexT
    FORTRAN_NAME(zdotc)(ComplexT * ret_val, int *n,
		    ComplexT *zx, int *incx, ComplexT *zy, int *incy);
    

/** \private */
extern ComplexT
    FORTRAN_NAME(zdotu)(ComplexT * ret_val, int *n,
		    ComplexT *zx, int *incx, ComplexT *zy, int *incy);
    

/** \private */
extern void
    FORTRAN_NAME(zdrot)(int *n, ComplexT *zx, int *incx, ComplexT *zy,
		int *incy, double *c, double *s);
    

/** \private */
extern void
    FORTRAN_NAME(zdscal)(int *n, double *da, ComplexT *zx, int *incx);
    

/** \private */
extern void
    FORTRAN_NAME(zgbmv)(char *trans, int *m, int *n, int *kl,
		    int *ku, ComplexT *alpha, ComplexT *a, int *lda,
		    ComplexT *x, int *incx, ComplexT *beta, ComplexT *y,
		    int *incy);
    

/** \private */
extern void
    FORTRAN_NAME(zgemm)(const char *transa, const char *transb, const int *m,
		    const int *n, const int *k, const ComplexT *alpha,
		    const ComplexT *a, const int *lda,
		    const ComplexT *b, const int *ldb,
		    const ComplexT *beta, ComplexT *c, const int *ldc);
    

/** \private */
extern void
    FORTRAN_NAME(zgemv)(char *trans, int *m, int *n, ComplexT *alpha,
		    ComplexT *a, int *lda, ComplexT *x, int *incx,
		    ComplexT *beta, ComplexT *y, int * incy);
    

/** \private */
extern void
    FORTRAN_NAME(zgerc)(int *m, int *n, ComplexT *alpha, ComplexT *x,
		    int *incx, ComplexT *y, int *incy, ComplexT *a, int *lda);
    

/** \private */
extern void
    FORTRAN_NAME(zgeru)(int *m, int *n, ComplexT *alpha, ComplexT *x,
		    int *incx, ComplexT *y, int *incy, ComplexT *a, int *lda);
    

/** \private */
extern void
    FORTRAN_NAME(zhbmv)(char *uplo, int *n, int *k, ComplexT *alpha,
		    ComplexT *a, int *lda, ComplexT *x, int *incx,
		    ComplexT *beta, ComplexT *y, int *incy);
    

/** \private */
extern void
    FORTRAN_NAME(zhemm)(char *side, char *uplo, int *m, int *n,
		    ComplexT *alpha, ComplexT *a, int *lda, ComplexT *b,
		    int *ldb, ComplexT *beta, ComplexT *c, int *ldc);
    

/** \private */
extern void
    FORTRAN_NAME(zhemv)(char *uplo, int *n, ComplexT *alpha, ComplexT *a,
		    int *lda, ComplexT *x, int *incx, ComplexT *beta,
		    ComplexT *y, int *incy);
    

/** \private */
extern void
    FORTRAN_NAME(zher)(char *uplo, int *n, double *alpha, ComplexT *x,
		   int *incx, ComplexT *a, int *lda);
    

/** \private */
extern void
    FORTRAN_NAME(zher2)(char *uplo, int *n, ComplexT *alpha, ComplexT *x,
		    int *incx, ComplexT *y, int *incy, ComplexT *a, int *lda);
    

/** \private */
extern void
    FORTRAN_NAME(zher2k)(char *uplo, char *trans, int *n, int *k,
		     ComplexT *alpha, ComplexT *a, int *lda, ComplexT *b,
		     int *ldb, double *beta, ComplexT *c, int *ldc);
    

/** \private */
extern void
    FORTRAN_NAME(zherk)(char *uplo, char *trans, int *n, int *k,
		    double *alpha, ComplexT *a, int *lda, double *beta,
		    ComplexT *c, int *ldc);
    

/** \private */
extern void
    FORTRAN_NAME(zhpmv)(char *uplo, int *n, ComplexT *alpha, ComplexT *ap,
		    ComplexT *x, int *incx, ComplexT * beta, ComplexT *y,
		    int *incy);
    

/** \private */
extern void
    FORTRAN_NAME(zhpr)(char *uplo, int *n, double *alpha,
		   ComplexT *x, int *incx, ComplexT *ap);
    

/** \private */
extern void
    FORTRAN_NAME(zhpr2)(char *uplo, int *n, ComplexT *alpha, ComplexT *x,
		    int *incx, ComplexT *y, int *incy, ComplexT *ap);
    

/** \private */
extern void
    FORTRAN_NAME(zrotg)(ComplexT *ca, ComplexT *cb, double *c, ComplexT *s);
    

/** \private */
extern void
    FORTRAN_NAME(zscal)(int *n, ComplexT *za, ComplexT *zx, int *incx);
    

/** \private */
extern void
    FORTRAN_NAME(zswap)(int *n, ComplexT *zx, int *incx, ComplexT *zy, int *incy);
    

/** \private */
extern void
    FORTRAN_NAME(zsymm)(char *side, char *uplo, int *m, int *n,
		    ComplexT *alpha, ComplexT *a, int *lda, ComplexT *b,
		    int *ldb, ComplexT *beta, ComplexT *c, int *ldc);
    

/** \private */
extern void
    FORTRAN_NAME(zsyr2k)(char *uplo, char *trans, int *n, int *k,
		     ComplexT *alpha, ComplexT *a, int *lda, ComplexT *b,
		     int *ldb, ComplexT *beta, ComplexT *c, int *ldc);
    

/** \private */
extern void
    FORTRAN_NAME(zsyrk)(char *uplo, char *trans, int *n, int *k,
		    ComplexT *alpha, ComplexT *a, int *lda,
		    ComplexT *beta, ComplexT *c, int *ldc);
    

/** \private */
extern void
    FORTRAN_NAME(ztbmv)(char *uplo, char *trans, char *diag, int *n, int *k,
		    ComplexT *a, int *lda, ComplexT *x, int *incx);
    

/** \private */
extern void
    FORTRAN_NAME(ztbsv)(char *uplo, char *trans, char *diag, int *n, int *k,
		    ComplexT *a, int *lda, ComplexT *x, int *incx);
    

/** \private */
extern void
    FORTRAN_NAME(ztpmv)(char *uplo, char *trans, char *diag, int *n,
		    ComplexT *ap, ComplexT *x, int *incx);
    

/** \private */
extern void
    FORTRAN_NAME(ztpsv)(char *uplo, char *trans, char *diag, int *n,
		    ComplexT *ap, ComplexT *x, int *incx);
    

/** \private */
extern void
    FORTRAN_NAME(ztrmm)(char *side, char *uplo, char *transa, char *diag,
		    int *m, int *n, ComplexT *alpha, ComplexT *a,
		    int *lda, ComplexT *b, int *ldb);
    

/** \private */
extern void
    FORTRAN_NAME(ztrmv)(char *uplo, char *trans, char *diag, int *n,
		    ComplexT *a, int *lda, ComplexT *x, int *incx);
    

/** \private */
extern void
    FORTRAN_NAME(ztrsm)(char *side, char *uplo, char *transa, char *diag,
		    int *m, int *n, ComplexT *alpha, ComplexT *a,
		    int *lda, ComplexT *b, int *ldb);
    

/** \private */
extern void
    FORTRAN_NAME(ztrsv)(char *uplo, char *trans, char *diag, int *n,
		    ComplexT *a, int *lda, ComplexT *x, int *incx);
/* #endif */

#ifdef  __cplusplus
}
#endif

#endif 
