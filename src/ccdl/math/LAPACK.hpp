#ifndef _LAPACK_FortranInterface_H
#define _LAPACK_FortranInterface_H

#include "BLAS.hpp"

#ifdef	__cplusplus
extern "C" {
#endif

/* Double precision BiDiagonal matrices */

/* DBDSQR - compute the singular value decomposition (SVD) of a real */
/* N-by-N (upper or lower) bidiagonal matrix B */


/** \private */
extern void
FORTRAN_NAME(dbdsqr)(const char* uplo, const int* n, const int* ncvt,
		     const int* nru, const int* ncc, double* d, double* e,
		     double* vt, const int* ldvt, double* u, const int* ldu,
		     double* c, const int* ldc, double* work, int* info);
  /* DDISNA - compute the reciprocal condition numbers for the */
  /* eigenvectors of a real symmetric or complex Hermitian matrix or */
  /* for the left or right singular vectors of a general m-by-n */
  /* matrix */


  /** \private */
  extern void
  FORTRAN_NAME(ddisna)(const char* job, const int* m, const int* n,
		       double* d, double* sep, int* info);
  
/* Double precision General Banded matrices */

/* DGBBRD - reduce a real general m-by-n band matrix A to upper */
/* bidiagonal form B by an orthogonal transformation  */


/** \private */
extern void
FORTRAN_NAME(dgbbrd)(const char* vect, const int* m, const int* n,
		 const int* ncc, const int* kl, const int* ku,
		 double* ab, const int* ldab,
		 double* d, double* e, double* q,
		 const int* ldq, double* pt, const int* ldpt,
		 double* c, const int* ldc,
		 double* work, int* info);
/* DGBCON - estimate the reciprocal of the condition number of a */
/* real general band matrix A, in either the 1-norm or the */
/* infinity-norm */


/** \private */
extern void
FORTRAN_NAME(dgbcon)(const char* norm, const int* n, const int* kl,
		 const int* ku, double* ab, const int* ldab,
		 int* ipiv, double const* anorm, double* rcond,
		 double* work, int* iwork, int* info);
/* DGBEQU - compute row and column scalings intended to equilibrate */
/* an M-by-N band matrix A and reduce its condition number */


/** \private */
extern void
FORTRAN_NAME(dgbequ)(const int* m, const int* n, const int* kl, const int* ku,
		 double* ab, const int* ldab, double* r, double* c,
		 double* rowcnd, double* colcnd, double* amax, int* info);
/* DGBRFS - improve the computed solution to a system of linear */
/* equations when the coefficient matrix is banded, and provides */
/* error bounds and backward error estimates for the solution */


/** \private */
extern void
FORTRAN_NAME(dgbrfs)(const char* trans, const int* n, const int* kl,
		 const int* ku, const int* nrhs, double* ab,
		 const int* ldab, double* afb, const int* ldafb,
		 int* ipiv, double* b, const int* ldb,
		 double* x, const int* ldx, double* ferr, double* berr,
		 double* work, int* iwork, int* info);
/* DGBSV - compute the solution to a real system of linear */
/* equations A * X = B, where A is a band matrix of order N with */
/* KL subdiagonals and KU superdiagonals, and X and B are */
/* N-by-NRHS matrices */


/** \private */
extern void
FORTRAN_NAME(dgbsv)(const int* n, const int* kl,const int* ku,
		const int* nrhs, double* ab, const int* ldab,
		int* ipiv, double* b, const int* ldb, int* info);
/* DGBSVX - use the LU factorization to compute the solution to a */
/* real system of linear equations A * X = B or A**T * X = B */


/** \private */
extern void
FORTRAN_NAME(dgbsvx)(const int* fact, const char* trans,
		 const int* n, const int* kl,const int* ku,
		 const int* nrhs, double* ab, const int* ldab,
		 double* afb, const int* ldafb, int* ipiv,
		 const char* equed, double* r, double* c,
		 double* b, const int* ldb,
		 double* x, const int* ldx,
		 double* rcond, double* ferr, double* berr,
		 double* work, int* iwork, int* info);
/* DGBTF2 - compute an LU factorization of a real m-by-n band */
/* matrix A using partial pivoting with row interchanges */


/** \private */
extern void
FORTRAN_NAME(dgbtf2)(const int* m, const int* n, const int* kl,const int* ku,
		 double* ab, const int* ldab, int* ipiv, int* info);
/* DGBTRF - compute an LU factorization of a real m-by-n band */
/* matrix A using partial pivoting with row interchanges */


/** \private */
extern void
FORTRAN_NAME(dgbtrf)(const int* m, const int* n, const int* kl,const int* ku,
		  double* ab, const int* ldab, int* ipiv, int* info);
/* DGBTRS - solve a system of linear equations	A * X = B or  */
/* A' * X = B with a general band matrix A using the LU */
/* factorization computed by DGBTRF */


/** \private */
extern void
FORTRAN_NAME(dgbtrs)(const char* trans, const int* n,
		 const int* kl, const int* ku, const int* nrhs,
		 double const* ab, const int* ldab, const int* ipiv,
		 double* b, const int* ldb, int* info);

/* Double precision GEneral matrices */

/* DGEBAK - form the right or left eigenvectors of a real general */
/* matrix by backward transformation on the computed eigenvectors */
/* of the balanced matrix output by DGEBAL  */


/** \private */
extern void
FORTRAN_NAME(dgebak)(const char* job, const char* side, const int* n,
		 const int* ilo, const int* ihi, double* scale,
		 const int* m, double* v, const int* ldv, int* info);
/* DGEBAL - balance a general real matrix A */


/** \private */
extern void
FORTRAN_NAME(dgebal)(const char* job, const int* n, double* a, const int* lda,
		  int* ilo, int* ihi, double* scale, int* info);
/* DGEBD2 - reduce a real general m by n matrix A to upper or */
/* lower bidiagonal form B by an orthogonal transformation */


/** \private */
extern void
FORTRAN_NAME(dgebd2)(const int* m, const int* n, double* a, const int* lda,
		 double* d, double* e, double* tauq, double* taup,
		 double* work, int* info);
/* DGEBRD - reduce a general real M-by-N matrix A to upper or */
/* lower bidiagonal form B by an orthogonal transformation */


/** \private */
extern void
FORTRAN_NAME(dgebrd)(const int* m, const int* n, double* a, const int* lda,
		 double* d, double* e, double* tauq, double* taup,
		 double* work, const int* lwork, int* info);
/* DGECON - estimate the reciprocal of the condition number of a */
/* general real matrix A, in either the 1-norm or the */
/* infinity-norm, using the LU factorization computed by DGETRF */


/** \private */
extern void
FORTRAN_NAME(dgecon)(const char* norm, const int* n,
		 double const* a, const int* lda,
		 double const* anorm, double* rcond,
		 double* work, int* iwork, int* info);
/* DGEEQU - compute row and column scalings intended to equilibrate */
/* an M-by-N matrix A and reduce its condition number */


/** \private */
extern void
FORTRAN_NAME(dgeequ)(const int* m, const int* n, double* a, const int* lda,
		 double* r, double* c, double* rowcnd, double* colcnd,
		 double* amax, int* info);
/* DGEES - compute for an N-by-N real nonsymmetric matrix A, the */
/* eigenvalues, the real Schur form T, and, optionally, the matrix */
/* of Schur vectors Z */


/** \private */
extern void
FORTRAN_NAME(dgees)(const char* jobvs, const char* sort,
		int (*select)(double const*, double const*),
		const int* n, double* a, const int* lda,
		int* sdim, double* wr, double* wi,
		double* vs, const int* ldvs,
		double* work, const int* lwork, int* bwork, int* info);
/* DGEESX - compute for an N-by-N real nonsymmetric matrix A, the */
/* eigenvalues, the real Schur form T, and, optionally, the matrix */
/* of Schur vectors Z */


/** \private */
extern void
FORTRAN_NAME(dgeesx)(const char* jobvs, const char* sort,
		 int (*select)(double const*, double const*),
		 const char* sense, const int* n, double* a,
		 const int* lda, int* sdim, double* wr, double* wi,
		 double* vs, const int* ldvs, double* rconde,
		 double* rcondv, double* work, const int* lwork,
		 int* iwork, const int* liwork, int* bwork, int* info);
/* DGEEV - compute for an N-by-N real nonsymmetric matrix A, the */
/* eigenvalues and, optionally, the left and/or right eigenvectors */


/** \private */
extern void
FORTRAN_NAME(dgeev)(const char* jobvl, const char* jobvr,
		const int* n, double* a, const int* lda,
		double* wr, double* wi, double* vl, const int* ldvl,
		double* vr, const int* ldvr,
		double* work, const int* lwork, int* info);
/* DGEEVX - compute for an N-by-N real nonsymmetric matrix A, the */
/* eigenvalues and, optionally, the left and/or right eigenvectors */


/** \private */
extern void
FORTRAN_NAME(dgeevx)(const char* balanc, const char* jobvl, const char* jobvr,
		 const char* sense, const int* n, double* a, const int* lda,
		 double* wr, double* wi, double* vl, const int* ldvl,
		 double* vr, const int* ldvr, int* ilo, int* ihi,
		 double* scale, double* abnrm, double* rconde, double* rcondv,
		 double* work, const int* lwork, int* iwork, int* info);
/* DGEGV - compute for a pair of n-by-n real nonsymmetric */
/* matrices A and B, the generalized eigenvalues (alphar +/- */
/* alphai*i, beta);, and optionally, the left and/or right */
/* generalized eigenvectors (VL and VR); */


/** \private */
extern void
FORTRAN_NAME(dgegv)(const char* jobvl, const char* jobvr,
		const int* n, double* a, const int* lda,
		double* b, const int* ldb,
		double* alphar, double* alphai,
		double const* beta, double* vl, const int* ldvl,
		double* vr, const int* ldvr,
		double* work, const int* lwork, int* info);
/* DGEHD2 - reduce a real general matrix A to upper Hessenberg */
/* form H by an orthogonal similarity transformation */


/** \private */
extern void
FORTRAN_NAME(dgehd2)(const int* n, const int* ilo, const int* ihi,
		 double* a, const int* lda, double* tau,
		 double* work, int* info);
/* DGEHRD - reduce a real general matrix A to upper Hessenberg */
/* form H by an orthogonal similarity transformation */


/** \private */
extern void
FORTRAN_NAME(dgehrd)(const int* n, const int* ilo, const int* ihi,
		 double* a, const int* lda, double* tau,
		 double* work, const int* lwork, int* info);
/* DGELQ2 - compute an LQ factorization of a real m by n matrix A */


/** \private */
extern void
FORTRAN_NAME(dgelq2)(const int* m, const int* n,
		 double* a, const int* lda, double* tau,
		 double* work, int* info);
/* DGELQF - compute an LQ factorization of a real M-by-N matrix A */


/** \private */
extern void
FORTRAN_NAME(dgelqf)(const int* m, const int* n,
		 double* a, const int* lda, double* tau,
		 double* work, const int* lwork, int* info);
/* DGELS - solve overdetermined or underdetermined real linear */
/* systems involving an M-by-N matrix A, or its transpose, using a */
/* QR or LQ factorization of A */


/** \private */
extern void
FORTRAN_NAME(dgels)(const char* trans, const int* m, const int* n,
		const int* nrhs, double* a, const int* lda,
		double* b, const int* ldb,
		double* work, const int* lwork, int* info);
/* DGELSS - compute the minimum norm solution to a real linear */
/* least squares problem */


/** \private */
extern void
FORTRAN_NAME(dgelss)(const int* m, const int* n, const int* nrhs,
		 double* a, const int* lda, double* b, const int* ldb,
		 double* s, double* rcond, int* rank,
		 double* work, const int* lwork, int* info);
/* DGELSY - compute the minimum-norm solution to a real linear */
/* least squares problem */


/** \private */
extern void
FORTRAN_NAME(dgelsy)(const int* m, const int* n, const int* nrhs,
		 double* a, const int* lda, double* b, const int* ldb,
		 int* jpvt, double const* rcond, int* rank,
		 double* work, const int* lwork, int* info);
/* DGEQL2 - compute a QL factorization of a real m by n matrix A */


/** \private */
extern void
FORTRAN_NAME(dgeql2)(const int* m, const int* n, double* a, const int* lda,
		 double* tau, double* work, int* info);
/* DGEQLF - compute a QL factorization of a real M-by-N matrix A */


/** \private */
extern void
FORTRAN_NAME(dgeqlf)(const int* m, const int* n,
		 double* a, const int* lda, double* tau,
		 double* work, const int* lwork, int* info);
/* DGEQP3 - compute a QR factorization with column pivoting of a */
/* real M-by-N matrix A using level 3 BLAS */


/** \private */
extern void
FORTRAN_NAME(dgeqp3)(const int* m, const int* n, double* a, const int* lda,
		 int* jpvt, double* tau, double* work, const int* lwork,
		 int* info);
/* DGEQPF - compute a QR factorization with column pivoting of a */
/* real M-by-N matrix A */


/** \private */
extern void
FORTRAN_NAME(dgeqpf)(const int* m, const int* n, double* a, const int* lda,
		 int* jpvt, double* tau, double* work, int* info);
/* DGEQR2 - compute a QR factorization of a real m by n matrix A */


/** \private */
extern void
FORTRAN_NAME(dgeqr2)(const int* m, const int* n, double* a, const int* lda,
		 double* tau, double* work, int* info);
/* DGEQRF - compute a QR factorization of a real M-by-N matrix A */


/** \private */
extern void
FORTRAN_NAME(dgeqrf)(const int* m, const int* n, double* a, const int* lda,
		 double* tau, double* work, const int* lwork, int* info);
/* DGERFS - improve the computed solution to a system of linear */
/* equations and provides error bounds and backward error */
/* estimates for the solution */


/** \private */
extern void
FORTRAN_NAME(dgerfs)(const char* trans, const int* n, const int* nrhs,
		 double* a, const int* lda, double* af, const int* ldaf,
		 int* ipiv, double* b, const int* ldb,
		 double* x, const int* ldx, double* ferr, double* berr,
		 double* work, int* iwork, int* info);
/* DGERQ2 - compute an RQ factorization of a real m by n matrix A */


/** \private */
extern void
FORTRAN_NAME(dgerq2)(const int* m, const int* n, double* a, const int* lda,
		 double* tau, double* work, int* info);
/* DGERQF - compute an RQ factorization of a real M-by-N matrix A */


/** \private */
extern void
FORTRAN_NAME(dgerqf)(const int* m, const int* n, double* a, const int* lda,
		 double* tau, double* work, const int* lwork, int* info);
/* DGESV - compute the solution to a real system of linear */
/* equations  A * X = B, */


/** \private */
extern void
FORTRAN_NAME(dgesv)(const int* n, const int* nrhs, double* a, const int* lda,
		int* ipiv, double* b, const int* ldb, int* info);
/* DGESVD - compute the singular value decomposition (SVD); of a */
/* real M-by-N matrix A, optionally computing the left and/or */
/* right singular vectors */


/** \private */
extern void
FORTRAN_NAME(dgesvd)(const char* jobu, const char* jobvt, const int* m,
		 const int* n, double* a, const int* lda, double* s,
		 double* u, const int* ldu, double* vt, const int* ldvt,
		 double* work, const int* lwork, int* info);
/* DGESVX - use the LU factorization to compute the solution to a */
/* real system of linear equations  A * X = B, */


/** \private */
extern void
FORTRAN_NAME(dgesvx)(const int* fact, const char* trans, const int* n,
		 const int* nrhs, double* a, const int* lda,
		 double* af, const int* ldaf, int* ipiv,
		 char *equed, double* r, double* c,
		 double* b, const int* ldb,
		 double* x, const int* ldx,
		 double* rcond, double* ferr, double* berr,
		 double* work, int* iwork, int* info);
/* DGETF2 - compute an LU factorization of a general m-by-n */
/* matrix A using partial pivoting with row interchanges */


/** \private */
extern void
FORTRAN_NAME(dgetf2)(const int* m, const int* n, double* a, const int* lda,
		 int* ipiv, int* info);
/* DGETRF - compute an LU factorization of a general M-by-N */
/* matrix A using partial pivoting with row interchanges */


/** \private */
extern void
FORTRAN_NAME(dgetrf)(const int* m, const int* n, double* a, const int* lda,
		 int* ipiv, int* info);
/* DGETRI - compute the inverse of a matrix using the LU */
/* factorization computed by DGETRF */


/** \private */
extern void
FORTRAN_NAME(dgetri)(const int* n, double* a, const int* lda,
		 int* ipiv, double* work, const int* lwork, int* info);
/* DGETRS - solve a system of linear equations	A * X = B or A' * */
/* X = B with a general N-by-N matrix A using the LU factorization */
/* computed by DGETRF */


/** \private */
extern void
FORTRAN_NAME(dgetrs)(const char* trans, const int* n, const int* nrhs,
		 double const* a, const int* lda, const int* ipiv,
		 double* b, const int* ldb, int* info);

/* Double precision General matrices Generalized problems */

/* DGGBAK - form the right or left eigenvectors of a real */
/* generalized eigenvalue problem A*x = lambda*B*x, by backward */
/* transformation on the computed eigenvectors of the balanced */
/* pair of matrices output by DGGBAL */


/** \private */
extern void
FORTRAN_NAME(dggbak)(const char* job, const char* side,
		 const int* n, const int* ilo, const int* ihi,
		 double* lscale, double* rscale, const int* m,
		 double* v, const int* ldv, int* info);
/* DGGBAL - balance a pair of general real matrices (A,B); */


/** \private */
extern void
FORTRAN_NAME(dggbal)(const char* job, const int* n, double* a, const int* lda,
		 double* b, const int* ldb, int* ilo, int* ihi,
		 double* lscale, double* rscale, double* work, int* info);
/* DGGES - compute for a pair of N-by-N real nonsymmetric */
/* matrices A, B the generalized eigenvalues, the generalized */
/* real Schur form (S,T), optionally, the left and/or right matrices */
/* of Schur vectors (VSL and VSR)*/


/** \private */
extern void
FORTRAN_NAME(dgges)(const char* jobvsl, const char* jobvsr, const char* sort,
		int (*delztg)(double*, double*, double*),
		const int* n, double* a, const int* lda,
		double* b, const int* ldb, double* alphar,
		double* alphai, double const* beta,
		double* vsl, const int* ldvsl,
		double* vsr, const int* ldvsr,
		double* work, const int* lwork, int* bwork, int* info);

/* DGGGLM - solve a general Gauss-Markov linear model (GLM) problem */


/** \private */
extern void
FORTRAN_NAME(dggglm)(const int* n, const int* m, const int* p,
		 double* a, const int* lda, double* b, const int* ldb,
		 double* d, double* x, double* y,
		 double* work, const int* lwork, int* info);
/* DGGHRD - reduce a pair of real matrices (A,B); to generalized */
/* upper Hessenberg form using orthogonal transformations, where A */
/* is a general matrix and B is upper triangular */


/** \private */
extern void
FORTRAN_NAME(dgghrd)(const char* compq, const char* compz, const int* n,
		 const int* ilo, const int* ihi, double* a, const int* lda,
		 double* b, const int* ldb, double* q, const int* ldq,
		 double* z, const int* ldz, int* info);
/* DGGLSE - solve the linear equality-constrained least squares */
/* (LSE) problem */


/** \private */
extern void
FORTRAN_NAME(dgglse)(const int* m, const int* n, const int* p,
		 double* a, const int* lda,
		 double* b, const int* ldb,
		 double* c, double* d, double* x,
		 double* work, const int* lwork, int* info);
/* DGGQRF - compute a generalized QR factorization of an N-by-M */
/* matrix A and an N-by-P matrix B */


/** \private */
extern void
FORTRAN_NAME(dggqrf)(const int* n, const int* m, const int* p,
		 double* a, const int* lda, double* taua,
		 double* b, const int* ldb, double* taub,
		 double* work, const int* lwork, int* info);
/* DGGRQF - compute a generalized RQ factorization of an M-by-N */
/* matrix A and a P-by-N matrix B */


/** \private */
extern void
FORTRAN_NAME(dggrqf)(const int* m, const int* p, const int* n,
		 double* a, const int* lda, double* taua,
		 double* b, const int* ldb, double* taub,
		 double* work, const int* lwork, int* info);
/* DGGSVD - compute the generalized singular value decomposition */
/* (GSVD) of an M-by-N real matrix A and P-by-N real matrix B */


/** \private */
extern void
FORTRAN_NAME(dggsvd)(const char* jobu, const char* jobv, const char* jobq,
		 const int* m, const int* n, const int* p,
		 const int* k, const int* l,
		 double* a, const int* lda,
		 double* b, const int* ldb,
		 double const* alpha, double const* beta,
		 double* u, const int* ldu,
		 double* v, const int* ldv,
		 double* q, const int* ldq,
		 double* work, int* iwork, int* info);

/* Double precision General Tridiagonal matrices */

/* DGTCON - estimate the reciprocal of the condition number of a real */
/* tridiagonal matrix A using the LU factorization as computed by DGTTRF */


/** \private */
extern void
FORTRAN_NAME(dgtcon)(const char* norm, const int* n, double* dl, double* d,
		 double* du, double* du2, int* ipiv, double const* anorm,
		 double* rcond, double* work, int* iwork, int* info);
/* DGTRFS - improve the computed solution to a system of linear equations */
/* when the coefficient matrix is tridiagonal, and provides error bounds */
/* and backward error estimates for the solution */


/** \private */
extern void
FORTRAN_NAME(dgtrfs)(const char* trans, const int* n, const int* nrhs,
		 double* dl, double* d, double* du, double* dlf,
		 double* df, double* duf, double* du2,
		 int* ipiv, double* b, const int* ldb,
		 double* x, const int* ldx,
		 double* ferr, double* berr,
		 double* work, int* iwork, int* info);
/* DGTSV - solve the equation	A*X = B, */


/** \private */
extern void
FORTRAN_NAME(dgtsv)(const int* n, const int* nrhs,
		double* dl, double* d, double* du,
		double* b, const int* ldb, int* info);
/* DGTSVX - use the LU factorization to compute the solution to a */
/* real system of linear equations A * X = B or A**T * X = B, */


/** \private */
extern void
FORTRAN_NAME(dgtsvx)(const int* fact, const char* trans,
		 const int* n, const int* nrhs,
		 double* dl, double* d, double* du,
		 double* dlf, double* df, double* duf,
		 double* du2, int* ipiv,
		 double* b, const int* ldb,
		 double* x, const int* ldx,
		 double* rcond, double* ferr, double* berr,
		 double* work, int* iwork, int* info);
/* DGTTRF - compute an LU factorization of a real tridiagonal matrix */
/* A using elimination with partial pivoting and row interchanges */


/** \private */
extern void
FORTRAN_NAME(dgttrf)(const int* n, double* dl, double* d,
		 double* du, double* du2, int* ipiv, int* info);
/* DGTTRS - solve one of the systems of equations  A*X = B or */
/* A'*X = B, */


/** \private */
extern void
FORTRAN_NAME(dgttrs)(const char* trans, const int* n, const int* nrhs,
		 double* dl, double* d, double* du, double* du2,
		 int* ipiv, double* b, const int* ldb, int* info);

/* Double precision Orthogonal matrices */

/* DOPGTR - generate a real orthogonal matrix Q which is defined */
/* as the product of n-1 elementary reflectors H(i); of order n, */
/* as returned by DSPTRD using packed storage */


/** \private */
extern void
FORTRAN_NAME(dopgtr)(const char* uplo, const int* n,
		 double const* ap, double const* tau,
		 double* q, const int* ldq,
		 double* work, int* info);
/* DOPMTR - overwrite the general real M-by-N matrix C with */
/* SIDE = 'L' SIDE = 'R' TRANS = 'N' */


/** \private */
extern void
FORTRAN_NAME(dopmtr)(const char* side, const char* uplo,
		 const char* trans, const int* m, const int* n,
		 double const* ap, double const* tau,
		 double* c, const int* ldc,
		 double* work, int* info);
/* DORG2L - generate an m by n real matrix Q with orthonormal */
/* columns, */


/** \private */
extern void
FORTRAN_NAME(dorg2l)(const int* m, const int* n, const int* k,
		 double* a, const int* lda,
		 double const* tau, double* work, int* info);
/* DORG2R - generate an m by n real matrix Q with orthonormal */
/* columns, */


/** \private */
extern void
FORTRAN_NAME(dorg2r)(const int* m, const int* n, const int* k,
		 double* a, const int* lda,
		 double const* tau, double* work, int* info);
/* DORGBR - generate one of the real orthogonal matrices Q or */
/* P**T determined by DGEBRD when reducing a real matrix A to */
/* bidiagonal form */


/** \private */
extern void
FORTRAN_NAME(dorgbr)(const char* vect, const int* m,
		 const int* n, const int* k,
		 double* a, const int* lda,
		 double const* tau, double* work,
		 const int* lwork, int* info);
/* DORGHR - generate a real orthogonal matrix Q which is defined */
/* as the product of IHI-ILO elementary reflectors of order N, as */
/* returned by DGEHRD */


/** \private */
extern void
FORTRAN_NAME(dorghr)(const int* n, const int* ilo, const int* ihi,
		 double* a, const int* lda, double const* tau,
		 double* work, const int* lwork, int* info);
/* DORGL2 - generate an m by n real matrix Q with orthonormal */
/* rows, */


/** \private */
extern void
FORTRAN_NAME(dorgl2)(const int* m, const int* n, const int* k,
		 double* a, const int* lda, double const* tau,
		 double* work, int* info);
/* DORGLQ - generate an M-by-N real matrix Q with orthonormal */
/* rows, */


/** \private */
extern void
FORTRAN_NAME(dorglq)(const int* m, const int* n, const int* k,
		 double* a, const int* lda,
		 double const* tau, double* work,
		 const int* lwork, int* info);
/* DORGQL - generate an M-by-N real matrix Q with orthonormal */
/* columns, */


/** \private */
extern void
FORTRAN_NAME(dorgql)(const int* m, const int* n, const int* k,
		 double* a, const int* lda,
		 double const* tau, double* work,
		 const int* lwork, int* info);
/* DORGQR - generate an M-by-N real matrix Q with orthonormal */
/* columns, */


/** \private */
extern void
FORTRAN_NAME(dorgqr)(const int* m, const int* n, const int* k,
		 double* a, const int* lda, double const* tau,
		 double* work, const int* lwork, int* info);
/* DORGR2 - generate an m by n real matrix Q with orthonormal */
/* rows, */


/** \private */
extern void
FORTRAN_NAME(dorgr2)(const int* m, const int* n, const int* k,
		 double* a, const int* lda, double const* tau,
		 double* work, int* info);
/* DORGRQ - generate an M-by-N real matrix Q with orthonormal rows */


/** \private */
extern void
FORTRAN_NAME(dorgrq)(const int* m, const int* n, const int* k,
		 double* a, const int* lda, double const* tau,
		 double* work, const int* lwork, int* info);
/* DORGTR - generate a real orthogonal matrix Q which is defined */
/* as the product of n-1 elementary reflectors of order const int* n, as */
/* returned by DSYTRD */


/** \private */
extern void
FORTRAN_NAME(dorgtr)(const char* uplo, const int* n,
		 double* a, const int* lda, double const* tau,
		 double* work, const int* lwork, int* info);
/* DORM2L - overwrite the general real m by n matrix C with   Q * */
/* C if SIDE = 'L' and TRANS = 'N', or	 Q'* C if SIDE = 'L' and */
/* TRANS = 'T', or   C * Q if SIDE = 'R' and TRANS = 'N', or   C * */
/* Q' if SIDE = 'R' and TRANS = 'T', */


/** \private */
extern void
FORTRAN_NAME(dorm2l)(const char* side, const char* trans,
		 const int* m, const int* n, const int* k,
		 double const* a, const int* lda,
		 double const* tau, double* c, const int* ldc,
		 double* work, int* info);
/* DORM2R - overwrite the general real m by n matrix C with   Q * C */
/* if SIDE = 'L' and TRANS = 'N', or   Q'* C if SIDE = 'L' and */
/* TRANS = 'T', or   C * Q if SIDE = 'R' and TRANS = 'N', or   C * */
/* Q' if SIDE = 'R' and TRANS = 'T', */


/** \private */
extern void
FORTRAN_NAME(dorm2r)(const char* side, const char* trans,
		 const int* m, const int* n, const int* k,
		 double const* a, const int* lda, double const* tau,
		 double* c, const int* ldc, double* work, int* info);
/* DORMBR - VECT = 'Q', DORMBR overwrites the general real M-by-N */
/* matrix C with  SIDE = 'L' SIDE = 'R' TRANS = 'N' */


/** \private */
extern void
FORTRAN_NAME(dormbr)(const char* vect, const char* side, const char* trans,
		 const int* m, const int* n, const int* k,
		 double const* a, const int* lda, double const* tau,
		 double* c, const int* ldc,
		 double* work, const int* lwork, int* info);
/* DORMHR - overwrite the general real M-by-N matrix C with */
/* SIDE = 'L' SIDE = 'R' TRANS = 'N' */


/** \private */
extern void
FORTRAN_NAME(dormhr)(const char* side, const char* trans, const int* m,
		 const int* n, const int* ilo, const int* ihi,
		 double const* a, const int* lda, double const* tau,
		 double* c, const int* ldc,
		 double* work, const int* lwork, int* info);
/* DORML2 - overwrite the general real m by n matrix C with   Q * */
/* C if SIDE = 'L' and TRANS = 'N', or	 Q'* C if SIDE = 'L' and */
/* TRANS = 'T', or   C * Q if SIDE = 'R' and TRANS = 'N', or   C * */
/* Q' if SIDE = 'R' and TRANS = 'T', */


/** \private */
extern void
FORTRAN_NAME(dorml2)(const char* side, const char* trans,
		 const int* m, const int* n, const int* k,
		 double const* a, const int* lda, double const* tau,
		 double* c, const int* ldc, double* work, int* info);
/* DORMLQ - overwrite the general real M-by-N matrix C with */
/* SIDE = 'L' SIDE = 'R' TRANS = 'N'  */


/** \private */
extern void
FORTRAN_NAME(dormlq)(const char* side, const char* trans,
		 const int* m, const int* n, const int* k,
		 double const* a, const int* lda,
		 double const* tau, double* c, const int* ldc,
		 double* work, const int* lwork, int* info);
/* DORMQL - overwrite the general real M-by-N matrix C with */
/* SIDE = 'L' SIDE = 'R' TRANS = 'N' */


/** \private */
extern void
FORTRAN_NAME(dormql)(const char* side, const char* trans,
		 const int* m, const int* n, const int* k,
		 double const* a, const int* lda,
		 double const* tau, double* c, const int* ldc,
		 double* work, const int* lwork, int* info);
/* DORMQR - overwrite the general real M-by-N matrix C with   SIDE = */
/* 'L' SIDE = 'R' TRANS = 'N' */


/** \private */
extern void
FORTRAN_NAME(dormqr)(const char* side, const char* trans,
		 const int* m, const int* n, const int* k,
		 double const* a, const int* lda,
		 double const* tau, double* c, const int* ldc,
		 double* work, const int* lwork, int* info);
/* DORMR2 - overwrite the general real m by n matrix C with   Q * */
/* C if SIDE = 'L' and TRANS = 'N', or	 Q'* C if SIDE = 'L' and */
/* TRANS = 'T', or   C * Q if SIDE = 'R' and TRANS = 'N', or   C * */
/* Q' if SIDE = 'R' and TRANS = 'T', */


/** \private */
extern void
FORTRAN_NAME(dormr2)(const char* side, const char* trans,
		 const int* m, const int* n, const int* k,
		 double const* a, const int* lda,
		 double const* tau, double* c, const int* ldc,
		 double* work, int* info);
/* DORMRQ - overwrite the general real M-by-N matrix C with */
/* SIDE = 'L' SIDE = 'R' TRANS = 'N' */


/** \private */
extern void
FORTRAN_NAME(dormrq)(const char* side, const char* trans,
		 const int* m, const int* n, const int* k,
		 double const* a, const int* lda,
		 double const* tau, double* c, const int* ldc,
		 double* work, const int* lwork, int* info);
/* DORMTR - overwrite the general real M-by-N matrix C with */
/* SIDE = 'L' SIDE = 'R' TRANS = 'N' */


/** \private */
extern void
FORTRAN_NAME(dormtr)(const char* side, const char* uplo,
		 const char* trans, const int* m, const int* n,
		 double const* a, const int* lda,
		 double const* tau, double* c, const int* ldc,
		 double* work, const int* lwork, int* info);

/* Double precision Positive definite Band matrices */

/* DPBCON - estimate the reciprocal of the condition number (in */
/* the 1-norm); of a real symmetric positive definite band matrix */
/* using the Cholesky factorization A = U**T*U or A = L*L**T */
/* computed by DPBTRF */


/** \private */
extern void
FORTRAN_NAME(dpbcon)(const char* uplo, const int* n, const int* kd,
		 double const* ab, const int* ldab,
		 double const* anorm, double* rcond,
		 double* work, int* iwork, int* info);
/* DPBEQU - compute row and column scalings intended to */
/* equilibrate a symmetric positive definite band matrix A and */
/* reduce its condition number (with respect to the two-norm); */


/** \private */
extern void
FORTRAN_NAME(dpbequ)(const char* uplo, const int* n, const int* kd,
		 double const* ab, const int* ldab,
		 double* s, double* scond, double* amax, int* info);
/* DPBRFS - improve the computed solution to a system of linear */
/* equations when the coefficient matrix is symmetric positive */
/* definite and banded, and provides error bounds and backward */
/* error estimates for the solution */


/** \private */
extern void
FORTRAN_NAME(dpbrfs)(const char* uplo, const int* n,
		 const int* kd, const int* nrhs,
		 double const* ab, const int* ldab,
		 double const* afb, const int* ldafb,
		 double const* b, const int* ldb,
		 double* x, const int* ldx,
		 double* ferr, double* berr,
		 double* work, int* iwork, int* info);
/* DPBSTF - compute a split Cholesky factorization of a real */
/* symmetric positive definite band matrix A */


/** \private */
extern void
FORTRAN_NAME(dpbstf)(const char* uplo, const int* n, const int* kd,
		 double* ab, const int* ldab, int* info);
/* DPBSV - compute the solution to a real system of linear */
/* equations  A * X = B, */


/** \private */
extern void
FORTRAN_NAME(dpbsv)(const char* uplo, const int* n,
		const int* kd, const int* nrhs,
		double* ab, const int* ldab,
		double* b, const int* ldb, int* info);
/* DPBSVX - use the Cholesky factorization A = U**T*U or A = */
/* L*L**T to compute the solution to a real system of linear */
/* equations  A * X = B, */


/** \private */
extern void
FORTRAN_NAME(dpbsvx)(const int* fact, const char* uplo, const int* n,
		 const int* kd, const int* nrhs,
		 double* ab, const int* ldab,
		 double* afb, const int* ldafb,
		 char* equed, double* s,
		 double* b, const int* ldb,
		 double* x, const int* ldx, double* rcond,
		 double* ferr, double* berr,
		 double* work, int* iwork, int* info);
/* DPBTF2 - compute the Cholesky factorization of a real */
/* symmetric positive definite band matrix A */


/** \private */
extern void
FORTRAN_NAME(dpbtf2)(const char* uplo, const int* n, const int* kd,
		 double* ab, const int* ldab, int* info);
/* DPBTRF - compute the Cholesky factorization of a real */
/* symmetric positive definite band matrix A */


/** \private */
extern void
FORTRAN_NAME(dpbtrf)(const char* uplo, const int* n, const int* kd,
		 double* ab, const int* ldab, int* info);
/* DPBTRS - solve a system of linear equations A*X = B with a */
/* symmetric positive definite band matrix A using the Cholesky */
/* factorization A = U**T*U or A = L*L**T computed by DPBTRF */


/** \private */
extern void
FORTRAN_NAME(dpbtrs)(const char* uplo, const int* n,
		 const int* kd, const int* nrhs,
		 double const* ab, const int* ldab,
		 double* b, const int* ldb, int* info);

/* Double precision Positive definite matrices */

/* DPOCON - estimate the reciprocal of the condition number (in */
/* the 1-norm); of a real symmetric positive definite matrix using */
/* the Cholesky factorization A = U**T*U or A = L*L**T computed by */
/* DPOTRF */


/** \private */
extern void
FORTRAN_NAME(dpocon)(const char* uplo, const int* n,
		 double const* a, const int* lda,
		 double const* anorm, double* rcond,
		 double* work, int* iwork, int* info);
/* DPOEQU - compute row and column scalings intended to */
/* equilibrate a symmetric positive definite matrix A and reduce */
/* its condition number (with respect to the two-norm); */


/** \private */
extern void
FORTRAN_NAME(dpoequ)(const int* n, double const* a, const int* lda,
		 double* s, double* scond, double* amax, int* info);
/* DPORFS - improve the computed solution to a system of linear */
/* equations when the coefficient matrix is symmetric positive */
/* definite, */


/** \private */
extern void
FORTRAN_NAME(dporfs)(const char* uplo, const int* n, const int* nrhs,
		 double const* a, const int* lda,
		 double const* af, const int* ldaf,
		 double const* b, const int* ldb,
		 double* x, const int* ldx,
		 double* ferr, double* berr,
		 double* work, int* iwork, int* info);
/* DPOSV - compute the solution to a real system of linear */
/* equations  A * X = B, */


/** \private */
extern void
FORTRAN_NAME(dposv)(const char* uplo, const int* n, const int* nrhs,
		double* a, const int* lda,
		double* b, const int* ldb, int* info);
/* DPOSVX - use the Cholesky factorization A = U**T*U or A = */
/* L*L**T to compute the solution to a real system of linear */
/* equations  A * X = B, */


/** \private */
extern void
FORTRAN_NAME(dposvx)(const int* fact, const char* uplo,
		 const int* n, const int* nrhs,
		 double* a, const int* lda,
		 double* af, const int* ldaf, char* equed,
		 double* s, double* b, const int* ldb,
		 double* x, const int* ldx, double* rcond,
		 double* ferr, double* berr, double* work,
		 int* iwork, int* info);
/* DPOTF2 - compute the Cholesky factorization of a real */
/* symmetric positive definite matrix A */


/** \private */
extern void
FORTRAN_NAME(dpotf2)(const char* uplo, const int* n,
		 double* a, const int* lda, int* info);
/* DPOTRF - compute the Cholesky factorization of a real */
/* symmetric positive definite matrix A */


/** \private */
extern void
FORTRAN_NAME(dpotrf)(const char* uplo, const int* n,
		 double* a, const int* lda, int* info);
/* DPOTRI - compute the inverse of a real symmetric positive */
/* definite matrix A using the Cholesky factorization A = U**T*U */
/* or A = L*L**T computed by DPOTRF */


/** \private */
extern void
FORTRAN_NAME(dpotri)(const char* uplo, const int* n,
		 double* a, const int* lda, int* info);
/* DPOTRS - solve a system of linear equations A*X = B with a */
/* symmetric positive definite matrix A using the Cholesky */
/* factorization A = U**T*U or A = L*L**T computed by DPOTRF */


/** \private */
extern void
FORTRAN_NAME(dpotrs)(const char* uplo, const int* n,
		 const int* nrhs, double const* a, const int* lda,
		 double* b, const int* ldb, int* info);
/* DPPCON - estimate the reciprocal of the condition number (in */
/* the 1-norm); of a real symmetric positive definite packed */
/* matrix using the Cholesky factorization A = U**T*U or A = */
/* L*L**T computed by DPPTRF */


/** \private */
extern void
FORTRAN_NAME(dppcon)(const char* uplo, const int* n,
		 double const* ap, double const* anorm, double* rcond,
		 double* work, int* iwork, int* info);
/* DPPEQU - compute row and column scalings intended to */
/* equilibrate a symmetric positive definite matrix A in packed */
/* storage and reduce its condition number (with respect to the */
/* two-norm); */


/** \private */
extern void
FORTRAN_NAME(dppequ)(const char* uplo, const int* n,
		 double const* ap, double* s, double* scond,
		 double* amax, int* info);

/* Double precision Positive definite matrices in Packed storage */

/* DPPRFS - improve the computed solution to a system of linear */
/* equations when the coefficient matrix is symmetric positive */
/* definite and packed, and provides error bounds and backward */
/* error estimates for the solution */


/** \private */
extern void
FORTRAN_NAME(dpprfs)(const char* uplo, const int* n, const int* nrhs,
		 double const* ap, double const* afp,
		 double const* b, const int* ldb,
		 double* x, const int* ldx,
		 double* ferr, double* berr,
		 double* work, int* iwork, int* info);
/* DPPSV - compute the solution to a real system of linear */
/* equations  A * X = B, */


/** \private */
extern void
FORTRAN_NAME(dppsv)(const char* uplo, const int* n,
		const int* nrhs, double const* ap,
		double* b, const int* ldb, int* info);
/* DPPSVX - use the Cholesky factorization A = U**T*U or A = */
/* L*L**T to compute the solution to a real system of linear */
/* equations  A * X = B, */


/** \private */
extern void
FORTRAN_NAME(dppsvx)(const int* fact, const char* uplo,
		 const int* n, const int* nrhs, double* ap,
		 double* afp, char* equed, double* s,
		 double* b, const int* ldb,
		 double* x, const int* ldx,
		 double* rcond, double* ferr, double* berr,
		 double* work, int* iwork, int* info);
/* DPPTRF - compute the Cholesky factorization of a real */
/* symmetric positive definite matrix A stored in packed format */


/** \private */
extern void
FORTRAN_NAME(dpptrf)(const char* uplo, const int* n, double* ap, int* info);
/* DPPTRI - compute the inverse of a real symmetric positive */
/* definite matrix A using the Cholesky factorization A = U**T*U */
/* or A = L*L**T computed by DPPTRF  */


/** \private */
extern void
FORTRAN_NAME(dpptri)(const char* uplo, const int* n, double* ap, int* info);
/* DPPTRS - solve a system of linear equations A*X = B with a */
/* symmetric positive definite matrix A in packed storage using */
/* the Cholesky factorization A = U**T*U or A = L*L**T computed by */
/* DPPTRF */


/** \private */
extern void
FORTRAN_NAME(dpptrs)(const char* uplo, const int* n,
		 const int* nrhs, double const* ap,
		 double* b, const int* ldb, int* info);

/* Double precision symmetric Positive definite Tridiagonal matrices */

/* DPTCON - compute the reciprocal of the condition number (in */
/* the 1-norm); of a real symmetric positive definite tridiagonal */
/* matrix using the factorization A = L*D*L**T or A = U**T*D*U */
/* computed by DPTTRF */


/** \private */
extern void
FORTRAN_NAME(dptcon)(const int* n,
		 double const* d, double const* e,
		 double const* anorm, double* rcond,
		 double* work, int* info);
/* DPTEQR - compute all eigenvalues and, optionally, eigenvectors */
/* of a symmetric positive definite tridiagonal matrix by first */
/* factoring the matrix using DPTTRF, and then calling DBDSQR to */
/* compute the singular values of the bidiagonal factor */


/** \private */
extern void
FORTRAN_NAME(dpteqr)(const char* compz, const int* n, double* d,
		 double* e, double* z, const int* ldz,
		 double* work, int* info);
/* DPTRFS - improve the computed solution to a system of linear */
/* equations when the coefficient matrix is symmetric positive */
/* definite and tridiagonal, and provides error bounds and */
/* backward error estimates for the solution */


/** \private */
extern void
FORTRAN_NAME(dptrfs)(const int* n, const int* nrhs,
		 double const* d, double const* e,
		 double const* df, double const* ef,
		 double const* b, const int* ldb,
		 double* x, const int* ldx,
		 double* ferr, double* berr,
		 double* work, int* info);
/* DPTSV - compute the solution to a real system of linear */
/* equations A*X = B, where A is an N-by-N symmetric positive */
/* definite tridiagonal matrix, and X and B are N-by-NRHS matrices */


/** \private */
extern void
FORTRAN_NAME(dptsv)(const int* n, const int* nrhs, double* d,
		double* e, double* b, const int* ldb, int* info);
/* DPTSVX - use the factorization A = L*D*L**T to compute the */
/* solution to a real system of linear equations A*X = B, where A */
/* is an N-by-N symmetric positive definite tridiagonal matrix and */
/* X and B are N-by-NRHS matrices */


/** \private */
extern void
FORTRAN_NAME(dptsvx)(const int* fact, const int* n,
		 const int* nrhs,
		 double const* d, double const* e,
		 double* df, double* ef,
		 double const* b, const int* ldb,
		 double* x, const int* ldx, double* rcond,
		 double* ferr, double* berr,
		 double* work, int* info);
/* DPTTRF - compute the factorization of a real symmetric */
/* positive definite tridiagonal matrix A */


/** \private */
extern void
FORTRAN_NAME(dpttrf)(const int* n, double* d, double* e, int* info);
/* DPTTRS - solve a system of linear equations A * X = B with a */
/* symmetric positive definite tridiagonal matrix A using the */
/* factorization A = L*D*L**T or A = U**T*D*U computed by DPTTRF */


/** \private */
extern void
FORTRAN_NAME(dpttrs)(const int* n, const int* nrhs,
		 double const* d, double const* e,
		 double* b, const int* ldb, int* info);
/* DRSCL - multiply an n-element real vector x by the real scalar */
/* 1/a */


/** \private */
extern void
FORTRAN_NAME(drscl)(const int* n, double const* da,
		double* x, const int* incx);

/* Double precision Symmetric Band matrices */

/* DSBEV - compute all the eigenvalues and, optionally, */
/* eigenvectors of a real symmetric band matrix A */


/** \private */
extern void
FORTRAN_NAME(dsbev)(const char* jobz, const char* uplo,
		const int* n, const int* kd,
		double* ab, const int* ldab,
		double* w, double* z, const int* ldz,
		double* work, int* info);
/* DSBEVD - compute all the eigenvalues and, optionally, */
/* eigenvectors of a real symmetric band matrix A */


/** \private */
extern void
FORTRAN_NAME(dsbevd)(const char* jobz, const char* uplo,
		 const int* n, const int* kd,
		 double* ab, const int* ldab,
		 double* w, double* z, const int* ldz,
		 double* work, const int* lwork,
		 int* iwork, const int* liwork, int* info);
/* DSBEVX - compute selected eigenvalues and, optionally, */
/* eigenvectors of a real symmetric band matrix A */


/** \private */
extern void
FORTRAN_NAME(dsbevx)(const char* jobz, const char* range,
		 const char* uplo, const int* n, const int* kd,
		 double* ab, const int* ldab,
		 double* q, const int* ldq,
		 double const* vl, double const* vu,
		 const int* il, const int* iu,
		 double const* abstol,
		 int* m, double* w,
		 double* z, const int* ldz,
		 double* work, int* iwork,
		 int* ifail, int* info);
/* DSBGST - reduce a real symmetric-definite banded generalized */
/* eigenproblem A*x = lambda*B*x to standard form C*y = lambda*y, */


/** \private */
extern void
FORTRAN_NAME(dsbgst)(const char* vect, const char* uplo,
		 const int* n, const int* ka, const int* kb,
		 double* ab, const int* ldab,
		 double* bb, const int* ldbb,
		 double* x, const int* ldx,
		 double* work, int* info);
/* DSBGV - compute all the eigenvalues, and optionally, the */
/* eigenvectors of a real generalized symmetric-definite banded */
/* eigenproblem, of the form A*x=(lambda);*B*x */


/** \private */
extern void
FORTRAN_NAME(dsbgv)(const char* jobz, const char* uplo,
		const int* n, const int* ka, const int* kb,
		double* ab, const int* ldab,
		double* bb, const int* ldbb,
		double* w, double* z, const int* ldz,
		double* work, int* info);
/* DSBTRD - reduce a real symmetric band matrix A to symmetric */
/* tridiagonal form T by an orthogonal similarity transformation */


/** \private */
extern void
FORTRAN_NAME(dsbtrd)(const char* vect, const char* uplo,
		 const int* n, const int* kd,
		 double* ab, const int* ldab,
		 double* d, double* e,
		 double* q, const int* ldq,
		 double* work, int* info);

/* Double precision Symmetric Packed matrices */

/* DSPCON - estimate the reciprocal of the condition number (in */
/* the 1-norm); of a real symmetric packed matrix A using the */
/* factorization A = U*D*U**T or A = L*D*L**T computed by DSPTRF */


/** \private */
extern void
FORTRAN_NAME(dspcon)(const char* uplo, const int* n,
		 double const* ap, const int* ipiv,
		 double const* anorm, double* rcond,
		 double* work, int* iwork, int* info);
/* DSPEV - compute all the eigenvalues and, optionally, */
/* eigenvectors of a real symmetric matrix A in packed storage */


/** \private */
extern void
FORTRAN_NAME(dspev)(const char* jobz, const char* uplo, const int* n,
		double* ap, double* w, double* z, const int* ldz,
		double* work, int* info);
/* DSPEVD - compute all the eigenvalues and, optionally, */
/* eigenvectors of a real symmetric matrix A in packed storage */


/** \private */
extern void
FORTRAN_NAME(dspevd)(const char* jobz, const char* uplo,
		 const int* n, double* ap, double* w,
		 double* z, const int* ldz,
		 double* work, const int* lwork,
		 int* iwork, const int* liwork, int* info);
/* DSPEVX - compute selected eigenvalues and, optionally, */
/* eigenvectors of a real symmetric matrix A in packed storage */


/** \private */
extern void
FORTRAN_NAME(dspevx)(const char* jobz, const char* range,
		 const char* uplo, const int* n, double* ap,
		 double const* vl, double const* vu,
		 const int* il, const int* iu,
		 double const* abstol,
		 int* m, double* w,
		 double* z, const int* ldz,
		 double* work, int* iwork,
		 int* ifail, int* info);
/* DSPGST - reduce a real symmetric-definite generalized */
/* eigenproblem to standard form, using packed storage */


/** \private */
extern void
FORTRAN_NAME(dspgst)(const int* itype, const char* uplo,
		 const int* n, double* ap, double* bp, int* info);
/* DSPGV - compute all the eigenvalues and, optionally, the */
/* eigenvectors of a real generalized symmetric-definite */
/* eigenproblem, of the form A*x=(lambda)*B*x, A*Bx=(lambda)*x, */
/* or B*A*x=(lambda)*x */


/** \private */
extern void
FORTRAN_NAME(dspgv)(const int* itype, const char* jobz,
		const char* uplo, const int* n,
		double* ap, double* bp, double* w,
		double* z, const int* ldz,
		double* work, int* info);

/* DSPRFS - improve the computed solution to a system of linear */
/* equations when the coefficient matrix is symmetric indefinite */
/* and packed, and provides error bounds and backward error */
/* estimates for the solution */


/** \private */
extern void
FORTRAN_NAME(dsprfs)(const char* uplo, const int* n,
		 const int* nrhs, double const* ap,
		 double const* afp, const int* ipiv,
		 double const* b, const int* ldb,
		 double* x, const int* ldx,
		 double* ferr, double* berr,
		 double* work, int* iwork, int* info);

/* DSPSV - compute the solution to a real system of linear */
/* equations  A * X = B, */


/** \private */
extern void
FORTRAN_NAME(dspsv)(const char* uplo, const int* n,
		const int* nrhs, double* ap, int* ipiv,
		double* b, const int* ldb, int* info);

/* DSPSVX - use the diagonal pivoting factorization A = U*D*U**T */
/* or A = L*D*L**T to compute the solution to a real system of */
/* linear equations A * X = B, where A is an N-by-N symmetric */
/* matrix stored in packed format and X and B are N-by-NRHS */
/* matrices */


/** \private */
extern void
FORTRAN_NAME(dspsvx)(const int* fact, const char* uplo,
		 const int* n, const int* nrhs,
		 double const* ap, double* afp, int* ipiv,
		 double const* b, const int* ldb,
		 double* x, const int* ldx,
		 double* rcond, double* ferr, double* berr,
		 double* work, int* iwork, int* info);

/* DSPTRD - reduce a real symmetric matrix A stored in packed */
/* form to symmetric tridiagonal form T by an orthogonal */
/* similarity transformation */


/** \private */
extern void
FORTRAN_NAME(dsptrd)(const char* uplo, const int* n,
		 double* ap, double* d, double* e,
		 double* tau, int* info);

/* DSPTRF - compute the factorization of a real symmetric matrix */
/* A stored in packed format using the Bunch-Kaufman diagonal */
/* pivoting method */


/** \private */
extern void
FORTRAN_NAME(dsptrf)(const char* uplo, const int* n,
		 double* ap, int* ipiv, int* info);

/* DSPTRI - compute the inverse of a real symmetric indefinite */
/* matrix A in packed storage using the factorization A = U*D*U**T */
/* or A = L*D*L**T computed by DSPTRF */


/** \private */
extern void
FORTRAN_NAME(dsptri)(const char* uplo, const int* n,
		 double* ap, const int* ipiv,
		 double* work, int* info);

/* DSPTRS - solve a system of linear equations A*X = B with a */
/* real symmetric matrix A stored in packed format using the */
/* factorization A = U*D*U**T or A = L*D*L**T computed by DSPTRF */


/** \private */
extern void
FORTRAN_NAME(dsptrs)(const char* uplo, const int* n,
		 const int* nrhs, double const* ap,
		 const int* ipiv, double* b, const int* ldb, int* info);

/* Double precision Symmetric Tridiagonal matrices */

/* DSTEBZ - compute the eigenvalues of a symmetric tridiagonal */
/* matrix T */


/** \private */
extern void
FORTRAN_NAME(dstebz)(const char* range, const char* order, const int* n,
		 double const* vl, double const* vu,
		 const int* il, const int* iu,
		 double const *abstol,
		 double const* d, double const* e,
		 int* m, int* nsplit, double* w,
		 int* iblock, int* isplit,
		 double* work, int* iwork,
		 int* info);
/* DSTEDC - compute all eigenvalues and, optionally, eigenvectors */
/* of a symmetric tridiagonal matrix using the divide and conquer */
/* method */


/** \private */
extern void
FORTRAN_NAME(dstedc)(const char* compz, const int* n,
		 double* d, double* e,
		 double* z, const int* ldz,
		 double* work, const int* lwork,
		 int* iwork, const int* liwork, int* info);
/* DSTEIN - compute the eigenvectors of a real symmetric */
/* tridiagonal matrix T corresponding to specified eigenvalues, */
/* using inverse iteration */


/** \private */
extern void
FORTRAN_NAME(dstein)(const int* n, double const* d, double const* e,
		 const int* m, double const* w,
		 const int* iblock, const int* isplit,
		 double* z, const int* ldz,
		 double* work, int* iwork,
		 int* ifail, int* info);
/* DSTEQR - compute all eigenvalues and, optionally, eigenvectors */
/* of a symmetric tridiagonal matrix using the implicit QL or QR */
/* method */


/** \private */
extern void
FORTRAN_NAME(dsteqr)(const char* compz, const int* n, double* d, double* e,
		 double* z, const int* ldz, double* work, int* info);
/* DSTERF - compute all eigenvalues of a symmetric tridiagonal */
/* matrix using the Pal-Walker-Kahan variant of the QL or QR */
/* algorithm */


/** \private */
extern void
FORTRAN_NAME(dsterf)(const int* n, double* d, double* e, int* info);
/* DSTEV - compute all eigenvalues and, optionally, eigenvectors */
/* of a real symmetric tridiagonal matrix A */


/** \private */
extern void
FORTRAN_NAME(dstev)(const char* jobz, const int* n,
		double* d, double* e,
		double* z, const int* ldz,
		double* work, int* info);
/* DSTEVD - compute all eigenvalues and, optionally, eigenvectors */
/* of a real symmetric tridiagonal matrix */


/** \private */
extern void
FORTRAN_NAME(dstevd)(const char* jobz, const int* n,
		 double* d, double* e,
		 double* z, const int* ldz,
		 double* work, const int* lwork,
		 int* iwork, const int* liwork, int* info);
/* DSTEVX - compute selected eigenvalues and, optionally, */
/* eigenvectors of a real symmetric tridiagonal matrix A */


/** \private */
extern void
FORTRAN_NAME(dstevx)(const char* jobz, const char* range,
		 const int* n, double* d, double* e,
		 double const* vl, double const* vu,
		 const int* il, const int* iu,
		 double const* abstol,
		 int* m, double* w,
		 double* z, const int* ldz,
		 double* work, int* iwork,
		 int* ifail, int* info);

/* Double precision SYmmetric matrices */

/* DSYCON - estimate the reciprocal of the condition number (in */
/* the 1-norm); of a real symmetric matrix A using the */
/* factorization A = U*D*U**T or A = L*D*L**T computed by DSYTRF */


/** \private */
extern void
FORTRAN_NAME(dsycon)(const char* uplo, const int* n,
		 double const* a, const int* lda,
		 const int* ipiv,
		 double const* anorm, double* rcond,
		 double* work, int* iwork, int* info);
/* DSYEV - compute all eigenvalues and, optionally, eigenvectors */
/* of a real symmetric matrix A */


/** \private */
extern void
FORTRAN_NAME(dsyev)(const char* jobz, const char* uplo,
		const int* n, double* a, const int* lda,
		double* w, double* work, const int* lwork, int* info);
/* DSYEVD - compute all eigenvalues and, optionally, eigenvectors */
/* of a real symmetric matrix A */


/** \private */
extern void
FORTRAN_NAME(dsyevd)(const char* jobz, const char* uplo,
		 const int* n, double* a, const int* lda,
		 double* w, double* work, const int* lwork,
		 int* iwork, const int* liwork, int* info);
/* DSYEVX - compute selected eigenvalues and, optionally, */
/* eigenvectors of a real symmetric matrix A */


/** \private */
extern void
FORTRAN_NAME(dsyevx)(const char* jobz, const char* range,
		 const char* uplo, const int* n,
		 double* a, const int* lda,
		 double const* vl, double const* vu,
		 const int* il, const int* iu,
		 double const* abstol,
		 int* m, double* w,
		 double* z, const int* ldz,
		 double* work, const int* lwork, int* iwork,
		 int* ifail, int* info);
/* DSYEVR - compute all eigenvalues and, optionally, eigenvectors   */
/* of a real symmetric matrix A					   */


/** \private */
extern void
FORTRAN_NAME(dsyevr)(const char *jobz, const char *range, const char *uplo,
		 const int *n, double *a, const int *lda,
		 double const *vl, double const *vu,
		 const int *il, const int *iu,
		 double const *abstol, int *m, double *w,
		 double *z, const int *ldz, int *isuppz,
		 double *work, const int *lwork,
		 int *iwork, const int *liwork,
		 int *info);
/* DSYGS2 - reduce a real symmetric-definite generalized */
/* eigenproblem to standard form */


/** \private */
extern void
FORTRAN_NAME(dsygs2)(const int* itype, const char* uplo,
		 const int* n, double* a, const int* lda,
		 double const* b, const int* ldb, int* info);
/* DSYGST - reduce a real symmetric-definite generalized */
/* eigenproblem to standard form */


/** \private */
extern void
FORTRAN_NAME(dsygst)(const int* itype, const char* uplo,
		 const int* n, double* a, const int* lda,
		 double const* b, const int* ldb, int* info);
/* DSYGV - compute all the eigenvalues, and optionally, the */
/* eigenvectors of a real generalized symmetric-definite */
/* eigenproblem, of the form A*x=(lambda);*B*x, A*Bx=(lambda);*x, */
/* or B*A*x=(lambda);*x */


/** \private */
extern void
FORTRAN_NAME(dsygv)(const int* itype, const char* jobz,
		const char* uplo, const int* n,
		double* a, const int* lda,
		double* b, const int* ldb,
		double* w, double* work, const int* lwork,
		int* info);
/* DSYRFS - improve the computed solution to a system of linear */
/* equations when the coefficient matrix is symmetric indefinite, */
/* and provides error bounds and backward error estimates for the */
/* solution */


/** \private */
extern void
FORTRAN_NAME(dsyrfs)(const char* uplo, const int* n,
		 const int* nrhs,
		 double const* a, const int* lda,
		 double const* af, const int* ldaf,
		 const int* ipiv,
		 double const* b, const int* ldb,
		 double* x, const int* ldx,
		 double* ferr, double* berr,
		 double* work, int* iwork, int* info);

/* DSYSV - compute the solution to a real system of linear */
/* equations  A * X = B, */


/** \private */
extern void
FORTRAN_NAME(dsysv)(const char* uplo, const int* n,
		const int* nrhs,
		double* a, const int* lda, int* ipiv,
		double* b, const int* ldb,
		double* work, const int* lwork, int* info);

/* DSYSVX - use the diagonal pivoting factorization to compute */
/* the solution to a real system of linear equations A * X = B, */


/** \private */
extern void
FORTRAN_NAME(dsysvx)(const int* fact, const char* uplo,
		 const int* n, const int* nrhs,
		 double const* a, const int* lda,
		 double* af, const int* ldaf, int* ipiv,
		 double const* b, const int* ldb,
		 double* x, const int* ldx, double* rcond,
		 double* ferr, double* berr,
		 double* work, const int* lwork,
		 int* iwork, int* info);

/* DSYTD2 - reduce a real symmetric matrix A to symmetric */
/* tridiagonal form T by an orthogonal similarity transformation */


/** \private */
extern void
FORTRAN_NAME(dsytd2)(const char* uplo, const int* n,
		 double* a, const int* lda,
		 double* d, double* e, double* tau,
		 int* info);

/* DSYTF2 - compute the factorization of a real symmetric matrix */
/* A using the Bunch-Kaufman diagonal pivoting method */


/** \private */
extern void
FORTRAN_NAME(dsytf2)(const char* uplo, const int* n,
		 double* a, const int* lda,
		 int* ipiv, int* info);

/* DSYTRD - reduce a real symmetric matrix A to real symmetric */
/* tridiagonal form T by an orthogonal similarity transformation */


/** \private */
extern void
FORTRAN_NAME(dsytrd)(const char* uplo, const int* n,
		 double* a, const int* lda,
		 double* d, double* e, double* tau,
		 double* work, const int* lwork, int* info);

/* DSYTRF - compute the factorization of a real symmetric matrix */
/* A using the Bunch-Kaufman diagonal pivoting method */


/** \private */
extern void
FORTRAN_NAME(dsytrf)(const char* uplo, const int* n,
		 double* a, const int* lda, int* ipiv,
		 double* work, const int* lwork, int* info);

/* DSYTRI - compute the inverse of a real symmetric indefinite */
/* matrix A using the factorization A = U*D*U**T or A = L*D*L**T */
/* computed by DSYTRF */


/** \private */
extern void
FORTRAN_NAME(dsytri)(const char* uplo, const int* n,
		 double* a, const int* lda, const int* ipiv,
		 double* work, int* info);

/* DSYTRS - solve a system of linear equations A*X = B with a */
/* real symmetric matrix A using the factorization A = U*D*U**T or */
/* A = L*D*L**T computed by DSYTRF */


/** \private */
extern void
FORTRAN_NAME(dsytrs)(const char* uplo, const int* n,
		 const int* nrhs,
		 double const* a, const int* lda,
		 const int* ipiv,
		 double* b, const int* ldb, int* info);

/* Double precision Triangular Band matrices */

/* DTBCON - estimate the reciprocal of the condition number of a */
/* triangular band matrix A, in either the 1-norm or the */
/* infinity-norm */


/** \private */
extern void
FORTRAN_NAME(dtbcon)(const char* norm, const char* uplo,
		 const char* diag, const int* n, const int* kd,
		 double const* ab, const int* ldab,
		 double* rcond, double* work,
		 int* iwork, int* info);
/* DTBRFS - provide error bounds and backward error estimates for */
/* the solution to a system of linear equations with a triangular */
/* band coefficient matrix */


/** \private */
extern void
FORTRAN_NAME(dtbrfs)(const char* uplo, const char* trans,
		 const char* diag, const int* n, const int* kd,
		 const int* nrhs,
		 double const* ab, const int* ldab,
		 double const* b, const int* ldb,
		 double* x, const int* ldx,
		 double* ferr, double* berr,
		 double* work, int* iwork, int* info);
/* DTBTRS - solve a triangular system of the form   A * X = B or */
/* A**T * X = B,  */


/** \private */
extern void
FORTRAN_NAME(dtbtrs)(const char* uplo, const char* trans,
		 const char* diag, const int* n,
		 const int* kd, const int* nrhs,
		 double const* ab, const int* ldab,
		 double* b, const int* ldb, int* info);

/* Double precision Triangular matrices Generalized problems */

/* DTGEVC - compute some or all of the right and/or left */
/* generalized eigenvectors of a pair of real upper triangular */
/* matrices (A,B); */


/** \private */
extern void
FORTRAN_NAME(dtgevc)(const char* side, const char* howmny,
		 const int* select, const int* n,
		 double const* a, const int* lda,
		 double const* b, const int* ldb,
		 double* vl, const int* ldvl,
		 double* vr, const int* ldvr,
		 const int* mm, int* m, double* work, int* info);

/* DTGSJA - compute the generalized singular value decomposition */
/* (GSVD); of two real upper triangular (or trapezoidal); matrices */
/* A and B */


/** \private */
extern void
FORTRAN_NAME(dtgsja)(const char* jobu, const char* jobv, const char* jobq,
		 const int* m, const int* p, const int* n,
		 const int* k, const int* l,
		 double* a, const int* lda,
		 double* b, const int* ldb,
		 double const* tola, double const* tolb,
		 double* alpha, double* beta,
		 double* u, const int* ldu,
		 double* v, const int* ldv,
		 double* q, const int* ldq,
		 double* work, int* ncycle, int* info);

/* Double precision Triangular matrices Packed storage */

/* DTPCON - estimate the reciprocal of the condition number of a */
/* packed triangular matrix A, in either the 1-norm or the */
/* infinity-norm */


/** \private */
extern void
FORTRAN_NAME(dtpcon)(const char* norm, const char* uplo,
		 const char* diag, const int* n,
		 double const* ap, double* rcond,
		 double* work, int* iwork, int* info);

/* DTPRFS - provide error bounds and backward error estimates for */
/* the solution to a system of linear equations with a triangular */
/* packed coefficient matrix */


/** \private */
extern void
FORTRAN_NAME(dtprfs)(const char* uplo, const char* trans,
		 const char* diag, const int* n,
		 const int* nrhs, double const* ap,
		 double const* b, const int* ldb,
		 double* x, const int* ldx,
		 double* ferr, double* berr,
		 double* work, int* iwork, int* info);

/* Double precision TRiangular matrices */

/* DTPTRI - compute the inverse of a real upper or lower */
/* triangular matrix A stored in packed format */


/** \private */
extern void
FORTRAN_NAME(dtptri)(const char* uplo, const char* diag,
		 const int* n, double* ap, int* info);

/* DTPTRS - solve a triangular system of the form   A * X = B or */
/* A**T * X = B, */


/** \private */
extern void
FORTRAN_NAME(dtptrs)(const char* uplo, const char* trans,
		 const char* diag, const int* n,
		 const int* nrhs, double const* ap,
		 double* b, const int* ldb, int* info);

/* DTRCON - estimate the reciprocal of the condition number of a */
/* triangular matrix A, in either the 1-norm or the infinity-norm */


/** \private */
extern void
FORTRAN_NAME(dtrcon)(const char* norm, const char* uplo,
		 const char* diag, const int* n,
		 double const* a, const int* lda,
		 double* rcond, double* work,
		 int* iwork, int* info);

/* DTREVC - compute some or all of the right and/or left */
/* eigenvectors of a real upper quasi-triangular matrix T */


/** \private */
extern void
FORTRAN_NAME(dtrevc)(const char* side, const char* howmny,
		 const int* select, const int* n,
		 double const* t, const int* ldt,
		 double* vl, const int* ldvl,
		 double* vr, const int* ldvr,
		 const int* mm, int* m, double* work, int* info);

/* DTREXC - reorder the real Schur factorization of a real matrix */
/* A = Q*T*Q**T, so that the diagonal block of T with row index */
/* IFST is moved to row ILST */


/** \private */
extern void
FORTRAN_NAME(dtrexc)(const char* compq, const int* n,
		 double* t, const int* ldt,
		 double* q, const int* ldq,
		 int* ifst, int* ILST,
		 double* work, int* info);

/* DTRRFS - provide error bounds and backward error estimates for */
/* the solution to a system of linear equations with a triangular */
/* coefficient matrix */


/** \private */
extern void
FORTRAN_NAME(dtrrfs)(const char* uplo, const char* trans,
		 const char* diag, const int* n, const int* nrhs,
		 double const* a, const int* lda,
		 double const* b, const int* ldb,
		 double* x, const int* ldx,
		 double* ferr, double* berr,
		 double* work, int* iwork, int* info);

/* DTRSEN - reorder the real Schur factorization of a real matrix */
/* A = Q*T*Q**T, so that a selected cluster of eigenvalues appears */
/* in the leading diagonal blocks of the upper quasi-triangular */
/* matrix T, */


/** \private */
extern void
FORTRAN_NAME(dtrsen)(const char* job, const char* compq,
		 const int* select, const int* n,
		 double* t, const int* ldt,
		 double* q, const int* ldq,
		 double* wr, double* wi,
		 int* m, double* s, double* sep,
		 double* work, const int* lwork,
		 int* iwork, const int* liwork, int* info);

/* DTRSNA - estimate reciprocal condition numbers for specified */
/* eigenvalues and/or right eigenvectors of a real upper */
/* quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q */
/* orthogonal); */


/** \private */
extern void
FORTRAN_NAME(dtrsna)(const char* job, const char* howmny,
		 const int* select, const int* n,
		 double const* t, const int* ldt,
		 double const* vl, const int* ldvl,
		 double const* vr, const int* ldvr,
		 double* s, double* sep, const int* mm,
		 int* m, double* work, const int* lwork,
		 int* iwork, int* info);

/* DTRSYL - solve the real Sylvester matrix equation */


/** \private */
extern void
FORTRAN_NAME(dtrsyl)(const char* trana, const char* tranb,
		 const int* isgn, const int* m, const int* n,
		 double const* a, const int* lda,
		 double const* b, const int* ldb,
		 double* c, const int* ldc,
		 double* scale, int* info);

/* DTRTI2 - compute the inverse of a real upper or lower */
/* triangular matrix */


/** \private */
extern void
FORTRAN_NAME(dtrti2)(const char* uplo, const char* diag,
		 const int* n, double* a, const int* lda,
		 int* info);

/* DTRTRI - compute the inverse of a real upper or lower */
/* triangular matrix A */


/** \private */
extern void
FORTRAN_NAME(dtrtri)(const char* uplo, const char* diag,
		 const int* n, double* a, const int* lda,
		 int* info);

/* DTRTRS - solve a triangular system of the form   A * X = B or */
/* A**T * X = B	 */


/** \private */
extern void
FORTRAN_NAME(dtrtrs)(const char* uplo, const char* trans,
		 const char* diag, const int* n, const int* nrhs,
		 double const* a, const int* lda,
		 double* b, const int* ldb, int* info);

/* DTZRQF - reduce the M-by-N ( M<=N ); real upper trapezoidal */
/* matrix A to upper triangular form by means of orthogonal */
/* transformations  */


/** \private */
extern void
FORTRAN_NAME(dtzrqf)(const int* m, const int* n,
		 double* a, const int* lda,
		 double* tau, int* info);



/* Double precision utilties in Lapack */
/* DHGEQZ - implement a single-/double-shift version of the QZ */
/* method for finding the generalized eigenvalues */
/* w(j);=(ALPHAR(j); + i*ALPHAI(j););/BETAR(j); of the equation */
/* det( A - w(i); B ); = 0  In addition, the pair A,B may be */
/* reduced to generalized Schur form */


/** \private */
extern void
FORTRAN_NAME(dhgeqz)(const char* job, const char* compq, const char* compz,
		 const int* n, const int *ILO, const int* IHI,
		 double* a, const int* lda,
		 double* b, const int* ldb,
		 double* alphar, double* alphai, double const* beta,
		 double* q, const int* ldq,
		 double* z, const int* ldz,
		 double* work, const int* lwork, int* info);
/* DHSEIN - use inverse iteration to find specified right and/or */
/* left eigenvectors of a real upper Hessenberg matrix H */


/** \private */
extern void
FORTRAN_NAME(dhsein)(const char* side, const char* eigsrc,
		 const char* initv, int* select,
		 const int* n, double* h, const int* ldh,
		 double* wr, double* wi,
		 double* vl, const int* ldvl,
		 double* vr, const int* ldvr,
		 const int* mm, int* m, double* work,
		 int* ifaill, int* ifailr, int* info);
/* DHSEQR - compute the eigenvalues of a real upper Hessenberg */
/* matrix H and, optionally, the matrices T and Z from the Schur */
/* decomposition H = Z T Z**T, where T is an upper */
/* quasi-triangular matrix (the Schur form);, and Z is the */
/* orthogonal matrix of Schur vectors */


/** \private */
extern void
FORTRAN_NAME(dhseqr)(const char* job, const char* compz, const int* n,
		 const int* ilo, const int* ihi,
		 double* h, const int* ldh,
		 double* wr, double* wi,
		 double* z, const int* ldz,
		 double* work, const int* lwork, int* info);
/* DLABAD - take as input the values computed by SLAMCH for */
/* underflow and overflow, and returns the square root of each of */
/* these values if the log of LARGE is sufficiently large */


/** \private */
extern void
FORTRAN_NAME(dlabad)(double* small, double* large);
/* DLABRD - reduce the first NB rows and columns of a real */
/* general m by n matrix A to upper or lower bidiagonal form by an */
/* orthogonal transformation Q' * A * P, and returns the matrices */
/* X and Y which are needed to apply the transformation to the */
/* unreduced part of A */


/** \private */
extern void
FORTRAN_NAME(dlabrd)(const int* m, const int* n, const int* nb,
		 double* a, const int* lda, double* d, double* e,
		 double* tauq, double* taup,
		 double* x, const int* ldx, double* y, const int* ldy);
/* DLACON - estimate the 1-norm of a square, real matrix A */


/** \private */
extern void
FORTRAN_NAME(dlacon)(const int* n, double* v, double* x,
		 int* isgn, double* est, int* kase);
/* DLACPY - copy all or part of a two-dimensional matrix A to */
/* another matrix B */


/** \private */
extern void
FORTRAN_NAME(dlacpy)(const char* uplo, const int* m, const int* n,
		 double const* a, const int* lda,
		 double* b, const int* ldb);
/* DLADIV - perform complex division in real arithmetic	 */


/** \private */
extern void
FORTRAN_NAME(dladiv)(double const* a, double const* b,
		 double const* c, double const* d,
		 double* p, double* q);
/* DLAE2 - compute the eigenvalues of a 2-by-2 symmetric matrix [ A B ] */
/*								[ B C ] */


/** \private */
extern void
FORTRAN_NAME(dlae2)(double const* a, double const* b, double const* c,
		double* rt1, double* rt2);
/* DLAEBZ - contain the iteration loops which compute and use the */
/* function N(w);, which is the count of eigenvalues of a */
/* symmetric tridiagonal matrix T less than or equal to its */
/* argument w  */


/** \private */
extern void
FORTRAN_NAME(dlaebz)(const int* ijob, const int* nitmax, const int* n,
		 const int* mmax, const int* minp, const int* nbmin,
		 double const* abstol, double const* reltol,
		 double const* pivmin, double* d, double* e,
		 double* e2, int* nval, double* ab, double* c,
		 int* mout, int* nab, double* work, int* iwork,
		 int* info);
/* DLAED0 - compute all eigenvalues and corresponding */
/* eigenvectors of a symmetric tridiagonal matrix using the divide */
/* and conquer method */


/** \private */
extern void
FORTRAN_NAME(dlaed0)(const int* icompq, const int* qsiz, const int* n,
		 double* d, double* e, double* q, const int* ldq,
		 double* qstore, const int* ldqs,
		 double* work, int* iwork, int* info);
/* DLAED1 - compute the updated eigensystem of a diagonal matrix */
/* after modification by a rank-one symmetric matrix */


/** \private */
extern void
FORTRAN_NAME(dlaed1)(const int* n, double* d, double* q, const int* ldq,
		 int* indxq, double const* rho, const int* cutpnt,
		 double* work, int* iwork, int* info);
/* DLAED2 - merge the two sets of eigenvalues together into a */
/* single sorted set */


/** \private */
extern void
FORTRAN_NAME(dlaed2)(const int* k, const int* n, double* d,
		 double* q, const int* ldq, int* indxq,
		 double* rho, const int* cutpnt, double* z,
		 double* dlamda, double* q2, const int *ldq2,
		 int* indxc, int* w, int* indxp, int* indx,
		 int* coltyp, int* info);
/* DLAED3 - find the roots of the secular equation, as defined by */
/* the values in double* d, W, and RHO, between KSTART and KSTOP */


/** \private */
extern void
FORTRAN_NAME(dlaed3)(const int* k, const int* kstart,
		 const int *kstop, const int* n,
		 double* d, double* q, const int* ldq,
		 double const* rho, const int* cutpnt,
		 double* dlamda, int* q2, const int* ldq2,
		 int* indxc, int* ctot, double* w,
		 double* s, const int* lds, int* info);
/* DLAED4 - subroutine computes the I-th updated eigenvalue of a */
/* symmetric rank-one modification to a diagonal matrix whose */
/* elements are given in the array d, and that	 D(i); < D(j); for */
/* i < j  and that RHO > 0 */


/** \private */
extern void
FORTRAN_NAME(dlaed4)(const int* n, const int* i, double const* d,
		 double const* z, double const* delta,
		 double const* rho, double* dlam, int* info);
/* DLAED5 - subroutine computes the I-th eigenvalue of a */
/* symmetric rank-one modification of a 2-by-2 diagonal matrix */
/* diag( D ); + RHO  The diagonal elements in the array D are */
/* assumed to satisfy	D(i); < D(j); for i < j	 */


/** \private */
extern void
FORTRAN_NAME(dlaed5)(const int* i, double const* d, double const* z,
		 double* delta, double const* rho, double* dlam);
/* DLAED6 - compute the positive or negative root (closest to the */
/* origin); of	z(1); z(2); z(3); f(x); = rho + --------- + */
/* ---------- + ---------  d(1);-x d(2);-x d(3);-x  It is assumed */
/* that	  if ORGATI = .true  */


/** \private */
extern void
FORTRAN_NAME(dlaed6)(const int* kniter, const int* orgati,
		 double const* rho, double const* d,
		 double const* z, double const* finit,
		 double* tau, int* info);
/* DLAED7 - compute the updated eigensystem of a diagonal matrix */
/* after modification by a rank-one symmetric matrix */


/** \private */
extern void
FORTRAN_NAME(dlaed7)(const int* icompq, const int* n,
		 const int* qsiz, const int* tlvls,
		 const int* curlvl, const int* curpbm,
		 double* d, double* q, const int* ldq,
		 int* indxq, double const* rho, const int* cutpnt,
		 double* qstore, double* qptr, const int* prmptr,
		 const int* perm, const int* givptr,
		 const int* givcol, double const* givnum,
		 double* work, int* iwork, int* info);
/* DLAED8 - merge the two sets of eigenvalues together into a */
/* single sorted set */


/** \private */
extern void
FORTRAN_NAME(dlaed8)(const int* icompq, const int* k,
		 const int* n, const int* qsiz,
		 double* d, double* q, const int* ldq,
		 const int* indxq, double* rho,
		 const int* cutpnt, double const* z,
		 double* dlamda, double* q2, const int* ldq2,
		 double* w, int* perm, int* givptr,
		 int* givcol, double* givnum, int* indxp,
		 int* indx, int* info);
/* DLAED9 - find the roots of the secular equation, as defined by */
/* the values in double* d, Z, and RHO, between KSTART and KSTOP */


/** \private */
extern void
FORTRAN_NAME(dlaed9)(const int* k, const int* kstart, const int* kstop,
		 const int* n, double* d, double* q, const int* ldq,
		 double const* rho, double const* dlamda,
		 double const* w, double* s, const int* lds, int* info);
/* DLAEDA - compute the Z vector corresponding to the merge step */
/* in the CURLVLth step of the merge process with TLVLS steps for */
/* the CURPBMth problem */


/** \private */
extern void
FORTRAN_NAME(dlaeda)(const int* n, const int* tlvls, const int* curlvl,
		 const int* curpbm, const int* prmptr, const int* perm,
		 const int* givptr, const int* givcol,
		 double const* givnum, double const* q,
		 const int* qptr, double* z, double* ztemp, int* info);
/* DLAEIN - use inverse iteration to find a right or left */
/* eigenvector corresponding to the eigenvalue (WR,WI); of a real */
/* upper Hessenberg matrix H */


/** \private */
extern void
FORTRAN_NAME(dlaein)(const int* rightv, const int* noinit, const int* n,
		 double const* h, const int* ldh,
		 double const* wr, double const* wi,
		 double* vr, double* vi,
		 double* b, const int* ldb, double* work,
		 double const* eps3, double const* smlnum,
		 double const* bignum, int* info);
/* DLAEV2 - compute the eigendecomposition of a 2-by-2 symmetric */
/* matrix  [ A B ]  [ B C ] */


/** \private */
extern void
FORTRAN_NAME(dlaev2)(double const* a, double const* b, double const* c,
		 double* rt1, double* rt2, double* cs1, double *sn1);
/* DLAEXC - swap adjacent diagonal blocks T11 and T22 of order 1 */
/* or 2 in an upper quasi-triangular matrix T by an orthogonal */
/* similarity transformation */


/** \private */
extern void
FORTRAN_NAME(dlaexc)(const int* wantq, const int* n, double* t, const int* ldt,
		  double* q, const int* ldq, const int* j1,
		 const int* n1, const int* n2, double* work, int* info);
/* DLAG2 - compute the eigenvalues of a 2 x 2 generalized */
/* eigenvalue problem A - w B, with scaling as necessary to aextern void */
/* over-/underflow */


/** \private */
extern void
FORTRAN_NAME(dlag2)(double const* a, const int* lda, double const* b,
		const int* ldb, double const* safmin,
		double* scale1, double* scale2,
		double* wr1, double* wr2, double* wi);
/* DLAGS2 - compute 2-by-2 orthogonal matrices U, V and Q, such */
/* that if ( UPPER ); then   U'*A*Q = U'*( A1 A2 );*Q = ( x 0 ); */
/* ( 0 A3 ); ( x x ); and  V'*B*Q = V'*( B1 B2 );*Q = ( x 0 );	( */
/* 0 B3 ); ( x x );  or if ( .NOT.UPPER ); then	  U'*A*Q = U'*( A1 */
/* 0 );*Q = ( x x );  ( A2 A3 ); ( 0 x ); and  V'*B*Q = V'*( B1 0 */
/* );*Q = ( x x );  ( B2 B3 ); ( 0 x );	 The rows of the */
/* transformed A and B are parallel, where   U = ( CSU SNU );, V = */
/* ( CSV SNV );, Q = ( CSQ SNQ );  ( -SNU CSU ); ( -SNV CSV ); ( */
/* -SNQ CSQ );	Z' denotes the transpose of Z */


/** \private */
extern void
FORTRAN_NAME(dlags2)(const int* upper,
		 double const* a1, double const* a2, double const* a3,
		 double const* b1, double const* b2, double const* b3,
		 double* csu, double* snu,
		 double* csv, double* snv, double *csq, double *snq);
/* DLAGTF - factorize the matrix (T - lambda*I);, where T is an n */
/* by n tridiagonal matrix and lambda is a scalar, as	T - */
/* lambda*I = PLU, */


/** \private */
extern void
FORTRAN_NAME(dlagtf)(const int* n, double* a, double const* lambda,
		 double* b, double* c, double const *tol,
		 double* d, int* in, int* info);
/* DLAGTM - perform a matrix-vector product of the form	  B := */
/* alpha * A * X + beta * B  where A is a tridiagonal matrix of */
/* order N, B and X are N by NRHS matrices, and alpha and beta are */
/* real scalars, each of which may be 0., 1., or -1 */


/** \private */
extern void
FORTRAN_NAME(dlagtm)(const char* trans, const int* n, const int* nrhs,
		 double const* alpha, double const* dl,
		 double const* d, double const* du,
		 double const* x, const int* ldx, double const* beta,
		 double* b, const int* ldb);
/* DLAGTS - may be used to solve one of the systems of equations */
/* (T - lambda*I);*x = y or (T - lambda*I);'*x = y, */


/** \private */
extern void
FORTRAN_NAME(dlagts)(const int* job, const int* n,
		 double const* a, double const* b,
		 double const* c, double const* d,
		 const int* in, double* y, double* tol, int* info);
/* DLAHQR - an auxiliary routine called by DHSEQR to update the */
/* eigenvalues and Schur decomposition already computed by DHSEQR, */
/* by dealing with the Hessenberg submatrix in rows and columns */
/* ILO to IHI */


/** \private */
extern void
FORTRAN_NAME(dlahqr)(const int* wantt, const int* wantz, const int* n,
		 const int* ilo, const int* ihi,
		 double* H, const int* ldh, double* wr, double* wi,
		 const int* iloz, const int* ihiz,
		 double* z, const int* ldz, int* info);
/* DLAHRD - reduce the first NB columns of a real general */
/* n-by-(n-k+1); matrix A so that elements below the k-th */
/* subdiagonal are zero */


/** \private */
extern void
FORTRAN_NAME(dlahrd)(const int* n, const int* k, const int* nb,
		 double* a, const int* lda,
		 double* tau, double* t, const int* ldt,
		 double* y, const int* ldy);
/* DLAIC1 - apply one step of incremental condition estimation in */
/* its simplest version */


/** \private */
extern void
FORTRAN_NAME(dlaic1)(const int* job, const int* j, double const* x,
		 double const* sest, double const* w,
		 double const* gamma, double* sestpr,
		 double* s, double* c);
/* DLALN2 - solve a system of the form (ca A - w D ); X = s B or */
/* (ca A' - w D); X = s B with possible scaling ("s"); and */
/* perturbation of A */


/** \private */
extern void
FORTRAN_NAME(dlaln2)(const int* ltrans, const int* na, const int* nw,
		 double const* smin, double const* ca,
		 double const* a, const int* lda,
		 double const* d1, double const* d2,
		 double const* b, const int* ldb,
		 double const* wr, double const* wi,
		 double* x, const int* ldx, double* scale,
		 double* xnorm, int* info);
/* DLAMCH - determine double precision machine parameters */


/** \private */
extern double
FORTRAN_NAME(dlamch)(const char* cmach);
/* DLAMRG - will create a permutation list which will merge the */
/* elements of A (which is composed of two independently sorted */
/* sets); into a single set which is sorted in ascending order */


/** \private */
extern void
FORTRAN_NAME(dlamrg)(const int* n1, const int* n2, double const* a,
		 const int* dtrd1, const int* dtrd2, int* index);
/* DLANGB - return the value of the one norm, or the Frobenius */
/* norm, or the infinity norm, or the element of largest absolute */
/* value of an n by n band matrix A, with kl sub-diagonals and ku */
/* super-diagonals */


/** \private */
extern double
FORTRAN_NAME(dlangb)(const char* norm, const int* n,
		 const int* kl, const int* ku, double const* ab,
		 const int* ldab, double* work);
/* DLANGE - return the value of the one norm, or the Frobenius */
/* norm, or the infinity norm, or the element of largest absolute */
/* value of a real matrix A */


/** \private */
extern double
FORTRAN_NAME(dlange)(const char* norm, const int* m, const int* n,
		 double const* a, const int* lda, double* work);
/* DLANGT - return the value of the one norm, or the Frobenius */
/* norm, or the infinity norm, or the element of largest absolute */
/* value of a real tridiagonal matrix A */


/** \private */
extern double
FORTRAN_NAME(dlangt)(const char* norm, const int* n,
		 double const* dl, double const* d,
		 double const* du);
/* DLANHS - return the value of the one norm, or the Frobenius */
/* norm, or the infinity norm, or the element of largest absolute */
/* value of a Hessenberg matrix A */


/** \private */
extern double
FORTRAN_NAME(dlanhs)(const char* norm, const int* n,
		 double const* a, const int* lda, double* work);
/* DLANSB - return the value of the one norm, or the Frobenius */
/* norm, or the infinity norm, or the element of largest absolute */
/* value of an n by n symmetric band matrix A, with k */
/* super-diagonals */


/** \private */
extern double
FORTRAN_NAME(dlansb)(const char* norm, const char* uplo,
		 const int* n, const int* k,
		 double const* ab, const int* ldab, double* work);
/* DLANSP - return the value of the one norm, or the Frobenius */
/* norm, or the infinity norm, or the element of largest absolute */
/* value of a real symmetric matrix A, supplied in packed form */


/** \private */
extern double
FORTRAN_NAME(dlansp)(const char* norm, const char* uplo,
		 const int* n, double const* ap, double* work);
/* DLANST - return the value of the one norm, or the Frobenius */
/* norm, or the infinity norm, or the element of largest absolute */
/* value of a real symmetric tridiagonal matrix A */


/** \private */
extern double
FORTRAN_NAME(dlanst)(const char* norm, const int* n,
		 double const* d, double const* e);
/* DLANSY - return the value of the one norm, or the Frobenius */
/* norm, or the infinity norm, or the element of largest absolute */
/* value of a real symmetric matrix A */


/** \private */
extern double
FORTRAN_NAME(dlansy)(const char* norm, const char* uplo, const int* n,
		 double const* a, const int* lda, double* work);
/* DLANTB - return the value of the one norm, or the Frobenius */
/* norm, or the infinity norm, or the element of largest absolute */
/* value of an n by n triangular band matrix A, with ( k + 1 ) diagonals */


/** \private */
extern double
FORTRAN_NAME(dlantb)(const char* norm, const char* uplo,
		 const char* diag, const int* n, const int* k,
		 double const* ab, const int* ldab, double* work);
/* DLANTP - return the value of the one norm, or the Frobenius */
/* norm, or the infinity norm, or the element of largest absolute */
/* value of a triangular matrix A, supplied in packed form */


/** \private */
extern double
FORTRAN_NAME(dlantp)(const char* norm, const char* uplo, const char* diag,
		 const int* n, double const* ap, double* work);
/* DLANTR - return the value of the one norm, or the Frobenius */
/* norm, or the infinity norm, or the element of largest absolute */
/* value of a trapezoidal or triangular matrix A */


/** \private */
extern double
FORTRAN_NAME(dlantr)(const char* norm, const char* uplo,
		 const char* diag, const int* m, const int* n,
		 double const* a, const int* lda, double* work);
/* DLANV2 - compute the Schur factorization of a real 2-by-2 */
/* nonsymmetric matrix in standard form */


/** \private */
extern void
FORTRAN_NAME(dlanv2)(double* a, double* b, double* c, double* d,
		 double* rt1r, double* rt1i, double* rt2r, double* rt2i,
		 double* cs, double *sn);
/* DLAPLL - two column vectors X and Y, let A = ( X Y ); */


/** \private */
extern void
FORTRAN_NAME(dlapll)(const int* n, double* x, const int* incx,
		 double* y, const int* incy, double* ssmin);
/* DLAPMT - rearrange the columns of the M by N matrix X as */
/* specified by the permutation K(1);,K(2);,...,K(N); of the */
/* integers 1,...,N */


/** \private */
extern void
FORTRAN_NAME(dlapmt)(const int* forwrd, const int* m, const int* n,
		 double* x, const int* ldx, const int* k);
/* DLAPY2 - return sqrt(x**2+y**2);, taking care not to cause */
/* unnecessary overflow */


/** \private */
extern double
FORTRAN_NAME(dlapy2)(double const* x, double const* y);
/* DLAPY3 - return sqrt(x**2+y**2+z**2);, taking care not to */
/* cause unnecessary overflow */


/** \private */
extern double
FORTRAN_NAME(dlapy3)(double const* x, double const* y, double const* z);
/* DLAQGB - equilibrate a general M by N band matrix A with KL */
/* subdiagonals and KU superdiagonals using the row and scaling */
/* factors in the vectors R and C */


/** \private */
extern void
FORTRAN_NAME(dlaqgb)(const int* m, const int* n,
		 const int* kl, const int* ku,
		 double* ab, const int* ldab,
		 double* r, double* c,
		 double* rowcnd, double* colcnd,
		 double const* amax, char* equed);
/* DLAQGE - equilibrate a general M by N matrix A using the row */
/* and scaling factors in the vectors R and C */


/** \private */
extern void
FORTRAN_NAME(dlaqge)(const int* m, const int* n,
		 double* a, const int* lda,
		 double* r, double* c,
		 double* rowcnd, double* colcnd,
		 double const* amax, char* equed);
/* DLAQSB - equilibrate a symmetric band matrix A using the */
/* scaling factors in the vector S */


/** \private */
extern void
FORTRAN_NAME(dlaqsb)(const char* uplo, const int* n, const int* kd,
		 double* ab, const int* ldab, double const* s,
		 double const* scond, double const* amax, char* equed);
/* DLAQSP - equilibrate a symmetric matrix A using the scaling */
/* factors in the vector S */


/** \private */
extern void
FORTRAN_NAME(dlaqsp)(const char* uplo, const int* n,
		 double* ap, double const* s, double const* scond,
		 double const* amax, int* equed);
/* DLAQSY - equilibrate a symmetric matrix A using the scaling */
/* factors in the vector S */


/** \private */
extern void
FORTRAN_NAME(dlaqsy)(const char* uplo, const int* n,
		 double* a, const int* lda,
		 double const* s, double const* scond,
		 double const* amax, int* equed);
/* DLAQTR - solve the real quasi-triangular system   */
/* op(T) * p = scale*c */


/** \private */
extern void
FORTRAN_NAME(dlaqtr)(const int* ltran, const int* lreal, const int* n,
		 double const* t, const int* ldt,
		 double const* b, double const* w,
		 double* scale, double* x, double* work, int* info);
/* DLAR2V - apply a vector of real plane rotations from both */
/* sides to a sequence of 2-by-2 real symmetric matrices, defined */
/* by the elements of the vectors x, y and z  */


/** \private */
extern void
FORTRAN_NAME(dlar2v)(const int* n, double* x, double* y,
		 double* z, const int* incx,
		 double const* c, double const* s,
		 const int* incc);
/* DLARF - apply a real elementary reflector H to a real m by n */
/* matrix C, from either the left or the right */


/** \private */
extern void
FORTRAN_NAME(dlarf)(const char* side, const int* m, const int* n,
		double const* v, const int* incv, double const* tau,
		double* c, const int* ldc, double* work);
/* DLARFB - apply a real block reflector H or its transpose H' */
/* to a real m by n matrix C, from either the left or the right */


/** \private */
extern void
FORTRAN_NAME(dlarfb)(const char* side, const char* trans,
		 const char* direct, const char* storev,
		 const int* m, const int* n, const int* k,
		 double const* v, const int* ldv,
		 double const* t, const int* ldt,
		 double* c, const int* ldc,
		 double* work, const int* lwork);
/* DLARFG - generate a real elementary reflector H of order n, */
/* such that   H * ( alpha ) = ( beta ), H' * H = I */


/** \private */
extern void
FORTRAN_NAME(dlarfg)(const int* n, double const* alpha,
		 double* x, const int* incx, double* tau);
/* DLARFT - form the triangular factor T of a real block */
/* reflector H of order n, which is defined as a product of k */
/* elementary reflectors */


/** \private */
extern void
FORTRAN_NAME(dlarft)(const char* direct, const char* storev,
		 const int* n, const int* k, double* v, const int* ldv,
		 double const* tau, double* t, const int* ldt);
/* DLARFX - apply a real elementary reflector H to a real m by n */
/* matrix C, from either the left or the right */


/** \private */
extern void
FORTRAN_NAME(dlarfx)(const char* side, const int* m, const int* n,
		 double const* v, double const* tau,
		 double* c, const int* ldc, double* work);
/* DLARGV - generate a vector of real plane rotations, determined */
/* by elements of the real vectors x and y */


/** \private */
extern void
FORTRAN_NAME(dlargv)(const int* n, double* x, const int* incx,
		 double* y, const int* incy, double* c, const int* incc);
/* DLARNV - return a vector of n random real numbers from a */
/* uniform or normal distribution */


/** \private */
extern void
FORTRAN_NAME(dlarnv)(const int* idist, int* iseed, const int* n, double* x);
/* DLARTG - generate a plane rotation so that	[ CS SN ]  */


/** \private */
extern void
FORTRAN_NAME(dlartg)(double const* f, double const* g, double* cs,
		 double* sn, double *r);
/* DLARTV - apply a vector of real plane rotations to elements of */
/* the real vectors x and y */


/** \private */
extern void
FORTRAN_NAME(dlartv)(const int* n, double* x, const int* incx,
		 double* y, const int* incy,
		 double const* c, double const* s,
		 const int* incc);
/* DLARUV - return a vector of n random real numbers from a */
/* uniform (0,1); */


/** \private */
extern void
FORTRAN_NAME(dlaruv)(int* iseed, const int* n, double* x);

/* DLAS2 - compute the singular values of the 2-by-2 matrix */
/* [ F G ]  [ 0 H ] */


/** \private */
extern void
FORTRAN_NAME(dlas2)(double const* f, double const* g, double const* h,
		 double* ssmin, double* ssmax);

/* DLASCL - multiply the M by N real matrix A by the real scalar */
/* CTO/CFROM */


/** \private */
extern void
FORTRAN_NAME(dlascl)(const char* type,
		 const int* kl,const int* ku,
		 double* cfrom, double* cto,
		 const int* m, const int* n,
		 double* a, const int* lda, int* info);

/* DLASET - initialize an m-by-n matrix A to BETA on the diagonal */
/* and ALPHA on the offdiagonals */


/** \private */
extern void
FORTRAN_NAME(dlaset)(const char* uplo, const int* m, const int* n,
		 double const* alpha, double const* beta,
		 double* a, const int* lda);
/* DLASQ1 - DLASQ1 computes the singular values of a real N-by-N */
/* bidiagonal  matrix with diagonal D and off-diagonal E */


/** \private */
extern void
FORTRAN_NAME(dlasq1)(const int* n, double* d, double* e,
		 double* work, int* info);
/* DLASQ2 - DLASQ2 computes the singular values of a real N-by-N */
/* unreduced  bidiagonal matrix with squared diagonal elements in */
/* Q and  squared off-diagonal elements in E */


/** \private */
extern void
FORTRAN_NAME(dlasq2)(const int* m, double* q, double* e,
		 double* qq, double* ee, double const* eps,
		 double const* tol2, double const* small2,
		 double* sup, int* kend, int* info);
/* DLASQ3 - DLASQ3 is the workhorse of the whole bidiagonal SVD */
/* algorithm */


/** \private */
extern void
FORTRAN_NAME(dlasq3)(int* n, double* q, double* e, double* qq,
		 double* ee, double* sup, double *sigma,
		 int* kend, int* off, int* iphase,
		 const int* iconv, double const* eps,
		 double const* tol2, double const* small2);
/* DLASQ4 - DLASQ4 estimates TAU, the smallest eigenvalue of a */
/* matrix */


/** \private */
extern void
FORTRAN_NAME(dlasq4)(const int* n, double const* q, double const* e,
		 double* tau, double* sup);
/* DLASR - perform the transformation	A := P*A, when SIDE = 'L' */
/* or 'l' ( Left-hand side );	A := A*P', when SIDE = 'R' or 'r' */
/* ( Right-hand side );	 where A is an m by n real matrix and P is */
/* an orthogonal matrix, */


/** \private */
extern void
FORTRAN_NAME(dlasr)(const char* side, const char* pivot,
		const char* direct, const int* m, const int* n,
		double const* c, double const* s,
		double* a, const int* lda);
/* DLASRT - the numbers in D in increasing order (if ID = 'I'); */
/* or in decreasing order (if ID = 'D' ); */


/** \private */
extern void
FORTRAN_NAME(dlasrt)(const char* id, const int* n, double* d, int* info);
/* DLASSQ - return the values scl and smsq such that   ( scl**2 */
/* );*smsq = x( 1 );**2 +...+ x( n );**2 + ( scale**2 );*sumsq, */


/** \private */
extern void
FORTRAN_NAME(dlassq)(const int* n, double const* x, const int* incx,
		 double* scale, double* sumsq);
/* DLASV2 - compute the singular value decomposition of a 2-by-2 */
/* triangular matrix  [ F G ]  [ 0 H ] */


/** \private */
extern void
FORTRAN_NAME(dlasv2)(double const* f, double const* g, double const* h,
		 double* ssmin, double* ssmax, double* snr, double* csr,
		 double* snl, double* csl);
/* DLASWP - perform a series of row interchanges on the matrix A */


/** \private */
extern void
FORTRAN_NAME(dlaswp)(const int* n, double* a, const int* lda,
		 const int* k1, const int* k2,
		 const int* ipiv, const int* incx);
/* DLASY2 - solve for the N1 by N2 matrix double* x, 1 <= N1,N2 <= 2, in */
/* op(TL);*X + ISGN*X*op(TR); = SCALE*B, */


/** \private */
extern void
FORTRAN_NAME(dlasy2)(const int* ltranl, const int* ltranr,
		 const int* isgn, const int* n1, const int* n2,
		 double const* tl, const int* ldtl,
		 double const* tr, const int* ldtr,
		 double const* b, const int* ldb,
		 double* scale, double* x, const int* ldx,
		 double* xnorm, int* info);
/* DLASYF - compute a partial factorization of a real symmetric */
/* matrix A using the Bunch-Kaufman diagonal pivoting method */


/** \private */
extern void
FORTRAN_NAME(dlasyf)(const char* uplo, const int* n,
		 const int* nb, const int* kb,
		 double* a, const int* lda, int* ipiv,
		 double* w, const int* ldw, int* info);
/* DLATBS - solve one of the triangular systems	  A *x = s*b or */
/* A'*x = s*b  with scaling to prevent overflow, where A is an */
/* upper or lower triangular band matrix */


/** \private */
extern void
FORTRAN_NAME(dlatbs)(const char* uplo, const char* trans,
		 const char* diag, const char* normin,
		 const int* n, const int* kd,
		 double const* ab, const int* ldab,
		 double* x, double* scale, double* cnorm, int* info);
/* DLATPS - solve one of the triangular systems	  A *x = s*b or */
/* A'*x = s*b  with scaling to prevent overflow, where A is an */
/* upper or lower triangular matrix stored in packed form */


/** \private */
extern void
FORTRAN_NAME(dlatps)(const char* uplo, const char* trans,
		 const char* diag, const char* normin,
		 const int* n, double const* ap,
		 double* x, double* scale, double* cnorm, int* info);
/* DLATRD - reduce NB rows and columns of a real symmetric matrix */
/* A to symmetric tridiagonal form by an orthogonal similarity */
/* transformation Q' * A * Q, and returns the matrices V and W */
/* which are needed to apply the transformation to the unreduced */
/* part of A */


/** \private */
extern void
FORTRAN_NAME(dlatrd)(const char* uplo, const int* n, const int* nb,
		 double* a, const int* lda, double* e, double* tau,
		 double* w, const int* ldw);
/* DLATRS - solve one of the triangular systems	  A *x = s*b or */
/* A'*x = s*b  with scaling to prevent overflow */


/** \private */
extern void
FORTRAN_NAME(dlatrs)(const char* uplo, const char* trans,
		 const char* diag, const char* normin,
		 const int* n, double const* a, const int* lda,
		 double* x, double* scale, double* cnorm, int* info);
/* DLATZM - apply a Householder matrix generated by DTZRQF to a */
/* matrix */


/** \private */
extern void
FORTRAN_NAME(dlatzm)(const char* side, const int* m, const int* n,
		 double const* v, const int* incv,
		 double const* tau, double* c1, double* c2,
		 const int* ldc, double* work);
/* DLAUU2 - compute the product U * U' or L' * const int* l, where the */
/* triangular factor U or L is stored in the upper or lower */
/* triangular part of the array A */


/** \private */
extern void
FORTRAN_NAME(dlauu2)(const char* uplo, const int* n,
		 double* a, const int* lda, int* info);
/* DLAUUM - compute the product U * U' or L' * L, where the */
/* triangular factor U or L is stored in the upper or lower */
/* triangular part of the array A */


/** \private */
extern void
FORTRAN_NAME(dlauum)(const char* uplo, const int* n,
		 double* a, const int* lda, int* info);


/* ======================================================================== */

/* Selected Double Complex Lapack Routines
   ========
 */

/* IZMAX1 finds the index of the element whose real part has maximum
 * absolute value. */


/** \private */
extern int
FORTRAN_NAME(izmax1)(const int *n, ComplexT *cx, const int *incx);


/*  ZGECON estimates the reciprocal of the condition number of a general
 *  complex matrix A, in either the 1-norm or the infinity-norm, using
 *  the LU factorization computed by ZGETRF.
 */


/** \private */
extern void
FORTRAN_NAME(zgecon)(const char *norm, const int *n,
		 const ComplexT *a, const int *lda,
		 double const *anorm, double *rcond,
		 ComplexT *work, double *rwork, int *info);

/* ZGESV computes the solution to a complex system of linear equations */


/** \private */
extern void
FORTRAN_NAME(zgesv)(const int *n, const int *nrhs, ComplexT *a,
		const int *lda, int *ipiv, ComplexT *b,
		const int *ldb, int *info);

/*  ZGEQP3 computes a QR factorization with column pivoting */


/** \private */
extern void
FORTRAN_NAME(zgeqp3)(const int *m, const int *n,
		 ComplexT *a, const int *lda,
		 int *jpvt, ComplexT *tau,
		 ComplexT *work, const int *lwork,
		 double *rwork, int *info);

/* ZUNMQR applies Q or Q**H from the Left or Right */


/** \private */
extern void
FORTRAN_NAME(zunmqr)(const char *side, const char *trans,
		 const int *m, const int *n, const int *k,
		 ComplexT *a, const int *lda,
		 ComplexT *tau,
		 ComplexT *c, const int *ldc,
		 ComplexT *work, const int *lwork, int *info);

/*  ZTRTRS solves triangular systems */


/** \private */
extern void
FORTRAN_NAME(ztrtrs)(const char *uplo, const char *trans, const char *diag,
		 const int *n, const int *nrhs,
		 ComplexT *a, const int *lda,
		 ComplexT *b, const int *ldb, int *info);
/* ZGESVD - compute the singular value decomposition (SVD); of a   */
/* real M-by-N matrix A, optionally computing the left and/or	   */
/* right singular vectors					   */


/** \private */
extern void
FORTRAN_NAME(zgesvd)(const char *jobu, const char *jobvt,
		 const int *m, const int *n,
		 ComplexT *a, const int *lda, double *s,
		 ComplexT *u, const int *ldu,
		 ComplexT *vt, const int *ldvt,
		 ComplexT *work, const int *lwork, double *rwork,
		 int *info);

/* ZGHEEV - compute all eigenvalues and, optionally, eigenvectors */
/* of a Hermitian matrix A */


/** \private */
extern void
FORTRAN_NAME(zheev)(const char *jobz, const char *uplo,
		const int *n, ComplexT *a, const int *lda,
		double *w, ComplexT *work, const int *lwork,
		double *rwork, int *info);

/* ZGGEEV - compute all eigenvalues and, optionally, eigenvectors */
/* of a complex non-symmetric matrix A */


/** \private */
extern void
FORTRAN_NAME(zgeev)(const char *jobvl, const char *jobvr,
		const int *n, ComplexT *a, const int *lda,
		ComplexT *wr, ComplexT *vl, const int *ldvl,
		ComplexT *vr, const int *ldvr,
		ComplexT *work, const int *lwork,
		double *rwork, int *info);


/* NOTE: The following entry points were traditionally in this file,
   but are not provided by R's libRlapack */

/* DZSUM1 - take the sum of the absolute values of a complex */
/* vector and returns a double precision result	 */


/** \private */
extern double
FORTRAN_NAME(dzsum1)(const int *n, ComplexT *CX, const int *incx);

/*  ZLACN2 estimates the 1-norm of a square, complex matrix A.
 *  Reverse communication is used for evaluating matrix-vector products.
*/


/** \private */
extern void
FORTRAN_NAME(zlacn2)(const int *n, ComplexT *v, ComplexT *x,
                 double *est, int *kase, int *isave);

/* ZLANTR  -  return the value of the one norm, or the Frobenius norm, */
/* or the infinity norm, or the element of largest absolute value of */
/* a trapezoidal or triangular matrix A */


/** \private */
extern double
FORTRAN_NAME(zlantr)(const char *norm, const char *uplo, const char *diag,
		 const int *m, const int *n, ComplexT *a,
		 const int *lda, double *work);

/* ======================================================================== */

/* Other double precision and double complex Lapack routines
   provided by libRlapack.

   These are extracted from the CLAPACK headers.
*/



/** \private */
extern void
FORTRAN_NAME(dbdsdc)(char *uplo, char *compq, int *n, double *
	d, double *e, double *u, int *ldu, double *vt,
	int *ldvt, double *q, int *iq, double *work, int * iwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dgegs)(char *jobvsl, char *jobvsr, int *n,
	double *a, int *lda, double *b, int *ldb, double *
	alphar, double *alphai, double *beta, double *vsl,
	int *ldvsl, double *vsr, int *ldvsr, double *work,
	int *lwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dgelsd)(int *m, int *n, int *nrhs,
	double *a, int *lda, double *b, int *ldb, double *
	s, double *rcond, int *rank, double *work, int *lwork,
	 int *iwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dgelsx)(int *m, int *n, int *nrhs,
	double *a, int *lda, double *b, int *ldb, int *
	jpvt, double *rcond, int *rank, double *work, int *
	info);



/** \private */
extern void
FORTRAN_NAME(dgesc2)(int *n, double *a, int *lda,
	double *rhs, int *ipiv, int *jpiv, double *scale);

/* DGESDD - compute the singular value decomposition (SVD); of a   */
/* real M-by-N matrix A, optionally computing the left and/or	   */
/* right singular vectors.  If singular vectors are desired, it uses a */
/* divide-and-conquer algorithm.				   */


/** \private */
extern void
FORTRAN_NAME(dgesdd)(const char *jobz,
		 const int *m, const int *n,
		 double *a, const int *lda, double *s,
		 double *u, const int *ldu,
		 double *vt, const int *ldvt,
		 double *work, const int *lwork, int *iwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dgetc2)(int *n, double *a, int *lda, int
	*ipiv, int *jpiv, int *info);

typedef int (*L_fp)();


/** \private */
extern void
FORTRAN_NAME(dggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp
	delctg, char *sense, int *n, double *a, int *lda,
	double *b, int *ldb, int *sdim, double *alphar,
	double *alphai, double *beta, double *vsl, int *ldvsl,
	 double *vsr, int *ldvsr, double *rconde, double *
	rcondv, double *work, int *lwork, int *iwork, int *
	liwork, int *bwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dggev)(char *jobvl, char *jobvr, int *n, double *
	a, int *lda, double *b, int *ldb, double *alphar,
	double *alphai, double *beta, double *vl, int *ldvl,
	double *vr, int *ldvr, double *work, int *lwork,
	int *info);



/** \private */
extern void
FORTRAN_NAME(dggevx)(char *balanc, char *jobvl, char *jobvr, char *
	sense, int *n, double *a, int *lda, double *b,
	int *ldb, double *alphar, double *alphai, double *
	beta, double *vl, int *ldvl, double *vr, int *ldvr,
	int *ilo, int *ihi, double *lscale, double *rscale,
	double *abnrm, double *bbnrm, double *rconde, double *
	rcondv, double *work, int *lwork, int *iwork, int *
	bwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dggsvp)(char *jobu, char *jobv, char *jobq, int *m,
	int *p, int *n, double *a, int *lda, double *b,
	int *ldb, double *tola, double *tolb, int *k, int
	*l, double *u, int *ldu, double *v, int *ldv,
	double *q, int *ldq, int *iwork, double *tau,
	double *work, int *info);



/** \private */
extern void
FORTRAN_NAME(dgtts2)(int *itrans, int *n, int *nrhs,
	double *dl, double *d, double *du, double *du2,
	int *ipiv, double *b, int *ldb);


/** \private */
extern void
FORTRAN_NAME(dlagv2)(double *a, int *lda, double *b, int *ldb, double *alphar,
		 double *alphai, double * beta, double *csl, double *snl,
		 double *csr, double * snr);



/** \private */
extern void
FORTRAN_NAME(dlals0)(int *icompq, int *nl, int *nr,
	int *sqre, int *nrhs, double *b, int *ldb, double
	*bx, int *ldbx, int *perm, int *givptr, int *givcol,
	int *ldgcol, double *givnum, int *ldgnum, double *
	poles, double *difl, double *difr, double *z, int *
	k, double *c, double *s, double *work, int *info);



/** \private */
extern void
FORTRAN_NAME(dlalsa)(int *icompq, int *smlsiz, int *n,
	int *nrhs, double *b, int *ldb, double *bx, int *
	ldbx, double *u, int *ldu, double *vt, int *k,
	double *difl, double *difr, double *z, double *
	poles, int *givptr, int *givcol, int *ldgcol, int *
	perm, double *givnum, double *c, double *s, double *
	work, int *iwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dlalsd)(char *uplo, int *smlsiz, int *n, int
	*nrhs, double *d, double *e, double *b, int *ldb,
	double *rcond, int *rank, double *work, int *iwork,
	int *info);



/** \private */
extern void
FORTRAN_NAME(dlamc1)(int *beta, int *t, int *rnd, int
	*ieee1);



/** \private */
extern void
FORTRAN_NAME(dlamc2)(int *beta, int *t, int *rnd,
	double *eps, int *emin, double *rmin, int *emax,
	double *rmax);



/** \private */
extern double
FORTRAN_NAME(dlamc3)(double *a, double *b);



/** \private */
extern void
FORTRAN_NAME(dlamc4)(int *emin, double *start, int *base);



/** \private */
extern void
FORTRAN_NAME(dlamc5)(int *beta, int *p, int *emin,
	int *ieee, int *emax, double *rmax);



/** \private */
extern void
FORTRAN_NAME(dlaqp2)(int *m, int *n, int *offset,
	double *a, int *lda, int *jpvt, double *tau,
	double *vn1, double *vn2, double *work);



/** \private */
extern void
FORTRAN_NAME(dlaqps)(int *m, int *n, int *offset, int
	*nb, int *kb, double *a, int *lda, int *jpvt,
	double *tau, double *vn1, double *vn2, double *auxv,
	double *f, int *ldf);



/** \private */
extern void
FORTRAN_NAME(dlar1v)(int *n, int *b1, int *bn, double
	*sigma, double *d, double *l, double *ld, double *
	lld, double *gersch, double *z, double *ztz, double
	*mingma, int *r, int *isuppz, double *work);



/** \private */
extern void
FORTRAN_NAME(dlarrb)(int *n, double *d, double *l,
	double *ld, double *lld, int *ifirst, int *ilast,
	double *sigma, double *reltol, double *w, double *
	wgap, double *werr, double *work, int *iwork, int *
	info);



/** \private */
extern void
FORTRAN_NAME(dlarre)(int *n, double *d, double *e,
	double *tol, int *nsplit, int *isplit, int *m,
	double *w, double *woff, double *gersch, double *work,
	 int *info);



/** \private */
extern void
FORTRAN_NAME(dlarrf)(int *n, double *d, double *l,
	double *ld, double *lld, int *ifirst, int *ilast,
	double *w, double *dplus, double *lplus, double *work,
	 int *iwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dlarrv)(int *n, double *d, double *l,
	int *isplit, int *m, double *w, int *iblock,
	double *gersch, double *tol, double *z, int *ldz,
	int *isuppz, double *work, int *iwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dlarz)(char *side, int *m, int *n, int *l,
	double *v, int *incv, double *tau, double *c,
	int *ldc, double *work);



/** \private */
extern void
FORTRAN_NAME(dlarzb)(char *side, char *trans, char *direct, char *
	storev, int *m, int *n, int *k, int *l, double *v,
	 int *ldv, double *t, int *ldt, double *c, int *
	ldc, double *work, int *ldwork);



/** \private */
extern void
FORTRAN_NAME(dlarzt)(char *direct, char *storev, int *n, int *
	k, double *v, int *ldv, double *tau, double *t,
	int *ldt);



/** \private */
extern void
FORTRAN_NAME(dlasd0)(int *n, int *sqre, double *d,
	double *e, double *u, int *ldu, double *vt, int *
	ldvt, int *smlsiz, int *iwork, double *work, int *
	info);



/** \private */
extern void
FORTRAN_NAME(dlasd1)(int *nl, int *nr, int *sqre,
	double *d, double *alpha, double *beta, double *u,
	int *ldu, double *vt, int *ldvt, int *idxq, int *
	iwork, double *work, int *info);



/** \private */
extern void
FORTRAN_NAME(dlasd2)(int *nl, int *nr, int *sqre, int
	*k, double *d, double *z, double *alpha, double *
	beta, double *u, int *ldu, double *vt, int *ldvt,
	double *dsigma, double *u2, int *ldu2, double *vt2,
	int *ldvt2, int *idxp, int *idx, int *idxc, int *
	idxq, int *coltyp, int *info);



/** \private */
extern void
FORTRAN_NAME(dlasd3)(int *nl, int *nr, int *sqre, int
	*k, double *d, double *q, int *ldq, double *dsigma,
	double *u, int *ldu, double *u2, int *ldu2,
	double *vt, int *ldvt, double *vt2, int *ldvt2,
	int *idxc, int *ctot, double *z, int *info);



/** \private */
extern void
FORTRAN_NAME(dlasd4)(int *n, int *i, double *d,
	double *z, double *delta, double *rho, double *
	sigma, double *work, int *info);



/** \private */
extern void
FORTRAN_NAME(dlasd5)(int *i, double *d, double *z,
	double *delta, double *rho, double *dsigma, double *
	work);



/** \private */
extern void
FORTRAN_NAME(dlasd6)(int *icompq, int *nl, int *nr,
	int *sqre, double *d, double *vf, double *vl,
	double *alpha, double *beta, int *idxq, int *perm,
	int *givptr, int *givcol, int *ldgcol, double *givnum,
	 int *ldgnum, double *poles, double *difl, double *
	difr, double *z, int *k, double *c, double *s,
	double *work, int *iwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dlasd7)(int *icompq, int *nl, int *nr,
	int *sqre, int *k, double *d, double *z,
	double *zw, double *vf, double *vfw, double *vl,
	double *vlw, double *alpha, double *beta, double *
	dsigma, int *idx, int *idxp, int *idxq, int *perm,
	int *givptr, int *givcol, int *ldgcol, double *givnum,
	 int *ldgnum, double *c, double *s, int *info);



/** \private */
extern void
FORTRAN_NAME(dlasd8)(int *icompq, int *k, double *d,
	double *z, double *vf, double *vl, double *difl,
	double *difr, int *lddifr, double *dsigma, double *
	work, int *info);



/** \private */
extern void
FORTRAN_NAME(dlasd9)(int *icompq, int *ldu, int *k,
	double *d, double *z, double *vf, double *vl,
	double *difl, double *difr, double *dsigma, double *
	work, int *info);



/** \private */
extern void
FORTRAN_NAME(dlasda)(int *icompq, int *smlsiz, int *n,
	int *sqre, double *d, double *e, double *u, int
	*ldu, double *vt, int *k, double *difl, double *difr,
	double *z, double *poles, int *givptr, int *givcol,
	int *ldgcol, int *perm, double *givnum, double *c,
	double *s, double *work, int *iwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dlasdq)(char *uplo, int *sqre, int *n, int *
	ncvt, int *nru, int *ncc, double *d, double *e,
	double *vt, int *ldvt, double *u, int *ldu,
	double *c, int *ldc, double *work, int *info);



/** \private */
extern void
FORTRAN_NAME(dlasdt)(int *n, int *lvl, int *nd, int *
	inode, int *ndiml, int *ndimr, int *msub);



/** \private */
extern void
FORTRAN_NAME(dlasq5)(int *i0, int *n0, double *z,
	int *pp, double *tau, double *dmin, double *dmin1,
	double *dmin2, double *dn, double *dnm1, double *dnm2,
	 int *ieee);



/** \private */
extern void
FORTRAN_NAME(dlasq6)(int *i0, int *n0, double *z,
	int *pp, double *dmin, double *dmin1, double *dmin2,
	 double *dn, double *dnm1, double *dnm2);



/** \private */
extern void
FORTRAN_NAME(dlatdf)(int *ijob, int *n, double *z,
	int *ldz, double *rhs, double *rdsum, double *rdscal,
	int *ipiv, int *jpiv);



/** \private */
extern void
FORTRAN_NAME(dlatrz)(int *m, int *n, int *l, double *
	a, int *lda, double *tau, double *work);



/** \private */
extern void
FORTRAN_NAME(dormr3)(char *side, char *trans, int *m, int *n,
	int *k, int *l, double *a, int *lda, double *tau,
	double *c, int *ldc, double *work, int *info);



/** \private */
extern void
FORTRAN_NAME(dormrz)(char *side, char *trans, int *m, int *n,
	int *k, int *l, double *a, int *lda, double *tau,
	double *c, int *ldc, double *work, int *lwork,
	int *info);



/** \private */
extern void
FORTRAN_NAME(dptts2)(int *n, int *nrhs, double *d,
	double *e, double *b, int *ldb);



/** \private */
extern void
FORTRAN_NAME(dsbgvd)(char *jobz, char *uplo, int *n, int *ka,
	int *kb, double *ab, int *ldab, double *bb, int *
	ldbb, double *w, double *z, int *ldz, double *work,
	int *lwork, int *iwork, int *liwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dsbgvx)(char *jobz, char *range, char *uplo, int *n,
	int *ka, int *kb, double *ab, int *ldab, double *
	bb, int *ldbb, double *q, int *ldq, double *vl,
	double *vu, int *il, int *iu, double *abstol, int
	*m, double *w, double *z, int *ldz, double *work,
	int *iwork, int *ifail, int *info);



/** \private */
extern void
FORTRAN_NAME(dspgvd)(int *itype, char *jobz, char *uplo, int *
	n, double *ap, double *bp, double *w, double *z,
	int *ldz, double *work, int *lwork, int *iwork,
	int *liwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dspgvx)(int *itype, char *jobz, char *range, char *
	uplo, int *n, double *ap, double *bp, double *vl,
	double *vu, int *il, int *iu, double *abstol, int
	*m, double *w, double *z, int *ldz, double *work,
	int *iwork, int *ifail, int *info);



/** \private */
extern void
FORTRAN_NAME(dstegr)(char *jobz, char *range, int *n, double *
	d, double *e, double *vl, double *vu, int *il,
	int *iu, double *abstol, int *m, double *w,
	double *z, int *ldz, int *isuppz, double *work,
	int *lwork, int *iwork, int *liwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dstevr)(char *jobz, char *range, int *n, double *
	d, double *e, double *vl, double *vu, int *il,
	int *iu, double *abstol, int *m, double *w,
	double *z, int *ldz, int *isuppz, double *work,
	int *lwork, int *iwork, int *liwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dsygvd)(int *itype, char *jobz, char *uplo, int *
	n, double *a, int *lda, double *b, int *ldb,
	double *w, double *work, int *lwork, int *iwork,
	int *liwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dsygvx)(int *itype, char *jobz, char *range, char *
	uplo, int *n, double *a, int *lda, double *b, int
	*ldb, double *vl, double *vu, int *il, int *iu,
	double *abstol, int *m, double *w, double *z,
	int *ldz, double *work, int *lwork, int *iwork,
	int *ifail, int *info);



/** \private */
extern void
FORTRAN_NAME(dtgex2)(int *wantq, int *wantz, int *n,
	double *a, int *lda, double *b, int *ldb, double *
	q, int *ldq, double *z, int *ldz, int *j1, int *
	n1, int *n2, double *work, int *lwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dtgexc)(int *wantq, int *wantz, int *n,
	double *a, int *lda, double *b, int *ldb, double *
	q, int *ldq, double *z, int *ldz, int *ifst,
	int *ilst, double *work, int *lwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dtgsen)(int *ijob, int *wantq, int *wantz,
	int *select, int *n, double *a, int *lda, double *
	b, int *ldb, double *alphar, double *alphai, double *
	beta, double *q, int *ldq, double *z, int *ldz,
	int *m, double *pl, double *pr, double *dif,
	double *work, int *lwork, int *iwork, int *liwork,
	int *info);



/** \private */
extern void
FORTRAN_NAME(dtgsna)(char *job, char *howmny, int *select,
	int *n, double *a, int *lda, double *b, int *ldb,
	double *vl, int *ldvl, double *vr, int *ldvr,
	double *s, double *dif, int *mm, int *m, double *
	work, int *lwork, int *iwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dtgsy2)(char *trans, int *ijob, int *m, int *
	n, double *a, int *lda, double *b, int *ldb,
	double *c, int *ldc, double *d, int *ldd,
	double *e, int *lde, double *f, int *ldf, double *
	scale, double *rdsum, double *rdscal, int *iwork, int
	*pq, int *info);



/** \private */
extern void
FORTRAN_NAME(dtgsyl)(char *trans, int *ijob, int *m, int *
	n, double *a, int *lda, double *b, int *ldb,
	double *c, int *ldc, double *d, int *ldd,
	double *e, int *lde, double *f, int *ldf, double *
	scale, double *dif, double *work, int *lwork, int *
	iwork, int *info);



/** \private */
extern void
FORTRAN_NAME(dtzrzf)(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *lwork, int *info);




/** \private */
extern int
FORTRAN_NAME(lsame)(char *ca, char *cb);



/** \private */
extern void
FORTRAN_NAME(zbdsqr)(char *uplo, int *n, int *ncvt, int *
	nru, int *ncc, double *d, double *e, ComplexT *vt,
	int *ldvt, ComplexT *u, int *ldu, ComplexT *c,
	int *ldc, double *rwork, int *info);



/** \private */
extern void
FORTRAN_NAME(zdrot)(int *n, ComplexT *cx, int *incx,
	ComplexT *cy, int *incy, double *c, double *s);



/** \private */
extern void
FORTRAN_NAME(zgebak)(char *job, char *side, int *n, int *ilo,
	int *ihi, double *scale, int *m, ComplexT *v,
	int *ldv, int *info);



/** \private */
extern void
FORTRAN_NAME(zgebal)(char *job, int *n, ComplexT *a, int
	*lda, int *ilo, int *ihi, double *scale, int *info);



/** \private */
extern void
FORTRAN_NAME(zgebd2)(int *m, int *n, ComplexT *a,
	int *lda, double *d, double *e, ComplexT *tauq,
	ComplexT *taup, ComplexT *work, int *info);



/** \private */
extern void
FORTRAN_NAME(zgebrd)(int *m, int *n, ComplexT *a,
	int *lda, double *d, double *e, ComplexT *tauq,
	ComplexT *taup, ComplexT *work, int *lwork, int *
	info);


/** \private */
extern void
FORTRAN_NAME(zgehd2)(int *n, int *ilo, int *ihi,
	ComplexT *a, int *lda, ComplexT *tau, ComplexT *
	work, int *info);



/** \private */
extern void
FORTRAN_NAME(zgehrd)(int *n, int *ilo, int *ihi,
	ComplexT *a, int *lda, ComplexT *tau, ComplexT *
	work, int *lwork, int *info);



/** \private */
extern void
FORTRAN_NAME(zgelq2)(int *m, int *n, ComplexT *a,
	int *lda, ComplexT *tau, ComplexT *work, int *info);



/** \private */
extern void
FORTRAN_NAME(zgelqf)(int *m, int *n, ComplexT *a,
	int *lda, ComplexT *tau, ComplexT *work, int *lwork,
	 int *info);



/** \private */
extern void
FORTRAN_NAME(zgeqr2)(int *m, int *n, ComplexT *a,
	int *lda, ComplexT *tau, ComplexT *work, int *info);



/** \private */
extern void
FORTRAN_NAME(zgeqrf)(int *m, int *n, ComplexT *a,
		 int *lda, ComplexT *tau, ComplexT *work, int *lwork,
		 int *info);



/** \private */
extern void
FORTRAN_NAME(zgetf2)(int *m, int *n, ComplexT *a,
	int *lda, int *ipiv, int *info);



/** \private */
extern void
FORTRAN_NAME(zgetrf)(int *m, int *n, ComplexT *a,
	int *lda, int *ipiv, int *info);



/** \private */
extern void
FORTRAN_NAME(zgetrs)(char *trans, int *n, int *nrhs,
	ComplexT *a, int *lda, int *ipiv, ComplexT *b,
	int *ldb, int *info);




/** \private */
extern void
FORTRAN_NAME(zhetd2)(char *uplo, int *n, ComplexT *a, int *lda, double *d,
		 double *e, ComplexT *tau, int *info);



/** \private */
extern void
FORTRAN_NAME(zhetrd)(char *uplo, int *n, ComplexT *a,
	int *lda, double *d, double *e, ComplexT *tau,
	ComplexT *work, int *lwork, int *info);



/** \private */
extern void
FORTRAN_NAME(zhseqr)(char *job, char *compz, int *n, int *ilo,
	 int *ihi, ComplexT *h, int *ldh, ComplexT *w,
	ComplexT *z, int *ldz, ComplexT *work, int *lwork,
	 int *info);



/** \private */
extern void
FORTRAN_NAME(zlabrd)(int *m, int *n, int *nb,
	ComplexT *a, int *lda, double *d, double *e,
	ComplexT *tauq, ComplexT *taup, ComplexT *x, int *
	ldx, ComplexT *y, int *ldy);



/** \private */
extern void
FORTRAN_NAME(zlacgv)(int *n, ComplexT *x, int *incx);



/** \private */
extern void
FORTRAN_NAME(zlacpy)(char *uplo, int *m, int *n,
	ComplexT *a, int *lda, ComplexT *b, int *ldb);



/** \private */
extern void
FORTRAN_NAME(zlahqr)(int *wantt, int *wantz, int *n,
	int *ilo, int *ihi, ComplexT *h, int *ldh,
	ComplexT *w, int *iloz, int *ihiz, ComplexT *z,
	int *ldz, int *info);



/** \private */
extern void
FORTRAN_NAME(zlahrd)(int *n, int *k, int *nb,
	ComplexT *a, int *lda, ComplexT *tau, ComplexT *t,
	int *ldt, ComplexT *y, int *ldy);



/** \private */
extern double
FORTRAN_NAME(zlange)(char *norm, int *m, int *n, ComplexT *a, int *lda,
		 double *work);



/** \private */
extern double
FORTRAN_NAME(zlanhe)(char *norm,  char *uplo, int *n, ComplexT *a,
		 int *lda, double *work);



/** \private */
extern double
FORTRAN_NAME(zlanhs)(char *norm, int *n, ComplexT *a, int *lda, double *work);




/** \private */
extern void
FORTRAN_NAME(zlaqp2)(int *m, int *n, int *offset,
	ComplexT *a, int *lda, int *jpvt, ComplexT *tau,
	double *vn1, double *vn2, ComplexT *work);



/** \private */
extern void
FORTRAN_NAME(zlaqps)(int *m, int *n, int *offset, int
	*nb, int *kb, ComplexT *a, int *lda, int *jpvt,
	ComplexT *tau, double *vn1, double *vn2, ComplexT *
	auxv, ComplexT *f, int *ldf);



/** \private */
extern void
FORTRAN_NAME(zlarf)(char *side, int *m, int *n, ComplexT
	*v, int *incv, ComplexT *tau, ComplexT *c, int *
	ldc, ComplexT *work);



/** \private */
extern void
FORTRAN_NAME(zlarfb)(char *side, char *trans, char *direct, char *
	storev, int *m, int *n, int *k, ComplexT *v, int
	*ldv, ComplexT *t, int *ldt, ComplexT *c, int *
	ldc, ComplexT *work, int *ldwork);



/** \private */
extern void
FORTRAN_NAME(zlarfg)(int *n, ComplexT *alpha, ComplexT *
	x, int *incx, ComplexT *tau);



/** \private */
extern void
FORTRAN_NAME(zlarft)(char *direct, char *storev, int *n, int *
	k, ComplexT *v, int *ldv, ComplexT *tau, ComplexT *
	t, int *ldt);



/** \private */
extern void
FORTRAN_NAME(zlarfx)(char *side, int *m, int *n,
	ComplexT *v, ComplexT *tau, ComplexT *c, int *
	ldc, ComplexT *work);



/** \private */
extern void
FORTRAN_NAME(zlascl)(char *type, int *kl, int *ku,
	double *cfrom, double *cto, int *m, int *n,
	ComplexT *a, int *lda, int *info);



/** \private */
extern void
FORTRAN_NAME(zlaset)(char *uplo, int *m, int *n,
	ComplexT *alpha, ComplexT *beta, ComplexT *a, int *
	lda);



/** \private */
extern void
FORTRAN_NAME(zlasr)(char *side, char *pivot, char *direct, int *m,
	 int *n, double *c, double *s, ComplexT *a,
	int *lda);



/** \private */
extern void
FORTRAN_NAME(zlassq)(int *n, ComplexT *x, int *incx,
	double *scale, double *sumsq);



/** \private */
extern void
FORTRAN_NAME(zlaswp)(int *n, ComplexT *a, int *lda,
	int *k1, int *k2, int *ipiv, int *incx);



/** \private */
extern void
FORTRAN_NAME(zlatrd)(char *uplo, int *n, int *nb,
	ComplexT *a, int *lda, double *e, ComplexT *tau,
	ComplexT *w, int *ldw);



/** \private */
extern void
FORTRAN_NAME(zlatrs)(char *uplo, char *trans, char *diag, char *
	normin, int *n, ComplexT *a, int *lda, ComplexT *x,
	double *scale, double *cnorm, int *info);



/** \private */
extern void
FORTRAN_NAME(zsteqr)(char *compz, int *n, double *d,
	double *e, ComplexT *z, int *ldz, double *work,
	int *info);

/* ZTRCON estimates the reciprocal of the condition number of a
 * triangular matrix A, in either the 1-norm or the infinity-norm.
 */


/** \private */
extern void
FORTRAN_NAME(ztrcon)(const char *norm, const char *uplo, const char *diag,
                 const int *n, const ComplexT *a, const int *lda,
		 double *rcond, ComplexT *work, double *rwork, int *info);



/** \private */
extern void
FORTRAN_NAME(ztrevc)(char *side, char *howmny, int *select,
	int *n, ComplexT *t, int *ldt, ComplexT *vl,
	int *ldvl, ComplexT *vr, int *ldvr, int *mm, int
	*m, ComplexT *work, double *rwork, int *info);



/** \private */
extern void
FORTRAN_NAME(zung2l)(int *m, int *n, int *k,
	ComplexT *a, int *lda, ComplexT *tau, ComplexT *
	work, int *info);



/** \private */
extern void
FORTRAN_NAME(zung2r)(int *m, int *n, int *k,
	ComplexT *a, int *lda, ComplexT *tau, ComplexT *
	work, int *info);



/** \private */
extern void
FORTRAN_NAME(zungbr)(char *vect, int *m, int *n, int *k,
	ComplexT *a, int *lda, ComplexT *tau, ComplexT *
	work, int *lwork, int *info);



/** \private */
extern void
FORTRAN_NAME(zunghr)(int *n, int *ilo, int *ihi,
	ComplexT *a, int *lda, ComplexT *tau, ComplexT *
	work, int *lwork, int *info);



/** \private */
extern void
FORTRAN_NAME(zungl2)(int *m, int *n, int *k,
	ComplexT *a, int *lda, ComplexT *tau, ComplexT *
	work, int *info);



/** \private */
extern void
FORTRAN_NAME(zunglq)(int *m, int *n, int *k,
	ComplexT *a, int *lda, ComplexT *tau, ComplexT *
	work, int *lwork, int *info);



/** \private */
extern void
FORTRAN_NAME(zungql)(int *m, int *n, int *k,
	ComplexT *a, int *lda, ComplexT *tau, ComplexT *
	work, int *lwork, int *info);



/** \private */
extern void
FORTRAN_NAME(zungqr)(int *m, int *n, int *k,
	ComplexT *a, int *lda, ComplexT *tau, ComplexT *
	work, int *lwork, int *info);



/** \private */
extern void
FORTRAN_NAME(zungr2)(int *m, int *n, int *k,
	ComplexT *a, int *lda, ComplexT *tau, ComplexT *
	work, int *info);



/** \private */
extern void
FORTRAN_NAME(zungrq)(int *m, int *n, int *k,
	ComplexT *a, int *lda, ComplexT *tau, ComplexT *
	work, int *lwork, int *info);



/** \private */
extern void
FORTRAN_NAME(zungtr)(char *uplo, int *n, ComplexT *a,
	int *lda, ComplexT *tau, ComplexT *work, int *lwork,
	 int *info);



/** \private */
extern void
FORTRAN_NAME(zunm2r)(char *side, char *trans, int *m, int *n,
	int *k, ComplexT *a, int *lda, ComplexT *tau,
	ComplexT *c, int *ldc, ComplexT *work, int *info);



/** \private */
extern void
FORTRAN_NAME(zunmbr)(char *vect, char *side, char *trans, int *m,
	int *n, int *k, ComplexT *a, int *lda, ComplexT
	*tau, ComplexT *c, int *ldc, ComplexT *work, int *
	lwork, int *info);



/** \private */
extern void
FORTRAN_NAME(zunml2)(char *side, char *trans, int *m, int *n,
	int *k, ComplexT *a, int *lda, ComplexT *tau,
	ComplexT *c, int *ldc, ComplexT *work, int *info);



/** \private */
extern void
FORTRAN_NAME(zunmlq)(char *side, char *trans, int *m, int *n,
	int *k, ComplexT *a, int *lda, ComplexT *tau,
	ComplexT *c, int *ldc, ComplexT *work, int *lwork,
	 int *info);

#ifdef	__cplusplus
}
#endif

#endif /* R_LAPACK_H */
