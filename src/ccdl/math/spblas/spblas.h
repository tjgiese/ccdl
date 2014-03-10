#ifndef SPBLAS_H
#define SPBLAS_H

#ifdef __cplusplus
extern "C" {
#endif

/*+++++++++++++++++++++++++++++++++++*/
/*  NIST SPARSE BLAS                 */
/*  C/C++ Toolkit layer interfaces:  */
/*+++++++++++++++++++++++++++++++++++*/

/* void  dcoomm( */
/*              const int transa, const int m, const int n, const int k, */
/*              const double alpha, const int descra[], const double val[], */
/*              const int indx[], const int jndx[], const int nnz, */
/*              const double b[], const int ldb, */
/*              const double beta, double c[], const int ldc, */
/*              double work[], const int lwork); */

/* void dcoosm(const int transa, const int m, const int n, */
/*    const int unitd, const double dv[], const double alpha, */
/*    const int descra[], const double val[], */
/*    const int indx[], const int jndx[], const int nnz, const double b[], const int ldb, */
/*    const double beta, double c[], const int ldc, */
/*    double work[], const int lwork); */

void  dcscmm(
             const int transa, const int m, const int n, const int k,
             const double alpha, const int descra[], const double val[],
             const int indx[], const int pntrb[], const int pntre[],
             const double b[], const int ldb,
             const double beta, double c[], const int ldc);
/* void  dcscsm( */
/*              const int transa, const int m, const int n,  */
/*              const int unitd, const double dv[],  */
/*              const double alpha, const int descra[], const double val[], */
/*              const int indx[], const int pntrb[], const int pntre[], */
/*              const double b[], const int ldb, */
/*              const double beta, double c[], const int ldc, */
/*              double work[], const int lwork); */
void  dcsrmm(
             const int transa, const int m, const int n, const int k,
             const double alpha, const int descra[], const double val[],
             const int indx[], const int pntrb[], const int pntre[],
             const double b[], const int ldb,
             const double beta, double c[], const int ldc);
void  dcsrsm(
             const int transa, const int m, const int n,
             const int unitd, const double dv[],
             const double alpha, const int descra[], const double val[],
             const int indx[], const int pntrb[], const int pntre[],
             const double b[], const int ldb,
             const double beta, double c[], const int ldc,
             double work[], const int lwork);

/* void  dbcomm( */
/*              const int transa, const int mb, const int n, const int kb, */
/*              const double alpha, const int descra[], const double val[], */
/*              const int bindx[], const int bjndx[], const int bnnz, */
/*              const int lb, const double b[], const int ldb, */
/*              const double beta, double c[], const int ldc, */
/*              double work[], const int lwork); */
/* void  dbscmm( */
/*              const int transa, const int mb, const int n, const int kb, */
/*              const double alpha, const int descra[], const double val[], */
/*              const int bindx[], const int bpntrb[], const int bpntre[], */
/*              const int lb, const double b[], const int ldb, */
/*              const double beta, double c[], const int ldc, */
/*              double work[], const int lwork); */
/* void  dbscsm( */
/*              const int transa, const int mb, const int n,  */
/*              const int unitd, const double dv[],  */
/*              const double alpha, const int descra[], const double val[], */
/*              const int bindx[], const int bpntrb[], const int bpntre[], */
/*              const int lb, const double b[], const int ldb, */
/*              const double beta, double c[], const int ldc, */
/*              double work[], const int lwork); */

/* void  dbsrmm( */
/*              const int transa, const int mb, const int n, const int kb, */
/*              const double alpha, const int descra[], const double val[], */
/*              const int bindx[], const int bpntrb[], const int bpntre[], */
/*              const int lb, const double b[], const int ldb, */
/*              const double beta, double c[], const int ldc, */
/*              double work[], const int lwork); */
/* void  dbsrsm( */
/*              const int transa, const int mb, const int n,  */
/*              const int unitd, const double dv[],  */
/*              const double alpha, const int descra[], const double val[], */
/*              const int bindx[], const int bpntrb[], const int bpntre[], */
/*              const int lb, const double b[], const int ldb, */
/*              const double beta, double c[], const int ldc, */
/*              double work[], const int lwork); */

/* void  dvbrmm( */
/*              const int transa, const int mb, const int n, const int kb, */
/*              const double alpha, const int descra[], const double val[], */
/*              const int indx[], const int bindx[],  */
/*              const int rpntr[], const int cpntr[],  */
/*              const int bpntrb[], const int bpntre[], */
/*              const double b[], const int ldb, */
/*              const double beta, double c[], const int ldc, */
/*              double work[], const int lwork); */
/* void  dvbrsm( */
/*              const int transa, const int mb, const int n,  */
/*              const int unitd, const double dv[],  */
/*              const double alpha, const int descra[], const double val[], */
/*              const int indx[], const int bindx[], const int rpntr[], */
/*              const int cpntr[], const int bpntrb[], const int bpntre[], */
/*              const double b[], const int ldb, */
/*              const double beta, double c[], const int ldc, */
/*              double work[], const int lwork); */

/* /\*+++++++++++++++++++++++++++++++++++++*\/ */
/* /\*  NIST SPARSE BLAS                   *\/ */
/* /\*  Fortran Toolkit layer interfaces:  *\/ */
/* /\*+++++++++++++++++++++++++++++++++++++*\/ */

/* void dcoomm_( */
/*              const int *transa, const int *m, const int *n, const int *k, */
/*              const double *alpha, const int descra[], const double val[], */
/*              const int indx[], const int jndx[], const int *nnz, */
/*              const double b[], const int *ldb, */
/*              const double *beta, double c[], const int *ldc, */
/*              double work[], const int *lwork); */
/* void dcoosm_(const int *transa, const int *m, const int *n, */
/*              const int *unitd, const double dv[], const double *alpha, */
/*              const int descra[], const double val[], */
/*              const int indx[], const int jndx[], const int *nnz, const double b[], const int *ldb, */
/*              const double *beta, double c[], const int *ldc, */
/*              double work[], const int *lwork); */
/* void dcscmm_( */
/*              const int *transa, const int *m, const int *n, const int *k, */
/*              const double *alpha, const int descra[], const double val[], */
/*              const int indx[], const int pntrb[], const int pntre[], */
/*              const double b[], const int *ldb, */
/*              const double *beta, double c[], const int *ldc, */
/*              double work[], const int *lwork); */
/* void dcscsm_( */
/*              const int *transa, const int *m, const int *n,  */
/*              const int *unitd, const double dv[],  */
/*              const double *alpha, const int descra[], const double val[], */
/*              const int indx[], const int pntrb[], const int pntre[], */
/*              const double b[], const int *ldb, */
/*              const double *beta, double c[], const int *ldc, */
/*              double work[], const int *lwork); */
/* void dcsrmm_( */
/*              const int *transa, const int *m, const int *n, const int *k, */
/*              const double *alpha, const int descra[], const double val[], */
/*              const int indx[], const int pntrb[], const int pntre[], */
/*              const double b[], const int *ldb, */
/*              const double *beta, double c[], const int *ldc, */
/*              double work[], const int *lwork); */
/* void dcsrsm_( */
/*              const int *transa, const int *m, const int *n,  */
/*              const int *unitd, const double dv[],  */
/*              const double *alpha, const int descra[], const double val[], */
/*              const int indx[], const int pntrb[], const int pntre[], */
/*              const double b[], const int *ldb, */
/*              const double *beta, double c[], const int *ldc, */
/*              double work[], const int *lwork); */
/* void dbcomm_( */
/*              const int *transa, const int *mb, const int *n, const int *kb, */
/*              const double *alpha, const int descra[], const double val[], */
/*              const int bindx[], const int bjndx[], const int *bnnz, */
/*              const int *lb, const double b[], const int *ldb, */
/*              const double *beta, double c[], const int *ldc, */
/*              double work[], const int *lwork); */
/* void dbscmm_( */
/*              const int *transa, const int *mb, const int *n, const int *kb, */
/*              const double *alpha, const int descra[], const double val[], */
/*              const int bindx[], const int bpntrb[], const int bpntre[], */
/*              const int *lb, const double b[], const int *ldb, */
/*              const double *beta, double c[], const int *ldc, */
/*              double work[], const int *lwork); */
/* void dbscsm_( */
/*              const int *transa, const int *mb, const int *n,  */
/*              const int *unitd, const double dv[],  */
/*              const double *alpha, const int descra[], const double val[], */
/*              const int bindx[], const int bpntrb[], const int bpntre[], */
/*              const int *lb, const double b[], const int *ldb, */
/*              const double *beta, double c[], const int *ldc, */
/*              double work[], const int *lwork); */
/* void dbsrmm_( */
/*              const int *transa, const int *mb, const int *n, const int *kb, */
/*              const double *alpha, const int descra[], const double val[], */
/*              const int bindx[], const int bpntrb[], const int bpntre[], */
/*              const int *lb, const double b[], const int *ldb, */
/*              const double *beta, double c[], const int *ldc, */
/*              double work[], const int *lwork); */
/* void dbsrsm_( */
/*              const int *transa, const int *mb, const int *n,  */
/*              const int *unitd, const double dv[],  */
/*              const double *alpha, const int descra[], const double val[], */
/*              const int bindx[], const int bpntrb[], const int bpntre[], */
/*              const int *lb, const double b[], const int *ldb, */
/*              const double *beta, double c[], const int *ldc, */
/*              double work[], const int *lwork); */
/* void dvbrmm_( */
/*              const int *transa, const int *mb, const int *n, const int *kb, */
/*              const double *alpha, const int descra[], const double val[], */
/*              const int indx[], const int bindx[], const int rpntr[], */
/*              const int cpntr[], const int bpntrb[], const int bpntre[], */
/*              const double b[], const int *ldb, */
/*              const double *beta, double c[], const int *ldc, */
/*              double work[], const int *lwork); */
/* void dvbrsm_( */
/*              const int *transa, const int *mb, const int *n,  */
/*              const int *unitd, const double dv[],  */
/*              const double *alpha, const int descra[], const double val[], */
/*              const int indx[], const int bindx[], const int rpntr[], */
/*              const int cpntr[], const int bpntrb[], const int bpntre[], */
/*              const double b[], const int *ldb, */
/*              const double *beta, double c[], const int *ldc, */
/*              double work[], const int *lwork); */

/*+++++++++++++++++++*/
/* Utility routines: */
/*+++++++++++++++++++*/
void ScaleArray_double(int m, int n, double *c, int ldc, const double beta);

/* void ScaleArray_float(int m, int n, float *c, int ldc, const float beta); */

/* void ScaleArray_complex(int m, int n, float *c, int ldc, const float *beta); */

/* void ScaleArray_Z(int m, int n, double *c, int ldc, const double *beta); */


#ifdef __cplusplus
}
#endif

#endif
