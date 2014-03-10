/*-------------------------------------------------------|
|  NIST SPARSE BLAS v. 0.9 (Sat Jul 6 14:27:21 EDT 1996) |
|                                                        |
|  Authors:                                              |
|     Karin A. Remington and Roldan Pozo                 |
|     National Institute of Standards and Technology     |
|                                                        |
|  Based on the interface standard proposed in:          | 
|   "A Revised Proposal for a Sparse BLAS Toolkit" by    |
|    S. Carney and K. Wu -- University of Minnesota      |
|    M. Heroux and G. Li -- Cray Research                |  
|    R. Pozo and K.A. Remington -- NIST                  |
|                                                        |
|  Contact:                                              |
|     Karin A. Remington, email: kremington@nist.gov     |
--------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include "spblas.h"
#include "dcscmtsl.h"
#include "dcscvtsl.h"
#include "dcsrmtsl.h"
#include "dcsrvtsl.h"


/* Sparse BLAS Toolkit interface routine: */
void  dcsrsm(
             const int transa, const int m, const int n, 
             const int unitd, const double dv[], 
             const double alpha, const int descra[], const double val[],
             const int indx[], const int pntrb[], const int pntre[],
             const double b[], const int ldb,
             const double beta, double c[], const int ldc,
             double work[], const int lwork)
{
/* ------------ begin interface description ------------
   Toolkit interface:
   dcsrsm -- compressed sparse row format triangular solve
  
   C <- alpha D inv(A) B + beta C    C <- alpha D inv(A') B + beta C
   C <- alpha inv(A) D B + beta C    C <- alpha inv(A') D B + beta C
   
                                      ( ' indicates matrix transpose)
  
   Arguments:
  
   int transa	Indicates how to operate with the sparse matrix
  		0 : operate with matrix
  		1 : operate with transpose matrix
  
   int m	Number of rows in matrix c
  
   int n	Number of columns in matrix c
  
   int unitd	Type of scaling:
                        1 : Identity matrix (argument dv[] is ignored)
                        2 : Scale on left (row scaling)
                        3 : Scale on right (column scaling)
  
   double alpha	Scalar parameter
  
   double beta 	Scalar parameter
  
   int descra[]	Descriptor argument.  Nine element integer array
  		descra[0] matrix structure
  			0 : general
  			1 : symmetric
  			2 : Hermitian
  			3 : Triangular
  			4 : Skew(Anti-Symmetric
  			5 : Diagonal
  		descra[1] upper/lower triangular indicator
  			1 : lower
  			2 : upper
  		descra[2] main diagonal type
  			0 : non-unit
  			1 : unit
  		descra[3] Array base 
  			0 : C/C++ compatible
  			1 : Fortran compatible
  		descra[4] repeated indices?
  			0 : unknown
  			1 : no repeated indices
  
  
   double *val	scalar array of length nnz containing matrix entries
  
   int *indx    integer array of length nnz containing column indices
  
   int *pntrb   integer array of length k such that pntrb(j)-pntrb(1)
                points to location in val of the first nonzero element in row j
  
   int *pntre   integer array of length k such that pntre(j)-pntrb(1)
                points to location in val of the last nonzero element in row j
  
   double *b	rectangular array with first dimension ldb
  
   double *c	rectangular array with first dimension ldc
  
   double *work	scratch array of length lwork.  
                lwork should be at least (n*m)
  
   ------------ end interface description --------------*/
int ind_base = descra[3];

  if (lwork < m*n ){
    printf("Insufficient work space for dcsrsm.\n");
    printf("   lwork must be at least (n*m) = %d \n",m*n);
    return;
  }

if (alpha == 0.0) {
   ScaleArray_double(m, n, c, ldc, beta);
   return;
} 
switch (descra[0]) {
case 3:  /* Matrix MUST be triangular, of course */
  switch (transa) {
  case 0:
    switch  ( 10*descra[1]+descra[2] ) {
    case 10:  /* Lower triangular, non-unit diagonal */
       switch (n) {
       case 1: /* Vec Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLD_CABC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLD_CDABC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLD_CADBC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLD_CABmC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLD_CDABmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLD_CADBmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLD_CAB_double(m,  val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLD_CDAB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLD_CADB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLD_CABbC_double(m,  val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLD_CDABbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLD_CADBbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLD_CaABC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLD_CaDABC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLD_CaADBC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLD_CaABmC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLD_CaDABmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLD_CaADBmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLD_CaAB_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLD_CaDAB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLD_CaADB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLD_CaABbC_double(m,  alpha, val, indx, pntrb, pntre, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLD_CaDABbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLD_CaADBbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       default: /* Mat Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLD_CABC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLD_CDABC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLD_CADBC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLD_CABmC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLD_CDABmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLD_CADBmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLD_CAB_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLD_CDAB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLD_CADB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLD_CABbC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLD_CDABbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLD_CADBbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLD_CaABC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLD_CaDABC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLD_CaADBC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLD_CaABmC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLD_CaDABmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLD_CaADBmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLD_CaAB_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLD_CaDAB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLD_CaADB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLD_CaABbC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLD_CaDABbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLD_CaADBbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       } /* end switch on n */
       break;
    case 11:  /* Lower triangular, Unit diagonal */
       switch (n) {
       case 1: /* Vec Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLU_CABC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLU_CDABC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLU_CADBC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLU_CABmC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLU_CDABmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLU_CADBmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLU_CAB_double(m,  val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLU_CDAB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLU_CADB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLU_CABbC_double(m,  val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLU_CDABbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLU_CADBbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLU_CaABC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLU_CaDABC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLU_CaADBC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLU_CaABmC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLU_CaDABmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLU_CaADBmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLU_CaAB_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLU_CaDAB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLU_CaADB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvLU_CaABbC_double(m,  alpha, val, indx, pntrb, pntre, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvLU_CaDABbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvLU_CaADBbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       default: /* Mat Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLU_CABC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLU_CDABC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLU_CADBC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLU_CABmC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLU_CDABmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLU_CADBmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLU_CAB_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLU_CDAB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLU_CADB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLU_CABbC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLU_CDABbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLU_CADBbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLU_CaABC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLU_CaDABC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLU_CaADBC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLU_CaABmC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLU_CaDABmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLU_CaADBmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLU_CaAB_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLU_CaDAB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLU_CaADB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvLU_CaABbC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvLU_CaDABbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvLU_CaADBbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       } /* end switch on n */
       break;
    case 20:  /* Upper triangular, non-unit diagonal */
       switch (n) {
       case 1: /* Vec Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUD_CABC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUD_CDABC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUD_CADBC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUD_CABmC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUD_CDABmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUD_CADBmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUD_CAB_double(m,  val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUD_CDAB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUD_CADB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUD_CABbC_double(m,  val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUD_CDABbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUD_CADBbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUD_CaABC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUD_CaDABC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUD_CaADBC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUD_CaABmC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUD_CaDABmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUD_CaADBmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUD_CaAB_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUD_CaDAB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUD_CaADB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUD_CaABbC_double(m,  alpha, val, indx, pntrb, pntre, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUD_CaDABbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUD_CaADBbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       default: /* Mat Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUD_CABC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUD_CDABC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUD_CADBC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUD_CABmC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUD_CDABmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUD_CADBmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUD_CAB_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUD_CDAB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUD_CADB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUD_CABbC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUD_CDABbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUD_CADBbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUD_CaABC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUD_CaDABC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUD_CaADBC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUD_CaABmC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUD_CaDABmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUD_CaADBmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUD_CaAB_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUD_CaDAB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUD_CaADB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUD_CaABbC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUD_CaDABbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUD_CaADBbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       } /* end switch on n */
       break;
    case 21:  /* Upper triangular, Unit diagonal */
       switch (n) {
       case 1: /* Vec Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUU_CABC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUU_CDABC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUU_CADBC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUU_CABmC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUU_CDABmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUU_CADBmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUU_CAB_double(m,  val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUU_CDAB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUU_CADB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUU_CABbC_double(m,  val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUU_CDABbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUU_CADBbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUU_CaABC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUU_CaDABC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUU_CaADBC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUU_CaABmC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUU_CaDABmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUU_CaADBmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUU_CaAB_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUU_CaDAB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUU_CaADB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_VecTriangSlvUU_CaABbC_double(m,  alpha, val, indx, pntrb, pntre, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_VecTriangSlvUU_CaDABbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_VecTriangSlvUU_CaADBbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       default: /* Mat Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUU_CABC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUU_CDABC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUU_CADBC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUU_CABmC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUU_CDABmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUU_CADBmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUU_CAB_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUU_CDAB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUU_CADB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUU_CABbC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUU_CDABbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUU_CADBbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUU_CaABC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUU_CaDABC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUU_CaADBC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUU_CaABmC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUU_CaDABmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUU_CaADBmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUU_CaAB_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUU_CaDAB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUU_CaADB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSR_MatTriangSlvUU_CaABbC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSR_MatTriangSlvUU_CaDABbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSR_MatTriangSlvUU_CaADBbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       } /* end switch on n */
       break;
     default: /*  Invalid descra[1] or descra[2] */
       printf("Check values of descra[1] and descra[2] (descra(2) and descra(3) in Fortran.\n");
       printf("Valid values of descra[1] (descra(2)) are 1 and 2.\n");
       printf("Valid values of descra[2] (descra(3)) are 0 and 1.\n");
       break;
     } /* end switch on descra[1] and descra[2] */
     break;
  case 1:
    switch  ( 10*descra[1]+descra[2] ) {
    case 20:  /* Lower triangular, non-unit diagonal */
       switch (n) {
       case 1: /* Vec Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLD_CABC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLD_CDABC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLD_CADBC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLD_CABmC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLD_CDABmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLD_CADBmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLD_CAB_double(m,  val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLD_CDAB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLD_CADB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLD_CABbC_double(m,  val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLD_CDABbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLD_CADBbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLD_CaABC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLD_CaDABC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLD_CaADBC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLD_CaABmC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLD_CaDABmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLD_CaADBmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLD_CaAB_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLD_CaDAB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLD_CaADB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLD_CaABbC_double(m,  alpha, val, indx, pntrb, pntre, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLD_CaDABbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLD_CaADBbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       default: /* Mat Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLD_CABC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLD_CDABC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLD_CADBC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLD_CABmC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLD_CDABmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLD_CADBmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLD_CAB_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLD_CDAB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLD_CADB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLD_CABbC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLD_CDABbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLD_CADBbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLD_CaABC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLD_CaDABC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLD_CaADBC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLD_CaABmC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLD_CaDABmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLD_CaADBmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLD_CaAB_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLD_CaDAB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLD_CaADB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLD_CaABbC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLD_CaDABbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLD_CaADBbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       } /* end switch on n */
       break;
    case 21:  /* Lower triangular, Unit diagonal */
       switch (n) {
       case 1: /* Vec Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLU_CABC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLU_CDABC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLU_CADBC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLU_CABmC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLU_CDABmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLU_CADBmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLU_CAB_double(m,  val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLU_CDAB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLU_CADB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLU_CABbC_double(m,  val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLU_CDABbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLU_CADBbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLU_CaABC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLU_CaDABC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLU_CaADBC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLU_CaABmC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLU_CaDABmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLU_CaADBmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLU_CaAB_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLU_CaDAB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLU_CaADB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvLU_CaABbC_double(m,  alpha, val, indx, pntrb, pntre, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvLU_CaDABbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvLU_CaADBbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       default: /* Mat Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLU_CABC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLU_CDABC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLU_CADBC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLU_CABmC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLU_CDABmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLU_CADBmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLU_CAB_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLU_CDAB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLU_CADB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLU_CABbC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLU_CDABbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLU_CADBbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLU_CaABC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLU_CaDABC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLU_CaADBC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLU_CaABmC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLU_CaDABmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLU_CaADBmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLU_CaAB_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLU_CaDAB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLU_CaADB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvLU_CaABbC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvLU_CaDABbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvLU_CaADBbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       } /* end switch on n */
       break;
    case 10:  /* Upper triangular, non-unit diagonal */
       switch (n) {
       case 1: /* Vec Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUD_CABC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUD_CDABC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUD_CADBC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUD_CABmC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUD_CDABmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUD_CADBmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUD_CAB_double(m,  val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUD_CDAB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUD_CADB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUD_CABbC_double(m,  val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUD_CDABbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUD_CADBbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUD_CaABC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUD_CaDABC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUD_CaADBC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUD_CaABmC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUD_CaDABmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUD_CaADBmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUD_CaAB_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUD_CaDAB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUD_CaADB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUD_CaABbC_double(m,  alpha, val, indx, pntrb, pntre, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUD_CaDABbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUD_CaADBbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       default: /* Mat Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUD_CABC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUD_CDABC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUD_CADBC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUD_CABmC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUD_CDABmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUD_CADBmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUD_CAB_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUD_CDAB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUD_CADB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUD_CABbC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUD_CDABbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUD_CADBbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUD_CaABC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUD_CaDABC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUD_CaADBC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUD_CaABmC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUD_CaDABmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUD_CaADBmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUD_CaAB_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUD_CaDAB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUD_CaADB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUD_CaABbC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUD_CaDABbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUD_CaADBbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       } /* end switch on n */
       break;
    case 11:  /* Upper triangular, Unit diagonal */
       switch (n) {
       case 1: /* Vec Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUU_CABC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUU_CDABC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUU_CADBC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUU_CABmC_double(m,  val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUU_CDABmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUU_CADBmC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUU_CAB_double(m,  val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUU_CDAB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUU_CADB_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else  { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUU_CABbC_double(m,  val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUU_CDABbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUU_CADBbC_double(m,  dv, val, indx, pntrb, pntre, b, 
                                     beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else { /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUU_CaABC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUU_CaDABC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUU_CaADBC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUU_CaABmC_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUU_CaDABmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUU_CaADBmC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUU_CaAB_double(m,  alpha, val, indx, pntrb, pntre, b, 
                                     c, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUU_CaDAB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUU_CaADB_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, c, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_VecTriangSlvUU_CaABbC_double(m,  alpha, val, indx, pntrb, pntre, 
                                     b, beta, c, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_VecTriangSlvUU_CaDABbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_VecTriangSlvUU_CaADBbC_double(m,  dv, alpha, val, indx, pntrb, 
                                     pntre, b, beta, c, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       default: /* Mat Mult */
         if (alpha == 1) {
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUU_CABC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUU_CDABC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUU_CADBC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUU_CABmC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUU_CDABmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUU_CADBmC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUU_CAB_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUU_CDAB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUU_CADB_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUU_CABbC_double(m, n,  val, indx, pntrb, pntre, b, 
                                     ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUU_CDABbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUU_CADBbC_double(m, n,  dv, val, indx, pntrb, pntre, 
                                     b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         } else {   /* alpha is general nonzero */
           if (beta == 1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUU_CaABC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUU_CaDABC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUU_CaADBC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#if (0)
           } else if (beta == -1) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUU_CaABmC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUU_CaDABmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUU_CaADBmC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
#endif
           } else if (beta == 0) {
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUU_CaAB_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, c, ldc, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUU_CaDAB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUU_CaADB_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, c, ldc, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           } else { /*  beta is general nonzero */
             switch (unitd) {
             case 1: /* No scaling */
                CSC_MatTriangSlvUU_CaABbC_double(m, n,  alpha, val, indx, pntrb, 
                                     pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 2: /* Left scaling */
                CSC_MatTriangSlvUU_CaDABbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             case 3: /* Right scaling */
                CSC_MatTriangSlvUU_CaADBbC_double(m, n,  dv, alpha, val, indx, 
                                     pntrb, pntre, b, ldb, beta, c, ldc, work, ind_base);
                break;
             default:
                printf("Invalid argument unitd in dcsrsm. Use 1-3. \n");
                break;
             } /* end switch on unitd */
           }
         }
         break;
       } /* end switch on n */
       break;
     default: /*  Invalid descra[1] or descra[2] */
       printf("Check values of descra[1] and descra[2] (descra(2) and descra(3) in Fortran.\n");
       printf("Valid values of descra[1] (descra(2)) are 1 and 2.\n");
       printf("Valid values of descra[2] (descra(3)) are 0 and 1.\n");
       break;
     } /* end switch on descra[1] and descra[2] */
     break;
  default: 
     printf("Invalid argument transa in dcsrsm. Use 0 or 1. \n");
     break;
  } /* end switch on transa */
  break;
default:
  printf("Invalid argument descra[0] in dcsrsm. Must be 3 (triangular matrix). \n");
  break;
} /* end switch on descra[0] */

  
}
   
  
  
   
  
