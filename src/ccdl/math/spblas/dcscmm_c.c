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
#include "dcscmml.h"
#include "dcscvml.h"
#include "dcsrmml.h"
#include "dcsrvml.h"


/* Sparse BLAS Toolkit interface routine: */
void  dcscmm(
             const int transa, const int m, const int n, const int k,
             const double alpha, const int descra[], const double val[],
             const int indx[], const int pntrb[], const int pntre[],
             const double b[], const int ldb,
             const double beta, double c[], const int ldc )
{
/* ------------ begin interface description ------------
   Toolkit interface:
   dcscmm -- compressed sparse column format matrix-matrix multiply
  
   C <- alpha A B + beta C
  
   Arguments:
  
   int transa	Indicates how to operate with the sparse matrix
  		0 : operate with matrix
  		1 : operate with transpose matrix
  
   int m	Number of rows in matrix c
  
   int n	Number of columns in matrix c
  
   int k	Number of rows in matrix b
  
   double alpha Scalar parameter
  
   double beta  Scalar parameter
  
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
  
   int *indx    integer array of length nnz containing row indices
  
   int *pntrb   integer array of length k such that pntrb(j)-pntrb(1)
                points to location in val of the first nonzero element in column j
  
   int *pntre   integer array of length k such that pntre(j)-pntrb(1)
                points to location in val of the last nonzero element in column j
  
   double *b	rectangular array with first dimension ldb
  
   double *c	rectangular array with first dimension ldc
  
   double *work	scratch array of length lwork.  lwork should be at least
  		max(m,n)
  
   ------------ end interface description --------------*/
int ind_base = descra[3];

if (alpha == 0.0) {
   ScaleArray_double(m, n, c, ldc, beta);
   return;
} 
switch (descra[0]) {
case 1: /* Symmetric */
case 2: /* Hermitian (for real same as symmetric) */
  if ( m != k ) {
    printf("In dcscmm: inconsistant dimensions for a symmetric matrix");
    printf("m = %d  k = %d\nExiting...\n",m,k);
    exit(-1); 
  }
  switch (descra[1]) {
  case 1:  /* Lower triangular stored */
    switch (n) {
    case 1: /* Vec Mult */
      if (alpha == 1) {
        if (beta == 1) {
          CSCsymm_VecMult_CABC_double(m, k, 
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          CSCsymm_VecMult_CAB_double(m, k, 
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          CSCsymm_VecMult_CABbC_double(m, k, 
                                  val, indx, pntrb, pntre, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          CSCsymm_VecMult_CaABC_double(m, k, alpha,
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          CSCsymm_VecMult_CaAB_double(m, k, alpha,
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          CSCsymm_VecMult_CaABbC_double(m, k, alpha,
                                  val, indx, pntrb, pntre, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default: /* n is greater than 1 -- doing Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          CSCsymm_MatMult_CABC_double(m, n, k, 
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          CSCsymm_MatMult_CAB_double(m, n, k, 
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          CSCsymm_MatMult_CABbC_double(m, n, k, 
                                  val, indx, pntrb, pntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          CSCsymm_MatMult_CaABC_double(m, n, k, alpha,
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          CSCsymm_MatMult_CaAB_double(m, n, k, alpha,
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          CSCsymm_MatMult_CaABbC_double(m, n, k, alpha,
                                  val, indx, pntrb, pntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    } /* end switch on n */
    break;
 case 2: /* Upper triangular part stored */
    switch (n) {
    case 1: /* Vec Mult */
      if (alpha == 1) {
        if (beta == 1) {
          CSRsymm_VecMult_CABC_double(m, k, 
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          CSRsymm_VecMult_CAB_double(m, k, 
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          CSRsymm_VecMult_CABbC_double(m, k, 
                                  val, indx, pntrb, pntre, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          CSRsymm_VecMult_CaABC_double(m, k, alpha,
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          CSRsymm_VecMult_CaAB_double(m, k, alpha,
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          CSRsymm_VecMult_CaABbC_double(m, k, alpha,
                                  val, indx, pntrb, pntre, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default: /* Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          CSRsymm_MatMult_CABC_double(m, n, k, 
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          CSRsymm_MatMult_CAB_double(m, n, k, 
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          CSRsymm_MatMult_CABbC_double(m, n, k, 
                                  val, indx, pntrb, pntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          CSRsymm_MatMult_CaABC_double(m, n, k, alpha,
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          CSRsymm_MatMult_CaAB_double(m, n, k, alpha,
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          CSRsymm_MatMult_CaABbC_double(m, n, k, alpha,
                                  val, indx, pntrb, pntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    } /* end switch on n */
    break;
  default:
      printf("Invalid argument descra[1] in dcscmm. Use 1 or 2. \n");
      break;
  } /* end switch on descra[1] */
  break;
case 4: /* Skew-symmetric */
  if ( m != k ) {
    printf("In dcscmm: inconsistant dimensions for a skew-symmetric matrix");
    printf("m = %d  k = %d\nExiting...\n",m,k);
    exit(-1); 
  }
  switch ( transa ) {
  case 0:
    switch (n) {
    case 1: /* Vec Mult */
      if (alpha == 1) {
        if (beta == 1) {
          CSCskew_VecMult_CABC_double(m, k, 
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          CSCskew_VecMult_CAB_double(m, k, 
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          CSCskew_VecMult_CABbC_double(m, k, 
                                  val, indx, pntrb, pntre, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          CSCskew_VecMult_CaABC_double(m, k, alpha,
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          CSCskew_VecMult_CaAB_double(m, k, alpha,
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          CSCskew_VecMult_CaABbC_double(m, k, alpha,
                                  val, indx, pntrb, pntre, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default: /* Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          CSCskew_MatMult_CABC_double(m, n, k, 
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          CSCskew_MatMult_CAB_double(m, n, k, 
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          CSCskew_MatMult_CABbC_double(m, n, k, 
                                  val, indx, pntrb, pntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          CSCskew_MatMult_CaABC_double(m, n, k, alpha,
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          CSCskew_MatMult_CaAB_double(m, n, k, alpha,
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          CSCskew_MatMult_CaABbC_double(m, n, k, alpha,
                                  val, indx, pntrb, pntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    } /* end switch on n */
    break;
  case 1:    /* call CSR to get proper sign  */
    switch (n) {
    case 1:  /* Vec Mult */
      if (alpha == 1) {
        if (beta == 1) {
          CSRskew_VecMult_CABC_double(k, m, 
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          CSRskew_VecMult_CAB_double(k, m, 
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          CSRskew_VecMult_CABbC_double(k, m, 
                                  val, indx, pntrb, pntre, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          CSRskew_VecMult_CaABC_double(k, m, alpha,
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          CSRskew_VecMult_CaAB_double(k, m, alpha,
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          CSRskew_VecMult_CaABbC_double(k, m, alpha,
                                  val, indx, pntrb, pntre, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default: /* Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          CSRskew_MatMult_CABC_double(k, n, m, 
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          CSRskew_MatMult_CAB_double(k, n, m, 
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          CSRskew_MatMult_CABbC_double(k, n, m, 
                                  val, indx, pntrb, pntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          CSRskew_MatMult_CaABC_double(k, n, m, alpha,
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          CSRskew_MatMult_CaAB_double(k, n, m, alpha,
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          CSRskew_MatMult_CaABbC_double(k, n, m, alpha,
                                  val, indx, pntrb, pntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    } /* end switch on n */
  } /* end switch on transa */
  break;
case 0: case 3: case 5:
  switch (transa) {
  case 0:
    switch (n) {
    case 1: /* Vec Mult */
      if (alpha == 1) {
        if (beta == 1) {
          CSC_VecMult_CABC_double(m, k, 
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          CSC_VecMult_CAB_double(m, k, 
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          CSC_VecMult_CABbC_double(m, k, 
                                  val, indx, pntrb, pntre, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          CSC_VecMult_CaABC_double(m, k, alpha,
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          CSC_VecMult_CaAB_double(m, k, alpha,
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          CSC_VecMult_CaABbC_double(m, k, alpha,
                                  val, indx, pntrb, pntre, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default: /* Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          CSC_MatMult_CABC_double(m, n, k, 
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          CSC_MatMult_CAB_double(m, n, k, 
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          CSC_MatMult_CABbC_double(m, n, k, 
                                  val, indx, pntrb, pntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          CSC_MatMult_CaABC_double(m, n, k, alpha,
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          CSC_MatMult_CaAB_double(m, n, k, alpha,
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          CSC_MatMult_CaABbC_double(m, n, k, alpha,
                                  val, indx, pntrb, pntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    } /* end switch on n */
    break;
  case 1: 
    switch (n) {
    case 1: /* Vec Mult */
      if (alpha == 1) {
        if (beta == 1) {
          CSR_VecMult_CABC_double(k, m, 
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          CSR_VecMult_CAB_double(k, m, 
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else  { /*  beta is general nonzero */
          CSR_VecMult_CABbC_double(k, m, 
                                  val, indx, pntrb, pntre, b, beta,
                                  c, ind_base);
        }
      } else { /* alpha is general nonzero */
        if (beta == 1) {
          CSR_VecMult_CaABC_double(k, m, alpha,
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else if (beta == 0) {
          CSR_VecMult_CaAB_double(k, m, alpha,
                                  val, indx, pntrb, pntre, b, 
                                  c, ind_base);
        } else { /*  beta is general nonzero */
          CSR_VecMult_CaABbC_double(k, m, alpha,
                                  val, indx, pntrb, pntre, b, beta,
                                  c, ind_base);
        }
      }
      break;
    default: /* Mat Mult */
      if (alpha == 1) {
        if (beta == 1) {
          CSR_MatMult_CABC_double(k, n, m, 
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          CSR_MatMult_CAB_double(k, n, m, 
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          CSR_MatMult_CABbC_double(k, n, m, 
                                  val, indx, pntrb, pntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      } else {   /* alpha is general nonzero */
        if (beta == 1) {
          CSR_MatMult_CaABC_double(k, n, m, alpha,
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else if (beta == 0) {
          CSR_MatMult_CaAB_double(k, n, m, alpha,
                                  val, indx, pntrb, pntre, b, ldb, 
                                  c, ldc, ind_base);
        } else { /*  beta is general nonzero */
          CSR_MatMult_CaABbC_double(k, n, m, alpha,
                                  val, indx, pntrb, pntre, b, ldb, beta,
                                  c, ldc, ind_base);
        }
      }
      break;
    } /* end switch on n */
    break;
  default: 
    printf("Invalid argument transa in dcscmm. Use 0 or 1. \n");
    break;
  } /* end switch on transa */
  break;
default:
  printf("Invalid argument descar[0] in dcscmm. Use 0-5. \n");
  break;
} /* end switch on descra[0] */

  
}
   
  
  
   
  
