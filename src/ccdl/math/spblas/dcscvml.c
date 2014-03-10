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


/* Created:  Sat Jul 6 14:30:42 EDT 1996 */

#include "dcscvml.h"



void CSC_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  val-=ind_base;
  indx-=ind_base;
  c-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                 

    for (i=0;i!=k;i++){
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        c[indx[j]] +=  b[i] * val[j];
    }
}

void CSCsymm_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                 

    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        c[j] +=  b[j] * val[jb];
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] +=  b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] +=  b[j] * val[i];
      }
      c+=ind_base;
    }
}

void CSCskew_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                 

    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] -=  b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] +=  b[j] * val[i];
      }
      c+=ind_base;
    }
}




void CSC_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  val-=ind_base;
  indx-=ind_base;
  c-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                 

    for (i=0;i!=k;i++){
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        c[indx[j]] += alpha * b[i] * val[j];
    }
}

void CSCsymm_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                 

    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        c[j] += alpha * b[j] * val[jb];
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] += alpha * b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] += alpha * b[j] * val[i];
      }
      c+=ind_base;
    }
}

void CSCskew_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                 

    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] -= alpha * b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] += alpha * b[j] * val[i];
      }
      c+=ind_base;
    }
}




void CSC_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  val-=ind_base;
  indx-=ind_base;
  c-=ind_base;


    for (i=0;i!=k;i++){
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        c[indx[j]] +=  b[i] * val[j];
    }
}

void CSCsymm_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  val-=ind_base;
  indx-=ind_base;


    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        c[j] +=  b[j] * val[jb];
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] +=  b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] +=  b[j] * val[i];
      }
      c+=ind_base;
    }
}

void CSCskew_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  val-=ind_base;
  indx-=ind_base;


    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] -=  b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] +=  b[j] * val[i];
      }
      c+=ind_base;
    }
}




void CSC_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  val-=ind_base;
  indx-=ind_base;
  c-=ind_base;


    for (i=0;i!=k;i++){
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        c[indx[j]] += alpha * b[i] * val[j];
    }
}

void CSCsymm_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  val-=ind_base;
  indx-=ind_base;


    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        c[j] += alpha * b[j] * val[jb];
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] += alpha * b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] += alpha * b[j] * val[i];
      }
      c+=ind_base;
    }
}

void CSCskew_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  val-=ind_base;
  indx-=ind_base;


    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] -= alpha * b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] += alpha * b[j] * val[i];
      }
      c+=ind_base;
    }
}




void CSC_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  double *pc=c;                                       
  val-=ind_base;
  indx-=ind_base;
  c-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                 

    for (i=0;i!=k;i++){
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        c[indx[j]] +=  b[i] * val[j];
    }
}

void CSCsymm_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  double *pc=c;                                       
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                 

    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        c[j] +=  b[j] * val[jb];
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] +=  b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] +=  b[j] * val[i];
      }
      c+=ind_base;
    }
}

void CSCskew_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  double *pc=c;                                       
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                 

    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] -=  b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] +=  b[j] * val[i];
      }
      c+=ind_base;
    }
}




void CSC_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  double *pc=c;                                       
  val-=ind_base;
  indx-=ind_base;
  c-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                 

    for (i=0;i!=k;i++){
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        c[indx[j]] += alpha * b[i] * val[j];
    }
}

void CSCsymm_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  double *pc=c;                                       
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                 

    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        c[j] += alpha * b[j] * val[jb];
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] += alpha * b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] += alpha * b[j] * val[i];
      }
      c+=ind_base;
    }
}

void CSCskew_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,   
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  double *pc=c;                                       
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                 

    for (j=0;j!=k;j++){
      jb = pntrb[j];
      je = pntre[j];
      if ( indx[jb] == j+ind_base ) {
        jb++;
      }
      b-=ind_base;
      for (i=jb;i!=je;i++){
        c[j] -= alpha * b[indx[i]] * val[i];
      }
      c-=ind_base;
      b+=ind_base;
      for (i=jb;i!=je;i++){
        c[indx[i]] += alpha * b[j] * val[i];
      }
      c+=ind_base;
    }
}

