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


/* Created:  Sat Jul 6 14:30:41 EDT 1996 */

#include "dcscmml.h"



void CSC_MatMult_CAB_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  int l;                                               
  val-=ind_base;
  indx-=ind_base;
  c-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ = 0;                 

  for (l=0;l!=n;l++) {                                 
    for (i=0;i!=k;i++){
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        c[indx[j]] +=  b[i] * val[j];
    }
    c += ldc; b += ldb;                                
  }                                                    
}

void CSCsymm_MatMult_CAB_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  int l;                                               
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ = 0;                 

  for (l=0;l!=n;l++) {                                 
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
    c += ldc; b += ldb;                                
  }                                                    
}

void CSCskew_MatMult_CAB_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  int l;                                               
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ = 0;                 

  for (l=0;l!=n;l++) {                                 
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
    c += ldc; b += ldb;                                
  }                                                    
}




void CSC_MatMult_CaAB_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  int l;                                               
  val-=ind_base;
  indx-=ind_base;
  c-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ = 0;                 

  for (l=0;l!=n;l++) {                                 
    for (i=0;i!=k;i++){
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        c[indx[j]] += alpha * b[i] * val[j];
    }
    c += ldc; b += ldb;                                
  }                                                    
}

void CSCsymm_MatMult_CaAB_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  int l;                                               
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ = 0;                 

  for (l=0;l!=n;l++) {                                 
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
    c += ldc; b += ldb;                                
  }                                                    
}

void CSCskew_MatMult_CaAB_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  int l;                                               
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ = 0;                 

  for (l=0;l!=n;l++) {                                 
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
    c += ldc; b += ldb;                                
  }                                                    
}




void CSC_MatMult_CABC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  int l;                                               
  val-=ind_base;
  indx-=ind_base;
  c-=ind_base;


  for (l=0;l!=n;l++) {                                 
    for (i=0;i!=k;i++){
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        c[indx[j]] +=  b[i] * val[j];
    }
    c += ldc; b += ldb;                                
  }                                                    
}

void CSCsymm_MatMult_CABC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  int l;                                               
  val-=ind_base;
  indx-=ind_base;


  for (l=0;l!=n;l++) {                                 
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
    c += ldc; b += ldb;                                
  }                                                    
}

void CSCskew_MatMult_CABC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  int l;                                               
  val-=ind_base;
  indx-=ind_base;


  for (l=0;l!=n;l++) {                                 
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
    c += ldc; b += ldb;                                
  }                                                    
}




void CSC_MatMult_CaABC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  int l;                                               
  val-=ind_base;
  indx-=ind_base;
  c-=ind_base;


  for (l=0;l!=n;l++) {                                 
    for (i=0;i!=k;i++){
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        c[indx[j]] += alpha * b[i] * val[j];
    }
    c += ldc; b += ldb;                                
  }                                                    
}

void CSCsymm_MatMult_CaABC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  int l;                                               
  val-=ind_base;
  indx-=ind_base;


  for (l=0;l!=n;l++) {                                 
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
    c += ldc; b += ldb;                                
  }                                                    
}

void CSCskew_MatMult_CaABC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  i=m;
  int l;                                               
  val-=ind_base;
  indx-=ind_base;


  for (l=0;l!=n;l++) {                                 
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
    c += ldc; b += ldb;                                
  }                                                    
}




void CSC_MatMult_CABbC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  int l;                                               
  val-=ind_base;
  indx-=ind_base;
  c-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ *= beta;                 

  for (l=0;l!=n;l++) {                                 
    for (i=0;i!=k;i++){
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        c[indx[j]] +=  b[i] * val[j];
    }
    c += ldc; b += ldb;                                
  }                                                    
}

void CSCsymm_MatMult_CABbC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  int l;                                               
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ *= beta;                 

  for (l=0;l!=n;l++) {                                 
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
    c += ldc; b += ldb;                                
  }                                                    
}

void CSCskew_MatMult_CABbC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  int l;                                               
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ *= beta;                 

  for (l=0;l!=n;l++) {                                 
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
    c += ldc; b += ldb;                                
  }                                                    
}




void CSC_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  int l;                                               
  val-=ind_base;
  indx-=ind_base;
  c-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ *= beta;                 

  for (l=0;l!=n;l++) {                                 
    for (i=0;i!=k;i++){
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        c[indx[j]] += alpha * b[i] * val[j];
    }
    c += ldc; b += ldb;                                
  }                                                    
}

void CSCsymm_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  int l;                                               
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ *= beta;                 

  for (l=0;l!=n;l++) {                                 
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
    c += ldc; b += ldb;                                
  }                                                    
}

void CSCskew_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base)
{
  int i,j,jb,je;
  double *pc=c;                                       
  int l;                                               
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ *= beta;                 

  for (l=0;l!=n;l++) {                                 
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
    c += ldc; b += ldb;                                
  }                                                    
}

