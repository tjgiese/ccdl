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


/* Created:  Sat Jul 6 14:28:31 EDT 1996 */

#include "dcsrmml.h"



void CSR_MatMult_CAB_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base)
{
  double t;
  const double *pval;
  double *pc=c;                                       
  int i,j,jb,je;
  i=k;
  int l;                                               
  b-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ = 0;                 

  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (i=0;i!=m;i++) {
      t = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        t +=  b[indx[j]] * (*pval++);
      c[i] += t;
    }
    c += ldc; b += ldb;                                
  }                                                    
}

void CSRsymm_MatMult_CAB_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int l;                                               
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ = 0;                 

  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          c[j] +=  b[j] * (*pval++);
          continue;
        }
        c-=ind_base;
        c[index] +=  b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] +=  b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                
  }                                                    
}



void CSRskew_MatMult_CAB_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int l;                                               
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ = 0;                 
  
  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          pval++;
          continue;
        }
        c-=ind_base;
        c[index] -=  b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] +=  b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                 
  }                                                     
}



void CSR_MatMult_CaAB_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base)
{
  double t;
  const double *pval;
  double *pc=c;                                       
  int i,j,jb,je;
  i=k;
  int l;                                               
  b-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ = 0;                 

  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (i=0;i!=m;i++) {
      t = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        t += alpha * b[indx[j]] * (*pval++);
      c[i] += t;
    }
    c += ldc; b += ldb;                                
  }                                                    
}

void CSRsymm_MatMult_CaAB_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int l;                                               
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ = 0;                 

  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          c[j] += alpha * b[j] * (*pval++);
          continue;
        }
        c-=ind_base;
        c[index] += alpha * b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] += alpha * b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                
  }                                                    
}



void CSRskew_MatMult_CaAB_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int l;                                               
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ = 0;                 
  
  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          pval++;
          continue;
        }
        c-=ind_base;
        c[index] -= alpha * b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] += alpha * b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                 
  }                                                     
}



void CSR_MatMult_CABC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base)
{
  double t;
  const double *pval;
  int i,j,jb,je;
  i=k;
  int l;                                               
  b-=ind_base;
  indx-=ind_base;


  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (i=0;i!=m;i++) {
      t = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        t +=  b[indx[j]] * (*pval++);
      c[i] += t;
    }
    c += ldc; b += ldb;                                
  }                                                    
}

void CSRsymm_MatMult_CABC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  int j;
  j=m;
  int l;                                               
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      

  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          c[j] +=  b[j] * (*pval++);
          continue;
        }
        c-=ind_base;
        c[index] +=  b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] +=  b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                
  }                                                    
}



void CSRskew_MatMult_CABC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  int j;
  j=m;
  int l;                                               
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
  
  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          pval++;
          continue;
        }
        c-=ind_base;
        c[index] -=  b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] +=  b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                 
  }                                                     
}



void CSR_MatMult_CaABC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base)
{
  double t;
  const double *pval;
  int i,j,jb,je;
  i=k;
  int l;                                               
  b-=ind_base;
  indx-=ind_base;


  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (i=0;i!=m;i++) {
      t = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        t += alpha * b[indx[j]] * (*pval++);
      c[i] += t;
    }
    c += ldc; b += ldb;                                
  }                                                    
}

void CSRsymm_MatMult_CaABC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  int j;
  j=m;
  int l;                                               
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      

  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          c[j] += alpha * b[j] * (*pval++);
          continue;
        }
        c-=ind_base;
        c[index] += alpha * b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] += alpha * b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                
  }                                                    
}



void CSRskew_MatMult_CaABC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  int j;
  j=m;
  int l;                                               
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
  
  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          pval++;
          continue;
        }
        c-=ind_base;
        c[index] -= alpha * b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] += alpha * b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                 
  }                                                     
}



void CSR_MatMult_CABbC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base)
{
  double t;
  const double *pval;
  double *pc=c;                                       
  int i,j,jb,je;
  i=k;
  int l;                                               
  b-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ *= beta;                 

  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (i=0;i!=m;i++) {
      t = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        t +=  b[indx[j]] * (*pval++);
      c[i] += t;
    }
    c += ldc; b += ldb;                                
  }                                                    
}

void CSRsymm_MatMult_CABbC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int l;                                               
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ *= beta;                 

  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          c[j] +=  b[j] * (*pval++);
          continue;
        }
        c-=ind_base;
        c[index] +=  b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] +=  b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                
  }                                                    
}



void CSRskew_MatMult_CABbC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int l;                                               
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ *= beta;                 
  
  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          pval++;
          continue;
        }
        c-=ind_base;
        c[index] -=  b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] +=  b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                 
  }                                                     
}



void CSR_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base)
{
  double t;
  const double *pval;
  double *pc=c;                                       
  int i,j,jb,je;
  i=k;
  int l;                                               
  b-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ *= beta;                 

  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (i=0;i!=m;i++) {
      t = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        t += alpha * b[indx[j]] * (*pval++);
      c[i] += t;
    }
    c += ldc; b += ldb;                                
  }                                                    
}

void CSRsymm_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int l;                                               
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ *= beta;                 

  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          c[j] += alpha * b[j] * (*pval++);
          continue;
        }
        c-=ind_base;
        c[index] += alpha * b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] += alpha * b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                
  }                                                    
}



void CSRskew_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int l;                                               
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
  for (l=0;l!=n;l++)                           
    for (i=0;i!=m;i++) *pc++ *= beta;                 
  
  for (l=0;l!=n;l++) {                                 
    pval = val;
    for (j=0;j!=k;j++){
      rpntrb = pntrb[j];
      rpntre = pntre[j];
      for (jj=rpntrb;jj!=rpntre;jj++) {
        index = indx[jj];
        if ( index == j+ind_base ) {
          pval++;
          continue;
        }
        c-=ind_base;
        c[index] -= alpha * b[j] * (*pval);
        c+=ind_base;
        b-=ind_base;
        c[j] += alpha * b[index] * (*pval++);
        b+=ind_base;
      }
    }
    c += ldc; b += ldb;                                 
  }                                                     
}
