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


/* Created:  Sat Jul 6 14:28:32 EDT 1996 */

#include "dcsrvml.h"



void CSR_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  double t;
  const double *pval;
  double *pc=c;                                       
  int i,j,jb,je;
  i=k;
  b-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                 

    pval = val;
    for (i=0;i!=m;i++) {
      t = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        t +=  b[indx[j]] * (*pval++);
      c[i] += t;
    }
}

void CSRsymm_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
    for (i=0;i!=m;i++) *pc++ = 0;                 

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
}



void CSRskew_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
    for (i=0;i!=m;i++) *pc++ = 0;                 
  
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
}



void CSR_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  double t;
  const double *pval;
  double *pc=c;                                       
  int i,j,jb,je;
  i=k;
  b-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pc++ = 0;                 

    pval = val;
    for (i=0;i!=m;i++) {
      t = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        t += alpha * b[indx[j]] * (*pval++);
      c[i] += t;
    }
}

void CSRsymm_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
    for (i=0;i!=m;i++) *pc++ = 0;                 

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
}



void CSRskew_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
    for (i=0;i!=m;i++) *pc++ = 0;                 
  
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
}



void CSR_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  double t;
  const double *pval;
  int i,j,jb,je;
  i=k;
  b-=ind_base;
  indx-=ind_base;


    pval = val;
    for (i=0;i!=m;i++) {
      t = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        t +=  b[indx[j]] * (*pval++);
      c[i] += t;
    }
}

void CSRsymm_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  const double *pval;
  int j;
  j=m;
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      

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
}



void CSRskew_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  const double *pval;
  int j;
  j=m;
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
  
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
}



void CSR_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  double t;
  const double *pval;
  int i,j,jb,je;
  i=k;
  b-=ind_base;
  indx-=ind_base;


    pval = val;
    for (i=0;i!=m;i++) {
      t = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        t += alpha * b[indx[j]] * (*pval++);
      c[i] += t;
    }
}

void CSRsymm_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  const double *pval;
  int j;
  j=m;
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      

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
}



void CSRskew_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base)
{
  const double *pval;
  int j;
  j=m;
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
  
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
}



void CSR_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base)
{
  double t;
  const double *pval;
  double *pc=c;                                       
  int i,j,jb,je;
  j=k;
  b-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                 

    pval = val;
    for (i=0;i!=m;i++) {
      t = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        t +=  b[indx[j]] * (*pval++);
      c[i] += t;
    }
}

void CSRsymm_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
    for (i=0;i!=m;i++) *pc++ *= beta;                 

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
}



void CSRskew_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
    for (i=0;i!=m;i++) *pc++ *= beta;                 
  
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
}



void CSR_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base)
{
  double t;
  const double *pval;
  double *pc=c;                                       
  int i,j,jb,je;
  i=k;
  b-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pc++ *= beta;                 

    pval = val;
    for (i=0;i!=m;i++) {
      t = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j!=je;j++)
        t += alpha * b[indx[j]] * (*pval++);
      c[i] += t;
    }
}

void CSRsymm_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
    for (i=0;i!=m;i++) *pc++ *= beta;                 

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
}



void CSRskew_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base)
{
  const double *pval;
  double *pc=c;                                       
  int i,j;
  int jj;
  int rpntrb, rpntre;
  int index;
  indx-=ind_base;
      
    for (i=0;i!=m;i++) *pc++ *= beta;                 
  
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
}
