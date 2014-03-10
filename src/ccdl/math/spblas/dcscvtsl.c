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


/* Created:  Sat Jul 6 14:30:43 EDT 1996 */

#include "dcscvtsl.h"



void CSC_VecTriangSlvLD_CAB_double( 
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
}

void CSC_VecTriangSlvLU_CAB_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
}

void CSC_VecTriangSlvUD_CAB_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
}

void CSC_VecTriangSlvUU_CAB_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
}




void CSC_VecTriangSlvLD_CaAB_double( 
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CaAB_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CaAB_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CaAB_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CABC_double( 
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CABC_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CABC_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CABC_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CaABC_double( 
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CaABC_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CaABC_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CaABC_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CABbC_double( 
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CABbC_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CABbC_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CABbC_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CaABbC_double( 
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CaABbC_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CaABbC_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CaABbC_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CDAB_double( 
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CDAB_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CDAB_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CDAB_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CaDAB_double( 
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CaDAB_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CaDAB_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CaDAB_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CDABC_double( 
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CDABC_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CDABC_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CDABC_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CaDABC_double( 
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CaDABC_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CaDABC_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CaDABC_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CDABbC_double( 
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CDABbC_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CDABbC_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CDABbC_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CaDABbC_double( 
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CaDABbC_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CaDABbC_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CaDABbC_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CADB_double( 
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
}

void CSC_VecTriangSlvLU_CADB_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
}

void CSC_VecTriangSlvUD_CADB_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
}

void CSC_VecTriangSlvUU_CADB_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
}




void CSC_VecTriangSlvLD_CaADB_double( 
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CaADB_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CaADB_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CaADB_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CADBC_double( 
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CADBC_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CADBC_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CADBC_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CaADBC_double( 
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CaADBC_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CaADBC_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CaADBC_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CADBbC_double( 
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CADBbC_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CADBbC_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CADBbC_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CaADBbC_double( 
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CaADBbC_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CaADBbC_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CaADBbC_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CDADB_double( 
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CDADB_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CDADB_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CDADB_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CaDADB_double( 
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CaDADB_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CaDADB_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CaDADB_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CDADBC_double( 
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CDADBC_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CDADBC_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CDADBC_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CaDADBC_double( 
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CaDADBC_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CaDADBC_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CaDADBC_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ =  (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CDADBbC_double( 
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CDADBbC_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CDADBbC_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CDADBbC_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSC_VecTriangSlvLD_CaDADBbC_double( 
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      jb = pntrb[i];
      je = pntre[i];
      z =  pc[i] / val[jb];
      pc[i] = z;
      for (j=jb+1;j<je;j++) {
        c[indx[j]] -= z*val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvLU_CaDADBbC_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
       jb = pntrb[i];
       je = pntre[i];
       z = pc[i];
       for (j=jb;j<je;j++)
         c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUD_CaDADBbC_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i]-1;
      z = pc[i] /val[je];
      pc[i] = z;
      for (j=jb;j<je;j++) {
          c[indx[j]] -= z * val[j];
      }
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSC_VecTriangSlvUU_CaDADBbC_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++){
      *pwork++ = beta * (*pc);                    
      *pc = dvr[i]*b[i];
      pc++;
    }                                     

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      jb = pntrb[i];
      je = pntre[i];
      z = pc[i];
      pc[i] = z;
      for (j=jb;j<je;j++)
        c[indx[j]] -= z * val[j];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

