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


/* Created:  Sat Jul 6 14:28:33 EDT 1996 */

#include "dcsrvtsl.h"



void CSR_VecTriangSlvLD_CAB_double(
                 const int m,   
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,   const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (b[i] - z) / valtmp;
    }
}

void CSR_VecTriangSlvLU_CAB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
}

void CSR_VecTriangSlvUD_CAB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (b[i] - z) / val[jb];
    }
}

void CSR_VecTriangSlvUU_CAB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
}




void CSR_VecTriangSlvLD_CaAB_double(
                 const int m,   
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,   const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CaAB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CaAB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CaAB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CABC_double(
                 const int m,   
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CABC_double(
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

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CABC_double(
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

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CABC_double(
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

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CaABC_double(
                 const int m,   
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CaABC_double(
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

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CaABC_double(
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

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CaABC_double(
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

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CABbC_double(
                 const int m,   
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CABbC_double(
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

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CABbC_double(
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

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CABbC_double(
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

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CaABbC_double(
                 const int m,   
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CaABbC_double(
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

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CaABbC_double(
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

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CaABbC_double(
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

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CDAB_double(
                 const int m,  const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CDAB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CDAB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CDAB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CaDAB_double(
                 const int m,  const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CaDAB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CaDAB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CaDAB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CDABC_double(
                 const int m,  const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CDABC_double(
                 const int m,  const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CDABC_double(
                 const int m,  const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CDABC_double(
                 const int m,  const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CaDABC_double(
                 const int m,  const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CaDABC_double(
                 const int m,  const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CaDABC_double(
                 const int m,  const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CaDABC_double(
                 const int m,  const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CDABbC_double(
                 const int m,  const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CDABbC_double(
                 const int m,  const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CDABbC_double(
                 const int m,  const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CDABbC_double(
                 const int m,  const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CaDABbC_double(
                 const int m,  const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CaDABbC_double(
                 const int m,  const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CaDABbC_double(
                 const int m,  const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CaDABbC_double(
                 const int m,  const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CADB_double(
                 const int m,   
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,   const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (dvr[i]*b[i] - z) / valtmp;
    }
}

void CSR_VecTriangSlvLU_CADB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
}

void CSR_VecTriangSlvUD_CADB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (dvr[i]*b[i] - z) / val[jb];
    }
}

void CSR_VecTriangSlvUU_CADB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
}




void CSR_VecTriangSlvLD_CaADB_double(
                 const int m,   
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,   const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (dvr[i]*b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CaADB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CaADB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (dvr[i]*b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CaADB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CADBC_double(
                 const int m,   
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (dvr[i]*b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CADBC_double(
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

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CADBC_double(
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

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (dvr[i]*b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CADBC_double(
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

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CaADBC_double(
                 const int m,   
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (dvr[i]*b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CaADBC_double(
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

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CaADBC_double(
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

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (dvr[i]*b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CaADBC_double(
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

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CADBbC_double(
                 const int m,   
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (dvr[i]*b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CADBbC_double(
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

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CADBbC_double(
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

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (dvr[i]*b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CADBbC_double(
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

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CaADBbC_double(
                 const int m,   
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (dvr[i]*b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CaADBbC_double(
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

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CaADBbC_double(
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

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (dvr[i]*b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CaADBbC_double(
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

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha ;            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CDADB_double(
                 const int m,  const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (dvr[i]*b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CDADB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CDADB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (dvr[i]*b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CDADB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CaDADB_double(
                 const int m,  const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (dvr[i]*b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CaDADB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CaDADB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (dvr[i]*b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CaDADB_double(
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


  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CDADBC_double(
                 const int m,  const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (dvr[i]*b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CDADBC_double(
                 const int m,  const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CDADBC_double(
                 const int m,  const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (dvr[i]*b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CDADBC_double(
                 const int m,  const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CaDADBC_double(
                 const int m,  const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (dvr[i]*b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CaDADBC_double(
                 const int m,  const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CaDADBC_double(
                 const int m,  const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (dvr[i]*b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CaDADBC_double(
                 const int m,  const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,   
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CDADBbC_double(
                 const int m,  const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (dvr[i]*b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CDADBbC_double(
                 const int m,  const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CDADBbC_double(
                 const int m,  const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (dvr[i]*b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CDADBbC_double(
                 const int m,  const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *=  dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}




void CSR_VecTriangSlvLD_CaDADBbC_double(
                 const int m,  const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        index = indx[j];
        if ( index == i+ind_base ) {
          valtmp = val[j];
        } else {
          z += c[index] * val[j];
        }
      }
      pc[i] = (dvr[i]*b[i] - z) / valtmp;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvLU_CaDADBbC_double(
                 const int m,  const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=0;i!=m;i++) {
      z = 0;
      jb = pntrb[i];
      je = pntre[i];
      for (j=jb;j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUD_CaDADBbC_double(
                 const int m,  const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1; i--) {
      z = 0;
      jb = pntrb[i];
      je =  pntre[i];
      for (j=jb+1; j<je; j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = (dvr[i]*b[i] - z) / val[jb];
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

void CSR_VecTriangSlvUU_CaDADBbC_double(
                 const int m,  const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b,  const double beta, 
                 double *c,  double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  double z; 
  val-=ind_base;
  indx-=ind_base;

    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
    pc=c+ind_base;
    for (i=m-1;i!=-1;i--) {
      z = 0;
      jb =  pntrb[i];
      je =  pntre[i];
      for (j=jb; j<je;j++) {
        z += c[indx[j]] * val[j];
      }
      pc[i] = dvr[i]*b[i] - z;
    }
    for (i=0;i!=m;i++) {                  
        *pc *= alpha * dvl[i];            
        *pc += *pwork;                    
        pwork++;                          
        pc++;                             
    }                                     
}

