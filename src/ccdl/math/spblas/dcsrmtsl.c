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

#include "dcsrmtsl.h"



void CSR_MatTriangSlvLD_CAB_double(
                 const int m, const int n,  
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CAB_double(
                 const int m, const int n,  
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CAB_double(
                 const int m, const int n,  
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CAB_double(
                 const int m, const int n,  
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CaAB_double(
                 const int m, const int n,  
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CaAB_double(
                 const int m, const int n,  
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CaAB_double(
                 const int m, const int n,  
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CaAB_double(
                 const int m, const int n,  
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CABC_double(
                 const int m, const int n,  
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CABC_double(
                 const int m, const int n,  
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CABC_double(
                 const int m, const int n,  
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CABC_double(
                 const int m, const int n,  
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CaABC_double(
                 const int m, const int n,  
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CaABC_double(
                 const int m, const int n,  
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CaABC_double(
                 const int m, const int n,  
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CaABC_double(
                 const int m, const int n,  
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CABbC_double(
                 const int m, const int n,  
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CABbC_double(
                 const int m, const int n,  
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CABbC_double(
                 const int m, const int n,  
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CABbC_double(
                 const int m, const int n,  
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CaABbC_double(
                 const int m, const int n,  
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CaABbC_double(
                 const int m, const int n,  
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CaABbC_double(
                 const int m, const int n,  
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CaABbC_double(
                 const int m, const int n,  
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CDAB_double(
                 const int m, const int n, const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CDAB_double(
                 const int m, const int n, const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CDAB_double(
                 const int m, const int n, const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CDAB_double(
                 const int m, const int n, const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CaDAB_double(
                 const int m, const int n, const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CaDAB_double(
                 const int m, const int n, const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CaDAB_double(
                 const int m, const int n, const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CaDAB_double(
                 const int m, const int n, const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CDABC_double(
                 const int m, const int n, const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CDABC_double(
                 const int m, const int n, const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CDABC_double(
                 const int m, const int n, const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CDABC_double(
                 const int m, const int n, const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CaDABC_double(
                 const int m, const int n, const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CaDABC_double(
                 const int m, const int n, const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CaDABC_double(
                 const int m, const int n, const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CaDABC_double(
                 const int m, const int n, const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CDABbC_double(
                 const int m, const int n, const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CDABbC_double(
                 const int m, const int n, const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CDABbC_double(
                 const int m, const int n, const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CDABbC_double(
                 const int m, const int n, const double *dvl, 
                   const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CaDABbC_double(
                 const int m, const int n, const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CaDABbC_double(
                 const int m, const int n, const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CaDABbC_double(
                 const int m, const int n, const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CaDABbC_double(
                 const int m, const int n, const double *dvl, 
                  const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CADB_double(
                 const int m, const int n,  
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CADB_double(
                 const int m, const int n,  
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CADB_double(
                 const int m, const int n,  
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CADB_double(
                 const int m, const int n,  
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CaADB_double(
                 const int m, const int n,  
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CaADB_double(
                 const int m, const int n,  
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CaADB_double(
                 const int m, const int n,  
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CaADB_double(
                 const int m, const int n,  
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc,  const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CADBC_double(
                 const int m, const int n,  
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CADBC_double(
                 const int m, const int n,  
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CADBC_double(
                 const int m, const int n,  
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CADBC_double(
                 const int m, const int n,  
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CaADBC_double(
                 const int m, const int n,  
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CaADBC_double(
                 const int m, const int n,  
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CaADBC_double(
                 const int m, const int n,  
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CaADBC_double(
                 const int m, const int n,  
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CADBbC_double(
                 const int m, const int n,  
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CADBbC_double(
                 const int m, const int n,  
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CADBbC_double(
                 const int m, const int n,  
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CADBbC_double(
                 const int m, const int n,  
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CaADBbC_double(
                 const int m, const int n,  
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CaADBbC_double(
                 const int m, const int n,  
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CaADBbC_double(
                 const int m, const int n,  
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CaADBbC_double(
                 const int m, const int n,  
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CDADB_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CDADB_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CDADB_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CDADB_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CaDADB_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je, index;
  double *pc=c;
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CaDADB_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CaDADB_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CaDADB_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  *work=0.;
  int i, j, jb, je;
  double *pc=c;
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;


  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CDADBC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CDADBC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CDADBC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CDADBC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CaDADBC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ =  (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CaDADBC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CaDADBC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CaDADBC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb,  
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ =  (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CDADBbC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CDADBbC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CDADBbC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CDADBbC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr,  const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}




void CSR_MatTriangSlvLD_CaDADBbC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je, index;
  double *pc=c;
  double *pwork=work;                             
  int l;                                              
  double z; 
  double valtmp=0.;
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                       
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++); 

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvLU_CaDADBbC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUD_CaDADBbC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

void CSR_MatTriangSlvUU_CaDADBbC_double(
                 const int m, const int n, const double *dvl, 
                 const double *dvr, const double alpha, const double *val, 
                 const int *indx, const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta, 
                 double *c, const int ldc, double *work, const int ind_base)
{
  int i, j, jb, je;
  double *pc=c;
  double *pwork=work;                                
  int l;                                              
  double z; 
  val-=ind_base;
  indx-=ind_base;

  for (l=0;l!=n;l++)                          
    for (i=0;i!=m;i++) *pwork++ = beta * (*pc++);    

  pwork=work;                                        
  c-=ind_base;
  for (l=0;l!=n;l++) {                                
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
    c += ldc; b += ldb;                               
  }                                                   
}

