#ifndef dcscvml_DOUBLE_H
#define dcscvml_DOUBLE_H
void CSC_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base);

void CSCsymm_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base);

void CSCskew_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base);

void CSC_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base);

void CSCsymm_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base);

void CSCskew_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base);

void CSC_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base);

void CSCsymm_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base);

void CSCskew_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base);

void CSC_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base);

void CSCsymm_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base);

void CSCskew_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,   
                 const int ind_base);

void CSC_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,   
                 const int ind_base);

void CSCsymm_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,   
                 const int ind_base);

void CSCskew_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,   
                 const int ind_base);

void CSC_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,   
                 const int ind_base);

void CSCsymm_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,   
                 const int ind_base);

void CSCskew_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,   
                 const int ind_base);
#endif
