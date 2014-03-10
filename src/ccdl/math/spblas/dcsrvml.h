#ifndef dcsrvml_DOUBLE_H
#define dcsrvml_DOUBLE_H
void CSR_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);

void CSRsymm_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);

void CSRskew_VecMult_CAB_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);

void CSR_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);

void CSRsymm_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);

void CSRskew_VecMult_CaAB_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);

void CSR_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);

void CSRsymm_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);

void CSRskew_VecMult_CABC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);

void CSR_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);

void CSRsymm_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);

void CSRskew_VecMult_CaABC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  
                 double *c,  
                 const int ind_base);

void CSR_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base);

void CSRsymm_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base);

void CSRskew_VecMult_CABbC_double(
                 const int m,  const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base);

void CSR_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base);

void CSRsymm_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base);

void CSRskew_VecMult_CaABbC_double(
                 const int m,  const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b,  const double beta,
                 double *c,  
                 const int ind_base);
#endif
