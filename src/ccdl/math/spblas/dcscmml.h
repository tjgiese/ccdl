#ifndef dcscmml_DOUBLE_H
#define dcscmml_DOUBLE_H
void CSC_MatMult_CAB_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base);

void CSCsymm_MatMult_CAB_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base);

void CSCskew_MatMult_CAB_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base);

void CSC_MatMult_CaAB_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base);

void CSCsymm_MatMult_CaAB_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base);

void CSCskew_MatMult_CaAB_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base);

void CSC_MatMult_CABC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base);

void CSCsymm_MatMult_CABC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base);

void CSCskew_MatMult_CABC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base);

void CSC_MatMult_CaABC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base);

void CSCsymm_MatMult_CaABC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base);

void CSCskew_MatMult_CaABC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc,  
                 const int ind_base);

void CSC_MatMult_CABbC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base);

void CSCsymm_MatMult_CABbC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base);

void CSCskew_MatMult_CABbC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base);

void CSC_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base);

void CSCsymm_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base);

void CSCskew_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc,  
                 const int ind_base);
#endif
