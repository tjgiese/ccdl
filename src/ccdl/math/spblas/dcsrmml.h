#ifndef dcsrmml_DOUBLE_H
#define dcsrmml_DOUBLE_H
void CSR_MatMult_CAB_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base);

void CSRsymm_MatMult_CAB_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base);

void CSRskew_MatMult_CAB_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base);

void CSR_MatMult_CaAB_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base);

void CSRsymm_MatMult_CaAB_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base);

void CSRskew_MatMult_CaAB_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base);

void CSR_MatMult_CABC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base);

void CSRsymm_MatMult_CABC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base);

void CSRskew_MatMult_CABC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base);

void CSR_MatMult_CaABC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base);

void CSRsymm_MatMult_CaABC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base);

void CSRskew_MatMult_CaABC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, 
                 double *c, const int ldc, 
                 const int ind_base);

void CSR_MatMult_CABbC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base);

void CSRsymm_MatMult_CABbC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base);

void CSRskew_MatMult_CABbC_double(
                 const int m, const int n, const int k, 
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base);

void CSR_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base);

void CSRsymm_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base);

void CSRskew_MatMult_CaABbC_double(
                 const int m, const int n, const int k, const double alpha,
                 const double *val, const int *indx, 
                 const int *pntrb, const int *pntre, 
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, 
                 const int ind_base);
#endif
