#ifndef dcscmtsl_DOUBLE_H
#define dcscmtsl_DOUBLE_H
void CSC_MatTriangSlvLD_CAB_double( 
                 const int m, const int n, 
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvLU_CAB_double(
                 const int m, const int n, 
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvUD_CAB_double(
                 const int m, const int n, 
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvUU_CAB_double(
                 const int m, const int n, 
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvLD_CaAB_double( 
                 const int m, const int n, 
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvLU_CaAB_double(
                 const int m, const int n, 
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvUD_CaAB_double(
                 const int m, const int n, 
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvUU_CaAB_double(
                 const int m, const int n, 
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvLD_CABC_double( 
                 const int m, const int n, 
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CABC_double(
                 const int m, const int n, 
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CABC_double(
                 const int m, const int n, 
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CABC_double(
                 const int m, const int n, 
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CaABC_double( 
                 const int m, const int n, 
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CaABC_double(
                 const int m, const int n, 
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CaABC_double(
                 const int m, const int n, 
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CaABC_double(
                 const int m, const int n, 
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CABbC_double( 
                 const int m, const int n, 
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CABbC_double(
                 const int m, const int n, 
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CABbC_double(
                 const int m, const int n, 
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CABbC_double(
                 const int m, const int n, 
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CaABbC_double( 
                 const int m, const int n, 
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CaABbC_double(
                 const int m, const int n, 
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CaABbC_double(
                 const int m, const int n, 
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CaABbC_double(
                 const int m, const int n, 
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CDAB_double( 
                 const int m, const int n, const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CDAB_double(
                 const int m, const int n, const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CDAB_double(
                 const int m, const int n, const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CDAB_double(
                 const int m, const int n, const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CaDAB_double( 
                 const int m, const int n, const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CaDAB_double(
                 const int m, const int n, const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CaDAB_double(
                 const int m, const int n, const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CaDAB_double(
                 const int m, const int n, const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CDABC_double( 
                 const int m, const int n, const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CDABC_double(
                 const int m, const int n, const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CDABC_double(
                 const int m, const int n, const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CDABC_double(
                 const int m, const int n, const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CaDABC_double( 
                 const int m, const int n, const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CaDABC_double(
                 const int m, const int n, const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CaDABC_double(
                 const int m, const int n, const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CaDABC_double(
                 const int m, const int n, const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CDABbC_double( 
                 const int m, const int n, const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CDABbC_double(
                 const int m, const int n, const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CDABbC_double(
                 const int m, const int n, const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CDABbC_double(
                 const int m, const int n, const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CaDABbC_double( 
                 const int m, const int n, const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CaDABbC_double(
                 const int m, const int n, const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CaDABbC_double(
                 const int m, const int n, const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CaDABbC_double(
                 const int m, const int n, const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CADB_double( 
                 const int m, const int n, 
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvLU_CADB_double(
                 const int m, const int n, 
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvUD_CADB_double(
                 const int m, const int n, 
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvUU_CADB_double(
                 const int m, const int n, 
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvLD_CaADB_double( 
                 const int m, const int n, 
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvLU_CaADB_double(
                 const int m, const int n, 
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvUD_CaADB_double(
                 const int m, const int n, 
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvUU_CaADB_double(
                 const int m, const int n, 
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc,  const int ind_base);

void CSC_MatTriangSlvLD_CADBC_double( 
                 const int m, const int n, 
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CADBC_double(
                 const int m, const int n, 
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CADBC_double(
                 const int m, const int n, 
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CADBC_double(
                 const int m, const int n, 
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CaADBC_double( 
                 const int m, const int n, 
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CaADBC_double(
                 const int m, const int n, 
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CaADBC_double(
                 const int m, const int n, 
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CaADBC_double(
                 const int m, const int n, 
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CADBbC_double( 
                 const int m, const int n, 
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CADBbC_double(
                 const int m, const int n, 
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CADBbC_double(
                 const int m, const int n, 
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CADBbC_double(
                 const int m, const int n, 
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CaADBbC_double( 
                 const int m, const int n, 
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CaADBbC_double(
                 const int m, const int n, 
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CaADBbC_double(
                 const int m, const int n, 
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CaADBbC_double(
                 const int m, const int n, 
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CDADB_double( 
                 const int m, const int n, const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CDADB_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CDADB_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CDADB_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CaDADB_double( 
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CaDADB_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CaDADB_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CaDADB_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CDADBC_double( 
                 const int m, const int n, const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CDADBC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CDADBC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CDADBC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CaDADBC_double( 
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CaDADBC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CaDADBC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CaDADBC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, 
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CDADBbC_double( 
                 const int m, const int n, const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CDADBbC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CDADBbC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CDADBbC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLD_CaDADBbC_double( 
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvLU_CaDADBbC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUD_CaDADBbC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);

void CSC_MatTriangSlvUU_CaDADBbC_double(
                 const int m, const int n, const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b, const int ldb, const double beta,
                 double *c, const int ldc, double *work, const int ind_base);
#endif
