#ifndef dcscvtsl_DOUBLE_H
#define dcscvtsl_DOUBLE_H
void CSC_VecTriangSlvLD_CAB_double( 
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvLU_CAB_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvUD_CAB_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvUU_CAB_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvLD_CaAB_double( 
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvLU_CaAB_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvUD_CaAB_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvUU_CaAB_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvLD_CABC_double( 
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CABC_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CABC_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CABC_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CaABC_double( 
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CaABC_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CaABC_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CaABC_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CABbC_double( 
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CABbC_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CABbC_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CABbC_double(
                 const int m,  
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CaABbC_double( 
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CaABbC_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CaABbC_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CaABbC_double(
                 const int m,  
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CDAB_double( 
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CDAB_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CDAB_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CDAB_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CaDAB_double( 
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CaDAB_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CaDAB_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CaDAB_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CDABC_double( 
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CDABC_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CDABC_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CDABC_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CaDABC_double( 
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CaDABC_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CaDABC_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CaDABC_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CDABbC_double( 
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CDABbC_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CDABbC_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CDABbC_double(
                 const int m,  const double *dvl,
                   const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CaDABbC_double( 
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CaDABbC_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CaDABbC_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CaDABbC_double(
                 const int m,  const double *dvl,
                  const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CADB_double( 
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvLU_CADB_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvUD_CADB_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvUU_CADB_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvLD_CaADB_double( 
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvLU_CaADB_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvUD_CaADB_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvUU_CaADB_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,   const int ind_base);

void CSC_VecTriangSlvLD_CADBC_double( 
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CADBC_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CADBC_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CADBC_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CaADBC_double( 
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CaADBC_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CaADBC_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CaADBC_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CADBbC_double( 
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CADBbC_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CADBbC_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CADBbC_double(
                 const int m,  
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CaADBbC_double( 
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CaADBbC_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CaADBbC_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CaADBbC_double(
                 const int m,  
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CDADB_double( 
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CDADB_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CDADB_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CDADB_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CaDADB_double( 
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CaDADB_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CaDADB_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CaDADB_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CDADBC_double( 
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CDADBC_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CDADBC_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CDADBC_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CaDADBC_double( 
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CaDADBC_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CaDADBC_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CaDADBC_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CDADBbC_double( 
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CDADBbC_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CDADBbC_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CDADBbC_double(
                 const int m,  const double *dvl,
                 const double *dvr,  const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLD_CaDADBbC_double( 
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvLU_CaDADBbC_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUD_CaDADBbC_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);

void CSC_VecTriangSlvUU_CaDADBbC_double(
                 const int m,  const double *dvl,
                 const double *dvr, const double alpha, const double *val,
                 const int *indx, const int *pntrb, const int *pntre,
                 const double *b,  const double beta,
                 double *c,  double *work, const int ind_base);
#endif
