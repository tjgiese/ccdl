#ifndef _ccdl_linalg_h_
#define _ccdl_linalg_h_


namespace ccdl
{

  double v_dot_v( int n, double const * A, double const * B );

  void ge_dot_v( int nfA, int nsA, double const * A, double const * x,double * Ax,
		 double alpha = 1., double beta  = 0. );

  void gt_dot_v( int nfA, int nsA, double const * A, double const * x, double * Ax,
		 double alpha = 1., double beta  = 0. );

  void sy_dot_v( int nfA, int nsA, double const * A, double const * x, double * Ax,
		 double alpha = 1., double beta  = 0. );

  void ge_dot_ge( int nfA, int nsA, double const * A, int nfB, int nsB, double const * B, double * AB,
		  double alpha = 1., double beta  = 0. );

  void gt_dot_ge( int nfA, int nsA, double const * A, int nfB, int nsB, double const * B, double * AB,
		  double alpha = 1., double beta  = 0. );

  void ge_dot_gt( int nfA, int nsA, double const * A, int nfB, int nsB, double const * B, double * AB,
		  double alpha = 1., double beta  = 0. );

  void di_dot_gt( int nfA, int nsA, double const * A, int nfB, int nsB, double const * B,double * AB,
		  double alpha = 1.,double beta  = 0. );

  void ge_dot_di( int nfA, int nsA, double const * A, int nfB, int nsB, double const * B,double * AB,
		  double alpha = 1.,double beta  = 0. );

  void di_dot_ge( int nfA, int nsA, double const * A, int nfB, int nsB, double const * B,double * AB,
		  double alpha = 1.,double beta  = 0. );


  void sy_dot_ge( int nfA, int nsA, double const * A, int nfB, int nsB, double const * B, double * AB,
		  double alpha = 1., double beta  = 0. );

  void ge_dot_sy( int nfA, int nsA, double const * A, int nfB, int nsB, double const * B, double * AB,
		  double alpha = 1., double beta  = 0. );


  int query_sdd( int const M, int const N );

  int sdd_decomp( int const nf, int const ns, double const * A,
		  double * U, double * w, double * VT );

  int sdd_power( double const power, 
		  int const nf, int const ns, double * A, 
		  double TOL = 1.e-20 );

  int sdd_sym_power( double const power, 
		      int const nf, int const ns, double * A, 
		      double TOL = 1.e-20 );

  int sym_inverse( int const nf, int const ns, double * A );

}


#endif

