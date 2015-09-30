#ifndef _ccdl_bspline_hpp_
#define _ccdl_bspline_hpp_


namespace ccdl
{

  void bspline_one_pass( double * c, double const w, int const n );

  void bspline_diff( int const order, double const * array, double * diff );
  
  void bspline_eval( double const w, int const order, double * array );
  
  void bspline_eval( double const w, int const order, double * array, double * darray );
  
  void bspline_eval( double const w, int const order, int const nder, double * array );
  


  void bspline_periodic( double x, double const lenx, int const nx, int const order, double * w );
  void bspline_periodic_deriv( double x, double const lenx, int const nx, int const order, double * w, double * dw );

  void bspline_periodic( double x, double const lenx, int const nx, int const order, double * w, int * gidx );
  
  void bspline_periodic( double x, double const lenx, int const nx, int const order, double * w, int * gidx, double * delta );
  
  
  void bspline_periodic_deriv( double x, double const lenx, int const nx, int const order, double * w, double * dw, int * gidx );
  
  void bspline_periodic_deriv( double x, double const lenx, int const nx, int const order, double * w, double * dw, int * gidx, double * delta );
  


  void bspline_aperiodic( double x, double const lenx, int const nx, int const order, double * w, int * gidx );
  
  void bspline_aperiodic( double x, double const lenx, int const nx, int const order, double * w, int * gidx, double * delta );
  
  void bspline_aperiodic_deriv( double x, double const lenx, int const nx, int const order, double * w, double * dw, int * gidx );
  
  void bspline_aperiodic_deriv( double x, double const lenx, int const nx, int const order, double * w, double * dw, int * gidx, double * delta );




  void bspline_periodic_nderiv( double x, double const lenx, int const nx, int const order, int const nder, double * dwdx, int * gidx, double * delta );

  void bspline_aperiodic_nderiv( double x, double const lenx, int const nx, int const order, int const nder, double * dwdx, int * gidx, double * delta );

  void bspline_dct( int const N, int const order, double const * in, double * out );

  void bspline_hdct( int const N, int const order, double const * in, double * out );

  void bspline_fourier_full( int const N, int const order, double * fc );

  void bspline_fourier_half( int const N, int const order, double * fc );

  void bspline_renormalize( double const lx, int const nx, int const order, double * data );
  void bspline_renormalize( double const lx, double const ly, int const nx, int const ny, int const order, double * data );

}



#endif

