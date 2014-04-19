#ifndef _ccdl_bspline_hpp_
#define _ccdl_bspline_hpp_


namespace ccdl
{

  void one_pass_bspline( double * c, double const w, int const n );

  void diff_bspline( int const order, double const * array, double * diff );
  
  void fill_bspline_0( double const w, int const order, double * array );
  
  void fill_bspline_1( double const w, int const order, double * array, double * darray );
  
  void fill_bspline_n( double const w, int const order, int const nder, double * array );
  


  void periodicBspline( double x, double const lenx, int const nx, int const order, double * w );
  void periodicBsplineDeriv( double x, double const lenx, int const nx, int const order, double * w, double * dw );

  void periodicBspline( double x, double const lenx, int const nx, int const order, double * w, int * gidx );
  
  void periodicBspline( double x, double const lenx, int const nx, int const order, double * w, int * gidx, double * delta );
  
  
  void periodicBsplineDeriv( double x, double const lenx, int const nx, int const order, double * w, double * dw, int * gidx );
  
  void periodicBsplineDeriv( double x, double const lenx, int const nx, int const order, double * w, double * dw, int * gidx, double * delta );
  


  void aperiodicBspline( double x, double const lenx, int const nx, int const order, double * w, int * gidx );
  
  void aperiodicBspline( double x, double const lenx, int const nx, int const order, double * w, int * gidx, double * delta );
  
  void aperiodicBsplineDeriv( double x, double const lenx, int const nx, int const order, double * w, double * dw, int * gidx );
  
  void aperiodicBsplineDeriv( double x, double const lenx, int const nx, int const order, double * w, double * dw, int * gidx, double * delta );




  void periodicBsplineNthDeriv( double x, double const lenx, int const nx, int const order, int const nder, double * dwdx, int * gidx, double * delta );

  void aperiodicBsplineNthDeriv( double x, double const lenx, int const nx, int const order, int const nder, double * dwdx, int * gidx, double * delta );

}



#endif

