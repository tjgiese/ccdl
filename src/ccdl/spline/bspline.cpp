#include "bspline.hpp"

#include <iostream>
#include <cmath>
#include <cassert>

#include "../constants.hpp"

#ifndef alloca
#define alloca __builtin_alloca
#endif


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


void ccdl::bspline_one_pass( double * c, double const w, int const n )
{
  int nm1 = n-1;
  double const div = 1. / nm1;
  c[nm1] = div * w * c[nm1 - 1];
  for ( int j=1; j<nm1; ++j )
    c[nm1 - j] = div * ((w + j) * c[nm1 - j - 1] + (n - j - w) * c[nm1 - j]);
  c[0] = div * (1 - w) * c[0];
}
void ccdl::bspline_diff( int const order, double const * array, double * diff )
{
  assert( order > 1 );
  int const nm1 = order-1;
  diff[0] = -array[0];
  for ( int j=1; j<nm1; ++j )
    diff[j] = array[j-1] - array[j];
  diff[nm1] = array[nm1-1];
}

void ccdl::bspline_eval( double const w, int const order, double * array )
{
  array[0] = 1. - w;
  array[1] = w;
  if (order > 2)
    {
      // One pass to order 3:
      array[2] = 0.5 * w * array[1];
      array[1] = 0.5 * ((w + 1.) * array[0] + (2. - w) * array[1]);
      array[0] = 0.5 * (1. - w) * array[0];
      if (order > 3)
        {
          // One pass to order 4:         
          double const div = 1./3.;
          array[3] = div * w * array[2];
          array[2] = div * ((w + 1.) * array[1] + (3. - w) * array[2]);
          array[1] = div * ((w + 2.) * array[0] + (2. - w) * array[1]);
          array[0] = div * (1. - w) * array[0];
          // and the rest
          for ( int k = 5; k < order+1; ++k )
            ccdl::bspline_one_pass(array, w, k);
        };
    };
}

void ccdl::bspline_eval( double const w, int const order, double * array, double * darray )
{
  assert( order > 2 );
  double const div = 1./3.;

  array[0] = 1. - w;
  array[1] = w;

  if (order == 4)
    {
      // One pass to order 3:
      array[2] = 0.5 * w * array[1];
      array[1] = 0.5 * ((w + 1.) * array[0] + (2. - w) * array[1]);
      array[0] = 0.5 * (1. - w) * array[0];
      
      darray[0] = -array[0];
      darray[1] = array[0]-array[1];
      darray[2] = array[1]-array[2];
      darray[3] = array[2];

      // One pass to order 4:     
      array[3] = div * w * array[2];
      array[2] = div * ((w + 1.) * array[1] + (3. - w) * array[2]);
      array[1] = div * ((w + 2.) * array[0] + (2. - w) * array[1]);
      array[0] = div * (1. - w) * array[0];
      
    }
  else if ( order > 4 )
    {
      array[2] = 0.5 * w * array[1];
      array[1] = 0.5 * ((w + 1.) * array[0] + (2. - w) * array[1]);
      array[0] = 0.5 * (1. - w) * array[0];

      array[3] = div * w * array[2];
      array[2] = div * ((w + 1.) * array[1] + (3. - w) * array[2]);
      array[1] = div * ((w + 2.) * array[0] + (2. - w) * array[1]);
      array[0] = div * (1. - w) * array[0];

      // and the rest
      for ( int k = 5; k < order; ++k ) // don't do k==order
        ccdl::bspline_one_pass(array, w, k);

      ccdl::bspline_diff(order,array,darray);

      // One more recursion: // do the k==order
      ccdl::bspline_one_pass(array, w, order);

    }
  else // order == 3
    {
      darray[0] = -array[0];
      darray[1] = array[0]-array[1];
      darray[2] = array[1];

      // One pass to order 3:
      array[2] = 0.5 * w * array[1];
      array[1] = 0.5 * ((w + 1.) * array[0] + (2. - w) * array[1]);
      array[0] = 0.5 * (1. - w) * array[0];
    };
}

void ccdl::bspline_eval( double const w, int const order, int const nder, double * array )
{
  assert( order > 2 );
  double ** darray = (double**)alloca ( sizeof(double**) * (nder+1) );
  //std::vector<double *> darray(nder+1);
  for ( int ider=0; ider<nder+1; ++ider )
    darray[ider] = array+ider*order;
  for ( int i=0; i<(nder+1)*order; ++i )
    array[i] = 0.;
  array[0] = 1.;
  for ( int d=0, nd=std::min(nder,order-1); d<nd; ++d )
    ccdl::bspline_diff(order,darray[d],darray[d+1]);
  array[0] = 1. - w;
  array[1] = w;
  for ( int o=3; o<=order; ++o )
    {
      for ( int d=0, nd=std::min(order-(o-1),nder); d<nd; ++d )
        bspline_diff(order,darray[d],darray[d+1]);
      ccdl::bspline_one_pass(array, w, o);
    };
}



///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


void ccdl::bspline_periodic
( double x, 
  double const lenx, 
  int const nx,
  int const order,
  double * w )
{
  x = x*nx/lenx - (order%2)*0.5;
  ccdl::bspline_eval(x-std::floor(x),order,w);
}


void ccdl::bspline_periodic
( double x, 
  double const lenx, 
  int const nx,
  int const order,
  double * w,
  int * gidx )
{
  x = x*nx/lenx - (order%2)*0.5;
  int ilo = std::floor(x);
  x -= ilo;
  ilo = ((ilo+1-order/2)%nx+nx)%nx;
  for ( int b=0; b<order; ++b, ++ilo )
    gidx[b] = ilo%nx;
  ///////////////////////////
  ccdl::bspline_eval(x,order,w);
}


void ccdl::bspline_periodic
( double x, 
  double const lenx, 
  int const nx,
  int const order,
  double * w,
  int * gidx, 
  double * delta )
{
  double del = lenx/nx;
  double q = x/del - (order%2)*0.5;
  int ilo = std::floor(q);
  q -= ilo;
  ilo += 1 - order/2;
  delta[0] = del*ilo-x;
  gidx[0]  = ((ilo++)%nx+nx)%nx;
  for ( int b=1; b<order; ++b, ++ilo )
    {
      delta[b] = delta[b-1] + del;
      gidx[b] = ilo%nx;
    };
  ///////////////////////////
  ccdl::bspline_eval(q,order,w);
}





void ccdl::bspline_periodic_deriv
( double x, 
  double const lenx, 
  int const nx,
  int const order,
  double * w,
  double * dw )
{
  x = x*nx/lenx - (order%2)*0.5;
  int ilo = std::floor(x);
  x -= ilo;
  ///////////////////////////
  ccdl::bspline_eval(x,order,w,dw);
  x = nx/lenx;
  for ( int b=0; b<order; ++b )
    dw[b] *= x;
}


void ccdl::bspline_periodic_deriv
( double x, 
  double const lenx, 
  int const nx,
  int const order,
  double * w,
  double * dw,
  int * gidx )
{
  x = x*nx/lenx - (order%2)*0.5;
  int ilo = std::floor(x);
  x -= ilo;
  ilo = ((ilo+1-order/2)%nx+nx)%nx;
  ccdl::bspline_eval(x,order,w,dw);
  x = nx/lenx;
  for ( int b=0; b<order; ++b, ++ilo )
    {
      gidx[b] = ilo%nx;
      dw[b] *= x;
    };
}


void ccdl::bspline_periodic_deriv
( double x, 
  double const lenx, 
  int const nx,
  int const order,
  double * w,
  double * dw,
  int * gidx, 
  double * delta )
{
  double del = lenx/nx;
  double q = x/del - (order%2)*0.5;
  int ilo = std::floor(q);
  ccdl::bspline_eval(q-ilo,order,w,dw);
  ilo += 1 - order/2;
  dw[0] /= del;
  delta[0] = del*ilo-x;
  gidx[0]  = ((ilo++)%nx+nx)%nx;
  //std::cout << order << " : " << gidx[0] << " ";
  for ( int b=1; b<order; ++b, ++ilo )
    {
      dw[b] /= del;
      delta[b] = delta[b-1] + del;
      gidx[b] = ilo%nx;
      //std::cout << gidx[b] << " ";
    };
  //std::cout << "\n";
}


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


void ccdl::bspline_aperiodic
( double x, 
  double const lenx, 
  int const nx,
  int const order, 
  double * w,
  int * gidx )
{
  x = x*nx/lenx - (order%2)*0.5;
  int ilo = std::floor(x);
  ccdl::bspline_eval(x-ilo,order,w);
  ilo += 1 - order/2;
  for ( int b=0; b<order; ++b, ++ilo )
    {
      if ( ilo > -1 and ilo < nx )
	{
	  gidx[b] = ilo;
	}
      else
	{
	  gidx[b] = 0;
	  w[b] = 0.;
	};
    };
}

void ccdl::bspline_aperiodic
( double x, 
  double const lenx, 
  int const nx,
  int const order,
  double * w,
  int * gidx, 
  double * delta )
{
  double del = lenx/nx;
  double q = x/del - (order%2)*0.5;
  int ilo = std::floor(q);
  ccdl::bspline_eval(q-ilo,order,w);
  ilo += 1 - order/2;
  for ( int b=0; b<order; ++b, ++ilo )
    {
      delta[b] = delta[b-1] + del;
      if ( ilo > -1 and ilo < nx )
	{
	  delta[b] = del*ilo-x;
	  gidx[b] = ilo;
	}
      else
	{
	  delta[b] = 0.;
	  gidx[b] = 0;
	  w[b] = 0.;
	};
    };
}

void ccdl::bspline_aperiodic_deriv
( double x, 
  double const lenx, 
  int const nx,
  int const order,
  double * w,
  double * dw,
  int * gidx )
{
  x = x*nx/lenx - (order%2)*0.5;
  int ilo = std::floor(x);
  ccdl::bspline_eval(x-ilo,order,w,dw);
  ilo += 1-order/2;
  x = nx/lenx;
  for ( int b=0; b<order; ++b, ++ilo )
    {
      if ( ilo > -1 and ilo < nx )
	{
	  gidx[b] = ilo;
	  dw[b] *= x;
	}
      else
	{
	  gidx[b] = 0;
	  w[b] = 0.;
	  dw[b] = 0.;
	};
    };
}



void ccdl::bspline_aperiodic_deriv
( double x, 
  double const lenx, 
  int const nx,
  int const order,
  double * w,
  double * dw,
  int * gidx, 
  double * delta )
{
  double del = lenx/nx;
  double q = x/del - (order%2)*0.5;
  int ilo = std::floor(q);
  ccdl::bspline_eval(q-ilo,order,w,dw);
  ilo += 1 - order/2;
  for ( int b=0; b<order; ++b, ++ilo )
    {
      if ( ilo > -1 and ilo < nx )
	{
	  dw[b] /= del;
	  gidx[b] = ilo;
	  delta[b] = del*ilo-x;
	}
      else
	{
	  w[b] = 0.;
	  dw[b] = 0.;
	  gidx[b] = 0;
	  delta[b] = 0.;
	};
    };
}








void ccdl::bspline_periodic_nderiv
( double x, 
  double const lenx, 
  int const nx,
  int const order,
  int const nder,
  double * dwdx,
  int * gidx, 
  double * delta )
{
  double del = lenx/nx;
  double q = x/del - (order%2)*0.5;
  int ilo = std::floor(q);
  ccdl::bspline_eval(q-ilo,order,nder,dwdx);
  ilo += 1 - order/2;
  delta[0] = del*ilo-x;
  gidx[0]  = ((ilo++)%nx+nx)%nx;
  //std::cout << order << " : " << gidx[0] << " ";
  for ( int b=1; b<order; ++b, ++ilo )
    {
      delta[b] = delta[b-1] + del;
      gidx[b] = ilo%nx;
      //std::cout << gidx[b] << " ";
    };
  //std::cout << "\n";
  q = 1.;
  for ( int k=1; k<nder+1; ++k )
    {
      q /= del;
      for ( int i=0; i<order; ++i )
	dwdx[i+k*order] *= q;
    };
}


void ccdl::bspline_aperiodic_nderiv
( double x, 
  double const lenx, 
  int const nx,
  int const order,
  int const nder,
  double * dwdx,
  int * gidx, 
  double * delta )
{
  double del = lenx/nx;
  double q = x/del - (order%2)*0.5;
  int ilo = std::floor(q);
  ccdl::bspline_eval(q-ilo,order,nder,dwdx);
  ilo += 1 - order/2;
  for ( int b=0; b<order; ++b, ++ilo )
    {
      if ( ilo > -1 and ilo < nx )
	{
	  for ( int k=0; k<nder; ++k )
	    dwdx[b+k*order] /= del;
	  gidx[b] = ilo;
	  delta[b] = del*ilo-x;
	}
      else
	{
	  for ( int k=0; k<nder; ++k )
	    dwdx[b+k*order] = 0.;
	  gidx[b] = 0;
	  delta[b] = 0.;
	};
    };
}






////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////




#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>

////////////////////////////////////////////////////////////////////
// discrete cosine expansion; computes the fourier coefficients
// of an even function positioned according to the fftw3 conventions
////////////////////////////////////////////////////////////////////
void ccdl::bspline_dct
( int const N, int const order, double const * in, double * out )
{
  for ( int i=0; i<N; ++i ) out[i] = 0.;
  
  int const o = (order-1)/2;
  int const Np = (order)/2+1;
  int const Nh = N-o;
  
  double const tpion = TWO_PI/N;
  for ( int m=0; m<N; ++m )
    {
      double const k = m * tpion;
      for ( int n=0; n<Np; ++n )
        out[m] += in[n+o] * std::cos(k*n);
      for ( int i=0; i<o; ++i )
        out[m] += in[i] * std::cos(k*(Nh+i));
    };
}


////////////////////////////////////////////////////////////////////
// discrete cosine expansion; computes the fourier coefficients
// of an even function positioned according to the fftw3 conventions
// This version computes only half of the spectrum on account
// of the symmetry of the complex numbers in the fast-loop
////////////////////////////////////////////////////////////////////
void ccdl::bspline_hdct
( int const N, int const order, double const * in, double * out )
{
  for ( int i=0; i<N/2+1; ++i ) out[i] = 0.;
  
  int const Nmax = N/2+1;
  int const o = (order-1)/2;
  int const Np = (order)/2+1;
  int const Nh = N-o;
  
  double const tpion = TWO_PI/N;
  for ( int m=0; m<Nmax; ++m )
    {
      double const k = m * tpion;
      for ( int n=0; n<Np; ++n )
        out[m] += in[n+o] * std::cos(k*n);
      for ( int i=0; i<o; ++i )
        out[m] += in[i] * std::cos(k*(Nh+i));
    };
}


////////////////////////////////////////////////////////////////////
// Computes the fourier transform of a b-spline
////////////////////////////////////////////////////////////////////
void ccdl::bspline_fourier_full
( int const N, int const order, double * fc )
{
  double * wts = (double*)alloca (sizeof(double*)*order);
    //std::vector<double> wts( order, 0. );
  ccdl::bspline_eval( ( order%2 == 0 ? 0.0 : 0.5 ), order, wts );
  ccdl::bspline_dct( N, order, wts, fc );
}

////////////////////////////////////////////////////////////////////
// Computes the fourier transform of a b-spline
// This version only computes half of the fourier spectrum on
// account of the complex number symmetry within the fast loop
////////////////////////////////////////////////////////////////////
void ccdl::bspline_fourier_half
( int const N, int const order, double * fc )
{
  double * wts = (double*)alloca (sizeof(double*)*order);
  //std::vector<double> wts( order, 0. );
  ccdl::bspline_eval( ( order%2 == 0 ? 0.0 : 0.5 ), order, wts );
  ccdl::bspline_hdct( N, order, wts, fc );
}



///////////////////////////////////////////////////////////////////
// renormalizes a uniform grid of data so that b-spline
// interpolation exactly passes through the data
///////////////////////////////////////////////////////////////////
void ccdl::bspline_renormalize
( double const lx,
  int const nx,
  int const order,
  double * data )
{
  int const nfourier      = nx/2+1;
  
  //
  // allocate fft data and plans
  //
  double * value = fftw_alloc_real( nx );

  std::complex<double> * fourier 
    = reinterpret_cast< std::complex<double> * >
    ( fftw_alloc_complex( nfourier ) );

  fftw_plan * fplan = new fftw_plan
    ( fftw_plan_dft_r2c_1d
      ( nx,
        value, 
        reinterpret_cast< fftw_complex *>( fourier ), 
        FFTW_ESTIMATE ) );

  fftw_plan * rplan = new fftw_plan
    ( fftw_plan_dft_c2r_1d
      ( nx,
        reinterpret_cast< fftw_complex *>( fourier ), 
        value, 
        FFTW_ESTIMATE ) );
  //
  // forward fft of data
  //
  for ( int i=0; i<nx; ++i )
    value[i] = data[i];
  fftw_execute( *fplan );
  double const volElement = lx/nx;
  for ( int i=0; i<nfourier; ++i )
    fourier[i] *= volElement;
  //
  // forward fft of b-splines
  //
  std::vector<double> fx( nx, 0. );
  bspline_fourier_half( nx, order, fx.data() );
  //
  // scale the fourier coefs
  //
  for ( int i=0; i < nfourier; ++i )
    fourier[i] /= ( fx[i] );
  //
  // reverse fft
  //
  fftw_execute( *rplan ); 
  double const ooV = 1./lx;
  for ( int i=0; i<nx; ++i )
    data[i] = ooV * value[i];
  //
  // delete
  //
  fftw_destroy_plan( *fplan );
  delete fplan;
  fplan = NULL;
  fftw_destroy_plan( *rplan );
  delete rplan;
  rplan = NULL;
  fftw_free(reinterpret_cast< fftw_complex * &>(fourier));
  fourier = NULL;
  fftw_free( value );
  value = NULL;
}

