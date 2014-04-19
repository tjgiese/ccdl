#include "bspline.hpp"

#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>

#ifndef alloca
#define alloca __builtin_alloca
#endif


///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


void ccdl::one_pass_bspline( double * c, double const w, int const n )
{
  int nm1 = n-1;
  double const div = 1. / nm1;
  c[nm1] = div * w * c[nm1 - 1];
  for ( int j=1; j<nm1; ++j )
    c[nm1 - j] = div * ((w + j) * c[nm1 - j - 1] + (n - j - w) * c[nm1 - j]);
  c[0] = div * (1 - w) * c[0];
}
void ccdl::diff_bspline( int const order, double const * array, double * diff )
{
  assert( order > 1 );
  int const nm1 = order-1;
  diff[0] = -array[0];
  for ( int j=1; j<nm1; ++j )
    diff[j] = array[j-1] - array[j];
  diff[nm1] = array[nm1-1];
}

void ccdl::fill_bspline_0( double const w, int const order, double * array )
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
            ccdl::one_pass_bspline(array, w, k);
        };
    };
}

void ccdl::fill_bspline_1( double const w, int const order, double * array, double * darray )
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
        ccdl::one_pass_bspline(array, w, k);

      ccdl::diff_bspline(order,array,darray);

      // One more recursion: // do the k==order
      ccdl::one_pass_bspline(array, w, order);

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

void ccdl::fill_bspline_n( double const w, int const order, int const nder, double * array )
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
    ccdl::diff_bspline(order,darray[d],darray[d+1]);
  array[0] = 1. - w;
  array[1] = w;
  for ( int o=3; o<=order; ++o )
    {
      for ( int d=0, nd=std::min(order-(o-1),nder); d<nd; ++d )
        diff_bspline(order,darray[d],darray[d+1]);
      ccdl::one_pass_bspline(array, w, o);
    };
}



///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


void ccdl::periodicBspline
( double x, 
  double const lenx, 
  int const nx,
  int const order,
  double * w )
{
  x = x*nx/lenx - (order%2)*0.5;
  ccdl::fill_bspline_0(x-std::floor(x),order,w);
}


void ccdl::periodicBspline
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
  ccdl::fill_bspline_0(x,order,w);
}


void ccdl::periodicBspline
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
  ccdl::fill_bspline_0(q,order,w);
}





void ccdl::periodicBsplineDeriv
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
  ccdl::fill_bspline_1(x,order,w,dw);
  x = nx/lenx;
  for ( int b=0; b<order; ++b )
    dw[b] *= x;
}


void ccdl::periodicBsplineDeriv
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
  ccdl::fill_bspline_1(x,order,w,dw);
  x = nx/lenx;
  for ( int b=0; b<order; ++b, ++ilo )
    {
      gidx[b] = ilo%nx;
      dw[b] *= x;
    };
}


void ccdl::periodicBsplineDeriv
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
  ccdl::fill_bspline_1(q-ilo,order,w,dw);
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


void ccdl::aperiodicBspline
( double x, 
  double const lenx, 
  int const nx,
  int const order, 
  double * w,
  int * gidx )
{
  x = x*nx/lenx - (order%2)*0.5;
  int ilo = std::floor(x);
  ccdl::fill_bspline_0(x-ilo,order,w);
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

void ccdl::aperiodicBspline
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
  ccdl::fill_bspline_0(q-ilo,order,w);
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

void ccdl::aperiodicBsplineDeriv
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
  ccdl::fill_bspline_1(x-ilo,order,w,dw);
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



void ccdl::aperiodicBsplineDeriv
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
  ccdl::fill_bspline_1(q-ilo,order,w,dw);
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








void ccdl::periodicBsplineNthDeriv
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
  ccdl::fill_bspline_n(q-ilo,order,nder,dwdx);
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


void ccdl::aperiodicBsplineNthDeriv
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
  ccdl::fill_bspline_n(q-ilo,order,nder,dwdx);
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

