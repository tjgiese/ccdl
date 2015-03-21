#ifndef _CCDL_GieseBoysFcn_H_
#define _CCDL_GieseBoysFcn_H_

#include <cmath>

extern const double GlobalGieseBoysFcnData[];

namespace ccdl
{
  void GieseBoysFcn( int const Lin, double const T, double *__restrict__ Fm );

  template <int Lin>
  void GieseBoysFcnTemplate( double const T, double *__restrict__ Fm );
}





template <int Lin>
void ccdl::GieseBoysFcnTemplate
( double const T, 
  double *__restrict__ Fm ) 
{
  const double Rupper = 30.;
  const double FmDEL = 6.250000000000000e-02;
  const int FmM = 18;
  double const PI = 3.141592653589793238462643383279502884197;
  double const EPS = 1.e-13;
  int const MAXFACT = 300;
  double const DELTAS[] = { FmDEL, 2.*FmDEL, 4*FmDEL };
  int const IOFFS[] = { 0, 80, 160 };
 
  //int const L = Lin+1;


  if ( T >= Rupper ) 
    {
      
      // F0(T) ~ (1/2) SQRT(PI/T)
      // Fn(T) ~ (1/2) SQRT(PI/T) (2n-1)!!/2^n (1/T)^n
      
      double df = 0.5;
      double oot = 1.;
      double const ooT = 1./T;
      double const fact = sqrt( PI * ooT );
      for ( int m=0; m<Lin+1; ++m )
	{
	  Fm[m] = df * fact * oot;
	  df   *= (m+0.5);
	  oot  *= ooT;
	};

    }
  else if ( T < 1.e-4 )
    {
      // Fn(T) ~ 1/(2n+1)
      for ( int m=0; m < Lin+1; ++m )
	Fm[m] = 1. / (2.*m+1.);

    }
  else if ( Lin+6 <= FmM )
    {

      // 6 TERM TAYLOR EXPANSION

      // 3 ranges:
      //  0-10  ; 160 pts ; delta =     0.0625
      // 10-20  ;  80 pts ; delta = 2 * 0.0625
      // 20-30  ;  40 pts ; delta = 4 * 0.0625

      // What range are we in? (0-2)
      int const nr = (int)(T/10.);

      // What is the delta for this range?
      double const del = DELTAS[nr];
      // What is the dT from the nearest point?
      int ipt = (int)( T/del + 0.5 );
      double const dT = T-ipt*del;
      // We actually have more data pretabulated than the delta implies
      // Add an extra offset to account for the extra data in the
      // previous regions
      //ipt = ipt + IOFFS[nr];
      // Get the offset in the master data array for this point
      //int const o = Lin + ipt * FmM;

      // 6 TERM TAYLOR EXPANSION
      double const dT2 = dT*dT;
      double const dT3 = dT2*dT;
      double const dT4 = dT2*dT2;
      double const dT5 = dT2*dT3;
      /*
	Fm[Lin]
	= GlobalGieseBoysFcnData[o]                   
	- GlobalGieseBoysFcnData[o+1] * dT            
	+ GlobalGieseBoysFcnData[o+2] * dT2 / 2. 
	- GlobalGieseBoysFcnData[o+3] * dT3 / 6.
	+ GlobalGieseBoysFcnData[o+4] * dT4 / 24.
	- GlobalGieseBoysFcnData[o+5] * dT5 / 120.;
      */

      double const * f = GlobalGieseBoysFcnData + Lin + (ipt+IOFFS[nr]) * FmM;
      //double const * f = &GlobalGieseBoysFcnData[Lin + (ipt+IOFFS[nr]) * FmM];
      //printf("%p\n",(void*)f);
      Fm[Lin]  = *f;
      Fm[Lin] -= *(++f) * dT;
      Fm[Lin] += *(++f) * dT2 / 2.;
      Fm[Lin] -= *(++f) * dT3 / 6.;
      Fm[Lin] += *(++f) * dT4 / 24.;
      Fm[Lin] -= *(++f) * dT5 / 120.;


      // DOWNWARD RECURSION
      double const TwoT = 2. * T;
      double const ExpT = std::exp(-T);
      //for ( int m=Lin-1; m >= 0; --m )
      //Fm[m]=(TwoT*Fm[m+1]+ExpT)/(2.*m+1.);
      
      for ( int m=Lin; m-- > 0; )
	Fm[m]=(TwoT*Fm[m+1]+ExpT)/(2*m+1);
      
    }
  else
    {

      // Fn(T) = Exp(-T) ( Exp(T) * Fn(T) )
      //       ~ Exp(-T) PowerExpansion( Exp(T) * Fn(T) )
      //       ~ Exp(-T) ( sum_i (2n-1)!! (2T)^i / (2n+2i+1)!! )

      double const ExpT  = std::exp(-T);
      double const TwoT  = 2. * T;
       
      // num = (2L-1)!!
      double num = 1.;
      for ( int i=2; i < Lin+1; ++i )
	num *= ( 2.*i-1. );
      
      
      double dfi = 2*Lin+1;
      
      // df = (2L+1)!!
      double df = num * dfi;
      
      double sum = 1.0/dfi;
      double term1 = 0.0;
      for ( int i=0; i < MAXFACT; ++i )
	{
	  dfi   = dfi+2;
	  df    = df*dfi;
	  num   = num*TwoT;
	  term1 = num/df;
	  sum   = sum+term1;
	  if ( (term1 < EPS) || (num > 1.e+290) ) break;
	};

      Fm[Lin] = sum*ExpT;
      for ( int m=Lin; m-- > 0; )
	Fm[m]=(TwoT*Fm[m+1]+ExpT)/(2*m+1);

    };

}



#endif
