#include <cmath>
#include <vector>
#include "NumericalBasis.hpp"
#include "OrthogonalPolynomials.hpp"
#include "../constants.hpp"








//  FUNCTION SHGaussianL2RadialNorm(z,l) RESULT(norm) 
//****f* NumericalBasisMod/SHGaussianL2RadialNorm
//
// NAME
//     SHGaussianL2RadialNorm - Spherical Harmonic Gaussian Normalization Constant
//
// USAGE
//     norm = SHGaussianL2RadialNorm(z,l)
//
// DESCRIPTION
//     Given a spherical harmonic Gaussian function of the form
//      G(r,th,phi) = N * r**l * std::exp(-z*r**2) * Ylm(th,phi),
//      (Note that this is a solid harmonic, actually - people refer to it, incorrectly, as
//       a spherical harmonic Gaussian function.)
//       This routine computes the constant "N" needed for L2 normalization.
//       Specifically, the result of this routine is
//       result = SQRT(1/blah),
//       where
//       blah = INT_0^inf  (r^l std::exp(-z*r**2))*(r^l std::exp(-z*r**2)) r^2 dr
//       Note that only the radial portion is integrated because the spherical harmonics
//       defined in this module ALREADY have the appropriate angular normalization.
//
// INPUTS
//     REAL(SP),INTENT(IN) :: z
//     INTEGER(I4B),INTENT(IN) :: l
//
// OUTPUTS
//     double norm ! FUNCTION RESULT
//
//***
    // THIS FUNCTION RETURNS THE SQUARE ROOT INVERSE OF THE
    // 1-DIMENSIONAL INTEGRAL
    // blah = INT_0^inf  r^2  (r^l std::exp(-z*r**2))*(r^l std::exp(-z*r**2)) dr
    // result = SQRT(1/blah)
    //
    // Note that
    // blah = N[l] = Gamma[3/2+l] * z**(-3/2-l) * 2**(-5/2-l)
    //      = Gamma[3/2+l] * z**(-3/2) * z**(-l) * 2**(-5/2) * 2**(-l)
    //      = Gamma[3/2+l] * (2*(2z)**(-3/2)) * (2z)**(-l)
    // N[0]      = SQRT(PI/(2*z**3))/8
    // N[1]/N[0] = 3/(4z) ; N[1] = (N[1]/N[0])*N[0]                         = N[0] * 3/(4z)
    // N[2]/N[1] = 5/(4z) ; N[2] = (N[2]/N[1])*(N[1]/N[0])*N[0]             = N[0] * 3/(4z) * 5/(4z)
    // N[3]/N[2] = 7/(4z) ; N[3] = (N[3]/N[2])*(N[2]/N[1])*(N[1]/N[0])*N[0] = N[0] * 3/(4z) * 5/(4z) * 7/(4z)
    // ...
    //
    // THIS IS USEFUL if YOU ARE USING SPHERICAL HARMONIC
    // GAUSSIAN FUNCTIONS OF THE FORM
    // r^l std::exp(-Z*r**2) Ylm(theta,phi)
    // Because CheapSphericalHarmonic is really fast and
    // returns Ylm(theta,phi) such that
    // INT INT Ylm(theta,phi)*Ylm(theta,phi) Cos(theta) dtheta dphi = 1
    // So, that part of the normalization is taken care of by that
    // routine.  The only other part of the normalization to consider
    // is the 1-dimensional radial integral.
    //
double SHGaussianL2RadialNorm
( double const z,
  int const l )
{
  double norm;
  double fourz;
  int cnt;

  norm  = std::sqrt(ccdl::PI/(2.0*z*z*z))/8.0;
  fourz = 4.0 * z;
  cnt = 1;
  for ( int i=1; i <= l; ++i )
    {
      cnt+=2;
      norm *= cnt/fourz;
    };
  return std::sqrt(1.0/norm);
}



//  FUNCTION SHGaussianL1RadialNorm(z,l) RESULT(norm) 
//****f* NumericalBasisMod/SHGaussianL2RadialNorm
//
// NAME
//     SHGaussianL1RadialNorm - Spherical Harmonic Gaussian Normalization Constant
//
// USAGE
//     norm = SHGaussianL1RadialNorm(z,l)
//
// DESCRIPTION
//
// INPUTS
//     REAL(SP),INTENT(IN) :: z
//     INTEGER(I4B),INTENT(IN) :: l
//
// OUTPUTS
//     double norm ! FUNCTION RESULT
//
//***
    //
    //   result = 1/blah
    //   blah   = Int[ r**2 r**L (4*PI/(2*L+1)) (2*z*r)**L * EXP(-z*r**2) , {r,0,inf} ]
    //   blah = N[L] = (PI*2**L)/(z**1.5)  * Gamma[L+1/2]
    //   N[0] = (PI/z)**1.5
    //   N[1] = N[0] * (N[1]/N[0]) = N[0] * 1
    //   N[2] = N[0] * (N[1]/N[0]) * (N[2]/N[1]) = N[0] * 1 * 3
    //   N[3] = N[0] * (N[1]/N[0]) * (N[2]/N[1]) * (N[3]/N[2]) = N[0] * 1 * 3 * 5
    //   ...
    //   N[L] = N[0] * (2*L-1)////  
double SHGaussianL1RadialNorm
( double const z,
  int const l )
{
  double norm;
  int cnt;
  norm  = 1.0;
  cnt = -1;
  for ( int i=1; i<=l; ++i )
    {
      cnt += 2;
      norm *= cnt;
    };
  norm = std::pow( (z/ccdl::PI), 1.5 ) * std::sqrt(ccdl::FOUR_PI/(2*l+1)) / norm;
  if ( l > 0 ) norm *= std::pow( (2*z),l );
  return norm;
}





//    Radial Hermite basis function
//       Rh(z,n,x) = Nc * H(n,z*x) * EXP(-0.5 * (z*x)**2)
//    where Nc = SQRT( z / SQRT_PI * n! * 2**n )
//
//    Note that for a perfect simple harmonic oscillator,
//    this function returns the exact wavefunction when
//    z = SQRT( SQRT(m*k) / hbar )
//    where m is the reduced mass (in atomic units),
//    k is the force constant (k = d^2 V(x)/dx^2 in a.u.)
//    and hbar is the reduced Planck's constant (which is unity in a.u.).
//    Recall that a SHO has a potential of
//      V(x) = 1/2 * k * x^2
//
//    This basis holds the property
//       \int_{-inf}^{inf} Rh(z,n,x) Rh(z,n',x) dx = delta_{n,n'}
//
//    The derivatives are
//       D[Rh(z,n,x),x]      = Nc*z*( 2*n*H(n-1,z*x)-z*x*H(n,z*x) ) * EXP(-0.5 * (z*x)**2)
//       D[D[Rh(z,n,x),x],x] = Nc*z**2*( 4*(n-1)*n*H(n-2,z*x) - 4*n*z*x*H(n-1,z*x) + ((z*x)**2-1)*H(n,z*x) )
//

  // Radial Hermite basis function
  //    Rh(z,n,x) = Nc * H(n,z*x) * EXP(-0.5 * (z*x)**2)
  // where Nc = SQRT( z / SQRT_PI * n! * 2**n )
  // It holds the property
  //    Int_{-inf}^{inf} Rh(z,n,x) Rh(z,n',x) dx = delta_{n,n'}
  // The derivatives are
  //    D[Rh(z,n,x),x]      = Nc*z*( 2*n*H(n-1,z*x)-z*x*H(n,z*x) ) * EXP(-0.5 * (z*x)**2)
  //    D[D[Rh(z,n,x),x],x] = Nc*z**2*( 4*(n-1)*n*H(n-2,z*x) - 4*n*z*x*H(n-1,z*x) + ((z*x)**2-1)*H(n,z*x) )

inline void fill_factrl( int const n, double * a )
{
  a[0] = 1.;
  for ( int i=1; i<n; ++i )
    a[i] *= i*a[i-1];
}

inline double factorial( int const n )
{
  double a = 1.;
  for ( int i=2; i<n; ++i )
    a *= i;
  return a;
}

//SUBROUTINE l2_hermite_basis_3(z,npoly,npts,RadPts,f,fp,fpp,fppp) BIND(C, name="l2_hermite_basis_3_")
  
void ccdl::HermiteBasis
( double const z,
  int const npoly,
  int const npts,
  double const * RadPts,
  double *  f,
  double *  fp,
  double *  fpp,
  double *  fppp )
{

  std::vector<double> H(npoly+1,0.);
  double H0,H1,H2,H3,r,c,zz,zr,zzr,zzrr,e;
  double dedr,d2edr2,d3edr3;
  double dHdr,d2Hdr2,d3Hdr3;

  c = std::sqrt( z / ( ccdl::SQRT_PI * factorial(npoly) * std::pow( 2.0, npoly ) ) );
  zz = z*z;
  
  for ( int i=0; i<npts; ++i )
    {
      r    = RadPts[i];
      zr   = z*r;
      zzr  = z*zr;
      zzrr = zr*zr;
      e = c*std::exp(-zzrr*0.5);

      ccdl::HermitePolynomial(zr,npoly+1,H.data());
      H0 = H.back();
      H1 = 0.;
      H2 = 0.;
      H3 = 0.;
      if ( npoly > 0 ) H1 = H[npoly-2];
      if ( npoly > 1 ) H2 = H[npoly-3];
      if ( npoly > 2 ) H3 = H[npoly-4];

      dedr   = -zzr*e;
      d2edr2 =  -zz*e   - zzr*dedr;
      d3edr3 = -2*zz*dedr - zzr*d2edr2;
      dHdr   = 2*z*npoly* H1;
      d2Hdr2 = 2*z*npoly*( 2*z*(npoly-1) )* H2;
      d3Hdr3 = 2*z*npoly*( 2*z*(npoly-1) )*( 2*z*(npoly-2) )* H3;

      f[i]    = e*H0;
      fp[i]   = dedr*H0 + e*dHdr;
      fpp[i]  = d2edr2*H0 + 2*dedr*dHdr + e*d2Hdr2;
      fppp[i] = d3edr3*H0 + 3*d2edr2*dHdr + 3*dedr*d2Hdr2 + e*d3Hdr3;

    };
}





  // The basis function is
  //   B(r;n,l,m,z) = f(r;n,l,z) * Nc(n,l,z) * Clm(r)
  //   for  n >=0 and L >=0
  //
  // where 
  //   Nc(n,l,z) = SQRT( (2*z)**3  * n! / (n+2*l+2)! ) / SQRT( FOUR_PI / (2*L+1) )
  // and
  //
  // f(r)  =  (2*z)**L * std::exp(-z*r) * Lag(n;2*L+2;2*z*r)
  // f'(r) = -(2*z)**L * std::exp(-z*r) * 
  //           (   z   * Lag(n;2*L+2;2*z*r) 
  //           +   2*z * Lag(n-1;2*L+3;2*z*r)  )
  // f"(r) =  (2*z)**L * std::exp(-z*r) *
  //           (   z*z * Lag(n;2*L+2;2*z*r) 
  //           + 4*z*z * Lag(n-1;2*L+3;2*z*r)
  //           + 4*z*z * Lag(n-2;2*L+4;2*z*r)  )
  //
  // Lag(n;m;0) = (n+m)! / ( n! * m! ) ; if n >= 0
  //              0                    ; otherwise
  //
  //
  // It holds the properties
  // \int B(r;n,l,m,z) * B(r;n',l',m',z) d3r = \delta_{n,n'} \delta_{l,l'} \delta_{m,m'}
  // \int f(r;n,l,z)   * Nc(n,l,z)    *  
  //      f(r;n',l',z) * Nc(n',l',z)  * 
  //      ( FOUR_PI / (2*L+1) )       * r**2 dr = \delta_{n,n'} \delta_{l,l'}



  // // Notice that
  // //    n L(n,m,x) = (n+m) L(n-1,m,x) - x L(n-1,m+1,x)
  // // or alternatively
  // //    L(n-1,m+1,x) = ( (n+m) L(n-1,m,x) - n L(n,m,x) ) / x
  // // Also notice that we need
  // //    L(n,2l+2,x), L(n-1,2l+3,x) and L(n-2,2l+4,x)
  // // We have all n for L(:,2l+2,x), so
  // //    L(n-1,2l+3,x) = ( (n+2l+2) L(n-1,2l+2,x) - n     * L(n,2l+2,x) ) / x
  // //    L(n-2,2l+4,x) = ( (n+2l+1) L(n-2,2l+3,x) - (n-1) * L(n-1,2l+3,x) ) / x
  // // where
  // //    L(n-2,2l+3,x) = ( (n+2l+1) L(n-2,2l+2,x) - (n-1) * L(n-1,2l+2,x) ) / x

  // IF ( n == 0 ) THEN
  //    L0 = 1.
  //    L1 = 0.
  //    L2 = 0.
  // ELSE IF ( n == 1 ) THEN
  //    L0 = -tzx + (2*l+2) + 1.
  //    L1 = 1.
  //    L2 = 0.
  // ELSE
  //    CALL laguerre_polynomial( tzx, 2*l+2, n+1, Lag(0) )
  //    L0  = Lag(n)
  //    L1  = ( (n+2*l+2) * Lag(n-1) - n     * L(n)   ) / tzx
  //    L11 = ( (n+2*l+1) * Lag(n-2) - (n-1) * L(n-1) ) / tzx
  //    L2  = ( (n+2*l+1) * L11      - (n-1) * L1     ) / tzx
  // END IF

void ccdl::ClmLaguerreBasis
( double const z,
  int const n,
  int const l,
  int const npt,
  double const * xvec,
  double * f,
  double * fp,
  double * fpp,
  double * fppp )
{
  std::vector<double> factrl(n+2*l+2,1.);
  std::vector<double> Lag(n+1,0.);
  double Nc,x,tz,tzl,tzx;
  double L0,L1,L2,L3;
  int tlt;
  double dLdx, d2Ldx2, d3Ldx3, dpdx, d2pdx2, d3pdx3, p;

  tz  = 2*z;
  tzl = std::pow(tz,l);
  tlt = 2*l+2;
  fill_factrl( n+tlt+1, factrl.data() );

  Nc = std::sqrt( std::pow(tz,3) * factrl[n] / factrl[n+tlt] ) 
    / std::sqrt( ccdl::FOUR_PI / (2*l+1) );

  for ( int i=0; i<npt; ++i )
    {
      x = xvec[i];

      if ( x <= 1.E-25 )
	{
          L0 = 0.;
          L1 = 0.;
          L2 = 0.;
          L3 = 0.;

          L0 = factrl[n+tlt] / ( factrl[n] * factrl[tlt] );
          if ( n > 0 ) L1 = factrl[n+tlt] / ( factrl[n-1] * factrl[tlt+1] );
          if ( n > 1 ) L2 = factrl[n+tlt] / ( factrl[n-2] * factrl[tlt+2] );
          if ( n > 2 ) L3 = factrl[n+tlt] / ( factrl[n-3] * factrl[tlt+3] );

          p = Nc * tzl;
          dpdx   = -z*p;
          d2pdx2 = -z*dpdx;
          d3pdx3 = -z*d2pdx2;

          dLdx   = (-tz)*L1;
          d2Ldx2 = (-tz)*(-tz)*L2;
          d3Ldx3 = (-tz)*(-tz)*(-tz)*L3;

          f[i]    = p*L0;
          fp[i]   = dpdx*L0 + p*dLdx;
          fpp[i]  = d2pdx2*L0 + 2*dpdx*dLdx + p*d2Ldx2;
          fppp[i] = d3pdx3*L0 + 3*d2pdx2*dLdx + 3*dpdx*d2Ldx2 + p*d3Ldx3;
	}
      else
	{
          tzx = tz*x;
          p   = Nc * tzl * std::exp(-z*x);

          L0 = 0.;
          L1 = 0.;
          L2 = 0.;
          L3 = 0.;

	  ccdl::LaguerrePolynomial( tzx, tlt, n+1, Lag.data() );
          L0 = Lag[n];
          if ( n > 0 )
	    {
	      ccdl::LaguerrePolynomial( tzx, tlt+1, n, Lag.data() );
	      L1 = Lag[n-1];
	    }
          if ( n > 1 )
	    {
	      ccdl::LaguerrePolynomial( tzx, tlt+2, n-1, Lag.data() );
	      L2 = Lag[n-2];
	    }
          if ( n > 2 )
	    {
	      ccdl::LaguerrePolynomial( tzx, tlt+3, n-2, Lag.data() );
	      L3 = Lag[n-3];
	    };

          dpdx   = -z*p;
          d2pdx2 = -z*dpdx;
          d3pdx3 = -z*d2pdx2;

          dLdx   = (-tz)*L1;
          d2Ldx2 = (-tz)*(-tz)*L2;
          d3Ldx3 = (-tz)*(-tz)*(-tz)*L3;

          f[i]    = p*L0;
          fp[i]   = dpdx*L0 + p*dLdx;
          fpp[i]  = d2pdx2*L0 + 2*dpdx*dLdx + p*d2Ldx2;
          fppp[i] = d3pdx3*L0 + 3*d2pdx2*dLdx + 3*dpdx*d2Ldx2 + p*d3Ldx3;
	};
    };
}



 
void ccdl::YlmLaguerreBasis
( double const z,
  int const n,
  int const l,
  int const npt,
  double const * xvec,
  double * f,
  double * fp,
  double * fpp,
  double * fppp )
{
  std::vector<double> factrl( n+2*l+5, 1. );
  std::vector<double> Lag(n+1,0.);
  double Nc,x,tz,tzl,tzx,L0,L1,L2,L3;
  int tlt;

  double X0,X1,X2,X3;
  double A,dAdx,d2Adx2,d3Adx3;
  double p,dpdx,d2pdx2,d3pdx3;
  double dLdx,d2Ldx2,d3Ldx3;

  tz  = 2*z;
  tzl = std::pow(tz,l);
  tlt = 2*l+2;

  fill_factrl( n+tlt+1, factrl.data() );

  Nc = std::sqrt( std::pow(tz,3)  * factrl[n] / factrl[n+tlt] );

  for ( int i=0; i<npt; ++i )
    {
      x = xvec[i];
      if ( x <= 1.E-25 ) 
	{
	  L0 = 0.;
	  L1 = 0.;
	  L2 = 0.;
          L3 = 0.;

          L0 = factrl[n+tlt] / ( factrl[n] * factrl[tlt] );
          if ( n > 0 ) L1 = factrl[n+tlt] / ( factrl[n-1] * factrl[tlt+1] );
          if ( n > 1 ) L2 = factrl[n+tlt] / ( factrl[n-2] * factrl[tlt+2] );
          if ( n > 2 ) L3 = factrl[n+tlt] / ( factrl[n-3] * factrl[tlt+3] );

          X0=0.;
          X1=0.;
          X2=0.;
          X3=0.;

          if ( l == 0 ) X0 = 1.;
          if ( l == 1 ) X1 = l;
          if ( l == 2 ) X2 = l*(l-1);
          if ( l == 3 ) X3 = l*(l-1)*(l-2);

          p      = Nc * tzl;
          dpdx   = -z*p;
          d2pdx2 = -z*dpdx;
          d3pdx3 = -z*d2pdx2;

          A      = p*X0;
          dAdx   = dpdx*X0 + p*X1;
          d2Adx2 = d2pdx2*X0 + 2*dpdx*X1 + p*X2;
          d3Adx3 = d3pdx3*X0 + 3*d2pdx2*X1 + 3*dpdx*X2 + p*X3;

          dLdx   = (-tz)*L1;
          d2Ldx2 = (-tz)*(-tz)*L2;
          d3Ldx3 = (-tz)*(-tz)*(-tz)*L3;

          f[i]    = A*L0;
          fp[i]   = dAdx*L0 + A*dLdx;
          fpp[i]  = d2Adx2*L0 + 2*dAdx*dLdx + A*d2Ldx2;
          fppp[i] = d3Adx3*L0 + 3*d2Adx2*dLdx + 3*dAdx*d2Ldx2 + A*d3Ldx3;
	}
      else
	{
          tzx = tz*x;

          L0 = 0.;
          L1 = 0.;
          L2 = 0.;
          L3 = 0.;

	  ccdl::LaguerrePolynomial( tzx, tlt, n+1, Lag.data() );
          L0 = Lag[n];
          if ( n > 0 )
	    { 
	      ccdl::LaguerrePolynomial( tzx, tlt+2, n, Lag.data() );
	      L1 = Lag[n-1];
	    };
          if ( n > 1 )
	    { 
	      ccdl::LaguerrePolynomial( tzx, tlt+3, n-1, Lag.data() );
	      L2 = Lag[n-2];
	    };
          if ( n > 2 )
	    { 
	      ccdl::LaguerrePolynomial( tzx, tlt+4, n-2, Lag.data() );
	      L3 = Lag[n-3];
	    };
	  
          X0 = std::pow(x,l);
          X1 = l*std::pow(x,(l-1));
          X2 = l*(l-1)*std::pow(x,(l-2));
          X3 = l*(l-1)*(l-2)*std::pow(x,(l-3));
	  
          p      = Nc * tzl * std::exp(-z*x);
          dpdx   = -z*p;
          d2pdx2 = -z*dpdx;
          d3pdx3 = -z*d2pdx2;
	  
          A      = p*X0;
          dAdx   = dpdx*X0 + p*X1;
          d2Adx2 = d2pdx2*X0 + 2*dpdx*X1 + p*X2;
          d3Adx3 = d3pdx3*X0 + 3*d2pdx2*X1 + 3*dpdx*X2 + p*X3;

          dLdx   = (-tz)*L1;
          d2Ldx2 = (-tz)*(-tz)*L2;
          d3Ldx3 = (-tz)*(-tz)*(-tz)*L3;


          f[i]    = A*L0;
          fp[i]   = dAdx*L0 + A*dLdx;
          fpp[i]  = d2Adx2*L0 + 2*dAdx*dLdx + A*d2Ldx2;
          fppp[i] = d3Adx3*L0 + 3*d2Adx2*dLdx + 3*dAdx*d2Ldx2 + A*d3Ldx3;
	};
    };
}







void ccdl::ClmGaussianBasis
( double const z,
  int const l,
  int const npt,
  double const * XVec,
  double * f,
  double * fp,
  double * fpp,
  double * fppp )
{
  double Nc,x;
  double e,dedx,d2edx2,d3edx3;

  Nc = SHGaussianL2RadialNorm(z,l) / std::sqrt( ccdl::FOUR_PI/(2*l+1) );

  for ( int i=0; i<npt; ++i )
    {
      x = XVec[i];
       if ( x < 1.E-25 )
	 { 
	   f[i]    = Nc;
	   fp[i]   = 0.;
	   fpp[i]  = -2*z*Nc;
	   fppp[i] = 0.;
	 }
       else
	 {
	   e      = Nc * std::exp(-z*x*x);
	   dedx   = -2*z*x*e;
	   d2edx2 = -2*z*e - 2*z*x*dedx;
	   d3edx3 = -4*z*dedx - 2*z*x*d2edx2;
	   f[i]    = e;
	   fp[i]   = dedx;
	   fpp[i]  = d2edx2;
	   fppp[i] = d3edx3;
	 };
    };
}







void ccdl::YlmGaussianBasis
( double const z,
  int const l,
  int const npt,
  double const * XVec,
  double * f,
  double * fp,
  double * fpp,
  double * fppp )
{

  double Nc,x;
  double e,dedx,d2edx2,d3edx3;
  double A,dAdx,d2Adx2,d3Adx3;
  double X0,X1,X2,X3;

//    double SHGaussianL2RadialNorm

  Nc = SHGaussianL2RadialNorm(z,l);

  for ( int i=0; i<npt; ++i )
    {
      x = XVec[i];
      if ( x < 1.E-25 )
	{ 
          X0 = 0.;
          X1 = 0.;
          X2 = 0.;
          X3 = 0.;
          if ( l == 0 ) X0 = 1.;
          if ( l == 1 ) X1 = l;
          if ( l == 2 ) X2 = l*(l-1);
          if ( l == 3 ) X3 = l*(l-1)*(l-2);
          e      = Nc;
          dedx   = 0.;
          d2edx2 = -2*z*e ;
          d3edx3 = 0.;

          A      = e*X0;
          dAdx   = dedx*X0 + e*X1;
          d2Adx2 = d2edx2*X0 + 2*dedx*X1 + e*X2;
          d3Adx3 = d3edx3*X0 + 3*d2edx2*X1 + 3*dedx*X2 + e*X3;

          f[i]    = A;
          fp[i]   = dAdx;
          fpp[i]  = d2Adx2;
          fppp[i] = d3Adx3;
	}
      else
	{
          X0 = std::pow(x,l);
	  X1 = l*std::pow(x,(l-1));
          X2 = l*(l-1)*std::pow(x,(l-2));
          X3 = l*(l-1)*(l-2)*std::pow(x,(l-3));

          e      = Nc * std::exp(-z*x*x);
          dedx   = -2*z*x*e;
          d2edx2 = -2*z*e - 2*z*x*dedx;
          d3edx3 = -4*z*dedx - 2*z*x*d2edx2;

          A      = e*X0;
          dAdx   = dedx*X0 + e*X1;
          d2Adx2 = d2edx2*X0 + 2*dedx*X1 + e*X2;
          d3Adx3 = d3edx3*X0 + 3*d2edx2*X1 + 3*dedx*X2 + e*X3;

          f[i]    = A;
          fp[i]   = dAdx;
          fpp[i]  = d2Adx2;
          fppp[i] = d3Adx3;
	};
    };
}

















