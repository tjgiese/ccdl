#ifndef _ccdl_wignerd_hpp_
#define _ccdl_wignerd_hpp_


namespace ccdl
{
  //
  // The little d Wigner matrix
  //    d^j_{m,k} = d[ (j+m) + (j+k)*(2*j+1) ]
  // where 
  //   -j <= m <= +j 
  //   -j <= k <= +j
  //
  // That is, the "zero"-component is stored in the middle
  //
  // This matrix can be reproduced in mathematica by
  //   WignerD[ {j, m, k}, 0, \[Beta], 0 ]
  //
  void WignerLittleD( int j, double Y, double * d );

  //
  // The real and imaginary components of the complex-valued
  // Wigner-D matrix
  //    D^j_{m,k} = e^{-i m Zp} d^j_{m,k}(Y) e^{-i k Z}
  //
  // The real values are stored in c
  // The imaginary values are stored in s
  //
  // Specifically,
  //   c[ (j+m) + (j+k)*(2*j+1) ] = Re[ D^j_{m,k} ]
  //   s[ (j+m) + (j+k)*(2*j+1) ] = Im[ D^j_{m,k} ]
  //
  // That is, the "zero"-component is stored in the middle
  //
  // These matrices can be reproduced in mathematica by
  //   c=ComplexExpand[ Re[ WignerD[ {j, m, k}, \[Alpha], \[Beta], \[Gamma] ]]]
  //   s=ComplexExpand[ Im[ WignerD[ {j, m, k}, \[Alpha], \[Beta], \[Gamma] ]]]
  //
  void WignerBigD( int j, double Z, double Y, double Zp, 
		   double * c, double * s );

  // The real-valued Wigner-D matrix used to rotate
  // real-valued spherical harmonics
  //
  // The matrix has the form:
  //   D^{cc}_{00} D^{cc}_{01} D^{cs}_{01} D^{cc}_{02} D^{cs}_{02} ...
  //   D^{cc}_{10} D^{cc}_{11} D^{cs}_{11} D^{cc}_{12} D^{cs}_{12} ...
  //   D^{sc}_{10} D^{sc}_{11} D^{ss}_{11} D^{sc}_{12} D^{ss}_{12} ...
  //   D^{cc}_{20} D^{cc}_{21} D^{cs}_{21} D^{cc}_{22} D^{cs}_{22} ...
  //   D^{sc}_{20} D^{sc}_{21} D^{ss}_{21} D^{sc}_{22} D^{ss}_{22} ...
  //
  // That is, the "zero"-component is stored as the first element.
  //
  // The Euler angles Z, Y, Zp correspond to the rotation
  //   r' = R(Zp).R(Y).R(Z).r
  //   r' = R.r
  //
  // Y(R.r) = D(Z,Y,Zp) . Y(r)
  //
  // D_{lm,jk}(R) = \int  Y_{lm}( R.r ) Y_{jk}( r ) dOmega
  //
  void RealWignerBigD( int j, double Z, double Y, double Zp, double * D );


  //
  // This is equivalent to the function above, but the rotation matrix
  // is passed-in, instead of the euler angles Z,Y,Zp
  //
  // Y(R.r) = D(Z,Y,Zp) . Y(r)
  //
  void RealWignerBigD( int j, double const * R, double * D );

}


#endif
