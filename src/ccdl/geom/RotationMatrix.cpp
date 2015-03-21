#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>

#include "RotationMatrix.hpp"
#include "../math.hpp"

/*

void ZYZConsistentAngles( double & Z, double & Y, double & Zp, double * R )
{
  double c1 = std::cos(Z);
  double c2 = std::cos(Y);
  double c3 = std::cos(Zp);
  double s1 = std::sin(Z);
  double s2 = std::sin(Y);
  double s3 = std::sin(Zp);
  R[0+0*3] =  c1*c2*c3 - s1*s3;
  R[1+0*3] =  c1*s3    + c2*c3*s1;
  R[2+0*3] = -c3*s2;
  R[0+1*3] = -c3*s1    - c1*c2*s3;
  R[1+1*3] =  c1*c3    - c2*s1*s3;
  R[2+1*3] =  s2*s3;
  R[0+2*3] =  c1*s2;
  R[1+2*3] =  s1*s2;
  R[2+2*3] =  c2;


  if ( R[2+2*3] < 1. )
    {
      if ( R[2+2*3] > -1. )
	{
	  Z  = std::atan2( R[1+2*3],  R[0+2*3] );
	  Y  = std::acos(  R[2+2*3] );
	  Zp = std::atan2( R[2+1*3], -R[2+0*3] );
	}
      else
	{
	  // not a unique solution
	  Z = -std::atan2( R[1+0*3], R[1+1*3] );
	  Y = M_PI;
	  Zp = 0.;
	};
    }
  else
    {
      // not a unique solution
      Z = std::atan2( R[1+0*3], R[1+1*3] );
      Y = 0.;
      Zp = 0.;
    };

}

*/


//
// http://www.geometrictools.com/Documentation/EulerAngles.pdf
//
void ccdl::RotationMatrixFromZYZ( double Z, double Y, double Zp, double * R )
{
  //ZYZConsistentAngles( Z, Y, Zp, R );
  double c1 = std::cos(Z);
  double c2 = std::cos(Y);
  double c3 = std::cos(Zp);
  double s1 = std::sin(Z);
  double s2 = std::sin(Y);
  double s3 = std::sin(Zp);
  R[0+0*3] =  c1*c2*c3 - s1*s3;
  R[1+0*3] =  c1*s3    + c2*c3*s1;
  R[2+0*3] = -c3*s2;
  R[0+1*3] = -c3*s1    - c1*c2*s3;
  R[1+1*3] =  c1*c3    - c2*s1*s3;
  R[2+1*3] =  s2*s3;
  R[0+2*3] =  c1*s2;
  R[1+2*3] =  s1*s2;
  R[2+2*3] =  c2;
}


//
// http://www.geometrictools.com/Documentation/EulerAngles.pdf
//
void ccdl::ZYZFromRotationMatrix( double const * R, double & Z, double & Y, double & Zp )
{
  // 22 12 02 21 20 10 11

  if ( R[2+2*3] < 1. )
    {
      if ( R[2+2*3] > -1. )
	{
	  Z  = std::atan2( R[1+2*3],  R[0+2*3] );
	  Y  = std::acos(  R[2+2*3] );
	  Zp = std::atan2( R[2+1*3], -R[2+0*3] );
	}
      else
	{
	  // not a unique solution
	  Z = -std::atan2( R[1+0*3], R[1+1*3] );
	  Y = M_PI;
	  Zp = 0.;
	};
    }
  else
    {
      // not a unique solution
      Z = std::atan2( R[1+0*3], R[1+1*3] );
      Y = 0.;
      Zp = 0.;
    };
}



void ccdl::RotationMatrixGradientsFromZYZ
( double  Z, double  Y, double Zp, 
  double * dRdZ,
  double * dRdY,
  double * dRdZp )
{
  //ZYZConsistentAngles( Z, Y, Zp, dRdZ );

  double c1 = std::cos(Z);
  double c2 = std::cos(Y);
  double c3 = std::cos(Zp);
  double s1 = std::sin(Z);
  double s2 = std::sin(Y);
  double s3 = std::sin(Zp);

  double dc1 = -s1;
  double dc2 = -s2;
  double dc3 = -s3;
  double ds1 =  c1;
  double ds2 =  c2;
  double ds3 =  c3;

  dRdZ[0+0*3] =  dc1*c2*c3 - ds1*s3;
  dRdZ[1+0*3] =  dc1*s3    + c2*c3*ds1;
  dRdZ[2+0*3] = 0.;
  dRdZ[0+1*3] = -c3*ds1    - dc1*c2*s3;
  dRdZ[1+1*3] =  dc1*c3    - c2*ds1*s3;
  dRdZ[2+1*3] =  0.;
  dRdZ[0+2*3] =  dc1*s2;
  dRdZ[1+2*3] =  ds1*s2;
  dRdZ[2+2*3] =  0.;



  dRdY[0+0*3] =  c1*dc2*c3;
  dRdY[1+0*3] =  dc2*c3*s1;
  dRdY[2+0*3] = -c3*ds2;
  dRdY[0+1*3] = - c1*dc2*s3;
  dRdY[1+1*3] = - dc2*s1*s3;
  dRdY[2+1*3] =  ds2*s3;
  dRdY[0+2*3] =  c1*ds2;
  dRdY[1+2*3] =  s1*ds2;
  dRdY[2+2*3] =  dc2;

  dRdZp[0+0*3] =  c1*c2*dc3 - s1*ds3;
  dRdZp[1+0*3] =  c1*ds3    + c2*dc3*s1;
  dRdZp[2+0*3] = -dc3*s2;
  dRdZp[0+1*3] = -dc3*s1    - c1*c2*ds3;
  dRdZp[1+1*3] =  c1*dc3    - c2*s1*ds3;
  dRdZp[2+1*3] =  s2*ds3;
  dRdZp[0+2*3] =  0.;
  dRdZp[1+2*3] =  0.;
  dRdZp[2+2*3] =  0.;
}



//
// https://www.ics.forth.gr/_publications/2000_eccv_SVD_jacobian.pdf
//
// A = U.D.VT
// Aij = Uik dk Vjk
// Aij = Uik dk VTkj
//
// d(dk)/d(Aij) = Uik Vjk
// d(Ukl)/d(Aij) = Ukl' O^{ij}_{l'j}
// d(Vkl)/d(Aij) = -Vkl' P^{ij}_{l'j}
//
// O^{ij}_{l'j} and P^{ij}_{l'j} are obtained from
// solving the set 2x2 linear equations
//   H . x = b
//  x = [ O^{ij}_{l'j} , P^{ij}_{l'j} ]
//  b = [ Uik Vjl'     , - Uil' Vjk   ]
//  H = [ dl'  dk  ]
//      [ dk   dl' ]
//  x = H^{-1}.b
//  H^{-1} = ( dl' * dl' - dk * dk )^{-1} [  dl'  -dk  ]
//                                        [ -dk    dl' ]
// O^{ij}_{l'j} =  (  (Uik Vjl') dl' + (Uil' Vjk) * dk  ) / det
// P^{ij}_{l'j} = -(  (Uik Vjl') dk  + (Uil' Vjk) * dl' ) / det



namespace ccdl
{
  namespace RotationMatrix
  {
    void SvdGradients
    ( int nf, int ns, 
      double const * U,  // (nf*nf)
      double const * W,  // min(nf,ns)
      double const * VT, // (ns*ns)
      double * dUdA,   // min(nf,ns) * (nf*ns)
      double * dWdA,   // (nf*nf) * (nf*ns)
      double * dVTdA ); // (ns*ns) * (nf*ns)
    //
    // The gradients are the slow index (NOTE: this is unusual for me)
    //
  }
}


void ccdl::RotationMatrix::SvdGradients
( int nf, int ns, 
  double const * U,  // (nf*nf)
  double const * W,  // min(nf,ns)
  double const * VT, // (ns*ns)
  double * dUdA,   // min(nf,ns) * (nf*ns)
  double * dWdA,   // (nf*nf) * (nf*ns)
  double * dVTdA ) // (ns*ns) * (nf*ns)
{
  int n = nf*ns;
  int nmin = std::min(nf,ns);
  for ( int i=0; i < n*nmin; ++i )
    dWdA[i] = 0.;

  for ( int j=0; j<ns; ++j )
    for ( int i=0; i<nf; ++i )
      {
	for ( int k=0; k<nmin; ++k )
	  dWdA[ k + (i+j*nf)*nmin ] = U[i+k*nf] * VT[k+j*ns];
      };

  for ( int j=0; j<ns; ++j )
    for ( int i=0; i<nf; ++i )
      {
	for ( int l=0; l<ns; ++l )
	  for ( int k=0; k<ns; ++k )
	    {
	      double dl = (l < nmin) ? W[l] : 0.;
	      double dk = (k < nmin) ? W[k] : 0.;
	      double det = dl*dl - dk*dk;
	      if ( std::abs( det ) < 1.e-10 ) continue;

	      double u = U[i+k*nf] * VT[l+j*ns];
	      double v = U[i+l*nf] * VT[k+j*ns];
	      double Okl =   ( u * dl + v * dk ) / det;
	      double Pkl = - ( u * dk + v * dl ) / det;
	      for ( int kp=0; kp<nf; ++kp )
		dUdA[ (kp+l*nf) + (i+j*nf)*(nf*nf) ] += U[kp+k*nf]*Okl;
	      for ( int kp=0; kp<ns; ++kp )
		dVTdA[ (l+kp*ns) + (i+j*nf)*(ns*ns) ] -= VT[k+kp*ns]*Pkl;
	    };
      };

}



double ccdl::RmsTransformData
( int n, 
  double const * W,
  double const * rcrd,
  double const * crd,
  double * R,
  double * rcom,
  double * com )
{
  double totw = 0.;
  for ( int i=0; i<n; ++i ) 
    totw += W[i];
  std::vector<double> w(W,W+n);
  if ( totw > 1.e-15 )
    for ( int i=0; i<n; ++i ) 
      w[i] /= totw;
  
  com[0] = 0.;
  com[1] = 0.;
  com[2] = 0.;
  for ( int i=0; i<n; ++i )
    for ( int k=0; k<3; ++k )
      com[k] += w[i] * crd[k+i*3];

  rcom[0] = 0.;
  rcom[1] = 0.;
  rcom[2] = 0.;
  for ( int i=0; i<n; ++i )
    for ( int k=0; k<3; ++k )
      rcom[k] += w[i] * rcrd[k+i*3];

  double sumsq = 0.;
  std::vector<double> A(9,0.);
  for ( int iat=0; iat<n; ++iat )
    for ( int i=0; i<3; ++i )
      {
	double c1 =  crd[i+iat*3] -  com[i];
	double c2 = rcrd[i+iat*3] - rcom[i];
	sumsq += w[iat] * ( c1*c1+c2*c2 );
	for ( int j=0; j<3; ++j )
	  A[j+i*3] += w[iat] * (crd[j+iat*3]-com[j]) * c2;
      };

  std::vector<double> U(9,0.), S(3,0.), VT(9,0.);
  ccdl::ge Umat(3,3,U.data());
  ccdl::di Smat(3,3,S.data());
  ccdl::ge VTmat(3,3,VT.data());
  ccdl::ge Amat(3,3,A.data());

  Amat.svd( Umat, Smat, VTmat );


  double const detU = 
      U[0] * (U[4]*U[8] - U[7]*U[5])
    - U[3] * (U[1]*U[8] - U[7]*U[2])
    + U[6] * (U[1]*U[5] - U[4]*U[2]);

  // double const detV = 
  //     VT[0] * (VT[4]*VT[8] - VT[7]*VT[5])
  //   - VT[3] * (VT[1]*VT[8] - VT[7]*VT[2])
  //   + VT[6] * (VT[1]*VT[5] - VT[4]*VT[2]);
  double const detV = 
      VT[0*3+0] * (VT[1*3+1]*VT[2*3+2] - VT[1*3+2]*VT[2*3+1])
    - VT[0*3+1] * (VT[1*3+0]*VT[2*3+2] - VT[1*3+2]*VT[2*3+0])
    + VT[0*3+2] * (VT[1*3+0]*VT[2*3+1] - VT[1*3+1]*VT[2*3+0]);


  double dot = S[0]+S[1]+S[2];
  if ( detU*detV < 0. )
    {
      int jmin = 0;
      if ( S[1] < S[jmin] ) jmin=1;
      if ( S[2] < S[jmin] ) jmin=2;
      dot -= 2.*S[jmin];
      U[0+jmin*3] = -U[0+jmin*3];
      U[1+jmin*3] = -U[1+jmin*3];
      U[2+jmin*3] = -U[2+jmin*3];
    };

  double const rmsPD = std::sqrt( std::max( sumsq - 2.0 * dot , 0. ) );

  //ccdl::ge Rmat(3,3,R);
  //Rmat.dot( VT.t(), U.t() );
  for ( int i=0; i<9; ++i )
    R[i] = 0.;
  for ( int j=0; j<3; ++j )
    for ( int i=0; i<3; ++i )
      for ( int k=0; k<3; ++k )
	R[i+j*3] += VT[k+i*3] * U[j+k*3];

  return rmsPD;
}


void ccdl::RotateCrds
( int nat, 
  double * crd,
  double const * R,
  double const * rcom,
  double const * com )
{
  for ( int a=0; a<nat; ++a )
    {
      double crdp[3] = {0.,0.,0.};
      double c[3] = { crd[0+a*3]-com[0],
		      crd[1+a*3]-com[1],
		      crd[2+a*3]-com[2] };
		      
      for ( int j=0; j<3; ++j )
	for ( int k=0; k<3; ++k )
	  crdp[k] += R[k+j*3] * c[j]; // regular...

      crd[0+a*3] = crdp[0] + rcom[0];
      crd[1+a*3] = crdp[1] + rcom[1];
      crd[2+a*3] = crdp[2] + rcom[2];
    };
}

void ccdl::UnrotateCrds
( int nat, 
  double * crd,
  double const * R,
  double const * rcom,
  double const * com )
{
  for ( int a=0; a<nat; ++a )
    {
      double crdp[3] = {0.,0.,0.};
      double c[3] = { crd[0+a*3]-rcom[0],
		      crd[1+a*3]-rcom[1],
		      crd[2+a*3]-rcom[2] };
		      
      for ( int j=0; j<3; ++j )
	for ( int k=0; k<3; ++k )
	  crdp[k] += R[j+k*3] * c[j]; // transposed...

      crd[0+a*3] = crdp[0] + com[0];
      crd[1+a*3] = crdp[1] + com[1];
      crd[2+a*3] = crdp[2] + com[2];
    };
}


void ccdl::RmsOverlay
( int nat, 
  double const * w,
  double const * rcrd,
  double * crd )
{
  double R[9];
  double rcom[3];
  double com[3];
  ccdl::RmsTransformData( nat, w, rcrd, crd, R, rcom, com );
  ccdl::RotateCrds( nat, crd, R, rcom, com );
}
		     
  
void ccdl::RmsOverlayGrd
( int n, 
  double const * W,
  double const * rcrd,
  double const * crd,
  double const * rotated_grd,
  double * unrotated_grd )
{

  double totw = 0.;
  for ( int i=0; i<n; ++i ) 
    totw += W[i];
  std::vector<double> w(W,W+n);
  if ( totw > 1.e-15 )
    for ( int i=0; i<n; ++i ) 
      w[i] /= totw;
  
  double com[3] = {0.,0.,0.};
  for ( int i=0; i<n; ++i )
    for ( int k=0; k<3; ++k )
      com[k] += w[i] * crd[k+i*3];

  double rcom[3] = {0.,0.,0.};
  for ( int i=0; i<n; ++i )
    for ( int k=0; k<3; ++k )
      rcom[k] += w[i] * rcrd[k+i*3];

  double sumsq = 0.;
  std::vector<double> A(9,0.);
  for ( int iat=0; iat<n; ++iat )
    for ( int i=0; i<3; ++i )
      {
	double c1 =  crd[i+iat*3] -  com[i];
	double c2 = rcrd[i+iat*3] - rcom[i];
	sumsq += w[iat] * ( c1*c1+c2*c2 );
	for ( int j=0; j<3; ++j )
	  A[j+i*3] += w[iat] * (crd[j+iat*3]-com[j]) * c2;
      };

  std::vector<double> U(9,0.), S(3,0.), VT(9,0.);
  ccdl::ge Umat(3,3,U.data());
  ccdl::di Smat(3,3,S.data());
  ccdl::ge VTmat(3,3,VT.data());
  ccdl::ge Amat(3,3,A.data());

  Amat.svd( Umat, Smat, VTmat );

  std::vector<double> dU(9*9,0.), dS(3*9,0.), dVT(9*9,0.);
  ccdl::RotationMatrix::SvdGradients( 3,3, U.data(), S.data(), VT.data(),
				      dU.data(), dS.data(), dVT.data() );


  double const detU = 
      U[0] * (U[4]*U[8] - U[7]*U[5])
    - U[3] * (U[1]*U[8] - U[7]*U[2])
    + U[6] * (U[1]*U[5] - U[4]*U[2]);

  // double const detV = 
  //     VT[0] * (VT[4]*VT[8] - VT[7]*VT[5])
  //   - VT[3] * (VT[1]*VT[8] - VT[7]*VT[2])
  //   + VT[6] * (VT[1]*VT[5] - VT[4]*VT[2]);
  double const detV = 
      VT[0*3+0] * (VT[1*3+1]*VT[2*3+2] - VT[1*3+2]*VT[2*3+1])
    - VT[0*3+1] * (VT[1*3+0]*VT[2*3+2] - VT[1*3+2]*VT[2*3+0])
    + VT[0*3+2] * (VT[1*3+0]*VT[2*3+1] - VT[1*3+1]*VT[2*3+0]);


  double dot = S[0]+S[1]+S[2];
  if ( detU*detV < 0. )
    {
      int jmin = 0;
      if ( S[1] < S[jmin] ) jmin=1;
      if ( S[2] < S[jmin] ) jmin=2;
      dot -= 2.*S[jmin];
      for ( int k=0; k<3; ++k )
	U[k+jmin*3] *= -1.;
      for ( int ij=0; ij<9; ++ij )
	for ( int k=0; k<3; ++k )
	  dU[ (k+jmin*3) + ij*9 ] *= -1.;
    };


  //ccdl::ge Rmat(3,3,R);
  //Rmat.dot( VT.t(), U.t() );
  double R[9];
  for ( int i=0; i<9; ++i )
    R[i] = 0.;

  double dRdA[9*9];
  for ( int i=0; i<9*9; ++i )
    dRdA[i] = 0.;

  for ( int j=0; j<3; ++j )
    for ( int i=0; i<3; ++i )
      for ( int k=0; k<3; ++k )
	{
	  R[i+j*3] += VT[k+i*3] * U[j+k*3];
	  for ( int aij=0; aij<9; ++aij )
	    dRdA[(i+j*3)+aij*9] += 
	      dVT[(k+i*3)+aij*9] *  U[j+k*3] +
	      VT[k+i*3]          * dU[(j+k*3)+aij*9];
	};



  for ( int i=0; i<3*n; ++i )
    unrotated_grd[i] = 0.;


  std::vector<double> dRdXa(9,0.);
  std::vector<double> dXdXa(3*n,0.);
  for ( int a=0; a<n; ++a )
    for ( int k=0; k<3; ++k )
      {
	dRdXa.assign( dRdXa.size(), 0. );
	dXdXa.assign( dXdXa.size(), 0. );
	if ( w[a] > 1.e-13 )
	  {
	    for ( int l=0; l<3; ++l )
	      {
		double f = w[a] * (rcrd[l+a*3]-rcom[l]);
		for ( int j=0; j<3; ++j )
		  for ( int i=0; i<3; ++i )
		    dRdXa[ i+j*3 ] += dRdA[(i+j*3)+(k+l*3)*9] * f;
	      };
	    
	    for ( int b=0; b<n; ++b )
	      for ( int kp=0; kp<3; ++kp )
		{
		  for ( int kpp=0; kpp<3; ++kpp )
		    dXdXa[kpp+b*3] += dRdXa[kpp+kp*3] * (crd[kp+b*3]-com[kp]);
		  dXdXa[kp+b*3] -= R[kp+k*3] * w[a];
		};
	  }; // w[a] is non-zero

	for ( int kp=0; kp<3; ++kp ) 
	  dXdXa[kp+a*3] += R[kp+k*3];

	double dEdXa = 0.;
	for ( int b=0; b<n; ++b )
	  for ( int kp=0; kp<3; ++kp )
	    dEdXa += rotated_grd[kp+b*3] * dXdXa[kp+b*3];
	unrotated_grd[k+a*3] += dEdXa;
      };

  /*
  std::vector<double> dRdX(3*n*9,0.);
  for ( int a=0; a<n; ++a )
    {
      if ( std::abs(w[a]) < 1.e-12 ) continue;
      for ( int k=0; k<3; ++k )
	{
	  for ( int j=0; j<3; ++j )
	    for ( int i=0; i<3; ++i )
	      for ( int l=0; l<3; ++l )
		dRdX[ (k+a*3) + (i+j*3)*(3*n) ] += 
		  dRdA[(i+j*3)+(k+l*3)*9] * ( w[a] * ( rcrd[l+a*3]-rcom[l] ) );
	};
    };

 
  for ( int i=0; i<3*n; ++i )
    unrotated_grd[i] = 0.;

  for ( int a=0; a<n; ++a )
    for ( int k=0; k<3; ++k )
      {
	std::vector<double> dXdXa(3*n,0.);
	for ( int b=0; b<n; ++b )
	  for ( int kp=0; kp<3; ++kp )
	    {
	      for ( int kpp=0; kpp<3; ++kpp )
		dXdXa[kpp+b*3] += 
		  dRdX[(k+a*3)+(kpp+kp*3)*3*n] * (crd[kp+b*3]-com[kp]);
	      dXdXa[kp+b*3] -= R[kp+k*3] * (w[a]);
	    };
	for ( int kp=0; kp<3; ++kp ) 
	  dXdXa[kp+a*3] += R[kp+k*3];

	double dEdXa = 0.;
	for ( int b=0; b<n; ++b )
	  for ( int kp=0; kp<3; ++kp )
	    dEdXa += rotated_grd[kp+b*3] * dXdXa[kp+b*3];
	unrotated_grd[k+a*3] += dEdXa;
      };
  */

  /* 
  std::vector<double> dXdX(3*n*3*n,0.);
  ccdl::Rms_dXdX( n, W, rcrd, crd, dXdX.data() );
 
  for ( int i=0; i<3*n; ++i )
    unrotated_grd[i] = 0.;

  for ( int a=0; a<n; ++a )
    for ( int k=0; k<3; ++k )
      {
	double dE = 0.;
	for ( int b=0; b<n; ++b )
	  for ( int kp=0; kp<3; ++kp )
	    dE += rotated_grd[kp+b*3] * dXdX[ (k+a*3) + (kp+b*3)*(3*n) ];
	unrotated_grd[k+a*3] += dE;
      };
  */
}












///////////////////////////////////////
///////////////////////////////////////
// TESTING JUNK
///////////////////////////////////////
///////////////////////////////////////



void ccdl::RotationMatrix::test_RmsOverlayGrd
( int n, 
  double const * W,
  double const * rcrd,
  double * crd )
{

  std::vector<double> ocrd(crd,crd+3*n);
  std::vector<double> pgrd(3*n,0.);
  std::vector<double> ogrd(3*n,0.);

  ccdl::RmsOverlay( n, W, rcrd, crd );

  // e = \sum_a Xa*Xa + Ya*Ya + Za*Za
  for ( int a=0; a<n; ++a )
    for ( int k=0; k<3; ++k )
      pgrd[k+a*3] = 2. * crd[k+a*3];

  //for ( int k=0; k<3; ++k ) pgrd[k+0*3] = 2. * crd[k+0*3];


  ccdl::RmsOverlayGrd( n, W, rcrd, ocrd.data(), pgrd.data(), ogrd.data() );

  double DEL = 1.e-6;

  for ( int a=0; a<n; ++a )
    for ( int k=0; k<3; ++k )
      {
	for ( int j=0; j<3*n; ++j )
	  crd[j] = ocrd[j];
	crd[k+a*3] += DEL;

	// for ( int b=0; b<n; ++b )
	//   printf("A %12.7f %12.7f %12.7f    %12.7f %12.7f %12.7f\n",
	// 	 crd[0+b*3],crd[1+b*3],crd[2+b*3],
	// 	 ocrd[0+b*3],ocrd[1+b*3],ocrd[2+b*3] );

	ccdl::RmsOverlay( n, W, rcrd, crd );

	// for ( int b=0; b<n; ++b )
	//   printf("B %12.7f %12.7f %12.7f    %12.7f %12.7f %12.7f\n",
	// 	 crd[0+b*3],crd[1+b*3],crd[2+b*3],
	// 	 ocrd[0+b*3],ocrd[1+b*3],ocrd[2+b*3] );

	double Ehi = 0.;
	for ( int j=0; j<3*n; ++j ) Ehi += crd[j]*crd[j];
	//Ehi = crd[0]*crd[0] + crd[1]*crd[1] + crd[2]*crd[2];


	for ( int j=0; j<3*n; ++j )
	  crd[j] = ocrd[j];
	crd[k+a*3] -= DEL;
	ccdl::RmsOverlay( n, W, rcrd, crd );

	// for ( int b=0; b<n; ++b )
	//   printf("C %12.7f %12.7f %12.7f\n",
	// 	 crd[0+b*3],crd[1+b*3],crd[2+b*3]);

	double Elo = 0.;
	for ( int j=0; j<3*n; ++j ) Elo += crd[j]*crd[j];
	//Elo = crd[0]*crd[0] + crd[1]*crd[1] + crd[2]*crd[2];



	double num = (Ehi-Elo) / (2.*DEL);
	double ana = ogrd[k+a*3];
	double d = num-ana;

	printf("%i %i  %12.4e %12.4e %12.4e  %12.4e %12.4e\n",
	       a,k, num, ana, d, Ehi, Elo );
      };
}




// namespace ccdl
// {
//   namespace RotationMatrix
//   {

    
//     void test_SvdGradients( int nf, int ns, double * A );
    
//     void test_dRdX
//     ( int n, double const * W,
//       double const * rcrd, double const * crd );
    
    
//     void test_dAdX
//     ( int n, double const * W,
//       double const * rcrd, double const * crd );
    
//     void test_dXdX
//     ( int n, double const * W,
//       double const * rcrd, double const * crd );
    
    
//     void Rms_dAdX
//     ( int n,  double const * W,
//       double const * rcrd, double const * crd,
//       double * dAdX );
    
//     void Rms_dRdX
//     ( int n, double const * W,
//       double const * rcrd, double const * crd,
//       double * dRdX );
    
//     void Rms_A
//     ( int n, double const * W,
//       double const * rcrd, double const * crd,
//       double * A );
    
//     void Rms_dXdX
//     ( int n, double const * W,
//       double const * rcrd, double const * crd,
//       double * dXdX );

//   }
// }










// void ccdl::RotationMatrix::Rms_A
// ( int n, 
//   double const * W,
//   double const * rcrd,
//   double const * crd,
//   double * A )
// {
//   double totw = 0.;
//   for ( int i=0; i<n; ++i ) 
//     totw += W[i];
//   std::vector<double> w(W,W+n);
//   if ( totw > 1.e-15 )
//     for ( int i=0; i<n; ++i ) 
//       w[i] /= totw;
  
//   double com[3];
//   com[0] = 0.;
//   com[1] = 0.;
//   com[2] = 0.;
//   for ( int i=0; i<n; ++i )
//     for ( int k=0; k<3; ++k )
//       com[k] += w[i] * crd[k+i*3];

//   double rcom[3];
//   rcom[0] = 0.;
//   rcom[1] = 0.;
//   rcom[2] = 0.;
//   for ( int i=0; i<n; ++i )
//     for ( int k=0; k<3; ++k )
//       rcom[k] += w[i] * rcrd[k+i*3];

//   for ( int i=0; i<9; ++i )
//     A[i] = 0.;
//   for ( int iat=0; iat<n; ++iat )
//     for ( int i=0; i<3; ++i )
//       {
// 	double c2 = rcrd[i+iat*3] - rcom[i];
// 	for ( int j=0; j<3; ++j )
// 	  A[j+i*3] += w[iat] * (crd[j+iat*3]-com[j]) * c2;
//       };
// }


// void ccdl::RotationMatrix::Rms_dRdX
// ( int n, 
//   double const * W,
//   double const * rcrd,
//   double const * crd,
//   double * dRdX )
// {

//   double totw = 0.;
//   for ( int i=0; i<n; ++i ) 
//     totw += W[i];
//   std::vector<double> w(W,W+n);
//   if ( totw > 1.e-15 )
//     for ( int i=0; i<n; ++i ) 
//       w[i] /= totw;
  
//   double com[3] = {0.,0.,0.};
//   for ( int i=0; i<n; ++i )
//     for ( int k=0; k<3; ++k )
//       com[k] += w[i] * crd[k+i*3];

//   double rcom[3] = {0.,0.,0.};
//   for ( int i=0; i<n; ++i )
//     for ( int k=0; k<3; ++k )
//       rcom[k] += w[i] * rcrd[k+i*3];

//   double sumsq = 0.;
//   std::vector<double> A(9,0.);
//   for ( int iat=0; iat<n; ++iat )
//     for ( int i=0; i<3; ++i )
//       {
// 	double c1 =  crd[i+iat*3] -  com[i];
// 	double c2 = rcrd[i+iat*3] - rcom[i];
// 	sumsq += w[iat] * ( c1*c1+c2*c2 );
// 	for ( int j=0; j<3; ++j )
// 	  A[j+i*3] += w[iat] * (crd[j+iat*3]-com[j]) * c2;
//       };

//   std::vector<double> U(9,0.), S(3,0.), VT(9,0.);
//   ccdl::ge Umat(3,3,U.data());
//   ccdl::di Smat(3,3,S.data());
//   ccdl::ge VTmat(3,3,VT.data());
//   ccdl::ge Amat(3,3,A.data());

//   Amat.svd( Umat, Smat, VTmat );

//   std::vector<double> dU(9*9,0.), dS(3*9,0.), dVT(9*9,0.);
//   ccdl::RotationMatrix::SvdGradients( 3,3, U.data(), S.data(), VT.data(),
// 				      dU.data(), dS.data(), dVT.data() );


//   double const detU = 
//       U[0] * (U[4]*U[8] - U[7]*U[5])
//     - U[3] * (U[1]*U[8] - U[7]*U[2])
//     + U[6] * (U[1]*U[5] - U[4]*U[2]);

//   double const detV = 
//       VT[0*3+0] * (VT[1*3+1]*VT[2*3+2] - VT[1*3+2]*VT[2*3+1])
//     - VT[0*3+1] * (VT[1*3+0]*VT[2*3+2] - VT[1*3+2]*VT[2*3+0])
//     + VT[0*3+2] * (VT[1*3+0]*VT[2*3+1] - VT[1*3+1]*VT[2*3+0]);


//   double dot = S[0]+S[1]+S[2];
//   if ( detU*detV < 0. )
//     {
//       int jmin = 0;
//       if ( S[1] < S[jmin] ) jmin=1;
//       if ( S[2] < S[jmin] ) jmin=2;
//       dot -= 2.*S[jmin];
//       for ( int k=0; k<3; ++k )
// 	U[k+jmin*3] *= -1.;
//       for ( int ij=0; ij<9; ++ij )
// 	for ( int k=0; k<3; ++k )
// 	  dU[ (k+jmin*3) + ij*9 ] *= -1.;
//     };

//   double R[9];
//   for ( int i=0; i<9; ++i )
//     R[i] = 0.;

//   double dRdA[9*9];
//   for ( int i=0; i<9*9; ++i )
//     dRdA[i] = 0.;


//   for ( int j=0; j<3; ++j )
//     for ( int i=0; i<3; ++i )
//       for ( int k=0; k<3; ++k )
// 	{
// 	  R[i+j*3] += VT[k+i*3] * U[j+k*3];
// 	  for ( int aij=0; aij<9; ++aij )
// 	    dRdA[(i+j*3)+aij*9] += 
// 	      dVT[(k+i*3)+aij*9] *  U[j+k*3] +
// 	      VT[k+i*3]          * dU[(j+k*3)+aij*9];
// 	};



//   for ( int i=0; i<9*3*n; ++i )
//     dRdX[i] = 0.;

//   for ( int j=0; j<3; ++j )
//     for ( int i=0; i<3; ++i )
//       for ( int a=0; a<n; ++a )
// 	for ( int k=0; k<3; ++k )
// 	  for ( int l=0; l<3; ++l )
// 	    dRdX[ (k+a*3) + (i+j*3)*(3*n) ] += 
// 	      dRdA[(i+j*3)+(k+l*3)*9] * ( w[a] * ( rcrd[l+a*3]-rcom[l] ) );
// }
		



// void ccdl::RotationMatrix::Rms_dXdX
// ( int n, 
//   double const * W,
//   double const * rcrd,
//   double const * crd,
//   double * dXdX )
// {

//   double totw = 0.;
//   for ( int i=0; i<n; ++i ) 
//     totw += W[i];
//   std::vector<double> w(W,W+n);
//   if ( totw > 1.e-15 )
//     for ( int i=0; i<n; ++i ) 
//       w[i] /= totw;
  
//   double com[3] = {0.,0.,0.};
//   for ( int i=0; i<n; ++i )
//     for ( int k=0; k<3; ++k )
//       com[k] += w[i] * crd[k+i*3];

//   double rcom[3] = {0.,0.,0.};
//   for ( int i=0; i<n; ++i )
//     for ( int k=0; k<3; ++k )
//       rcom[k] += w[i] * rcrd[k+i*3];

//   double sumsq = 0.;
//   std::vector<double> A(9,0.);
//   for ( int iat=0; iat<n; ++iat )
//     for ( int i=0; i<3; ++i )
//       {
// 	double c1 =  crd[i+iat*3] -  com[i];
// 	double c2 = rcrd[i+iat*3] - rcom[i];
// 	sumsq += w[iat] * ( c1*c1+c2*c2 );
// 	for ( int j=0; j<3; ++j )
// 	  A[j+i*3] += w[iat] * (crd[j+iat*3]-com[j]) * c2;
//       };

//   std::vector<double> U(9,0.), S(3,0.), VT(9,0.);
//   ccdl::ge Umat(3,3,U.data());
//   ccdl::di Smat(3,3,S.data());
//   ccdl::ge VTmat(3,3,VT.data());
//   ccdl::ge Amat(3,3,A.data());

//   Amat.svd( Umat, Smat, VTmat );

//   std::vector<double> dU(9*9,0.), dS(3*9,0.), dVT(9*9,0.);
//   ccdl::RotationMatrix::SvdGradients( 3,3, U.data(), S.data(), VT.data(),
// 				      dU.data(), dS.data(), dVT.data() );


//   double const detU = 
//       U[0] * (U[4]*U[8] - U[7]*U[5])
//     - U[3] * (U[1]*U[8] - U[7]*U[2])
//     + U[6] * (U[1]*U[5] - U[4]*U[2]);

//   double const detV = 
//       VT[0*3+0] * (VT[1*3+1]*VT[2*3+2] - VT[1*3+2]*VT[2*3+1])
//     - VT[0*3+1] * (VT[1*3+0]*VT[2*3+2] - VT[1*3+2]*VT[2*3+0])
//     + VT[0*3+2] * (VT[1*3+0]*VT[2*3+1] - VT[1*3+1]*VT[2*3+0]);


//   double dot = S[0]+S[1]+S[2];
//   if ( detU*detV < 0. )
//     {
//       int jmin = 0;
//       if ( S[1] < S[jmin] ) jmin=1;
//       if ( S[2] < S[jmin] ) jmin=2;
//       dot -= 2.*S[jmin];
//       for ( int k=0; k<3; ++k )
// 	U[k+jmin*3] *= -1.;
//       for ( int ij=0; ij<9; ++ij )
// 	for ( int k=0; k<3; ++k )
// 	  dU[ (k+jmin*3) + ij*9 ] *= -1.;
//     };

//   double R[9];
//   for ( int i=0; i<9; ++i )
//     R[i] = 0.;


//   for ( int j=0; j<3; ++j )
//     for ( int i=0; i<3; ++i )
//       for ( int k=0; k<3; ++k )
// 	{
// 	  R[i+j*3] += VT[k+i*3] * U[j+k*3];
// 	};

//   std::vector<double> dRdX(9*3*n,0.);
//   ccdl::RotationMatrix::Rms_dRdX( n, W, rcrd, crd, dRdX.data() );
//   for ( int i=0; i<3*n*3*n; ++i )
//     dXdX[i] = 0.;


//   for ( int a=0; a<n; ++a )
//     for ( int k=0; k<3; ++k )
//       {
// 	for ( int b=0; b<n; ++b )
// 	  for ( int kp=0; kp<3; ++kp )
// 	    {
// 	      for ( int kpp=0; kpp<3; ++kpp )
// 		dXdX[(k+a*3) + (kpp+b*3)*3*n] += 
// 		  dRdX[(k+a*3)+(kpp+kp*3)*3*n] * (crd[kp+b*3]-com[kp]);
// 	      dXdX[(k+a*3) + (kp+b*3)*3*n] -= R[kp+k*3] * (w[a]);
// 	    };
// 	for ( int kp=0; kp<3; ++kp ) dXdX[(k+a*3) + (kp+a*3)*3*n] += R[kp+k*3] ;
//       };

// }









// void ccdl::RotationMatrix::test_SvdGradients( int nf, int ns, double * A )
// {
//   int nmin = std::min(nf,ns);
//   std::vector<double> U( nf*nf, 0. );
//   std::vector<double> W( nmin, 0. );
//   std::vector<double> VT( ns*ns, 0. );
//   std::vector<double> dUdA(  nf*nf * nf*ns, 0. );
//   std::vector<double> dWdA(  nmin  * nf*ns, 0. );
//   std::vector<double> dVTdA( ns*ns * nf*ns, 0. );

//   {
//     ccdl::ge Umat( nf,nf,U.data() );
//     ccdl::di Wmat( nf,ns,W.data() );
//     ccdl::ge VTmat( ns,ns,VT.data() );
//     ccdl::ge( nf,ns,A ).svd( Umat, Wmat, VTmat );
//   };

//   ccdl::RotationMatrix::SvdGradients( nf,ns,
// 				      U.data(), W.data(), VT.data(),
// 				      dUdA.data(),
// 				      dWdA.data(),
// 				      dVTdA.data() );
  
//   double DEL = 1.e-6;
//   for ( int j=0; j<ns; ++j )
//     for ( int i=0; i<nf; ++i )
//       {
// 	std::vector<double> Uhi( nf*nf, 0. );
// 	std::vector<double> Whi( nmin, 0. );
// 	std::vector<double> VThi( ns*ns, 0. );
// 	std::vector<double> Ulo( nf*nf, 0. );
// 	std::vector<double> Wlo( nmin, 0. );
// 	std::vector<double> VTlo( ns*ns, 0. );

// 	A[i+j*nf] += DEL;
// 	{
// 	  ccdl::ge Umat( nf,nf,Uhi.data() );
// 	  ccdl::di Wmat( nf,ns,Whi.data() );
// 	  ccdl::ge VTmat( ns,ns,VThi.data() );
// 	  ccdl::ge( nf,ns,A ).svd( Umat, Wmat, VTmat );
// 	};


// 	A[i+j*nf] -= 2.*DEL;
// 	{
// 	  ccdl::ge Umat( nf,nf,Ulo.data() );
// 	  ccdl::di Wmat( nf,ns,Wlo.data() );
// 	  ccdl::ge VTmat( ns,ns,VTlo.data() );
// 	  ccdl::ge( nf,ns,A ).svd( Umat, Wmat, VTmat );
// 	};

// 	A[i+j*nf] += DEL;

// 	for ( int k=0; k<nmin; ++k )
// 	  {
// 	    double d = (Whi[k]-Wlo[k])/(2.*DEL);
// 	    double a = dWdA[ k + (i+j*nf)*nmin ];
// 	    printf("dW[%i]/dA[%i,%i] : %12.4e %12.4e %12.4e\n",
// 		   k,i,j,d,a,d-a);
// 	  }
// 	printf("\n");
// 	for ( int l=0; l<nf; ++l )
// 	  for ( int k=0; k<nf; ++k )
// 	    {
// 	      double d = (Uhi[k+l*nf]-Ulo[k+l*nf])/(2.*DEL);
// 	      double a = dUdA[ (k+l*nf) + (i+j*nf)*nf*nf ];
// 	      printf("dU[%i,%i]/dA[%i,%i] : %12.4e %12.4e %12.4e\n",
// 		     k,l,i,j,d,a,d-a);
// 	  }
// 	printf("\n");
// 	for ( int l=0; l<ns; ++l )
// 	  for ( int k=0; k<ns; ++k )
// 	    {
// 	      double d = (VThi[k+l*ns]-VTlo[k+l*ns])/(2.*DEL);
// 	      double a = dVTdA[ (k+l*ns) + (i+j*nf)*ns*ns ];
// 	      printf("dV[%i,%i]/dA[%i,%i] : %12.4e %12.4e %12.4e\n",
// 		     k,l,i,j,d,a,d-a);
// 	  }

// 	printf("\n\n");

//       };

// }



			

// void ccdl::RotationMatrix::test_dRdX
// ( int n, 
//   double const * W,
//   double const * rcrd,
//   double const * crd )
// {
//   std::vector<double> dRdX( 3*n*9, 0. );
//   ccdl::RotationMatrix::Rms_dRdX( n, W, rcrd, crd, dRdX.data() );

//   double DEL = 1.e-7;

//   for ( int a=0; a<n; ++a )
//     for ( int k=0; k<3; ++k )
//       {
//       std::vector<double> chi(crd,crd+3*n);
//       chi[k+a*3] += DEL;
//       double Rhi[9];
//       double rcomhi[3];
//       double comhi[3];
//       ccdl::RmsTransformData( n, W, rcrd, chi.data(), Rhi, rcomhi, comhi );

//       std::vector<double> clo(crd,crd+3*n);
//       clo[k+a*3] -= DEL;
//       double Rlo[9];
//       double rcomlo[3];
//       double comlo[3];
//       ccdl::RmsTransformData( n, W, rcrd, clo.data(), Rlo, rcomlo, comlo );


//       for ( int i=0; i<3; ++i )
// 	for ( int j=0; j<3; ++j )
// 	  {
// 	    double num = (Rhi[i+j*3]-Rlo[i+j*3])/(2.*DEL);
// 	    double ana = dRdX[ (k+a*3) + (i+j*3)*(3*n) ];
// 	    double d = num-ana;
// 	    if ( std::abs(d) > 5.e-9 )
// 	      printf("%i %i  %i %i %13.5e %13.5e %13.5e\n",
// 		     a,k,i,j,num,ana,d);
// 	  };

//     };

// }



// void ccdl::RotationMatrix::Rms_dAdX
// ( int n, 
//   double const * W,
//   double const * rcrd,
//   double const * crd,
//   double * dAdX )
// {
//   double totw = 0.;
//   for ( int i=0; i<n; ++i ) 
//     totw += W[i];
//   std::vector<double> w(W,W+n);
//   if ( totw > 1.e-15 )
//     for ( int i=0; i<n; ++i ) 
//       w[i] /= totw;
  
//   double com[3];
//   com[0] = 0.;
//   com[1] = 0.;
//   com[2] = 0.;
//   for ( int i=0; i<n; ++i )
//     for ( int k=0; k<3; ++k )
//       com[k] += w[i] * crd[k+i*3];

//   double rcom[3];
//   rcom[0] = 0.;
//   rcom[1] = 0.;
//   rcom[2] = 0.;
//   for ( int i=0; i<n; ++i )
//     for ( int k=0; k<3; ++k )
//       rcom[k] += w[i] * rcrd[k+i*3];

//   for ( int i=0; i<9*3*n; ++i )
//     dAdX[i] = 0.;
//   for ( int iat=0; iat<n; ++iat )
//     for ( int i=0; i<3; ++i )
//       {
// 	double c2 = rcrd[i+iat*3] - rcom[i];
// 	for ( int j=0; j<3; ++j )
// 	  dAdX[ (j+iat*3) + (j+i*3)*(3*n) ] +=  w[iat] * c2;
//       };
// }



// void ccdl::RotationMatrix::test_dAdX
// ( int n, 
//   double const * W,
//   double const * rcrd,
//   double const * crd )
// {
//   double DEL = 1.e-5;
//   std::vector<double> dAdX( 3*n*9, 0. );
//   ccdl::RotationMatrix::Rms_dAdX( n, W, rcrd, crd, dAdX.data() );
//   for ( int a=0; a<n; ++a )
//     for ( int k=0; k<3; ++k )
//       {
// 	std::vector<double> chi(crd,crd+3*n);
// 	chi[k+a*3] += DEL;
// 	std::vector<double> Ahi(9,0.);
// 	ccdl::RotationMatrix::Rms_A(n,W,rcrd,chi.data(),Ahi.data());

// 	std::vector<double> clo(crd,crd+3*n);
// 	clo[k+a*3] -= DEL;
// 	std::vector<double> Alo(9,0.);
// 	ccdl::RotationMatrix::Rms_A(n,W,rcrd,clo.data(),Alo.data());

// 	for ( int i=0; i<3; ++i )
// 	  for ( int j=0; j<3; ++j )
// 	    {
// 	      double ana = dAdX[ (k+a*3) + (i+j*3)*(3*n) ];
// 	      double num = (Ahi[i+j*3]-Alo[i+j*3])/(2.*DEL);
// 	      double d = num-ana;
// 	      printf("%i %i %i %i %12.4e %12.4e %12.4e\n",
// 		     a,k,i,j,num,ana,d);
// 	    };


//       };
// }





// void ccdl::RotationMatrix::test_dXdX
// ( int n, 
//   double const * W,
//   double const * rcrd,
//   double const * crd )
// {
//   std::vector<double> dXdX(3*n*3*n,0.);
//   ccdl::RotationMatrix::Rms_dXdX( n, W, rcrd, crd, dXdX.data() );

//   double DEL = 2.e-6;
//   for ( int a=0; a<n; ++a )
//     for ( int k=0; k<3; ++k )
//       {
// 	std::vector<double> chi(crd,crd+3*n);
// 	chi[k+a*3] += DEL;
// 	ccdl::RmsOverlay( n, W, rcrd, chi.data() );

// 	std::vector<double> clo(crd,crd+3*n);
// 	clo[k+a*3] -= DEL;
// 	ccdl::RmsOverlay( n, W, rcrd, clo.data() );

// 	for ( int b=0; b<n; ++b )
// 	  {
// 	    for ( int kp=0; kp<3; ++kp )
// 	      printf("%12.4e",rcrd[kp+3*b]);
// 	    printf("  ");
// 	    for ( int kp=0; kp<3; ++kp )
// 	      printf("%12.4e",chi[kp+3*b]);
// 	    printf("  ");
// 	    for ( int kp=0; kp<3; ++kp )
// 	      printf("%12.4e",clo[kp+3*b]);
// 	    printf("\n");
// 	  };


// 	for ( int b=0; b<n; ++b )
// 	  for ( int kp=0; kp<3; ++kp )
// 	    {
// 	      double ana = dXdX[ (k+a*3) + (kp+b*3)*(3*n) ];
// 	      double num = ( chi[kp+b*3]-clo[kp+b*3] ) / (2.*DEL);
// 	      double d = num-ana;
// 	      if ( std::abs(d) < 1.e-8 ) d = 0.;
// 	      printf("%i %i  %i %i  %12.4e %12.4e %12.4e\n",
// 		     a,k,b,kp,num,ana,d);
// 	    };
//       }
// }






