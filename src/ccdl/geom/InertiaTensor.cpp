#include "InertiaTensor.hpp"
#include "EuclideanGeometry.hpp"
#include "RotationMatrix.hpp"
#include "../constants.hpp"
#include <algorithm>
// # include <vector>


extern "C"
{
  void dsyev_(const char* jobz, const char* uplo,
	      const int* n, double* a, const int* lda,
	      double* w, double* work, const int* lwork, int* info);
}

// -----------------------------------------------------------------
// 3x3 moment of inertia matrix
// -----------------------------------------------------------------
void ccdl::InertiaTensor( int const nat, double const * m, double const * c, double * I )
{
  std::fill( I, I+9, 0. );
  double Rcom[3] = {0.,0.,0.};
  double mtot = 0.;
  for ( int a=0; a<nat; ++a )
    {
      mtot += m[a];
      Rcom[0] += m[a] * c[0+a*3];
      Rcom[1] += m[a] * c[1+a*3];
      Rcom[2] += m[a] * c[2+a*3];
    };
  Rcom[0] /= mtot;
  Rcom[1] /= mtot;
  Rcom[2] /= mtot;
  for ( int a=0; a<nat; ++a )
    {
      double x = c[0+a*3]-Rcom[0];
      double y = c[1+a*3]-Rcom[1];
      double z = c[2+a*3]-Rcom[2];

      I[0] += m[a]*( y*y + z*z );
      I[1] -= m[a]*x*y;
      I[2] -= m[a]*x*z;
      I[4] += m[a]*( x*x + z*z );
      I[5] -= m[a]*y*z;
      I[8] += m[a]*( x*x + y*y );
    };
  I[3] = I[1];
  I[6] = I[2];
  I[7] = I[5];
}


// -----------------------------------------------------------------
// inertia tensor eigenvalues
// -----------------------------------------------------------------
void ccdl::InertiaMoments( int const nat, double const * m, double const * c, double * E )
{
  double U[9];
  ccdl::InertiaTensor(nat,m,c,U);
  int N=3, LWORK=16, INFO=0;
  double WORK[16];
  dsyev_("N","U",&N,U,&N,E,WORK,&LWORK,&INFO);
}

// -----------------------------------------------------------------
// inertia tensor eigenvalues and gradients wrt c
// -----------------------------------------------------------------
void ccdl::InertiaMoments( int const nat, double const * m, double const * c, double * E, double * G )
{
  double Rcom[3] = {0.,0.,0.};
  double mtot = 0.;
  for ( int a=0; a<nat; ++a )
    {
      mtot += m[a];
      Rcom[0] += m[a] * c[0+a*3];
      Rcom[1] += m[a] * c[1+a*3];
      Rcom[2] += m[a] * c[2+a*3];
    };
  Rcom[0] /= mtot;
  Rcom[1] /= mtot;
  Rcom[2] /= mtot;

  double U[9];
  ccdl::InertiaTensor(nat,m,c,U);
  int N=3, LWORK=16, INFO=0;
  double WORK[16];
  dsyev_("V","U",&N,U,&N,E,WORK,&LWORK,&INFO);
  std::fill( G, G+3*nat*3, 0. );
  double dEdA[9];
  for ( int k=0; k<3; ++k ) // foreach eigval
    {
      double * g = G + 3*nat*k;
      for ( int j=0; j<3; ++j )
	for ( int i=0; i<3; ++i )
	  dEdA[i+j*3] = U[i+k*3]*U[j+k*3]; // dEigval/dAij
      for ( int a=0; a<nat; ++a )
	{
	  double x = c[0+a*3]-Rcom[0];
	  double y = c[1+a*3]-Rcom[1];
	  double z = c[2+a*3]-Rcom[2];

	  // dEigval/dx = dEigval/dAij * dAij/dx
	  g[0+a*3] = m[a]*( -y*(dEdA[1]+dEdA[3])
			    -z*(dEdA[2]+dEdA[6])
			    +2*x*(dEdA[4]+dEdA[8]) ) ;

	  g[1+a*3] = m[a]*( 2*y*(dEdA[0]+dEdA[8])
			    -x*(dEdA[1]+dEdA[3])
			    -z*(dEdA[5]+dEdA[7]) ) ;

	  g[2+a*3] = m[a]*( 2*z*(dEdA[0]+dEdA[4])
			    -x*(dEdA[2]+dEdA[6])
			    -y*(dEdA[5]+dEdA[7]) ) ;

	}
    }
}

// -----------------------------------------------------------------
// inertia matrix eigenvectors
// -----------------------------------------------------------------
// Rotating the coordinates by 
//   c' = U^t . (c-com) + com
// aligns c' to the principle moments of inertia
void ccdl::InertiaAxis( int const nat, double const * m, double const * c, double * U )
{
  ccdl::InertiaTensor(nat,m,c,U);
  int N=3, LWORK=16, INFO=0;
  double WORK[16],E[3];
  dsyev_("V","U",&N,U,&N,E,WORK,&LWORK,&INFO);
}

// -----------------------------------------------------------------
// computes the 3 euler angles describing the orientation of the
// principle moments of inertia
// These angles are computed from the eigenvector matrix of the
// inertia tensor 
// -----------------------------------------------------------------
void ccdl::InertiaOrientation( int const nat, double const * m, double const * c, double * euler )
{
  double U[9];
  ccdl::InertiaAxis(nat,m,c,U);
  ccdl::ZYZFromRotationMatrix( U, euler[0], euler[1], euler[2] );
}

// -----------------------------------------------------------------
// computes the moments of inertia about the axis described by the
// input euler angles
// -----------------------------------------------------------------
void ccdl::InertiaMomentsAroundOrientation( int const nat, double const * m, double const * c, double const * euler, double * E )
{
  double I[9];
  ccdl::InertiaTensor(nat,m,c,I);
  double R[9];
  ccdl::RotationMatrixFromZYZ( euler[0],euler[1],euler[2], R );
  std::fill(E,E+3,0.);
  for ( int k=0; k<3; ++k )
    for ( int j=0; j<3; ++j )
      for ( int i=0; i<3; ++i )
	E[k] += R[i+k*3] * I[i+j*3] * R[j+k*3];
}

// -----------------------------------------------------------------
// computes the moments of inertia about the axis described by the
// input euler angles, and computes the derivative of the inertia
// wrt the cartesian position
// E is len 3
// G is len (3*nat)*3 -- the slow index runs over axis
// -----------------------------------------------------------------
void ccdl::InertiaMomentsAroundOrientation( int const nat, double const * m, double const * c, double const * euler, double * E, double * G )
{
  std::fill(E,E+3,0.);
  std::fill(G,G+3*nat*3,0.);
  double R[9];
  ccdl::RotationMatrixFromZYZ( euler[0],euler[1],euler[2], R );

  double Rcom[3] = {0.,0.,0.};
  double mtot = 0.;
  for ( int a=0; a<nat; ++a )
    {
      mtot += m[a];
      Rcom[0] += m[a] * c[0+a*3];
      Rcom[1] += m[a] * c[1+a*3];
      Rcom[2] += m[a] * c[2+a*3];
    };
  Rcom[0] /= mtot;
  Rcom[1] /= mtot;
  Rcom[2] /= mtot;

  for ( int k=0; k<3; ++k )
    {
      double * g = G + 3*nat*k;
      double rx = R[0+k*3];
      double ry = R[1+k*3];
      double rz = R[2+k*3];
      for ( int a=0; a<nat; ++a )
	{
	  double x = c[0+a*3]-Rcom[0];
	  double y = c[1+a*3]-Rcom[1];
	  double z = c[2+a*3]-Rcom[2];

	  g[0+a*3] = m[a] * 
	    (  rx * ( (  0)*rx+( -y)*ry+( -z)*rz )
	       +
	       ry * ( ( -y)*rx+(2*x)*ry+(  0)*rz )
	       +
	       rz * ( ( -z)*rx+(  0)*ry+(2*x)*rz ) );

	  g[1+a*3] = m[a] * 
	    (  rx * ( (2*y)*rx+( -x)*ry+(  0)*rz )
	       +
	       ry * ( ( -x)*rx+(  0)*ry+( -z)*rz )
	       +
	       rz * ( (  0)*rx+( -z)*ry+(2*y)*rz ) );

	  g[2+a*3] = m[a] * 
	    (  rx * ( (2*z)*rx+(  0)*ry+( -x)*rz )
	       +
	       ry * ( (  0)*rx+(2*z)*ry+( -y)*rz )
	       +
	       rz * ( ( -x)*rx+( -y)*ry+(  0)*rz ) );

	};
    };
  double I[9];
  std::fill(I,I+9,0.);
  for ( int a=0; a<nat; ++a )
    {
      double x = c[0+a*3]-Rcom[0];
      double y = c[1+a*3]-Rcom[1];
      double z = c[2+a*3]-Rcom[2];

      I[0] += m[a]*( y*y + z*z );
      I[1] -= m[a]*x*y;
      I[2] -= m[a]*x*z;
      I[4] += m[a]*( x*x + z*z );
      I[5] -= m[a]*y*z;
      I[8] += m[a]*( x*x + y*y );
    };
  I[3] = I[1];
  I[6] = I[2];
  I[7] = I[5];
  for ( int k=0; k<3; ++k )
    for ( int j=0; j<3; ++j )
      for ( int i=0; i<3; ++i )
	E[k] += R[i+k*3] * I[i+j*3] * R[j+k*3];
}



// -----------------------------------------------------------------
// computes the moments of inertia about the axis described by the
// input euler angles, and computes the derivative and hessian of
// the inertia wrt the cartesian position
// E is len 3
// G is len (3*nat)*3 -- the slow index runs over axis
// H in len (3*nat)*(3*nat)*3 -- each axis has a 3nat X 3nat Hessian
// -----------------------------------------------------------------
void ccdl::InertiaMomentsAroundOrientation( int const nat, double const * m, double const * c, double const * euler,  double * E, double * G, double * H )
{
  std::fill(E,E+3,0.);
  std::fill(G,G+3*nat*3,0.);
  std::fill(H,H+3*nat*3*nat*3,0.);
  double R[9];
  ccdl::RotationMatrixFromZYZ( euler[0],euler[1],euler[2], R );

  double Rcom[3] = {0.,0.,0.};
  double mtot = 0.;
  for ( int a=0; a<nat; ++a )
    {
      mtot += m[a];
      Rcom[0] += m[a] * c[0+a*3];
      Rcom[1] += m[a] * c[1+a*3];
      Rcom[2] += m[a] * c[2+a*3];
    };
  Rcom[0] /= mtot;
  Rcom[1] /= mtot;
  Rcom[2] /= mtot;

  for ( int k=0; k<3; ++k )
    {
      double * g = G + (3*nat)*k;
      double * h = H + (3*nat)*(3*nat)*k;
      double rx = R[0+k*3];
      double ry = R[1+k*3];
      double rz = R[2+k*3];
      for ( int a=0; a<nat; ++a )
	{
	  double x = c[0+a*3]-Rcom[0];
	  double y = c[1+a*3]-Rcom[1];
	  double z = c[2+a*3]-Rcom[2];

	  g[0+a*3] = m[a] * 
	    (  rx * ( (  0)*rx+( -y)*ry+( -z)*rz )
	       +
	       ry * ( ( -y)*rx+(2*x)*ry+(  0)*rz )
	       +
	       rz * ( ( -z)*rx+(  0)*ry+(2*x)*rz ) );

	  g[1+a*3] = m[a] * 
	    (  rx * ( (2*y)*rx+( -x)*ry+(  0)*rz )
	       +
	       ry * ( ( -x)*rx+(  0)*ry+( -z)*rz )
	       +
	       rz * ( (  0)*rx+( -z)*ry+(2*y)*rz ) );

	  g[2+a*3] = m[a] * 
	    (  rx * ( (2*z)*rx+(  0)*ry+( -x)*rz )
	       +
	       ry * ( (  0)*rx+(2*z)*ry+( -y)*rz )
	       +
	       rz * ( ( -x)*rx+( -y)*ry+(  0)*rz ) );


	  double hxx = (  rx * ( (  0)*rx+(  0)*ry+(  0)*rz )
			  +
			  ry * ( (  0)*rx+(  2)*ry+(  0)*rz )
			  +
			  rz * ( (  0)*rx+(  0)*ry+(  2)*rz ) );

	  double hyx = (  rx * ( (  0)*rx+( -1)*ry+(  0)*rz )
	       +
	       ry * ( ( -1)*rx+(  0)*ry+(  0)*rz )
	       +
	       rz * ( (  0)*rx+(  0)*ry+(  0)*rz ) );
	  double hzx = (  rx * ( (  0)*rx+(  0)*ry+( -1)*rz )
	       +
	       ry * ( (  0)*rx+(  0)*ry+(  0)*rz )
	       +
	       rz * ( ( -1)*rx+(  0)*ry+(  0)*rz ) );
	  double hxy = (  rx * ( (  0)*rx+( -1)*ry+(  0)*rz )
	       +
	       ry * ( ( -1)*rx+(  0)*ry+(  0)*rz )
	       +
	       rz * ( (  0)*rx+(  0)*ry+(  0)*rz ) );
	  double hyy = (  rx * ( (  2)*rx+(  0)*ry+(  0)*rz )
	       +
	       ry * ( (  0)*rx+(  0)*ry+(  0)*rz )
	       +
	       rz * ( (  0)*rx+(  0)*ry+(  2)*rz ) );
	  double hzy = (  rx * ( (  0)*rx+(  0)*ry+(  0)*rz )
	       +
	       ry * ( (  0)*rx+(  0)*ry+( -1)*rz )
	       +
	       rz * ( (  0)*rx+( -1)*ry+(  0)*rz ) );
	  double hxz = (  rx * ( (  0)*rx+(  0)*ry+( -1)*rz )
	       +
	       ry * ( (  0)*rx+(  0)*ry+( -0)*rz )
	       +
	       rz * ( ( -1)*rx+(  0)*ry+(  0)*rz ) );
	  double hyz = (  rx * ( (  0)*rx+(  0)*ry+(  0)*rz )
	       +
	       ry * ( (  0)*rx+(  0)*ry+( -1)*rz )
	       +
	       rz * ( (  0)*rx+( -1)*ry+(  0)*rz ) );
	  double hzz = (  rx * ( (  2)*rx+(  0)*ry+(  0)*rz )
	       +
	       ry * ( (  0)*rx+(  2)*ry+(  0)*rz )
	       +
	       rz * ( (  0)*rx+(  0)*ry+(  0)*rz ) );


	  for ( int b=0; b < nat; ++b )
	    {
	      double mm = (a==b ? m[a]*(1.-m[b]/mtot) : -m[a]*m[b]/mtot );
	      h[ (0+b*3) + (0+a*3)*(3*nat) ] = mm * hxx;
	      h[ (1+b*3) + (0+a*3)*(3*nat) ] = mm * hyx;
	      h[ (2+b*3) + (0+a*3)*(3*nat) ] = mm * hzx;

	      h[ (0+b*3) + (1+a*3)*(3*nat) ] = mm * hxy;
	      h[ (1+b*3) + (1+a*3)*(3*nat) ] = mm * hyy;
	      h[ (2+b*3) + (1+a*3)*(3*nat) ] = mm * hzy;

	      h[ (0+b*3) + (2+a*3)*(3*nat) ] = mm * hxz;
	      h[ (1+b*3) + (2+a*3)*(3*nat) ] = mm * hyz;
	      h[ (2+b*3) + (2+a*3)*(3*nat) ] = mm * hzz;
	    };
	}
    };
  double I[9];
  std::fill(I,I+9,0.);
  for ( int a=0; a<nat; ++a )
    {
      double x = c[0+a*3]-Rcom[0];
      double y = c[1+a*3]-Rcom[1];
      double z = c[2+a*3]-Rcom[2];

      I[0] += m[a]*( y*y + z*z );
      I[1] -= m[a]*x*y;
      I[2] -= m[a]*x*z;
      I[4] += m[a]*( x*x + z*z );
      I[5] -= m[a]*y*z;
      I[8] += m[a]*( x*x + y*y );
    };
  I[3] = I[1];
  I[6] = I[2];
  I[7] = I[5];
  for ( int k=0; k<3; ++k )
    for ( int j=0; j<3; ++j )
      for ( int i=0; i<3; ++i )
	E[k] += R[i+k*3] * I[i+j*3] * R[j+k*3];
}



