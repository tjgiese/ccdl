#ifndef _ccdl_InertiaTensor_hpp_

namespace ccdl
{
  // -----------------------------------------------------------------
  // 3x3 moment of inertia matrix
  // -----------------------------------------------------------------
  void InertiaTensor( int const nat, double const * m, double const * c, double * I );
  // -----------------------------------------------------------------
  // inertia tensor eigenvalues
  // -----------------------------------------------------------------
  void InertiaMoments( int const nat, double const * m, double const * c, double * E );
  // -----------------------------------------------------------------
  // inertia tensor eigenvalues and gradients wrt c
  // -----------------------------------------------------------------
  void InertiaMoments( int const nat, double const * m, double const * c, double * E, double * G );
  // -----------------------------------------------------------------
  // inertia matrix eigenvectors
  // -----------------------------------------------------------------
  // Rotating the coordinates by 
  //   c' = U^t . (c-com) + com
  // aligns c' to the principle moments of inertia
  void InertiaAxis( int const nat, double const * m, double const * c, double * U );
  // -----------------------------------------------------------------
  // computes the 3 euler angles describing the orientation of the
  // principle moments of inertia
  // These angles are computed from the eigenvector matrix of the
  // inertia tensor 
  // -----------------------------------------------------------------
  void InertiaOrientation( int const nat, double const * m, double const * c, double * euler );
  // -----------------------------------------------------------------
  // computes the moments of inertia about the axis described by the
  // input euler angles
  // -----------------------------------------------------------------
  void InertiaMomentsAroundOrientation( int const nat, double const * m, double const * c, double const * euler, double * E );
  // -----------------------------------------------------------------
  // computes the moments of inertia about the axis described by the
  // input euler angles, and computes the derivative of the inertia
  // wrt the cartesian position
  // E is len 3
  // G is len (3*nat)*3 -- the slow index runs over axis
  // -----------------------------------------------------------------
  void InertiaMomentsAroundOrientation( int const nat, double const * m, double const * c, double const * euler, double * E, double * G );
  // -----------------------------------------------------------------
  // computes the moments of inertia about the axis described by the
  // input euler angles, and computes the derivative and hessian of
  // the inertia wrt the cartesian position
  // E is len 3
  // G is len (3*nat)*3 -- the slow index runs over axis
  // H in len (3*nat)*(3*nat)*3 -- each axis has a 3nat X 3nat Hessian
  // -----------------------------------------------------------------
  void InertiaMomentsAroundOrientation( int const nat, double const * m, double const * c, double const * euler,  double * E, double * G, double * H );

}

#endif

