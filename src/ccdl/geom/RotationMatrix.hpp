#ifndef _ccdl_rotationmatrix_hpp_
#define _ccdl_rotationmatrix_hpp_

namespace ccdl
{
  //
  // Given Euler angles Z,Y,Zp; this routine constructs
  // a rotation matrix 
  //    R(Z,Y,Zp) = R(Zp).R(Y).R(Z)
  // So that a rotated vector, r' is obtained by right-multiplication
  //    r' = R.r
  //
  // That is, the first rotation is about Z and the last rotation
  // is about Zp
  //
  void RotationMatrixFromZYZ
  ( double Z, double Y, double Zp, double * R );

  //
  // The gradients dR/dZ, dR/dY and dR/dZp
  //
  void RotationMatrixGradientsFromZYZ
  ( double  Z, double  Y, double Zp, 
    double * dRdZ, double * dRdY, double * dRdZp );

  //
  // Obtains the z-y-z Euler angles from a rotation matrix.
  // Note that the output Euler angles can be used to
  // reconstruct the rotation matrix EXACTLY; however,
  // these Euler angles are not unique!  There is another set
  // of z-y-z Euler angles that could also reproduce
  // this rotation matrix.
  //
  void ZYZFromRotationMatrix
  ( double const * R, double & Z, double & Y, double & Zp );


  //
  // Returns the rotation and translations required to
  // overlay crd onto rcrd, the reference coordinates.
  //
  // W    = weights
  // rcrd = reference crds
  // crd  = unrotated crds
  // R    = a 3x3 rotation matrix
  // rcom = the "center of mass" of rcrd
  // com  = the "center of mass" of crd
  //
  // You can rotate crd using R, rcom, and com, as follows:
  //
  // new_crd[k+a*3] = rcom[k] + \sum_kp R[k+kp*3] * ( crd[kp+a*3]-com[kp] )
  //
  double RmsTransformData
  ( int n, double const * W,
    double const * rcrd, double const * crd,
    double * R, double * rcom, double * com );
    
  //
  // Given the RmsTransformData, perform the rotation on crd,
  // so that, on output:
  //
  // crd[k+a*3] => rcom[k] + \sum_kp R[k+kp*3] * ( crd[kp+a*3]-com[kp] )
  //
  void RotateCrds
  ( int nat, double * crd, 
    double const * R, double const * rcom, double const * com );

  //
  // Given the RmsTransformData, perform the reverse rotation on
  // the rotated crds to return the original set of crds
  //
  void UnrotateCrds
  ( int nat, double * crd,
    double const * R, double const * rcom, double const * com );

  //
  // RMS overlay a set of crds onto a set of reference crds (rcrd).
  // The crd array is overwritten with the rotated crds,
  // and you are NOT provided the RmsTransformData;
  // therefore, you should not use this if you need to rotate back!
  //
  void RmsOverlay
  ( int nat, double const * w,
    double const * rcrd, double * crd );

  //
  // Unrotates a set of gradients that were computed in the
  // rotated frame (rotated_grd).
  //
  //   w             = the weights used to perform the rms fit
  //   crd           = the set of unrotated coordinates.
  //   rcrd          = the set of reference coordinates; the overlay of crd
  //                   onto rcrd defines the rotation
  //   rotated_grd   = the gradients computed in the rotated frame
  //   unrotated_grd = the gradients consistent with the
  //                   orientation of crd
  //
  void RmsOverlayGrd
  ( int nat, double const * w,
    double const * rcrd, double const * crd,
    double const * rotated_grd, double * unrotated_grd );


}





//
// debug / testing routines
//
namespace ccdl
{
  namespace RotationMatrix
  {
    void test_RmsOverlayGrd
    ( int n, double const * W,
      double const * rcrd, double * crd );
  }
}


#endif
