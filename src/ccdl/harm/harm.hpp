#ifndef _ccdl_harm_hpp_
#define _ccdl_harm_hpp_

namespace ccdl
{
  void CptAlm( int const Lmax, double * Alm );
  void SolidHarm_Rlm( int const Lmax, double const * R, double * Y );
  void SolidHarm_dRlm( int const Lmax, double const * Y, double * dY );
  void SolidHarm_d2Rlm( int const Lmax, double const * dY, double * d2Y );
  void RlmTranslation( int const Lto, int const Lfrom, double const * Rft, double * W );
  void RlmTranslationGrd( int const Lto, int const Lfrom, double const * Rft, double * W, double * dW );
  void ClmTranslation( int const Lto, int const Lfrom, double const * Rft, double * W );


  void SolidHarm_Ilm( int const Lmax, double const * R, double * Y );
  void SolidHarm_dIlm( int const Lmax, double const * Y, double * dY );
  void SolidHarm_d2Ilm( int const Lmax, double const * dY, double * d2Y );
  void IlmInteraction( int const La, int const Lb, double const * Rab, double * T );
  void IlmInteractionGrd( int const La, int const Lb, double const * Rab, double * T, double * dT );

  // void ExpIntMat( int const La, int const Lb, double const * O, double * T, double * WORK );
  // void ExpGrdMat( int const La, int const Lb, double const * O, double * T, double * D, double * WORK );
  // void SegIntMat( int const La, int const Lb, double const * O, double * T, double * WORK );
  // void SegGrdMat( int const La, int const Lb, double const * O, double * T, double * D, double * WORK );

  
}

#endif
