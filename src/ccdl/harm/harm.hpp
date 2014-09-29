#ifndef _ccdl_harm_hpp_
#define _ccdl_harm_hpp_

namespace ccdl
{
  void SolidHarmRlm_v5( int const Lmax, double const * R, double * Y );
  void SolidHarmIlm_v5( int const Lmax, double const * R, double * Y );
  void RlmTranslation( int const Lto, int const Lfrom, 
		       double const * Rtf, double * Wtf );
  void IlmInteraction( int const La, int const Lb, 
		       double const * Rab, double * T );

}

#endif
