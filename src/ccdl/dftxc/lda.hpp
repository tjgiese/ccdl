#ifndef _XCF_lda_H_
#define _XCF_lda_H_

extern "C"
{

  void uks_xf_lda_spw92
  ( double const * pa, double const * pb,
    double * F );
  
  void uks_cf_lda_spw92
  ( double const * pa, double const * pb,
    double * F );
  
  void uks_xcf_lda_spw92_potential
  ( double const * pa, double const * pb,
    double * F );
  
  void rks_xcf_lda_spw92_potential
  ( double const * pa, 
    double * F );

}

#endif

