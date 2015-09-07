#ifndef _ccdl_goptprimitives_hpp_
#define _ccdl_goptprimitives_hpp_

namespace ccdl
{
  namespace gopt
  {

    double CptDeltaX_RFO
    ( int const n, double const * H, double const * g,
      double const trust_radius, double * dx );

    double CptDeltaX_TrustRadius
    ( int const n, double const * H, double const * g,
      double const trust_radius, double * dx, double * evals );

    double CptDeltaX_EigenFollow
    ( int const n, double const * H, double const * g, 
      double trust_radius, double * dx, double * evals );


    double PredictEnergyChange
    ( int const n, double const * H, double const * g, 
      double const * dx );

    void UpdateHessian_BFGS
    ( int const n, double const * H, double const * g, 
      double const * dx, double * Hnew );

    void UpdateHessian_MS
    ( int const n, double const * H, double const * g, 
      double const * dx, double * Hnew );

    void UpdateHessian_Bofill_BFGS_MS
    ( int const n, double const * H, double const * g, 
      double const * dx, double * Hnew );

    void UpdateHessian_PSB
    ( int const n, double const * H, double const * g, 
      double const * dx, double * Hnew );

    void UpdateHessian_Bofill_PSB_MS
    ( int const n, double const * H, double const * g, 
      double const * dx, double * Hnew );

    // recommended
    int CptDeltaX_dsysv
    ( int const n, double const * H, double const * g, 
      double const lambda, double * dx );

    // slower
    int CptDeltaX_dposv
    ( int const n, double const * H, double const * g, 
      double const lambda, double * dx );

    // this is super-slow
    int CptDeltaX_dpotrf
    ( int const n, double const * H, double const * g, 
      double const lambda, double * dx );



    // this is a debug version -- don't use
    void CptDeltaX_eigen
    ( int const n, double const * H, double const * g, 
      double const lambda, double * dx );

    // this is a debug version -- don't use
    double CptDeltaXNorm
    ( int const n, double const * H, double const * g, 
      double const lambda  );

  }
}
#endif
