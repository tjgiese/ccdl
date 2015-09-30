#ifndef _ccdl_lsfit_hpp_
#define _ccdl_lsfit_hpp_

#include "exceptions.hpp"

namespace ccdl
{

  int LeastSquaresFit
  ( int const nobs, int const nparam,
    double const * A_obs_by_param,
    double * x_param,
    double const * b_obs,
    double relative_accuracy_of_the_obs = 1.e-32 );

  int WeightedLeastSquaresFit
  ( int const nobs, int const nparam,
    double const * A_obs_by_param,
    double * x_param,
    double const * b_obs,
    double const * w_obs );

  int ConstrainedLeastSquaresFit
  ( int const nobs, int const nparam,
    double const * A_obs_by_param,
    double * x_param,
    double const * b_obs,
    int const ncon,
    double const * D_con_by_param,
    double const * c_con );

  int ConstrainedWeightedLeastSquaresFit
  ( int const nobs, int const nparam,
    double const * A_obs_by_param,
    double * x_param,
    double const * b_obs,
    double const * w_obs,
    int const ncon,
    double const * D_con_by_param,
    double const * c_con );

}

#endif

