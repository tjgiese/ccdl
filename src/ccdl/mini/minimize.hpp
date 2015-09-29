#ifndef _ccdl_minimizer_hpp_
#define _ccdl_minimizer_hpp_

#include <iostream>
#include "minifcn.hpp"

// see minimize.cpp for examples / test functions


namespace ccdl
{

  namespace mini
  {
    struct options
    {
      options()
	: maxiter(1000),
	  nsd(1),
	  ncg(3),
	  gmax_tol( 1.e-4 ),
	  grms_tol( 5.e-5 ),
	  ostr( &std::cout )
      {}

      int maxiter,nsd,ncg;
      double gmax_tol,grms_tol;
      std::ostream * ostr;
    };
  }

  int minimize( ccdl::minifcn * pfcn, ccdl::mini::options opts = ccdl::mini::options() );
}




#endif
