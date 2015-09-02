#include "GeomOpt.hpp"


ccdl::gopt::Step::Step( int n )
  : n(n),
    e(0.), 
    predicted_de(0.),
    x(n,0.),
    g(n,0.), 
    h(n*n,0.)
{
}


ccdl::OptOptions::OptOptions()
  : type( ccdl::gopt::MIN ),
    update( ccdl::gopt::BFGS ),
    delta( 1.e-4 ),
    calcfc( false ),
    calcall( false ),
    maxiter( 100 ),
    maxstep( 0.5 ),
    varmaxstep( true ),
    eigvec( -1 ),
    ener_tol( 1.e-8 ),
    gmax_tol( 1.e-4 ),
    xmax_tol( 1.e-4 ),
    grms_tol( 1.e-4 ),
    xrms_tol( 1.e-4 ),
    ostr( &std::cout )
{}

