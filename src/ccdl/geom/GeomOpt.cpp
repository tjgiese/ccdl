#include "GeomOpt.hpp"
#include "../bmath.hpp"

ccdl::gopt::StepInfo::StepInfo ()
  : de(0.), de_abs(0.),
    gc_rms(0.), gc_max( -1.e+30 ),
    dxc_rms(0.), dxc_max( -1.e+30 ), dxc_len(0.)
{}

void ccdl::gopt::StepInfo::CptInfo
( ccdl::gopt::Step & step, 
  ccdl::gopt::Step & prevstep)
{
  de      = 0.; 
  de_abs  = 0.;
  gc_rms  = 0.; 
  gc_max  = -1.e+30;
  dxc_rms = 0.; 
  dxc_max = -1.e+30;
  dxc_len = 0.;

  int nat = step.nat;

  de = step.e - prevstep.e;
  de_abs = std::abs(de);

  step.dxc.resize( 3*nat );
  step.dgc.resize( 3*nat );

  step.GrdTransform( true );

  for ( int a=0; a<nat; ++a )
    {
      double xcnrm = 0.;
      double gcnrm = 0.;
      for ( int k=0; k<3; ++k )
	{
	  int i=k+a*3;
	  step.dxc[i] = step.x[i]-prevstep.x[i];
	  step.dgc[i] = step.g[i]-prevstep.g[i];
	  xcnrm += step.dxc[i]*step.dxc[i];
	  gcnrm += step.g[i]*step.g[i];
	};
      dxc_len += xcnrm;
      dxc_rms += xcnrm;
      gc_rms  += gcnrm;
      dxc_max  = std::max( dxc_max, std::sqrt( xcnrm ) );
      gc_max   = std::max( gc_max,  std::sqrt( gcnrm ) );
    };
  dxc_rms = std::sqrt( dxc_rms / step.n );
  dxc_len = std::sqrt( dxc_len );
  gc_rms  = std::sqrt( gc_rms / step.n );

  if ( step.dlc != NULL )
    {
      step.dxq.resize( step.nq );
      step.dgq.resize( step.nq );
      for ( int i=0; i<step.nq; ++i )
	{
	  step.dgq[i] = step.gq[i] - prevstep.gq[i];
	  step.dxq[i] = step.xq[i] - prevstep.xq[i];
	}
      //step.dlc->CptDifference( step.xq.data(), prevstep.xq.data(), step.dxq.data() );
    }
}

// ccdl::gopt::Step::Step
// ( int nat, double const * crd, 
//   ccdl::OptOptions & options )
//   : nat(nat), n(nat*3), nq(0), nmax( nat*3 ),
//     maxstep( options.maxstep ),
//     e(0.), 
//     predicted_de(0.),
//     x(crd,crd+n),
//     g(n,0.), 
//     h(n*n,0.),
//     eigvals( nmax, 0. ),
//     opts( &options ),
//     dlc( NULL ),
//     prevstep( NULL )
// {
// }

ccdl::gopt::Step::Step
( int nat, double const * crd,
  ccdl::OptOptions & options, ccdl::RedundantIC * dlcobj )
  : nat( nat ), 
    n( nat*3 ), 
    nq( 0 ),
    nmax( nat*3 ),
    maxstep( options.maxstep ),
    e(0.), 
    predicted_de(0.),
    x(crd,crd+n),
    g(n,0.), 
    h(n*n,0.),
    opts( &options ),
    dlc( dlcobj ),
    prevstep( NULL )
{
  if ( dlc != NULL )
    {
      nq = dlc->GetNumInternalCrds();
      nmax = std::max(n,nq);
      xq.assign(nq,0.);
      gq.assign(nq,0.);
      hq.assign(nq*nq,0.);
      dxq.assign(nq,0.);
      dgq.assign(nq,0.);
      dlc->DisplaceByDeltaQ( dxq.data(), x.data(), 1.e-8 );
      dlc->CptInternalCrds( x.data(), xq.data() );
      for ( int i=0; i<nq; ++i )
	hq[i+i*nq] = 1.;
    }
  eigvals.assign( nmax, 0. );
  for ( int i=0; i<n; ++i )
    h[i+i*n] = 1.;
}

void ccdl::gopt::Step::GrdTransform( bool q2c )
{
  if ( dlc == NULL ) return;

  // std::cerr << "Pre transform gradients\n";
  // std::vector<int> z(nat,1);
  // ccdl::WriteXyz( std::cerr, nat, z.data(), g.data() );

  dlc->GrdTransform( x.data(), g.data(), 
		     gq.data(), q2c );

  //std::cerr << "Post transform gradients\n";
  //ccdl::WriteXyz( std::cerr, nat, z.data(), g.data() );

}

void ccdl::gopt::Step::GrdAndHesTransform( bool q2c )
{
  if ( dlc == NULL ) return;
  dlc->GrdAndHesTransform( x.data(), g.data(), h.data(), 
			   gq.data(), hq.data(), q2c );
}

void ccdl::gopt::Step::CptInfo( ccdl::gopt::Step & prev )
{
  prevstep = &prev;
  info.CptInfo( *this, *prevstep );
}


void ccdl::gopt::Step::Print
( std::ostream & cout, 
  int const iter )
{

#define FMTF(a,b) std::fixed      << std::setw(a) << std::setprecision(b)
#define FMTE(a,b) std::scientific << std::setw(a) << std::setprecision(b)
#define FMT1(str,a,b) (str) << FMTF(14,8) << (a) << FMTE(14,4) << (b) << ( std::abs((a)) > (b) ? " F\n" : " T\n" )
#define FMT2(str,a)   (str) << FMTE(14,4) << (a) << "\n"
#define FMT3(str,a,b) (str) << FMTE(14,4) << (a) << FMTE(14,4) << (b) << ( std::abs((a)) > (b) ? " F\n" : " T\n" )


  cout << "\n\n";
  cout << "GEOMOPT ITER " << std::setw(4) << iter << "\n";

  if ( iter > 0 )
    cout << "GEOMOPT SUMMARY "
	 << std::setw(4) << iter
	 << FMTE(16,8) << e
	 << FMTE(13,4) << maxstep
	 << FMTE(13,4) << info.gc_max
	 << FMTE(13,4) << info.gc_rms
	 << "\n";
  
  if ( dlc != NULL )
    dlc->PrintReport( cout, x.data() );
  
  cout << "GEOMOPT Energy " << FMTF(20,8) << e << "\n";
  cout << FMT1("GEOMOPT MAX Force    ",info.gc_max,opts->gmax_tol);
  cout << FMT1("GEOMOPT RMS Force    ",info.gc_rms,opts->grms_tol);
  if ( iter > 0 )
    {
      cout << FMT1("GEOMOPT MAX Displ    ",info.dxc_max,opts->xmax_tol);
      cout << FMT1("GEOMOPT RMS Displ    ",info.dxc_rms,opts->xrms_tol);
      cout << FMT3("GEOMOPT Actual dE    ",info.de,opts->ener_tol);
      cout << "GEOMOPT Step Length  "
	   << FMTE(14,4) << info.dxc_len
	   << FMTE(14,4) << maxstep << " max\n";    
    };

}


bool ccdl::gopt::Step::CheckConvergence()
{
  bool CONVERGED = false;
  CONVERGED = ( info.gc_max <= opts->gmax_tol and
		info.gc_rms <= opts->grms_tol and
		info.dxc_max <= opts->xmax_tol and
		info.dxc_rms <= opts->xrms_tol and
		info.de_abs  <= opts->ener_tol );

  if ( ! CONVERGED )
    CONVERGED = ( ( info.gc_max < 0.05 * opts->gmax_tol and 
		    info.de_abs  < opts->ener_tol ) or
		  ( info.gc_max < opts->gmax_tol and 
		    info.de_abs < 1.e-12 ) );

  if ( ! CONVERGED )
    CONVERGED = ( info.gc_max < 0.01 * opts->gmax_tol );

  return CONVERGED;
}




void ccdl::gopt::Step::UpdateHessian()
{
  int ncrd = n;

  double * Hold = prevstep->h.data();
  double * Hnew = h.data();
  double * DG = dgc.data();
  double * DX = dxc.data();

  if ( dlc != NULL )
    {
      ncrd = nq;
      Hold = prevstep->hq.data();
      Hnew = hq.data();
      DG = dgq.data();
      DX = dxq.data();
    };

  double dxnrm = ccdl::v_dot_v( ncrd, DX, DX );
  if ( dxnrm > 1.e-15 )
    {
      switch ( opts->update )
	{
	case ccdl::gopt::BFGS:
	  // MIN SEARCH
	  ccdl::gopt::UpdateHessian_BFGS
	    ( ncrd, Hold, DG, DX, Hnew );
	  break;
	  
	case ccdl::gopt::MS:
	  // MIN SEARCH
	  ccdl::gopt::UpdateHessian_MS
	    ( ncrd, Hold, DG, DX, Hnew );
	  break;
	  
	case ccdl::gopt::BFGSMS:
	  // MIN SEARCH
	  ccdl::gopt::UpdateHessian_Bofill_BFGS_MS
	    ( ncrd, Hold, DG, DX, Hnew );
	  break;
	  
	case ccdl::gopt::PSB:
	  // TS SEARCH
	  ccdl::gopt::UpdateHessian_PSB
	    ( ncrd, Hold, DG, DX, Hnew );
	  break;
	  
	case ccdl::gopt::PSBMS:
	  // TS SEARCH
	  ccdl::gopt::UpdateHessian_Bofill_PSB_MS
	    ( ncrd, Hold, DG, DX, Hnew );
	  break;
	  
	default:
	  break;
	}
    }

  if ( dlc != NULL ) dlc->HesBackTransform( x.data(), gq.data(), hq.data(), h.data() );


  std::printf("Hold\n");
  for ( int i=0; i<ncrd; ++i )
    {
      for ( int j=0; j<ncrd; ++j )
  	std::printf("%12.3e",Hold[i+j*ncrd]);
      std::printf("\n");
    };
  std::printf("DG\n");
  for ( int i=0; i<ncrd; ++i )
    std::printf("%12.3e",DG[i]);
  std::printf("\n");

  std::printf("DX\n");
  for ( int i=0; i<ncrd; ++i )
    std::printf("%12.3e",DX[i]);
  std::printf("\n");

  std::printf("Hnew\n");
  for ( int i=0; i<ncrd; ++i )
    {
      for ( int j=0; j<ncrd; ++j )
  	std::printf("%12.3e",Hnew[i+j*ncrd]);
      std::printf("\n");
    };


}



void ccdl::gopt::Step::Move( double & maxstep )
{
  maxstep = maxstep;

  int N = n;
  double * H = h.data();
  double * G = g.data();
  double * DX = dxc.data();
  if ( dlc != NULL )
    {
      N = nq;
      H = hq.data();
      G = gq.data();
      DX = dxq.data();
    }


  switch ( opts->type )
    {
    case ccdl::gopt::MIN:
      {
	
	// for ( int i=0; i<N; ++i )
	//   {
	//     for ( int j=0; j<N; ++j )
	//       std::printf("%12.3e",H[i+j*N]);
	//     std::printf("\n");
	//   };
	    
	ccdl::gopt::CptDeltaX_TrustRadius
	  ( N, H, G, maxstep, DX, eigvals.data() );
      };
      break;
      
    case ccdl::gopt::TS: // need to add following options
      ccdl::gopt::CptDeltaX_EigenFollow
	( N, H, G, maxstep, DX, eigvals.data() );
      break;
      
    default:
      break;
    }

  x = prevstep->x;
  xq = prevstep->xq;
  if ( dlc == NULL )
    {
      for ( int i=0; i<n; ++i ) x[i] += dxc[i];
    }
  else
    {
      dlc->DisplaceByDeltaQ( dxq.data(), x.data(), 1.e-8 );
      dlc->CptInternalCrds( x.data(), xq.data() );
    };
}






ccdl::OptOptions::OptOptions()
  : type( ccdl::gopt::MIN ),
    update( ccdl::gopt::BFGS ),
    delta( 1.e-4 ),
    calcfc( false ),
    calcall( false ),
    maxiter( 100 ),
    maxstep( 0.5 ),
    limstep( 1.0 ),
    varmaxstep( true ),
    eigvec( -1 ),
    ener_tol( 1.e-8 ),
    gmax_tol( 1.e-4 ),
    xmax_tol( 1.e-4 ),
    grms_tol( 1.e-4 ),
    xrms_tol( 1.e-4 ),
    ostr( &std::cout )
{}

