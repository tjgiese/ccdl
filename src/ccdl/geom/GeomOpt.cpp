#include "GeomOpt.hpp"
#include "../bmath.hpp"
#include <algorithm>


namespace ccdl
{
  namespace gopt
  {

    struct steppt;
    bool steps_rfo_min( ccdl::gopt::steppt const & a,
			ccdl::gopt::steppt const & b );


    struct steppt
    {

      int n;
      double const * h;
      double const * evals;
      double const * evecs;
      double const * g;
      std::vector<double> numer;
      double const * q0;
      double const * x0;
      std::vector<double> dq;
      std::vector<double> dx;
      ccdl::RedundantIC * dlc;

      double lambda;
      double de,rfo;
      double dqlen,dxlen;
      double dxmax;

      steppt() : n(0),h(NULL),evals(NULL),g(NULL),
		 q0(NULL),x0(NULL),dlc(NULL),
		 lambda(0.),de(0.),rfo(0.),
		 dqlen(0.),dxlen(0.),dxmax(0.) 
      {}

      steppt( int N, 
	      double const * H, 
	      double const * EVALS, 
	      double const * EVECS, 
	      double const * G, 
	      double const * Q0,
	      double const * X0,
	      ccdl::RedundantIC * DLC = NULL )
      {
	reset( N,H,EVALS,EVECS,G,Q0,X0,DLC );
      }

      void reset( int N, 
		  double const * H, 
		  double const * EVALS, 
		  double const * EVECS, 
		  double const * G, 
		  double const * Q0,
		  double const * X0,
		  ccdl::RedundantIC * DLC = NULL )
      {
	n=N;
	h=H; 
	evals=EVALS;
	evecs=EVECS;
	g=G;
	numer.assign(N,0.);
	q0=Q0;
	x0=X0;
	dq.assign(N,0.); 
	dlc=DLC;
	lambda=0.; 
	de=0.;
	rfo=0.; 
	dqlen=0.; 
	dxlen=0.; 
	dxmax=0.;
      
	if ( dlc == NULL )
	  dx=dq;
	else
	  dx.assign(3*dlc->GetNumAtoms(),0.);
	ccdl::gt_dot_v( n, n, evecs, g, numer.data() );
      }

      bool SetLambda( double L )
      {
	
	lambda = L;

	std::fill(dq.data(),dq.data()+n,0.);
	for ( int i=0; i<n; ++i )
	  {
	    double alpha = 0.;
	    if ( std::abs( evals[i]-lambda ) > 1.e-30 )
	      alpha = numer[i]/(lambda-evals[i]);
	    ccdl::axpy( alpha, n, evecs+i*n, dq.data() );
	  };


	double dq2 = ccdl::v_dot_v(n,dq.data(),dq.data());
	if ( dq2 > 1000. ) return false;

	std::vector<double> t(g,g+n);
	ccdl::sy_dot_v( n, n, h, dq.data(), t.data(), 0.5, 1. );
	de = ccdl::v_dot_v( n, dq.data(), t.data() );
	dqlen = std::sqrt( ccdl::v_dot_v( n, dq.data(), dq.data() ) );
	rfo = de / ( 1. + dqlen );

	int nat = n/3;
	if ( dlc == NULL )
	  {
	    dxlen = dqlen;
	    dx = dq;
	  }
	else
	  {
	    nat = dlc->GetNumAtoms();
	    t.resize( 3*nat );
	    std::copy( x0, x0 + 3*nat, t.data() );
	    std::vector<double> ddq(dq);
	    dlc->DisplaceByDeltaQ( ddq.data(), t.data() );
	    nat = dlc->GetNumAtoms();

	    // {
	    //   std::printf("%i\n",nat);
	    //   std::printf("%20.10f\n",lambda);
	    //   for ( int i=0; i<nat; ++i )
	    // 	{
	    // 	  std::printf("X ");
	    // 	  for ( int k=0; k<3; ++k )
	    // 	    std::printf("%20.10f",t[k+i*3]/ccdl::AU_PER_ANGSTROM);
	    // 	  std::printf("\n");
	    // 	}
	    // };

	    dxlen = 0.;
	    dxmax = 0.;
	    for ( int i=0; i<nat*3; ++i )
	      {
		dx[i] = t[i] - x0[i];
		dxlen += dx[i]*dx[i];		
	      }
	    dxlen = std::sqrt( dxlen );
	  }
	//std::printf("nat = %i\n",nat);
	dxmax = 0.;
	int i=0;
	for ( int a=0; a<nat; ++a )
	  {
	    double len2 = 0.;
	    for ( int k=0; k<3; ++k, ++i )
	      len2 += dx[i]*dx[i];
	    dxmax = std::max( dxmax, len2 );
	  };
	dxmax = std::sqrt( dxmax );
	return true;
      }

      double MinSeparation()
      {
	int nat = n/3;
	if ( dlc != NULL ) nat = dlc->GetNumAtoms();
	double minr2 = 1.e+30;
	for ( int a=1; a<nat; ++a )
	  {
	    double dxa = x0[0+a*3]+dx[0+a*3];
	    double dya = x0[1+a*3]+dx[1+a*3];
	    double dza = x0[2+a*3]+dx[2+a*3];
	    for ( int b=0; b<a; ++b )
	      {
		double DX = dxa-(x0[0+b*3]+dx[0+b*3]);
		double DY = dya-(x0[1+b*3]+dx[1+b*3]);
		double DZ = dza-(x0[2+b*3]+dx[2+b*3]);
		double r2 = DX*DX + DY*DY + DZ*DZ;
		minr2 = std::min( minr2, r2 );
	      }
	  }
	if ( nat < 2 )
	  minr2 = 100.;
	else
	  minr2 = std::sqrt( minr2 );
	return minr2;
      }

    };


    double FindStepLength
    ( int const ng, double const * H, double const * g, double const * x0,
      double const trust_radius, 
      double * dx, double * evals,
      ccdl::RedundantIC * dlc = NULL )
    {
      std::vector<double> q0;
      //int nat = ng/3;
      if ( dlc == NULL )
	q0.assign( x0, x0 + ng );
      else
	{
	  q0.assign( ng, 0. );
	  dlc->CptInternalCrds( x0, q0.data() );
	};
      std::vector<double> evecs(ng*ng,0.);
      ccdl::eigen( ng, H, evals, evecs.data() );
      std::vector< ccdl::gopt::steppt > steps(3);

      std::vector<double> tmpH;
      // if ( evals[0] < -100. )
      // 	{
      // 	  tmpH.assign(H,H+ng*ng);
      // 	  for ( int i=0; i<ng; ++i )
      // 	    tmpH[i+i*ng] += 10.;
      // 	  ccdl::eigen( ng, tmpH.data(), evals, evecs.data() );
      // 	  steps[0].reset( ng, tmpH.data(), evals, evecs.data(), g, q0.data(), x0, dlc );
      // 	}
      // else
      steps[0].reset( ng, H, evals, evecs.data(), g, q0.data(), x0, dlc );
      steps[1] = steps[0];
      steps[2] = steps[0];

     
      // forward search from the lowest eigenvalue
      // to find some lambda that does not exceed the trust_radius
      double lambda_never_gt = evals[0] - 1.e-8;
      if ( evals[0] < -100. )
	lambda_never_gt -= std::max( 1.e-8*std::abs(evals[0]), 1.e-6 );

      double lambda_hi = lambda_never_gt;
      double lambda_lo = lambda_hi;
      double delta =  std::max( std::abs(lambda_hi) , 5.e-5 );
      while ( true ) // first find a lambda that doesn't cause lapack to explode
	{
	  //std::printf("lambda %13.4e ",lambda_hi);
	  if ( steps[1].SetLambda( lambda_hi ) )
	    {
	      //std::printf("(%13.4e)\n",steps[1].dxmax);
	      if ( steps[1].dxmax < trust_radius ) 
		{ // we went too far -- do a backtrack
		  lambda_lo = lambda_hi;
		  delta /= 1.5;
		  if ( lambda_hi + delta >= lambda_never_gt )
		    {
		      lambda_hi = (lambda_hi + lambda_never_gt)/2.;
		      delta = 0.5 * std::abs(lambda_hi-lambda_lo);
		    }
		  else
		    lambda_hi += delta;
		}
	      else // lambda_hi is safe, but trust_radius is too high
		break;
	      if ( std::abs( lambda_lo - lambda_hi ) < 1.e-8 )
		break;
	    }
	  else
	    {
	      //std::printf("\n");
	      delta *= 2.;
	      lambda_lo = lambda_hi;
	      lambda_hi -= delta;
	    };
	}
      lambda_lo = steps[1].lambda;

      //std::printf("dxmax @ evals[0]  %13.4e (%13.4e)\n",lambda_lo,steps[1].dxmax);

      steps[0] = steps[1];
      if ( steps[0].dxmax > trust_radius )
	{
	  double dxprev = steps[0].dxmax;
	  //for ( int i=0; i<20; ++i )
	  while ( true )
	    {
	      steps[0] = steps[1];
	      if ( std::abs( lambda_lo ) > 0.1 ) 
		lambda_lo -= std::max(0.1,std::abs(lambda_lo));
	      else 
		lambda_lo -= 0.2;
	      steps[1].SetLambda( lambda_lo );	  

	      //std::printf("lambda %13.4e %13.4e\n",steps[1].lambda, steps[1].dxmax);

	      if ( steps[1].dxmax < trust_radius )
		break;
	      if ( steps[1].dxmax > dxprev )
		break;
	      dxprev = steps[1].dxmax;
	    }
	}
      if ( steps[1].dxmax > trust_radius )
	{
	  if ( dlc == NULL )
	    std::copy( steps[1].dx.data(), steps[1].dx.data() + ng, dx );
	  else
	    {
	      //int nat = dlc->GetNumAtoms();
	      //int nat3 = nat*3;
	      steps[0] = steps[1];
	      std::vector<double> ddq( steps[1].dq );
	      double scl = 0.9 * trust_radius / steps[0].dxmax;
	      {
		for ( std::size_t k=0; k<ddq.size(); ++k )
		  ddq[k] = scl * steps[0].dq[k];
		steps[1].dq = ddq;
		//steps[1].x.assign( x0, x0+nat3 );
		//dlc->DisplaceByDeltaQ( ddq, steps[1].x.data() );
	      }
	      std::copy( steps[1].dq.data(), steps[1].dq.data() + ng, dx );
	    }
	  return 1;
	  //std::printf("DXMAX TOO BIG %13.4e %13.4e\n",steps[1].dxmax, trust_radius);
	}

      //std::printf("trust between     %13.4e %13.4e (%13.4e %13.4e)  tr %13.4e\n",
      //	  steps[0].lambda, steps[1].lambda,
      //	  steps[0].dxmax, steps[1].dxmax,
      //	  trust_radius);
      // now bisect to the trust radius
      steps[2]  = steps[1];      
      if ( std::abs(steps[2].dxmax - trust_radius) > 1.e-6 )
	{
	  for ( int i=0; i<30; ++i )
	    {
	      steps[1].SetLambda( 0.5 * ( steps[0].lambda + steps[2].lambda ) );
	      if ( steps[1].dxmax > trust_radius )
		steps[0] = steps[1];
	      else
		steps[2] = steps[1];
	      if ( std::abs(steps[2].dxmax - trust_radius) < 1.e-6 )
		break;
	    }
	};

      //std::printf("trust @           %13.4e (%13.4e)\n",
      //	  steps[2].lambda, steps[2].dxmax );

      // steps[2] is now at the trust radius
      // We can now search for a lambda to minimize the rfo, as long as
      // lambda is lower than what steps[2] currently is 
      lambda_hi = steps[2].lambda;
      steps[1]  = steps[2];
      steps[0]  = steps[1];

      // double delta = 0.05 * std::abs( steps[1].lambda );
      // std::printf("lambda rfo %13.4e %13.4e\n",steps[1].lambda, steps[1].rfo );
      // while ( true )
      // 	{
      // 	  steps[1].SetLambda( steps[1].lambda - delta );
      // 	  std::printf("lambda rfo %13.4e %13.4e\n",steps[1].lambda, steps[1].rfo );
      // 	  if ( steps[1].rfo < steps[0].rfo )
      // 	    steps[0] = steps[1];
      // 	  else
      // 	    break;
      // 	}

      
      // forward search for an rfo that's less than the current rfo
      delta = 0.05 * std::abs(lambda_hi);
      lambda_hi -= delta;
      while ( true )
	{
	  steps[1] = steps[2];
	  steps[2].SetLambda( lambda_hi );
	  if ( steps[2].rfo >= steps[1].rfo )
	    { // we went too far
	      break;
	    }
	  else
	    {
	      delta *= 2.;
	      lambda_hi -= delta;
	    };
	}
      ccdl::gopt::steppt trust_step( steps[0] );
      // the rfo minimum is between steps[1] and steps[2]
      //std::printf("rfo between       %13.4e %13.4e (%13.4e %13.4e) [%13.4e %13.4e]\n",
      //	  steps[1].lambda, steps[2].lambda,
      //	  steps[1].dxmax, steps[2].dxmax,
      //	  steps[1].rfo,    steps[2].rfo  );

      if ( steps[1].dxmax < trust_radius and steps[2].dxmax < trust_radius )
	{
	  steps[0] = steps[1];
	  steps[1] = steps[2];
	  while ( true )
	    {
	      double lnew = 0.5 * ( steps[0].lambda + steps[1].lambda );
	      steps[2].SetLambda( lnew );
	      std::sort( steps.begin(), steps.end(), &ccdl::gopt::steps_rfo_min );
	      //      std::printf("rfo %13.5e %13.5e %13.5e (%13.5e %13.5e %13.5e)\n",
	      //	  steps[0].lambda,steps[1].lambda,steps[2].lambda,
	      //	  steps[0].rfo,steps[1].rfo,steps[2].rfo);
	      //if ( steps[2].SetLambda( lnew ) );
	      
	      if ( std::abs( steps[1].rfo-steps[0].rfo ) < 1.e-9 or
		   std::abs( steps[1].lambda-steps[0].lambda ) < 1.e-8 ) break;
	    }
	  
	  //std::printf("rfo   @           %13.4e\n",
	  //	      steps[0].lambda );
	}
      else
	{
	  //std::printf("rfo search skipped\n");
	}
      // steps[0] now has the minimum rfo below the trust_radius
      // Did we make any close contacts?
      // Decrease lambda until the close contacts are removed
      double min0 = steps[0].MinSeparation();
      if ( min0 < 0.6 * ccdl::AU_PER_ANGSTROM )
	{
	  steps[1] = steps[0];
	  for ( int i=0; i<10; ++i )
	    {
	      steps[1].SetLambda( steps[1].lambda + delta );
	      min0 = steps[0].MinSeparation();
	      if ( min0 > 0.6 * ccdl::AU_PER_ANGSTROM )
		break;
	    }
	  steps[0] = steps[1];
	}
      // steps[0] now avoids excessively close contacts
      // interpret the trust_radius as the dxmax


      //std::printf("close @           %13.4e (%13.4e) [%13.4e]\n",
      //	  steps[0].lambda,steps[0].dxmax,min0/ccdl::AU_PER_ANGSTROM );


      if ( dlc == NULL )
	std::copy( steps[0].dx.data(), steps[0].dx.data() + ng, dx );
      else
	std::copy( steps[0].dq.data(), steps[0].dq.data() + ng, dx );


      {
	

      }

      return steps[0].dxmax;
    }



  }
}

bool ccdl::gopt::steps_rfo_min( ccdl::gopt::steppt const & a,
				ccdl::gopt::steppt const & b )
{
  return a.rfo < b.rfo;
}




ccdl::gopt::StepInfo::StepInfo ()
  : de(0.), de_abs(0.), de_pred(0.),
    gc_rms(0.), gc_max( -1.e+30 ),
    dxc_rms(0.), dxc_max( -1.e+30 ), dxc_len(0.)
{}

void ccdl::gopt::StepInfo::CptInfo
( ccdl::gopt::Step & step, 
  ccdl::gopt::Step & prevstep)
{
  de      = 0.; 
  de_abs  = 0.;
  de_pred = 0.;
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

  //std::printf("cpt info\n");
  for ( int a=0; a<nat; ++a )
    {
      double xcnrm = 0.;
      double gcnrm = 0.;
      for ( int k=0; k<3; ++k )
	{
	  int i=k+a*3;
	  //std::printf("%20.10f (%20.10f)",step.x[i],prevstep.x[i]);
	  step.dxc[i] = step.x[i]-prevstep.x[i];
	  step.dgc[i] = step.g[i]-prevstep.g[i];
	  xcnrm += step.dxc[i]*step.dxc[i];
	  gcnrm += step.g[i]*step.g[i];
	};
      //std::printf("\n");
      dxc_len += xcnrm;
      dxc_rms += xcnrm;
      gc_rms  += gcnrm;
      //std::printf("max -> %12.8f\n",std::sqrt( xcnrm ) );
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

      de_pred = ccdl::gopt::PredictEnergyChange
	( step.nq, step.hq.data(), prevstep.gq.data(), step.dxq.data() );
    }
  else
    {
      de_pred = ccdl::gopt::PredictEnergyChange
	( step.n, step.h.data(), prevstep.g.data(), step.dxc.data() );
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

  // std::cerr << "Post transform gradients\n";
  // ccdl::WriteXyz( std::cerr, nat, z.data(), g.data() );

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
	 << FMTE(11,2) << maxstep
	 << FMTE(11,2) << info.dxc_max 
	 << FMTE(11,2) << info.dxc_rms
	 << FMTE(11,2) << info.gc_max
	 << FMTE(11,2) << info.gc_rms
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
      cout <<      "GEOMOPT Predicted dE " << FMTE(14,4) << info.de_pred << "\n";
      cout <<      "GEOMOPT Actual/Pred  " << FMTE(14,4) << info.de/info.de_pred << "\n";

      cout << "GEOMOPT Step Length  "
	   << FMTE(14,4) << info.dxc_max //info.dxc_len
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


  // std::printf("Hold\n");
  // for ( int i=0; i<ncrd; ++i )
  //   {
  //     for ( int j=0; j<ncrd; ++j )
  // 	std::printf("%12.3e",Hold[i+j*ncrd]);
  //     std::printf("\n");
  //   };
  // std::printf("DG\n");
  // for ( int i=0; i<ncrd; ++i )
  //   std::printf("%12.3e",DG[i]);
  // std::printf("\n");

  // std::printf("DX\n");
  // for ( int i=0; i<ncrd; ++i )
  //   std::printf("%12.3e",DX[i]);
  // std::printf("\n");

  // std::printf("Hnew\n");
  // for ( int i=0; i<ncrd; ++i )
  //   {
  //     for ( int j=0; j<ncrd; ++j )
  // 	std::printf("%12.3e",Hnew[i+j*ncrd]);
  //     std::printf("\n");
  //   };


}



void ccdl::gopt::Step::Move( double & mxstep )
{
  maxstep = mxstep;

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
	    
	//ccdl::gopt::CptDeltaX_TrustRadius( N, H, G, maxstep, DX, eigvals.data() );
	//ccdl::gopt::CptDeltaX_RFO( N, H, G, maxstep, DX );
	ccdl::gopt::FindStepLength( N, H, G, x.data(), maxstep, DX, eigvals.data(), dlc );

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


      // {
      // 	std::printf("%i\n",nat);
      // 	std::printf("MOVE\n");
      // 	for ( int i=0; i<nat; ++i )
      // 	  {
      // 	    std::printf("X ");
      // 	    for ( int k=0; k<3; ++k )
      // 	      std::printf("%20.10f",x[k+i*3]/ccdl::AU_PER_ANGSTROM);
      // 	    std::printf("\n");
      // 	  }
      // };

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

