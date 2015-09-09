#include "GeomOpt.hpp"
#include "../bmath.hpp"
#include <algorithm>


namespace ccdl
{
  namespace gopt
  {


    double MinSeparation( int nat, double const * c );
    double MaxDisplacement( int nat, double const * newc, double const * oldc );
    double MaxDisplacement( int nat, double const * dc );

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
      std::tr1::shared_ptr< std::vector<double> > BGinv;
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
	  {
	    int nq = dlc->GetNumInternalCrds();
	    int nat = dlc->GetNumAtoms();
	    int nat3 = 3 * nat;
	    dx.assign(nat3,0.);
	    BGinv.reset( new std::vector<double>( nat3*nq, 0. ) );
	    dlc->CptQ2CEstimator( X0, BGinv->data() );
	  };
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
	


	/*
	double delta = lambda-evals[0];
	std::fill(dq.data(),dq.data()+n,0.);
	for ( int i=0; i<n; ++i )
	  {
	    double alpha = 0.;

	    // ? hmm why does this  work... doesn't seem like it should
	    
	    //double denom = i == 0 ? delta : delta-evals[i]; 

	    // This one does not work at all
	    //double lambda_i = evals[i]-evals[0] + lambda;
	    //double denom = (lambda_i - evals[i]);

	    // This is the regular one
	    double denom = (lambda - evals[i]);

	    if ( std::abs( denom ) > 1.e-30 )
	      alpha = numer[i]/denom;
	    ccdl::axpy( alpha, n, evecs+i*n, dq.data() );
	  };
	*/



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
	  // {
	  //   nat = dlc->GetNumAtoms();
	  //   t.resize( 3*nat );
	  //   std::copy( x0, x0 + 3*nat, t.data() );
	  //   std::vector<double> ddq(dq);
	  //   dlc->DisplaceByDeltaQ( ddq.data(), t.data() );
	  //   nat = dlc->GetNumAtoms();
	  //   dxlen = 0.;
	  //   dxmax = 0.;
	  //   for ( int i=0; i<nat*3; ++i )
	  //     {
	  // 	dx[i] = t[i] - x0[i];
	  // 	dxlen += dx[i]*dx[i];		
	  //     }
	  //   dxlen = std::sqrt( dxlen );
	  // }
	  {
	    nat = dlc->GetNumAtoms();
	    int nat3 = nat*3;
	    dlc->EstimateQ2C( BGinv->data(), dq.data(), dx.data() );
	    dxlen = 0.;
	    dxmax = 0.;
	    for ( int i=0; i<nat3; ++i )
	      dxlen += dx[i]*dx[i]; 
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

      // find the first non-negligble eigenvalue
      // {
      // 	int istart = 0;
      // 	int ineg = ng;
      // 	for ( int i=0; i<ng; ++i )
      // 	  {
      // 	    if ( evals[i] < -1.e-9 )
      // 	      ineg = i;
      // 	    if ( evals[i] > 1.e-9 )
      // 	      {
      // 		istart = i;
      // 		break;
      // 	      }
      // 	  };
      // 	if ( ineg == ng and evals[istart] > 0.1 and std::abs(evals[istart]-1.) > 1 )
      // 	  {
      // 	    steps[0].reset( ng, H, evals, evecs.data(), g, q0.data(), x0, dlc );
      // 	    steps[0].SetLambda( 0. );
	  
      // 	    if ( dlc == NULL )
      // 	      std::copy( steps[0].dx.data(), steps[0].dx.data() + ng, dx );
      // 	    else
      // 	      std::copy( steps[0].dq.data(), steps[0].dq.data() + ng, dx );
	    
      // 	    return steps[0].dxmax;
      // 	  };
	
      // }


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
      int ctr = 0;
      while ( true ) // first find a lambda that doesn't cause lapack to explode
	{
	  ctr++;
	  if ( ctr > 1000 )
	    {
	      std::cerr << "ccdl::gopt::steppt::FindStepLength "
			<< "infinite loop encountered while finding "
			<< "a reasonable lambda\n"
			<< "Expect your optimization to fail"
			<< std::endl;
	      break;
	    };
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

      //std::printf("dxmax @ evals[0] %5i  %13.4e (%13.4e)\n",ctr,lambda_lo,steps[1].dxmax);

      steps[0] = steps[1];
      if ( steps[0].dxmax > trust_radius )
	{
	  double dxprev = steps[0].dxmax;
	  //for ( int i=0; i<20; ++i )
	  ctr = 0;
	  while ( true )
	    {
	      ctr++;
	      if ( ctr > 1000 )
		{
		  std::cerr << "ccdl::gopt::steppt::FindStepLength "
			    << "infinite loop encountered while finding "
			    << "a lambda less than trust_radius\n"
			    << "Expect your optimization to fail"
			    << std::endl;
		  break;
		};

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

      //std::printf("dxmax lbound %5i  %13.4e (%13.4e)\n",ctr,lambda_lo,steps[1].dxmax);

      if ( steps[1].dxmax > trust_radius )
	{


	  //std::printf("DXMAX TOO BIG %13.4e %13.4e\n",steps[1].dxmax, trust_radius);
	  if ( dlc == NULL )
	    {
	      double s   = trust_radius / steps[0].dxmax;
	      for ( int i=0; i<ng; ++i )
		steps[1].dx[i] *= s;
	      std::copy( steps[1].dx.data(), steps[1].dx.data() + ng, dx );
	    }
	  else
	    {
	      int nat = dlc->GetNumAtoms();
	      int nat3 = nat*3;
	      steps[0] = steps[1];
	      std::vector<double> ddq( steps[1].dq );
	      std::vector<double> xv( nat3 );
	      double shi = 1.;
	      double slo = 1.e-12;
	      double s   = trust_radius / steps[0].dxmax;
	      for ( int iter=0; iter < 100; ++iter )
		{
		  for ( std::size_t k=0; k<ddq.size(); ++k )
		    ddq[k] = s * steps[0].dq[k];
		  steps[1].dq = ddq;
		  double maxdx = 0.;

		  // {
		  //   std::copy( x0, x0+nat3, xv.data() );
		  //   dlc->DisplaceByDeltaQ( ddq.data(), xv.data() );
		  //   for ( int a=0; a<nat; ++a )
		  //     {
		  // 	double x = xv[0+a*3]-x0[0+a*3];
		  // 	double y = xv[1+a*3]-x0[1+a*3];
		  // 	double z = xv[2+a*3]-x0[2+a*3];
		  // 	double r2 = x*x+y*y+z*z;
		  // 	maxdx = std::max( maxdx, std::sqrt(r2) );
		  //     };
		  // }
		  {
		    dlc->EstimateQ2C( steps[1].BGinv->data(), ddq.data(), xv.data() );
		    for ( int a=0; a<nat; ++a )
		      {
			double x = xv[0+a*3];
			double y = xv[1+a*3];
			double z = xv[2+a*3];
			double r2 = x*x+y*y+z*z;
			maxdx = std::max( maxdx, std::sqrt(r2) );
		      };
		  }


		  steps[1].dxmax = maxdx;
		  // std::printf("s %13.4e < %13.4e < %13.4e (%13.4e)\n",
		  // 	      slo,s,shi,maxdx);
		  if ( maxdx < trust_radius )
		    {
		      slo = s;
		      if ( std::abs(maxdx-trust_radius) < 1.e-5 )
			{ 
			  //std::printf("dxmax error %5i  %13.4e\n",iter,maxdx);
			  break;
			};
		    }
		  else if ( maxdx > trust_radius )
		    shi = s;
		  s = (slo+shi)/2.;
		}
	      std::copy( steps[1].dq.data(), steps[1].dq.data() + ng, dx );
	    }
	  //std::printf("DXMAX REDUCED TO %13.4e %13.4e\n",steps[1].dxmax, trust_radius);
	  return 1;
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
		{
		  //std::printf("bisection took %i\n",i);
		  break;
		};
	    }
	};

      //std::printf("trust @           %13.4e (%13.4e)\n",steps[2].lambda, steps[2].dxmax );

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
      ctr=0;
      while ( true )
	{
	  ctr++;
	  if ( ctr > 100 )
	    {
	      std::cerr << "ccdl::gopt::steppt::FindStepLength "
			<< "infinite loop encountered while bracketing "
			<< "the rfo\n"
			<< "Expect your optimization to fail"
			<< std::endl;
	      break;
	    };
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
      //std::printf("rfo ubound took %i\n",ctr);
      // the rfo minimum is between steps[1] and steps[2]
      //std::printf("rfo between       %13.4e %13.4e (%13.4e %13.4e) [%13.4e %13.4e]\n",
      //	  steps[1].lambda, steps[2].lambda,
      //	  steps[1].dxmax, steps[2].dxmax,
      //	  steps[1].rfo,    steps[2].rfo  );

      if ( steps[1].dxmax < trust_radius and steps[2].dxmax < trust_radius )
	{
	  steps[0] = steps[1];
	  steps[1] = steps[2];
	  double lold = steps[0].lambda;
	  ctr=0;
	  while ( true )
	    {
	      ctr++;
	      if ( ctr > 100 )
		{
		  std::cerr << "ccdl::gopt::steppt::FindStepLength "
			    << "infinite loop encountered while bisecting "
			    << "the rfo\n"
			    << "Expect your optimization to fail"
			    << std::endl;
		  break;
		};
	      double lnew = 0.5 * ( steps[0].lambda + steps[1].lambda );
	      steps[2].SetLambda( lnew );
	      std::sort( steps.begin(), steps.end(), &ccdl::gopt::steps_rfo_min );

	      /// screw this =============================
	      if ( steps[1].rfo > steps[0].rfo )
		{
		  steps[1] = steps[0];
		  break;
		}
	      else if ( steps[1].rfo > steps[2].rfo )
		{
		  steps[1] = steps[2];
		  break;
		}
	      else
		break;
	      /// screw this =============================

	      // std::printf("rfo %13.5e %13.5e %13.5e (%13.5e %13.5e %13.5e)\n",
	      // 		  steps[0].lambda,steps[1].lambda,steps[2].lambda,
	      // 		  steps[0].rfo,steps[1].rfo,steps[2].rfo);
	      //if ( steps[2].SetLambda( lnew ) );
	      
	      if ( std::abs( steps[1].rfo-steps[0].rfo ) < 1.e-9 or
		   std::abs( steps[1].lambda-steps[0].lambda ) < 1.e-8 or
		   lnew == lold ) break;
	      lold = lnew;

	    }
	  //std::printf("rfo search took %i\n",ctr);
	  //std::printf("rfo   @           %13.4e\n",steps[0].lambda );



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
	  delta = 0.05 * std::abs(steps[1].lambda);
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


      //std::printf("close @           %13.4e (%13.4e) [%13.4e]\n",steps[0].lambda,steps[0].dxmax,min0/ccdl::AU_PER_ANGSTROM );


      if ( dlc == NULL )
	std::copy( steps[0].dx.data(), steps[0].dx.data() + ng, dx );
      else
	std::copy( steps[0].dq.data(), steps[0].dq.data() + ng, dx );


      return steps[0].dxmax;
    }



  }
}

bool ccdl::gopt::steps_rfo_min( ccdl::gopt::steppt const & a,
				ccdl::gopt::steppt const & b )
{
  return a.rfo < b.rfo;
}














double ccdl::gopt::MinSeparation( int nat, double const * t )
{
  double minsep = 10000000.;
  for ( int i=0; i<nat; ++i )
    {
      for ( int j=0; j<i; ++j )
	{
	  double x = t[0+i*3]-t[0+j*3];
	  double y = t[1+i*3]-t[1+j*3];
	  double z = t[2+i*3]-t[2+j*3];
	  double r = std::sqrt( x*x+y*y+z*z );
	  minsep = std::min( minsep, r );
	};
    };
  return minsep;
}

double ccdl::gopt::MaxDisplacement( int nat, double const * newc, double const * oldc )
{
  double maxdx = -1.;
  for ( int i=0; i<nat; ++i )
    {
      double x = newc[0+i*3]-oldc[0+i*3];
      double y = newc[1+i*3]-oldc[1+i*3];
      double z = newc[2+i*3]-oldc[2+i*3];
      double r = std::sqrt( x*x+y*y+z*z );
      maxdx = std::max( maxdx, r );
    };
  return maxdx;
}

double ccdl::gopt::MaxDisplacement( int nat, double const * dc )
{
  double maxdx = -1.;
  for ( int i=0; i<nat; ++i )
    {
      double x = dc[0+i*3];
      double y = dc[1+i*3];
      double z = dc[2+i*3];
      double r = std::sqrt( x*x+y*y+z*z );
      maxdx = std::max( maxdx, r );
    };
  return maxdx;
}










ccdl::gopt::StepInfo::StepInfo()
  : de(0.), de_abs(0.), de_pred(0.), de_ratio(1.0),
    pgc_rms(0.), pgc_max( -1.e+30 ),
    dxc_rms(0.), dxc_max( -1.e+30 ), dxc_len(0.),
    pgc_max_atm(0), dxc_max_atm(0),
    pgc_max_atm_dxc_pgc_cosangle(1.),
    dxc_max_atm_dxc_pgc_cosangle(1.)
{}

void ccdl::gopt::StepInfo::CptInfo
( ccdl::gopt::Step & step, 
  ccdl::gopt::Step & prevstep)
{
  de      = 0.; 
  de_abs  = 0.;
  de_pred = 0.;
  de_ratio = 1.;
  pgc_rms  = 0.; 
  pgc_max  = -1.e+30;
  dxc_rms = 0.; 
  dxc_max = -1.e+30;
  dxc_len = 0.;
  pgc_max_atm = 0;
  dxc_max_atm = 0;
  std::fill( pgc_max_atm_pgc_uvec.begin(), pgc_max_atm_pgc_uvec.end(), 0. );
  std::fill( dxc_max_atm_pgc_uvec.begin(), dxc_max_atm_pgc_uvec.end(), 0. );
  std::fill( pgc_max_atm_dxc_uvec.begin(), pgc_max_atm_dxc_uvec.end(), 0. );
  std::fill( dxc_max_atm_dxc_uvec.begin(), dxc_max_atm_dxc_uvec.end(), 0. );
  pgc_max_atm_dxc_pgc_cosangle=1.;
  dxc_max_atm_dxc_pgc_cosangle=1.;


  int nat = step.nat;

  de = step.e - prevstep.e;
  de_abs = std::abs(de);

  step.dxc.resize( 3*nat );
  step.dgc.resize( 3*nat );
  step.dpgc.resize( 3*nat );

  step.GrdTransform();

  //std::printf("cpt info\n");
  for ( int a=0; a<nat; ++a )
    {
      double xcnrm = 0.;
      double gcnrm = 0.;
      for ( int k=0; k<3; ++k )
	{
	  int i=k+a*3;
	  //std::printf("%20.10f (%20.10f)",step.x[i],prevstep.x[i]);
	  step.dxc[i]  = step.xc[i]-prevstep.xc[i];
	  step.dgc[i]  = step.gc[i]-prevstep.gc[i];
	  step.dpgc[i] = step.pgc[i]-prevstep.pgc[i];
	  xcnrm += step.dxc[i]*step.dxc[i];
	  gcnrm += step.pgc[i]*step.pgc[i];
	};
      //std::printf("\n");
      dxc_len += xcnrm;
      dxc_rms += xcnrm;
      pgc_rms += gcnrm;
      //std::printf("max -> %12.8f\n",std::sqrt( xcnrm ) );
      xcnrm = std::sqrt( xcnrm );
      if ( xcnrm > dxc_max )
	{
	  dxc_max = xcnrm;
	  dxc_max_atm = a;
	  dxc_max_atm_dxc_uvec[0] = step.dxc[0+a*3] / xcnrm;
	  dxc_max_atm_dxc_uvec[1] = step.dxc[1+a*3] / xcnrm;
	  dxc_max_atm_dxc_uvec[2] = step.dxc[2+a*3] / xcnrm;
	  dxc_max_atm_dxc_uvec[3] = xcnrm;
	  double tmp = 0.;
	  tmp = std::sqrt( step.pgc[0+a*3]*step.pgc[0+a*3] +
			   step.pgc[1+a*3]*step.pgc[1+a*3] +
			   step.pgc[2+a*3]*step.pgc[2+a*3] );
	  dxc_max_atm_pgc_uvec[0] = step.pgc[0+a*3] / tmp;
	  dxc_max_atm_pgc_uvec[1] = step.pgc[1+a*3] / tmp;
	  dxc_max_atm_pgc_uvec[2] = step.pgc[2+a*3] / tmp;
	  dxc_max_atm_pgc_uvec[3] = tmp;
	}
      //dxc_max  = std::max( dxc_max, std::sqrt( xcnrm ) );
      gcnrm = std::sqrt( gcnrm );
      if ( gcnrm > pgc_max )
	{
	  pgc_max = gcnrm;
	  pgc_max_atm = a;
	  pgc_max_atm_pgc_uvec[0] = step.pgc[0+a*3] / gcnrm;
	  pgc_max_atm_pgc_uvec[1] = step.pgc[1+a*3] / gcnrm;
	  pgc_max_atm_pgc_uvec[2] = step.pgc[2+a*3] / gcnrm;
	  pgc_max_atm_pgc_uvec[3] = gcnrm;
	  double tmp = 0.;
	  tmp = std::sqrt( step.dxc[0+a*3]*step.dxc[0+a*3] +
			   step.dxc[1+a*3]*step.dxc[1+a*3] +
			   step.dxc[2+a*3]*step.dxc[2+a*3] );
	  pgc_max_atm_dxc_uvec[0] = step.dxc[0+a*3] / tmp;
	  pgc_max_atm_dxc_uvec[1] = step.dxc[1+a*3] / tmp;
	  pgc_max_atm_dxc_uvec[2] = step.dxc[2+a*3] / tmp;
	  pgc_max_atm_dxc_uvec[3] = tmp;
	};
      //pgc_max  = std::max( pgc_max, std::sqrt( gcnrm ) );

    };
  dxc_rms = std::sqrt( dxc_rms / step.nc );
  dxc_len = std::sqrt( dxc_len );
  pgc_rms = std::sqrt( pgc_rms / step.nc );

  pgc_max_atm_dxc_pgc_cosangle = 
    pgc_max_atm_dxc_uvec[0]*pgc_max_atm_pgc_uvec[0] +
    pgc_max_atm_dxc_uvec[1]*pgc_max_atm_pgc_uvec[1] +
    pgc_max_atm_dxc_uvec[2]*pgc_max_atm_pgc_uvec[2];

  dxc_max_atm_dxc_pgc_cosangle = 
    dxc_max_atm_dxc_uvec[0]*dxc_max_atm_pgc_uvec[0] +
    dxc_max_atm_dxc_uvec[1]*dxc_max_atm_pgc_uvec[1] +
    dxc_max_atm_dxc_uvec[2]*dxc_max_atm_pgc_uvec[2];


  pgc_max_atm_dxc_pgc_cosangle = std::min(1.,std::max(-1.,pgc_max_atm_dxc_pgc_cosangle));
  dxc_max_atm_dxc_pgc_cosangle = std::min(1.,std::max(-1.,dxc_max_atm_dxc_pgc_cosangle));


  if ( step.dlc != NULL )
    {
      step.dxq.resize( step.nq );
      step.dgq.resize( step.nq );
      for ( int i=0; i<step.nq; ++i )
	{
	  step.dgq[i] = step.gq[i] - prevstep.gq[i];
	  step.dxq[i] = step.xq[i] - prevstep.xq[i];
	}
      step.dlc->CptDifference( step.xq.data(), prevstep.xq.data(), step.dxq.data() );

      //de_pred = ccdl::gopt::PredictEnergyChange
      //( step.nq, step.hq.data(), prevstep.gq.data(), step.dxq.data() );
    };
  //  else
    {
      de_pred = ccdl::gopt::PredictEnergyChange
	( step.nc, step.phc.data(), prevstep.pgc.data(), step.dxc.data() );
    }
    if ( std::abs( de_pred ) > 1.e-30 )
      de_ratio = de / de_pred;
    else
      de_ratio = 0.;
}



void ccdl::gopt::StepInfo::CptDeltaCrds
( ccdl::gopt::Step & step, 
  ccdl::gopt::Step & prevstep)
{
  dxc_rms = 0.; 
  dxc_max = -1.e+30;
  dxc_len = 0.;
  dxc_max_atm = 0;
  std::fill( dxc_max_atm_dxc_uvec.begin(), dxc_max_atm_dxc_uvec.end(), 0. );
  std::fill( dxc_max_atm_pgc_uvec.begin(), dxc_max_atm_pgc_uvec.end(), 0. );
  std::fill( pgc_max_atm_dxc_uvec.begin(), pgc_max_atm_dxc_uvec.end(), 0. );

  int nat = step.nat;
  step.dxc.resize( 3*nat );

  for ( int a=0; a<nat; ++a )
    {
      double xcnrm = 0.;
      for ( int k=0; k<3; ++k )
	{
	  int i=k+a*3;
	  step.dxc[i]  = step.xc[i]-prevstep.xc[i];
	  xcnrm += step.dxc[i]*step.dxc[i];
	};
      dxc_len += xcnrm;
      dxc_rms += xcnrm;
      xcnrm = std::sqrt( xcnrm );
      if ( xcnrm > dxc_max )
	{
	  dxc_max = xcnrm;
	  dxc_max_atm = a;
	  dxc_max_atm_dxc_uvec[0] = step.dxc[0+a*3] / xcnrm;
	  dxc_max_atm_dxc_uvec[1] = step.dxc[1+a*3] / xcnrm;
	  dxc_max_atm_dxc_uvec[2] = step.dxc[2+a*3] / xcnrm;
	  dxc_max_atm_dxc_uvec[3] = xcnrm;

	  double tmp = 0.;
	  tmp = std::sqrt( step.pgc[0+a*3]*step.pgc[0+a*3] +
			   step.pgc[1+a*3]*step.pgc[1+a*3] +
			   step.pgc[2+a*3]*step.pgc[2+a*3] );
	  dxc_max_atm_pgc_uvec[0] = step.pgc[0+a*3] / tmp;
	  dxc_max_atm_pgc_uvec[1] = step.pgc[1+a*3] / tmp;
	  dxc_max_atm_pgc_uvec[2] = step.pgc[2+a*3] / tmp;
	  dxc_max_atm_pgc_uvec[3] = tmp;
	}
    };

  {
    int a = pgc_max_atm;
    double tmp = 0.;
    tmp = std::sqrt( step.dxc[0+a*3]*step.dxc[0+a*3] +
		     step.dxc[1+a*3]*step.dxc[1+a*3] +
		     step.dxc[2+a*3]*step.dxc[2+a*3] );
    pgc_max_atm_dxc_uvec[0] = step.dxc[0+a*3] / tmp;
    pgc_max_atm_dxc_uvec[1] = step.dxc[1+a*3] / tmp;
    pgc_max_atm_dxc_uvec[2] = step.dxc[2+a*3] / tmp;
    pgc_max_atm_dxc_uvec[3] = tmp;
  }

  dxc_rms = std::sqrt( dxc_rms / step.nc );
  dxc_len = std::sqrt( dxc_len );

  pgc_max_atm_dxc_pgc_cosangle = 
    pgc_max_atm_dxc_uvec[0]*pgc_max_atm_pgc_uvec[0] +
    pgc_max_atm_dxc_uvec[1]*pgc_max_atm_pgc_uvec[1] +
    pgc_max_atm_dxc_uvec[2]*pgc_max_atm_pgc_uvec[2];

  dxc_max_atm_dxc_pgc_cosangle = 
    dxc_max_atm_dxc_uvec[0]*dxc_max_atm_pgc_uvec[0] +
    dxc_max_atm_dxc_uvec[1]*dxc_max_atm_pgc_uvec[1] +
    dxc_max_atm_dxc_uvec[2]*dxc_max_atm_pgc_uvec[2];

  pgc_max_atm_dxc_pgc_cosangle = std::min(1.,std::max(-1.,pgc_max_atm_dxc_pgc_cosangle));
  dxc_max_atm_dxc_pgc_cosangle = std::min(1.,std::max(-1.,dxc_max_atm_dxc_pgc_cosangle));


  if ( step.dlc != NULL )
    {
      step.dxq.resize( step.nq );
      for ( int i=0; i<step.nq; ++i )
	{
	  step.dxq[i] = step.xq[i] - prevstep.xq[i];
	}
      step.dlc->CptDifference( step.xq.data(), prevstep.xq.data(), step.dxq.data() );
    };
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
    nc( nat*3 ), 
    nq( 0 ),
    nmax( nat*3 ),
    maxstep( options.maxstep ),
    e(0.), 
    predicted_de(0.),
    xc(crd,crd+nc),
    gc(nc,0.), 
    hc(nc*nc,0.),
    pgc(nc,0.),
    phc(nc*nc,0.),
    opts( &options ),
    dlc( dlcobj ),
    prevstep( NULL )
{
  if ( dlc != NULL )
    {
      nq = dlc->GetNumInternalCrds();
      nmax = std::max(nc,nq);
      xq.assign(nq,0.);
      gq.assign(nq,0.);
      hq.assign(nq*nq,0.);
      dxq.assign(nq,0.);
      dgq.assign(nq,0.);
      dlc->DisplaceByDeltaQ( dxq.data(), xc.data(), 1.e-7 );
      dlc->CptInternalCrds( xc.data(), xq.data() );
      for ( int i=0; i<nq; ++i )
	hq[i+i*nq] = 1.;
    }
  eigvals.assign( nmax, 0. );
  for ( int i=0; i<nc; ++i )
    hc[i+i*nc] = 1.;
}


void ccdl::gopt::Step::GrdTransform()
{
  pgc = gc;
  if ( dlc != NULL )
    dlc->GrdTransform( xc.data(), pgc.data(), 
		       gq.data(), true );

  // std::cerr << "Pre transform gradients\n";
  // std::vector<int> z(nat,1);
  // ccdl::WriteXyz( std::cerr, nat, z.data(), gc.data() );
  
  // std::cerr << "Post transform gradients\n";
  // ccdl::WriteXyz( std::cerr, nat, z.data(), pgc.data() );

}


void ccdl::gopt::Step::GrdAndHesTransform()
{
  pgc = gc;
  phc = hc;
  if ( dlc != NULL ) 
    dlc->GrdAndHesTransform( xc.data(), pgc.data(), phc.data(), 
			     gq.data(), hq.data(), true );
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
	 << FMTE(11,2) << info.pgc_max
	 << FMTE(11,2) << info.pgc_rms
	 << "\n";
  
  if ( dlc != NULL )
    dlc->PrintReport( cout, xc.data() );
  
  cout << "GEOMOPT Energy " << FMTF(20,8) << e << "\n";
  cout << FMT1("GEOMOPT MAX Force    ",info.pgc_max,opts->gmax_tol);
  cout << FMT1("GEOMOPT RMS Force    ",info.pgc_rms,opts->grms_tol);
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
  CONVERGED = ( info.pgc_max <= opts->gmax_tol and
		info.pgc_rms <= opts->grms_tol and
		info.dxc_max <= opts->xmax_tol and
		info.dxc_rms <= opts->xrms_tol and
		info.de_abs  <= opts->ener_tol );

  if ( ! CONVERGED )
    CONVERGED = ( ( info.pgc_max < 0.05 * opts->gmax_tol and 
		    info.de_abs  < opts->ener_tol ) or
		  ( info.pgc_max < opts->gmax_tol and 
		    info.de_abs < 1.e-12 ) );

  if ( ! CONVERGED )
    CONVERGED = ( info.pgc_max < 0.01 * opts->gmax_tol );

  return CONVERGED;
}




void ccdl::gopt::Step::UpdateHessian()
{
  int ncrd = nc;

  double * Hold = prevstep->hc.data();
  double * Hnew = hc.data();
  double * DG = dgc.data();
  double * DX = dxc.data();

  // if ( dlc != NULL )
  //   {
  //     ncrd = nq;
  //     Hold = prevstep->hq.data();
  //     Hnew = hq.data();
  //     DG = dgq.data();
  //     DX = dxq.data();
  //   };

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

  GrdAndHesTransform();

  //if ( dlc != NULL ) dlc->HesBackTransform( x.data(), gq.data(), hq.data(), h.data() );


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

  std::vector<double> xc0( xc );

  int N = nc;
  double * H = hc.data();
  double * G = gc.data();
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
      //ccdl::gopt::CptDeltaX_TrustRadius( N, H, G, maxstep, DX, eigvals.data() );
      ccdl::gopt::FindStepLength( N, H, G, xc.data(), maxstep, DX, eigvals.data(), dlc );
      break;
    case ccdl::gopt::TS: // need to add following options
      ccdl::gopt::CptDeltaX_EigenFollow
	( N, H, G, maxstep, DX, eigvals.data() );
      break;
    default:
      break;
    }

  xc = prevstep->xc;
  xq = prevstep->xq;
  if ( dlc == NULL )
    {
      for ( int i=0; i<nc; ++i ) xc[i] += dxc[i];
    }
  else
    {
      double disp_tol = std::max( 1.e-9, std::min( 1.e-6, info.dxc_max/10. ) );
      dlc->DisplaceByDeltaQ( dxq.data(), xc.data(), disp_tol );

      double minsep = ccdl::gopt::MinSeparation(nat,xc.data());
      double maxdx  = ccdl::gopt::MaxDisplacement(nat,xc.data(),xc0.data());

      bool bad_step =  minsep < 0.6 * ccdl::AU_PER_ANGSTROM or maxdx > 3. * ccdl::AU_PER_ANGSTROM;

      if ( bad_step )
	{
	  *(opts->ostr) << "GEOMOPT DLC step failed; performing step in cartesians\n";
	  xc = xc0;
	  N = nc;
	  H = phc.data();
	  G = pgc.data();
	  DX = dxc.data();
	  switch ( opts->type )
	    {
	    case ccdl::gopt::MIN:
	      //ccdl::gopt::CptDeltaX_TrustRadius( N, H, G, maxstep, DX, eigvals.data() );
	      ccdl::gopt::FindStepLength( N, H, G, xc.data(), maxstep, DX, eigvals.data(), NULL );
	      break;
	    case ccdl::gopt::TS: // need to add following options
	      ccdl::gopt::CptDeltaX_EigenFollow
		( N, H, G, maxstep, DX, eigvals.data() );
	      break;
	    default:
	      break;
	    }

	  maxdx = ccdl::gopt::MaxDisplacement(nat,dxc.data());

	  double nrm = maxstep / maxdx;
	  for ( int i=0; i<nc; ++i ) 
	    {
	      dxc[i] *= nrm;
	      xc[i] += dxc[i];
	    };
	  // now re-enforce constraints
	  dxq.assign(nq,0.);
	  dlc->DisplaceByDeltaQ( dxq.data(), xc.data() );
	};

      for ( int i=0; i<nc; ++i )
	dxc[i] = xc[i]-xc0[i];

      double nrm = 1.;

      maxdx = ccdl::gopt::MaxDisplacement( nat, xc.data(), xc0.data() );
      if ( maxdx > maxstep )
	{
	  nrm *= maxstep / maxdx;
	  for ( int i=0; i<nc; ++i )
	    {
	      dxc[i] *= nrm;
	      xc[i] = xc0[i] + nrm * dxc[i];
	    };
	  dxq.assign(nq,0.);
	  dlc->DisplaceByDeltaQ( dxq.data(), xc.data() );
	}

      bad_step =  minsep < 0.6 * ccdl::AU_PER_ANGSTROM or maxdx > 3. * ccdl::AU_PER_ANGSTROM;

      if ( bad_step )
	{
	  *(opts->ostr) << "GEOMOPT Warning: constraints are not being enforced at this step\n";
	  maxdx = ccdl::gopt::MaxDisplacement( nat, xc.data(), xc0.data() );
	  nrm *= maxstep / maxdx;
	  for ( int i=0; i<nc; ++i )
	    {
	      dxc[i] *= nrm;
	      xc[i] = xc0[i] + nrm * dxc[i];
	    };
	}
    };

  info.CptDeltaCrds( *this, *prevstep );

}






ccdl::OptOptions::OptOptions()
  : type( ccdl::gopt::MIN ),
    update( ccdl::gopt::BFGS ),
    delta( 1.e-4 ),
    calcfc( false ),
    calcall( false ),
    maxiter( 100 ),
    maxstep( 1.0 ),
    limstep( 1.0 ),
    varmaxstep( true ),
    eigvec( -1 ),
    ener_tol( 1.e-7 ),
    gmax_tol( 1.e-4 ),
    xmax_tol( 1.e-3 ),
    grms_tol( 5.e-5 ),
    xrms_tol( 5.e-4 ),
    ostr( &std::cout )
{}

