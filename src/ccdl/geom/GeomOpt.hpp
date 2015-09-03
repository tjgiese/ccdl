#ifndef _GeomOpt_hpp_
#define _GeomOpt_hpp_

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "OptPrimitives.hpp"

namespace ccdl
{
  namespace gopt
  {

    struct Step
    {
      
      Step( int n );
      Step( ccdl::gopt::Step const & rhs );

      int n;
      double e,predicted_de;
      std::vector<double> x,g,h;
      
    };

  }
}


namespace ccdl
{

  namespace gopt
  {
    enum type { MIN, TS };
  }

  namespace gopt
  {
    enum update { BFGS, MS, BFGSMS, PSB, PSBMS };
  }

  struct OptOptions
  {
    OptOptions();

    ccdl::gopt::type type;
    ccdl::gopt::update update;
    double delta;
    bool   calcfc;
    bool   calcall;
    int    maxiter;
    double maxstep;
    bool   varmaxstep;
    int    eigvec;
    double ener_tol;
    double gmax_tol;
    double xmax_tol;
    double grms_tol;
    double xrms_tol;
    std::ostream * ostr;
  };
}


namespace ccdl
{
  namespace gopt
  {
    template <class T>
    int CptHessian( int nat, double * crd, double delta, double * h, T fcn );
  }
}



namespace ccdl
{
  template <class T>
  int CartOptMin( int nat, double * crd, ccdl::OptOptions opts, T fcn );
}






template <class T>
int ccdl::gopt::CptHessian
( int nat, double * crd, double DEL, double * h, T fcn )
{
  int n = 3*nat;
  std::vector<double> ghi(n,0.);
  std::vector<double> glo(n,0.);
  for ( int i=0; i<n; ++i )
    {
      crd[i] += DEL;
      fcn( nat, crd, ghi.data() );
      crd[i] -= 2*DEL;
      fcn( nat, crd, glo.data() );
      crd[i] += DEL;
      for ( int j=0; j<n; ++j )
	h[j+i*n] = (ghi[j]-glo[j])/(2.*DEL);
    };
  for ( int i=0; i<n; ++i )
    for ( int j=0; j<i; ++j )
      {
	double avg = 0.5*(h[i+j*n]+h[j+i*n]);
	h[i+j*n] = avg;
	h[j+i*n] = avg;
      };
  return 0;
}



#include "XyzFile.hpp"

template <class T>
int ccdl::CartOptMin
( int nat, double * xc, 
  ccdl::OptOptions opts,
  T fcn )
{
  int ERROR = 0;
  std::ostream & cout = *(opts.ostr);
  double maxstep      = opts.maxstep;

  if ( opts.type == ccdl::gopt::TS )
    {
      if ( opts.update != ccdl::gopt::PSB and
	   opts.update != ccdl::gopt::PSBMS )
	opts.update = ccdl::gopt::PSBMS;
      opts.calcfc = true;
    }

  int nxc = 3*nat;
  int nxh = nxc*nxc;
  ccdl::gopt::Step step( nxc ), prevstep( nxc );
  std::vector<double> DGC( nxc, 0. ), DXC( nxc, 0. ),EIGVALS( nxc, 0. );
  double *__restrict__ dgc = DGC.data();
  double *__restrict__ dxc = DXC.data();
  double *__restrict__ eigvals = EIGVALS.data();
  double de = 0.;
  std::copy( xc, xc+nxc, step.x.data() );


  bool CONVERGED = false;
  for ( int iter=0; iter < opts.maxiter; ++iter )
    {


      {
	std::vector<int> z( nat, 6 );
	ccdl::WriteXyz( std::cout, nat, z.data(), step.x.data() );
      }
      
      //-----------------------------------------------------
      step.e = fcn( nat, step.x.data(), step.g.data() );
      //-----------------------------------------------------
      de = step.e - prevstep.e;
      double grms = 0.;
      double xrms = 0.;
      double xlen = 0.;
      double gmax = -1.e+30;
      double xmax = -1.e+30;
      for ( int a=0; a<nat; ++a )
	{
	  double gnrm = 0.;
	  double xnrm = 0.;
	  for ( int k=0; k<3; ++k )
	    {
	      int i=k+a*3;
	      dxc[i] = step.x[i]-prevstep.x[i];
	      dgc[i] = step.g[i]-prevstep.g[i];
	      xnrm += dxc[i]*dxc[i];
	      gnrm += step.g[i]*step.g[i];
	    };
	  xlen += xnrm;
	  xrms += xnrm;
	  grms += gnrm;
	  xmax = std::max( xmax, std::sqrt( xnrm ) );
	  gmax = std::max( gmax, std::sqrt( gnrm ) );
	};
      xrms = std::sqrt( xrms / nxc );
      grms = std::sqrt( grms / nxc );
      xlen = std::sqrt( xlen );

      double ade = std::abs(de);
      double dde = std::abs(de-step.predicted_de);

#define FMTF(a,b) std::fixed      << std::setw(a) << std::setprecision(b)
#define FMTE(a,b) std::scientific << std::setw(a) << std::setprecision(b)
#define FMT1(str,a,b) (str) << FMTF(14,8) << (a) << FMTE(14,4) << (b) << ( std::abs((a)) > (b) ? " F\n" : " T\n" )
#define FMT2(str,a)   (str) << FMTE(14,4) << (a) << "\n"
#define FMT3(str,a,b) (str) << FMTE(14,4) << (a) << FMTE(14,4) << (b) << ( std::abs((a)) > (b) ? " F\n" : " T\n" )

      cout << "\n\n";
      cout << "GEOMOPT ITER " << std::setw(4) << iter << "\n";
      cout << "GEOMOPT Energy " << FMTF(20,8) << step.e << "\n";
      cout << FMT1("GEOMOPT MAX Force    ",gmax,opts.gmax_tol);
      cout << FMT1("GEOMOPT RMS Force    ",grms,opts.grms_tol);
      if ( iter > 0 )
	{
	  cout << FMT1("GEOMOPT MAX Displ    ",xmax,opts.xmax_tol);
	  cout << FMT1("GEOMOPT RMS Displ    ",xrms,opts.xrms_tol);
	  //cout << FMT2("GEOMOPT Predicted dE ",step.predicted_de);
	  cout << FMT3("GEOMOPT Actual dE    ",de,opts.ener_tol);
	  // if ( std::abs( step.predicted_de ) > 1.e-8 )
	  //   cout << FMT2("GEOMOPT Act/Pred  dE ",de/step.predicted_de);
	  // else
	  //   cout << FMT2("GEOMOPT Act/Pred  dE ",1.234567890123456789);

	  cout << "GEOMOPT Step Length  "
	       << FMTE(14,4) << xlen
	       << FMTE(14,4) << maxstep << " max\n";
	};

      CONVERGED = ( gmax <= opts.gmax_tol and
		    grms <= opts.grms_tol and
		    xmax <= opts.xmax_tol and
		    xrms <= opts.xrms_tol and
		    ade  <= opts.ener_tol );
      if ( ! CONVERGED )
	CONVERGED = ( ( gmax < 0.05 * opts.gmax_tol and 
			ade  < opts.ener_tol ) or
		      ( gmax < opts.gmax_tol and 
			ade < 1.e-12 ) );
      if ( ! CONVERGED )
	CONVERGED = ( gmax < 0.01 * opts.gmax_tol );
      if ( CONVERGED )
	{
	  std::copy( step.x.data(), step.x.data() + nxc, xc );
	  break;
	};

      // std::printf(" x:");
      // for ( int i=0; i<n; ++i ) std::printf("%12.4e",step.x[i]);
      // std::printf(" g:");
      // for ( int i=0; i<n; ++i ) std::printf("%10.2e",step.g[i]);
      // std::printf(" h:");
      // for ( int i=0; i<n*n; ++i ) std::printf("%9.2e",step.h[i]);
      // std::printf("\n");

      /*
      if ( ade > 1.e-30 and iter and opts.varmaxstep )
	{
	  //double ratio = dde/ade;
	  if ( std::abs( step.predicted_de ) > 1.e-8 )
	    {
	      double ratio = de / step.predicted_de;
	      if ( ratio > 0.75 and xlen < 1. )
	        maxstep *= 1.25;
	      else if ( ratio < 0. and de > 0. and xlen > opts.xmax_tol )
		maxstep /= 1.3333333333333333333333333333;
	    }
	  // double rho = de / step.predicted_de;
	  // if ( rho > 0.75 and 5./(4. * xlen) > maxstep )
	  // //if ( rho > 0.75 and 5./4. * xlen > maxstep )
	  //   maxstep *= 2.;
	  // else if ( rho < 0.25 )
	  //   maxstep = 1./(4. * xlen);
	  // //maxstep = 1./4. * xlen;
	};
      */

      bool redo = false;
      {
	std::vector<double> go(prevstep.g), gn(step.g);
	double onrm = 0., nnrm = 0., ond = 0.;
	for ( int i=0; i<nxc; ++i )
	  {
	    onrm += go[i]*go[i];
	    nnrm += gn[i]*gn[i];
	    ond  += go[i]*gn[i];
	  };
	if ( onrm > 1.e-8 and nnrm > 1.e-8 )
	  ond = ond / std::sqrt( onrm*nnrm );
	nnrm = std::sqrt(nnrm/nxc);
	onrm = std::sqrt(onrm/nxc);
	if ( nnrm > 1.5 * onrm and iter and de > 0. )
	  {
	    cout << "GEOMOPT Reject step  " 
		 << FMTF(14,8) << nnrm 
		 << FMTF(14,8) << onrm << " treat as line search" << "\n"; 
	    redo = true;
	    maxstep /= ( 1.5 * 0.8 + std::min(4.,nnrm/onrm) * 0.2 );
	  }
	else if ( nnrm > onrm and iter and de > 0. )
	  {
	    cout << "GEOMOPT Reject step  " 
		 << FMTF(14,8) << nnrm 
		 << FMTF(14,8) << onrm << " treat as line search" << "\n"; 
	    redo = true;
	    maxstep /= 1.5;
	  }		 
	else if ( onrm > 1.e-8 and nnrm > 1.e-8 and de <= 0. and ond > 0. )
	  maxstep *= ( 1. + std::exp(-onrm) * ond );
	//else if ( de < 0. and iter ) maxstep *= 1.05;
	else if ( de > 0. and iter ) 
	  maxstep /= 1.25;
      }

      bool update_hessian = false;
      if ( ! redo )
	{
	  bool cpt_hessian = opts.calcall;
	  if ( iter == 0 and opts.calcfc ) 
	    cpt_hessian = true;
	  //if ( iter > 0 and de > 0. ) cpt_hessian = true;
	  if ( cpt_hessian ) 
	    {
	      ccdl::gopt::CptHessian( nat, step.x.data(), opts.delta,
				      step.h.data(), fcn );
	      prevstep = step;
	    }
	  else
	    {
	      prevstep = step;
	      update_hessian = true;

	      switch ( opts.update )
		{
		case ccdl::gopt::BFGS:
		  // MIN SEARCH
		  ccdl::gopt::UpdateHessian_BFGS
		    ( nxc, prevstep.h.data(), dgc, dxc, step.h.data() );
		  break;
		  
		case ccdl::gopt::MS:
		  // MIN SEARCH
		  ccdl::gopt::UpdateHessian_MS
		    ( nxc, prevstep.h.data(), dgc, dxc, step.h.data() );
		  break;
		  
		case ccdl::gopt::BFGSMS:
		  // MIN SEARCH
		  ccdl::gopt::UpdateHessian_Bofill_BFGS_MS
		    ( nxc, prevstep.h.data(), dgc, dxc, step.h.data() );
		  break;
		  
		case ccdl::gopt::PSB:
		  // TS SEARCH
		  ccdl::gopt::UpdateHessian_PSB
		    ( nxc, prevstep.h.data(), dgc, dxc, step.h.data() );
		  break;
		  
		case ccdl::gopt::PSBMS:
		  // TS SEARCH
		  ccdl::gopt::UpdateHessian_Bofill_PSB_MS
		    ( nxc, prevstep.h.data(), dgc, dxc, step.h.data() );
		  break;
		  
		default:
		  break;
		}


	      // if ( posdef_h and eigvals[0] < -1.e-8 )
	      // 	{
	      // 	  redo = true;
	      // 	  step = prevstep;
	      // 	  maxstep /= 2.;
	      // 	};

	    };
	};

      double steplen = 0.;
      bool posdef_h = eigvals[0] > -1.e-8 and iter > 1 and opts.calcfc;
      while ( true )
	{
	  switch ( opts.type )
	    {
	    case ccdl::gopt::MIN:
	      steplen = ccdl::gopt::CptDeltaX_TrustRadius
		( nxc, prevstep.h.data(), prevstep.g.data(), 
		  maxstep, dxc, eigvals );
	      break;
	      
	    case ccdl::gopt::TS: // need to add following options
	      steplen = ccdl::gopt::CptDeltaX_EigenFollow
		( nxc, prevstep.h.data(), prevstep.g.data(), 
		  maxstep, dxc, eigvals );
	      break;
	      
	    default:
	      break;
	    }
	  if ( ! update_hessian or ! posdef_h ) break;
	  if ( eigvals[0] >= -1.e-8 ) break;
	  //cout << "Scaling step because Hessian is no longer positive\n";// << FMTE(12,3) << maxstep * 0.25 << "\n";
	  //maxstep *= 0.1;
	  posdef_h = false;
	  update_hessian = false;
	  cout << "Recomputing Hessian because it is no longer positive definite\n";
	  ccdl::gopt::CptHessian( nat, prevstep.x.data(), opts.delta,
				  prevstep.h.data(), fcn );
	};

      cout << "GEOMOPT Eigv";
      int neig = std::min( 6, nxc );
      for ( int i=0; i<neig; ++i )
	cout << FMTE(11,2) << eigvals[i];
      cout << "\n";


      if ( opts.type == ccdl::gopt::TS and eigvals[0] >= -1.e-4 )
	cout << "WARNING The lowest eigenvalue not suited for a TS search. This won't work.\n";

      for ( int i=0; i<nxc; ++i ) step.x[i] = prevstep.x[i] + dxc[i];
      std::copy( step.x.data(), step.x.data() + nxc, xc );

      step.predicted_de = ccdl::gopt::PredictEnergyChange
	( nxc, prevstep.h.data(), prevstep.g.data(), dxc );
      

    }
  if ( ! CONVERGED and ! ERROR ) ERROR = 1;
  return ERROR;
}

#undef FMTF
#undef FMTE
#undef FMT1
#undef FMT2

#endif

