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





template <class T>
int ccdl::CartOptMin
( int nat, double * xc, 
  ccdl::OptOptions opts,
  T fcn )
{
  int ERROR = 0;
  std::ostream & cout = *(opts.ostr);
  double maxstep      = opts.maxstep;

  int nxc = 3*nat;
  int nxh = nxc*nxc;
  ccdl::gopt::Step step( nxc ), prevstep( nxc );
  std::vector<double> DGC( nxc, 0. ), DXC( nxc, 0. );
  double *__restrict__ dgc = DGC.data();
  double *__restrict__ dxc = DXC.data();

  double de = 0.;
  std::copy( xc, xc+nxc, step.x.data() );

  bool CONVERGED = false;
  for ( int iter=0; iter < opts.maxiter; ++iter )
    {
      bool cpt_hessian = opts.calcall;
      if ( iter == 0 and opts.calcfc ) 
	cpt_hessian = true;
      if ( cpt_hessian ) 
	ccdl::gopt::CptHessian( nat, step.x.data(), opts.delta,
				step.h.data(), fcn );

      //-----------------------------------------------------
      step.e = fcn( nat, step.x.data(), step.g.data() );
      //-----------------------------------------------------
      de = step.e - prevstep.e;
      double grms = 0.;
      double xrms = 0.;
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
	  xrms += xnrm;
	  grms += gnrm;
	  xmax = std::max( xmax, std::sqrt( xnrm ) );
	  gmax = std::max( gmax, std::sqrt( gnrm ) );
	};
      xrms = std::sqrt( xrms / nxc );
      grms = std::sqrt( grms / nxc );

      double ade = std::abs(de);
      double dde = std::abs(de-step.predicted_de);

#define FMTF(a,b) std::fixed      << std::setw(a) << std::setprecision(b)
#define FMTE(a,b) std::scientific << std::setw(a) << std::setprecision(b)
#define FMT1(str,a,b) (str) << FMTF(14,8) << (a) << FMTE(11,1) << (b) << ( (a) > (b) ? " F\n" : " T\n" )
#define FMT2(str,a)   (str) << FMTF(14,8) << (a) << "\n"

      cout << FMT1("GEOMOPT MAX Force    ",gmax,opts.gmax_tol);
      cout << FMT1("GEOMOPT RMS Force    ",grms,opts.grms_tol);
      cout << FMT1("GEOMOPT MAX Displ    ",xmax,opts.xmax_tol);
      cout << FMT1("GEOMOPT RMS Displ    ",xrms,opts.xrms_tol);
      cout << FMT2("GEOMOPT Predicted dE ",step.predicted_de);
      cout << FMT1("GEOMOPT Actual dE    ",ade,opts.ener_tol);
      cout << FMT2("GEOMOPT Pred-Act ddE ",dde);

#undef FMTF
#undef FMTE
#undef FMT1
#undef FMT2

      CONVERGED = ( gmax <= opts.gmax_tol and
		    grms <= opts.grms_tol and
		    xmax <= opts.xmax_tol and
		    xrms <= opts.xrms_tol and
		    ade  <= opts.ener_tol );
      if ( CONVERGED )
	{
	  std::copy( step.x.data(), step.x.data() + nxc, xc );
	};

      // std::printf(" x:");
      // for ( int i=0; i<n; ++i ) std::printf("%12.4e",step.x[i]);
      // std::printf(" g:");
      // for ( int i=0; i<n; ++i ) std::printf("%10.2e",step.g[i]);
      // std::printf(" h:");
      // for ( int i=0; i<n*n; ++i ) std::printf("%9.2e",step.h[i]);
      // std::printf("\n");

      if ( ade > 1.e-30 and iter and opts.varmaxstep )
	{
	  double ratio = dde/ade;
	  if ( ratio < 0.3 ) maxstep *= 1.1;
	  else if ( ratio > 0.9 ) maxstep /= 1.1;
	};


      prevstep = step;

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

      double steplen = 0.;
      switch ( opts.type )
	{
	case ccdl::gopt::MIN:
	  steplen = ccdl::gopt::CptDeltaX_TrustRadius
	    ( nxc, prevstep.h.data(), prevstep.g.data(), 
	      maxstep, dxc );
	  break;
	  
	case ccdl::gopt::TS: // need to add following options
	  steplen = ccdl::gopt::CptDeltaX_EigenFollow
	    ( nxc, prevstep.h.data(), prevstep.g.data(), 
	      maxstep, dxc );
	  break;
	  
	default:
	  break;
	}

      
      for ( int i=0; i<nxc; ++i ) step.x[i] = prevstep.x[i] + dxc[i];
      std::copy( step.x.data(), step.x.data() + nxc, xc );

      step.predicted_de = ccdl::gopt::PredictEnergyChange
	( nxc, prevstep.h.data(), prevstep.g.data(), dxc );
      

    }
  if ( ! CONVERGED and ! ERROR ) ERROR = 1;
  return ERROR;
}

#endif

