#ifndef _GeomOpt_hpp_
#define _GeomOpt_hpp_

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "../constants.hpp"
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
    double limstep; // variable step cannot exceed this
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
  if ( DEL > 0. )
    {
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
    }
  else
    {
      std::vector<double> ghi(n,0.);
      fcn( nat, crd, ghi.data() );

      std::vector<double> glo(n,0.);
      for ( int i=0; i<n; ++i )
	{
	  crd[i] -= DEL;
	  fcn( nat, crd, glo.data() );
	  crd[i] += DEL;
	  for ( int j=0; j<n; ++j )
	    h[j+i*n] = (ghi[j]-glo[j])/(DEL);
	};
      for ( int i=0; i<n; ++i )
	for ( int j=0; j<i; ++j )
	  {
	    double avg = 0.5*(h[i+j*n]+h[j+i*n]);
	    h[i+j*n] = avg;
	    h[j+i*n] = avg;
	  };
    }
  return 0;
}


#include <cstdio>
#include "XyzFile.hpp"

extern "C"
{
  double ddot_(const int *n, double const *dx, const int *incx,
	       double const *dy, const int *incy);
}
namespace ccdl
{
  namespace gopt
  {
    double ddot( int n, double const * A, double const * B )
    {
      int inc=1;
      return ddot_( &n, A, &inc, B, &inc );
    }
  }
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

  //std::vector<double> BTUC( nxc, 0. ); // backtracking direction
  //double *__restrict__ btuc = BTUC.data();
  bool backtracking = false;

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
	double gogn = ccdl::gopt::ddot( nxc, prevstep.g.data(), step.g.data() );
	double gogo = ccdl::gopt::ddot( nxc, prevstep.g.data(), prevstep.g.data() );
	double gngn = ccdl::gopt::ddot( nxc, step.g.data(), step.g.data() );
	double gndx = ccdl::gopt::ddot( nxc, step.g.data(), dxc );
	double godx = ccdl::gopt::ddot( nxc, prevstep.g.data(), dxc );

	double ugg = 0.;
	if ( gogo > 1.e-8 and gngn > 1.e-8 )
	  ugg = gogn / std::sqrt(gogo*gngn);
	double ogrms = std::sqrt( gogo/nxc );
	double ngrms = std::sqrt( gngn/nxc );

	double WolfeAlpha = 0.1;
	double WolfeBeta  = 0.5;
	double Wolfe1 = WolfeAlpha * gndx;
	double Wolfe2 = WolfeBeta  * godx;
	// cout << "GEOMOPT Wolfe dE<g.dx" 
	//      << FMTF(14,8) << de
	//      << FMTF(14,8) << Wolfe1 
	//      << (de<Wolfe1 ? " T" : " F") << "\n"; 
	// cout << "GEOMOPT Wolfe ndx<odx" 
	//      << FMTF(14,8) << gndx
	//      << FMTF(14,8) << Wolfe2
	//      << (gndx<Wolfe2 ? " T" : " F") << "\n";
	 

	// if ( gndx > Wolfe2 and de < 0. )
	//   {
	//     cout << "GEOMOPT Reject step  Wolfe 2 -- increase trust radius\n";
	//     redo = true;
	//     maxstep *= 1.5;
	//   }
	// else 
	  if ( ngrms > 1.5 * ogrms and iter and de > 0. )
	  {
	    cout << "GEOMOPT Reject step  " 
		 << FMTF(14,8) << ngrms
		 << FMTF(14,8) << ogrms << " perform a backtrack" << "\n"; 
	    redo = true;
	    backtracking = true;
	    maxstep /= ( 1.5 * 0.8 + std::min(4.,ngrms/ogrms) * 0.2 );
	  }
	else if ( ngrms > ogrms and iter and de > 0. )
	  {
	    cout << "GEOMOPT Reject step  " 
		 << FMTF(14,8) << ngrms
		 << FMTF(14,8) << ogrms << " perform a backtrack" << "\n"; 
	    redo = true;
	    backtracking = true;
	    maxstep /= 1.5;
	  }		 
	else if ( ogrms > 1.e-8 and ngrms > 1.e-8 and de <= 0. and ugg > 0. )
	  {
	    if ( ! backtracking )
	      {
		maxstep *= ( 1. + 0.3*std::exp(-ogrms) * ugg );
	      }
	      else
	      {
		maxstep /= 1.3; //( 1. + std::exp(-ogrms) * ugg );
		cout << "GEOMOPT Decreasing step to avoid backtrack path\n";
	      };
	    backtracking = false;
	  }
	else if ( de > 0. and iter ) //and de > Wolfe1 ) 
	  {
	    maxstep /= 1.25;
	    backtracking = false;
	  }
	else
	  backtracking = false;
	  //else if ( gndx > Wolfe2 ) maxstep *= 1.001;
      }


      maxstep = std::min( maxstep, opts.limstep );
      //maxstep = std::max( maxstep, 0.001 );

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
	  {
	    std::vector<double> t(nxc);
	    for ( int i=0; i<nxc; ++i ) t[i] = prevstep.x[i] + dxc[i];
	    double minsep = 1000000.;
	    for ( int i=1; i<nat; ++i )
	      for ( int j=0; j<i; ++j )
		{
		  double x = t[0+i*3]-t[0+j*3];
		  double y = t[1+i*3]-t[1+j*3];
		  double z = t[2+i*3]-t[2+j*3];
		  double r = std::sqrt( x*x+y*y+z*z );
		  minsep = std::min( minsep, r );
		};
	    if ( minsep < 0.5 * ccdl::AU_PER_ANGSTROM )
	      {
		cout << "GEOMOPT Rescaling step to avoid close contact: " 
		     << FMTF(14,8) << minsep / ccdl::AU_PER_ANGSTROM << " A\n";
		maxstep /= 2.;
		continue;
	      };
	  }
	 
	  if ( ! update_hessian or ! posdef_h ) break;
	  if ( eigvals[0] >= -1.e-8 ) break;
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

