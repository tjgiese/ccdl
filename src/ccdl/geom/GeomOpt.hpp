#ifndef _GeomOpt_hpp_
#define _GeomOpt_hpp_

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "../constants.hpp"
#include "OptPrimitives.hpp"
#include "RedundantIC.hpp"




namespace ccdl
{

  struct OptOptions; // fwd decl
  namespace gopt
  {

    struct Step;  // fwd decl
    struct StepInfo
    {
      StepInfo();
      void CptInfo( ccdl::gopt::Step & step, ccdl::gopt::Step & prevstep );

      double de,de_abs,de_pred;
      double gc_rms, gc_max;
      double dxc_rms, dxc_max, dxc_len;
    };


    struct Step
    {
      
      //Step( int nat, double const * crd, ccdl::OptOptions & opts );
      Step( int nat, double const * crd, ccdl::OptOptions & opts, ccdl::RedundantIC * dlc = NULL );


      void CptInfo( ccdl::gopt::Step & prevstep );
      void Print( std::ostream & cout, int const iter );
      bool CheckConvergence();

      template <class T>
      void CptHessian( T fcn );
      void UpdateHessian();
      void Move( double & maxstep );

      void GrdTransform( bool q2c );
      void GrdAndHesTransform( bool q2c );


      int nat,n,nq,nmax;

      double maxstep;

      double e,predicted_de;
      std::vector<double> x,g,h;
      std::vector<double> xq,gq,hq;

      std::vector<double> eigvals;

      std::vector<double> dxc,dgc;
      std::vector<double> dxq,dgq;

      ccdl::gopt::StepInfo info;
      ccdl::OptOptions * opts;
      ccdl::RedundantIC * dlc;
      ccdl::gopt::Step * prevstep;
      
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

  template <class T>
  int DLCOptMin( int nat, double * crd, ccdl::OptOptions opts, T fcn, ccdl::RedundantIC * dlc = NULL );

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




template <class T>
void ccdl::gopt::Step::CptHessian( T fcn )
{
  ccdl::gopt::CptHessian( nat, x.data(), opts->delta,
			  h.data(), fcn );
  GrdAndHesTransform( true );
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

  ccdl::gopt::Step step( nat, xc, opts ), prevstep( nat, xc, opts );
  step.prevstep = &prevstep;

  std::vector<double> DGC( nxc, 0. ), DXC( nxc, 0. ),EIGVALS( nxc, 0. );
  double *__restrict__ dgc = DGC.data();
  double *__restrict__ dxc = DXC.data();
  double *__restrict__ eigvals = EIGVALS.data();
  double de = 0.;
  std::copy( xc, xc+nxc, step.x.data() );

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
	  cout << FMT3("GEOMOPT Actual dE    ",de,opts.ener_tol);
	  cout << "GEOMOPT Step Length  "
	       << FMTE(14,4) << xlen
	       << FMTE(14,4) << maxstep << " max\n";

	  cout << "GEOMOPT SUMMARY "
	       << std::setw(4) << iter
	       << FMTE(16,8) << step.e
	       << FMTE(13,4) << maxstep
	       << FMTE(13,4) << gmax
	       << FMTE(13,4) << grms
	       << "\n";

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

      bool update_hessian = false;
      if ( ! redo )
	{
	  bool cpt_hessian = opts.calcall;
	  if ( iter == 0 and opts.calcfc ) 
	    cpt_hessian = true;
	  if ( cpt_hessian ) 
	    {
	      ccdl::gopt::CptHessian( nat, step.x.data(), opts.delta,
				      step.h.data(), fcn );
	      prevstep = step; step.prevstep = &prevstep;
	    }
	  else
	    {
	      prevstep = step; step.prevstep = &prevstep;
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





















namespace ccdl
{
  template <class T>
  int DLCOptMin_1
  ( int nat, double * xc, 
    ccdl::RedundantIC & dlc,
    ccdl::OptOptions opts,
    T fcn );
}



template <class T>
int ccdl::DLCOptMin_1
( int nat, double * xc, 
  ccdl::RedundantIC & dlc,
  ccdl::OptOptions opts,
  T fcn )
{

  {
    int nxq = dlc.GetNumInternalCrds();
    std::vector<double> dxq(nxq,0.);
    dlc.DisplaceByDeltaQ( dxq.data(), xc, 1.e-8 );
  }


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
  ccdl::gopt::Step step( nat, xc, opts ), prevstep( nat, xc, opts );
  std::vector<double> DGC( nxc, 0. ), DXC( nxc, 0. ),EIGVALS( nxc, 0. );
  double *__restrict__ dgc = DGC.data();
  double *__restrict__ dxc = DXC.data();
  double *__restrict__ eigvals = EIGVALS.data();
  double de = 0.;
  std::copy( xc, xc+nxc, step.x.data() );

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
      dlc.PrintReport( cout, step.x.data() );
      cout << "GEOMOPT Energy " << FMTF(20,8) << step.e << "\n";
      cout << FMT1("GEOMOPT MAX Force    ",gmax,opts.gmax_tol);
      cout << FMT1("GEOMOPT RMS Force    ",grms,opts.grms_tol);
      if ( iter > 0 )
	{
	  cout << FMT1("GEOMOPT MAX Displ    ",xmax,opts.xmax_tol);
	  cout << FMT1("GEOMOPT RMS Displ    ",xrms,opts.xrms_tol);
	  cout << FMT3("GEOMOPT Actual dE    ",de,opts.ener_tol);
	  cout << "GEOMOPT Step Length  "
	       << FMTE(14,4) << xlen
	       << FMTE(14,4) << maxstep << " max\n";

	  cout << "GEOMOPT SUMMARY "
	       << std::setw(4) << iter
	       << FMTE(16,8) << step.e
	       << FMTE(13,4) << maxstep
	       << FMTE(13,4) << gmax
	       << FMTE(13,4) << grms
	       << "\n";

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

      bool update_hessian = false;
      if ( ! redo )
	{
	  bool cpt_hessian = opts.calcall;
	  if ( iter == 0 and opts.calcfc ) 
	    cpt_hessian = true;
	  if ( cpt_hessian ) 
	    {
	      ccdl::gopt::CptHessian( nat, step.x.data(), opts.delta,
				      step.h.data(), fcn );
	      prevstep = step; step.prevstep = &prevstep;
	    }
	  else
	    {
	      prevstep = step; step.prevstep = &prevstep;
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



      int nxq = dlc.GetNumInternalCrds();
      std::vector<double> gq(nxq,0.),hq(nxq*nxq,0.),dxq(nxq,0.);
      { // cart -> dlc
	dlc.GrdAndHesTransform( prevstep.x.data(), prevstep.g.data(), 
				prevstep.h.data(), gq.data(), hq.data(), true );
      }


      while ( true )
	{

	  switch ( opts.type )
	    {
	    case ccdl::gopt::MIN:
	      steplen = ccdl::gopt::CptDeltaX_TrustRadius
		( nxq, hq.data(), gq.data(), 
		  maxstep, dxq.data(), eigvals );
	      break;
	      
	    case ccdl::gopt::TS: // need to add following options
	      steplen = ccdl::gopt::CptDeltaX_EigenFollow
		( nxq, hq.data(), gq.data(), 
		  maxstep, dxq.data(), eigvals );
	      break;
	      
	    default:
	      break;
	    }

	  { // dlc -> cart
	    step.x = prevstep.x;
	    dlc.DisplaceByDeltaQ( dxq.data(), step.x.data(), 1.e-8 );
	    // cout << "GEOMOPT DLC Grd (Hartree-per-unit) and Displ (unit)\n";
	    // for ( int i=0; i< nxq; ++i )
	    //   {
	    // 	cout << "GEOMOPT "; 
	    // 	dlc.PrintCrd(cout,i);
	    // 	dlc.PrintPrettyValue(cout,i,gq.data());
	    // 	dlc.PrintPrettyValue(cout,i,dxq.data());
	    // 	cout << "\n";
	    //   };
	  }


	  {
	    double const * t = step.x.data();
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
	 
	  break; //
	  if ( ! update_hessian or ! posdef_h ) break;
	  if ( eigvals[0] >= -1.e-8 ) break;
	  posdef_h = false;
	  update_hessian = false;
	  cout << "Recomputing Hessian because it is no longer positive definite\n";
	  ccdl::gopt::CptHessian( nat, prevstep.x.data(), opts.delta,
				  prevstep.h.data(), fcn );

	  { // cart -> dlc
	    dlc.GrdAndHesTransform( prevstep.x.data(), prevstep.g.data(), 
				    prevstep.h.data(), gq.data(), hq.data(), true );
	  }
	};

      cout << "GEOMOPT Eigv";
      int neig = std::min( 6, nxc );
      for ( int i=0; i<neig; ++i )
	cout << FMTE(11,2) << eigvals[i];
      cout << "\n";


      if ( opts.type == ccdl::gopt::TS and eigvals[0] >= -1.e-4 )
	cout << "WARNING The lowest eigenvalue not suited for a TS search. This won't work.\n";

      // this needs to change
      //for ( int i=0; i<nxc; ++i ) step.x[i] = prevstep.x[i] + dxc[i];
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









template <class T>
int ccdl::DLCOptMin
( int nat, double * xc, 
  ccdl::OptOptions opts,
  T fcn,
  ccdl::RedundantIC * dlc )
{

  if ( dlc != NULL )
    {
      int nxq = dlc->GetNumInternalCrds();
      std::vector<double> dxq(nxq,0.);
      dlc->DisplaceByDeltaQ( dxq.data(), xc, 1.e-8 );
    }

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
  ccdl::gopt::Step step( nat, xc, opts, dlc );
  ccdl::gopt::Step prevstep( nat, xc, opts, dlc );
  step.prevstep = &prevstep;

  bool backtracking = false;
  int nbacktracks = 0;
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
      step.CptInfo( prevstep );
      step.Print( cout, iter );
      CONVERGED = step.CheckConvergence();
      if ( CONVERGED )
	{
	  std::copy( step.x.data(), step.x.data() + nxc, xc );
	  break;
	};


#define FMTF(a,b) std::fixed      << std::setw(a) << std::setprecision(b)
#define FMTE(a,b) std::scientific << std::setw(a) << std::setprecision(b)
#define FMT1(str,a,b) (str) << FMTF(14,8) << (a) << FMTE(14,4) << (b) << ( std::abs((a)) > (b) ? " F\n" : " T\n" )
#define FMT2(str,a)   (str) << FMTE(14,4) << (a) << "\n"
#define FMT3(str,a,b) (str) << FMTE(14,4) << (a) << FMTE(14,4) << (b) << ( std::abs((a)) > (b) ? " F\n" : " T\n" )




      bool redo = false;
      {
	double gogn = ccdl::gopt::ddot( nxc, prevstep.g.data(), step.g.data() );
	double gogo = ccdl::gopt::ddot( nxc, prevstep.g.data(), prevstep.g.data() );
	double gngn = ccdl::gopt::ddot( nxc, step.g.data(), step.g.data() );
	double gndx = ccdl::gopt::ddot( nxc, step.g.data(), step.dxc.data() );
	double godx = ccdl::gopt::ddot( nxc, prevstep.g.data(), step.dxc.data() );

	double ugg = 0.;
	if ( gogo > 1.e-8 and gngn > 1.e-8 )
	  ugg = gogn / std::sqrt(gogo*gngn);
	double ogrms = std::sqrt( gogo/nxc );
	double ngrms = std::sqrt( gngn/nxc );

	// double WolfeAlpha = 0.1;
	// double WolfeBeta  = 0.5;
	// double Wolfe1 = WolfeAlpha * gndx;
	// double Wolfe2 = WolfeBeta  * godx;

	 

	  if ( ngrms > 1.5 * ogrms and iter and step.info.de > 0. )
	  {
	    cout << "GEOMOPT Reject step  " 
		 << FMTF(14,8) << ngrms
		 << FMTF(14,8) << ogrms << " perform a backtrack" << "\n"; 
	    redo = true;
	    backtracking = true;
	    maxstep /= ( 1.5 * 0.8 + std::min(4.,ngrms/ogrms) * 0.2 );
	  }
	else if ( ngrms > ogrms and iter and step.info.de > 0. )
	  {
	    cout << "GEOMOPT Reject step  " 
		 << FMTF(14,8) << ngrms
		 << FMTF(14,8) << ogrms << " perform a backtrack" << "\n"; 
	    redo = true;
	    backtracking = true;
	    maxstep /= 1.5;
	  }		 
	else if ( ogrms > 1.e-8 and ngrms > 1.e-8 and step.info.de <= 0. and ugg > 0. )
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
	else if ( step.info.de > 0. and iter ) //and de > Wolfe1 ) 
	  {
	    maxstep /= 1.25;
	    backtracking = false;
	  }
	else
	  backtracking = false;
	  //else if ( gndx > Wolfe2 ) maxstep *= 1.001;
      }


      
      if ( step.info.dxc_max > 1.e-20 )
	maxstep = std::min( 3 * step.info.dxc_max, maxstep );



      maxstep = std::min( maxstep, opts.limstep );
      if ( redo )
	{
	  step = prevstep;
	}
      else if ( step.info.gc_rms < 1. )
	{
	  // //////////////////////
	  if ( step.info.dxc_max > 1.e-20 && not backtracking )
	    {
	      maxstep = 2 * step.info.dxc_max;
	      if ( step.info.gc_rms > 0.1 )maxstep = std::min( maxstep, step.info.gc_rms / 3. );
	    };
	  // //////////////////////
	  maxstep = std::max( maxstep, step.info.gc_rms / 100. );
	  maxstep = std::min( maxstep, opts.limstep );
	}
      step.maxstep = maxstep;
      std::printf("maxstep = %13.4e\n",maxstep);


      /*
      if ( backtracking )
	{
	  nbacktracks++;
	  if ( nbacktracks == 4 )
	    {
	      cout << "GEOMOPT Too many backtracks. Resetting Hessian to diagonal.\n";
	      // if ( dlc != NULL )
	      // 	for ( int i=0; i<step.nq; ++i )
	      // 	  for ( int j=0; j < step.nq; ++j )
	      // 	    if ( i != j )
	      // 	      step.hq[j+i*step.nq] = 0.;
	      // 	    else
	      // 	      step.hq[j+i*step.nq] += 1.;
	      // int n3 = step.nat*3;
	      // for ( int i=0; i<n3; ++i )
	      // 	for ( int j=0; j < n3; ++j )
	      // 	  if ( i != j )
	      // 	    step.h[j+i*n3] = 0.;
	      // 	  else
	      // 	    step.h[j+i*n3] += 1.;
	    }
	};
      */

      bool update_hessian = false;
      if ( ! redo )
	{
	  bool cpt_hessian = opts.calcall;
	  if ( iter == 0 and opts.calcfc ) 
	    cpt_hessian = true;

	  if ( cpt_hessian ) 
	    {
	      //std::printf("cpt hessian\n");
	      step.CptHessian( fcn );
	    }
	  else
	    {
	      update_hessian = true;
	      //std::printf("update hessian\n");
	      step.UpdateHessian();
	    };
	  //std::printf("prevstep = step\n");
	  prevstep = step; step.prevstep = &prevstep;
	};

      double steplen = 0.;
      bool posdef_h = step.eigvals[0] > -1.e-8 and iter > 1 and opts.calcfc;

      while ( true )
	{
	  //std::printf("move\n");
	  step.Move( maxstep );

	  {
	    double const * t = step.x.data();
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
	 
	  break; //
	  if ( ! update_hessian or ! posdef_h ) break;
	  if ( step.eigvals[0] >= -1.e-8 ) break;

	  posdef_h = false;
	  update_hessian = false;
	  cout << "Recomputing Hessian because it is no longer positive definite\n";

	  //std::printf("cpthessian");
	  step.CptHessian( fcn );
	};


      cout << "GEOMOPT Eigv";
      int neig = std::min( 6, std::min(step.n,step.nq) );
      for ( int i=0; i<neig; ++i )
	cout << FMTE(11,2) << step.eigvals[i];
      cout << "\n";


      if ( opts.type == ccdl::gopt::TS and step.eigvals[0] >= -1.e-4 )
	cout << "WARNING The lowest eigenvalue not suited for a TS search. This won't work.\n";

      // this needs to change
      //for ( int i=0; i<nxc; ++i ) step.x[i] = prevstep.x[i] + dxc[i];
      std::copy( step.x.data(), step.x.data() + nxc, xc );

      // step.predicted_de = ccdl::gopt::PredictEnergyChange
      // 	( nxc, prevstep.h.data(), prevstep.g.data(), dxc );
      

    }
  if ( ! CONVERGED and ! ERROR ) ERROR = 1;
  return ERROR;
}

#undef FMTF
#undef FMTE
#undef FMT1
#undef FMT2

#endif

