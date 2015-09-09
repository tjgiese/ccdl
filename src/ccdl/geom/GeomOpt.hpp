#ifndef _GeomOpt_hpp_
#define _GeomOpt_hpp_

#include <vector>
#include <tr1/array>
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
      void CptDeltaCrds( ccdl::gopt::Step & step, ccdl::gopt::Step & prevstep );

      double de,de_abs,de_pred,de_ratio;
      double pgc_rms, pgc_max;
      double dxc_rms, dxc_max, dxc_len;
      int pgc_max_atm;
      int dxc_max_atm;
      double pgc_max_atm_dxc_pgc_cosangle;
      double dxc_max_atm_dxc_pgc_cosangle;

      std::tr1::array<double,4> pgc_max_atm_pgc_uvec;
      std::tr1::array<double,4> dxc_max_atm_pgc_uvec;
      std::tr1::array<double,4> pgc_max_atm_dxc_uvec;
      std::tr1::array<double,4> dxc_max_atm_dxc_uvec;
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

      void GrdTransform();
      void GrdAndHesTransform();


      int nat,nc,nq,nmax;

      double maxstep;

      double e,predicted_de;
      std::vector<double> xc,gc,hc,pgc,phc;
      std::vector<double> xq,gq,hq;

      std::vector<double> eigvals;

      std::vector<double> dxc,dgc,dpgc;
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
















namespace ccdl
{
  namespace gopt
  {


    double MinSeparation( int nat, double const * c );
    double MaxDisplacement( int nat, double const * newc, double const * oldc );
    double MaxDisplacement( int nat, double const * dc );

  }
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
      DEL = std::abs(DEL);

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
  ccdl::gopt::CptHessian( nat, xc.data(), opts->delta,
			  hc.data(), fcn );
  GrdAndHesTransform();
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
  double limstep      = opts.limstep;
  if ( opts.type == ccdl::gopt::TS )
    {
      if ( opts.update != ccdl::gopt::PSB and
	   opts.update != ccdl::gopt::PSBMS )
	opts.update = ccdl::gopt::PSBMS;
      opts.calcfc = true;
    }

  int nc = 3*nat;
  int nhc = nc*nc;
  ccdl::gopt::Step step( nat, xc, opts, dlc );
  ccdl::gopt::Step prevstep( nat, xc, opts, dlc );
  step.prevstep = &prevstep;

  bool backtracking = false;
  int  nbacktracks = 0;
  bool CONVERGED = false;
  bool calculated_fc = false;

  
  for ( int iter=0; iter < opts.maxiter; ++iter )
    {

      {
	std::vector<int> z( nat, 6 );
	ccdl::WriteXyz( std::cout, nat, z.data(), step.xc.data() );
      }
      
      //-----------------------------------------------------
      step.e = fcn( nat, step.xc.data(), step.gc.data() );
      //-----------------------------------------------------
      step.CptInfo( prevstep );
      step.Print( cout, iter );
      CONVERGED = step.CheckConvergence();

      if ( CONVERGED )
	{
	  std::copy( step.xc.data(), step.xc.data() + nc, xc );
	  break;
	};




#define FMTF(a,b) std::fixed      << std::setw(a) << std::setprecision(b)
#define FMTE(a,b) std::scientific << std::setw(a) << std::setprecision(b)
#define FMT1(str,a,b) (str) << FMTF(14,8) << (a) << FMTE(14,4) << (b) << ( std::abs((a)) > (b) ? " F\n" : " T\n" )
#define FMT2(str,a)   (str) << FMTE(14,4) << (a) << "\n"
#define FMT3(str,a,b) (str) << FMTE(14,4) << (a) << FMTE(14,4) << (b) << ( std::abs((a)) > (b) ? " F\n" : " T\n" )




      bool redo = false;
      {
	double gogn = ccdl::gopt::ddot( nc, prevstep.pgc.data(), step.pgc.data() );
	double gogo = ccdl::gopt::ddot( nc, prevstep.pgc.data(), prevstep.pgc.data() );
	double gngn = ccdl::gopt::ddot( nc, step.pgc.data(), step.pgc.data() );
	double gndx = ccdl::gopt::ddot( nc, step.pgc.data(), step.dxc.data() );
	double godx = ccdl::gopt::ddot( nc, prevstep.pgc.data(), step.dxc.data() );

	double ugg = 0.;
	if ( gogo > 1.e-8 and gngn > 1.e-8 )
	  ugg = gogn / std::sqrt(gogo*gngn);
	double ogrms = std::sqrt( gogo/nc );
	double ngrms = std::sqrt( gngn/nc );

	double WolfeAlpha = 0.1;
	double WolfeBeta  = 0.5;
	double Wolfe1 = WolfeAlpha * gndx;
	double Wolfe2 = WolfeBeta  * godx;


	if ( ngrms > 1.5 * ogrms and iter and step.info.de > 0. )
	  {
	    cout << "GEOMOPT Reject step  " 
		 << FMTF(14,8) << ngrms
		 << FMTF(14,8) << ogrms << " perform a backtrack" << "\n"; 
	    redo = true;
	    backtracking = true;
	    limstep = std::min(limstep,maxstep*1.5);
	    maxstep /= ( 1.5 * 0.8 + std::min(4.,ngrms/ogrms) * 0.2 );
	  }
	else if ( ngrms > ogrms and iter and step.info.de > 0. )
	  {
	    cout << "GEOMOPT Reject step  " 
		 << FMTF(14,8) << ngrms
		 << FMTF(14,8) << ogrms << " perform a backtrack" << "\n"; 
	    redo = true;
	    backtracking = true;
	    limstep = std::min(limstep,maxstep*1.5);
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
      }

      
      if ( step.info.dxc_max > 1.e-20 )
       	maxstep = std::min( 3 * step.info.dxc_max, maxstep );

      
      if ( ! backtracking )
	{

	  if ( step.info.de_ratio < 0.5 and iter > 0 )
	    maxstep *= 0.5;
	  else if ( step.info.de_ratio < 0.9 and iter > 0 )
	    maxstep *= 0.9;
	  else if ( step.info.de_ratio < 0.95 and iter > 0 ) 
	    maxstep *= 1.00;
	  else if ( step.info.de_ratio < 1.1 and step.info.de_ratio > 0.95 )
	    {
	      maxstep *= 2.0;
	      limstep  = std::min( opts.limstep, 1.25 * limstep );
	    }
	  maxstep = std::max( 0.0001 * step.info.pgc_max, maxstep );
	};
      

      maxstep = std::min( maxstep, limstep );
      if ( redo )
	{
	  step = prevstep;
	}
      else if ( step.info.pgc_rms < 1. )
	{
	  // //////////////////////
	  if ( step.info.dxc_max > 1.e-20 && not backtracking )
	    {
	      maxstep = 2 * step.info.dxc_max;
	      if ( step.info.pgc_rms > 0.1 )maxstep = std::min( maxstep, step.info.pgc_rms / 3. );
	    };
	  // //////////////////////
	  maxstep = std::max( maxstep, step.info.pgc_rms / 100. );
	  maxstep = std::min( maxstep, limstep );
	}
      step.maxstep = maxstep;


      

      bool update_hessian = false;
      if ( ! redo )
	{
	  bool cpt_hessian = opts.calcall;
	  if ( opts.calcfc and not calculated_fc )
	    { 
	      if ( step.info.pgc_max > 50. )
		cout << "GEOMOPT Postponing calcfc because the gradient is so large\n";
	      else
		{
		  cpt_hessian = true;
		  calculated_fc = true;
		}
	    };
	  
	  if ( cpt_hessian ) 
	    {
	      step.CptHessian( fcn );
	    }
	  else
	    {
	      update_hessian = true;
	      step.UpdateHessian();
	    };
	  //std::printf("prevstep = step\n");
	  prevstep = step; step.prevstep = &prevstep;
	}
      else if ( backtracking and not opts.calcall )
	{
	  nbacktracks++;
	  // if ( nbacktracks > 9 and nbacktracks % 10 == 0 )
	  //   {
	  //     cout << "GEOMOPT Too many backtracks. Resetting Hessian to diagonal.\n";
	  //     if ( dlc != NULL )
	  //     	for ( int i=0; i<step.nq; ++i )
	  //     	  for ( int j=0; j < step.nq; ++j )
	  //     	    if ( i != j )
	  //     	      step.hq[j+i*step.nq] = 0.;
	  //     	    else
	  //     	      step.hq[j+i*step.nq] += 1.;
	  //     int n3 = step.nat*3;
	  //     for ( int i=0; i<n3; ++i )
	  //     	for ( int j=0; j < n3; ++j )
	  //     	  if ( i != j )
	  //     	    step.h[j+i*n3] = 0.;
	  //     	  else
	  //     	    step.h[j+i*n3] += 1.;
	  //     prevstep = step; step.prevstep = &prevstep;
	  //   }
	};
      

      step.Move( maxstep );


      // double steplen = 0.;
      // bool posdef_h = step.eigvals[0] > -1.e-8 and iter > 1 and opts.calcfc;

      // while ( true )
      // 	{
      // 	  step.Move( maxstep );


      // 	  {
      // 	    double const * t = step.xc.data();
      // 	    double minsep = 1000000.;
      // 	    for ( int i=1; i<nat; ++i )
      // 	      for ( int j=0; j<i; ++j )
      // 		{
      // 		  double x = t[0+i*3]-t[0+j*3];
      // 		  double y = t[1+i*3]-t[1+j*3];
      // 		  double z = t[2+i*3]-t[2+j*3];
      // 		  double r = std::sqrt( x*x+y*y+z*z );
      // 		  minsep = std::min( minsep, r );
      // 		};
      // 	    if ( minsep < 0.5 * ccdl::AU_PER_ANGSTROM )
      // 	      {
      // 		cout << "GEOMOPT Rescaling step to avoid close contact: " 
      // 		     << FMTF(14,8) << minsep / ccdl::AU_PER_ANGSTROM << " A";
      // 		maxstep /= 2.;
      // 		continue;
      // 	      };
      // 	  }
	 
      // 	  break; //
      // 	  if ( ! update_hessian or ! posdef_h ) break;
      // 	  if ( step.eigvals[0] >= -1.e-8 ) break;

      // 	  posdef_h = false;
      // 	  update_hessian = false;
      // 	  cout << "Recomputing Hessian because it is no longer positive definite\n";

      // 	  //std::printf("cpthessian");
      // 	  step.CptHessian( fcn );
      // 	};


      cout << "GEOMOPT Max dx atom " 
	   << std::setw(4) << step.info.dxc_max_atm+1
	   << " dx/|dx| " 
	   << FMTF(8,4) << step.info.dxc_max_atm_dxc_uvec[0]
	   << FMTF(8,4) << step.info.dxc_max_atm_dxc_uvec[1]
	   << FMTF(8,4) << step.info.dxc_max_atm_dxc_uvec[2]
	   << " |dx| "
	   << FMTF(12,4) << step.info.dxc_max_atm_dxc_uvec[3]
	   << " proposed for the next iter\n";

      cout << "GEOMOPT Max dx atom " 
	   << std::setw(4) << step.info.dxc_max_atm+1
	   << " -g/|g|  " 
	   << FMTF(8,4) << -step.info.dxc_max_atm_pgc_uvec[0]
	   << FMTF(8,4) << -step.info.dxc_max_atm_pgc_uvec[1]
	   << FMTF(8,4) << -step.info.dxc_max_atm_pgc_uvec[2]
	   << " |g|  "
	   << FMTF(12,4) <<  step.info.dxc_max_atm_pgc_uvec[3]
	   << " in the current iter\n";
      
      cout << "GEOMOPT Max dx atom "
	   << std::setw(4) << step.info.dxc_max_atm+1
	   << " steepest descent angle  "
	   << FMTF(8,2) << std::acos(-step.info.dxc_max_atm_dxc_pgc_cosangle) * 180./ccdl::PI
	   << " deg\n";


      cout << "GEOMOPT Max  g atom " 
	   << std::setw(4) << step.info.pgc_max_atm+1
	   << " -g/|g|  " 
	   << FMTF(8,4) << -step.info.pgc_max_atm_pgc_uvec[0]
	   << FMTF(8,4) << -step.info.pgc_max_atm_pgc_uvec[1]
	   << FMTF(8,4) << -step.info.pgc_max_atm_pgc_uvec[2]
	   << " |g|  "
	   << FMTF(12,4) <<  step.info.pgc_max_atm_pgc_uvec[3]
	   << " in the current iter\n";

      cout << "GEOMOPT Max  g atom " 
	   << std::setw(4) << step.info.pgc_max_atm+1
	   << " dx/|dx| " 
	   << FMTF(8,4) << step.info.pgc_max_atm_dxc_uvec[0]
	   << FMTF(8,4) << step.info.pgc_max_atm_dxc_uvec[1]
	   << FMTF(8,4) << step.info.pgc_max_atm_dxc_uvec[2]
	   << " |dx| "
	   << FMTF(12,4) << step.info.pgc_max_atm_dxc_uvec[3]
	   << " proposed for the next iter\n";

      double pgc_angle = std::acos(-step.info.pgc_max_atm_dxc_pgc_cosangle) * 180./ccdl::PI;
      cout << "GEOMOPT Max  g atom "
	   << std::setw(4) << step.info.pgc_max_atm+1
	   << " steepest descent angle  "
	   << FMTF(8,2) << pgc_angle
	   << " deg\n";


      if ( step.info.dxc_max_atm != step.info.pgc_max_atm )
	cout << "GEOMOPT The atom with the largest force is not the one with proposed greatest change\n";

      if ( pgc_angle > 90. )
	cout << "GEOMOPT The atom with the largest force is not moving in a downhill direction\n";

      
      if ( (pgc_angle > 75.) )
	{
	  cout << "GEOMOPT The steepest descent angle of the atom with the largest gradient\n"
	       << "GEOMOPT is excessive. I will use a steepest descent direction instead.\n";
	  step.xc = prevstep.xc;
	  std::vector<double> s( step.pgc );
	  double maxdisp = 0.;
	  for ( int a=0; a<step.nat; ++a )
	    {
	      double disp = std::sqrt( s[0+a*3]*s[0+a*3]+
				       s[1+a*3]*s[1+a*3]+
				       s[2+a*3]*s[2+a*3] );
	      s[0+a*3] *= -1;
	      s[1+a*3] *= -1;
	      s[2+a*3] *= -1;
	      maxdisp = std::max( maxdisp, disp );
	    };
	  if ( step.info.pgc_rms < 10. )
	    for ( int a=0; a<step.nat; ++a )
	      {
		double disp = std::sqrt( s[0+a*3]*s[0+a*3]+
					 s[1+a*3]*s[1+a*3]+
					 s[2+a*3]*s[2+a*3] );
		
		if ( disp < 0.5 * maxdisp )
		  {
		    s[0+a*3] *= 0.5;
		    s[1+a*3] *= 0.5;
		    s[2+a*3] *= 0.5;
		  }
	      };

	  double scl = std::min( 0.5*maxstep, step.info.pgc_max_atm_dxc_uvec[3]) / maxdisp;
	  scl = 0.5*maxstep/maxdisp;
	  for ( int i=0; i<step.nc; ++i ) 
	    s[i] *= scl;
	  for ( int i=0; i<step.nc; ++i )
	    step.xc[i] += s[i];

	  if ( step.dlc != NULL )
	    {
	      std::vector<double> xbak( step.xc );
	      std::vector<double> dq0( step.nq, 0. );
	      step.dlc->DisplaceByDeltaQ( dq0.data(), step.xc.data() );
	      double minsep = ccdl::gopt::MinSeparation( step.nat, step.xc.data() );
	      double maxdx = ccdl::gopt::MaxDisplacement( step.nat, step.xc.data(), xbak.data() );
	      if ( minsep < 0.6 * ccdl::AU_PER_ANGSTROM or maxdx > 3. * ccdl::AU_PER_ANGSTROM )
		{
		  cout << "GEOMOPT Warning: constraints are not being enforced at this step\n";
		  step.xc=xbak;
		}
	    };

	  step.info.CptDeltaCrds( step, prevstep );

	  cout << "GEOMOPT Max  g atom " 
	       << std::setw(4) << step.info.pgc_max_atm+1
	       << " dx/|dx| " 
	       << FMTF(8,4) << step.info.pgc_max_atm_dxc_uvec[0]
	       << FMTF(8,4) << step.info.pgc_max_atm_dxc_uvec[1]
	       << FMTF(8,4) << step.info.pgc_max_atm_dxc_uvec[2]
	       << " |dx| "
	       << FMTF(12,4) << step.info.pgc_max_atm_dxc_uvec[3]
	       << " proposed for the next iter\n";
	  
	  double pgc_angle = std::acos(-step.info.pgc_max_atm_dxc_pgc_cosangle) * 180./ccdl::PI;
	  cout << "GEOMOPT Max  g atom "
	       << std::setw(4) << step.info.pgc_max_atm+1
	       << " steepest descent angle  "
	       << FMTF(8,2) << pgc_angle
	       << " deg\n";

	  // if ( step.dlc != NULL )
	  //   {
	  //     std::vector<double> dq( step.nq, 0. );
	  //     step.dlc->
	  //   }
	}
      

      cout << "GEOMOPT Backtrack count " << nbacktracks << "\n";
      cout << "GEOMOPT Eigv";
      int neig = std::min( 6, std::min(step.nc,step.nq) );
      for ( int i=0; i<neig; ++i )
	cout << FMTE(11,2) << step.eigvals[i];
      cout << "\n";


      if ( opts.type == ccdl::gopt::TS and step.eigvals[0] >= -1.e-4 )
	cout << "WARNING The lowest eigenvalue not suited for a TS search. This won't work.\n";

      std::copy( step.xc.data(), step.xc.data() + nc, xc );

      

    }
  if ( ! CONVERGED and ! ERROR ) ERROR = 1;
  return ERROR;
}

#undef FMTF
#undef FMTE
#undef FMT1
#undef FMT2

#endif

