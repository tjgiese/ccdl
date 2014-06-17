#include "lda.hpp"
#include <cmath>
#include <algorithm>
#include "../constants.hpp"


void uks_xf_lda_spw92
( double const * pa, double const * pb,
  double * F )
{

  for ( int i=0; i<11; ++i ) F[i] = 0.;

  double const mpa = std::max( pa[0], 0. );
  double const mpb = std::max( pb[0], 0. );

  if ( mpa + mpb < 1.e-100 ) return;

  /*
  double const C13    = 1./3.;
  double const C23    = 2./3.;
  double const C43    = 4./3.;
  double const C213   = std::pow(2.,C13);
  double const Cx     = 0.75 * std::pow(3./ccdl::PI,C13);
  double const ALPHA = C23;

  double const Ea = -Cx * 0.75 * ALPHA * std::pow(2.*mpa,C43);
  double const Eb = -Cx * 0.75 * ALPHA * std::pow(2.*mpb,C43);

  F[0] = Ea+Eb;
  F[1] = - C43 * Cx * 1.5 * ALPHA * std::pow(2.*mpa,C13);
  F[2] = - C43 * Cx * 1.5 * ALPHA * std::pow(2.*mpb,C13);
  */

  //double const a = -Cx * 0.75 * ALPHA * std::pow(2.,C43);
			       
  double const a = -9.305257363491000e-01;
  //double const b = -C43 * Cx * 1.5 * ALPHA * std::pow(2.,C13);
  double const b = -1.240700981798800e+00;
  double const C13 = 3.333333333333333e-01;
  double const pa13 = std::pow(mpa,C13);
  double const pb13 = std::pow(mpb,C13);
  F[0] = a * ( pa13*mpa + pb13*mpb );
  F[1] = b * pa13;
  F[2] = b * pb13;
}


void uks_cf_lda_spw92
( double const * pa, double const * pb,
  double * F )
{
  for ( int i=0; i<11; ++i ) F[i] = 0.;
  double const p = pa[0]+pb[0];
  if ( p < 1.e-100 ) return;


  // PARAMAGENTIC 
  double const A0  = 0.0621814 / 2.;
  double const A00 = 0.21370;
  double const B01 = 7.5957;
  double const B02 = 3.5876;
  double const B03 = 1.6382;
  double const B04 = 0.49294;

  // FEROMAGNETIC 
  double const A1  = A0 / 2.;
  double const A10 =  0.20548;
  double const B11 = 14.1189;
  double const B12 =  6.1977;
  double const B13 =  3.3662;
  double const B14 =  0.62517;

  // PADE APPROX ALPHA
  double const Aa  =  1./(3.*ccdl::PI*ccdl::PI) / 2.;
  double const Aa0 =  0.11125;
  double const Ba1 = 10.357;
  double const Ba2 =  3.6231;
  double const Ba3 =  0.88026;
  double const Ba4 =  0.49671;


  double const C13 = 1./3.;
  double const C43 = 4./3.;
  double const C213 = std::pow(2.,C13);
  double const CDENOM = 2. * (C213-1.);
  double const fpp0 = 4. / (9.*(C213-1.));

  double const z = std::max( std::min(  (pa[0]-pb[0])/p,  1. ), -1. );
  double const rs = std::pow( 0.75/(ccdl::PI*p), C13 );
  double const x = std::sqrt(rs);
  double const f = ( std::pow(1.+z,C43) + std::pow(1.-z,C43) - 2. ) / CDENOM;
  double const z2 = z*z;
  double const z3 = z2*z;
  double const z4 = z3*z;


  double e0,e1,a;
  double D0Drs,D1Drs,DaDrs;
  {
    double Q0,Q1,Q1p,lterm;
    double const r12 = x;
    double const r1  = rs;
    double const r32 = r12*r1;
    double const r2  = r1*r1;

    Q0  = -2.*A0*(1.+A00*r1);
    Q1  =  2.*A0*( B01*r12 + B02*r1 + B03*r32 + B04*r2  );
    Q1p =  A0*( B01/r12 + 2.*B02 + 3.*B03*r12 + 4.*B04*r1 );
    lterm = std::log(1.+1./Q1);

    e0     = Q0 * lterm;
    D0Drs  = -2*A0*A00*lterm - Q0 * Q1p / ( Q1*Q1+Q1 );

    Q0  = -2.*A1*(1.+A10*r1);
    Q1  =  2.*A1*( B11*r12 + B12*r1 + B13*r32 + B14*r2  );
    Q1p =  A1*( B11/r12 + 2.*B12 + 3.*B13*r12 + 4.*B14*r1 );
    lterm = std::log(1.+1./Q1);

    e1     = Q0 * lterm;
    D1Drs  = -2*A1*A10*lterm - Q0 * Q1p / ( Q1*Q1+Q1 );


    Q0  = -2.*Aa*(1.+Aa0*r1);
    Q1  =  2.*Aa*( Ba1*r12 + Ba2*r1 + Ba3*r32 + Ba4*r2  );
    Q1p =  Aa*( Ba1/r12 + 2.*Ba2 + 3.*Ba3*r12 + 4.*Ba4*r1 );
    lterm = std::log(1.+1./Q1);

    a      = - ( Q0 * lterm );
    DaDrs  = - ( -2*Aa*Aa0*lterm - Q0 * Q1p / ( Q1*Q1+Q1 ) );
  };			       

  double const k = a/fpp0*(1.-z4) + (e1-e0)*z4;
  
  F[0] = p * ( e0 + f*k );
  
  double const p2 = p*p;
  double const DzDpa =  2.*pb[0]/p2;
  double const DrsDpa = -rs*C13/p;
  double const DfDz  = C43 * ( std::pow(1.+z,C13) - std::pow(1.-z,C13) ) / CDENOM;
  double const DfDpa = DfDz*DzDpa;
  double const D0Dpa = D0Drs*DrsDpa;
  double const D1Dpa = D1Drs*DrsDpa;
  double const DaDpa = DaDrs*DrsDpa;
  double const PkPpa = DaDpa/fpp0*(1.-z4) + (D1Dpa-D0Dpa)*z4;
  double const DkDz = a/fpp0*(-4.*z3) + (e1-e0)*(4.*z3);
  double const DkDpa = PkPpa + DkDz*DzDpa;

  F[1] = e0 + p*D0Dpa + f*k + p*DfDpa*k + p*f*DkDpa;


  double const DzDpb = -2.*pa[0]/p2;
  double const DrsDpb = DrsDpa;
  double const DfDpb = DfDz*DzDpb;
  double const D0Dpb = D0Drs*DrsDpb;
  double const D1Dpb = D1Drs*DrsDpb;
  double const DaDpb = DaDrs*DrsDpb;
  double const PkPpb = DaDpb/fpp0*(1.-z4) + (D1Dpb-D0Dpb)*z4;
  double const DkDpb = PkPpb + DkDz*DzDpb;

  F[2] = e0 + p*D0Dpb + f*k + p*DfDpb*k + p*f*DkDpb;

}


// F dFdpa dFdpb dFdga dFdgb dFdta dFdtb

void uks_xcf_lda_spw92_potential( double const * pa, double const * pb, double * F )
{
  for ( int i=0; i<11; ++i ) F[i] = 0.;

  double const mpa = std::max( pa[0], 0. );
  double const mpb = std::max( pb[0], 0. );
  double const p = mpa+mpb;

  if ( p < 1.e-100 ) return;

  //double const a = -Cx * 0.75 * ALPHA * std::pow(2.,C43);	       
  double const xfa = -9.305257363491000e-01;
  //double const b = -C43 * Cx * 1.5 * ALPHA * std::pow(2.,C13);
  double const xfb = -1.240700981798800e+00;
  double const C13 = 3.333333333333333e-01;
  double const pa13 = std::pow(mpa,C13);
  double const pb13 = std::pow(mpb,C13);

  // PARAMAGENTIC 
  double const A0  = 0.0621814 / 2.;
  double const A00 = 0.21370;
  double const B01 = 7.5957;
  double const B02 = 3.5876;
  double const B03 = 1.6382;
  double const B04 = 0.49294;

  // FEROMAGNETIC 
  double const A1  = A0 / 2.;
  double const A10 =  0.20548;
  double const B11 = 14.1189;
  double const B12 =  6.1977;
  double const B13 =  3.3662;
  double const B14 =  0.62517;

  // PADE APPROX ALPHA
  double const Aa  =  1./(3.*ccdl::PI*ccdl::PI) / 2.;
  double const Aa0 =  0.11125;
  double const Ba1 = 10.357;
  double const Ba2 =  3.6231;
  double const Ba3 =  0.88026;
  double const Ba4 =  0.49671;


  //double const C13 = 1./3.;
  double const C43 = 4./3.;
  double const C213 = std::pow(2.,C13);
  double const CDENOM = 2. * (C213-1.);
  double const fpp0 = 4. / (9.*(C213-1.));

  double const z = std::max( std::min(  (pa[0]-pb[0])/p,  1. ), -1. );
  double const rs = std::pow( 0.75/(ccdl::PI*p), C13 );
  double const x = std::sqrt(rs);
  double const f = ( std::pow(1.+z,C43) + std::pow(1.-z,C43) - 2. ) / CDENOM;
  double const z2 = z*z;
  double const z3 = z2*z;
  double const z4 = z3*z;


  double e0,e1,a;
  double D0Drs,D1Drs,DaDrs;
  {
    double Q0,Q1,Q1p,lterm;
    double const r12 = x;
    double const r1  = rs;
    double const r32 = r12*r1;
    double const r2  = r1*r1;

    Q0  = -2.*A0*(1.+A00*r1);
    Q1  =  2.*A0*( B01*r12 + B02*r1 + B03*r32 + B04*r2  );
    Q1p =  A0*( B01/r12 + 2.*B02 + 3.*B03*r12 + 4.*B04*r1 );
    lterm = std::log(1.+1./Q1);

    e0     = Q0 * lterm;
    D0Drs  = -2*A0*A00*lterm - Q0 * Q1p / ( Q1*Q1+Q1 );

    Q0  = -2.*A1*(1.+A10*r1);
    Q1  =  2.*A1*( B11*r12 + B12*r1 + B13*r32 + B14*r2  );
    Q1p =  A1*( B11/r12 + 2.*B12 + 3.*B13*r12 + 4.*B14*r1 );
    lterm = std::log(1.+1./Q1);

    e1     = Q0 * lterm;
    D1Drs  = -2*A1*A10*lterm - Q0 * Q1p / ( Q1*Q1+Q1 );


    Q0  = -2.*Aa*(1.+Aa0*r1);
    Q1  =  2.*Aa*( Ba1*r12 + Ba2*r1 + Ba3*r32 + Ba4*r2  );
    Q1p =  Aa*( Ba1/r12 + 2.*Ba2 + 3.*Ba3*r12 + 4.*Ba4*r1 );
    lterm = std::log(1.+1./Q1);

    a      = - ( Q0 * lterm );
    DaDrs  = - ( -2*Aa*Aa0*lterm - Q0 * Q1p / ( Q1*Q1+Q1 ) );
  };			       

  double const k = a/fpp0*(1.-z4) + (e1-e0)*z4;
  
  F[0] = xfa * ( pa13*mpa + pb13*mpb ) + p * ( e0 + f*k );
  
  double const p2 = p*p;
  double const DzDpa =  2.*pb[0]/p2;
  double const DrsDpa = -rs*C13/p;
  double const DfDz  = C43 * ( std::pow(1.+z,C13) - std::pow(1.-z,C13) ) / CDENOM;
  double const DfDpa = DfDz*DzDpa;
  double const D0Dpa = D0Drs*DrsDpa;
  double const D1Dpa = D1Drs*DrsDpa;
  double const DaDpa = DaDrs*DrsDpa;
  double const PkPpa = DaDpa/fpp0*(1.-z4) + (D1Dpa-D0Dpa)*z4;
  double const DkDz = a/fpp0*(-4.*z3) + (e1-e0)*(4.*z3);
  double const DkDpa = PkPpa + DkDz*DzDpa;

  F[1] = xfb * pa13 + e0 + p*D0Dpa + f*k + p*DfDpa*k + p*f*DkDpa;


  double const DzDpb = -2.*pa[0]/p2;
  double const DrsDpb = DrsDpa;
  double const DfDpb = DfDz*DzDpb;
  double const D0Dpb = D0Drs*DrsDpb;
  double const D1Dpb = D1Drs*DrsDpb;
  double const DaDpb = DaDrs*DrsDpb;
  double const PkPpb = DaDpb/fpp0*(1.-z4) + (D1Dpb-D0Dpb)*z4;
  double const DkDpb = PkPpb + DkDz*DzDpb;

  F[2] = xfb * pb13 + e0 + p*D0Dpb + f*k + p*DfDpb*k + p*f*DkDpb;

}



void rks_xcf_lda_spw92_potential( double const * pa, double * F )
{
  for ( int i=0; i<11; ++i ) F[i] = 0.;

  double const mpa = std::max( pa[0], 0. );
  double const p = 2.*mpa;

  if ( p < 1.e-100 ) return;

  //double const a = -Cx * 0.75 * ALPHA * std::pow(2.,C43);	       
  double const xfa = -9.305257363491000e-01;
  //double const b = -C43 * Cx * 1.5 * ALPHA * std::pow(2.,C13);
  double const xfb = -1.240700981798800e+00;
  double const C13 = 3.333333333333333e-01;
  double const pa13 = std::pow(mpa,C13);

  // PARAMAGENTIC 
  double const A0  = 0.0621814 / 2.;
  double const A00 = 0.21370;
  double const B01 = 7.5957;
  double const B02 = 3.5876;
  double const B03 = 1.6382;
  double const B04 = 0.49294;
  double const rs = std::pow( 0.75/(ccdl::PI*p), C13 );
  double const x = std::sqrt(rs);
  double e0;
  double D0Drs;
  {
    double Q0,Q1,Q1p,lterm;
    double const r12 = x;
    double const r1  = rs;
    double const r32 = r12*r1;
    double const r2  = r1*r1;
    Q0  = -2.*A0*(1.+A00*r1);
    Q1  =  2.*A0*( B01*r12 + B02*r1 + B03*r32 + B04*r2  );
    Q1p =  A0*( B01/r12 + 2.*B02 + 3.*B03*r12 + 4.*B04*r1 );
    lterm = std::log(1.+1./Q1);
    e0     = Q0 * lterm;
    D0Drs  = -2*A0*A00*lterm - Q0 * Q1p / ( Q1*Q1+Q1 );
  };			       
  F[0] = xfa * ( 2. * pa13*mpa ) + p * e0;
  double const DrsDpa = -rs*C13/p;
  double const D0Dpa = D0Drs*DrsDpa;
  F[1] = xfb * pa13 + e0 + p*D0Dpa;
  F[2] = F[1];
}


