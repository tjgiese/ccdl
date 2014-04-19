#include "harm.hpp"
#include <algorithm>
#include <cmath>

#ifndef alloca
#define alloca __builtin_alloca
#endif


void ccdl::SolidHarmRlm_v5
( int const Lmax, double const * R, double *__restrict__ Y )
{
  Y[0] = 1.;
  if ( Lmax > 0 )
    {	
      Y[1] =      R[2];
      Y[2] = -0.5*R[0];
      Y[3] = -0.5*R[1];
      double const r2 = R[0]*R[0]+R[1]*R[1]+R[2]*R[2];
      double *__restrict__ p = Y+3;
      double const *__restrict__  pa = Y;
      double const *__restrict__  pb = Y-1;
      double const tz = 2.*R[2];
      double tlm1 = R[2];
#ifdef __GNUG__
      double *__restrict__ ooi = (double *__restrict__) alloca ( sizeof(double)*(Lmax+Lmax) );
      ooi[1]=1.;
      for ( int i=2; i<Lmax+Lmax; ++i ) ooi[i] = 1./i;
      for ( int L=2; L<Lmax+1; ++L )
	{
	  tlm1 += tz;
	  *++p = ( tlm1*(*(++pa)) - r2*(*(++pb)) ) *ooi[L]*ooi[L];
	  for ( int M=1; M<L-1; ++M )
	    {
	      double mm = ooi[L+M]*ooi[L-M];
	      *++p = (tlm1*(*++pa)-r2*(*++pb))*mm;
	      *++p = (tlm1*(*++pa)-r2*(*++pb))*mm;
	    }
	  *++p = R[2]*(*++pa);
	  *++p = R[2]*(*++pa);
	  *++p = (Y[2]*(*(pa-1))-Y[3]*(*pa))*ooi[L];
	  *++p = (Y[2]*(*pa)+Y[3]*(*(pa-1)))*ooi[L];
	}
#else
      for ( int L=2; L<Lmax+1; ++L )
	{
	  double ool = 1./L;
	  tlm1 += tz;
	  *++p = ( tlm1*(*(++pa)) - r2*(*(++pb)) ) *ool*ool;
	  for ( int M=1; M<L-1; ++M )
	    {
	      double mm = 1./( (L+M) * (L-M));
	      *++p = (tlm1*(*++pa)-r2*(*++pb))*mm;
	      *++p = (tlm1*(*++pa)-r2*(*++pb))*mm;
	    }
	  *++p = R[2]*(*++pa);
	  *++p = R[2]*(*++pa);
	  *++p = (Y[2]*(*(pa-1))-Y[3]*(*pa))*ool;
	  *++p = (Y[2]*(*pa)+Y[3]*(*(pa-1)))*ool;
	}
#endif
    };   
}







void ccdl::SolidHarmIlm_v5
( int const Lmax, double const * R, double *__restrict__ Y )
{
  double const oor2 = 1. / ( R[0]*R[0] + R[1]*R[1] + R[2]*R[2] );
  Y[0] = std::sqrt(oor2);
  if ( Lmax > 0 )
    {
      double Q[3] = { -R[0]*oor2, -R[1]*oor2, R[2]*oor2 };
      Y[1] = Q[2]*Y[0];
      Y[2] = Q[0]*Y[0];
      Y[3] = Q[1]*Y[0];
      double *__restrict__ p = Y+3;
      double const *__restrict__  pa = Y;
      double const *__restrict__  pb = Y-1;
      for ( int L=2; L<Lmax+1; ++L )
	{
	  int tlm1 = L+L-1;
	  int l1 = (L-1)*(L-1);
	  double t = tlm1*Q[2];
	  *++p = *++pa*t - *++pb*l1*oor2;
	  for ( int M=1; M<L-1; ++M )
	    {
	      double h=(l1-M*M)*oor2;
	      *++p = *++pa*t - *++pb*h;
	      *++p = *++pa*t - *++pb*h;
	    };
	  *++p = *++pa*t;
	  *++p = *++pa*t;
	  
	  *++p = tlm1*(Q[0]*(*(pa-1))-Q[1]*(*pa));
	  *++p = tlm1*(Q[0]*(*pa)+Q[1]*(*(pa-1)));
	};
    };
}



void ccdl::RlmTranslation
( int const LtoMax, 
  int const LfromMax, 
  double const *__restrict__ Rtf,
  double *__restrict__ W )
{
  int const nt = ( LtoMax  +1 ) * ( LtoMax + 1 );
  int const nf = ( LfromMax+1 ) * ( LfromMax+1 );
  int const Lmin = std::min( LtoMax, LfromMax );
  for ( int i=0, n=nt*nf; i<n; ++i )
    W[i]=0.;
  for ( int i=0, n=(Lmin+1)*(Lmin+1); i<n; ++i )
    W[i+i*nt] = 1.;
  ccdl::SolidHarmRlm_v5( LtoMax, Rtf, W );
  
  int jend = (LtoMax <= LfromMax) ? LtoMax : LfromMax+1; 
  
  double s = 1.;
  for ( int j=1; j<jend; ++j )
    {
      s *= -1;
      // all but k=j
      for ( int k=0; k<2*j-1; ++k )
	{
	  int jk = j*j+k;
	  int pjk = (j-1)*(j-1)+k;
	  for ( int l=j+1; l<LtoMax+1; ++l )
	    {
	      int lm=l*l;
	      int plm=(l-1)*(l-1);
	      for ( int m=0; m<2*l-1; ++m, ++lm, ++plm )
		W[lm+jk*nt] = W[plm+pjk*nt];
	    };
	};
      // do k=j
      for ( int l=j+1; l<LtoMax+1; ++l )
	{
	  int k0c = l*l + (j*j+2*j-1)*nt;
	  int k0s = k0c+nt;
	  int o   = (l-j)*(l-j);
	  
	  if ( l >= j+j )
	    {
	      W[k0c] = s * 2. * W[o+2*j-1];
	      W[k0s] = s * 2. * W[o+2*j];
	    };
	  
	  for ( int m=1; m<l+1; ++m )
	    {
	      int cc = (2*m-1)+k0c;
	      int sc = (2*m  )+k0c;
	      int cs = (2*m-1)+k0s;
	      int ss = (2*m  )+k0s;
	      
	      int mj = m-j;
	      int amj = mj;
	      double csgn = 1.;
	      double ssgn = 1.;
	      if ( mj < 0 )
		{
		  amj = -mj;
		  csgn = (amj%2==1)? -1.: 1.;
		  ssgn = -csgn;
		};
	      if ( amj <= l-j )
		{
		  if ( amj > 0 )
		    {
		      // cc
		      W[cc] =  csgn*W[o+2*amj-1];
		      // sc
		      W[sc] =  ssgn*W[o+2*amj  ];
		      // cs
		      W[cs] = -ssgn*W[o+2*amj  ];
		      // ss
		      W[ss] =  csgn*W[o+2*amj-1];
		    }
		  else
		    {
		      // cc
		      W[cc] =  W[o];
		      // ss
		      W[ss] =  W[o];
		    };
		};
	      if ( m+j <= l-j )
		{
		  int O=o+2*(m+j);
		  // cc
		  W[cc] += s*W[O-1];
		  // sc
		  W[sc] += s*W[O];
		  // cs
		  W[cs] += s*W[O];
		  // ss
		  W[ss] -= s*W[O-1];
		};
	    };
	  
	};
    }
}




/*

void ccdl::SolidHarmIlm_v1
( int const Lmax, double const *__restrict__ R, double *__restrict__ Y )
{
  int const nl1 = Lmax+1;
  double x = R[0];
  double y = R[1];
  double z = R[2];
  double r2     = x*x + y*y + z*z;
  double oor2   = 1.0/r2;
  Y[0] = std::sqrt(oor2);
  double co  = Y[0];
  double coo = 0.0;  
  int TLm1=-1;
  int LPMtLMM=0;
  for (int L=1; L<nl1; ++L )
    {
      TLm1 += 2;
      double c  = oor2 * (TLm1*z*co-LPMtLMM*coo);
      coo= co;
      co = c;
      LPMtLMM += TLm1;
      Y[LPMtLMM] = c;
    };

  double cmmo=0.0;
  double cpmo=Y[0];

  for ( int M=1; M<nl1; ++M )
    {
      int TM = M+M;
      int TMm1=TM-1;
      int M2 = M*M;
      double cpm =  ( y*cmmo - x*cpmo ) * TMm1 * oor2;
      double cmm = -( y*cpmo + x*cmmo ) * TMm1 * oor2;
      cpmo=cpm;
      cmmo=cmm;
      int i  = M2+TM;
      Y[i-1] = cpm;
      Y[i] = cmm;
      double cpm2o = cpm;
      double cmm2o = cmm;
      double cpm2oo= 0.0;
      double cmm2oo= 0.0;
      LPMtLMM = -TMm1;
      for ( int L=M+1; L<nl1; ++L )
	{
	  int TL = L+L;
	  TLm1 = TL-1;
	  LPMtLMM += TLm1 - 2;
	  double a = TLm1 * z * oor2;
	  double b = LPMtLMM * oor2;
	  double cpm2 = a * cpm2o - b * cpm2oo;
	  double cmm2 = a * cmm2o - b * cmm2oo;
	  i += TLm1;
	  Y[i-1]   = cpm2;
	  Y[i] = cmm2;
	  cpm2oo = cpm2o;
	  cmm2oo = cmm2o;
	  cpm2o  = cpm2;
	  cmm2o  = cmm2;
	  };
      };
}




  void SolidHarmIlm_v2
  ( int const Lmax, double const * R, double *__restrict__ Y )
  {
    double const oor2 = 1. / ( R[0]*R[0] + R[1]*R[1] + R[2]*R[2] );
    Y[0] = std::sqrt(oor2);
    if ( Lmax > 0 )
      {
	double Q[3] = { -R[0]*oor2, -R[1]*oor2, R[2]*oor2 };
	Y[1] = Q[2]*Y[0];
	Y[2] = Q[0]*Y[0];
	Y[3] = Q[1]*Y[0];
	for ( int L=2; L<Lmax+1; ++L )
	  {
	    int tlm1 = L+L-1;
	    int l0 = L*L;
	    int l1 = l0-tlm1;
	    int l2 = l1-tlm1+2;
	    double t = tlm1*Q[2];
	    Y[l0] = t*Y[l1] - l1*Y[l2]*oor2;
	    for ( int M=1; M<L-1; ++M )
	      {
		Y[l0+2*M-1] = t*Y[l1+2*M-1]-(l1-M*M)*Y[l2+2*M-1]*oor2;
		Y[l0+2*M]   = t*Y[l1+2*M]  -(l1-M*M)*Y[l2+2*M]  *oor2;
	      };
	    Y[l0+tlm1-2] = t*Y[l1+tlm1-2];
	    Y[l0+tlm1-1] = t*Y[l1+tlm1-1];
	    
	    Y[l0+tlm1]   = tlm1*(Q[0]*Y[l1+tlm1-2]-Q[1]*Y[l1+tlm1-1]);
 	    Y[l0+tlm1+1] = tlm1*(Q[0]*Y[l1+tlm1-1]+Q[1]*Y[l1+tlm1-2]);
	  };
      };
  }


  void SolidHarmIlm_v3
  ( int const Lmax, double const * R, double *__restrict__ Y )
  {
    double const oor2 = 1. / ( R[0]*R[0] + R[1]*R[1] + R[2]*R[2] );
    Y[0] = std::sqrt(oor2);
    if ( Lmax > 0 )
      {
	double Q[3] = { -R[0]*oor2, -R[1]*oor2, R[2]*oor2 };
	Y[1] = Q[2]*Y[0];
	Y[2] = Q[0]*Y[0];
	Y[3] = Q[1]*Y[0];
	double *__restrict__ p = Y+3;
	double const *__restrict__  pa = Y;
	double const *__restrict__  pb = Y-1;
	for ( int L=2; L<Lmax+1; ++L )
	  {
	    int tlm1 = L+L-1;
	    //int l0 = L*L;
	    int l1 = (L-1)*(L-1);
	    //int l2 = l1-tlm1+2;
	    double t = tlm1*Q[2];
	    *++p = *++pa*t - *++pb*l1*oor2;
	    for ( int M=1; M<L-1; ++M )
	      {
		*++p = *++pa*t - *++pb*(l1-M*M)*oor2;
		*++p = *++pa*t - *++pb*(l1-M*M)*oor2;
	      };
	    *++p = *++pa*t;
	    *++p = *++pa*t;
	    
	    *++p = tlm1*(Q[0]*(*(pa-1))-Q[1]*(*pa));
 	    *++p = tlm1*(Q[0]*(*pa)+Q[1]*(*(pa-1)));
	  };
      };
  }

*/




/*

void ccdl::SolidHarmRlm_v1
( int const Lmax, double const *__restrict__ R, double *__restrict__ Y )
{
  Y[0] = 1.;

  int const nl1 = Lmax+1;
  double const x = R[0];
  double const y = R[1];
  double const z = R[2];
  double const r2 = x*x+y*y+z*z;

  double c;
  double co  = 1.0;
  double coo = 0.0;
  for ( int L=1; L < nl1; ++L )
    {
      int L2=L*L;
      c      = ( (2*L-1)*z*co-r2*coo ) / L2;
      coo   = co;
      co    = c;
      Y[L2] = c;
    };

  double cmmo=0.0;
  double cpmo=1.0;

  for ( int M=1; M < nl1; ++M )
    {
      int TM = 2*M;
      int M2 = M*M;
      double cpm =  ( y*cmmo - x*cpmo ) / TM;
      double cmm = -( y*cpmo + x*cmmo ) / TM;
      cpmo=cpm;
      cmmo=cmm;
      int i  = M2+TM;
      Y[i-1] = cpm;
      Y[i]   = cmm;

      double cpm2o = cpm;
      double cmm2o = cmm;
      double cpm2oo= 0.0;
      double cmm2oo= 0.0;
      int LPMtLMM = 0;
      for ( int L=M+1; L < nl1; ++L )
	{
	  int TL       = L+L;
	  int TLm1     = TL-1;
	  LPMtLMM += TLm1;
	  double a    = TLm1 * z / LPMtLMM;
	  double b    = r2 / LPMtLMM;
	  double cpm2 = a * cpm2o - b * cpm2oo;
	  double cmm2 = a * cmm2o - b * cmm2oo;
	  i += TLm1;
	  Y[i-1] = cpm2;
	  Y[i]   = cmm2;
	  cpm2oo   = cpm2o;
	  cmm2oo   = cmm2o;
	  cpm2o    = cpm2;
	  cmm2o    = cmm2;
	};
    };
  
}

void ccdl::SolidHarmRlm_v3
( int const Lmax, double const * R, double * Y )
{
  Y[0] = 1.;


  //0 1 2 3 4 5 6 7 8
  //0 1 1 1 2 2 2 2 2
  //0-1 0 1-2-1 0 1 2
  
  if ( Lmax > 0 )
    {
      double const x = R[0];
      double const y = R[1];
      double const z = R[2];

      Y[1] =      z;
      Y[2] = -0.5*x;
      Y[3] = -0.5*y;

      int const nl1 = Lmax+1;
      double const r2 = x*x+y*y+z*z;
      for ( int L=2; L<nl1; ++L )
	{
	  int l2 = (L-2)*(L-2);
	  int l1 = (L-1)*(L-1);
	  int l0 = L*L;
	  int tl = L+L;
	  double tlm1 = tl-1;
	  Y[l0] = ( tlm1*z*Y[l1] - r2*Y[l2] ) / (L*L);
	  int M,m;
	  double a,b;
	  
	  for ( M=1; M<L-1; ++M )
	    {
	      m = M+M;
	      double mm = (L+M)*(L-M);
	      a = tlm1*z/mm;
	      b = r2/mm;
	      Y[l0+m-1] = a*Y[l1+m-1]-b*Y[l2+m-1];
	      Y[l0+m] = a*Y[l1+m]-b*Y[l2+m];
	    }
	  M=L-1;
	  m=M+M;
	  a = tlm1*z / ( (L+M)*(L-M) );
	  Y[l0+m-1]=a*Y[l1+m-1];
	  Y[l0+m]=a*Y[l1+m];
	  M=L;
	  m=M+M;
	  Y[l0+m-1]=-(x*Y[l1+(m-2)-1]-y*Y[l1+(m-2)])/tl;
	  Y[l0+m]=-(x*Y[l1+(m-2)]+y*Y[l1+(m-2)-1])/tl;
	}
    };
}

void ccdl::SolidHarmRlm_v4
( int const Lmax, double const * R, double *__restrict__ Y )
{
  Y[0] = 1.;
  if ( Lmax > 0 )
    {
      double const x = R[0];
      double const y = R[1];
      double const z = R[2];
	
      Y[1] =      z;
      Y[2] = -0.5*x;
      Y[3] = -0.5*y;
	
      int const nl1 = Lmax+1;
      double const r2 = x*x+y*y+z*z;
      int M;
      //double a,b;
      double *__restrict__ p = Y+3;
	
      for ( int L=2; L<nl1; ++L )
	{
	  double ool0 = 1./(L*L);
	  double ootl = 0.5/L;
	  int o = L+L-1;
	  double const tlm1=o*z;
	  double const *__restrict__  pa = p-(--o);
	  double const *__restrict__  pb = pa-(--o);
	  *++p = ( tlm1*(*pa) - r2*(*pb) ) * ool0;
	  for ( M=1; M<L-1; ++M )
	    {
	      double mm = 1./((L+M)*(L-M));
	      *++p = (tlm1*(*++pa)-r2*(*++pb))*mm;
	      *++p = (tlm1*(*++pa)-r2*(*++pb))*mm;
	      // a = tlm1*mm;
	      // b = r2*mm;
	      // *++p = a*(*++pa)-b*(*++pb);
	      // *++p = a*(*++pa)-b*(*++pb);
	    }
	  *++p = z*(*++pa);
	  *++p = z*(*++pa);
	  pb=pa-1;
	  *++p = -(x*(*(pb))-y*(*pa))*ootl;
	  *++p = -(x*(*pa)+y*(*(pb)))*ootl;
	}
    };   
}

*/
