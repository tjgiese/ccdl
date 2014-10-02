#include "harm.hpp"
#include <algorithm>
#include <cmath>

#ifndef alloca
#define alloca __builtin_alloca
#endif


void ccdl::SolidHarm_Rlm
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



void ccdl::SolidHarm_dRlm
( int const lmax,
  double const *__restrict__ Y,
  double *__restrict__ dY )
{
  // l = 0
  dY[0] = 0.; dY[1] = 0.; dY[2] = 0.;

  if ( lmax > 0 )
    {
      // l=1
      dY[3] =  0.;       dY[4]  =  0.;      dY[5]  = Y[0];
      dY[6] = -0.5*Y[0]; dY[7]  =  0.;      dY[8]  = 0.;
      dY[9] =  0.;      dY[10] = -0.5*Y[0]; dY[11] = 0.;
      if ( lmax > 1 )
	{
	  // l=2
	  dY[12] =      Y[2];  dY[13]  =  Y[3];      dY[14]  = Y[1];
	  dY[15] = -0.5*Y[1];  dY[16]  =  0.;        dY[17]  = Y[2];
	  dY[18] =  0.;        dY[19]  = -0.5*Y[1];  dY[20]  = Y[3];
	  dY[21] = -0.5*Y[2];  dY[22]  =  0.5*Y[3];  dY[23]  = 0.;
	  dY[24] = -0.5*Y[3];  dY[25]  = -0.5*Y[2];  dY[26]  = 0.;

	  double *__restrict__ p = dY+26;
	  for ( int l=3; l<lmax+1; ++l )
	    {
	      int lo = (l-1)*(l-1);
		
	      // M=0
	      *(++p) = Y[lo+1];
	      *(++p) = Y[lo+2];
	      *(++p) = Y[lo  ];
		
	      // M=1 (no l-1,m-1s; l-1,m-1c is 0 case)
	      *(++p) =  0.5 * ( Y[lo+3]-Y[lo] );
	      *(++p) =  0.5 * Y[lo+4];
	      *(++p) =  Y[lo+1];
	      *(++p) =  0.5 * Y[lo+4];
	      *(++p) = -0.5 * ( Y[lo+3]+Y[lo] );
	      *(++p) =  Y[lo+2];
		
	      // 2 <= M < L-2
	      int ms = 4;
	      int mc = 3;
	      for ( int  m=2; m < l-1; ++m )
		{
		  //int ms = 2*m;
		  //int mc = ms-1;
		  *(++p) =  0.5 * ( Y[lo+mc+2]-Y[lo+mc-2] );
		  *(++p) =  0.5 * ( Y[lo+ms+2]+Y[lo+ms-2] );
		  *(++p) = Y[lo+mc];
		  *(++p) =  0.5 * ( Y[lo+ms+2]-Y[lo+ms-2] );
		  *(++p) = -0.5 * ( Y[lo+mc+2]+Y[lo+mc-2] );
		  *(++p) = Y[lo+ms];
		  ms += 2;
		  mc += 2;
		}
		
	      // M = L-1 (no l-1,m+1)
	      //int ms = 2*(l-1);
	      //int mc = ms-1;
	      *(++p) = -0.5 * Y[lo+mc-2];
	      *(++p) =  0.5 * Y[lo+ms-2];
	      *(++p) = Y[lo+mc];
	      *(++p) = -0.5 * Y[lo+ms-2];
	      *(++p) = -0.5 * Y[lo+mc-2];
	      *(++p) = Y[lo+ms];
		
	      // M = L  (no l-1,m+1; no l-1,m)
	      //ms += 2;
	      //mc += 2;
	      *(++p) = -0.5 * Y[lo+mc];
	      *(++p) =  0.5 * Y[lo+ms];
	      *(++p) =  0.;
	      *(++p) = -0.5 * Y[lo+ms];
	      *(++p) = -0.5 * Y[lo+mc];
	      *(++p) =  0.;
	    };
	};
    };
}





void ccdl::SolidHarm_d2Rlm
( int const lmax,
  double const *__restrict__ dY,
  double *__restrict__ d2Y )
{
  --d2Y;
  for ( int i=0; i<6; ++i ) *(++d2Y) = 0.;

  if ( lmax > 0 )
    {
      for ( int i=6; i<6*4; ++i ) *(++d2Y) = 0.;

      if ( lmax > 1 )
	{

	  *(++d2Y) = dY[6];
	  *(++d2Y) = dY[9];
	  *(++d2Y) = dY[3];
	  *(++d2Y) = dY[10];
	  *(++d2Y) = dY[4];
	  *(++d2Y) = dY[5];

	  *(++d2Y) = -0.5 * dY[3];
	  *(++d2Y) = 0.;
	  *(++d2Y) = dY[6];
	  *(++d2Y) = 0.;
	  *(++d2Y) = dY[7];
	  *(++d2Y) = dY[8];

	  *(++d2Y) = 0.;
	  *(++d2Y) = -0.5 * dY[3];
	  *(++d2Y) = dY[9];
	  *(++d2Y) = -0.5 * dY[4];
	  *(++d2Y) = dY[10];
	  *(++d2Y) = dY[11];

	  *(++d2Y) = -0.5 * dY[6];
	  *(++d2Y) =  0.5 * dY[9];
	  *(++d2Y) = 0.;
	  *(++d2Y) =  0.5 * dY[10];
	  *(++d2Y) = 0.;
	  *(++d2Y) = 0.;

	  *(++d2Y) = -0.5 * dY[9];
	  *(++d2Y) = -0.5 * dY[6];
	  *(++d2Y) = 0.;
	  *(++d2Y) = -0.5 * dY[7];
	  *(++d2Y) = 0.;
	  *(++d2Y) = 0.;

	  // here

	  for ( int L=3; L <= lmax; ++L )
	    {
	      int lo = (L-1)*(L-1);
	      double const *__restrict__ q = dY + lo*3;
	      // M=0
	      *(++d2Y) = q[3];
	      *(++d2Y) = q[6];
	      *(++d2Y) = q[0];
	      *(++d2Y) = q[7];
	      *(++d2Y) = q[1];
	      *(++d2Y) = q[2];



	      // M=1 (no l-1,m-1s; l-1,m-1c is 0 case)
	      *(++d2Y) = 0.5 * ( q[9]-q[0] );
	      *(++d2Y) = 0.5 * q[12];
	      *(++d2Y) = q[3];
	      *(++d2Y) = 0.5 * q[13];
	      *(++d2Y) = q[4];
	      *(++d2Y) = q[5];

	      *(++d2Y) =  0.5 * q[12];
	      *(++d2Y) = -0.5 * ( q[9]+q[0] );
	      *(++d2Y) = q[6];
	      *(++d2Y) = -0.5 * ( q[10]+q[1] );
	      *(++d2Y) = q[7];
	      *(++d2Y) = q[8];



	      // 2 <= M < L-2
	      //int mc = 3;
	      int tmc = 9;
	      for ( int m=2; m < L-1; ++m )
		{
		  *(++d2Y) = 0.5 * ( q[tmc+6]-q[tmc-6] );
		  *(++d2Y) = 0.5 * ( q[tmc+9]+q[tmc-3] );
		  *(++d2Y) = q[tmc];
		  *(++d2Y) = 0.5 * ( q[tmc+10]+q[tmc-2] );
		  *(++d2Y) = q[tmc+1];
		  *(++d2Y) = q[tmc+2];

		  *(++d2Y) =  0.5 * ( q[tmc+9]-q[tmc-3] );
		  *(++d2Y) = -0.5 * ( q[tmc+6]+q[tmc-6] );
		  *(++d2Y) = q[tmc+3];
		  *(++d2Y) = -0.5 * ( q[tmc+7]+q[tmc-5] );
		  *(++d2Y) = q[tmc+4];
		  *(++d2Y) = q[tmc+5];

		  //mc += 2;
		  tmc += 6;
		}

	      // M = L-1 (no l-1,m+1)
	      //ms = 2*(L-1);
	      //mc = ms-1;

	      *(++d2Y) = -0.5 * q[tmc-6];
	      *(++d2Y) =  0.5 * q[tmc-3];
	      *(++d2Y) = q[tmc];
	      *(++d2Y) =  0.5 * q[tmc-2];
	      *(++d2Y) = q[tmc+1];
	      *(++d2Y) = q[tmc+2];

	      *(++d2Y) = -0.5 * q[tmc-3];
	      *(++d2Y) = -0.5 * q[tmc-6];
	      *(++d2Y) = q[tmc+3];
	      *(++d2Y) = -0.5 * q[tmc-5];
	      *(++d2Y) = q[tmc+4];
	      *(++d2Y) = q[tmc+5];


	      // M = L  (no l-1,m+1; no l-1,m)
	      //ms = 2*L;
	      //mc = ms-1;

	      *(++d2Y) = -0.5 * q[tmc];
	      *(++d2Y) =  0.5 * q[tmc+3];
	      *(++d2Y) = 0.;
	      *(++d2Y) =  0.5 * q[tmc+4];
	      *(++d2Y) = 0.;
	      *(++d2Y) = 0.;

	      *(++d2Y) = -0.5 * q[tmc+3];
	      *(++d2Y) = -0.5 * q[tmc];
	      *(++d2Y) = 0.;
	      *(++d2Y) = -0.5 * q[tmc+1];
	      *(++d2Y) = 0.;
	      *(++d2Y) = 0.;
	    }


	}
    }


}






void ccdl::RlmTranslation
( int const LtoMax, 
  int const LfromMax,
  double const *__restrict__ Rft,
  double *__restrict__ W )
{
  int const nt = ( LtoMax  +1 ) * ( LtoMax + 1 );
  int const nf = ( LfromMax+1 ) * ( LfromMax+1 );
  int const Lmin = std::min( LtoMax, LfromMax );
  for ( int i=nt, n=nt*nf; i<n; ++i )
    W[i]=0.;
  for ( int i=0, n=(Lmin+1)*(Lmin+1); i<n; ++i )
    W[i+i*nt] = 1.;
  ccdl::SolidHarm_Rlm( LtoMax, Rft, W );
  
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



// void ccdl::RlmTranslationFromCrd
// ( int const LtoMax, 
//   int const LfromMax, 
//   double const *__restrict__ Rft,
//   double *__restrict__ W )
// {
//   ccdl::SolidHarm_Rlm( LtoMax, Rft, W );
//   ccdl::RlmTranslationFromRlm( LtoMax, LfromMax, W );
// }





void ccdl::RlmTranslationGrd
( int const lam, int const lbm, 
  double const *__restrict__ Rab,
  double *__restrict__ W,
  double *__restrict__ D )
{

  ccdl::RlmTranslation(lam,lbm,Rab,W);

  // l = 0
  int na = (lam+1)*(lam+1);
  int nb = (lbm+1)*(lbm+1);
  int nd = na*nb*3;
  for ( int i=0; i<nd; ++i ) D[i] = 0.;

  if ( lam < 1 ) return;

  int lmin = std::min(lam,lbm);


  // l=1
  D[3] =  0.;       D[4]  =  0.;      D[5]  = W[0];
  D[6] = -0.5*W[0]; D[7]  =  0.;      D[8]  = 0.;
  D[9] =  0.;      D[10] = -0.5*W[0]; D[11] = 0.;
  if ( lam > 1 )
    {
      // l=2
      D[12] =      W[2];  D[13]  =  W[3];      D[14]  = W[1];
      D[15] = -0.5*W[1];  D[16]  =  0.;        D[17]  = W[2];
      D[18] =  0.;        D[19]  = -0.5*W[1];  D[20]  = W[3];
      D[21] = -0.5*W[2];  D[22]  =  0.5*W[3];  D[23]  = 0.;
      D[24] = -0.5*W[3];  D[25]  = -0.5*W[2];  D[26]  = 0.;
	
      double *__restrict__ p = D+26;
      for ( int l=3; l<lam+1; ++l )
	{
	  int lo = (l-1)*(l-1);
            
	  // M=0
	  *(++p) = W[lo+1];
	  *(++p) = W[lo+2];
	  *(++p) = W[lo  ];
            
	  // M=1 (no l-1,m-1s; l-1,m-1c is 0 case)
	  *(++p) =  0.5 * ( W[lo+3]-W[lo] );
	  *(++p) =  0.5 * W[lo+4];
	  *(++p) =  W[lo+1];
	  *(++p) =  0.5 * W[lo+4];
	  *(++p) = -0.5 * ( W[lo+3]+W[lo] );
	  *(++p) =  W[lo+2];
            
	  // 2 <= M < L-2
	  int ms = 4;
	  int mc = 3;
	  for ( int  m=2; m < l-1; ++m )
	    {
	      //int ms = 2*m;
	      //int mc = ms-1;
	      *(++p) =  0.5 * ( W[lo+mc+2]-W[lo+mc-2] );
	      *(++p) =  0.5 * ( W[lo+ms+2]+W[lo+ms-2] );
	      *(++p) = W[lo+mc];
	      *(++p) =  0.5 * ( W[lo+ms+2]-W[lo+ms-2] );
	      *(++p) = -0.5 * ( W[lo+mc+2]+W[lo+mc-2] );
	      *(++p) = W[lo+ms];
	      ms += 2;
	      mc += 2;
	    }
	    
	  // M = L-1 (no l-1,m+1)
	  //int ms = 2*(l-1);
	  //int mc = ms-1;
	  *(++p) = -0.5 * W[lo+mc-2];
	  *(++p) =  0.5 * W[lo+ms-2];
	  *(++p) = W[lo+mc];
	  *(++p) = -0.5 * W[lo+ms-2];
	  *(++p) = -0.5 * W[lo+mc-2];
	  *(++p) = W[lo+ms];
            
	  // M = L  (no l-1,m+1; no l-1,m)
	  //ms += 2;
	  //mc += 2;
	  *(++p) = -0.5 * W[lo+mc];
	  *(++p) =  0.5 * W[lo+ms];
	  *(++p) =  0.;
	  *(++p) = -0.5 * W[lo+ms];
	  *(++p) = -0.5 * W[lo+mc];
	  *(++p) =  0.;
	};
	
    };


	  

  if ( lam > 1 and lbm > 0 )
    {
      D[2+(4+(1)*na)*3] =  1.0;
      D[0+(5+(1)*na)*3] = -0.5;
      D[1+(6+(1)*na)*3] = -0.5;
	
      D[0+(4+(2)*na)*3] =  1.0;
      D[2+(5+(2)*na)*3] =  1.0;
      D[0+(7+(2)*na)*3] = -0.5;
      D[1+(8+(2)*na)*3] = -0.5;
	
      D[1+(4+(3)*na)*3] =  1.0;
      D[2+(6+(3)*na)*3] =  1.0;
      D[1+(7+(3)*na)*3] =  0.5;
      D[0+(8+(3)*na)*3] = -0.5;
    };
    
  for ( int l=3; l <= lam; ++l )
    {
      int l0 = l*l;
      int l1 = (l-1)*(l-1);

      int jmx = std::min(l-1,lmin);
      for ( int j=1; j <= jmx; ++j )
	{
	  int j0  = j*j;
	  int j1 = (j-1)*(j-1);
	  int j2 = (j+1)*(j+1);
	    
	  for ( int k=0; k < 2*j-1; ++k )
	    for ( int m=0; m < 2*l-1; ++m )
	      {
		D[0+(l0+m+(j0+k)*na)*3] = D[0+(l1+m+(j1+k)*na)*3];
		D[1+(l0+m+(j0+k)*na)*3] = D[1+(l1+m+(j1+k)*na)*3];
		D[2+(l0+m+(j0+k)*na)*3] = D[2+(l1+m+(j1+k)*na)*3];
	      };

	  // Now only do k=j cases;
	  int l0a = l0+(j2-2)*na;
	  int l1a = l1+(j2-2)*na;
	  for ( int jk= j2-2 ; jk < j2; ++jk, l0a += na, l1a += na )
	    {
	      //double *__restrict__ p = D+l0a*3-1;
	      // m=0;
	      D[0+(l0a)*3] = W[l1a+1];
	      D[1+(l0a)*3] = W[l1a+2];
	      D[2+(l0a)*3] = W[l1a+0];
		
	      // m=1c;
	      D[0+(l0a+1)*3] =  0.5 * W[l1a+3]- 0.5*W[l1a+0];
	      D[1+(l0a+1)*3] =  0.5 * W[l1a+4];
	      D[2+(l0a+1)*3] = W[l1a+1];

	      // m=1s;
	      D[0+(l0a+2)*3] =  0.5 * W[l1a+4];
	      D[1+(l0a+2)*3] = -0.5 * W[l1a+3]- 0.5*W[l1a+0];
	      D[2+(l0a+2)*3] = W[l1a+2];

	      // m>1, 0 < m+1 < l-1;
	      for ( int m=2; m < l-1; ++m )
		{
		  D[0+(l0a+2*m-1)*3] =  0.5 * W[l1a+2*(m+1)-1] - 0.5 * W[l1a+2*(m-1)-1];
		  D[1+(l0a+2*m-1)*3] =  0.5 * W[l1a+2*(m+1)  ] + 0.5 * W[l1a+2*(m-1)  ];
		  D[2+(l0a+2*m-1)*3] = W[l1a+2*m-1];
		  D[0+(l0a+2*m  )*3] =  0.5 * W[l1a+2*(m+1)  ] - 0.5 * W[l1a+2*(m-1)  ];
		  D[1+(l0a+2*m  )*3] = -0.5 * W[l1a+2*(m+1)-1] - 0.5 * W[l1a+2*(m-1)-1];
		  D[2+(l0a+2*m  )*3] = W[l1a+2*m  ];
		};

	      // m=l-1;
	      int m=l-1;
	      D[0+(l0a+2*m-1)*3] =  - 0.5 * W[l1a+2*(m-1)-1];
	      D[1+(l0a+2*m-1)*3] =    0.5 * W[l1a+2*(m-1)  ];
	      D[2+(l0a+2*m-1)*3] = W[l1a+2*m-1];

	      D[0+(l0a+2*m  )*3] = - 0.5 * W[l1a+2*(m-1)  ];
	      D[1+(l0a+2*m  )*3] = - 0.5 * W[l1a+2*(m-1)-1];
	      D[2+(l0a+2*m  )*3] = W[l1a+2*m  ];

	      m=l;
	      D[0+(l0a+2*m-1)*3] =  - 0.5 * W[l1a+2*(m-1)-1];
	      D[1+(l0a+2*m-1)*3] =    0.5 * W[l1a+2*(m-1)  ];

	      D[0+(l0a+2*m  )*3] = - 0.5 * W[l1a+2*(m-1)  ];
	      D[1+(l0a+2*m  )*3] = - 0.5 * W[l1a+2*(m-1)-1];

	    };
	};
    };
    
}










void ccdl::SolidHarm_Ilm
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





void ccdl::SolidHarm_dIlm
( int const lmax,
  double const *__restrict__ Y,
  double *__restrict__ dY )
{
  // L=0
  *(dY) =  Y[2];
  *(++dY) =  Y[3];
  *(++dY) = -Y[1];


  double const *__restrict__ q = Y+4;
  for ( int L=1; L < lmax+1; ++L )
    {
      //int l1 = (L+1)*(L+1);
      // double const *__restrict__ q = Y+l1;
      // m=0
      *(++dY) =  q[1];
      *(++dY) =  q[2];
      *(++dY) = -*q;

      // m=1
      *(++dY) =  0.5 * ( q[3] - *q );
      *(++dY) =  0.5 * q[4];
      *(++dY) = -q[1];

      *(++dY) =  0.5 * q[4];
      *(++dY) = -0.5 * ( q[3] + *q );
      *(++dY) = -q[2];

      // m>1
      // for ( int m=2; m < L+1; ++m )
      //   {
      //     *(++dY) =  0.5 * ( q[2*m+1] - q[2*m-3] );
      //     *(++dY) =  0.5 * ( q[2*m+2] + q[2*m-2] );
      //     *(++dY) = -q[2*m-1];
      //     *(++dY) =  0.5 * ( q[2*m+2] - q[2*m-2] );
      //     *(++dY) = -0.5 * ( q[2*m+1] + q[2*m-3] );
      //     *(++dY) = -q[2*m  ];
      //   }
      //q += 2*L+3;

      q += 4;
      for ( int m=2; m < L+1; ++m )
	{
	  *(++dY) =  0.5 * ( q[1] - q[-3] );
	  *(++dY) =  0.5 * ( q[2] + q[-2] );
	  *(++dY) = -q[-1];
	  *(++dY) =  0.5 * ( q[2] - q[-2] );
	  *(++dY) = -0.5 * ( q[1] + q[-3] );
	  *(++dY) = -*q;
	  q += 2;
	}
      ++q;

    }
}





void ccdl::SolidHarm_d2Ilm
( int const lmax,
  double const *__restrict__ dY,
  double *__restrict__ d2Y )
{

  // L=0
  *(d2Y)   =  dY[6];
  *(++d2Y) =  dY[9];
  *(++d2Y) = -dY[3];
  *(++d2Y) =  dY[10];
  *(++d2Y) = -dY[4];
  *(++d2Y) = -dY[5];

  for ( int L=1; L < lmax+1; ++L )
    {
      //int l0 = L*L;
      int l1 = (L+1)*(L+1);
      double const *__restrict__ q = dY+3*l1;
      // m=0
      *(++d2Y) =  q[3];
      *(++d2Y) =  q[6];
      *(++d2Y) = -*q;
      *(++d2Y) =  q[7];
      *(++d2Y) = -q[1];
      *(++d2Y) = -q[2];

      // m=1
      *(++d2Y) =  0.5 * ( q[9] - *q );
      *(++d2Y) =  0.5 * q[12];
      *(++d2Y) = -q[3];
      *(++d2Y) =  0.5 * q[13];
      *(++d2Y) = -q[4];
      *(++d2Y) = -q[5];

      *(++d2Y) =  0.5 * q[12];
      *(++d2Y) = -0.5 * ( q[9] + *q );
      *(++d2Y) = -q[6];
      *(++d2Y) = -0.5 * ( q[10] + q[1] );
      *(++d2Y) = -q[7];
      *(++d2Y) = -q[8];
      // m>1
      q += 12;
      for ( int m=2; m < L+1; ++m )
	{
	  // *(++d2Y) =  0.5 * ( q[6*m+3] - q[6*m-9] );
	  // *(++d2Y) =  0.5 * ( q[6*m+6] + q[6*m-6] );
	  // *(++d2Y) = -q[6*m-3];
	  // *(++d2Y) =  0.5 * ( q[6*m+7] + q[6*m-5] );
	  // *(++d2Y) = -q[6*m-2];
	  // *(++d2Y) = -q[6*m-1];


	  // *(++d2Y) =  0.5 * ( q[6*m+6] - q[6*m-6] );
	  // *(++d2Y) = -0.5 * ( q[6*m+3] + q[6*m-9] );
	  // *(++d2Y) = -q[6*m];
	  // *(++d2Y) = -0.5 * ( q[6*m+4] + q[6*m-8] );
	  // *(++d2Y) = -q[6*m+1];
	  // *(++d2Y) = -q[6*m+2];

	  *(++d2Y) =  0.5 * ( q[3] - q[-9] );
	  *(++d2Y) =  0.5 * ( q[6] + q[-6] );
	  *(++d2Y) = -q[-3];
	  *(++d2Y) =  0.5 * ( q[7] + q[-5] );
	  *(++d2Y) = -q[-2];
	  *(++d2Y) = -q[-1];


	  *(++d2Y) =  0.5 * ( q[6] - q[-6] );
	  *(++d2Y) = -0.5 * ( q[3] + q[-9] );
	  *(++d2Y) = -*q;
	  *(++d2Y) = -0.5 * ( q[4] + q[-8] );
	  *(++d2Y) = -q[1];
	  *(++d2Y) = -q[2];
	    
	  q += 6;
	}
    }
}



/*
void ccdl::IlmInteraction
( int const lam, int const lbm, 
  double const *__restrict__ Rab,
  double *__restrict__ T )
{
  int l,m,j,k;
  int l0,j0,lj0,l1,j1,jc,js;
  int jms,kms,mms;


  if ( lam == 0 and lbm == 0 )
    {
      T[0] = 1. / std::sqrt( Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2] );
      return;
    };

  int const na = (lam+1)*(lam+1);
#ifdef __GNUG__
  double *__restrict__ I = (double *__restrict__) alloca ( sizeof(double)*((lam+lbm+1)*(lam+lbm+1)) );
#else
  std::vector<double> I( (lam+lbm+1)*(lam+lbm+1) , 0. );
#endif

  // minus Rab !!
  double mrab[] = { -Rab[0], -Rab[1], -Rab[2] };
  ccdl::SolidHarm_Ilm(lam+lbm,mrab,I);





  if ( lam >= lbm )
    {

      for ( l=0; l <= lam; ++l )
	{
	  l0 = l*l;
	  T[l0] = I[l0];
	  for ( m=1 ; m < 2*l+1; ++m )
	    T[l0+m] = I[l0+m] * 2.;
	};

      jms = -1;
      for ( j=1; j <= lbm; ++j )
	{
	  j0 = j*j;
	  T[j0*na] = jms * I[j0];
	  for ( k=1; k < 2*j+1; ++k )
	    T[(j0+k)*na] = jms * I[j0+k] * 2.;
	  jms = -jms;
	};
	
      for ( j=1; j <= lbm; ++j )
	{
	  j0 = j*j;


	  jms = -std::pow(-1,j);
	  int lmx = std::min(j-1,lam);
	  for ( l=1; l <= lmx; ++l )
	    {
	      l0 = l*l;
	      for ( m=0; m < 2*l+1; ++m )
		for ( k=0; k < 2*j+1; ++k )
		  T[(l0+m)+(j0+k)*na] = jms * T[(j0+k)+(l0+m)*na];
	      jms = -jms;
	    };
          

	  //  we can build partially build l=j,lam-1
	  //  from the l+1,j-1 block... we will just need to compute the k=j cases for all m


	  jms = std::pow(-1,j);

	  j1 = (j-1)*(j-1);
	  for ( l=j; l <= lam-1; ++l )
	    {
	      l0 = l*l;
	      l1 = (l+1)*(l+1);
	      lj0 = (l+j)*(l+j);
		
	      //  k != j
	      for ( k=0; k < 2*j-1; ++k )
		for ( m=0; m < 2*l+1; ++m )
		  T[(l0+m)+(j0+k)*na] = -T[(l1+m)+(j1+k)*na];
		
	      //  k = j
		
	      jc = j0+2*j-1;
	      js = jc+1;
		
	      //  m=0, k=j
	      //  (m-k < 0)
	      T[l0+jc*na] = jms * 2. * I[lj0+2*j-1];
	      T[l0+js*na] = jms * 2. * I[lj0+2*j  ];
		
		
	      //  m=j, k=j
	      //  (m-k = 0)
	      //  (m = k)
	      T[(l0+2*j-1)+jc*na] = jms * 2. * ( I[lj0+4*j-1] + jms * I[lj0] );
	      T[(l0+2*j-1)+js*na] = jms * 2. * ( I[lj0+4*j  ] );
	      T[(l0+2*j  )+jc*na] = jms * 2. * ( I[lj0+4*j  ] );
	      T[(l0+2*j  )+js*na] = jms * 2. * (-I[lj0+4*j-1] + jms * I[lj0] );
		
		
	      //  m>0, k>0
	      //  m-k < 0
	      mms = -1;
	      int mmx = std::min(l,j-1);
	      for ( m=1; m <= mmx; ++m )
		{
		  T[(l0+2*m-1)+jc*na] = jms * 2. * ( I[lj0+2*(m+j)-1] + mms * I[lj0-2*(m-j)-1] );
		  T[(l0+2*m-1)+js*na] = jms * 2. * ( I[lj0+2*(m+j)  ] + mms * I[lj0-2*(m-j)  ] );
		  T[(l0+2*m  )+jc*na] = jms * 2. * ( I[lj0+2*(m+j)  ] - mms * I[lj0-2*(m-j)  ] );
		  T[(l0+2*m  )+js*na] = jms * 2. * (-I[lj0+2*(m+j)-1] + mms * I[lj0-2*(m-j)-1] );
		  mms = -mms;
		}
		
	      //  m>0, k>0
	      //  m-k > 0
	      for ( m=j+1; m <= l; ++m )
		{
		  T[(l0+2*m-1)+jc*na] = jms * 2. * ( I[lj0+2*(m+j)-1] + jms * I[lj0+2*(m-j)-1] );
		  T[(l0+2*m-1)+js*na] = jms * 2. * ( I[lj0+2*(m+j)  ] - jms * I[lj0+2*(m-j)  ] );
		  T[(l0+2*m  )+jc*na] = jms * 2. * ( I[lj0+2*(m+j)  ] + jms * I[lj0+2*(m-j)  ] );
		  T[(l0+2*m  )+js*na] = jms * 2. * (-I[lj0+2*(m+j)-1] + jms * I[lj0+2*(m-j)-1] );
		}
	    };
	    
	    
	    
	  //  we have to build the entire l=lam block
	    

	  l = lam;
	    
	  l0  = l*l;
	  lj0 = (l+j)*(l+j);
	    
	  //  m=0, k=0
	  T[l0+j0*na] = jms * I[lj0];
	    
	  //  m=0, k>0
	  //  (m-k < 0)
	  for ( k=1; k <= j; ++k )
	    {
	      T[l0+(j0+2*k-1)*na] = jms * 2. * I[lj0+2*k-1];
	      T[l0+(j0+2*k  )*na] = jms * 2. * I[lj0+2*k];
	    };
	  //  m>0, k=0
	  //  (m-k > 0)
	  for ( m=1; m <= l; ++m )
	    {
	      T[(l0+2*m-1)+j0*na] = jms * 2. * I[lj0+2*m-1];
	      T[(l0+2*m  )+j0*na] = jms * 2. * I[lj0+2*m  ];
	    }
	    
	  //  m>0, k>0
	  //  (m-k = 0)
	  //  (m = k)
	  mms = -1;
	  for ( m=1; m <= j; ++m )
	    {
	      T[(l0+2*m-1)+(j0+2*m-1)*na] = jms * 2. * ( I[lj0+4*m-1] + mms * I[lj0] );
	      T[(l0+2*m-1)+(j0+2*m  )*na] = jms * 2. * ( I[lj0+4*m  ]);
	      T[(l0+2*m  )+(j0+2*m-1)*na] = jms * 2. * ( I[lj0+4*m  ]);
	      T[(l0+2*m  )+(j0+2*m  )*na] = jms * 2. * (-I[lj0+4*m-1] + mms * I[lj0] );
	      mms = -mms;
	    };
	    
	    
	  //  m>0, k>0
	  //  m-k < 0
	  mms = -1;
	  for ( m=1; m <= l; ++m )
	    {
	      for ( k=m+1; k <= j; ++k )
		{
		  T[(l0+2*m-1)+(j0+2*k-1)*na] = jms * 2. * ( I[lj0+2*(m+k)-1] + mms * I[lj0-2*(m-k)-1] );
		  T[(l0+2*m-1)+(j0+2*k  )*na] = jms * 2. * ( I[lj0+2*(m+k)  ] + mms * I[lj0-2*(m-k)  ] );
		  T[(l0+2*m  )+(j0+2*k-1)*na] = jms * 2. * ( I[lj0+2*(m+k)  ] - mms * I[lj0-2*(m-k)  ] );
		  T[(l0+2*m  )+(j0+2*k  )*na] = jms * 2. * (-I[lj0+2*(m+k)-1] + mms * I[lj0-2*(m-k)-1] );
		};
	      mms = -mms;
	    }
	    
	  //  m>0, k>0
	  //  m-k > 0
	  kms = -1;
	  for ( k=1; k <= j; ++k )
	    {
	      for ( m=k+1; m <= l; ++m )
		{
		  T[(l0+2*m-1)+(j0+2*k-1)*na] = jms * 2. * ( I[lj0+2*(m+k)-1] + kms * I[lj0+2*(m-k)-1] );
		  T[(l0+2*m-1)+(j0+2*k  )*na] = jms * 2. * ( I[lj0+2*(m+k)  ] - kms * I[lj0+2*(m-k)  ] );
		  T[(l0+2*m  )+(j0+2*k-1)*na] = jms * 2. * ( I[lj0+2*(m+k)  ] + kms * I[lj0+2*(m-k)  ] );
		  T[(l0+2*m  )+(j0+2*k  )*na] = jms * 2. * (-I[lj0+2*(m+k)-1] + kms * I[lj0+2*(m-k)-1] );
		};
	      kms = -kms;
	    };
	    
	    
	};
    }


  else
    {


      for ( l=0; l <= lam; ++l )
	{
	  l0 = l*l;
	  T[l0] = I[l0];
	  for ( m=1; m < 2*l+1; ++m )
	    T[l0+m] = I[l0+m] * 2.;
	}

      jms = -1;
      for ( j=1; j <= lbm; ++j )
	{
	  j0 = j*j;
	  T[j0*na] = jms * I[j0];
	  for ( k=1; k < 2*j+1; ++k )
	    T[(j0+k)*na] = jms * I[j0+k] * 2.;
	  jms = -jms;
	}


      for ( l=1; l <= lam; ++l )
	{
	  l0 = l*l;

	  jms = -std::pow(-1,l);
	  int jmx = std::min(l-1,lbm);
	  for ( j=1; j <= jmx; ++j )
	    {
	      j0 = j*j;
	      for ( m=0; m < 2*l+1; ++m )
		for ( k=0; k < 2*j+1; ++k )
		  T[(l0+m)+(j0+k)*na] = jms * T[(j0+k)+(l0+m)*na];
	      jms = -jms;
	    }

	  //  we can build partially build j=l,lbm-1
	  //  from the l-1,j+1 block... we will just need to compute the m=l cases for all j

	  l1 = (l-1)*(l-1);
	  for ( j=l; j <= lbm-1; ++j )
	    {
	      jms = std::pow(-1,j);
	      j0  = j*j;
	      j1  = (j+1)*(j+1);
	      lj0 = (l+j)*(l+j);


	      //  l != m
	      for ( k=0; k < 2*j+1; ++k )
		for ( m=0; m < 2*l-1; ++m )
		  T[(l0+m)+(j0+k)*na] = -T[(l1+m)+(j1+k)*na];
		

	      //  m>0, k=0
	      //  (m-k > 0)
	      T[(l0+2*l-1)+j0*na] = jms * 2. * I[lj0+2*l-1];
	      T[(l0+2*l  )+j0*na] = jms * 2. * I[lj0+2*l  ];

	      // HERE *na

	      //  m=l, k=l
	      //  (m-k = 0)
	      //  (m = k)
	      mms = std::pow(-1,l);
	      T[(l0+2*l-1)+(j0+2*l-1)*na] = jms * 2. * ( I[lj0+4*l-1] + mms * I[lj0] );
	      T[(l0+2*l-1)+(j0+2*l  )*na] = jms * 2. * ( I[lj0+4*l  ] );
	      T[(l0+2*l  )+(j0+2*l-1)*na] = jms * 2. * ( I[lj0+4*l  ] );
	      T[(l0+2*l  )+(j0+2*l  )*na] = jms * 2. * (-I[lj0+4*l-1] + mms * I[lj0] );


	      //  m>0, k>0
	      //  m-k < 0
	      for ( k=l+1; k <= j; ++k )
		{
		  T[(l0+2*l-1)+(j0+2*k-1)*na] = jms * 2. * ( I[lj0+2*(l+k)-1] + mms * I[lj0-2*(l-k)-1] );
		  T[(l0+2*l-1)+(j0+2*k  )*na] = jms * 2. * ( I[lj0+2*(l+k)  ] + mms * I[lj0-2*(l-k)  ] );
		  T[(l0+2*l  )+(j0+2*k-1)*na] = jms * 2. * ( I[lj0+2*(l+k)  ] - mms * I[lj0-2*(l-k)  ] );
		  T[(l0+2*l  )+(j0+2*k  )*na] = jms * 2. * (-I[lj0+2*(l+k)-1] + mms * I[lj0-2*(l-k)-1] );
		}

	      //  m>0, k>0
	      //  m-k > 0
	      kms = -1;
	      for ( k=1; k <= l-1; ++k )
		{
		  T[(l0+2*l-1)+(j0+2*k-1)*na] = jms * 2. * ( I[lj0+2*(l+k)-1] + kms * I[lj0+2*(l-k)-1] );
		  T[(l0+2*l-1)+(j0+2*k  )*na] = jms * 2. * ( I[lj0+2*(l+k)  ] - kms * I[lj0+2*(l-k)  ] );
		  T[(l0+2*l  )+(j0+2*k-1)*na] = jms * 2. * ( I[lj0+2*(l+k)  ] + kms * I[lj0+2*(l-k)  ] );
		  T[(l0+2*l  )+(j0+2*k  )*na] = jms * 2. * (-I[lj0+2*(l+k)-1] + kms * I[lj0+2*(l-k)-1] );
		  kms = -kms;
		};
	    }

	  //  we have to build the entire j=lbm block

	  j = lbm;
	  jms = std::pow(-1,j);
	  j0  = j*j;
	  lj0 = (l+j)*(l+j);

	  //  m=0, k=0
	  T[l0+j0*na] = jms * I[lj0];

	  //  m=0, k>0
	  //  (m-k < 0)
	  for ( k=1; k <= j; ++k )
	    {
	      T[l0+(j0+2*k-1)*na] = jms * 2. * I[lj0+2*k-1];
	      T[l0+(j0+2*k  )*na] = jms * 2. * I[lj0+2*k  ];
	    };

	  //  m>0, k=0
	  //  (m-k > 0)
	  for ( m=1; m <= l; ++m )
	    {
	      T[(l0+2*m-1)+j0*na] = jms * 2. * I[lj0+2*m-1];
	      T[(l0+2*m  )+j0*na] = jms * 2. * I[lj0+2*m  ];
	    };

	  //  m>0, k>0
	  //  (m-k = 0)
	  //  (m = k)
	  mms = -1;
	  for ( m=1; m <= l; ++m )
	    {
	      T[(l0+2*m-1)+(j0+2*m-1)*na] = jms * 2. * ( I[lj0+4*m-1] + mms * I[lj0] );
	      T[(l0+2*m-1)+(j0+2*m  )*na] = jms * 2. * ( I[lj0+4*m  ] );
	      T[(l0+2*m  )+(j0+2*m-1)*na] = jms * 2. * ( I[lj0+4*m  ] );
	      T[(l0+2*m  )+(j0+2*m  )*na] = jms * 2. * (-I[lj0+4*m-1] + mms * I[lj0] );
	      mms = -mms;
	    }

	  //  m>0, k>0
	  //  m-k < 0
	  mms = -1;
	  for ( m=1; m <= l; ++m )
	    {
	      for ( k=m+1; k <= j; ++k )
		{
		  T[(l0+2*m-1)+(j0+2*k-1)*na] = jms * 2. * ( I[lj0+2*(m+k)-1] + mms * I[lj0-2*(m-k)-1] );
		  T[(l0+2*m-1)+(j0+2*k  )*na] = jms * 2. * ( I[lj0+2*(m+k)  ] + mms * I[lj0-2*(m-k)  ] );
		  T[(l0+2*m  )+(j0+2*k-1)*na] = jms * 2. * ( I[lj0+2*(m+k)  ] - mms * I[lj0-2*(m-k)  ] );
		  T[(l0+2*m  )+(j0+2*k  )*na] = jms * 2. * (-I[lj0+2*(m+k)-1] + mms * I[lj0-2*(m-k)-1] );
		};
	      mms = -mms;
	    };

	  //  m>0, k>0
	  //  m-k > 0
	  kms = -1;
	  for ( k=1; k <= j; ++k )
	    {
	      for ( m=k+1; m <= l; ++m )
		{
		  T[(l0+2*m-1)+(j0+2*k-1)*na] = jms * 2. * ( I[lj0+2*(m+k)-1] + kms * I[lj0+2*(m-k)-1] );
		  T[(l0+2*m-1)+(j0+2*k  )*na] = jms * 2. * ( I[lj0+2*(m+k)  ] - kms * I[lj0+2*(m-k)  ] );
		  T[(l0+2*m  )+(j0+2*k-1)*na] = jms * 2. * ( I[lj0+2*(m+k)  ] + kms * I[lj0+2*(m-k)  ] );
		  T[(l0+2*m  )+(j0+2*k  )*na] = jms * 2. * (-I[lj0+2*(m+k)-1] + kms * I[lj0+2*(m-k)-1] );
		};
	      kms = -kms;
	    };

	};
    };
 
}
*/




void ccdl::IlmInteraction
( int const lam, int const lbm, 
  double const *__restrict__ Rab,
  double *__restrict__ T )
{
  int l,m,j,k;
  int l0,j0,lj0,l1,j1,jc,js;
  int lms,jms,kms,mms;
  double tlms;


  if ( lam == 0 and lbm == 0 )
    {
      T[0] = 1. / std::sqrt( Rab[0]*Rab[0] + Rab[1]*Rab[1] + Rab[2]*Rab[2] );
      // double t = -T[0]*T[0]*T[0];
      // D[0] = t*Rab[0];
      // D[1] = t*Rab[1];
      // D[2] = t*Rab[2];
      return;
    };

  int const na = (lam+1)*(lam+1);
#ifdef __GNUG__
  double *__restrict__ I = (double *__restrict__) alloca ( sizeof(double)*((lam+lbm+2)*(lam+lbm+2)) );
  //double *__restrict__ dI = (double *__restrict__) alloca ( sizeof(double)*(3*(lam+lbm+1)*(lam+lbm+1)) );
#else
  std::vector<double> I( (lam+lbm+2)*(lam+lbm+2) , 0. );
  //std::vector<double> dI( 3*(lam+lbm+1)*(lam+lbm+1) , 0. );
#endif

  //double mrab[] = { Rab[0], Rab[1], -Rab[2] };
  ccdl::SolidHarm_Ilm(lam+lbm,Rab,I);
  //ccdl::SolidHarm_dIlm(lam+lbm,I,dI);





  if ( lam >= lbm )
    {

      lms = 1;
      for ( l=0; l <= lam; ++l )
	{
	  l0 = l*l;
	  T[l0] = lms * I[l0];
	  for ( m=1 ; m < 2*l+1; ++m )
	    T[l0+m] = lms * I[l0+m] * 2.;
	  lms=-lms;
	};

      for ( j=1; j <= lbm; ++j )
	{
	  j0 = j*j;
	  T[j0*na] = I[j0];
	  for ( k=1; k < 2*j+1; ++k )
	    T[(j0+k)*na] = I[j0+k] * 2.;
	};
	
      jms = 1;
      for ( j=1; j <= lbm; ++j )
	{
	  j0 = j*j;
	  jms = -jms;

	  lms = -jms;
	  int lmx = std::min(j-1,lam);
	  for ( l=1; l <= lmx; ++l )
	    {
	      l0 = l*l;
	      for ( m=0; m < 2*l+1; ++m )
		for ( k=0; k < 2*j+1; ++k )
		  T[(l0+m)+(j0+k)*na] = lms * T[(j0+k)+(l0+m)*na];
	      lms = -lms;
	    };
          

	  //  we can build partially build l=j,lam-1
	  //  from the l+1,j-1 block... we will just need to compute the k=j cases for all m


	  //jms = std::pow(-1,j);

	  lms = -jms;
	  j1 = (j-1)*(j-1);
	  for ( l=j; l <= lam-1; ++l )
	    {
	      lms = -lms;
	      tlms = 2.*lms;
	      l0 = l*l;
	      l1 = (l+1)*(l+1);
	      lj0 = (l+j)*(l+j);
		
	      //  k != j
	      for ( k=0; k < 2*j-1; ++k )
		for ( m=0; m < 2*l+1; ++m )
		  T[(l0+m)+(j0+k)*na] = -T[(l1+m)+(j1+k)*na];
		
	      //  k = j
		
	      jc = j0+2*j-1;
	      js = jc+1;
		
	      //  m=0, k=j
	      //  (m-k < 0)
	      T[l0+jc*na] = tlms * I[lj0+2*j-1];
	      T[l0+js*na] = tlms * I[lj0+2*j  ];
		
		
	      //  m=j, k=j
	      //  (m-k = 0)
	      //  (m = k)
	      T[(l0+2*j-1)+jc*na] = tlms * ( I[lj0+4*j-1] + jms * I[lj0] );
	      T[(l0+2*j  )+jc*na] = tlms * ( I[lj0+4*j  ] );
	      T[(l0+2*j-1)+js*na] = tlms * ( I[lj0+4*j  ] );
	      T[(l0+2*j  )+js*na] = tlms * (-I[lj0+4*j-1] + jms * I[lj0] );
		
		
	      //  m>0, k>0
	      //  m-k < 0
	      mms = -1;
	      int mmx = std::min(l,j-1);
	      for ( m=1; m <= mmx; ++m )
		{
		  T[(l0+2*m-1)+jc*na] = tlms * ( I[lj0+2*(m+j)-1] + mms * I[lj0-2*(m-j)-1] );
		  T[(l0+2*m  )+jc*na] = tlms * ( I[lj0+2*(m+j)  ] - mms * I[lj0-2*(m-j)  ] );
		  T[(l0+2*m-1)+js*na] = tlms * ( I[lj0+2*(m+j)  ] + mms * I[lj0-2*(m-j)  ] );
		  T[(l0+2*m  )+js*na] = tlms * (-I[lj0+2*(m+j)-1] + mms * I[lj0-2*(m-j)-1] );
		  mms = -mms;
		}
		
	      //  m>0, k>0
	      //  m-k > 0
	      for ( m=j+1; m <= l; ++m )
		{
		  T[(l0+2*m-1)+jc*na] = tlms * ( I[lj0+2*(m+j)-1] + jms * I[lj0+2*(m-j)-1] );
		  T[(l0+2*m  )+jc*na] = tlms * ( I[lj0+2*(m+j)  ] + jms * I[lj0+2*(m-j)  ] );
		  T[(l0+2*m-1)+js*na] = tlms * ( I[lj0+2*(m+j)  ] - jms * I[lj0+2*(m-j)  ] );
		  T[(l0+2*m  )+js*na] = tlms * (-I[lj0+2*(m+j)-1] + jms * I[lj0+2*(m-j)-1] );
		}
	    };
	    
	    
	    
	  //  we have to build the entire l=lam block
	    

	  l = lam;
	  lms = -lms;
	  tlms = 2.*lms;
	    
	  l0  = l*l;
	  lj0 = (l+j)*(l+j);
	    
	  //  m=0, k=0
	  T[l0+j0*na] = lms * I[lj0];
	    
	  //  m=0, k>0
	  //  (m-k < 0)
	  for ( k=1; k <= j; ++k )
	    {
	      T[l0+(j0+2*k-1)*na] = tlms * I[lj0+2*k-1];
	      T[l0+(j0+2*k  )*na] = tlms * I[lj0+2*k];
	    };
	  //  m>0, k=0
	  //  (m-k > 0)
	  for ( m=1; m <= l; ++m )
	    {
	      T[(l0+2*m-1)+j0*na] = tlms * I[lj0+2*m-1];
	      T[(l0+2*m  )+j0*na] = tlms * I[lj0+2*m  ];
	    }
	    
	  //  m>0, k>0
	  //  (m-k = 0)
	  //  (m = k)
	  mms = -1;
	  for ( m=1; m <= j; ++m )
	    {
	      T[(l0+2*m-1)+(j0+2*m-1)*na] = tlms * ( I[lj0+4*m-1] + mms * I[lj0] );
	      T[(l0+2*m  )+(j0+2*m-1)*na] = tlms * ( I[lj0+4*m  ]);
	      T[(l0+2*m-1)+(j0+2*m  )*na] = tlms * ( I[lj0+4*m  ]);
	      T[(l0+2*m  )+(j0+2*m  )*na] = tlms * (-I[lj0+4*m-1] + mms * I[lj0] );
	      mms = -mms;
	    };
	    
	    
	  //  m>0, k>0
	  //  m-k < 0
	  mms = -1;
	  for ( m=1; m <= l; ++m )
	    {
	      for ( k=m+1; k <= j; ++k )
		{
		  T[(l0+2*m-1)+(j0+2*k-1)*na] = tlms * ( I[lj0+2*(m+k)-1] + mms * I[lj0-2*(m-k)-1] );
		  T[(l0+2*m  )+(j0+2*k-1)*na] = tlms * ( I[lj0+2*(m+k)  ] - mms * I[lj0-2*(m-k)  ] );
		  T[(l0+2*m-1)+(j0+2*k  )*na] = tlms * ( I[lj0+2*(m+k)  ] + mms * I[lj0-2*(m-k)  ] );
		  T[(l0+2*m  )+(j0+2*k  )*na] = tlms * (-I[lj0+2*(m+k)-1] + mms * I[lj0-2*(m-k)-1] );
		};
	      mms = -mms;
	    }
	    
	  //  m>0, k>0
	  //  m-k > 0
	  kms = -1;
	  for ( k=1; k <= j; ++k )
	    {
	      for ( m=k+1; m <= l; ++m )
		{
		  T[(l0+2*m-1)+(j0+2*k-1)*na] = tlms * ( I[lj0+2*(m+k)-1] + kms * I[lj0+2*(m-k)-1] );
		  T[(l0+2*m  )+(j0+2*k-1)*na] = tlms * ( I[lj0+2*(m+k)  ] + kms * I[lj0+2*(m-k)  ] );
		  T[(l0+2*m-1)+(j0+2*k  )*na] = tlms * ( I[lj0+2*(m+k)  ] - kms * I[lj0+2*(m-k)  ] );
		  T[(l0+2*m  )+(j0+2*k  )*na] = tlms * (-I[lj0+2*(m+k)-1] + kms * I[lj0+2*(m-k)-1] );
		};
	      kms = -kms;
	    };
	    
	    
	};
    }

  /* */

  else
    {

      lms = 1;
      for ( l=0; l <= lam; ++l )
	{
	  l0 = l*l;
	  T[l0] = lms * I[l0];
	  for ( m=1; m < 2*l+1; ++m )
	    T[l0+m] = lms * I[l0+m] * 2.;
	  lms = -lms;
	}

      for ( j=1; j <= lbm; ++j )
	{
	  j0 = j*j;
	  T[j0*na] = I[j0];
	  for ( k=1; k < 2*j+1; ++k )
	    T[(j0+k)*na] = I[j0+k] * 2.;
	}


      lms = 1;
      for ( l=1; l <= lam; ++l )
	{
	  l0 = l*l;
	  lms = -lms;
	  tlms = 2.*lms;

	  jms = -lms;
	  int jmx = std::min(l-1,lbm);
	  for ( j=1; j <= jmx; ++j )
	    {
	      j0 = j*j;
	      for ( m=0; m < 2*l+1; ++m )
		for ( k=0; k < 2*j+1; ++k )
		  T[(l0+m)+(j0+k)*na] = jms * T[(j0+k)+(l0+m)*na];
	      jms = -jms;
	    }

	  //  we can build partially build j=l,lbm-1
	  //  from the l-1,j+1 block... we will just need to compute the m=l cases for all j

	  jms = -lms;
	  l1 = (l-1)*(l-1);
	  for ( j=l; j <= lbm-1; ++j )
	    {
	      jms = -jms;
	      j0  = j*j;
	      j1  = (j+1)*(j+1);
	      lj0 = (l+j)*(l+j);


	      //  l != m
	      for ( k=0; k < 2*j+1; ++k )
		for ( m=0; m < 2*l-1; ++m )
		  T[(l0+m)+(j0+k)*na] = -T[(l1+m)+(j1+k)*na];
		

	      //  m>0, k=0
	      //  (m-k > 0)
	      T[(l0+2*l-1)+j0*na] = tlms * I[lj0+2*l-1];
	      T[(l0+2*l  )+j0*na] = tlms * I[lj0+2*l  ];

	      // HERE *na

	      //  m=l, k=l
	      //  (m-k = 0)
	      //  (m = k)
	      mms = lms; // std::pow(-1,l);
	      T[(l0+2*l-1)+(j0+2*l-1)*na] = tlms * ( I[lj0+4*l-1] + mms * I[lj0] );
	      T[(l0+2*l  )+(j0+2*l-1)*na] = tlms * ( I[lj0+4*l  ] );
	      T[(l0+2*l-1)+(j0+2*l  )*na] = tlms * ( I[lj0+4*l  ] );
	      T[(l0+2*l  )+(j0+2*l  )*na] = tlms * (-I[lj0+4*l-1] + mms * I[lj0] );


	      //  m>0, k>0
	      //  m-k < 0
	      for ( k=l+1; k <= j; ++k )
		{
		  T[(l0+2*l-1)+(j0+2*k-1)*na] = tlms * ( I[lj0+2*(l+k)-1] + mms * I[lj0-2*(l-k)-1] );
		  T[(l0+2*l  )+(j0+2*k-1)*na] = tlms * ( I[lj0+2*(l+k)  ] - mms * I[lj0-2*(l-k)  ] );
		  T[(l0+2*l-1)+(j0+2*k  )*na] = tlms * ( I[lj0+2*(l+k)  ] + mms * I[lj0-2*(l-k)  ] );
		  T[(l0+2*l  )+(j0+2*k  )*na] = tlms * (-I[lj0+2*(l+k)-1] + mms * I[lj0-2*(l-k)-1] );
		}

	      //  m>0, k>0
	      //  m-k > 0
	      kms = -1;
	      for ( k=1; k <= l-1; ++k )
		{
		  T[(l0+2*l-1)+(j0+2*k-1)*na] = tlms * ( I[lj0+2*(l+k)-1] + kms * I[lj0+2*(l-k)-1] );
		  T[(l0+2*l  )+(j0+2*k-1)*na] = tlms * ( I[lj0+2*(l+k)  ] + kms * I[lj0+2*(l-k)  ] );
		  T[(l0+2*l-1)+(j0+2*k  )*na] = tlms * ( I[lj0+2*(l+k)  ] - kms * I[lj0+2*(l-k)  ] );
		  T[(l0+2*l  )+(j0+2*k  )*na] = tlms * (-I[lj0+2*(l+k)-1] + kms * I[lj0+2*(l-k)-1] );
		  kms = -kms;
		};
	    }

	  //  we have to build the entire j=lbm block

	  j = lbm;
	  jms = -jms;
	  j0  = j*j;
	  lj0 = (l+j)*(l+j);

	  //  m=0, k=0
	  T[l0+j0*na] = lms * I[lj0];

	  //  m=0, k>0
	  //  (m-k < 0)
	  for ( k=1; k <= j; ++k )
	    {
	      T[l0+(j0+2*k-1)*na] = tlms * I[lj0+2*k-1];
	      T[l0+(j0+2*k  )*na] = tlms * I[lj0+2*k  ];
	    };

	  //  m>0, k=0
	  //  (m-k > 0)
	  for ( m=1; m <= l; ++m )
	    {
	      T[(l0+2*m-1)+j0*na] = tlms * I[lj0+2*m-1];
	      T[(l0+2*m  )+j0*na] = tlms * I[lj0+2*m  ];
	    };

	  //  m>0, k>0
	  //  (m-k = 0)
	  //  (m = k)
	  mms = -1;
	  for ( m=1; m <= l; ++m )
	    {
	      T[(l0+2*m-1)+(j0+2*m-1)*na] = tlms * ( I[lj0+4*m-1] + mms * I[lj0] );
	      T[(l0+2*m  )+(j0+2*m-1)*na] = tlms * ( I[lj0+4*m  ] );
	      T[(l0+2*m-1)+(j0+2*m  )*na] = tlms * ( I[lj0+4*m  ] );
	      T[(l0+2*m  )+(j0+2*m  )*na] = tlms * (-I[lj0+4*m-1] + mms * I[lj0] );
	      mms = -mms;
	    }

	  //  m>0, k>0
	  //  m-k < 0
	  mms = -1;
	  for ( m=1; m <= l; ++m )
	    {
	      for ( k=m+1; k <= j; ++k )
		{
		  T[(l0+2*m-1)+(j0+2*k-1)*na] = tlms * ( I[lj0+2*(m+k)-1] + mms * I[lj0-2*(m-k)-1] );
		  T[(l0+2*m  )+(j0+2*k-1)*na] = tlms * ( I[lj0+2*(m+k)  ] - mms * I[lj0-2*(m-k)  ] );
		  T[(l0+2*m-1)+(j0+2*k  )*na] = tlms * ( I[lj0+2*(m+k)  ] + mms * I[lj0-2*(m-k)  ] );
		  T[(l0+2*m  )+(j0+2*k  )*na] = tlms * (-I[lj0+2*(m+k)-1] + mms * I[lj0-2*(m-k)-1] );
		};
	      mms = -mms;
	    };

	  //  m>0, k>0
	  //  m-k > 0
	  kms = -1;
	  for ( k=1; k <= j; ++k )
	    {
	      for ( m=k+1; m <= l; ++m )
		{
		  T[(l0+2*m-1)+(j0+2*k-1)*na] = tlms * ( I[lj0+2*(m+k)-1] + kms * I[lj0+2*(m-k)-1] );
		  T[(l0+2*m  )+(j0+2*k-1)*na] = tlms * ( I[lj0+2*(m+k)  ] + kms * I[lj0+2*(m-k)  ] );
		  T[(l0+2*m-1)+(j0+2*k  )*na] = tlms * ( I[lj0+2*(m+k)  ] - kms * I[lj0+2*(m-k)  ] );
		  T[(l0+2*m  )+(j0+2*k  )*na] = tlms * (-I[lj0+2*(m+k)-1] + kms * I[lj0+2*(m-k)-1] );
		};
	      kms = -kms;
	    };

	};
    };
 
}





void ccdl::IlmInteractionGrd
( int const lam, int const lbm, 
  double const *__restrict__ Rab,
  double *__restrict__ T,
  double *__restrict__ D )
{
  int l,l1,  m,jk,lm,  na,nb,na1;

  na1 = (lam+2)*(lam+2);
  na = (lam+1)*(lam+1);
  nb = (lbm+1)*(lbm+1);
  //std::vector<double> V( na1*nb , 0. );
  double *__restrict__ V = (double *)malloc(sizeof(double)*(na1*nb));

  ccdl::IlmInteraction(lam+1,lbm,Rab,V);

  for ( jk=0; jk < nb; ++jk )
    for ( lm=0; lm < na; ++lm )
      T[lm+jk*na] = V[lm+jk*na1];

  --D;
  //D = 0.;

  if ( lam == 0 )
    {
      double const *__restrict__ p = V;
      for ( jk=0; jk < nb; ++jk, p += na1 )
	{
	  *(++D) = -0.5 * p[2];
	  *(++D) = -0.5 * p[3];
	  *(++D) = p[1];
	}
    }
  else if ( lam == 1 )
    {
      double const *__restrict__ p = V;
      for ( jk=0; jk < nb; ++jk, p += na1 )
	{
          *(++D) = -0.5 * p[2];
          *(++D) = -0.5 * p[3];
          *(++D) = p[1];

          *(++D) = -0.5 * p[5];
          *(++D) = -0.5 * p[6];
          *(++D) = p[4];

          *(++D) = -0.5 * ( p[7] - 2. * p[4] );
          *(++D) = -0.5 * ( p[8] );
          *(++D) = p[5];

          *(++D) = -0.5 * ( p[8] );
          *(++D) =  0.5 * ( p[7] + 2. * p[4] );
          *(++D) = p[6];
	}
    }
  else
    {
      double const *__restrict__ p = V;
      for ( jk=0; jk < nb; ++jk, p += na1 )
	{
          *(++D) = -0.5 * p[2];
          *(++D) = -0.5 * p[3];
          *(++D) = p[1];

          *(++D) = -0.5 * p[5];
          *(++D) = -0.5 * p[6];
          *(++D) = p[4];

          *(++D) = -0.5 * ( p[7] - 2. * p[4] );
          *(++D) = -0.5 * ( p[8] );
          *(++D) = p[5];

          *(++D) = -0.5 * ( p[8] );
          *(++D) =  0.5 * ( p[7] + 2. * p[4] );
          *(++D) = p[6];

          for ( l=2; l <= lam; ++l )
	    {
	      //l0 = l*l;
	      l1 = (l+1)*(l+1);

	      // m=0;
	      *(++D) = -0.5 * p[l1+1];
	      *(++D) = -0.5 * p[l1+2];
	      *(++D) = p[l1];

	      // m=1;
	      *(++D) = -0.5 * ( p[l1+3] - 2. * p[l1] );
	      *(++D) = -0.5 * ( p[l1+4] );
	      *(++D) = p[l1+1];
	      *(++D) = -0.5 * ( p[l1+4] );
	      *(++D) =  0.5 * ( p[l1+3] + 2. * p[l1] );
	      *(++D) = p[l1+2];

	      // m>1;
	      for ( m=2; m <= l; ++m )
		{
		  *(++D) = -0.5 * ( p[l1+2*(m+1)-1] - p[l1+2*(m-1)-1] );
		  *(++D) = -0.5 * ( p[l1+2*(m+1)  ] + p[l1+2*(m-1)  ] );
		  *(++D) = p[l1+2*m-1];

		  *(++D) = -0.5 * ( p[l1+2*(m+1)  ] - p[l1+2*(m-1)  ] );
		  *(++D) =  0.5 * ( p[l1+2*(m+1)-1] + p[l1+2*(m-1)-1] );
		  *(++D) = p[l1+2*m  ];
		}
	    }
	};
    }

  free(V);
}
