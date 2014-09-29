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
  ccdl::SolidHarmIlm_v5(lam+lbm,mrab,I);





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

  /* */

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
