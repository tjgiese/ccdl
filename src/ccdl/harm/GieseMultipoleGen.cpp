#include <algorithm>
#include "harm.hpp"
#include "../constants.hpp"
#include "GieseMultipoleGen.hpp"
#include "GieseBoysFcn.hpp"
#include "../fmt.hpp"

#ifndef alloca
#define alloca __builtin_alloca
#endif

#ifdef __GNUG__
#define ALLOC(n) (double *__restrict__) alloca ( sizeof(double)*(n) )
#define FREE(n)
#else
#define ALLOC(n) (double *__restrict__) malloc ( sizeof(double)*(n) )
#define FREE(n) free(n)
#endif




 

namespace ccdl
{

  void AuxExpInt_Helper
  ( int const la, int const lb, 
    int const wna,
    double const *__restrict__ W,
    double const *__restrict__ O,
    double *__restrict__ T )
  {
    int const NA = (la+1)*(la+1);
    int const NB = (lb+1)*(lb+1);
    int const nt   = NA*NB;
    for ( int i=0; i<nt; ++i ) T[i] = 0.;

    int lmax,lmin;
    if ( la > lb )
      {
	lmax = la;
	lmin = lb;
      }
    else
      {
	lmax = lb;
	lmin = la;
      };
    //int const nmax = (lmax+1)*(lmax+1);
    double *__restrict__ TL = ALLOC(lmax+1);
    TL[0] = 1.;
    for ( int i=1; i <= lmax; ++i )
      TL[i] = (i-0.5) * TL[i-1];

    // version 3
    double ms = -1;
    for ( int Lb=0; Lb <= lb; ++Lb )
      {
	ms = -ms;
	int ob = Lb*Lb;
	int nb = ob + 2*Lb+1;

	int lstart = 0;
	if ( Lb <= la ) // this block is related by symmetry
	  {
	    lstart = Lb;
	    for ( int La=0; La < Lb; ++La )
	      {
		int oa = La*La;
		int na = oa + 2*La+1;
		double mm = (La+Lb) % 2 ? -1. : 1.;
		for ( int mb=ob; mb < nb; ++mb )
		  for ( int ma=oa; ma < na; ++ma )
		    T[ma+mb*NA] = mm *  T[mb+ma*NA];
	      }
	  };

	for ( int La=lstart; La <= la; ++La ) // compute the missing blocks in this col
	  {
	    int oa = La*La;
	    int na = oa + 2*La+1;
	    lmin = std::min( La, Lb );
	    for ( int j=0; j <= lmin; ++j )
	      {
		int nj = (j+1)*(j+1);
		double x = ms * O[La+Lb-j] * TL[j] / ( TL[La] * TL[Lb] );

		for ( int jk=j*j; jk < nj; ++jk )
		  for ( int mb=ob; mb < nb; ++mb )
		    for ( int ma=oa; ma < na; ++ma )
		      T[ma+mb*NA] += x *  W[ma+jk*wna] * W[mb+jk*wna];
	      }
	  }
      }
    FREE(TL);
  }

  void AuxExpGrd_Helper
  ( int const LaMax, int const LbMax, 
    double const *__restrict__ Rab,
    double const *__restrict__ X1,
    double *__restrict__ dX )
  {
    int const na = (LaMax+1)*(LaMax+1);
    //int const nb = (LbMax+1)*(LbMax+1);
    int lmax = std::max(LaMax,LbMax);
    int nmax = (lmax+1)*(lmax+1);
    
    double *__restrict__ Alm = ALLOC(nmax);
    ccdl::CptAlm(lmax,Alm);
    
    double tcrd[3] = { 2.*Rab[0], 2.*Rab[1], 2.*Rab[2] };
    int i;

    for ( int lb=0; lb<=LbMax; ++lb )
      {
	int lb0 = lb*lb;
	int lb1 = (lb-1)*(lb-1);
	double oo2lb = 1./(2*lb-1);
	for ( int la=0; la<=LaMax; ++la )
	  {
	    int la0 = la*la;
	    int la1 = (la-1)*(la-1);
	    double oo2la = 1./(2*la-1);


	    // symmetry
	    if ( lb > la and lb <= LaMax )
	      {
	    	double s = ((la+lb)%2) ? -1. : 1.;
	    	int nma = la0+2*la+1;
	    	int nmb = lb0+2*lb+1;
	    	int j;
	    	for ( int ib=lb0; ib<nmb; ++ib )
	    	  for ( int ia=la0; ia<nma; ++ia )
	    	    {
	    	      i = (ia+ib*na)*3;
	    	      j = (ib+ia*na)*3;
	    	      dX[0+i] = s*dX[0+j];
	    	      dX[1+i] = s*dX[1+j];
	    	      dX[2+i] = s*dX[2+j];
	    	    }
	    	continue;
	      };


	    int ib = lb0-1;
	    for ( int mb=0; mb<=lb; ++mb )
	      for ( int bcs=0; bcs<2; ++bcs )
		{
		  if ( mb == 0 and bcs > 0 ) continue;
		  ++ib;
		  double bf = oo2lb*Alm[ib];


		  i = lb1 + 2*mb;
		  int ib0  = i   + (    (mb>0)?(bcs-1):0);
		  int ibp1 = i+2 + (((mb+1)>0)?(bcs-1):0);
		  int ibm1 = i-2 + (((mb-1)>0)?(bcs-1):0);


		  int ia = la0-1;
		  for ( int ma=0; ma<=la; ++ma )
		    for ( int acs=0; acs<2; ++acs )
		      {
			if ( ma == 0 and acs > 0 ) continue;
			++ia;
			double af = oo2la*Alm[ia];
			
			i = la1 + 2*ma; // sine
			int ia0  = i   + (    (ma>0)?(acs-1):0);
			int iap1 = i+2 + (((ma+1)>0)?(acs-1):0);
			int iam1 = i-2 + (((ma-1)>0)?(acs-1):0);



			// TERM 1
			i = ia+(ib)*na;
			double dx = tcrd[0]*X1[i];
			double dy = tcrd[1]*X1[i];
			double dz = tcrd[2]*X1[i];


			//==============================================
			// d/dz
			//==============================================
			
			// TERM 2
			if ( mb < lb )
			  dz -= 2.*bf/Alm[ib0] * X1[ia+ib0*na];
			
			// TERM 3
			if ( ma < la ) 
			  dz += 2.*af/Alm[ia0] * X1[ia0+ib*na];


			// TERM 2
			if ( ma+1 < la ) 
			  {
			    //  get the index of (la-1,ma+1)
			    // cosine yields a cosine; sine yields a sine
			    dx += af/Alm[iap1] * X1[iap1+(ib)*na];
			    // cosine yields a sine; sine yields a cosine
			    if ( acs == 0 ) 
			      {
			    	dy += af/Alm[(iap1+1)] * X1[(iap1+1)+(ib)*na];
			      }
			    else 
			      {
			    	dy -= af/Alm[(iap1-1)] * X1[(iap1-1)+(ib)*na];
			      }
			  }

			// TERM 3
			if ( ma == 0 )
			  {
			    if ( la > 1 )
			      {
				dx += af/Alm[iap1] * X1[iap1+(ib)*na];
				dy += af/Alm[iap1] * X1[(iap1+1)+(ib)*na];
			      };
			  }
			else if ( ma == 1 )
			  {
			    if ( acs == 0 )
			      {
				dx -= af/Alm[iam1] * X1[iam1+(ib)*na];
			      }
			    else
			      {
				i = la1 + 2*(ma-1);
				dy -= af/Alm[i] * X1[i+(ib)*na];
			      };
			  }
			else
			  {
			    dx -= af/Alm[iam1] * X1[iam1+(ib)*na];
			    if  ( acs == 0 )
			      {
				dy += af/Alm[iam1+1] * X1[(iam1+1)+(ib)*na];
			      }
			    else
			      {
				i = la1 + 2*(ma-1)-1;
				dy -= af/Alm[i] * X1[i+(ib)*na];
			      };
			  };


			// TERM 4
			if ( mb+1 < lb ) 
			  {
			    //  get the index of (lb-1,mb+1)
			    // cosine yields a cosine; sine yields a sine
			    dx -= bf/Alm[ibp1] * X1[ia+(ibp1)*na];

			    //  get the index of (la-1,ma+1)
			    // cosine yields a sine; sine yields a cosine
			    if ( bcs == 0 ) 
			      {
				dy -= bf/Alm[ibp1] * X1[ia+(ibp1+1)*na];
			      }
			    else 
			      {
				dy += bf/Alm[ibp1] * X1[ia+(ibp1-1)*na];
			      }
			  }
			 
			// TERM 5
			if ( mb == 0 )
			  {
			    if ( lb > 1 )
			      {
				dx -= bf/Alm[ibp1] * X1[ia+(ibp1)*na];
				dy -= bf/Alm[ibp1] * X1[ia+(ibp1+1)*na];
			      };
			  }
			else if ( mb == 1 )
			  {
			    if ( bcs == 0 )
			      {
				dx += bf/Alm[ibm1] * X1[ia+(ibm1)*na];
			      }
			    else
			      {
				i = lb1 + 2*(mb-1);
				dy += bf/Alm[i] * X1[ia+(i)*na];
			      };
			  }
			else
			  {
			    dx += bf/Alm[ibm1] * X1[ia+(ibm1)*na];
			    if  ( bcs == 0 )
			      {
				dy -= bf/Alm[ibm1+1] * X1[ia+(ibm1+1)*na];
			      }
			    else
			      {
				i = lb1 + 2*(mb-1)-1;
				dy += bf/Alm[i] * X1[ia+(i)*na];
			      };
			  };

			 
			i = (ia+ib*na)*3;
			dX[0+i] = dx;
			dX[1+i] = dy;
			dX[2+i] = dz;
			
		      };
		};
	  };
      };
  }






  void AuxExpPotGrd_Helper
  ( int const LaMax, int const LbMax, 
    double const *__restrict__ Rab,
    double const *__restrict__ X1,
    double const * qa,
    double const * qb,
    double *__restrict__ g )
  {
    g[0] = 0.;
    g[1] = 0.;
    g[2] = 0.;

    int const na = (LaMax+1)*(LaMax+1);
    //int const nb = (LbMax+1)*(LbMax+1);
    int lmax = std::max(LaMax,LbMax);
    int nmax = (lmax+1)*(lmax+1);
    
    double *__restrict__ Alm = ALLOC(nmax);
    ccdl::CptAlm(lmax,Alm);
    
    double tcrd[3] = { 2.*Rab[0], 2.*Rab[1], 2.*Rab[2] };
    int i;

    for ( int lb=0; lb<=LbMax; ++lb )
      {
	int lb0 = lb*lb;
	int lb1 = (lb-1)*(lb-1);
	double oo2lb = 1./(2*lb-1);
	for ( int la=0; la<=LaMax; ++la )
	  {
	    int la0 = la*la;
	    int la1 = (la-1)*(la-1);
	    double oo2la = 1./(2*la-1);


	    // symmetry
	    if ( lb > la and lb <= LaMax ) continue;
	    bool sym = false;
	    double s = 1.;
	    if ( la > lb and la <= LbMax ) 
	      {
		sym = true;
		s = ((la+lb)%2) ? -1. : 1.;
	      };



	    int ib = lb0-1;
	    for ( int mb=0; mb<=lb; ++mb )
	      for ( int bcs=0; bcs<2; ++bcs )
		{
		  if ( mb == 0 and bcs > 0 ) continue;
		  ++ib;
		  double bf = oo2lb*Alm[ib];


		  i = lb1 + 2*mb;
		  int ib0  = i   + (    (mb>0)?(bcs-1):0);
		  int ibp1 = i+2 + (((mb+1)>0)?(bcs-1):0);
		  int ibm1 = i-2 + (((mb-1)>0)?(bcs-1):0);


		  int ia = la0-1;
		  for ( int ma=0; ma<=la; ++ma )
		    for ( int acs=0; acs<2; ++acs )
		      {
			if ( ma == 0 and acs > 0 ) continue;
			++ia;
			double af = oo2la*Alm[ia];
			
			i = la1 + 2*ma; // sine
			int ia0  = i   + (    (ma>0)?(acs-1):0);
			int iap1 = i+2 + (((ma+1)>0)?(acs-1):0);
			int iam1 = i-2 + (((ma-1)>0)?(acs-1):0);



			// TERM 1
			i = ia+(ib)*na;
			double dx = tcrd[0]*X1[i];
			double dy = tcrd[1]*X1[i];
			double dz = tcrd[2]*X1[i];


			//==============================================
			// d/dz
			//==============================================
			
			// TERM 2
			if ( mb < lb )
			  dz -= 2.*bf/Alm[ib0] * X1[ia+ib0*na];
			
			// TERM 3
			if ( ma < la ) 
			  dz += 2.*af/Alm[ia0] * X1[ia0+ib*na];


			// TERM 2
			if ( ma+1 < la ) 
			  {
			    //  get the index of (la-1,ma+1)
			    // cosine yields a cosine; sine yields a sine
			    dx += af/Alm[iap1] * X1[iap1+(ib)*na];
			    // cosine yields a sine; sine yields a cosine
			    if ( acs == 0 ) 
			      {
			    	dy += af/Alm[(iap1+1)] * X1[(iap1+1)+(ib)*na];
			      }
			    else 
			      {
			    	dy -= af/Alm[(iap1-1)] * X1[(iap1-1)+(ib)*na];
			      }
			  }

			// TERM 3
			if ( ma == 0 )
			  {
			    if ( la > 1 )
			      {
				dx += af/Alm[iap1] * X1[iap1+(ib)*na];
				dy += af/Alm[iap1] * X1[(iap1+1)+(ib)*na];
			      };
			  }
			else if ( ma == 1 )
			  {
			    if ( acs == 0 )
			      {
				dx -= af/Alm[iam1] * X1[iam1+(ib)*na];
			      }
			    else
			      {
				i = la1 + 2*(ma-1);
				dy -= af/Alm[i] * X1[i+(ib)*na];
			      };
			  }
			else
			  {
			    dx -= af/Alm[iam1] * X1[iam1+(ib)*na];
			    if  ( acs == 0 )
			      {
				dy += af/Alm[iam1+1] * X1[(iam1+1)+(ib)*na];
			      }
			    else
			      {
				i = la1 + 2*(ma-1)-1;
				dy -= af/Alm[i] * X1[i+(ib)*na];
			      };
			  };


			// TERM 4
			if ( mb+1 < lb ) 
			  {
			    //  get the index of (lb-1,mb+1)
			    // cosine yields a cosine; sine yields a sine
			    dx -= bf/Alm[ibp1] * X1[ia+(ibp1)*na];

			    //  get the index of (la-1,ma+1)
			    // cosine yields a sine; sine yields a cosine
			    if ( bcs == 0 ) 
			      {
				dy -= bf/Alm[ibp1] * X1[ia+(ibp1+1)*na];
			      }
			    else 
			      {
				dy += bf/Alm[ibp1] * X1[ia+(ibp1-1)*na];
			      }
			  }
			 
			// TERM 5
			if ( mb == 0 )
			  {
			    if ( lb > 1 )
			      {
				dx -= bf/Alm[ibp1] * X1[ia+(ibp1)*na];
				dy -= bf/Alm[ibp1] * X1[ia+(ibp1+1)*na];
			      };
			  }
			else if ( mb == 1 )
			  {
			    if ( bcs == 0 )
			      {
				dx += bf/Alm[ibm1] * X1[ia+(ibm1)*na];
			      }
			    else
			      {
				i = lb1 + 2*(mb-1);
				dy += bf/Alm[i] * X1[ia+(i)*na];
			      };
			  }
			else
			  {
			    dx += bf/Alm[ibm1] * X1[ia+(ibm1)*na];
			    if  ( bcs == 0 )
			      {
				dy -= bf/Alm[ibm1+1] * X1[ia+(ibm1+1)*na];
			      }
			    else
			      {
				i = lb1 + 2*(mb-1)-1;
				dy += bf/Alm[i] * X1[ia+(i)*na];
			      };
			  };

			 
			if ( sym )
			  {
			    double qq = qa[ia]*qb[ib] + s*qa[ib]*qb[ia];
			    g[0] += qq*dx;
			    g[1] += qq*dy;
			    g[2] += qq*dz;
			  }
			else
			  {
			    double qq = qa[ia]*qb[ib];
			    g[0] += qq*dx;
			    g[1] += qq*dy;
			    g[2] += qq*dz;
			  };

		      };
		};
	  };
      };
  }



}



void ccdl::AuxExpInt
( int const la, int const lb, 
  double const *__restrict__ crd,
  double const *__restrict__ O,
  double *__restrict__ T )
{
  int const NA = (la+1)*(la+1);
  int const NB = (lb+1)*(lb+1);
  int const nt   = NA*NB;
  for ( int i=0; i<nt; ++i ) T[i] = 0.;

  int lmax,lmin;
  if ( la > lb )
    {
      lmax = la;
      lmin = lb;
    }
  else
    {
      lmax = lb;
      lmin = la;
    };
  int const nmax = (lmax+1)*(lmax+1);

  double *__restrict__ W  = ALLOC(nt);
  double *__restrict__ TL = ALLOC(lmax+1);
  ccdl::ClmTranslation( lmax, lmin, crd, W );
  TL[0] = 1.;
  for ( int i=1; i <= lmax; ++i )
    TL[i] = (i-0.5) * TL[i-1];

  // version 3
  double ms = -1;
  for ( int Lb=0; Lb <= lb; ++Lb )
    {
      ms = -ms;
      int ob = Lb*Lb;
      int nb = ob + 2*Lb+1;

      int lstart = 0;
      if ( Lb <= la ) // this block is related by symmetry
	{
	  lstart = Lb;
	  for ( int La=0; La < Lb; ++La )
	    {
	      int oa = La*La;
	      int na = oa + 2*La+1;
	      double mm = (La+Lb) % 2 ? -1. : 1.;
	      for ( int mb=ob; mb < nb; ++mb )
		for ( int ma=oa; ma < na; ++ma )
		  T[ma+mb*NA] = mm *  T[mb+ma*NA];
	    }
	};

      for ( int La=lstart; La <= la; ++La ) // compute the missing blocks in this col
	{
          int oa = La*La;
	  int na = oa + 2*La+1;
	  lmin = std::min( La, Lb );
          for ( int j=0; j <= lmin; ++j )
	    {
	      int nj = (j+1)*(j+1);
	      double x = ms * O[La+Lb-j] * TL[j] / ( TL[La] * TL[Lb] );

	      for ( int jk=j*j; jk < nj; ++jk )
		for ( int mb=ob; mb < nb; ++mb )
		  for ( int ma=oa; ma < na; ++ma )
		    T[ma+mb*NA] += x *  W[ma+jk*nmax] * W[mb+jk*nmax];
	    }
	}
    }
  FREE(W);
  FREE(TL);
}






void ccdl::AuxExpGrd
( int const la, int const lb, 
  double const *__restrict__ crd,
  double const *__restrict__ O,
  double *__restrict__ T,
  double *__restrict__ dT )
{
  int na = (la+1)*(la+1);
  int nb = (lb+1)*(lb+1);
  int lmax,lmin,nmax;
  if ( na > nb )
    {
      lmax = la;
      lmin = lb;
      nmax = na;
    }
  else
    {
      lmax = lb;
      lmin = la;
      nmax = nb;
    };
  double *__restrict__ W = ALLOC(na*nb);
  ccdl::ClmTranslation( lmax, lmin, crd, W );
  ccdl::AuxExpInt_Helper(la,lb,nmax,W,O+1,T);
  ccdl::AuxExpGrd_Helper(la,lb,crd,T,dT);
  ccdl::AuxExpInt_Helper(la,lb,nmax,W,O,T);
  FREE(W);
}


void ccdl::AuxExpPot
( int const la, int const lb, 
  double const *__restrict__ crd,
  double const *__restrict__ O,
  double const * qa,
  double const * qb,
  double * pa,
  double * pb )
{
  int na = (la+1)*(la+1);
  int nb = (lb+1)*(lb+1);
  double *__restrict__ T = ALLOC(na*nb);
  ccdl::AuxExpInt(la,lb,crd,O,T);
  for ( int j=0; j<nb; ++j )
    for ( int i=0; i<na; ++i )
      {
	pa[i] += T[i+j*na]*qb[j];
	pb[j] += T[i+j*na]*qa[i];
      };
  FREE(T);
}



void ccdl::AuxExpPotGrd
( int const la, int const lb, 
  double const *__restrict__ crd,
  double const *__restrict__ O,
  double const * qa,
  double const * qb,
  double * pa,
  double * pb,
  double *__restrict__ g )
{
  int na = (la+1)*(la+1);
  int nb = (lb+1)*(lb+1);
  int lmax,lmin,nmax;
  if ( na > nb )
    {
      lmax = la;
      lmin = lb;
      nmax = na;
    }
  else
    {
      lmax = lb;
      lmin = la;
      nmax = nb;
    };
  double *__restrict__ T = ALLOC(2*na*nb);
  double *__restrict__ W = T + na*nb;
  ccdl::ClmTranslation( lmax, lmin, crd, W );
  ccdl::AuxExpInt_Helper(la,lb,nmax,W,O+1,T);
  ccdl::AuxExpPotGrd_Helper(la,lb,crd,T,qa,qb,g);
  ccdl::AuxExpInt_Helper(la,lb,nmax,W,O,T);
  for ( int j=0; j<nb; ++j )
    for ( int i=0; i<na; ++i )
      {
	pa[i] += T[i+j*na]*qb[j];
	pb[j] += T[i+j*na]*qa[i];
      };
  FREE(T);
}



/*
void ccdl::PrimGauExpGrd_Overlap
( int const la, int const lb, 
  double const zab,
  double const *__restrict__ crd,
  double const r2,
  double *__restrict__ T,
  double *__restrict__ dT )
{
  int const order = la+lb+1;
  int na = (la+1)*(la+1);
  int nb = (lb+1)*(lb+1);
  int lmax,lmin,nmax;
  if ( na > nb )
    {
      lmax = la;
      lmin = lb;
      nmax = na;
    }
  else
    {
      lmax = lb;
      lmin = la;
      nmax = nb;
    };

  double *__restrict__ O = ALLOC(order+1);
  double *__restrict__ W = ALLOC(na*nb);
  ccdl::ClmTranslation( lmax, lmin, crd, W );
  ccdl::PrimGauAuxVec_Overlap(order,zab,r2,O);
  ccdl::AuxExpInt_Helper(la,lb,nmax,W,O+1,T);
  ccdl::AuxExpGrd_Helper(la,lb,crd,T,dT);
  ccdl::AuxExpInt_Helper(la,lb,nmax,W,O,T);
  FREE(O);
  FREE(W);
}
*/

 /*
void ccdl::PrimGauExpGrd_Overlap
( int const la, int const lb, 
  double const zab,
  double const *__restrict__ crd,
  double const r2,
  double *__restrict__ T,
  double *__restrict__ dT )
{
  int const order = la+lb+1;
  double *__restrict__ O = ALLOC(order+1);
  ccdl::PrimGauAuxVec_Overlap(order,zab,r2,O);
  ccdl::AuxExpInt(la,lb,crd,O+1,T);
  ccdl::AuxExpGrd_Helper(la,lb,crd,T,dT);
  ccdl::AuxExpInt(la,lb,crd,O,T);
  FREE(O);
}
 */






void ccdl::PrimGauAuxVec_Coulomb
( int const order, 
  double const zab, 
  double const r2, 
  double *__restrict__ O )
{
  ccdl::GieseBoysFcn(order,zab*r2,O);
  double T = ccdl::TWO_OVER_SQRT_PI * std::sqrt(zab);
  O[0] *= T;
  double ms = T;
  for ( int k=1; k<order+1; ++k )
    {
      ms *= -zab;
      O[k] *= ms;
    };
}



void ccdl::PrimGauAuxVec_Overlap
( int const order, 
  double const zab, 
  double const r2, 
  double *__restrict__ O )
{
  O[0] = std::pow(zab/PI,1.5) * std::exp(-zab*r2);
  for ( int k=1; k < order+1; ++k )
    O[k] = -zab * O[k-1];
}



void ccdl::PrimGauAuxVec_Ewald
( int const order, 
  double const sqrt_za, 
  double const r2, 
  double *__restrict__ O )
{
  double oor = 0.;
  if ( r2 > 0. ) oor = 1. / std::sqrt(r2);
  double za = sqrt_za*sqrt_za;
  ccdl::GieseBoysFcn(order,za*r2,O);
  double t = ccdl::TWO_OVER_SQRT_PI * sqrt_za;
  double z  = 0.5 * oor * oor;
  *O = oor - *O*t;
  for ( int k=1; k<order+1; ++k )
    {
      oor *= -(2*k-1)*z;
      t *= -za;
      O[k] = oor-O[k]*t;
    };
}







void ccdl::PrimGauExpInt_Coulomb
( int const la, int const lb, 
  double const zab,
  double const * crd,
  double const r2,
  double * T )
{
  int const order = la+lb;
  double *__restrict__ O = ALLOC(order+1);
  ccdl::PrimGauAuxVec_Coulomb(order,zab,r2,O);
  ccdl::AuxExpInt(la,lb,crd,O,T);
  FREE(O);
}

void ccdl::PrimGauExpInt_Overlap
( int const la, int const lb, 
  double const zab,
  double const * crd,
  double const r2,
  double * T )
{
  int const order = la+lb;
  double *__restrict__ O = ALLOC(order+1);
  ccdl::PrimGauAuxVec_Overlap(order,zab,r2,O);
  ccdl::AuxExpInt(la,lb,crd,O,T);
  FREE(O);
}

void ccdl::PrimGauExpInt_Ewald
( int const la, int const lb, 
  double const sqrt_za,
  double const * crd,
  double const r2,
  double * T )
{
  int const order = la+lb;
  double *__restrict__ O = ALLOC(order+1);
  ccdl::PrimGauAuxVec_Ewald(order,sqrt_za,r2,O);
  ccdl::AuxExpInt(la,lb,crd,O,T);
  FREE(O);
}





void ccdl::PrimGauExpGrd_Coulomb
( int const la, int const lb, 
  double const zab,
  double const *__restrict__ crd,
  double const r2,
  double *__restrict__ T,
  double *__restrict__ dT )
{
  int const order = la+lb+1;
  double *__restrict__ O = ALLOC(order+1);
  ccdl::PrimGauAuxVec_Coulomb(order,zab,r2,O);
  ccdl::AuxExpGrd(la,lb,crd,O,T,dT);
  FREE(O);
}


void ccdl::PrimGauExpGrd_Overlap
( int const la, int const lb, 
  double const zab,
  double const *__restrict__ crd,
  double const r2,
  double *__restrict__ T,
  double *__restrict__ dT )
{
  int const order = la+lb+1;
  double *__restrict__ O = ALLOC(order+1);
  ccdl::PrimGauAuxVec_Overlap(order,zab,r2,O);
  ccdl::AuxExpGrd(la,lb,crd,O,T,dT);
  FREE(O);
}


void ccdl::PrimGauExpGrd_Ewald
( int const la, int const lb, 
  double const sqrt_za,
  double const *__restrict__ crd,
  double const r2,
  double *__restrict__ T,
  double *__restrict__ dT )
{
  int const order = la+lb+1;
  double *__restrict__ O = ALLOC(order+1);
  ccdl::PrimGauAuxVec_Ewald(order,sqrt_za,r2,O);
  ccdl::AuxExpGrd(la,lb,crd,O,T,dT);
  FREE(O);
}




void ccdl::PrimGauExpPot_Coulomb
( int const la, int const lb, 
  double const zab,
  double const *__restrict__ crd,
  double const r2,
  double const * qa,
  double const * qb,
  double * pa,
  double * pb )
{
  int const order = la+lb;
  double *__restrict__ O = ALLOC(order+1);
  ccdl::PrimGauAuxVec_Coulomb(order,zab,r2,O);
  ccdl::AuxExpPot(la,lb,crd,O,qa,qb,pa,pb);
  FREE(O);
}

void ccdl::PrimGauExpPot_Overlap
( int const la, int const lb, 
  double const zab,
  double const *__restrict__ crd,
  double const r2,
  double const * qa,
  double const * qb,
  double * pa,
  double * pb )
{
  int const order = la+lb;
  double *__restrict__ O = ALLOC(order+1);
  ccdl::PrimGauAuxVec_Overlap(order,zab,r2,O);
  ccdl::AuxExpPot(la,lb,crd,O,qa,qb,pa,pb);
  FREE(O);
}

void ccdl::PrimGauExpPot_Ewald
( int const la, int const lb, 
  double const sqrt_za,
  double const *__restrict__ crd,
  double const r2,
  double const * qa,
  double const * qb,
  double * pa,
  double * pb )
{
  int const order = la+lb;
  double *__restrict__ O = ALLOC(order+1);
  ccdl::PrimGauAuxVec_Ewald(order,sqrt_za,r2,O);
  ccdl::AuxExpPot(la,lb,crd,O,qa,qb,pa,pb);
  FREE(O);
}


void ccdl::PrimGauExpPotGrd_Coulomb
( int const la, int const lb, 
  double const zab,
  double const *__restrict__ crd,
  double const r2,
  double const * qa,
  double const * qb,
  double * pa,
  double * pb,
  double * g )
{
  int const order = la+lb+1;
  double *__restrict__ O = ALLOC(order+1);
  ccdl::PrimGauAuxVec_Coulomb(order,zab,r2,O);
  ccdl::AuxExpPotGrd(la,lb,crd,O,qa,qb,pa,pb,g);
  FREE(O);
}

void ccdl::PrimGauExpPotGrd_Overlap
( int const la, int const lb, 
  double const zab,
  double const *__restrict__ crd,
  double const r2,
  double const * qa,
  double const * qb,
  double * pa,
  double * pb,
  double * g )
{
  int const order = la+lb+1;
  double *__restrict__ O = ALLOC(order+1);
  ccdl::PrimGauAuxVec_Overlap(order,zab,r2,O);
  ccdl::AuxExpPotGrd(la,lb,crd,O,qa,qb,pa,pb,g);
  FREE(O);
}

void ccdl::PrimGauExpPotGrd_Ewald
( int const la, int const lb, 
  double const sqrt_za,
  double const *__restrict__ crd,
  double const r2,
  double const * qa,
  double const * qb,
  double * pa,
  double * pb,
  double * g )
{
  int const order = la+lb+1;
  double *__restrict__ O = ALLOC(order+1);
  ccdl::PrimGauAuxVec_Ewald(order,sqrt_za,r2,O);
  ccdl::AuxExpPotGrd(la,lb,crd,O,qa,qb,pa,pb,g);
  FREE(O);
}
