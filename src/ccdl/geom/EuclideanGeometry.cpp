#include "EuclideanGeometry.hpp"
#include "../constants.hpp"
#include <algorithm>

double ccdl::Bond
( double const * a, 
  double const * b )
{
  double const rab[] = { a[0]-b[0],a[1]-b[1],a[2]-b[2] };
  return std::sqrt( rab[0]*rab[0] + rab[1]*rab[1] + rab[2]*rab[2] );
}

double ccdl::Bond
( double const * a, 
  double const * b, 
  double * G )
{
  double const rab[] = { a[0]-b[0],a[1]-b[1],a[2]-b[2] };
  double const r = std::sqrt( rab[0]*rab[0] + rab[1]*rab[1] + rab[2]*rab[2] );
  for ( std::size_t k=0; k<3; ++k )
    {
      G[k+0] =  rab[k]/r;
      G[k+3] = -rab[k]/r;
    };
  return r;
}

double ccdl::Bond
( double const * a, 
  double const * b, 
  double * G,
  double * H )
{
  double const rab[] = { a[0]-b[0],a[1]-b[1],a[2]-b[2] };
  double const r2 =  rab[0]*rab[0] + rab[1]*rab[1] + rab[2]*rab[2];
  double const r = std::sqrt( r2 );
  for ( std::size_t k=0; k<3; ++k )
    {
      G[k+0] =  rab[k]/r;
      G[k+3] = -rab[k]/r;
    };
  for ( std::size_t B=0, ib=0; B<2; ++B )
    for ( std::size_t Bk=0; Bk < 3; ++Bk, ++ib )
      for ( std::size_t A=0, ia=0; A<2; ++A )
  	for ( std::size_t Ak=0; Ak < 3; ++Ak, ++ia )
	  if ( ia == ib )
	    H[ia+ib*6] = ( - G[ia]*G[ib] + 1. )/r;
	  else
	    H[ia+ib*6] = (Ak==Bk ? 2. : -1.) * G[ia]*G[ib] / r;
  return r;
}



double ccdl::CosineAngle
( double const * I, 
  double const * J, 
  double const * K )
{
  double const JI[] = { J[0]-I[0],J[1]-I[1],J[2]-I[2] };
  double const JK[] = { J[0]-K[0],J[1]-K[1],J[2]-K[2] };
  double const Lji = std::sqrt( JI[0]*JI[0] + JI[1]*JI[1] +  JI[2]*JI[2] );
  double const Ljk = std::sqrt( JK[0]*JK[0] + JK[1]*JK[1] +  JK[2]*JK[2] );
  double const ca = ( JI[0]*JK[0] + JI[1]*JK[1] + JI[2]*JK[2] )/(Lji*Ljk);
  return ca;
}



double ccdl::Angle
( double const * I, 
  double const * J, 
  double const * K )
{
  double ca = ccdl::CosineAngle(I,J,K);
  double a = 0.;
  if ( ca > 1.-1.e-9 )
    { a = 0.; }
  else if ( ca < -1.+1.e-9 )
    { a = ccdl::PI; }
  else
    { a = std::acos(ca); }
  return a;
}
#include <cstdio>

double ccdl::CosineAngle
( double const * I, 
  double const * J, 
  double const * K,
  double * G )
{
  double const JI[] = { J[0]-I[0],J[1]-I[1],J[2]-I[2] };
  double const JK[] = { J[0]-K[0],J[1]-K[1],J[2]-K[2] };
  double const Lji = std::sqrt( JI[0]*JI[0] + JI[1]*JI[1] +  JI[2]*JI[2] );
  double const Ljk = std::sqrt( JK[0]*JK[0] + JK[1]*JK[1] +  JK[2]*JK[2] );
  double const num =  JI[0]*JK[0] + JI[1]*JK[1] + JI[2]*JK[2];
  double const den = Lji*Ljk;
  double const ca = num/den;

  double const dcadnum =  1./den;
  double const dcadden = -ca/den;
  double const ddendLji = Ljk;
  double const ddendLjk = Lji;
  double const dLjidJ[] = { JI[0]/Lji, JI[1]/Lji, JI[2]/Lji};
  double const dLjidI[] = {-dLjidJ[0],-dLjidJ[1],-dLjidJ[2]};
  double const dLjkdJ[] = { JK[0]/Ljk, JK[1]/Ljk, JK[2]/Ljk};
  double const dLjkdK[] = {-dLjkdJ[0],-dLjkdJ[1],-dLjkdJ[2]};

  for ( std::size_t q=0; q<3; ++q )
    G[q+0] = dcadnum*(-JK[q]) + dcadden*(ddendLji*dLjidI[q]);
  for ( std::size_t q=0; q<3; ++q )
    G[q+3] = dcadnum*(JK[q]+JI[q]) + dcadden*(ddendLji*dLjidJ[q]+ddendLjk*dLjkdJ[q]);
  for ( std::size_t q=0; q<3; ++q )
    G[q+6] = dcadnum*(-JI[q]) + dcadden*(ddendLjk*dLjkdK[q]);

  return ca;
}
double ccdl::Angle
( double const * I, 
  double const * J, 
  double const * K,
  double * G )
{
  double dca[] = {0.,0.,0., 0.,0.,0., 0.,0.,0.};
  double const ca = ccdl::CosineAngle(I,J,K,dca);
  double a = 0.;
  double dadca = 0.;
  if ( ca > 1.-1.e-9 )
    { a = 0.; dadca = 0.; }
  else if ( ca < -1.+1.e-9 )
    { a = ccdl::PI; dadca = 0.; }
  else
    { a = std::acos(ca); dadca = -1./std::sin(a); }
  for ( std::size_t q=0; q<3; ++q )
    {
      G[q+0] = dadca * dca[q+0];
      G[q+3] = dadca * dca[q+3];
      G[q+6] = dadca * dca[q+6];
    };
  return a;
}

double ccdl::Angle
( double const * I, 
  double const * J, 
  double const * K,
  double * G,
  double * H )
{
  double u[3] = { I[0]-J[0],I[1]-J[1],I[2]-J[2] };
  double Lu = std::sqrt( u[0]*u[0] + u[1]*u[1] + u[2]*u[2] );
  u[0] /= Lu; u[1] /= Lu; u[2] /= Lu; 
  double v[3] = { K[0]-J[0],K[1]-J[1],K[2]-J[2] };
  double Lv = std::sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
  v[0] /= Lv; v[1] /= Lv; v[2] /= Lv; 
  double w[3] = {0.,0.,0.};

  double ca =  u[0]*v[0] + u[1]*v[1] +u[2]*v[2];
  double a = 0.;
  if ( ca > 1.-1.e-9 )
    { ca = 1.; a = 0.; }
  else if ( ca < -1.+1.e-9 )
    { ca =-1.; a = ccdl::PI; }
  else
    { a = std::acos(ca); }
  double sa = std::sqrt( 1. - ca*ca );

  ccdl::CrossProduct(u,v,w);
  double Lw = std::sqrt( w[0]*w[0] + w[1]*w[1] + w[2]*w[2] );
  std::fill( H, H+81, 0. );
  if ( Lw > 1.e-20 )
    { w[0] /= Lw; w[1] /= Lw; w[2] /= Lw; }
  else
    { w[0]  = 0.; w[1]  = 0.; w[2]  = 0.; }

  double uw[3] = {0.,0.,0.};
  double wv[3] = {0.,0.,0.};
  ccdl::CrossProduct(u,w,uw);
  ccdl::CrossProduct(w,v,wv);

#define ZETA(a,m,n) ((a)==(m)?1.:((a)==(n)?-1.:0.))
#define DELTA(i,j) ((i)==(j)?1.:0.)
  for ( int a=0; a<3; ++a )
    for ( int ak=0; ak<3; ++ak )
      G[ak+a*3] = ZETA(a,0,1)*uw[ak]/Lu + ZETA(a,2,1)*wv[ak]/Lv;

  if ( sa < 1.e-20 ) return a;

  for (int a=0; a<3; ++a)
    for (int ak=0; ak<3; ++ak)
      for (int b=0; b<3; ++b)
        for (int bk=0; bk<3; ++bk) 
	  {
	    double t = 0.;
	    t = ZETA(a,0,1)*ZETA(b,0,1)*
	      (u[ak]*v[bk]+u[bk]*v[ak]-3*u[ak]*u[bk]*ca+DELTA(ak,bk)*ca) / (Lu*Lu*sa);
	    t +=ZETA(a,2,1)*ZETA(b,2,1)*
	      (v[ak]*u[bk]+v[bk]*u[ak]-3*v[ak]*v[bk]*ca+DELTA(ak,bk)*ca) / (Lv*Lv*sa);
	    t +=ZETA(a,0,1)*ZETA(b,2,1)*
	      (u[ak]*u[bk]+v[bk]*v[ak]-u[ak]*v[bk]*ca-DELTA(ak,bk)) / (Lu*Lv*sa);
	    t +=ZETA(a,2,1)*ZETA(b,0,1)*
	      (v[ak]*v[bk]+u[bk]*u[ak]-v[ak]*u[bk]*ca-DELTA(ak,bk)) / (Lu*Lv*sa);
	    t -= ca / sa * G[ak+a*3] * G[bk+b*3];
	    H[(ak+a*3)+(bk+b*3)*9] = t;
	  }

#undef ZETA
#undef DELTA
  return a;
}

double ccdl::DihedralAngle( double const * Ra, double const * Rb, double const * Rc, double const * Rd, bool & linear )
{
  double d = 0.;

  long double const Rba[]={Rb[0]-Ra[0],Rb[1]-Ra[1],Rb[2]-Ra[2]};
  long double const Rcb[]={Rc[0]-Rb[0],Rc[1]-Rb[1],Rc[2]-Rb[2]};
  long double const Rdc[]={Rd[0]-Rc[0],Rd[1]-Rc[1],Rd[2]-Rc[2]};

  long double A[3],B[3];
  ccdl::CrossProduct(Rba,Rcb,A);
  ccdl::CrossProduct(Rcb,Rdc,B);
  long double const n2A = A[0]*A[0]+A[1]*A[1]+A[2]*A[2];
  long double const n2B = B[0]*B[0]+B[1]*B[1]+B[2]*B[2];
  long double const den = std::sqrt(n2A*n2B);

  linear = true;
  if ( den > 1.e-30 ) 
    {
      linear = false;
      long double x = -(A[0]*B[0] + A[1]*B[1] + A[2]*B[2])/den;
      // if ( x < -1. ) 
      // 	d = ccdl::PI;
      // else if ( x > 1. )
      // 	d = 0.;
      // else
      // 	d = std::acos( x );
      if ( x <= -1. ) 
	d = 0.;
      else if ( x >= 1. )
	d = -ccdl::PI;
      else
	d = std::acos( x )-ccdl::PI;

      long double BxA[3];
      ccdl::CrossProduct(B,A,BxA);

      //long double BondProj = Rcb[0]*BxA[0] + Rcb[1]*BxA[1] + Rcb[2]*BxA[2];
      //if ( BondProj > 0.0 ) d = -d;
      //d = ccdl::PI - d;
    }

  return d;
}


double ccdl::DihedralAngle( double const * Ra, double const * Rb, double const * Rc, double const * Rd, bool & linear, double * G )
{
  double * ddda = G+0;
  double * dddb = G+3;
  double * dddc = G+6;
  double * dddd = G+9;

  double d = 0.;

  long double const Rba[]={Rb[0]-Ra[0],Rb[1]-Ra[1],Rb[2]-Ra[2]};
  long double const Rcb[]={Rc[0]-Rb[0],Rc[1]-Rb[1],Rc[2]-Rb[2]};
  long double const Rdc[]={Rd[0]-Rc[0],Rd[1]-Rc[1],Rd[2]-Rc[2]};
  long double const Rca[]={Rc[0]-Ra[0],Rc[1]-Ra[1],Rc[2]-Ra[2]};
  long double const Rdb[]={Rd[0]-Rb[0],Rd[1]-Rb[1],Rd[2]-Rb[2]};

  long double A[3],B[3];
  ccdl::CrossProduct(Rba,Rcb,A);
  ccdl::CrossProduct(Rcb,Rdc,B);
  long double const n2A = A[0]*A[0]+A[1]*A[1]+A[2]*A[2];
  long double const n2B = B[0]*B[0]+B[1]*B[1]+B[2]*B[2];
  long double const den = std::sqrt(n2A*n2B);
  linear = true;
  if ( den > 1.e-30 ) 
    {
      linear = false;
      long double x = -(A[0]*B[0] + A[1]*B[1] + A[2]*B[2])/den;
      // if ( x <= -1. ) 
      // 	d = ccdl::PI;
      // else if ( x >= 1. )
      // 	d = 0.;
      // else
      // 	d = std::acos( x );
      if ( x <= -1. ) 
	d = 0.;
      else if ( x >= 1. )
	d = -ccdl::PI;
      else
	d = std::acos( x )-ccdl::PI;

      long double BxA[3];
      ccdl::CrossProduct(B,A,BxA);

      // long double BondProj = Rcb[0]*BxA[0] + Rcb[1]*BxA[1] + Rcb[2]*BxA[2];
      // if ( BondProj > 0.0 ) d = -d;
      // d = ccdl::PI - d;

      long double const ncb = std::sqrt(Rcb[0]*Rcb[0]+Rcb[1]*Rcb[1]+Rcb[2]*Rcb[2]);

      long double const sclA = 1. / ( n2A*ncb );
      A[0] *= sclA;
      A[1] *= sclA;
      A[2] *= sclA;
      long double AxRcb[3];
      ccdl::CrossProduct(A,Rcb,AxRcb);

      long double const sclB = 1. / ( n2B*ncb );
      B[0] *= sclB;
      B[1] *= sclB;
      B[2] *= sclB;

      long double RcbxB[3];
      ccdl::CrossProduct(Rcb,B,RcbxB);

      ccdl::CrossProduct(AxRcb,Rcb,ddda);
      ccdl::CrossProduct(Rca,AxRcb,dddb);
      ccdl::AccumulateCrossProduct(RcbxB,Rdc,dddb);
      ccdl::CrossProduct(AxRcb,Rba,dddc);
      ccdl::AccumulateCrossProduct(Rdb,RcbxB,dddc);
      ccdl::CrossProduct(RcbxB,Rcb,dddd);
    }

  return d;
}



#define ZETA(a,m,n) ((a)==(m)?1.l:((a)==(n)?-1.l:0.l))
#define DELTA(i,j) ((i)==(j)?1.l:0.l)

double ccdl::DihedralAngle( double const * Ra, double const * Rb, double const * Rc, double const * Rd, bool & linear, double * G, double * H )
{
  double d = 0.;
  linear = true;
  std::fill( G, G+4*3, 0. );
  std::fill( H, H+4*3*4*3, 0. );
  long double u[3]={Ra[0]-Rb[0],Ra[1]-Rb[1],Ra[2]-Rb[2]};
  long double v[3]={Rd[0]-Rc[0],Rd[1]-Rc[1],Rd[2]-Rc[2]};
  long double w[3]={Rc[0]-Rb[0],Rc[1]-Rb[1],Rc[2]-Rb[2]};
  long double Lu = std::sqrt( u[0]*u[0]+u[1]*u[1]+u[2]*u[2] );
  long double Lv = std::sqrt( v[0]*v[0]+v[1]*v[1]+v[2]*v[2] );
  long double Lw = std::sqrt( w[0]*w[0]+w[1]*w[1]+w[2]*w[2] );
  u[0] /= Lu; u[1] /= Lu; u[2] /= Lu;
  v[0] /= Lv; v[1] /= Lv; v[2] /= Lv;
  w[0] /= Lw; w[1] /= Lw; w[2] /= Lw;
  long double cos_u =  u[0]*w[0]+u[1]*w[1]+u[2]*w[2];
  long double cos_v = -v[0]*w[0]-v[1]*w[1]-v[2]*w[2];
  long double sin_u = std::sqrt(1.-cos_u*cos_u);
  long double sin_v = std::sqrt(1.-cos_v*cos_v);
  long double uXw[3]={0.l,0.l,0.l}, vXw[3]={0.l,0.l,0.l};
  ccdl::CrossProduct(u,w,uXw);
  ccdl::CrossProduct(v,w,vXw);
  long double const n1 = uXw[0]*uXw[0]+uXw[1]*uXw[1]+uXw[2]*uXw[2];
  long double const n2 = vXw[0]*vXw[0]+vXw[1]*vXw[1]+vXw[2]*vXw[2];
  long double const den = std::sqrt(n1*n2);
  if ( den > 1.e-30l )
    {
      linear = false;

      long double x = -(uXw[0]*vXw[0]+uXw[1]*vXw[1]+uXw[2]*vXw[2])/den;
      if ( x <= -1. ) 
	d = 0.;
      else if ( x >= 1. )
	d = -ccdl::PI;
      else
	d = std::acos( x )-ccdl::PI;

      long double tval, tval1, tval2, tval3, tval4;
      for ( int a=0; a<4; ++a )
	for ( int i=0; i<3; ++i ) 
	  { 
	    tval1 = tval2 = tval3 = tval4 = 0.0l;
	    
	    if ((a == 0) || (a == 1))
	      tval1 = ZETA(a,0,1) * uXw[i] / (Lu*sin_u*sin_u);
	    
	    if ((a == 2) || (a == 3))
	      tval2 = ZETA(a,2,3) * vXw[i] / (Lv*sin_v*sin_v);
	    
	    if ((a == 1) || (a == 2))
	      tval3 =   ZETA(a,1,2) * uXw[i]*cos_u/(Lw*sin_u*sin_u);
	    
	    if ((a == 1) || (a == 2)) 
	      tval4 = - ZETA(a,2,1) * vXw[i]*cos_v/(Lw*sin_v*sin_v);
	    
	    G[i+a*3] = tval1 + tval2 + tval3 + tval4;
	  }


      long double sinu4 = sin_u*sin_u*sin_u*sin_u;
      long double sinv4 = sin_v*sin_v*sin_v*sin_v;
      long double cosu3 = cos_u*cos_u*cos_u;
      long double cosv3 = cos_v*cos_v*cos_v;
      

      int k; // cartesian ; not i or j
      for (int a=0; a<4; ++a) 
	for (int b=0; b<=a; ++b)
	  for (int i=0; i<3; ++i)
	    for (int j=0; j<3; ++j) 
	      {
		tval = 0.l;

		if ((a==0 && b==0) || (a==1 && b==0) || (a==1 && b ==1))
		  tval +=  ZETA(a,0,1)*ZETA(b,0,1)*(uXw[i]*(w[j]*cos_u-u[j]) + uXw[j]*(w[i]*cos_u-u[i])) / (Lu*Lu*sinu4);
		
		if ((a==3 && b==3) || (a==3 && b==2) || (a==2 && b==2))
		  tval += ZETA(a,3,2)*ZETA(b,3,2)*(vXw[i]*(w[j]*cos_v+v[j]) + vXw[j]*(w[i]*cos_v+v[i])) / (Lv*Lv*sinv4);

		if ((a==1 && b==1) || (a==2 && b==1) || (a==2 && b==0) || (a==1 && b==0))
		  tval +=   (ZETA(a,0,1)*ZETA(b,1,2)+ZETA(a,2,1)*ZETA(b,1,0)) *
		    (uXw[i] * (w[j] - 2*u[j]*cos_u + w[j]*cos_u*cos_u) +
		     uXw[j] * (w[i] - 2*u[i]*cos_u + w[i]*cos_u*cos_u)) / (2*Lu*Lw*sinu4);
		
		if ((a==3 && b==2) || (a==3 && b==1) || (a==2 && b==2) || (a==2 && b==1))
		  tval +=  (ZETA(a,3,2)*ZETA(b,2,1)+ZETA(a,1,2)*ZETA(b,2,3)) *
		    (vXw[i] * (w[j] + 2*v[j]*cos_v + w[j]*cos_v*cos_v) +
		     vXw[j] * (w[i] + 2*v[i]*cos_v + w[i]*cos_v*cos_v)) / (2*Lv*Lw*sinv4);
		
		if ((a==1 && b==1) || (a==2 && b==2) || (a==2 && b==1))
		  tval +=   ZETA(a,1,2)*ZETA(b,2,1)*
		    (uXw[i]*(u[j] + u[j]*cos_u*cos_u - 3*w[j]*cos_u + w[j]*cosu3) +
		     uXw[j]*(u[i] + u[i]*cos_u*cos_u - 3*w[i]*cos_u + w[i]*cosu3)) / (2*Lw*Lw*sinu4);
		
		if ((a==2 && b==1) || (a==2 && b==2) || (a==1 && b==1))
		  tval +=  ZETA(a,2,1)*ZETA(b,1,2)*
		    (vXw[i]*(-v[j] - v[j]*cos_v*cos_v - 3*w[j]*cos_v + w[j]*cosv3) +
		     vXw[j]*(-v[i] - v[i]*cos_v*cos_v - 3*w[i]*cos_v + w[i]*cosv3)) / (2*Lw*Lw*sinv4);

		if ((a != b) && (i != j)) 
		  {
		    if ( i!=0 && j!=0 ) k = 0;
		    else if ( i!=1 && j!=1 ) k = 1;
		    else k = 2;
		    
		    if (a==1 && b==1)
		      tval +=  ZETA(a,0,1)*ZETA(b,1,2) * (j-i) *
			pow(-0.5, std::abs(j-i)) * (+w[k]*cos_u - u[k]) / (Lu*Lw*sin_u*sin_u);
		    
		    if ((a==3 && b==2) || (a==3 && b==1) || (a==2 && b==2) || (a==2 && b==1))
		      tval +=  ZETA(a,3,2)*ZETA(b,2,1) * (j-i) *
			pow(-0.5, std::abs(j-i)) * (-w[k]*cos_v - v[k]) / (Lv*Lw*sin_v*sin_v);
		    
		    if ((a==2 && b==1) || (a==2 && b==0) || (a==1 && b==1) || (a==1 && b==0))
		      tval +=  ZETA(a,2,1)*ZETA(b,1,0) * (j-i) *
			pow(-0.5, std::abs(j-i)) * (-w[k]*cos_u + u[k]) / (Lu*Lw*sin_u*sin_u);
		    
		    if (a==2 && b==2)
		      tval +=  ZETA(a,1,2)*ZETA(b,2,3) * (j-i) *
			pow(-0.5, std::abs(j-i)) * (+w[k]*cos_v + v[k]) / (Lv*Lw*sin_v*sin_v);
		  }
		H[ (i+a*3) + (j+b*3)*4*3 ] = tval;
		H[ (j+b*3) + (i+a*3)*4*3 ] = tval;
	      }       
    };

  return d;
}
