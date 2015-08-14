#include "EuclideanGeometry.hpp"
#include "../constants.hpp"


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
  double * dbdI, 
  double * dbdJ )
{
  double const rab[] = { a[0]-b[0],a[1]-b[1],a[2]-b[2] };
  double const r = std::sqrt( rab[0]*rab[0] + rab[1]*rab[1] + rab[2]*rab[2] );
  for ( std::size_t k=0; k<3; ++k )
    {
      dbdI[k] =  rab[k]/r;
      dbdJ[k] = -rab[k]/r;
    };
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
  return std::acos( ccdl::CosineAngle(I,J,K) );
}

double ccdl::CosineAngle
( double const * I, 
  double const * J, 
  double const * K,
  double * dcadI,
  double * dcadJ,
  double * dcadK )
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
    dcadI[q] = dcadnum*(-JK[q]) + dcadden*(ddendLji*dLjidI[q]);
  for ( std::size_t q=0; q<3; ++q )
    dcadJ[q] = dcadnum*(JK[q]+JI[q]) + dcadden*(ddendLji*dLjidJ[q]+ddendLjk*dLjkdJ[q]);
  for ( std::size_t q=0; q<3; ++q )
    dcadK[q] = dcadnum*(-JI[q]) + dcadden*(ddendLjk*dLjkdK[q]);

  return ca;
}
double ccdl::Angle
( double const * I, 
  double const * J, 
  double const * K,
  double * dadI,
  double * dadJ,
  double * dadK )
{
  double dcadI[] = {0.,0.,0.};
  double dcadJ[] = {0.,0.,0.};
  double dcadK[] = {0.,0.,0.};
  double const ca = ccdl::CosineAngle(I,J,K,dcadI,dcadJ,dcadK);
  double const a  = std::acos(ca);
  //double const dadca = -1./std::sqrt(1.-ca*ca);
  double const dadca = -1./std::sin(a);
  for ( std::size_t q=0; q<3; ++q )
    {
      dadI[q] = dadca * dcadI[q];
      dadJ[q] = dadca * dcadJ[q];
      dadK[q] = dadca * dcadK[q];
    };
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

  // std::cout << "den = " 
  //        << std::scientific << std::setw(12) << std::setprecision(3) 
  //        << den << "\n";
  linear = true;
  if ( den > 1.e-30 ) 
    {
      linear = false;
      // std::cout << "acos( " 
      //                << std::scientific << std::setw(12) << std::setprecision(3) 
      //                << -(A[0]*B[0] + A[1]*B[1] + A[2]*B[2])/den
      //                << "\n";
      long double x = -(A[0]*B[0] + A[1]*B[1] + A[2]*B[2])/den;
      if ( x < -1. ) 
        {
          d = ccdl::PI;
        }
      else if ( x > 1. )
        {
          d = 0.;
        }
      else
        {
          d = std::acos( x );
        };

      long double BxA[3];
      ccdl::CrossProduct(B,A,BxA);

      long double BondProj = Rcb[0]*BxA[0] + Rcb[1]*BxA[1] + Rcb[2]*BxA[2];
      if ( BondProj > 0.0 ) d = -d;
      d = ccdl::PI - d;
    }

  return d;
}


double ccdl::DihedralAngle( double const * Ra, double const * Rb, double const * Rc, double const * Rd, bool & linear, double * ddda, double * dddb, double * dddc, double * dddd )
{
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
      if ( x <= -1. ) 
        {
          d = ccdl::PI;
        }
      else if ( x >= 1. )
        {
          d = 0.;
        }
      else
        {
          d = std::acos( x );
        };


      //d = std::acos( -(A[0]*B[0] + A[1]*B[1] + A[2]*B[2])/den );

      long double BxA[3];
      ccdl::CrossProduct(B,A,BxA);

      long double BondProj = Rcb[0]*BxA[0] + Rcb[1]*BxA[1] + Rcb[2]*BxA[2];
      if ( BondProj > 0.0 ) d = -d;
      d = ccdl::PI - d;

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

      // std::cout << std::scientific << std::setw(12) << std::setprecision(3) 
      //                << sclA << " " 
      //                << std::scientific << std::setw(12) << std::setprecision(3)
      //                << sclB << "\n";

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
