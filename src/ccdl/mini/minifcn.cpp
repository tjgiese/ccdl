#include "minifcn.hpp"
#include <cstdio>

ccdl::minifcn::minifcn( int nvar )
  : nvar(nvar), 
    is_new_reference(true), 
    refx(nvar,0.)
{
}


ccdl::minifcn::minifcn( int nvar, double const * RefX )
  : nvar(nvar), 
    is_new_reference(true), 
    refx(RefX,RefX+nvar)
{
}

ccdl::minifcn::~minifcn() {}


void ccdl::minifcn::SetNewReferencePt( double const alp, double const * s )
{
  for ( int i=0; i<nvar; ++i )
    refx[i] += alp*s[i];
  is_new_reference = true;
}


bool ccdl::minifcn::IsNewReferencePt() const
{
  return is_new_reference;
}



double ccdl::minifcn::numgrd( double const alp, double const * sin, double * dydx )
{
  double const DEL = 1.e-4;
  std::vector<double> s( sin, sin + nvar );
  for ( int i=0; i<nvar; ++i )
    {
      s[i] += DEL;
      double fhi = (*this)( alp, s.data() );
      s[i] -= 2*DEL;
      double flo = (*this)( alp, s.data() );
      s[i] += DEL;
      dydx[i] = (fhi-flo)/(2.*DEL);
    }
  return (*this)( alp, sin );
}

void ccdl::minifcn::testgrd()
{
  std::vector<double> s( nvar, 0. );
  double alp = 1.;

  std::vector<double> ana( nvar, 0. ), num( nvar, 0. );
  double fana = (*this)( alp, s.data(), ana.data() );
  double fnum = numgrd( alp, s.data(), num.data() );

  std::printf("fnum,fana,diff      %20.10e %20.10e %13.4e\n",
	      fnum,fana,fnum-fana);

  for ( int i=0; i<nvar; ++i )
  std::printf("gnum,gana,diff %4i %20.10e %20.10e %13.4e\n",
	      i,num[i],ana[i],num[i]-ana[i]);

}

