#include <algorithm>
#include <iostream>
#include "GaussianBasis.hpp"
#include "GaussianBasis/GaussianBasisDatabase.hpp"


GaussianBasis::Basis::Basis
( char const * name, int const atNum )
  : mType(name), 
    mAtomId(atNum),
    mFockSize(0),
    mLmax(0)
{
  std::transform( mType.begin(), mType.end(), mType.begin(), ::toupper );
  int Nseg = 0;
  int const * Lvec = NULL;
  int const * Nprim = NULL;
  double const * Cmat;
  double const * Zmat;
  
  if ( mType.compare( "STO-3G" ) == 0 )
    {
      GaussianBasis::STO3G( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "STO-6G" ) == 0 )
    {
      GaussianBasis::STO6G( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "MG3S" ) == 0 )
    {
      GaussianBasis::MG3S( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "3-21G" ) == 0 )
    {
      GaussianBasis::Pople321G( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "3-21GSP" ) == 0 )
    {
      GaussianBasis::Pople321GSP( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "6-31G" ) == 0 )
    {
      GaussianBasis::Pople631G( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "6-31G*" ) == 0 )
    {
      GaussianBasis::Pople631Gs( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "6-31G**" ) == 0 )
    {
      GaussianBasis::Pople631Gss( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "6-31++G**" ) == 0 )
    {
      GaussianBasis::Pople631ppGss( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "6-311G" ) == 0 )
    {
      GaussianBasis::Pople6311G( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "6-311G**" ) == 0 )
    {
      GaussianBasis::Pople6311Gss( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "6-311++G**" ) == 0 )
    {
      GaussianBasis::Pople6311ppGss( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "CC-PVDZ" ) == 0 )
    {
      GaussianBasis::ccpVDZ( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "CC-PVTZ" ) == 0 )
    {
      GaussianBasis::ccpVTZ( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "CC-PVQZ" ) == 0 )
    {
      GaussianBasis::ccpVQZ( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "CC-PV5Z" ) == 0 )
    {
      GaussianBasis::ccpV5Z( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "AUG-CC-PVDZ" ) == 0 )
    {
      GaussianBasis::augccpVDZ( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "AUG-CC-PVTZ" ) == 0 )
    {
      GaussianBasis::augccpVTZ( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else if ( mType.compare( "AUG-CC-PVQZ" ) == 0 )
    {
      GaussianBasis::augccpVQZ( mAtomId, Nseg, Lvec, Nprim, Cmat, Zmat );
    }
  else
    {
       std::cerr << "Invalid basis set " << mType << std::endl;
    };

  if ( Lvec == NULL or Nseg == 0 )
    {
      std::cerr <<"Invalid element (" << atNum << ") for basis set '" << mType << "'" << std::endl;
      abort();
    }
  else
    {
      double const * C = Cmat;
      double const * Z = Zmat;
      for ( int iseg=0; iseg<Nseg; ++iseg )
	{
	  PushSegment( GaussianBasis::Segment( Lvec[iseg], Nprim[iseg], C, Z ) );
	  C += Nprim[iseg];
	  Z += Nprim[iseg];
	};
    };
}

