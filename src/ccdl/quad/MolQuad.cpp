#include "MolQuad.hpp"
#include "Quadrature.hpp"
#include <cassert>
#include <cstddef>

ccdl::AtomQuadParam::AtomQuadParam
( int const numShells, 
  int const numAngPts, 
  double const outerRadius, 
  double const atomRadius )
  : mNumAngPts(numShells,numAngPts),
    mAtomRadius(atomRadius),
    mRadialPts(numShells,0.),
    mRadialWts(numShells,0.),
    mNumPts(numShells*numAngPts),
    mShellOffsets(numShells+1)
{
  ccdl::GaussLaguerreRule( 0., numShells, mRadialPts.data(), mRadialWts.data() );
  double const lambda = mRadialPts.back() / outerRadius;
  for ( int i=0; i<numShells; ++i )
    {
      mRadialWts[i] /= ( lambda * std::exp( -mRadialPts[i] ) );
      mRadialPts[i] /= lambda;
      mRadialWts[i] *= mRadialPts[i] * mRadialPts[i];
    };
  BuildShellOffsets();
}


void ccdl::AtomQuadParam::Prune( int const irad, int const numAngPts )
{
  int const nrad = GetNumRadialShells();
  for ( int i=irad; i < nrad; ++i )
    mNumAngPts[i] = numAngPts;
  mNumPts = 0;
  for ( int i=0; i < nrad; ++i )
    mNumPts += mNumAngPts[i];
  BuildShellOffsets();
}


void ccdl::AtomQuadParam::BuildShellOffsets()
{
  int const nrad = GetNumRadialShells();
  mShellOffsets[0] = 0;
  for ( int irad=0; irad < nrad; ++irad )
    mShellOffsets[irad+1] = mShellOffsets[irad] + GetNumAngPts(irad);
}






ccdl::MolQuad::MolQuad
( ccdl::AtomQuadParams atomParam )
  : mAtomParam( 1, atomParam[0] ),
    mAtomCrd( 3*1 ),
    mAtomOffsets( 1+1, 0 ),
    mNumAtoms( 1 )
{
  mAtomCrd[0] = 0.;
  mAtomCrd[1] = 0.;
  mAtomCrd[2] = 0.;
  BuildGrid();
}


// ccdl::MolQuad::MolQuad
// ( std::tr1::shared_ptr< ccdl::AtomQuadParam > atomParam )
//   : mAtomParam( 1, atomParam ),
//     mAtomCrd( 3*1 ),
//     mAtomOffsets( 1+1, 0 ),
//     mNumAtoms( 1 )
// {
//   mAtomCrd[0] = 0.;
//   mAtomCrd[1] = 0.;
//   mAtomCrd[2] = 0.;
//   BuildGrid();
// }


ccdl::MolQuad::MolQuad
()
  : mAtomParam( 1 ),
    mAtomCrd( 3*1 ),
    mAtomOffsets( 1+1, 0 ),
    mNumAtoms( 1 )
{
  ccdl::AtomQuadParams mqp;
  mAtomParam[0] = mqp[0];
  mAtomCrd[0] = 0.;
  mAtomCrd[1] = 0.;
  mAtomCrd[2] = 0.;
  BuildGrid();
}




void ccdl::MolQuad::BuildGrid()
{
  assert( mNumAtoms == (int)mAtomParam.size() );
  assert( mNumAtoms+1 == (int)mAtomOffsets.size() );

  int nat = GetNumAtoms();
  mAtomOffsets[0] = 0;
  for ( int iat=0; iat < nat; ++iat )
    mAtomOffsets[iat+1] = mAtomOffsets[iat] + GetNumAtomPts(iat);


  mNumRadialShells = 0;
  for ( int iat=0; iat < nat; ++iat )
    mNumRadialShells += GetNumRadialShells( iat );
  mShellOffsets.resize( mNumRadialShells + 1 );
  mShellOwner.resize( mNumRadialShells );
  mShellOffsets[0] = 0;
  int shlidx = 0;
  for ( int iat=0; iat < nat; ++iat )
    {
      int n = GetNumRadialShells( iat );
      for ( int ishl = 0; ishl < n; ++ishl, ++shlidx )
	{
	  mShellOffsets[shlidx+1] = mShellOffsets[shlidx] + GetNumAngPts( iat, ishl );
	  mShellOwner[shlidx] = iat;
	};
    };

  mNumPts = mAtomOffsets.back();
  mSingleCenterWt.resize( GetNumPts() );
  mPartitionWt.resize( GetNumPts() );
  mTotalWt.resize( GetNumPts() );
  mQuadCrd.resize( 3*GetNumPts() );

  for ( int i=0; i<GetNumPts(); ++i )
    mPartitionWt[i] = 1.;
  
  for ( int iat=0; iat < nat; ++iat )
    {
      int const nrad = GetNumRadialShells(iat);
      for ( int irad=0; irad<nrad; ++irad )
	{
	  int nang = GetNumAngPts(iat,irad);
	  int o = GetRadialShellBegin(iat,irad);
	  double const radius = GetAtomParam(iat).GetRadialPt(irad);
	  double const radWt  = GetAtomParam(iat).GetRadialWt(irad);
	  double * angWts = mAngQuadDB.GetWts(nang)->data();
	  std::tr1::array<double,3> * angPts = mAngQuadDB.GetPts(nang)->data();

	  for ( int iang=0; iang < nang; ++iang, ++o )
	    {
	      mSingleCenterWt[o] = radWt * angWts[iang];
	      mQuadCrd[0+o*3] = mAtomCrd[0+iat*3] + radius * angPts[iang][0];
	      mQuadCrd[1+o*3] = mAtomCrd[1+iat*3] + radius * angPts[iang][1];
	      mQuadCrd[2+o*3] = mAtomCrd[2+iat*3] + radius * angPts[iang][2];
	    };
	};
    };

  mTotalWt = mSingleCenterWt;
}



namespace ccdl
{
  namespace backend
  {
    class pair
    {
    public:
      pair() : atom(0), rij(0.) {};
      pair( int const atom, double const rij ) : atom(atom), rij(rij) {};
      bool operator< ( pair const & rhs ) const { return rij<rhs.rij; }
      bool operator> ( pair const & rhs ) const { return rij>rhs.rij; }
      bool operator== ( pair const & rhs ) const { return rij==rhs.rij; }

      static bool RijLessThan( pair const & lhs, pair const & rhs ) { return lhs.atom<rhs.atom; }
      static bool AtomLessThan( pair const & lhs, pair const & rhs ) { return lhs.atom<rhs.atom; }

      int atom;
      double rij;
    };

  }
}



void ccdl::MolQuad::SetAtomCrd
( int const iat, double const * c )
{
  double const d[3] = { c[0]-mAtomCrd[0+iat*3], 
			c[1]-mAtomCrd[1+iat*3], 
			c[2]-mAtomCrd[2+iat*3] };
  mAtomCrd[0+iat*3] = c[0];
  mAtomCrd[1+iat*3] = c[1];
  mAtomCrd[2+iat*3] = c[2];
  int const begin = GetAtomBegin(iat);
  int const end   = GetAtomEnd(iat);
  for ( int ipt=begin; ipt<end; ++ipt )
    {
      mQuadCrd[0+ipt*3] += d[0];
      mQuadCrd[1+ipt*3] += d[1];
      mQuadCrd[2+ipt*3] += d[2];
    };
}





void ccdl::MolQuad::EvaluateBeckeWeights_WithRadiiBruteForce()
{
  int const nat = GetNumAtoms();


  std::vector<double> RIJM(nat*(nat-1),0.);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for ( int jat=1; jat<nat; ++jat )
    {
      double const * rb = GetAtomCrdPtr(jat);
      for ( int iat=0; iat<jat; ++iat )
	{
	  double const * ra = GetAtomCrdPtr(iat);
	  double t[3] = { ra[0]-rb[0], ra[1]-rb[1], ra[2]-rb[2] };
	  double const rij = std::sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
	  double x = GetAtomRadius(iat)/GetAtomRadius(jat);
	  double u = (x-1.)/(x+1.);
	  double a = u/(u*u-1.);
	  if ( a > 0.5 )
	    {
	      a = 0.5;
	    }
	  else if ( a < -0.5 )
	    {
	      a = -0.5;
	    };
	  int const k = 2*iat+jat*(jat-1);
	  RIJM[0+k] = rij;
	  RIJM[1+k] = a;
	};
    };

#ifdef _OPENMP
#pragma omp parallel
  {
#endif
  std::vector<double> w( nat, 0. );
  std::vector<double> d( nat, 0. );

#ifdef _OPENMP

  int const npt = GetNumPts();
#pragma omp for
  for ( int ipt=0; ipt<npt; ++ipt )
    {
      int const gat = GetAtomIdx(ipt);

#else

  for ( int gat=0; gat<nat; ++gat )
    {
      int const begin = GetAtomBegin(gat);
      int const end   = GetAtomEnd(gat);
      for ( int ipt=begin; ipt<end; ++ipt )
  	{

#endif
	  w.assign(nat,1.);
      
	  double const * rq = GetQuadCrdPtr(ipt);
	  
	  for ( int iat=0; iat<nat; ++iat )
	    {
	      double const * ra = GetAtomCrdPtr(iat);
	      double const t[3] = { rq[0]-ra[0],
				    rq[1]-ra[1],
				    rq[2]-ra[2] };
	      d[iat] = std::sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
	    };
	  
	  for ( int jat=1; jat<nat; ++jat )
	    for ( int iat=0; iat<jat; ++iat )
	      {
		int const k = 2*iat+jat*(jat-1);
		double const rij = RIJM[0+k];
		double const a   = RIJM[1+k];
		double mu = rij > 1.e-8 ? ( d[iat] - d[jat] ) / rij : 1.;
		mu += a - a * mu * mu;
		double p = ( 1.5 - 0.5 * mu*mu ) * mu;
		p *= ( 1.5 - 0.5 * p*p );
		p *= ( 1.5 - 0.5 * p*p );
		w[iat] *= ( 0.5 - 0.5*p );
		w[jat] *= ( 0.5 + 0.5*p );
	      };
	  
	  double wt = 0.;
	  for ( int iat=0; iat<nat; ++iat )
	    wt += w[iat];
	  wt = w[gat]/wt;
	  
	  SetPartitionWt(ipt,wt);
	  
	};
      
    };
}






void ccdl::MolQuad::EvaluateBeckeWeights_WithoutRadiiBruteForce()
{
  int const nat = GetNumAtoms();

  if ( nat < 1 ) { return; }

  std::vector<double> RIJM(nat*(nat-1)/2,0.);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for ( int jat=1; jat<nat; ++jat )
    {
      double const * rb = GetAtomCrdPtr(jat);
      for ( int iat=0; iat<jat; ++iat )
	{
	  double const * ra = GetAtomCrdPtr(iat);
	  double t[3] = { ra[0]-rb[0], ra[1]-rb[1], ra[2]-rb[2] };
	  RIJM[ iat+jat*(jat-1)/2 ] = std::sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
	};
    };

#ifdef _OPENMP
#pragma omp parallel
  {
#endif
  std::vector<double> w( nat, 0. );
  std::vector<double> d( nat, 0. );

#ifdef _OPENMP

  int const npt = GetNumPts();
#pragma omp for schedule(dynamic)
  for ( int ipt=0; ipt<npt; ++ipt )
    {
      int const gat = GetAtomIdx(ipt);

#else

  for ( int gat=0; gat<nat; ++gat )
    {
      int const begin = GetAtomBegin(gat);
      int const end   = GetAtomEnd(gat);
      for ( int ipt=begin; ipt<end; ++ipt )
  	{

#endif
	  
	  w.assign(nat,1.);
	  double const * rq = GetQuadCrdPtr(ipt);
	  for ( int iat=0; iat<nat; ++iat )
	    {
	      double const * ra = GetAtomCrdPtr(iat);
	      double const t[3] = { rq[0]-ra[0], rq[1]-ra[1], rq[2]-ra[2] };
	      d[iat] = std::sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
	    };
	  
	  for ( int jat=1; jat<nat; ++jat )
	    for ( int iat=0; iat<jat; ++iat )
	      {
		double const rij = RIJM[ iat+jat*(jat-1)/2 ];
		double mu = rij > 1.e-8 ? ( d[iat] - d[jat] ) / rij : 1.;
		double p = ( 1.5 - 0.5 * mu*mu ) * mu;
		p *= ( 1.5 - 0.5 * p*p );
		p *= ( 1.5 - 0.5 * p*p );
		w[iat] *= ( 0.5 - 0.5*p );
		w[jat] *= ( 0.5 + 0.5*p );
	      };
	  
	  double wt = 0.;
	  for ( int iat=0; iat<nat; ++iat )
	    wt += w[iat];
	  wt = w[gat]/wt;
	  
	  SetPartitionWt(ipt,wt);
	  
	};
      
    };
}

 

