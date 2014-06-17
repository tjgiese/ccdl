#include "MolQuad.hpp"
#include "Quadrature.hpp"
#include <cassert>


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
  ccdl::GaussLaguerreRule( 0., numShells, mRadialPts, mRadialWts );
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
  for ( int i=irad; i<GetNumRadialShells(); ++i )
    mNumAngPts[i] = numAngPts;
  mNumPts = 0;
  for ( int i=0; i<GetNumRadialShells(); ++i )
    mNumPts += mNumAngPts[i];
  BuildShellOffsets();
}


void ccdl::AtomQuadParam::BuildShellOffsets()
{
  mShellOffsets[0] = 0;
  for ( int irad=0; irad<GetNumRadialShells(); ++irad )
    mShellOffsets[irad+1] = mShellOffsets[irad] + GetNumAngPts(irad);
}






ccdl::MolQuad::MolQuad
( std::vector< std::tr1::shared_ptr< ccdl::AtomQuadParam > > const & atomParam, 
  std::vector< std::tr1::array<double,3> > const & atomCrd )
  : mAtomParam( atomParam ),
    mAtomCrd( atomCrd ),
    mAtomOffset( atomCrd.size()+1, 0 ),
    mNumAtoms( atomCrd.size() )
{
  BuildGrid();
}



// ccdl::MolQuad::MolQuad
// ( std::vector< std::tr1::shared_ptr< ccdl::AtomQuadParam > > const & atomParam, 
//   std::vector< std::tr1::array<double,3> > const & atomCrd,
//   double const bufferRadius )
//   : mAtomParam( atomParam ),
//     mAtomCrd( atomCrd ),
//     mAtomOffset( atomCrd.size()+1, 0 ),
//     mNumAtoms( atomCrd.size() ),
//     mBufferRadius( bufferRadius )
// {
//   BuildGrid();
//   BuildNeighborList();
// }







void ccdl::MolQuad::BuildGrid()
{
  assert( mNumAtoms == mAtomParam.size() );
  assert( mNumAtoms+1 == mAtomOffset.size() );

  mAtomOffset[0] = 0;
  for ( int iat=0; iat<GetNumAtoms(); ++iat )
    mAtomOffset[iat+1] = mAtomOffset[iat] + GetNumPts(iat);
  
  mNumPts = mAtomOffset.back();
  mSingleCenterWt.resize( GetNumPts() );
  mPartitionWt.resize( GetNumPts() );
  mTotalWt.resize( GetNumPts() );
  mQuadCrd.resize( GetNumPts() );

  for ( int i=0; i<GetNumPts(); ++i )
    mPartitionWt[i] = 1.;

  std::vector<double> angPts,angWts;
  ccdl::LebedevRule( GetNumPts(0,0), angPts, angWts );
  for ( int iat=0; iat < GetNumAtoms(); ++iat )
    {
      for ( int irad=0, nrad=GetNumRadialShells(iat); irad<nrad; ++irad )
	{
	  int const nang = GetNumPts(iat,irad);

	  if ( nang != angWts.size() )
	    ccdl::LebedevRule( nang, angPts, angWts );

	  int const o = GetRadialShellOffset(iat,irad);
	  double const radius = GetAtomParam(iat).GetRadialPt(irad);
	  double const radWt  = GetAtomParam(iat).GetRadialWt(irad);
	  for ( int iang=0; iang < nang; ++iang )
	    {
	      mSingleCenterWt[iang+o] = radWt * angWts[iang];
	      mQuadCrd[iang+o][0] = mAtomCrd[iat][0] + radius * angPts[0+iang*3];
	      mQuadCrd[iang+o][1] = mAtomCrd[iat][1] + radius * angPts[1+iang*3];
	      mQuadCrd[iang+o][2] = mAtomCrd[iat][2] + radius * angPts[2+iang*3];
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


// void ccdl::MolQuad::BuildNeighborList()
// {
//   int const nat = GetNumAtoms();
//   mNeighborList.resize( nat );

//   int const MAXNUM = 10;

//   std::vector< std::vector<ccdl::backend::pair> > pairs(nat);

//   for ( int iat=0; iat<nat; ++iat )
//     {
//       mNeighborList[iat].resize(0);
//       mNeighborList[iat].reserve(MAXNUM);
//       pairs[iat].resize(nat);
//     };
// #ifdef _OPENMP
// #pragma omp parallel for schedule(dynamic)
// #endif
//   for ( int jat=1; jat<nat; ++jat )
//     {
//       for ( int iat=0; iat<jat; ++iat )
// 	{
// 	  double const d[3] = { GetAtomCrd(iat,0)-GetAtomCrd(jat,0),
// 				GetAtomCrd(iat,1)-GetAtomCrd(jat,1),
// 				GetAtomCrd(iat,2)-GetAtomCrd(jat,2) };
// 	  double const rij2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
// 	  pairs[iat][jat].atom = jat;
// 	  pairs[iat][jat].rij  = rij2;
// 	  pairs[jat][iat].atom = iat;
// 	  pairs[jat][iat].rij  = rij2;
// 	};
//     };
// #ifdef _OPENMP
// #pragma omp parallel for schedule(dynamic)
// #endif
//   for ( int iat=0; iat<nat; ++iat )
//     {
//       pairs[iat].erase( pairs[iat].begin()+iat );
//       std::sort( pairs[iat].begin(), pairs[iat].end() );
//       int const n = std::min( pairs[iat].size(), MAXNUM );
//       for ( int i=0; i<n; ++i )
// 	mNeighborList[iat].push_back( pairs[iat][i].atom );
//       std::sort( mNeighborList[iat].begin(), mNeighborList[iat].end() );
//     };
// }




void ccdl::MolQuad::SetAtomCrd
( int const iat, double const * c )
{
  double const d[3] = { c[0]-mAtomCrd[iat][0], c[1]-mAtomCrd[iat][1], c[2]-mAtomCrd[iat][2] };
  mAtomCrd[iat][0] = c[0];
  mAtomCrd[iat][1] = c[1];
  mAtomCrd[iat][2] = c[2];

  for ( int ipt=GetAtomOffset(iat), last=GetAtomOffset(iat+1); ipt<last; ++ipt )
    {
      mQuadCrd[ipt][0] += d[0];
      mQuadCrd[ipt][1] += d[1];
      mQuadCrd[ipt][2] += d[2];
    };
  // if ( HasNeighborList() )
  //   BuildNeighborList();
}




// void ccdl::MolQuad::EvaluateBeckeWeights_WithoutRadiiBruteForce()
// {
//   int const nat = GetNumAtoms();

  
//   std::vector<double> w( nat, 0. );
//   for ( int gat=0; gat<nat; ++gat )
//     {
//       for ( int ipt=GetAtomOffset(gat); ipt<GetAtomOffset(gat+1); ++ipt )
//   	{

//   // int const npt = GetNumPts();
//   // std::vector<double> w( nat, 0. );
//   // for ( int ipt=0; ipt<npt; ++ipt )
//   //   {
//   //     int const gat = GetAtomIdx(ipt);


//       for ( int iat=0; iat<nat; ++iat )
// 	w[iat] = 1.;
      
//       for ( int iat=0; iat<nat; ++iat )
// 	for ( int jat=0; jat<iat; ++jat )
// 	  {
// 	    double t[3] = { GetAtomCrd(iat,0)-GetAtomCrd(jat,0),
// 			    GetAtomCrd(iat,1)-GetAtomCrd(jat,1),
// 			    GetAtomCrd(iat,2)-GetAtomCrd(jat,2) };
// 	    double const rij = std::sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
// 	    t[0] = GetQuadCrd(ipt,0)-GetAtomCrd(iat,0);
// 	    t[1] = GetQuadCrd(ipt,1)-GetAtomCrd(iat,1);
// 	    t[2] = GetQuadCrd(ipt,2)-GetAtomCrd(iat,2);
// 	    double const rpi = std::sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
// 	    t[0] = GetQuadCrd(ipt,0)-GetAtomCrd(jat,0);
// 	    t[1] = GetQuadCrd(ipt,1)-GetAtomCrd(jat,1);
// 	    t[2] = GetQuadCrd(ipt,2)-GetAtomCrd(jat,2);
// 	    double const rpj = std::sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
// 	    double const mu = ( rpi - rpj ) / rij;
// 	    double p = ( 1.5 - 0.5 * mu*mu ) * mu;
// 	    p *= ( 1.5 - 0.5 * p*p );
// 	    p *= ( 1.5 - 0.5 * p*p );
// 	    w[iat] *= ( 0.5 - 0.5*p );
// 	    w[jat] *= ( 0.5 + 0.5*p );
// 	  };
      
//       double wt = 0.;
//       for ( int iat=0; iat<nat; ++iat )
// 	wt += w[iat];
//       wt = w[gat]/wt;
      
//       SetPartitionWt(ipt,wt);
      
//     };

//     };
// }





void ccdl::MolQuad::EvaluateBeckeWeights_WithRadiiBruteForce()
{
  int const nat = GetNumAtoms();


  std::vector<double> RIJM(nat*(nat-1),0.);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for ( int jat=1; jat<nat; ++jat )
    for ( int iat=0; iat<jat; ++iat )
      {
	double t[3] = { GetAtomCrd(iat,0)-GetAtomCrd(jat,0),
			GetAtomCrd(iat,1)-GetAtomCrd(jat,1),
			GetAtomCrd(iat,2)-GetAtomCrd(jat,2) };
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
      for ( int ipt=GetAtomOffset(gat); ipt<GetAtomOffset(gat+1); ++ipt )
  	{

#endif

	  for ( int iat=0; iat<nat; ++iat )
	    w[iat] = 1.;
      
	  for ( int iat=0; iat<nat; ++iat )
	    {
	      double const t[3] = { GetQuadCrd(ipt,0)-GetAtomCrd(iat,0), 
				    GetQuadCrd(ipt,1)-GetAtomCrd(iat,1),
				    GetQuadCrd(ipt,2)-GetAtomCrd(iat,2) };
	      d[iat] = std::sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
	    };
	  
	  for ( int jat=1; jat<nat; ++jat )
	    for ( int iat=0; iat<jat; ++iat )
	      {
		int const k = 2*iat+jat*(jat-1);
		double const rij = RIJM[0+k];
		double const a   = RIJM[1+k];
		double mu = ( d[iat] - d[jat] ) / rij;
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
    for ( int iat=0; iat<jat; ++iat )
      {
	double t[3] = { GetAtomCrd(iat,0)-GetAtomCrd(jat,0),
			GetAtomCrd(iat,1)-GetAtomCrd(jat,1),
			GetAtomCrd(iat,2)-GetAtomCrd(jat,2) };
	RIJM[ iat+jat*(jat-1)/2 ] = std::sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
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
      for ( int ipt=GetAtomOffset(gat); ipt<GetAtomOffset(gat+1); ++ipt )
  	{

#endif

	  for ( int iat=0; iat<nat; ++iat )
	    w[iat] = 1.;
      
	  for ( int iat=0; iat<nat; ++iat )
	    {
	      double const t[3] = { GetQuadCrd(ipt,0)-GetAtomCrd(iat,0), 
				    GetQuadCrd(ipt,1)-GetAtomCrd(iat,1),
				    GetQuadCrd(ipt,2)-GetAtomCrd(iat,2) };
	      d[iat] = std::sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
	    };
	  
	  for ( int jat=1; jat<nat; ++jat )
	    for ( int iat=0; iat<jat; ++iat )
	      {
		double const rij = RIJM[ iat+jat*(jat-1)/2 ];
		double const mu = ( d[iat] - d[jat] ) / rij;
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

 



// void ccdl::MolQuad::EvaluateBeckeWeights_WithoutRadiiNearestN( int const maxN )
// {
//   int const nat = GetNumAtoms();

//   if ( nat < 1 ) { return; }

//   int const MAXN = std::min( nat, maxN );

//   std::vector<double> RIJM(nat*(nat-1)/2,0.);

// #ifdef _OPENMP
// #pragma omp parallel for schedule(dynamic)
// #endif
//   for ( int jat=1; jat<nat; ++jat )
//     for ( int iat=0; iat<jat; ++iat )
//       {
// 	double t[3] = { GetAtomCrd(iat,0)-GetAtomCrd(jat,0),
// 			GetAtomCrd(iat,1)-GetAtomCrd(jat,1),
// 			GetAtomCrd(iat,2)-GetAtomCrd(jat,2) };
// 	RIJM[ iat+jat*(jat-1)/2 ] = std::sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
//       };

// #ifdef _OPENMP
// #pragma omp parallel
//   {
// #endif
//   std::vector<double> w( nat, 0. );
//   std::vector< ccdl::backend::pair > d( nat );

// #ifdef _OPENMP

//   int const npt = GetNumPts();
// #pragma omp for schedule(dynamic)
//   for ( int ipt=0; ipt<npt; ++ipt )
//     {
//       int const gat = GetAtomIdx(ipt);

// #else

//   for ( int gat=0; gat<nat; ++gat )
//     {
//       for ( int ipt=GetAtomOffset(gat); ipt<GetAtomOffset(gat+1); ++ipt )
//   	{

// #endif

// 	  for ( int iat=0; iat<nat; ++iat )
// 	    w[iat] = 1.;
      
// 	  for ( int iat=0; iat<nat; ++iat )
// 	    {
// 	      double const t[3] = { GetQuadCrd(ipt,0)-GetAtomCrd(iat,0), 
// 				    GetQuadCrd(ipt,1)-GetAtomCrd(iat,1),
// 				    GetQuadCrd(ipt,2)-GetAtomCrd(iat,2) };
// 	      d[iat].atom = iat;
// 	      d[iat].rij = std::sqrt( t[0]*t[0] + t[1]*t[1] + t[2]*t[2] );
// 	    };
	  
// 	  std::sort( d.begin(), d.end() );
// 	  std::sort( d.begin(), d.begin()+MAXN, ccdl::backend::pair::AtomLessThan );


// 	  for ( int j=1; j<MAXN; ++j )
// 	    {
// 	      int const jat = d[j].atom;
// 	      int const o = jat*(jat-1)/2;
// 	      for ( int i=0; i<j; ++i )
// 		{
// 		  int const iat = d[i].atom;
// 		  double const rij = RIJM[ iat+o ];
// 		  double const mu = ( d[i].rij - d[j].rij ) / rij;
// 		  double p = ( 1.5 - 0.5 * mu*mu ) * mu;
// 		  p *= ( 1.5 - 0.5 * p*p );
// 		  p *= ( 1.5 - 0.5 * p*p );
// 		  w[iat] *= ( 0.5 - 0.5*p );
// 		  w[jat] *= ( 0.5 + 0.5*p );
// 		};
// 	    };
	  
// 	  double wt = 0.;
// 	  for ( int iat=0; iat<nat; ++iat )
// 	    wt += w[iat];
// 	  wt = w[gat]/wt;
	  
// 	  SetPartitionWt(ipt,wt);
	  
// 	};
      
//     };
// }

