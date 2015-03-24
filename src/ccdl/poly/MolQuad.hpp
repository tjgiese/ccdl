#ifndef _ccdl_MolQuad_H_
#define _ccdl_MolQuad_H_

#include <tr1/array>
#include <tr1/memory>
#include <vector>
#include <algorithm>

namespace ccdl
{


  class AtomQuadParam
  {
  public:

    AtomQuadParam( int const numShells, int const numAngPts, double const quadRadius, double const atomRadius );

    int GetNumPts() const { return mNumPts; }

    int GetNumRadialShells() const { return mNumAngPts.size(); }

    double GetAtomRadius() const { return mAtomRadius; }

    double GetQuadratureRadius() const { return mRadialPts.back(); }

    double GetRadialPt( int const irad ) const { return mRadialPts[irad]; }

    double const * GetRadialPts() const { return mRadialPts.data(); }

    double GetRadialWt( int const irad ) const { return mRadialWts[irad]; }

    int GetNumAngPts( int const irad ) const { return mNumAngPts[irad]; }

    void Prune( int const irad, int const numAngPts );

    int GetRadialShellBegin( int const irad ) const { return mShellOffsets[irad]; }
    int GetRadialShellEnd( int const irad ) const { return mShellOffsets[irad+1]; }

    //int GetAngRuleType( int irad ) const { return mAngRuleType[irad]; };

  private:
    
    void BuildShellOffsets();


  private:

    std::vector<int> mNumAngPts;
    //std::vector<int> mAngRuleType;
    double mAtomRadius;
    std::vector<double> mRadialPts;
    std::vector<double> mRadialWts;
    int mNumPts;
    std::vector<int> mShellOffsets;
  };




  class MolQuad
  {
  public:

    template <class T>
    MolQuad( std::vector< std::tr1::shared_ptr< ccdl::AtomQuadParam > > const & atomParam, 
	     T const & atomCrd );
   
    int GetNumAtoms() const { return mNumAtoms; }
    int GetNumPts() const { return mNumPts; }
    int GetNumAtomPts( int const iat ) const { return GetAtomParam(iat).GetNumPts(); }
    int GetNumRadialShells( int const iat ) const { return GetAtomParam(iat).GetNumRadialShells(); }
    int GetNumAngPts( int const iat, int const irad ) const { return GetAtomParam(iat).GetNumAngPts(irad); }

    int GetAtomBegin( int const iat ) const { return mAtomOffsets[iat]; }
    int GetAtomEnd( int const iat ) const { return mAtomOffsets[iat+1]; }

    int GetRadialShellBegin( int const iat, int const irad ) const { return mAtomOffsets[iat] + GetAtomParam(iat).GetRadialShellBegin(irad); }
    int GetRadialShellEnd( int const iat, int const irad ) const { return mAtomOffsets[iat] + GetAtomParam(iat).GetRadialShellEnd(irad); }


    double const * GetQuadCrdPtr( int const ipt ) const { return mQuadCrd.data() + ipt*3; }
    double const * GetAtomCrdPtr( int const iat ) const { return mAtomCrd.data() + iat*3; }
    void SetAtomCrd( int const iat, double const * crd );

    double GetQuadWt( int const ipt ) const { return mTotalWt[ipt]; }
    double GetSingleCenterWt( int const ipt ) const { return mSingleCenterWt[ipt]; }
    double GetPartitionWt( int const ipt ) const { return mPartitionWt[ipt]; }
    void SetPartitionWt( int const ipt, double const w ) { mPartitionWt[ipt] = w; mTotalWt[ipt] = mSingleCenterWt[ipt] * w; }
    void ClearPartitionWts() { for ( int i=0,n=GetNumPts(); i<n; ++i ) SetPartitionWt(i,1.); }
    void EvaluateBeckeWeights_WithRadiiBruteForce();
    void EvaluateBeckeWeights_WithoutRadiiBruteForce();

    double GetAtomRadius( int const iat ) const { return GetAtomParam(iat).GetAtomRadius(); }
    double GetQuadratureRadius( int const iat ) const { return GetAtomParam(iat).GetQuadratureRadius(); }

    int GetAtomIdx( int const ipt ) const { return std::distance( mAtomOffsets.begin()+1, std::upper_bound( mAtomOffsets.begin()+1, mAtomOffsets.end(), ipt ) ); }
    
  private:


    std::vector<double> mSingleCenterWt;
    std::vector<double> mPartitionWt;
    std::vector<double> mTotalWt;
    std::vector<double> mQuadCrd;
    //std::vector< std::tr1::array<double,3> > mQuadCrd;
    int mNumPts;

    void BuildGrid();


  private:

    ccdl::AtomQuadParam const & GetAtomParam( int const iat ) const { return *(mAtomParam[iat]); }

    std::vector< std::tr1::shared_ptr< ccdl::AtomQuadParam > > mAtomParam;
    std::vector<double> mAtomCrd;
    //std::vector< std::tr1::array<double,3> > mAtomCrd;
    std::vector<int> mAtomOffsets;
    int mNumAtoms;


  };


}


// inline void ccdl::MolQuad::GetCommonNeighbors
// ( int const iat, int const jat, 
//   std::vector<int> & atoms ) const 
// { 
//   assert( HasNeighborList() ); 
//   atoms.resize(0); 
//   std::set_intersection( mNeighborList[iat].begin(), mNeighborList[iat].end(), 
// 			 mNeighborList[jat].begin(), mNeighborList[jat].end(), 
// 			 std::back_inserter(atoms) ); 
// }



template <class T>
ccdl::MolQuad::MolQuad
( std::vector< std::tr1::shared_ptr< ccdl::AtomQuadParam > > const & atomParam, 
  T const & atomCrd )
  : mAtomParam( atomParam ),
    mAtomCrd( 3*atomCrd.size() ),
    mAtomOffsets( atomCrd.size()+1, 0 ),
    mNumAtoms( atomCrd.size() )
{
  int nat = atomCrd.size();
  for ( int a=0; a<nat; ++a )
    {
      mAtomCrd[0+a*3] = atomCrd[a][0];
      mAtomCrd[1+a*3] = atomCrd[a][1];
      mAtomCrd[2+a*3] = atomCrd[a][2];
    };
  BuildGrid();
}


#endif
