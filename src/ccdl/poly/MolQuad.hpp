#ifndef _ccdl_MolQuad_H_
#define _ccdl_MolQuad_H_

#include <tr1/array>
#include <tr1/memory>
#include <cstddef>
#include <vector>
#include <algorithm>
#include <cassert>

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

    double GetRadialWt( int const irad ) const { return mRadialWts[irad]; }

    int GetNumAngPts( int const irad ) const { return mNumAngPts[irad]; }

    void Prune( int const irad, int const numAngPts );

    int GetRadialShellOffset( int const irad ) const { return mShellOffsets[irad]; }


  private:
    
    void BuildShellOffsets();


  private:

    std::vector<int> mNumAngPts;
    double mAtomRadius;
    std::vector<double> mRadialPts;
    std::vector<double> mRadialWts;
    int mNumPts;
    std::vector<int> mShellOffsets;
  };




  class MolQuad
  {
  public:

    typedef std::vector< std::tr1::array<double,3> >::iterator iterator;
    typedef std::vector< std::tr1::array<double,3> >::const_iterator const_iterator;

    MolQuad( std::vector< std::tr1::shared_ptr< ccdl::AtomQuadParam > > const & atomParam, 
	     std::vector< std::tr1::array<double,3> > const & atomCrd );
   
    MolQuad( std::vector< std::tr1::shared_ptr< ccdl::AtomQuadParam > > const & atomParam, 
	     std::vector< std::tr1::array<double,3> > const & atomCrd,
	     double const bufferRadius );
    

    double GetQuadWt( int const ipt ) const { return mTotalWt[ipt]; }

    double const * GetQuadCrdPtr( int const ipt ) const { return mQuadCrd[ipt].data(); }

    double GetQuadCrd( int const ipt, int const k ) const { return mQuadCrd[ipt][k]; }

    const_iterator GetElementBegin() const { return mQuadCrd.begin(); }
    const_iterator GetElementEnd() const { return mQuadCrd.end(); }

    iterator GetElementBegin() { return mQuadCrd.begin(); }
    iterator GetElementEnd() { return mQuadCrd.end(); }


    double GetSingleCenterWt( int const ipt ) const { return mSingleCenterWt[ipt]; }

    double GetPartitionWt( int const ipt ) const { return mPartitionWt[ipt]; }

    void SetPartitionWt( int const ipt, double const w ) { mPartitionWt[ipt] = w; mTotalWt[ipt] = mSingleCenterWt[ipt] * w; }

    void ClearPartitionWts() { for ( int i=0,n=GetNumPts(); i<n; ++i ) SetPartitionWt(i,1.); }

    void EvaluateBeckeWeights_WithRadiiBruteForce();

    void EvaluateBeckeWeights_WithoutRadiiBruteForce();

    //void EvaluateBeckeWeights_WithoutRadiiNearestN( int const maxN );

    int GetNumAtoms() const { return mNumAtoms; }

    double const * GetAtomCrdPtr( int const iat ) const { return mAtomCrd[iat].data(); }

    double GetAtomCrd( int const iat, int const k ) const { return mAtomCrd[iat][k]; }

    void SetAtomCrd( int const iat, double const * crd );

    double GetAtomRadius( int const iat ) const { return GetAtomParam(iat).GetAtomRadius(); }

    double GetQuadratureRadius( int const iat ) const { return GetAtomParam(iat).GetQuadratureRadius(); }

    int GetAtomOffset( int const iat ) const { return mAtomOffset[iat]; }

    int GetNumRadialShells( int const iat ) const { return GetAtomParam(iat).GetNumRadialShells(); }

    int GetRadialShellOffset( int const iat, int const irad ) const { return mAtomOffset[iat] + GetAtomParam(iat).GetRadialShellOffset(irad); }

    int GetNumPts() const { return mNumPts; }

    int GetNumPts( int const iat ) const { return GetAtomParam(iat).GetNumPts(); }

    int GetNumPts( int const iat, int const irad ) const { return GetAtomParam(iat).GetNumAngPts(irad); }

    int GetAtomIdx( int const ipt ) const { return std::distance( mAtomOffset.begin()+1, std::upper_bound( mAtomOffset.begin()+1, mAtomOffset.end(), ipt ) ); }
    
    void GetDifferenceCrd( int const iat, int const ipt, double * crd ) 
    { crd[0] = GetQuadCrd(ipt,0)-GetAtomCrd(iat,0); 
      crd[1] = GetQuadCrd(ipt,1)-GetAtomCrd(iat,1); 
      crd[2] = GetQuadCrd(ipt,2)-GetAtomCrd(iat,2); }


    // bool HasNeighborList() const { return ( mNeighborList.size() > 0 ); }
    
    // int GetNumNeighbors( int const iat ) const { return mNeighborList[iat].size(); }

    // std::vector<int> const & GetNeighbors( int const iat ) const { return mNeighborList[iat]; }

    // void GetCommonNeighbors( int const iat, int const jat, std::vector<int> & atoms ) const;


  private:


    std::vector<double> mSingleCenterWt;
    std::vector<double> mPartitionWt;
    std::vector<double> mTotalWt;
    std::vector< std::tr1::array<double,3> > mQuadCrd;
    int mNumPts;

    void BuildGrid();


  private:

    ccdl::AtomQuadParam const & GetAtomParam( int const iat ) const { return *(mAtomParam[iat]); }

    std::vector< std::tr1::shared_ptr< ccdl::AtomQuadParam > > mAtomParam;
    std::vector< std::tr1::array<double,3> > mAtomCrd;
    std::vector<int> mAtomOffset;
    int mNumAtoms;


  // private:

  //   void BuildNeighborList();

  //   double GetBufferedRadius( int const iat ) const { return GetQuadratureRadius(iat)+mBufferRadius; }

  //   double mBufferRadius;
  //   std::vector< std::vector<int> > mNeighborList;

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




#endif
