#ifndef _GaussianBasis_H_
#define _GaussianBasis_H_

#include <vector>
#include <string>

namespace GaussianBasis
{

  struct Segment
  {
  public:

    Segment() : mL(0) {}

    Segment( int const L, 
	     std::vector<double> const & C, 
	     std::vector<double> const & Z ) 
      : mL(L), mC(C), mZ(Z) {};
    
    Segment( int const L, 
	     int const n,
	     double const * C, 
	     double const * Z ) 
      : mL(L), mC(C,C+n), mZ(Z,Z+n) {}

    int GetNumPrim() const { return mC.size(); }
    int GetL() const { return mL; }
    double const * GetC() const { return mC.data(); }
    double const * GetZ() const { return mZ.data(); }
    int GetFockSize() const { return 2*GetL()+1; }

  private:

    int mL;
    std::vector<double> mC;
    std::vector<double> mZ;

  };



  class Basis
  {
  public:

    Basis( char const * name, int const atNum );

    int GetFockSize() const { return mFockSize; }

    int GetNumSegments() const { return mSeg.size(); }

    int GetLmax() const { return mLmax; }

    int GetL( int const iseg ) const { return mSeg[iseg].GetL(); }

    double const * GetC( int const iseg ) const { return mSeg[iseg].GetC(); }

    double const * GetZ( int const iseg ) const { return mSeg[iseg].GetZ(); }

    GaussianBasis::Segment const & GetSegment( int const iseg ) const { return mSeg[iseg]; }

    GaussianBasis::Segment & GetSegment( int const iseg ) { return mSeg[iseg]; }

    void PushSegment( GaussianBasis::Segment const & s ) 
    { 
      mSeg.push_back(s); 
      mFockSize += s.GetFockSize();
      mLmax = std::max( mLmax, s.GetL() );
    }


  protected:
    
    std::string mType;
    int mAtomId;
    int mFockSize;
    int mLmax;
    std::vector< GaussianBasis::Segment > mSeg;

  };


 
}

#endif
