#ifndef _InternalCrd_hpp_
#define _InternalCrd_hpp_

#include <tr1/memory>
#include <vector>
#include <ostream>

namespace ccdl
{
  class InternalCrd
  {
  public:
    InternalCrd() : mFrozen(false), mConVal(0.) {};
    virtual ~InternalCrd() {};

    void Unfreeze() { mFrozen = false; }
    void Freeze() { mFrozen = true; }
    void Freeze( double const x );
    void Freeze( int const nat, double const * c );
    virtual bool Freeze( ccdl::InternalCrd const * p ) =0;

    virtual double CptValue( int const nat, double const * c ) const =0;

    virtual void CptWilson( int const nat, double const * c, double * B ) const =0;
    virtual void CptHessian( int const nat, double const * c, double * H ) const =0;

    virtual double CptDifference( double const x1, double const x0 ) const =0;
    virtual double MoveValueToValidRange( double x ) const =0;
    virtual std::string Print() const =0;
    virtual std::string PrintPrettyValue( double a ) const =0;

    double CptConstraintDelta( int const nat, double const * c ) const;
    bool IsFrozen() const { return mFrozen; }
    bool IsFree() const { return ( ! IsFrozen() ); }
    double GetConstraintValue() const { return mConVal; }
  protected:
    bool mFrozen;
    double mConVal;
  };
  typedef std::tr1::shared_ptr< ccdl::InternalCrd > pIC;

  class BondT : public InternalCrd
  {
  public:

    using InternalCrd::Freeze;

    BondT() : InternalCrd(), i(0), j(0) {};
    BondT( int const i, int const j ) : InternalCrd(), i(i), j(j) {};
    virtual ~BondT() {};
    virtual std::string Print() const;
    virtual std::string PrintPrettyValue( double a ) const;
    virtual bool Freeze( ccdl::InternalCrd const * p );
    virtual double CptValue( int const nat, double const * c ) const;
    virtual void CptWilson( int const nat, double const * c, double * B ) const;
    virtual void CptHessian( int const nat, double const * c, double * B ) const;
    virtual double CptDifference( double const x1, double const x0 ) const;
    virtual double MoveValueToValidRange( double x ) const;
    int i,j;
  };

  class AngleT : public InternalCrd
  {
  public:

    using InternalCrd::Freeze;

    AngleT() : InternalCrd(), i(0), j(0), k(0) {};
    AngleT( int const i, int const j, int const k ) : InternalCrd(), i(i), j(j), k(k) {};
    virtual ~AngleT() {};
    virtual std::string Print() const;
    virtual std::string PrintPrettyValue( double a ) const;
    virtual double CptValue( int const nat, double const * c ) const;
    virtual bool Freeze( ccdl::InternalCrd const * p );
    virtual void CptWilson( int const nat, double const * c, double * B ) const;
    virtual void CptHessian( int const nat, double const * c, double * B ) const;
    virtual double CptDifference( double const x1, double const x0 ) const;
    virtual double MoveValueToValidRange( double x ) const;
    int i,j,k;
  };

  class DihedralT : public InternalCrd
  {
  public:

    using InternalCrd::Freeze;

    DihedralT() : InternalCrd(), i(0), j(0), k(0), l(0) {};
    DihedralT( int const i, int const j, int const k, int const l ) : InternalCrd(), i(i), j(j), k(k), l(l) {};
    virtual ~DihedralT() {};
    virtual std::string Print() const;
    virtual std::string PrintPrettyValue( double a ) const;
    virtual bool Freeze( ccdl::InternalCrd const * p );
    virtual double CptValue( int const nat, double const * c ) const;
    virtual void CptWilson( int const nat, double const * c, double * B ) const;
    virtual void CptHessian( int const nat, double const * c, double * B ) const;
    virtual double CptDifference( double const x1, double const x0 ) const;
    virtual double MoveValueToValidRange( double x ) const;
    int i,j,k,l;
  };


  class RedundantIC
  {
  public:

    RedundantIC( int const nat, int const * z, double const * crd, bool merge_fragments = true );

    void reset( int const nat, int const * z, double const * crd, bool merge_fragments = true );


    void FreezeBond( int i, int j, double const v );
    void FreezeAngle( int i, int j, int k, double const v );
    void FreezeDihedral( int i, int j, int k, int l, double const v );

    void PrintReport( std::ostream & cout );
    void PrintReport( std::ostream & cout, double const * crd );
    void PrintCrd( std::ostream & cout, int i );
    void PrintPrettyValue( std::ostream & cout, int i, double const * crd );

    int  GetNumInternalCrds() const { return qs.size(); }
    int  GetNumAtoms() const { return nat; }

    void DisplaceByDeltaQ( double * dq, double * crd0, double const TOL=1.e-7 );

    void GrdTransform( double const * crd, double * cg, double * qg, bool c2q_AND_q2c = false );

    void GrdAndHesTransform( double const * crd, double * cg, double * ch, double * qg, double * qh, bool c2q_AND_q2c = false );

    void HesBackTransform( double const * crd,  double const * qg, double const * qh, double * ch );


    void CptDifference( double const * hi, double const * lo, double * diff );

    void CptInternalCrds( double const * crd, double * q );

    void CptTransformData( double const * crd, double * B, double * G, double * Ginv );

    int  GetConstraintMask( double * C );

    void CptProjector( double const * G, double const * Ginv, double * P );

    int nat;
    std::vector< ccdl::pIC > qs;
  };

}

#endif

