#ifndef _AutoBondFrag_hpp_
#define _AutoBondFrag_hpp_

#include <vector>
#include <algorithm>
#include <tr1/memory>


namespace ccdl
{
  enum CRDSYS { CART, CART_AND_BONDS, DLC, DLC_EXCEPT_BENDS };
}


namespace ccdl
{
  struct BondFragT
  {
    BondFragT();
    void push_back( int a ) { bond.push_back(a); }
    std::vector<int>::iterator begin() { return bond.begin(); }
    std::vector<int>::iterator end()   { return bond.end(); }
    int size() { return bond.size(); }
    void sort();
    bool HasHBond() const;
    
    int idx,Z,frag;
    double covR,vdwR,R;
    std::vector<int> bond;
    int hbond_donor, hbond_acceptor;
    double hbond_acceptor_R;
  };
  
  struct AutoBondFrag
  {
    AutoBondFrag
    ( int const nat, int const * z, double const * crd, 
      ccdl::CRDSYS crdsys,
      bool const merge_fragments );

    void reset
    ( int const nat, int const * z, double const * crd, 
      ccdl::CRDSYS crdsys,
      bool const merge_fragments );

    BondFragT & operator[] ( int i ) { return atom[i]; }


    int nat,nfrag;
    std::vector<BondFragT> atom;

    int nbonds;
    std::vector<int> bonds;
    int nangles;
    std::vector<int> angles;
    int ntorsions;
    std::vector<int> torsions;
    int nhbonds;
    std::vector<int> hbonds;
    int ncarts;
    std::vector<int> carts;

    void PrintReport( std::ostream & cout );
  };

  typedef std::tr1::shared_ptr<ccdl::AutoBondFrag> pAtomBondFrag;
}

#endif
