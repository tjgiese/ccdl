#include "EuclideanGeometry.hpp"
#include "AutoBondFrag.hpp"
#include "../constants.hpp"
#include "../periodictable.hpp"

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <tr1/array>
#include <cmath>


ccdl::BondFragT::BondFragT() 
  : idx(0),Z(0),frag(-1),
    covR(0.),vdwR(0.),R(0.),
    hbond_donor(-1),hbond_acceptor(-1),hbond_acceptor_R(1000.)
{}

void ccdl::BondFragT::sort()
{
  std::sort( begin(), end() );
  bond.resize( std::unique( begin(), end() ) - begin() );
}

bool ccdl::BondFragT::HasHBond() const
{
  return Z==1 and hbond_donor >= 0 and hbond_acceptor >= 0;
}


namespace ccdl
{
  namespace work
  {
    bool SetFragment( std::vector<BondFragT> & bf, int a, int frag )
    {
      bool set = false;
      typedef std::vector<int>::iterator iter;
      if ( bf[a].frag < 0 )
	{
	  set = true;
	  bf[a].frag = frag;
	  for ( iter p=bf[a].begin(), e=bf[a].end(); p!=e; ++p )
	    ccdl::work::SetFragment( bf, *p, frag );
	}
      return set;
    }
  }
}


ccdl::AutoBondFrag::AutoBondFrag
( int const natom, int const * z, double const * crd, 
  bool const merge_fragments )
{
  reset( natom,z,crd,merge_fragments );
}

void ccdl::AutoBondFrag::reset
( int const natom, int const * z, double const * crd, 
  bool const merge_fragments )
{
  nbonds=0;
  nangles=0;
  ntorsions=0;
  nhbonds=0;
  nfrag = 0;
  nat = natom;
  atom.resize(0);
  atom.resize(nat);
  for ( int a=0; a<nat; ++a )
    {
      atom[a].idx  = a;
      atom[a].Z    = z[a];
      atom[a].covR = ccdl::CovalentRadius( z[a] );
      atom[a].vdwR = ccdl::vdWRadius( z[a] );
      atom[a].R    = 0.8 * atom[a].covR + 0.2 * atom[a].vdwR;
    }
  for ( int a=1; a<nat; ++a )
    for ( int b=0; b<a; ++b )
      {
	double c[3] = { crd[0+a*3]-crd[0+b*3],
			crd[1+a*3]-crd[1+b*3],
			crd[2+a*3]-crd[2+b*3] };
	double r2 = c[0]*c[0]+c[1]*c[1]+c[2]*c[2];
	double r = std::sqrt(r2);
	double tol = 1. * ( atom[a].R + atom[b].R );
	bool check_hbond = ( ((atom[a].Z == 1) + (atom[b].Z == 1)) == 1 );
	if ( check_hbond )
	  {
	    if ( r < tol ) // this is covalent
	      {
		atom[a].push_back(b);
		atom[b].push_back(a);
	      }
	    else if ( r < tol + ccdl::AU_PER_ANGSTROM ) // this is noncovalent
	      {
		int i = (atom[a].Z == 1) ? a : b;
		int j = (atom[a].Z == 1) ? b : a;
		if ( atom[j].Z ==  7 or atom[j].Z ==  8 or 
		     atom[j].Z ==  9 or atom[j].Z == 16 or 
		     atom[j].Z == 17 or atom[j].Z == 35 )
		  if ( atom[i].hbond_acceptor < 0 or
		       ( atom[i].hbond_acceptor >= 0 and
			 r < atom[i].hbond_acceptor_R ) )
		    {
		      atom[i].hbond_acceptor = j;
		      atom[i].hbond_acceptor_R = r;
		    };
	      }
	  }
	else if ( r < tol ) // this is covalent
	  {
	    atom[a].push_back(b);
	    atom[b].push_back(a);
	  }
      }
  // find the hbond donors for each H that has an acceptor
  // throw away hbond if the hbond angle doesn't make sense
  for ( int a=0; a<nat; ++a ) 
    {
      atom[a].sort();
      if ( atom[a].Z == 1 and atom[a].hbond_acceptor >= 0 and atom[a].size() > 0 )
	{
	  atom[a].hbond_donor = atom[a].bond[0];
	  int i=atom[a].hbond_donor;
	  int j=a;
	  int k=atom[a].hbond_acceptor;
	  double ang = (180./ccdl::PI) * 
	    ccdl::Angle( crd+3*i, crd+3*j, crd+3*k );
	  if ( ang < 100. or atom[i].Z < 7 )
	    {
	      atom[a].hbond_donor = -1;
	      atom[a].hbond_acceptor = -1;
	    }
	};
    };
  // set the fragment index
  nfrag=0;
  for ( int a=0; a<nat; ++a )
    nfrag += ccdl::work::SetFragment( atom, a, nfrag );

  
  if ( merge_fragments )
    {
      // join each fragment to every other fragment
      // by the closest bond between them
      // if that bond is an hbond, then replace that hbond
      for ( int afrag=1; afrag < nfrag; ++afrag )
	{
	  std::vector<int> imin( 2*nfrag, -1 );
	  std::vector<double> rmin( nfrag, 100000. );
	  for ( int a=0; a<nat; ++a )
	    {
	      if ( atom[a].frag != afrag ) continue;
	      for ( int atomrag=0; atomrag < afrag; ++atomrag )
		for ( int b=0; b<nat; ++b )
		  {
		    if ( atom[b].frag != atomrag ) continue;

		    double c[3] = { crd[0+a*3]-crd[0+b*3],
				    crd[1+a*3]-crd[1+b*3],
				    crd[2+a*3]-crd[2+b*3] };
		    double r2 = c[0]*c[0]+c[1]*c[1]+c[2]*c[2];
		    double r = std::sqrt(r2);
		    if ( r < rmin[atomrag] )
		      {
			rmin[atomrag] = r;
			imin[0+atomrag*2] = a;
			imin[1+atomrag*2] = b;
		      };
		  }
	    };
	  for ( int atomrag=0; atomrag < afrag; ++atomrag )
	    {
	      
	      int a = imin[0+atomrag*2];
	      int b = imin[1+atomrag*2];
	      atom[a].push_back( b );
	      atom[b].push_back( a );
	      if ( atom[a].Z == 1 and atom[a].hbond_acceptor == b )
		{
		  atom[a].hbond_acceptor = -1;
		  atom[a].hbond_donor = -1;
		}
	      else if ( atom[b].Z == 1 and atom[b].hbond_acceptor == a )
		{
		  atom[b].hbond_acceptor = -1;
		  atom[b].hbond_donor = -1;
		}
	    };
	};
      for ( int a=0; a<nat; ++a )
	{
	  atom[a].frag = 0;
	  atom[a].sort();
	};
      nfrag = 1;



      // -- hack -- remove all h-bonds and insert them as regular bonds
      for ( int a=0; a<nat; ++a )
	if ( atom[a].HasHBond() )
	  {
	    atom[a].push_back( atom[a].hbond_acceptor );
	    atom[ atom[a].hbond_acceptor ].push_back( a );
	    atom[a].hbond_acceptor = -1;
	    atom[a].hbond_donor = -1;
	  }
      //
    }
  else
    { // -- hack -- delete the hbonds because we're going to use
      // com displacement and orientation to move the fragments around
      // relative to each other
      for ( int a=0; a<nat; ++a )
	if ( atom[a].HasHBond() )
	  {
	    int acc = atom[a].hbond_acceptor;
	    int don = atom[a].hbond_donor;
	    if ( atom[a].frag == atom[acc].frag and
		 atom[a].frag == atom[don].frag )
	      {
		atom[a].push_back( acc );
		atom[ acc ].push_back( a );
	      };
	    atom[a].hbond_acceptor = -1;
	    atom[a].hbond_donor = -1;
	  }
    };

  for ( int a=0; a<nat; ++a )
    atom[a].sort();

  typedef std::tr1::array<int,3> angt;
  typedef std::tr1::array<int,4> tort;
  typedef std::vector<int>::iterator iter;

  std::vector< angt > angs;
  for ( int a=0; a<nat; ++a )
    for ( iter p=atom[a].begin(), pend=atom[a].end();
	  p != pend; ++p )
      for ( iter q=atom[*p].begin(), qend=atom[*p].end();
	    q != qend; ++q )
	if ( *q > a ) 
	  {
	    angt tmp;
	    tmp[0] = a;
	    tmp[1] = *p;
	    tmp[2] = *q;
	    bool has = false;
	    for ( std::vector< angt >::iterator 
		    h=angs.begin(), hend=angs.end(); h != hend; ++h )
	      if ( (*h)[0] == tmp[2] and 
		   (*h)[1] == tmp[1] and 
		   (*h)[2] == tmp[0] )
		{
		  has = true;
		  break;
		};
	    if ( ! has )
	      angs.push_back( tmp );
	  };
  
  int nangs = angs.size();
  std::vector< tort > tors;
  for ( int iang=0; iang < nangs; ++iang )
    {
      int const i = angs[iang][0];
      int const j = angs[iang][1];
      int const k = angs[iang][2];
      for ( iter p=atom[k].begin(), pend=atom[k].end();
	      p != pend; ++p )
        {
          int const l = *p;
	  if ( l != i and l != j and l != k )
	    {
	      tort tmp;
	      tmp[0] = i;
	      tmp[1] = j;
	      tmp[2] = k;
	      tmp[3] = l;

	      bool linear = false;
	      ccdl::DihedralAngle
		( crd+i*3, crd+j*3, crd+k*3, crd+l*3, linear );
	      if ( linear ) continue;
	      double ang1 = (180./ccdl::PI)*ccdl::Angle
		( crd+i*3, crd+j*3, crd+k*3 );
	      double ang2 = (180./ccdl::PI)*ccdl::Angle
		( crd+j*3, crd+k*3, crd+l*3 );
	      if ( std::abs(ang1) > 178.5 or std::abs(ang2) > 178.5 )
		continue;

	      bool has = false;
	      for ( std::vector< tort >::iterator 
		      h=tors.begin(), hend=tors.end(); h != hend; ++h )
		if ( (*h)[0] == tmp[3] and 
		     (*h)[1] == tmp[2] and 
		     (*h)[2] == tmp[1] and 
		     (*h)[3] == tmp[0] )
		  {
		    has = true;
		    break;
		  };
	      if ( ! has )
		tors.push_back( tmp );
	    };
	};
    }
  int ntors = tors.size();

  nbonds = 0;
  bonds.resize(0);
  for ( int a=0; a<nat; ++a )
    for ( iter p=atom[a].begin(), pend=atom[a].end();
	  p != pend; ++p )
      if ( *p > a )
	{
	  bonds.push_back(a);
	  bonds.push_back(*p);
	  nbonds++;
	};

  nangles = angs.size();
  angles.resize(0);
  for ( int a=0; a<nangs; ++a )
    {
      angles.push_back( angs[a][0] );
      angles.push_back( angs[a][1] );
      angles.push_back( angs[a][2] );
    }

  ntorsions = tors.size();
  torsions.resize(0);
  for ( int a=0; a<ntors; ++a )
    {
      torsions.push_back( tors[a][0] );
      torsions.push_back( tors[a][1] );
      torsions.push_back( tors[a][2] );
      torsions.push_back( tors[a][3] );
    }

  nhbonds = 0;
  hbonds.resize(0);
  for ( int a=0; a<nat; ++a )
    if ( atom[a].HasHBond() )
      {
	hbonds.push_back( atom[a].hbond_donor );
	hbonds.push_back( a );
	hbonds.push_back( atom[a].hbond_acceptor );
	nhbonds++;
      }
}


void ccdl::AutoBondFrag::PrintReport( std::ostream & cout )
{
#define fmt std::setw(5)
  for ( int a=0; a<nbonds; ++a )
    cout << "B " 
	 << fmt << bonds[0+a*2]+1
	 << fmt << bonds[1+a*2]+1
	 << "\n";
  for ( int a=0; a<nangles; ++a )
    cout << "A " 
	 << fmt << angles[0+a*3]+1
	 << fmt << angles[1+a*3]+1
	 << fmt << angles[2+a*3]+1
	 << "\n";
  for ( int a=0; a<ntorsions; ++a )
    cout << "D " 
	 << fmt << torsions[0+a*4]+1
	 << fmt << torsions[1+a*4]+1
	 << fmt << torsions[2+a*4]+1
	 << fmt << torsions[3+a*4]+1
	 << "\n";
  for ( int a=0; a<nhbonds; ++a )
    cout << "H " 
	 << fmt << hbonds[0+a*3]+1
	 << fmt << hbonds[1+a*3]+1
	 << fmt << hbonds[2+a*3]+1
	 << "\n";
#undef fmt
}
