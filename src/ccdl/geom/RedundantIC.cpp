#include "RedundantIC.hpp"
#include "AutoBondFrag.hpp"
#include "EuclideanGeometry.hpp"
#include "../constants.hpp"
#include <sstream>
#include <iomanip>

void ccdl::InternalCrd::Freeze( double const x ) 
{  
  mFrozen = true; 
  mConVal = x; 
}
void ccdl::InternalCrd::Freeze( int const nat, double const * c )
{
  mFrozen = true;
  mConVal = CptValue(nat,c);
}
double ccdl::InternalCrd::CptConstraintDelta( int const nat, double const * c ) const
{
  double d = 0.;
  if ( IsFrozen() )
    d = CptDifference( CptValue(nat,c), GetConstraintValue() );
  return d;
}



std::string ccdl::BondT::Print() const
{
  std::stringstream msg;
  msg << "B " 
      << std::setw(5) << i+1 << " " 
      << std::setw(5) << j+1;
  if ( IsFrozen() )
    msg << " F " << std::fixed << std::setprecision(8) 
	<< std::setw(12) << mConVal;
  msg << "\n";
  return msg.str();
}
bool ccdl::BondT::Freeze( ccdl::InternalCrd const * q )
{
  bool ok = false;
  ccdl::BondT const * p = dynamic_cast<ccdl::BondT const *>(q);
  if ( p )
    {
      if ( i == p->i and j == p->j )
	{
	  mFrozen = true;
	  mConVal = p->GetConstraintValue();
	  ok = true;
	};
    };
  return ok;
}
double ccdl::BondT::CptValue( int const , double const * c ) const
{
  return ccdl::Bond( c+3*i, c+3*j );
}
void ccdl::BondT::CptWilson( int const nat, double const * c, double * B ) const
{
  double tmp[ 3*2 ];
  ccdl::Bond( c+3*i, c+3*j, tmp );
  std::fill( B, B + 3*nat, 0. );
  for ( int u=0; u<3; ++u )
    {
      B[ u+i*3 ] = tmp[ u+0*3 ];
      B[ u+j*3 ] = tmp[ u+1*3 ];
    };
}
void ccdl::BondT::CptHessian( int const nat, double const * c, double * B ) const
{
  double grd[ 3*2 ];
  double tmp[ (3*2)*(3*2) ];
  ccdl::Bond( c+3*i, c+3*j, grd, tmp );
  std::fill( B, B + 3*nat*3*nat, 0. );
  int ats[2] = { i,j };
  for ( int a=0; a<2; ++a )
    {
      int aa = ats[a];
      for ( int u=0; u<3; ++u )
	for ( int b=0; b<2; ++b )
	  {
	    int bb = ats[b];
	    for ( int v=0; v<3; ++v )
	      B[ (v+bb*3) + (u+aa*3)*(3*nat) ] 
		= tmp[ (v+b*3) + (u+a*3)*(3*2) ];
	  };
    };
}
double ccdl::BondT::CptDifference( double const x1, double const x0 ) const
{
  return x1-x0;
}



std::string ccdl::AngleT::Print() const
{
  std::stringstream msg;
  msg << "A   " 
      << std::setw(5) << i+1 << " " 
      << std::setw(5) << j+1 << " "
      << std::setw(5) << k+1;
  if ( IsFrozen() )
    msg << " F " << std::fixed << std::setprecision(8) 
	<< std::setw(12) << mConVal;
  msg << "\n";
  return msg.str();
}
bool ccdl::AngleT::Freeze( ccdl::InternalCrd const * q )
{
  bool ok = false;
  ccdl::AngleT const * p = dynamic_cast<ccdl::AngleT const *>(q);
  if ( p )
    {
      if ( ( i == p->i and j == p->j and k == p->k ) or
	   ( k == p->i and j == p->j and i == p->k ) )
	{
	  mFrozen = true;
	  mConVal = p->GetConstraintValue();
	  ok = true;
	};
    };
  return ok;
}
double ccdl::AngleT::CptValue( int const, double const * c ) const
{
  return ccdl::Angle( c+3*i, c+3*j, c+3*k );
}
void ccdl::AngleT::CptWilson( int const nat, double const * c, double * B ) const
{
  double tmp[ 3*3 ];
  ccdl::Angle( c+3*i, c+3*j,c+3*k, tmp );
  std::fill( B, B + 3*nat, 0. );
  for ( int u=0; u<3; ++u )
    {
      B[ u+i*3 ] = tmp[ u+0*3 ];
      B[ u+j*3 ] = tmp[ u+1*3 ];
      B[ u+k*3 ] = tmp[ u+2*3 ];
    };
}
void ccdl::AngleT::CptHessian( int const nat, double const * c, double * B ) const
{
  double grd[ 3*3 ];
  double tmp[ (3*3)*(3*3) ];
  ccdl::Angle( c+3*i, c+3*j,c+3*k, grd, tmp );
  int ats[3] = { i,j,k };
  for ( int a=0; a<3; ++a )
    {
      int aa = ats[a];
      for ( int u=0; u<3; ++u )
	for ( int b=0; b<3; ++b )
	  {
	    int bb = ats[b];
	    for ( int v=0; v<3; ++v )
	      B[ (v+bb*3) + (u+aa*3)*(3*nat) ] 
		= tmp[ (v+b*3) + (u+a*3)*(3*3) ];
	  };
    };
}
double ccdl::AngleT::CptDifference( double const x1, double const x0 ) const
{
  double d = x1-x0;
  if ( d >  ccdl::PI ) d -= ccdl::TWO_PI;
  if ( d < -ccdl::PI ) d += ccdl::TWO_PI;
  return d;
}




std::string ccdl::DihedralT::Print() const
{
  std::stringstream msg;
  msg << "D " 
      << std::setw(5) << i+1 << " " 
      << std::setw(5) << j+1 << " "
      << std::setw(5) << k+1 << " "
      << std::setw(5) << l+1;
  if ( IsFrozen() )
    msg << " F " << std::fixed << std::setprecision(8) 
	<< std::setw(12) << mConVal;
  msg << "\n";
  return msg.str();
}
double ccdl::DihedralT::CptValue( int const, double const * c ) const
{
  bool linear = false;
  return ccdl::DihedralAngle( c+3*i, c+3*j, c+3*k, c+3*l, linear );
}
bool ccdl::DihedralT::Freeze( ccdl::InternalCrd const * q )
{
  bool ok = false;
  ccdl::DihedralT const * p = dynamic_cast<ccdl::DihedralT const *>(q);
  if ( p )
    {
      if ( ( i == p->i and j == p->j and k == p->k and l == p->l ) or
	   ( l == p->i and k == p->j and j == p->k and i == p->l) )
	{
	  mFrozen = true;
	  mConVal = p->GetConstraintValue();
	  ok = true;
	};
    };
  return ok;
}
void ccdl::DihedralT::CptWilson( int const nat, double const * c, double * B ) const
{
  double tmp[ 3*4 ];  
  bool linear = false;
  ccdl::DihedralAngle( c+3*i, c+3*j,c+3*k, c+3*l, linear, tmp );
  std::fill( B, B + 3*nat, 0. );
  for ( int u=0; u<3; ++u )
    {
      B[ u+i*3 ] = tmp[ u+0*3 ];
      B[ u+j*3 ] = tmp[ u+1*3 ];
      B[ u+k*3 ] = tmp[ u+2*3 ];
      B[ u+l*3 ] = tmp[ u+3*3 ];
    };
}
void ccdl::DihedralT::CptHessian( int const nat, double const * c, double * B ) const
{
  double grd[ (3*4) ];
  double tmp[ (3*4)*(3*4) ];
  bool linear = false;
  ccdl::DihedralAngle( c+3*i, c+3*j,c+3*k, c+3*l, linear, grd, tmp );
  int ats[4] = { i,j,k,l };
  for ( int a=0; a<4; ++a )
    {
      int aa = ats[a];
      for ( int u=0; u<3; ++u )
	for ( int b=0; b<4; ++b )
	  {
	    int bb = ats[b];
	    for ( int v=0; v<3; ++v )
	      B[ (v+bb*3) + (u+aa*3)*(3*nat) ] 
		= tmp[ (v+b*3) + (u+a*3)*(3*4) ];
	  };
    };
}
double ccdl::DihedralT::CptDifference( double const x1, double const x0 ) const
{
  double d = x1-x0;
  if ( d >  ccdl::PI ) d -= ccdl::TWO_PI;
  if ( d < -ccdl::PI ) d += ccdl::TWO_PI;
  return d;
}


// struct ictype
// {
//   ictype() : i(-1),j(-1),k(-1),l(-1) {}
//   ictype( int i ) : i(i),j(-1),k(-1),l(-1) {}
//   ictype( int i, int j ) : i(i),j(j),k(-1),l(-1) 
//   { 
//     if (j>i) std::swap(i,j); 
//   }
//   ictype( int i, int j, int k ) : i(i),j(j),k(k),l(-1) 
//   { 
//     if (k>i) std::swap(i,k); 
//   }
//   ictype( int i, int j, int k, int l ) : i(i),j(j),k(k),l(l) 
//   { 
//     if (k>j) { std::swap(i,l); std::swap(j,k); }
//   }
//   int i,j,k,l;
// };
// bool operator==( ictype const & lhs, ictype const & rhs )
// {
//   return (lhs.i==rhs.i) and (lhs.j==rhs.j) and (lhs.k==rhs.k) and (lhs.l==rhs.l);
// }
// bool operator==( ictype const & lhs, ictype const & rhs )
// {
//   return (lhs.i < rhs.i) or 
//     (lhs.i==rhs.i and lhs.j<rhs.j) or
//     (lhs.i==rhs.i and lhs.j==rhs.j and lhs.k<rhs.k) or
//     (lhs.i==rhs.i and lhs.j==rhs.j and lhs.k==rhs.k and lhs.l<rhs.l );
// }



ccdl::RedundantIC::RedundantIC
( int const nat, int const * z, double const * crd, bool merge_fragments )
{
  reset( nat,z,crd,merge_fragments );
}

void ccdl::RedundantIC::reset
( int const nat, int const * z, double const * crd, bool merge_fragments )
{
  std::vector< ccdl::pIC > frozens;
  for ( std::vector< ccdl::pIC >::iterator q=qs.begin(),qend=qs.end(); q!=qend; ++q )
    if ( (*q)->IsFrozen() )
      frozens.push_back( *q );
  qs.resize( 0 );

  ccdl::AutoBondFrag abf( nat,z,crd,merge_fragments );
  for ( int a=0; a<abf.nbonds; ++a )
    {
      pIC t( new ccdl::BondT( abf.bonds[0+a*2], 
			      abf.bonds[1+a*2] ) );
      qs.push_back( t );
    };
  for ( int a=0; a<abf.nangles; ++a )
    {
      pIC t( new ccdl::AngleT( abf.angles[0+a*3], 
			       abf.angles[1+a*3], 
			       abf.angles[2+a*3] ) );
      qs.push_back( t );
    };
  for ( int a=0; a<abf.ntorsions; ++a )
    {
      pIC t( new ccdl::DihedralT( abf.torsions[0+a*4], 
				  abf.torsions[1+a*4], 
				  abf.torsions[2+a*4],
				  abf.torsions[3+a*4] ) );
      qs.push_back( t );
    };

  for ( std::vector< ccdl::pIC >::iterator 
	  f=frozens.begin(), fend=frozens.end(); 
	f!=fend; ++f )
    {
      bool ok = false;
      for ( std::vector< ccdl::pIC >::iterator 
	      q=qs.begin(), qend=qs.end(); 
	    q!=qend; ++q )
	{
	  ok = (*q)->Freeze( f->get() );
	  if ( ok ) break;
	}
      if ( ! ok )
	qs.push_back( *f );
    }
}


void ccdl::RedundantIC::FreezeBond( int i, int j, double const v )
{
  ccdl::pIC f( new ccdl::BondT( i, j ) );
  f->Freeze( v );
  bool ok = false;
  for ( std::vector< ccdl::pIC >::iterator 
	  q=qs.begin(), qend=qs.end(); 
	q!=qend; ++q )
    {
      ok = (*q)->Freeze( f.get() );
      if ( ok ) break;
    }
  if ( ! ok )
    qs.push_back( f );
}


void ccdl::RedundantIC::PrintReport( std::ostream & cout )
{
  for ( std::vector< ccdl::pIC >::iterator 
	  q=qs.begin(), qend=qs.end(); 
	q!=qend; ++q )
    cout << (*q)->Print();
}
