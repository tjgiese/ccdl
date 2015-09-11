#include "RedundantIC.hpp"
#include "AutoBondFrag.hpp"
#include "EuclideanGeometry.hpp"
#include "../constants.hpp"
#include "../bmath.hpp"
#include <sstream>
#include <iomanip>
#include <cstdio>

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
      << std::setw(5) << j+1 << " "
      << std::setw(5) << 0 << " "
      << std::setw(5) << 0;
  if ( IsFrozen() )
    msg << " F " << std::fixed << std::setprecision(8) 
	<< std::setw(13) << mConVal / ccdl::AU_PER_ANGSTROM << " A  ";
  else
    msg << " A " << std::setw(13) << " " << "    ";
  return msg.str();
}
std::string ccdl::BondT::PrintPrettyValue( double a ) const
{
  a /= ccdl::AU_PER_ANGSTROM;
  std::stringstream msg;
  msg << std::fixed << std::setprecision(8) 
      << std::setw(13) << a << " A  ";
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
double ccdl::BondT::MoveValueToValidRange( double d ) const
{
  return d;
}



std::string ccdl::AngleT::Print() const
{
  std::stringstream msg;
  msg << "A " 
      << std::setw(5) << i+1 << " " 
      << std::setw(5) << j+1 << " "
      << std::setw(5) << k+1 << " "
      << std::setw(5) << 0;
  if ( IsFrozen() )
    msg << " F " << std::fixed << std::setprecision(8) 
	<< std::setw(13) << mConVal * 180./ccdl::PI << " deg";
  else
    msg << " A " << std::setw(13) << " " << "    ";
  return msg.str();
}
std::string ccdl::AngleT::PrintPrettyValue( double a ) const
{
  a *= 180./ccdl::PI;
  std::stringstream msg;
  msg << std::fixed << std::setprecision(8) 
      << std::setw(13) << a << " deg";
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
  double dold = d;
  while ( true )
    {
      if ( d >  ccdl::PI ) d = ccdl::TWO_PI-d;
      if ( d < -ccdl::PI ) d = ccdl::TWO_PI+d;
      if ( d == dold ) break;
      else dold = d;
    };
  return d;
}
double ccdl::AngleT::MoveValueToValidRange( double d ) const
{
  double dold = d;
  while ( true )
    {
      if ( d > ccdl::PI ) d = ccdl::TWO_PI-d;
      if ( d < 0. ) d = ccdl::TWO_PI+d;
      if ( d == dold ) break;
      else dold = d;
    };
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
	<< std::setw(13) << mConVal * 180./ccdl::PI << " deg";
  else
    msg << " A " << std::setw(13) << " " << "    ";
  return msg.str();
}
std::string ccdl::DihedralT::PrintPrettyValue( double a ) const
{
  a *= 180./ccdl::PI;
  std::stringstream msg;
  msg << std::fixed << std::setprecision(8) 
      << std::setw(13) << a << " deg";
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
	   ( l == p->i and k == p->j and j == p->k and i == p->l ) )
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
  //std::printf("angle %13.4e %i\n",angle,(int)linear); 
  if ( ! linear )
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
  double dold = d;
  while ( true )
    {
      //if ( d >  ccdl::PI ) d = ccdl::TWO_PI - d;
      //if ( d < -ccdl::PI ) d = ccdl::TWO_PI + d;
      if ( d >  ccdl::PI ) d -= ccdl::TWO_PI;
      if ( d < -ccdl::PI ) d += ccdl::TWO_PI;
      if ( d == dold ) break;
      else dold = d;
    };
  return d;
}
double ccdl::DihedralT::MoveValueToValidRange( double d ) const
{
  double dold = d;
  while ( true )
    {
      if ( d >  ccdl::PI ) d -= ccdl::TWO_PI;
      if ( d < -ccdl::PI ) d += ccdl::TWO_PI;
      if ( d == dold ) break;
      else dold = d;
    };
  return d;
}












std::string ccdl::CartT::Print() const
{
  std::stringstream msg;
  if ( k == 0 )
    msg << "X ";
  else if ( k == 1 )
    msg << "Y ";
  else
    msg << "Z ";
  msg << std::setw(5) << a+1 << " " 
      << std::setw(5) << 0 << " "
      << std::setw(5) << 0 << " "
      << std::setw(5) << 0;
  if ( IsFrozen() )
    msg << " F " << std::fixed << std::setprecision(8) 
	<< std::setw(13) << mConVal / ccdl::AU_PER_ANGSTROM << " A  ";
  else
    msg << " A " << std::setw(13) << " " << "    ";
  return msg.str();
}
std::string ccdl::CartT::PrintPrettyValue( double a ) const
{
  a /= ccdl::AU_PER_ANGSTROM;
  std::stringstream msg;
  msg << std::fixed << std::setprecision(8) 
      << std::setw(13) << a << " A  ";
  return msg.str();
}

bool ccdl::CartT::Freeze( ccdl::InternalCrd const * q )
{
  bool ok = false;
  ccdl::CartT const * p = dynamic_cast<ccdl::CartT const *>(q);
  if ( p )
    {
      if ( a == p->a and k == p->k )
	{
	  mFrozen = true;
	  mConVal = p->GetConstraintValue();
	  ok = true;
	};
    };
  return ok;
}
double ccdl::CartT::CptValue( int const , double const * c ) const
{
  return c[k+a*3];
}
void ccdl::CartT::CptWilson( int const nat, double const * , double * B ) const
{
  std::fill( B, B + 3*nat, 0. );
  B[ k + a*3 ] = 1.;
}
void ccdl::CartT::CptHessian( int const nat, double const * , double * B ) const
{
  std::fill( B, B + 3*nat*3*nat, 0. );
}
double ccdl::CartT::CptDifference( double const x1, double const x0 ) const
{
  return x1-x0;
}
double ccdl::CartT::MoveValueToValidRange( double d ) const
{
  return d;
}











std::string ccdl::R12T::Print() const
{
  std::stringstream msg;
  msg << "B " 
      << std::setw(5) << i+1 << " " 
      << std::setw(5) << j+1 << " "
      << std::setw(5) << k+1 << " "
      << std::setw(5) << l+1;
  if ( IsFrozen() )
    msg << " F " << std::fixed << std::setprecision(8) 
	<< std::setw(13) << mConVal / ccdl::AU_PER_ANGSTROM << " A  ";
  else
    msg << " A " << std::setw(13) << " " << "    ";
  return msg.str();
}
std::string ccdl::R12T::PrintPrettyValue( double a ) const
{
  a /= ccdl::AU_PER_ANGSTROM;
  std::stringstream msg;
  msg << std::fixed << std::setprecision(8) 
      << std::setw(13) << a << " A  ";
  return msg.str();
}

bool ccdl::R12T::Freeze( ccdl::InternalCrd const * q )
{
  bool ok = false;
  ccdl::R12T const * p = dynamic_cast<ccdl::R12T const *>(q);
  if ( p )
    {
      if ( i == p->i and j == p->j and k == p->k and l == p->l )
	{
	  mFrozen = true;
	  mConVal = p->GetConstraintValue();
	  ok = true;
	};
    };
  return ok;
}
double ccdl::R12T::CptValue( int const , double const * c ) const
{
  return ccdl::Bond( c+3*i, c+3*j )-ccdl::Bond( c+3*k, c+3*l );
}
void ccdl::R12T::CptWilson( int const nat, double const * c, double * B ) const
{
  double tmp1[ 3*2 ],tmp2[ 3*2 ];
  ccdl::Bond( c+3*i, c+3*j, tmp1 );
  ccdl::Bond( c+3*k, c+3*l, tmp2 );
  std::fill( B, B + 3*nat, 0. );
  for ( int u=0; u<3; ++u )
    {
      B[ u+i*3 ] += tmp1[ u+0*3 ];
      B[ u+j*3 ] += tmp1[ u+1*3 ];
      B[ u+k*3 ] -= tmp2[ u+0*3 ];
      B[ u+l*3 ] -= tmp2[ u+1*3 ];
    };
}
void ccdl::R12T::CptHessian( int const nat, double const * c, double * B ) const
{
  std::fill( B, B + 3*nat*3*nat, 0. );
  double grd[ 3*2 ];
  double tmp[ (3*2)*(3*2) ];

  ccdl::Bond( c+3*i, c+3*j, grd, tmp );
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
  ccdl::Bond( c+3*k, c+3*l, grd, tmp );
  ats[0] = k;
  ats[1] = l;
  for ( int a=0; a<2; ++a )
    {
      int aa = ats[a];
      for ( int u=0; u<3; ++u )
	for ( int b=0; b<2; ++b )
	  {
	    int bb = ats[b];
	    for ( int v=0; v<3; ++v )
	      B[ (v+bb*3) + (u+aa*3)*(3*nat) ] 
		-= tmp[ (v+b*3) + (u+a*3)*(3*2) ];
	  };
    };
}
double ccdl::R12T::CptDifference( double const x1, double const x0 ) const
{
  return x1-x0;
}
double ccdl::R12T::MoveValueToValidRange( double d ) const
{
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
( int const nat, int const * z, double const * crd, ccdl::CRDSYS crdsys, bool merge_fragments )
{
  reset( nat,z,crd,crdsys,merge_fragments );
}

void ccdl::RedundantIC::reset
( int const num_atoms, int const * z, double const * crd,ccdl::CRDSYS crdsys, bool merge_fragments )
{
  nat = num_atoms;
  std::vector< ccdl::pIC > frozens;
  for ( std::vector< ccdl::pIC >::iterator q=qs.begin(),qend=qs.end(); q!=qend; ++q )
    if ( (*q)->IsFrozen() )
      frozens.push_back( *q );
  qs.resize( 0 );


  ccdl::AutoBondFrag abf( nat,z,crd, crdsys, merge_fragments );

  for ( int a=0; a<abf.ncarts; ++a )
    for ( int k=0; k<3; ++k )
      {
	pIC t( new ccdl::CartT( abf.carts[a], k ) );
	qs.push_back( t );
      };
  

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
void ccdl::RedundantIC::FreezeAngle( int i, int j, int k, double const v )
{
  ccdl::pIC f( new ccdl::AngleT( i, j, k ) );
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
void ccdl::RedundantIC::FreezeDihedral( int i, int j, int k, int l, double const v )
{
  ccdl::pIC f( new ccdl::DihedralT( i, j, k, l ) );
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

void ccdl::RedundantIC::FreezeR12( int i, int j, int k, int l, double const v )
{
  ccdl::pIC f( new ccdl::R12T( i, j, k, l ) );
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

void ccdl::RedundantIC::FreezeR12( int i, int j, int l, double const v )
{
  FreezeR12( i,j, j,l, v );
}

void ccdl::RedundantIC::FreezeAtom( int a, double const * atom_crd )
{
   for ( int k=0; k<3; ++k )
     {
        ccdl::pIC f( new ccdl::CartT( a, k ) );
        f->Freeze( atom_crd[k] );
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
     };
}

void ccdl::RedundantIC::Freeze( std::tr1::shared_ptr< ccdl::InternalCrd > f )
{
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
    cout << (*q)->Print() << "\n";
}


void ccdl::RedundantIC::PrintReport( std::ostream & cout, double const * crd )
{
  for ( std::vector< ccdl::pIC >::iterator 
	  q=qs.begin(), qend=qs.end(); 
	q!=qend; ++q )
    {
      double a = (*q)->CptValue( nat, crd );
      cout << (*q)->Print() << " @ " << (*q)->PrintPrettyValue(a) << "\n";
    };
}




void ccdl::RedundantIC::CptTransformData
( double const * crd, double * B, double * G, double * Ginv )
{
  int const nat3 = 3*nat;
  int const nq = GetNumInternalCrds();
  std::fill( B, B + nq*nat3, 0. );
  // Bxq = dq/dx
  for ( int i=0; i<nq; ++i )
    qs[i]->CptWilson( nat, crd, B + i*nat3 );
  // G = Bt . B
  ccdl::gt_dot_ge( nat3, nq, B, nat3, nq, B, G );
  std::copy( G, G+nq*nq, Ginv );
  // G^{-1}
  ccdl::sdd_sym_power( -1., nq, nq, Ginv, 1.e-9 );
}


int ccdl::RedundantIC::GetConstraintMask( double * C )
{
  int ncon = 0;
  int const nq = GetNumInternalCrds();
  std::fill( C, C+nq, 0. );
  for ( int i=0; i<nq; ++i )
    if ( qs[i]->IsFrozen() ) 
      {
	++ncon;
	C[i] = 1.;
      }
  return ncon;
}

void ccdl::RedundantIC::CptProjector
( double const * G, double const * Ginv, double * P )
{
  int const nq = GetNumInternalCrds();
  std::vector<double> C(nq,0.);
  int const ncon = GetConstraintMask( C.data() );
  if ( ncon == 0 )
    {
      ccdl::sy_dot_ge( nq,nq, G, nq,nq, Ginv, P );
    }
  else
    {
      std::vector<double> Pp(nq*nq,0.);
      // P' = G . G^{-1}
      ccdl::sy_dot_ge( nq,nq, G, nq,nq, Ginv, Pp.data() );
      std::copy( Pp.data(), Pp.data()+nq*nq, P );
      std::vector<double> PC(nq*nq,0.);
      ccdl::ge_dot_di( nq,nq, Pp.data(), nq,nq, C.data(), PC.data() );
      std::vector<double> CPC(PC);
      // C.P'.C
      ccdl::di_dot_ge( nq,nq, C.data(), nq,nq, PC.data(), CPC.data() );
      // (C.P'.C)^{-1}
      ccdl::sdd_sym_power( -1., nq, nq, CPC.data(), 1.e-9 );
      ccdl::ge_dot_sy( nq,nq, PC.data(), nq,nq, CPC.data(), Pp.data() );
      // P = P' - (P'C).(C.P'.C)^{-1}.(C.P')
      ccdl::ge_dot_gt( nq,nq, Pp.data(), nq,nq, PC.data(), P, -1., 1. );
    }
}



void ccdl::RedundantIC::CptInternalCrds
( double const * crd, double * q )
{
  int const nq = qs.size();
  for ( int i=0; i<nq; ++i )
    q[i] = qs[i]->CptValue( nat, crd );
}




void ccdl::RedundantIC::DisplaceByDeltaQ
( double * dq, double * crd0, double const TOL )
{
  int const nat3 = nat*3;

  std::vector<double> new_crd( nat3, 0. );
  double *__restrict__ crd1 = new_crd.data();

  int const nq = GetNumInternalCrds();
  std::vector<double> old_q(nq,0.), new_q(nq,0.), tmp_q(nq,0.);
  double *__restrict__ q0 = old_q.data();
  double *__restrict__ q1 = new_q.data();
  double *__restrict__ qt = tmp_q.data();
  
  std::vector<double> Bmat(nat3*nq,0.), Gmat(nq*nq,0.), Ginvmat(nq*nq,0.);
  double *__restrict__ B = Bmat.data();
  double *__restrict__ G = Gmat.data();
  double *__restrict__ Ginv = Ginvmat.data();

  std::vector<double> C(nq,0.);
  int ncon = GetConstraintMask( C.data() );

  CptInternalCrds( crd0, q0 );

  for ( int i=0; i<nq; ++i )
    {
      double final = qs[i]->MoveValueToValidRange( dq[i]+q0[i] );
      if ( C[i] > 0. )
	final = qs[i]->MoveValueToValidRange( qs[i]->GetConstraintValue() );
      dq[i] = qs[i]->CptDifference( final, q0[i] );
    };

  std::vector<double> backup;
  double dqnorm = 0.;
  for ( int i=0; i<nq; ++i ) dqnorm += dq[i]*dq[i];
  bool OK = false;
  if ( dqnorm > 1.e-30 )
    {
      for ( int iter=0; iter < 50; ++iter )
	{
	   // std::printf("dq  ");
	   // for ( int i=0; i<nq; ++i )
	   //   std::printf("%14.6e",dq[i]);
	   // std::printf("\n");

	  CptTransformData( crd0, B, G, Ginv );
	  ccdl::sy_dot_v(nq,nq,Ginv,dq,qt);
	  ccdl::ge_dot_v(nat3,nq,B,qt,crd1);

	  if ( iter == 0 ) backup.assign( crd1, crd1+nat3 );

	  double rms = 0.;
	  for ( int i=0; i<nat3; ++i ) rms += crd1[i]*crd1[i];
	  rms = std::sqrt(rms/nat3);
	  for ( int i=0; i<nat3; ++i ) crd1[i] += crd0[i];
	  CptInternalCrds( crd1, q1 );
	  for ( int i=0; i<nq; ++i ) 
	    {
	      dq[i] -= qs[i]->CptDifference(q1[i],q0[i]);
	      double final = qs[i]->MoveValueToValidRange( dq[i]+q1[i] );
	      dq[i]  = qs[i]->CptDifference( final, q1[i] );
	    }

	  std::copy( q1, q1+nq, q0 );
	  std::copy( crd1, crd1+nat3, crd0 );
	  if ( rms < TOL ) 
	    { 
	      OK = true;
	      break;
	    };
	};
      if ( ! OK )
	{
	  std::copy( backup.data(), backup.data()+nat3, crd0 );
	  return;
	};
    };

  if ( ncon > 0 and OK )
    {
      
      for ( int i=0; i<nq; ++i )
	{
	  //if ( C[i] < 0.9 ) C[i] = 0.5;
	  double final = qs[i]->MoveValueToValidRange( qs[i]->GetConstraintValue() );
	  dq[i] = C[i] * qs[i]->CptDifference( final, q0[i] );
	};
      dqnorm = 0.;
      for ( int i=0; i<nq; ++i ) dqnorm += dq[i]*dq[i];

      //if ( dqnorm > 0.1 )
      //std::cerr << "A constraint is redundant with an active coordinate in a manner that cannot acheive the desired displacement simultaneously with the constraint.\n";
	
      if ( dqnorm > 1.e-30 )
	for ( int iter=0; iter < 50; ++iter )
	  {

	    {
	      // double final = qs[nq-1]->MoveValueToValidRange( qs[nq-1]->GetConstraintValue() );
	      // std::printf("dq %3i ",iter);
	      // std::printf("%14.6f %14.6f %14.6f",
	      // 		  q0[nq-1]*180/ccdl::PI,
	      // 		  final*180/ccdl::PI,
	      // 		  dq[nq-1]*180/ccdl::PI);
	      // std::printf("\n");
	    }


	    CptTransformData( crd0, B, G, Ginv );


	    ccdl::sy_dot_v(nq,nq,Ginv,dq,qt);
	    ccdl::ge_dot_v(nat3,nq,B,qt,crd1);


	    double rms = 0.;
	    for ( int i=0; i<nat3; ++i ) rms += crd1[i]*crd1[i];
	    rms = std::sqrt(rms/nat3);
	    for ( int i=0; i<nat3; ++i ) crd1[i] += crd0[i];
	    CptInternalCrds( crd1, q1 );

	    // std::printf("rms %14.10f\n",rms);
	    // for ( int i=0; i<nq; ++i )
	    //   std::printf("q1 %3i %20.10f\n",i,q1[i]);

	    for ( int i=0; i<nq; ++i ) 
	      {
		dq[i] -= C[i] * qs[i]->CptDifference(q1[i],q0[i]);
		double final = qs[i]->MoveValueToValidRange( dq[i]+q1[i] );
		dq[i]  = C[i] * qs[i]->CptDifference( final, q1[i] );
	      }
	    std::copy( q1, q1+nq, q0 );
	    std::copy( crd1, crd1+nat3, crd0 );
	    if ( rms < TOL ) break;
	  };
    };
}







void ccdl::RedundantIC::GrdAndHesTransform
( double const *__restrict__ crd, 
  double *__restrict__ cg, 
  double *__restrict__ ch, 
  double *__restrict__ qg, 
  double *__restrict__ qh,
  bool c2q_AND_q2c )
{
  int ncg = nat*3;
  int nch = ncg*ncg;
  int nq  = qs.size();
  int nqh = nq*nq;

  int tsize = std::max( nch, std::max( nq*nq, ncg*nq ) );
  std::vector<double> Bmat(ncg*nq,0.),BmatGinv(tsize,0.),Ginvmat(nq*nq,0.);
  std::vector<double> Hxmat(ch,ch+nch),Htmat(tsize,0.);
  double *__restrict__ B      = Bmat.data();
  double *__restrict__ BGinv  = BmatGinv.data();
  double *__restrict__ Ginv   = Ginvmat.data();
  double *__restrict__ Hx     = Hxmat.data();
  double *__restrict__ Ht     = Htmat.data();
  CptTransformData( crd, B, BGinv, Ginv );

  std::vector<double> P(nq*nq,0.);
  //int ncon = GetConstraintMask( C.data() );
  CptProjector( BGinv, Ginv, P.data() );

  // Eq (3)
  ccdl::ge_dot_sy( ncg, nq, B, nq, nq, Ginv, BGinv );
  ccdl::gt_dot_v( ncg, nq, BGinv, cg, qg );

  //for ( int i=0; i<nq; ++i ) std::printf("%14.5e",qg[i]); std::printf("\n");

  // Eq (6)
  for ( int q=0; q<nq; ++q )
    {
      qs[q]->CptHessian( nat, crd, Ht );
      ccdl::axpy( -qg[q], nch, Ht, Hx );
    };
  ccdl::sy_dot_ge( ncg, ncg, Hx, ncg, nq, BGinv, Ht );
  ccdl::gt_dot_ge( ncg, nq, BGinv, ncg, nq, Ht, qh ); 

  // Eq (9)
  ccdl::sy_dot_v( nq, nq, P.data(), qg, Ht );
  std::copy( Ht, Ht+nq, qg );

  //for ( int i=0; i<nq; ++i ) std::printf("%14.5e",qg[i]); std::printf("\n");

  ccdl::sy_dot_ge( nq, nq, qh, nq, nq, P.data(), Ht );
  ccdl::sy_dot_ge( nq, nq, P.data(), nq, nq, Ht, qh );
  for ( int i=0; i<nq; ++i )
    P[i+i*nq] -= 1.;
  for ( int i=0; i<nqh; ++i )
    qh[i] -= 1000. * P[i];

  //std::fill( qh, qh+nqh, 0. ); for ( int i=0; i<nq; ++i ) qh[i+i*nq] = 1.;

  if ( c2q_AND_q2c )
    {
      // cg = B . qg
      ccdl::ge_dot_v( ncg, nq, B, qg, cg );

      // ch = B.qh.Bt + B'.qg
      ccdl::ge_dot_sy( ncg, nq, B,  nq, nq, qh,  Ht );
      ccdl::ge_dot_gt( ncg, nq, Ht, ncg, nq, B, ch );

      for ( int q=0; q<nq; ++q )
	{
	  qs[q]->CptHessian( nat, crd, Ht );
	  ccdl::axpy( qg[q], nch, Ht, ch );
	};
    }

}


void ccdl::RedundantIC::GrdTransform
( double const *__restrict__ crd, 
  double *__restrict__ cg, 
  double *__restrict__ qg, 
  bool c2q_AND_q2c )
{
  int ncg = nat*3;
  int nch = ncg*ncg;
  int nq  = qs.size();

  int tsize = std::max( nch, std::max( nq*nq, ncg*nq ) );
  std::vector<double> Bmat(ncg*nq,0.),BmatGinv(tsize,0.),Ginvmat(nq*nq,0.);
  std::vector<double> Htmat(tsize,0.);
  double *__restrict__ B      = Bmat.data();
  double *__restrict__ BGinv  = BmatGinv.data();
  double *__restrict__ Ginv   = Ginvmat.data();
  double *__restrict__ Ht     = Htmat.data();
  CptTransformData( crd, B, BGinv, Ginv );

  std::vector<double> P(nq*nq,0.);
  CptProjector( BGinv, Ginv, P.data() );

  // Eq (3)
  ccdl::ge_dot_sy( ncg, nq, B, nq, nq, Ginv, BGinv );
  ccdl::gt_dot_v( ncg, nq, BGinv, cg, qg );

  // Eq (9)
  ccdl::sy_dot_v( nq, nq, P.data(), qg, Ht );
  std::copy( Ht, Ht+nq, qg );

  if ( c2q_AND_q2c )
    {
      // Eq (2)
      // cg = B . qg
      ccdl::ge_dot_v( ncg, nq, B, qg, cg );

      // std::printf("qg ");
      // for ( int i=0; i<nq; ++i )
      // 	std::printf("%13.4e",qg[i]);
      // std::printf("\n");

      // for ( int i=0; i<ncg; ++i )
      // 	{
      // 	  for ( int j=0; j<nq; ++j )
      // 	    std::printf("%13.4e",B[i+j*ncg]);
      // 	  std::printf("\n");
      // 	}

    }

}



void ccdl::RedundantIC::HesBackTransform
( double const *__restrict__ crd, 
  double const *__restrict__ qg, 
  double const *__restrict__ qh,
  double *__restrict__ ch )
{
  int ncg = nat*3;
  int nch = ncg*ncg;
  int nq  = qs.size();

  int hsize = std::max( nch, std::max( nq*nq, nq*ncg ) );
  std::vector<double> Bmat(ncg*nq,0.);
  std::vector<double> Htmat(hsize,0.);
  double *__restrict__ B      = Bmat.data();
  double *__restrict__ Ht     = Htmat.data();

  std::fill( B, B + nq*ncg, 0. );
  // Bxq = dq/dx
  for ( int i=0; i<nq; ++i )
    qs[i]->CptWilson( nat, crd, B + i*ncg );

  // ch = B.qh.Bt + B'.qg
  ccdl::ge_dot_sy( ncg, nq, B,  nq, nq, qh,  Ht );
  ccdl::ge_dot_gt( ncg, nq, Ht, ncg, nq, B, ch );
  for ( int q=0; q<nq; ++q )
    {
      qs[q]->CptHessian( nat, crd, Ht );
      ccdl::axpy( qg[q], nch, Ht, ch );
    };
}


void ccdl::RedundantIC::PrintCrd( std::ostream & cout, int i )
{
  cout << qs[i]->Print();
} 

void ccdl::RedundantIC::PrintPrettyValue( std::ostream & cout, int i, double const * c )
{
  double a = qs[i]->CptValue( nat, c );
  cout << qs[i]->PrintPrettyValue(a);
} 


void ccdl::RedundantIC::CptDifference( double const * hi, double const * lo, double * diff )
{
  int nq = GetNumInternalCrds();
  for ( int i=0; i<nq; ++i ) 
    {
      double a = qs[i]->MoveValueToValidRange( hi[i] );
      double b = qs[i]->MoveValueToValidRange( lo[i] );
      diff[i]  = qs[i]->CptDifference( a, b );
    }
}






void ccdl::RedundantIC::CptQ2CEstimator
( double const * x0, double * BGinv )
{
  int const nat3 = nat*3;
  int const nq = GetNumInternalCrds();

  std::vector<double> Bmat(nat3*nq,0.), Gmat(nq*nq,0.), Ginvmat(nq*nq,0.);
  double *__restrict__ B = Bmat.data();
  double *__restrict__ G = Gmat.data();
  double *__restrict__ Ginv = Ginvmat.data();

  CptTransformData( x0, B, G, Ginv );
  ccdl::ge_dot_sy( nat3, nq, B, nq, nq, Ginv, BGinv );
}


void ccdl::RedundantIC::EstimateQ2C
( double const * BGinv, double const * DQ, 
  double * approx_dx )
{
  int const nat3 = nat*3;
  int const nq = GetNumInternalCrds();

  std::vector<double> dq( DQ, DQ+nq );
  for ( int i=0; i<nq; ++i )
    if ( qs[i]->IsFrozen() ) 
      dq[i] = 0.;
  ccdl::ge_dot_v( nat3, nq, BGinv, dq.data(), approx_dx );
}
