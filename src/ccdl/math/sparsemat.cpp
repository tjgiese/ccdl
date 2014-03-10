#include "sparsemat.hpp"
#include <map>
#include <algorithm>
#include <cmath>
#include <cassert>
#include "types.hpp"
#include "types_impl.hpp"


void ccdl::sparsemat::rows_containing_col( int const icol, std::vector<int> & rows ) const
{
  typedef std::map<int,double>::const_iterator iter;
  rows.resize(0);
  for ( int i=0; i<mFast; ++i )
    {
      iter p = mData[i].find( icol );
      if ( p != mData[i].end() )
	rows.push_back( i );
    };
}

void ccdl::sparsemat::extract( ccdl::ge & a ) const
{
  typedef std::map<int,double>::const_iterator iter;
  assert( a.nfast() == mFast );
  assert( a.nslow() == mSlow );
  a = 0.;
  for ( int i=0; i<mFast; ++i )
    for ( iter p=mData[i].begin(), e=mData[i].end(); p!=e; ++p )
      a[ i + p->first * mFast] = p->second;
}

void ccdl::sparsemat::erase_upper_triangle()
{
  typedef std::map<int,double>::iterator iter;
  for ( int i=0; i<mFast; ++i )
    for ( iter p=mData[i].begin(), e=mData[i].end(); p!=e; ++p )
      if ( p->first > i )
	{
	  mData[i].erase( p, mData[i].end() );
	  break;
	};
  csr_store();
}

void ccdl::sparsemat::erase()
{
  for ( int i=0; i<mFast; ++i )
    mData[i].erase( mData[i].begin(), mData[i].end() );
  mValues.resize(0);
  mCols.resize(0);
  mRows.resize(0);
}


void ccdl::sparsesym::cholesky_lower_triangle( ccdl::sparsemat & a ) const
{
  typedef std::map<int,double>::const_iterator citer;
  typedef std::vector<int>::const_iterator rowiter;

  a.erase();
  int const N = nfast();
  for ( int i=0; i<N; ++i )
    for ( citer p=mData[i].begin(), e=mData[i].end(); p != e; ++p )
      a[i][ p->first ] = p->second;

  std::vector<int> rows;
  rows.reserve(30);

  for ( int j=0; j<N; ++j )
    {
      a[j][j] = std::sqrt( a[j][j] );
      a.rows_containing_col( j, rows );

      // for i > j
      for ( rowiter pi=rows.begin(), ei=rows.end(); pi != ei; ++pi )
	{
	  int const i = *pi;
	  if ( i > j )
	    {
	      // for k < j
	      for ( citer pk=a[j].begin(), ek=a[j].end(); pk != ek; ++pk )
		{
		  int const k = pk->first;
		  if ( k < j )
		    a[i][j] -= a[i][k] * a[j][k];
		};
	    };
	};

      // for i > j
      for ( rowiter pi=rows.begin(), ei=rows.end(); pi != ei; ++pi )
	{
	  int const i = *pi;
	  if ( i > j )
	    {
	      a[i][j] /= a[j][j];
	      a[i][i] -= std::pow(a[i][j],2);
	    };
	};

    };

  a.erase_upper_triangle();
}


ccdl::sparsecholeskypreconditioner::sparsecholeskypreconditioner( ccdl::sparsesym const & A )
  : ccdl::sparsemat(A.nfast(),A.nslow())
{

  A.cholesky_lower_triangle( reinterpret_cast< ccdl::sparsemat & >( *this ) );

  /*

  typedef std::map<int,double>::const_iterator citer;
  typedef std::vector<int>::const_iterator rowiter;
  
  int const N = A.nfast();
  for ( int i=0; i<N; ++i )
    for ( citer p=A[i].begin(), e=A[i].end(); p != e; ++p )
      mData[i][ p->first ] = p->second;
  
  std::vector<int> rows;
  rows.reserve(30);
  
  for ( int j=0; j<N; ++j )
    {
      mData[j][j] = std::sqrt( mData[j][j] );
      rows_containing_col( j, rows );

      // for i > j
      for ( rowiter pi=rows.begin(), ei=rows.end(); pi != ei; ++pi )
	{
	  int const i = *pi;
	  if ( i > j )
	    {
	      // for k < j
	      for ( citer pk=mData[j].begin(), ek=mData[j].end(); pk != ek; ++pk )
		{
		  int const k = pk->first;
		  if ( k < j )
		    mData[i][j] -= mData[i][k] * mData[j][k];
		};
	    };
	};

      // for i > j
      for ( rowiter pi=rows.begin(), ei=rows.end(); pi != ei; ++pi )
	{
	  int const i = *pi;
	  if ( i > j )
	    {
	      mData[i][j] /= mData[j][j];
	      mData[i][i] -= std::pow(mData[i][j],2);
	    };
	};

    };

  erase_upper_triangle();
  */
}
