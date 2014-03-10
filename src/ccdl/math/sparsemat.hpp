#ifndef _ccdl_sparsemat_hpp_
#define _ccdl_sparsemat_hpp_

#include <vector>
#include <map>
#include <iostream>
#include <iomanip>

namespace ccdl
{
  class ge;
}

namespace ccdl
{

  class sparsemat
  {
  public:
    sparsemat( int const m, int const n );
    
    int nfast() const;
    int nslow() const;
    
    std::map<int,double> const & operator[] ( int const i ) const;
    std::map<int,double> & operator[] ( int const i );
    
    double & operator() ( int const i, int const j );
    double operator() ( int const i, int const j ) const;
    
    void rows_containing_col( int const icol, std::vector<int> & rows ) const; 
    void extract( ccdl::ge & a ) const;
    void erase_upper_triangle();
    void erase();

    int            csr_size() const;
    double *       csr_values();
    double const * csr_values() const;
    double &       csr_values( int const i );
    double         csr_values( int const i ) const;
    int const *    csr_rowoff() const;
    int *          csr_rowoff();
    int            csr_rowoff( int const i ) const;
    int const *    csr_colidx() const;
    int *          csr_colidx();
    int            csr_colidx( int const i ) const;
    void csr_store();
    bool treat_as_transpose() const;
    ccdl::sparsemat & t();

  protected:
    int mFast,mSlow;
    bool mTranspose;
    std::vector< std::map<int,double> > mData;
    std::vector<double> mValues;
    std::vector<int> mCols,mRows;
  };


  class sparsesym : public sparsemat
  {
  public:
    sparsesym( int const nfast, int const nslow );
    void csr_store();
    void cholesky_lower_triangle( ccdl::sparsemat & l ) const;

    ccdl::sparsesym & t();
  };


  class sparsecholeskypreconditioner : public sparsemat
  {
  public:
    sparsecholeskypreconditioner( sparsesym const & m );
  };

}

inline ccdl::sparsemat::sparsemat( int const m, int const n )
  : mFast(m), mSlow(n), mTranspose(false), mData(m), mValues(m), mRows(m+1,0)
{
}
inline ccdl::sparsesym::sparsesym( int const m, int const n )
  : ccdl::sparsemat(m,n)
{
}
inline int ccdl::sparsemat::nfast() const
{
  return mFast;
}
inline int ccdl::sparsemat::nslow() const
{
  return mSlow;
}
inline std::map<int,double> const & ccdl::sparsemat::operator[] ( int const i ) const
{
  return mData[i];
}
inline std::map<int,double> & ccdl::sparsemat::operator[] ( int const i )
{
  return mData[i];
}
inline double & ccdl::sparsemat::operator() ( int const i, int const j ) 
{ 
  return mData[i][j]; 
}
inline double ccdl::sparsemat::operator() ( int const i, int const j ) const 
{ 
  return mData[i].at(j);
}
inline int ccdl::sparsemat::csr_size() const
{
  return mValues.size();
}
inline double * ccdl::sparsemat::csr_values()
{
  return mValues.data();
}
inline double const * ccdl::sparsemat::csr_values() const
{
  return mValues.data();
}
inline double & ccdl::sparsemat::csr_values( int const i )
{
  return *(mValues.data()+i);
}
inline double ccdl::sparsemat::csr_values( int const i ) const
{
  return *(mValues.data()+i);
}
inline int const * ccdl::sparsemat::csr_rowoff() const
{
  return mRows.data();
}
inline int * ccdl::sparsemat::csr_rowoff()
{
  return mRows.data();
}
inline int ccdl::sparsemat::csr_rowoff( int const i ) const
{
  return *(mRows.data()+i);
}
inline int const * ccdl::sparsemat::csr_colidx() const
{
  return mCols.data();
}
inline int * ccdl::sparsemat::csr_colidx()
{
  return mCols.data();
}
inline int ccdl::sparsemat::csr_colidx( int const i ) const
{
  return *(mCols.data()+i);
}
inline void ccdl::sparsemat::csr_store()
{
  typedef std::map<int,double>::iterator iter;
  mValues.resize(0);
  mCols.resize(0);
  mRows[0] = 0;
  for ( int i=0; i<mFast; ++i )
    {
      mRows[i+1] = mRows[i] + mData[i].size();
      for ( iter p=mData[i].begin(), e=mData[i].end(); p!=e; ++p )
	{
	  mCols.push_back( p->first );
	  mValues.push_back( p->second );
	};
    };
}
inline void ccdl::sparsesym::csr_store()
{
  typedef std::map<int,double>::iterator iter;
  mValues.resize(0);
  mCols.resize(0);
  mRows[0] = 0;
  for ( int i=0; i<mFast; ++i )
    {
      //mRows[i+1] = mRows[i] + mData[i].size();
      int ncol = 0;
      for ( iter p=mData[i].begin(), e=mData[i].end(); p!=e; ++p )
	{
	  int icol = p->first;
	  if ( icol <= i )
	    {
	      mCols.push_back( icol );
	      mValues.push_back( p->second );
	      ncol++;
	    };
	};
      mRows[i+1] = mRows[i] + ncol;
    };
}
inline bool ccdl::sparsemat::treat_as_transpose() const
{
  return mTranspose;
}
inline ccdl::sparsemat & ccdl::sparsemat::t()
{
  mTranspose = ! mTranspose;
  return *this;
}
inline ccdl::sparsesym & ccdl::sparsesym::t()
{
  mTranspose = ! mTranspose;
  return *this;
}
inline std::ostream & operator<< ( std::ostream & cout, ccdl::sparsemat const & a )
{
  for ( int i=0; i<a.nfast(); ++i )
    {
      int const b = a.csr_rowoff(i);
      int const e = a.csr_rowoff(i+1);
      for ( int jj=b; jj<e; ++jj )
	{
	  int const j = a.csr_colidx(jj);
	  double v = a.csr_values(jj);
	  cout << std::setw(7) << i 
	       << std::setw(7) << j 
	       << std::setw(20) << std::setprecision(10) 
	       << std::scientific << v << "\n";
	};
    }
  return cout;
}


#endif
