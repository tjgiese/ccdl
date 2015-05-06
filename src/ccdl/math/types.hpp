#ifndef _ccdl_types_hpp_
#define _ccdl_types_hpp_

#include <vector>

namespace ccdl
{
  class v;
  class cv;

  class ge;
  class cge;

  class sy;
  class csy;

  class di;
  class cdi;

  class cgt;

}

#include "sparsemat.hpp"



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


namespace ccdl
{
  namespace base
  {
    class vector
    {
    public:
      vector( int const n, double * d ) : mSize(n), mData(d) {};
      vector( std::vector<double> & d ) : mSize(d.size()), mData(d.data()) {};
      int size() const { return mSize; }
      double const * data() const { return mData; }
      double const * begin() const { return mData; }
      double const * end() const { return mData+mSize; }
      double operator[] ( int const i ) const { return mData[i]; }
  
      double * data() { return mData; }
      double * begin() { return mData; }
      double * end() { return mData+mSize; }
      double & operator[] ( int const i ) { return mData[i]; }
        
      void operator+= ( double const a )
      {
	for ( int i=0; i<mSize; ++i )
	  mData[i] += a;
      }

      void operator-= ( double const a )
      {
	for ( int i=0; i<mSize; ++i )
	  mData[i] -= a;
      }

      void operator*= ( double const a )
      {
	for ( int i=0; i<mSize; ++i )
	  mData[i] *= a;
      }

    protected:
      int mSize;
      double * mData;
    };



    class const_vector
    {
    public:
      const_vector( int const n, double const * d ) : mSize(n), mData(d) {};
      const_vector( std::vector<double> const & d ) : mSize(d.size()), mData(d.data()) {};
      int size() const { return mSize; }
      double const * data() const { return mData; }
      double const * begin() const { return mData; }
      double const * end() const { return mData+mSize; }
      double operator[] ( int const i ) const { return mData[i]; }
  
    protected:
      int mSize;
      double const * mData;
    };



    class matrix
    {
    public:
      matrix( int const m, int const n, double * d ) : mFast(m), mSlow(n), mSize(m*n), mData(d) {};

      int size() const { return mSize; }
      int nfast() const { return mFast; }
      int nslow() const { return mSlow; }
      double const * data() const { return mData; }
      double const * begin() const { return mData; }
      double const * end() const { return mData+mSize; }
      double operator[] ( int const i ) const { return mData[i]; }
      double * data() { return mData; }
      double * begin() { return mData; }
      double * end() { return mData+mSize; }
      double & operator[] ( int const i ) { return mData[i]; }

      void operator+= ( double const a )
      {
	for ( int i=0; i<mSize; ++i )
	  mData[i] += a;
      }

      void operator-= ( double const a )
      {
	for ( int i=0; i<mSize; ++i )
	  mData[i] -= a;
      }

      void operator*= ( double const a )
      {
	for ( int i=0; i<mSize; ++i )
	  mData[i] *= a;
      }

    protected:
      int mFast,mSlow,mSize;
      double * mData;
    };



    class const_matrix
    {
    public:
      const_matrix( int const m, int const n, double const * d ) : mFast(m), mSlow(n), mSize(m*n), mData(d) {};

      int size() const { return mSize; }
      int nfast() const { return mFast; }
      int nslow() const { return mSlow; }
      double const * data() const { return mData; }
      double const * begin() const { return mData; }
      double const * end() const { return mData+mSize; }
      double operator[] ( int const i ) const { return mData[i]; }
    protected:
      int mFast,mSlow,mSize;
      double const * mData;
    };

  }
}




//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////




namespace ccdl
{

  class v : public ccdl::base::vector
  {
  public:
    
    v( int const n, double * d );
    v( std::vector<double> & a );

    using ccdl::base::vector::operator+=;
    using ccdl::base::vector::operator-=;
    using ccdl::base::vector::operator*=;


    ccdl::v & operator+= ( ccdl::cv const a );
    ccdl::v & operator-= ( ccdl::cv const a );
    ccdl::v & operator= ( ccdl::cv const a );
    ccdl::v & operator= ( ccdl::v const & a );
    ccdl::v & operator= ( double const a );

    ccdl::v & operator*= ( ccdl::cv const a );


    double dot( ccdl::cv const x ) const;
    double nrm2() const;    
    int iabsmax() const;
    ccdl::v & axpy( double const alpha, ccdl::cv const x );
    ccdl::v & axpy( double const alpha, ccdl::cv const x, ccdl::cv const y );

    //////////////////////////////////////////////////////////////
    // overload these for custom matrix-like object
    
    template <class MAT>
    ccdl::v & dot( double const alpha, MAT A, ccdl::cv const x, double const beta );

    template <class MAT>
    ccdl::v & dot( MAT A, ccdl::cv const x );
    
    //////////////////////////////////////////////////////////////

    
    ccdl::v & dot( double const alpha, ccdl::cge const A, ccdl::cv const x, double const beta );
    ccdl::v & dot( double const alpha, ccdl::ge  const A, ccdl::cv const x, double const beta );
    ccdl::v & dot( double const alpha, ccdl::cge const A, ccdl::v  const x, double const beta );
    ccdl::v & dot( double const alpha, ccdl::ge  const A, ccdl::v  const x, double const beta );

    ccdl::v & dot( double const alpha, ccdl::cgt const A, ccdl::cv const x, double const beta );
    ccdl::v & dot( double const alpha, ccdl::cgt const A, ccdl::v  const x, double const beta );

    ccdl::v & dot( double const alpha, ccdl::csy const A, ccdl::cv const x, double const beta );
    ccdl::v & dot( double const alpha, ccdl::sy  const A, ccdl::cv const x, double const beta );
    ccdl::v & dot( double const alpha, ccdl::csy const A, ccdl::v  const x, double const beta );
    ccdl::v & dot( double const alpha, ccdl::sy  const A, ccdl::v  const x, double const beta );

    ccdl::v & dot( ccdl::cge const A, ccdl::cv const x );
    ccdl::v & dot( ccdl::ge  const A, ccdl::cv const x );
    ccdl::v & dot( ccdl::cge const A, ccdl::v  const x );
    ccdl::v & dot( ccdl::ge  const A, ccdl::v  const x );

    ccdl::v & dot( ccdl::cgt const A, ccdl::cv const x );
    ccdl::v & dot( ccdl::cgt const A, ccdl::v  const x );

    ccdl::v & dot( ccdl::csy const A, ccdl::cv const x );
    ccdl::v & dot( ccdl::sy  const A, ccdl::cv const x );
    ccdl::v & dot( ccdl::csy const A, ccdl::v  const x );
    ccdl::v & dot( ccdl::sy  const A, ccdl::v  const x );


    ccdl::v & dot( double const alpha, ccdl::sparsemat const & A, ccdl::cv const x, double const beta );
    ccdl::v & dot( double const alpha, ccdl::sparsesym const & A, ccdl::cv const x, double const beta );
    ccdl::v & dot( ccdl::sparsemat const & A, ccdl::cv const x );
    ccdl::v & dot( ccdl::sparsesym const & A, ccdl::cv const x );

    //
    // solves A.x = b
    // min { xt.b - 0.5 xt.A.x }
    //
    static int query_solve( int const N );
    ccdl::v & solve( ccdl::csy const A, ccdl::cv const b, int const nscr, double * scr  );
    ccdl::v & solve( ccdl::csy const A, ccdl::cv const b );


    //
    // solves L.Lt.x = b for x
    // where L.Lt is a sparse representation of a dense matrix A.
    // L is known as the incomplete cholesky decomposition,
    // and is a lower triangular sparse matrix. 
    //
    ccdl::v & solve( ccdl::sparsecholeskypreconditioner const & L, ccdl::cv const b );

    //
    // solves A.x = b ; st dt.x = N
    // min { xt.b - 0.5 xt.A.x + mu ( xt.d - N ) }
    //
    ccdl::v & constrained_solve( ccdl::csy const A, ccdl::cv const b, ccdl::cv const d, double const N );

    static int query_svd_solve( int const M, int const N );
    ccdl::v & svd_solve( ccdl::cge const A, ccdl::cv const b, double const tol, int const nscr, double * scr  );
    ccdl::v & svd_solve( ccdl::cge const A, ccdl::cv const b, double const tol );
    ccdl::v & svd_solve( ccdl::csy const A, ccdl::cv const b, double const tol, int const nscr, double * scr  );
    ccdl::v & svd_solve( ccdl::csy const A, ccdl::cv const b, double const tol );


    /*
    //ccdl::v & apply_orthogonal_constraints( ccdl::ge & D, ccdl::v & v );
    //
    // On input, the vector does not necessarily satisfy
    // the constraints
    // At.x = constraint_values
    // On output, it does.
    // x += constraint_projector_of_A . ( constraint_values - At.x )
    // where constraint_projector_of_A is
    // ( constraint_projector_of_A = A ).constraint_projector()
    //
    int query_enforce_constraints( int const nconstraints );
    v & enforce_constraints( ccdl::ge & A, ccdl::ge const & constraint_projector_of_A, ccdl::v const & constraint_values, int const nscr, double * scr );
    v & enforce_constraints( ccdl::ge & A, ccdl::ge const & constraint_projector_of_A, ccdl::v const & constraint_values );
    v & enforce_constraints( ccdl::ge & A, ccdl::v const & constraint_values );
    */
  };

}





//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////





namespace ccdl
{

  class cv : public ccdl::base::const_vector
  {
  public:
    cv( int const n, double const * d );
    cv( std::vector<double> const & d ); 
    cv( ccdl::v const & a ); 
    double dot( ccdl::cv const x ) const;
    double nrm2() const;
    int iabsmax() const;
  };

}






//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////






namespace ccdl
{

  class ge : public ccdl::base::matrix
  {
  public:
    
    ge( int m, int n, double * d );
    ge( int m, int n, std::vector<double> & a );


    using ccdl::base::matrix::operator+=;
    using ccdl::base::matrix::operator-=;
    using ccdl::base::matrix::operator*=;

    void operator+= ( ccdl::cge const a );
    void operator-= ( ccdl::cge const a );
    ccdl::ge & operator= ( ccdl::cge const a );
    ccdl::ge & operator= ( ccdl::ge const & a );
    ccdl::ge & operator= ( double const a );
    
    ccdl::cgt t() const;
    ccdl::v col( int const i );
    ccdl::cv col( int const i ) const;


    /////////////////////////////////////////////////////////////////
    // overload these for custom matrix-like object
    template <class MA,class MB>
    ccdl::ge & dot( double const alpha, MA & A, MB & B, double const beta );
    
    template <class MA,class MB>
    ccdl::ge & dot( MA & A, MB & B );
    ////////////////////////////////////////////////////////////////


    ccdl::ge & dot( double const alpha, ccdl::cge const A, ccdl::cge const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::ge  const A, ccdl::cge const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::cge const A, ccdl::ge  const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::ge  const A, ccdl::ge  const B, double const beta );

    ccdl::ge & dot( double const alpha, ccdl::cgt const A, ccdl::cge const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::cgt const A, ccdl::ge  const B, double const beta );

    ccdl::ge & dot( double const alpha, ccdl::cge const A, ccdl::cgt const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::ge  const A, ccdl::cgt const B, double const beta );

    ccdl::ge & dot( double const alpha, ccdl::cgt const A, ccdl::cgt const B, double const beta );

    ccdl::ge & dot( double const alpha, ccdl::csy const A, ccdl::cge const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::sy  const A, ccdl::cge const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::csy const A, ccdl::ge  const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::sy  const A, ccdl::ge  const B, double const beta );

    ccdl::ge & dot( double const alpha, ccdl::cge const A, ccdl::csy const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::ge  const A, ccdl::csy const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::cge const A, ccdl::sy  const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::ge  const A, ccdl::sy  const B, double const beta );

    ccdl::ge & dot( double const alpha, ccdl::sy  const A, ccdl::sy  const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::csy const A, ccdl::sy  const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::sy  const A, ccdl::csy const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::csy const A, ccdl::csy const B, double const beta );


    ccdl::ge & dot( ccdl::cdi const A, ccdl::cge const B );
    ccdl::ge & dot( ccdl::di  const A, ccdl::cge const B );
    ccdl::ge & dot( ccdl::cdi const A, ccdl::ge  const B );
    ccdl::ge & dot( ccdl::di  const A, ccdl::ge  const B );

    ccdl::ge & dot( ccdl::cgt const A, ccdl::cdi const B );
    ccdl::ge & dot( ccdl::cgt const A, ccdl::di  const B );


    ccdl::ge & dot( ccdl::cge const A, ccdl::cdi const B );
    ccdl::ge & dot( ccdl::cge const A, ccdl::di  const B );
    ccdl::ge & dot( ccdl::ge const A, ccdl::cdi const B );
    ccdl::ge & dot( ccdl::ge const A, ccdl::di  const B );


    ccdl::ge & dot( ccdl::cge const A, ccdl::cge const B );
    ccdl::ge & dot( ccdl::ge  const A, ccdl::cge const B );
    ccdl::ge & dot( ccdl::cge const A, ccdl::ge  const B );
    ccdl::ge & dot( ccdl::ge  const A, ccdl::ge  const B );

    ccdl::ge & dot( ccdl::cgt const A, ccdl::cge const B );
    ccdl::ge & dot( ccdl::cgt const A, ccdl::ge  const B );

    ccdl::ge & dot( ccdl::cge const A, ccdl::cgt const B );
    ccdl::ge & dot( ccdl::ge  const A, ccdl::cgt const B );

    ccdl::ge & dot( ccdl::cgt const A, ccdl::cgt const B );

    ccdl::ge & dot( ccdl::csy const A, ccdl::cge const B );
    ccdl::ge & dot( ccdl::sy  const A, ccdl::cge const B );
    ccdl::ge & dot( ccdl::csy const A, ccdl::ge  const B );
    ccdl::ge & dot( ccdl::sy  const A, ccdl::ge  const B );

    ccdl::ge & dot( ccdl::cge const A, ccdl::csy const B );
    ccdl::ge & dot( ccdl::ge  const A, ccdl::csy const B );
    ccdl::ge & dot( ccdl::cge const A, ccdl::sy  const B );
    ccdl::ge & dot( ccdl::ge  const A, ccdl::sy  const B );

    ccdl::ge & dot( ccdl::sy  const A, ccdl::sy  const B );
    ccdl::ge & dot( ccdl::csy const A, ccdl::sy  const B );
    ccdl::ge & dot( ccdl::sy  const A, ccdl::csy const B );
    ccdl::ge & dot( ccdl::csy const A, ccdl::csy const B );


    ccdl::ge & dot( double const alpha, ccdl::sparsemat const & A, ccdl::cge const B, double const beta );
    ccdl::ge & dot( double const alpha, ccdl::sparsesym const & A, ccdl::cge const B, double const beta );
    ccdl::ge & dot( ccdl::sparsemat const & A, ccdl::cge const B );
    ccdl::ge & dot( ccdl::sparsesym const & A, ccdl::cge const B );


    ccdl::ge & graham_schmidt();
    //ccdl::ge & orthogonalize_constraints( ccdl::v & constraint_values );

    static int query_svd_inverse( int const M, int const N );
    ccdl::ge & svd_inverse( double const tol, int const nscr, double * scr );
    ccdl::ge & svd_inverse( double const tol );



    static int query_svd( int const M, int const N );
    void svd( ccdl::ge & U, ccdl::di & w, ccdl::ge & VT, int const nscr, double * scr ) const;
    void svd( ccdl::ge & U, ccdl::di & w, ccdl::ge & VT ) const;


    /*
    //
    // On entry, the matrix A has dimensions M x Nconstraints
    // On exit, the matrix is
    //  P = A.(At.A)^{-1}
    // which is still M x Nconstraints
    //
    // The resulting matrix projects out constraints via
    // x += P.(c - At.x)
    // where c is the vector of Nconstraint constraint values
    //
    static int query_constraint_projector( int const M, int const Nconstraints );
    ccdl::ge & constraint_projector( int const nscr, double * scr );
    ccdl::ge & constraint_projector();
    */

  };
  

}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



namespace ccdl
{
  class cge : public ccdl::base::const_matrix
  {
  public:    
    cge( int m, int n, double const * d );
    cge( ccdl::ge const & a );
    ccdl::cgt t() const;
  };
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



namespace ccdl
{
  class cgt : public ccdl::base::const_matrix
  {
  public:    
    cgt( int m, int n, double const * d );
  };
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////




namespace ccdl
{
  class sy : public ccdl::base::matrix
  {
  public:
    sy( int const n, double * d );
    sy( int const n, std::vector<double> & d );

    using ccdl::base::matrix::operator+=;
    using ccdl::base::matrix::operator-=;
    using ccdl::base::matrix::operator*=;

    void operator+= ( ccdl::csy const a );
    void operator-= ( ccdl::csy const a );
    ccdl::sy & operator= ( ccdl::csy const a );
    ccdl::sy & operator= ( ccdl::sy const & a );
    ccdl::sy & operator= ( double const a );

    ccdl::cge ge() const;
    ccdl::ge ge();

    ccdl::sy & dot( double const alpha, ccdl::cge const A, ccdl::cge const B, double const beta );
    ccdl::sy & dot( double const alpha, ccdl::cgt const A, ccdl::cge const B, double const beta );
    ccdl::sy & dot( ccdl::cge const A, ccdl::cge const B );
    ccdl::sy & dot( ccdl::cgt const A, ccdl::cge const B );

    ccdl::sy & inverse();

    static int query_svd_inverse( int const N );
    ccdl::sy & svd_inverse( double const tol, int const nscr, double * scr );
    ccdl::sy & svd_inverse( double const tol );

    // static int query_safe_inverse( int const N );
    // ccdl::sy & safe_inverse( double const tol, int const nscr, double * scr );
    // ccdl::sy & safe_inverse( double const tol );

    static int query_svd( int const N );
    void svd( ccdl::ge & U, ccdl::di & w, ccdl::ge & VT, int const nscr, double * scr ) const;
    void svd( ccdl::ge & U, ccdl::di & w, ccdl::ge & VT ) const;


    static int query_dsyev( int const N );
    static int query_dsyevd( int const N );
    static int query_dsyevr( int const N );
    static int iquery_dsyevd( int const N );
    static int iquery_dsyevr( int const N );

    void dsyev( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const;
    void dsyevd( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const;
    void dsyevr( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const;
    void dsyev( ccdl::di & E, ccdl::ge & U ) const;
    void dsyevd( ccdl::di & E, ccdl::ge & U ) const;
    void dsyevr( ccdl::di & E, ccdl::ge & U ) const;

    static int query_eigen( int const N );
    void eigen( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const;
    void eigen( ccdl::di & E, ccdl::ge & U ) const;

  };
}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



namespace ccdl
{
  class csy : public ccdl::base::const_matrix
  {
  public:    
    csy( int m, int n, double const * d );
    csy( int m, int n, std::vector<double> const & d );
    csy( ccdl::sy const & a );
    ccdl::cge ge() const;



    void svd( ccdl::ge & U, ccdl::di & w, ccdl::ge & VT, int const nscr, double * scr ) const;
    void svd( ccdl::ge & U, ccdl::di & w, ccdl::ge & VT ) const;

    void dsyev( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const;
    void dsyevd( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const;
    void dsyevr( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const;
    void dsyev( ccdl::di & E, ccdl::ge & U ) const;
    void dsyevd( ccdl::di & E, ccdl::ge & U ) const;
    void dsyevr( ccdl::di & E, ccdl::ge & U ) const;

    void eigen( ccdl::di & E, ccdl::ge & U, int const nscr, double * scr ) const;
    void eigen( ccdl::di & E, ccdl::ge & U ) const;
  };
}




//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



namespace ccdl
{

  class di : public ccdl::base::vector
  {
  public:

    di( int const nfast, int const nslow, double * data );
    di( int const nfast, int const nslow, std::vector<double> & data );

    int nfast() const;
    int nslow() const;


    using ccdl::base::vector::operator+=;
    using ccdl::base::vector::operator-=;
    using ccdl::base::vector::operator*=;

    void operator+= ( ccdl::cdi const a );
    void operator-= ( ccdl::cdi const a );
    ccdl::di & operator= ( ccdl::cdi const a );
    ccdl::di & operator= ( ccdl::di const & a );
    ccdl::di & operator= ( double const a );

    //ccdl::di  t();
    ccdl::cdi t() const;

    ccdl::di & inverse( double const tol );

  private:
    int mFast,mSlow;
  };

}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


namespace ccdl
{

  class cdi : public ccdl::base::const_vector
  {
  public:

    cdi( int const nfast, int const nslow, double const * data );
    cdi( int const nfast, int const nslow, std::vector<double> const & data );
    cdi( ccdl::di const & a );

    int nfast() const;
    int nslow() const;
    cdi t() const;

  private:
    int mFast,mSlow;
  };

}


#endif
