#ifndef _ccdl_minifcn_hpp_
#define _ccdl_minifcn_hpp_

// see minimize.cpp for examples / test functions

#include <vector>

namespace ccdl
{
  struct minifcn
  {
  public:
    minifcn( int nvar );
    minifcn( int nvar, double const * refx );
    virtual ~minifcn();

    virtual double operator() ( double const alp, double const * s ) =0;
    virtual double operator() ( double const alp, double const * s, double * dydx ) =0;
    double numgrd( double const alp, double const * s, double * dydx );
    void testgrd();

    virtual void SetNewReferencePt( double const alp, double const * s );

   
    bool IsNewReferencePt() const;

    int size() const { return nvar; }

    double       * data()       { return refx.data(); }
    double const * data() const { return refx.data(); }

    double & operator[] ( int i )       { return refx[i]; }
    double   operator[] ( int i ) const { return refx[i]; }

    int nvar;
    bool is_new_reference;
    std::vector<double> refx;
  };
}

#endif
