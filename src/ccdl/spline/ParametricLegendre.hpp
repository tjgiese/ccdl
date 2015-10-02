#ifndef _ParametricLegendre_hpp_
#define _ParametricLegendre_hpp_

#include <vector>

namespace ccdl
{
  struct ParametricLegendre
  {
    ParametricLegendre( int const nbas, int const npt, int const naux );

    int GetNauxpts() const { return naux; }
    int GetNpts() const { return npts; }
    int GetNbasis() const { return nbas; }

    // this go over npt
    double GetQuadPt( int const ipt ) const { return ts[ipt]; }
    // this go over naux
    double GetAuxPt( int const iaux ) const { return xts[iaux]; }

    // these arrays are npt
    double const * GetXpts() const { return xpts.data(); }
    double const * GetYpts() const { return ypts.data(); }


    void SetEndptValues( double const x0, double const y0, double const x1, double const y1 );
    void SetEndptValues();

    // these arrays should be npt
    void SetXY( double const * x, double const * y );
    void GetXY( double const t, double & x, double & y ) const;
    void GetXY( double const t, double & x, double & y, double & dxdt, double & dydt ) const;
    void GetXY( double const t, double & x, double & y, double & dxdt, double & dydt, double & d2xt2, double & d2ydt2 ) const;

    double PerpForceProj( double const t, double & dfdx, double & dfdy );
    double ParaForceProj( double const t, double & dfdx, double & dfdy );

    // the inputs are size naux, the outputs are size npts
    double PerpForceSumSq( double const * dfdx, double const * dfdy );
    double PerpForceSumSqCoefDer( double const * dfdx, double const * dfdy, 
				  double * dSqdcx, double * dSqdcy );
    double PerpForceSumSqPtDer( double const * dfdx, double const * dfdy, 
				double * dSqdX, double * dSqdY );

    double ParaForceSumSq( double const * dfdx, double const * dfdy );
    double ParaForceSumSqCoefDer( double const * dfdx, double const * dfdy, 
				  double * dSqdcx, double * dSqdcy );
    double ParaForceSumSqPtDer( double const * dfdx, double const * dfdy, 
				double * dSqdX, double * dSqdY );

    // double TotalForceSumSq( double const * dfdx, double const * dfdy );
    // double TotalForceSumSqCoefDer( double const * dfdx, double const * dfdy, 
    // 				   double * dSqdcx, double * dSqdcy );
    // double TotalForceSumSqPtDer( double const * dfdx, double const * dfdy, 
    // 				 double * dSqdX, double * dSqdY );

    int nbas,npts,naux;
    std::vector<double> ws,ts,xws,xts;
    std::vector<double> xpts,ypts;
    std::vector<double> cx,cy;
    bool with_endpoint_values;
    double x0,x1,y0,y1;
  };
}


#endif

