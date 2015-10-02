#include "ParametricLegendre.hpp"
#include "../quad.hpp"
#include "../bmath.hpp"
#include <cstdio>


ccdl::ParametricLegendre::ParametricLegendre
( int const nbas, int const npt, int const naux )
  : nbas(nbas), npts(npt), naux(naux),
    ws(npt,0.), ts(npt,0.), xws(naux,0.), xts(naux,0.),
    xpts(npt,0.), ypts(npt,0.),
    cx(nbas,0.), cy(nbas,0.),
    with_endpoint_values(false),
    x0(0),x1(0),y0(0),y1(0)
{
  ccdl::GaussLegendreRule( 0., 1., npts,  ts.data(),  ws.data() );
  ccdl::GaussLegendreRule( 0., 1., naux, xts.data(), xws.data() );
}



void ccdl::ParametricLegendre::SetEndptValues
( double const X0, double const Y0, double const X1, double const Y1 )
{
  with_endpoint_values = true;
  x0=X0;
  y0=Y0;
  x1=X1;
  y1=Y1;
}

void ccdl::ParametricLegendre::SetEndptValues()
{
  with_endpoint_values = false;
}

// c.b - 1/2 c.A.c - u^T . ( P.c - h )
//   b[k] = \sum_i w(ti) pk(ti) x(ti)
//   A = unit matrix
//   where h[0] = x(t=0) and h[1] = x(t=1)
//         P[0+k*2] = pk(t=0)
//         P[1+k*2] = pk(t=1)
//  b - A.c - P^T.u = 0
//  b-P^T.u = A.c
//  b-P^T.u = c
//  P.b - P.P^T.u = h
//  P.b-h = P.P^T.u
// (P.P^T)^{-1}.(P.b-h) = u
//  
// c = A^{-1}.x, but A is the unit matrix 
// void ccdl::ParametricLegendre::SetX( double const * x )
// {
//   std::copy( x, x+n, cx.data() );
// }
// void ccdl::ParametricLegendre::SetX( double const * y )
// {
//   std::copy( y, y+n, cy.data() );
// }

void ccdl::ParametricLegendre::SetXY( double const * x, double const * y )
{
  std::copy( x, x+npts, xpts.data() );
  std::copy( y, y+npts, ypts.data() );

  std::fill( cx.data(), cx.data()+nbas, 0. );
  std::fill( cy.data(), cy.data()+nbas, 0. );
  std::vector<double> p(nbas,0.);
  if ( ! with_endpoint_values )
    {
      for ( int i=0; i<npts; ++i )
	{
	  double const t = ts[i];
	  double const w = ws[i];
	  ccdl::OrthonormalLegendreBasis( 0., 1., t, nbas, p.data() );
	  for ( int k=0; k<nbas; ++k )
	    {
	      cx[k] += w * x[i]*p[k];
	      cy[k] += w * y[i]*p[k];
	    }
	};
    }
  else
    {
      std::vector<double> A(nbas*nbas,0.);
      for ( int i=0; i<nbas; ++i ) A[i+i*nbas] = 1.;

      std::vector<double> ConMat(2*nbas,0.);
      ccdl::OrthonormalLegendreBasis( 0., 1., 0., nbas, p.data() );
      for ( int i=0; i<nbas; ++i ) ConMat[0+i*2] = p[i];
      ccdl::OrthonormalLegendreBasis( 0., 1., 1., nbas, p.data() );
      for ( int i=0; i<nbas; ++i ) ConMat[1+i*2] = p[i];

      std::vector<double> b( 2*nbas, 0. );
      for ( int i=0; i<npts; ++i )
	{
	  double const t = ts[i];
	  double const w = ws[i];
	  ccdl::OrthonormalLegendreBasis( 0., 1., t, nbas, p.data() );
	  for ( int k=0; k<nbas; ++k )
	    {
	      b[k+0]    += w * x[i]*p[k];
	      b[k+nbas] += w * y[i]*p[k];
	    }
	};

      std::vector<double> ConVal(2,0.);
     
      ConVal[0] = x0;
      ConVal[1] = x1;
      ccdl::ConstrainedLeastSquaresFit( nbas, nbas, A.data(), cx.data(), b.data()+0, 
					2, ConMat.data(), ConVal.data() );
      ConVal[0] = y0;
      ConVal[1] = y1;
      ccdl::ConstrainedLeastSquaresFit( nbas, nbas, A.data(), cy.data(), b.data()+nbas, 
					2, ConMat.data(), ConVal.data() );
    }
}


void ccdl::ParametricLegendre::GetXY( double const t, double & x, double & y ) const
{
  std::vector<double> p(nbas,0.);
  ccdl::OrthonormalLegendreBasis( 0., 1., t, nbas, p.data() );
  x = ccdl::v_dot_v( nbas, cx.data(), p.data() );
  y = ccdl::v_dot_v( nbas, cy.data(), p.data() );
}


void ccdl::ParametricLegendre::GetXY
( double const t, 
  double & x, double & y, 
  double & dxdt, double & dydt ) const
{
  std::vector<double> p(nbas,0.),dp(nbas,0.);
  ccdl::OrthonormalLegendreBasis( 0., 1., t, nbas, p.data(), dp.data() );
  x = ccdl::v_dot_v( nbas, cx.data(), p.data() );
  y = ccdl::v_dot_v( nbas, cy.data(), p.data() );
  dxdt = ccdl::v_dot_v( nbas, cx.data(), dp.data() );
  dydt = ccdl::v_dot_v( nbas, cy.data(), dp.data() );
}

void ccdl::ParametricLegendre::GetXY
( double const t, 
  double & x, double & y, 
  double & dxdt, double & dydt,
  double & d2xdt2, double & d2ydt2 ) const
{
  std::vector<double> p(nbas,0.),dp(nbas,0.),d2p(nbas,0.);
  ccdl::OrthonormalLegendreBasis( 0., 1., t, nbas, p.data(), dp.data(), d2p.data() );
  x = ccdl::v_dot_v( nbas, cx.data(), p.data() );
  y = ccdl::v_dot_v( nbas, cy.data(), p.data() );
  dxdt = ccdl::v_dot_v( nbas, cx.data(), dp.data() );
  dydt = ccdl::v_dot_v( nbas, cy.data(), dp.data() );
  d2xdt2 = ccdl::v_dot_v( nbas, cx.data(), d2p.data() );
  d2ydt2 = ccdl::v_dot_v( nbas, cy.data(), d2p.data() );

}


double ccdl::ParametricLegendre::PerpForceProj( double const t, double & dfdx, double & dfdy )
{
  double x,y,dxdt,dydt;
  GetXY( t,x,y,dxdt,dydt );
  double den = std::sqrt( dxdt*dxdt + dydt*dydt );
  double vx =  dydt / den;
  double vy = -dxdt / den;
  double fperp_mag = dfdx * vx + dfdy * vy;
  dfdx = fperp_mag * vx;
  dfdy = fperp_mag * vy;
  return fperp_mag;
}



double ccdl::ParametricLegendre::ParaForceProj( double const t, double & dfdx, double & dfdy )
{
  double x,y,dxdt,dydt;
  GetXY( t,x,y,dxdt,dydt );
  double den = std::sqrt( dxdt*dxdt + dydt*dydt );
  double vx =  dxdt / den;
  double vy =  dydt / den;
  double fperp_mag = dfdx * vx + dfdy * vy;
  dfdx = fperp_mag * vx;
  dfdy = fperp_mag * vy;
  return fperp_mag;
}


// naux
double ccdl::ParametricLegendre::PerpForceSumSq( double const * dfdx, double const * dfdy )
{
  double x,y,dxdt,dydt;
  double df2 = 0.;
  for ( int i=0; i<naux; ++i )
    {
      GetXY( xts[i],x,y,dxdt,dydt );
      double den = std::sqrt( dxdt*dxdt + dydt*dydt );
      double vx =  dydt / den;
      double vy = -dxdt / den;
      double fperp_mag = dfdx[i] * vx + dfdy[i] * vy;
      df2 += xws[i] * fperp_mag * fperp_mag;
    };
  return df2;
}


double ccdl::ParametricLegendre::ParaForceSumSq( double const * dfdx, double const * dfdy )
{
  double x,y,dxdt,dydt;
  double df2 = 0.;
  for ( int i=0; i<naux; ++i )
    {
      GetXY( xts[i],x,y,dxdt,dydt );
      double den = std::sqrt( dxdt*dxdt + dydt*dydt );
      double vx =  dxdt / den;
      double vy =  dydt / den;
      double fperp_mag = dfdx[i] * vx + dfdy[i] * vy;
      df2 += xws[i] * fperp_mag * fperp_mag;
    };
  return df2;
}

// double ccdl::ParametricLegendre::TotalForceSumSq( double const * dfdx, double const * dfdy )
// {
//   return PerpForceSumSq( dfdx, dfdy ) + ParaForceSumSq( dfdx, dfdy );
// }


double ccdl::ParametricLegendre::PerpForceSumSqCoefDer
( double const * dfdx, double const * dfdy,
  double * dSqdcx, double * dSqdcy )
{
  std::vector<double> p(nbas,0.),dp(nbas,0.);
  std::fill( dSqdcx, dSqdcx + nbas, 0. );
  std::fill( dSqdcy, dSqdcy + nbas, 0. );

  double df2 = 0.;
  for ( int i=0; i<naux; ++i )
    {
      ccdl::OrthonormalLegendreBasis( 0., 1., xts[i], nbas, p.data(), dp.data() );
      //double x = ccdl::v_dot_v( nbas, cx.data(), p.data() );
      //double y = ccdl::v_dot_v( nbas, cy.data(), p.data() );
      double dxdt = ccdl::v_dot_v( nbas, cx.data(), dp.data() );
      double dydt = ccdl::v_dot_v( nbas, cy.data(), dp.data() );

      double den = std::sqrt( dxdt*dxdt + dydt*dydt );


      double vx =  dydt / den;
      double vy = -dxdt / den;
      double fperp_mag = dfdx[i] * vx + dfdy[i] * vy;
      df2 += xws[i] * fperp_mag * fperp_mag;
      
      for ( int j=0; j<nbas; ++j )
	{
	  double ooden_x = (-0.5 / (den*den*den)) * ( 2.*dxdt * dp[j] );
	  double ooden_y = (-0.5 / (den*den*den)) * ( 2.*dydt * dp[j] );
	  double vx_x =  dydt*ooden_x;
	  double vx_y =  dydt*ooden_y + dp[j]/den;
	  double vy_y = -dxdt*ooden_y;
	  double vy_x = -dxdt*ooden_x - dp[j]/den;
	  double f2_x = xws[i] * 2. * fperp_mag * ( dfdx[i] * vx_x + dfdy[i] * vy_x );
	  double f2_y = xws[i] * 2. * fperp_mag * ( dfdx[i] * vx_y + dfdy[i] * vy_y );
	  dSqdcx[j] += f2_x;
	  dSqdcy[j] += f2_y;
	}
    };
  return df2;
}


double ccdl::ParametricLegendre::ParaForceSumSqCoefDer
( double const * dfdx, double const * dfdy,
  double * dSqdcx, double * dSqdcy )
{
  std::vector<double> p(nbas,0.),dp(nbas,0.);
  std::fill( dSqdcx, dSqdcx + nbas, 0. );
  std::fill( dSqdcy, dSqdcy + nbas, 0. );

  double df2 = 0.;
  for ( int i=0; i<naux; ++i )
    {
      ccdl::OrthonormalLegendreBasis( 0., 1., xts[i], nbas, p.data(), dp.data() );
      //double x = ccdl::v_dot_v( nbas, cx.data(), p.data() );
      //double y = ccdl::v_dot_v( nbas, cy.data(), p.data() );
      double dxdt = ccdl::v_dot_v( nbas, cx.data(), dp.data() );
      double dydt = ccdl::v_dot_v( nbas, cy.data(), dp.data() );

      double den = std::sqrt( dxdt*dxdt + dydt*dydt );
      double vx =  dxdt / den;
      double vy =  dydt / den;
      double fperp_mag = dfdx[i] * vx + dfdy[i] * vy;
      df2 += xws[i] * fperp_mag * fperp_mag;
      
      for ( int j=0; j<nbas; ++j )
	{
	  double ooden_x = (-0.5 / (den*den*den)) * ( 2.*dxdt * dp[j] );
	  double ooden_y = (-0.5 / (den*den*den)) * ( 2.*dydt * dp[j] );
	  double vx_y =  dxdt*ooden_y;
	  double vx_x =  dxdt*ooden_x + dp[j]/den;
	  double vy_x =  dydt*ooden_x;
	  double vy_y =  dydt*ooden_y + dp[j]/den;
	  double f2_x = xws[i] * 2. * fperp_mag * ( dfdx[i] * vx_x + dfdy[i] * vy_x );
	  double f2_y = xws[i] * 2. * fperp_mag * ( dfdx[i] * vx_y + dfdy[i] * vy_y );
	  dSqdcx[j] += f2_x;
	  dSqdcy[j] += f2_y;
	}
    };
  return df2;
}




double ccdl::ParametricLegendre::PerpForceSumSqPtDer
( double const * dfdx, double const * dfdy, 
  double * dndX, double * dndY )
{
  std::fill( dndX, dndX+npts, 0. );
  std::fill( dndY, dndY+npts, 0. );

  std::vector<double> dndcx(nbas,0.), dndcy(nbas,0.);
  double fperp2 = PerpForceSumSqCoefDer( dfdx, dfdy, dndcx.data(), dndcy.data() );

  std::vector<double> p(nbas,0.);

  if ( ! with_endpoint_values )
    {
      for ( int i=0; i<npts; ++i )
	{
	  double const t = ts[i];
	  double const w = ws[i];
	  ccdl::OrthonormalLegendreBasis( 0., 1., t, nbas, p.data() );
	  for ( int k=0; k<nbas; ++k )
	    {
	      //cx[k] += w * x[i]*p[k];
	      //cy[k] += w * y[i]*p[k];
	      dndX[i] += dndcx[k] * w * p[k];
	      dndY[i] += dndcy[k] * w * p[k];
	    }
	}
    }
  else
    {
      std::vector<double> D(2*nbas,0.);
      ccdl::OrthonormalLegendreBasis( 0., 1., 0., nbas, p.data() );
      for ( int i=0; i<nbas; ++i ) D[0+i*2] = p[i];
      ccdl::OrthonormalLegendreBasis( 0., 1., 1., nbas, p.data() );
      for ( int i=0; i<nbas; ++i ) D[1+i*2] = p[i];
      std::vector<double> DDt(4,0.);
      ccdl::ge_dot_gt( 2, nbas, D.data(), 2, nbas, D.data(), DDt.data() );
      ccdl::sdd_sym_power( -1., 2, 2, DDt.data() );
      std::vector<double> ID(2*nbas,0.);
      ccdl::ge_dot_ge( 2,2, DDt.data(), 2,nbas, D.data(), ID.data() );
      std::vector<double> dcdb(nbas*nbas,0.);
      ccdl::gt_dot_ge( 2,nbas, D.data(), 2,nbas, ID.data(), dcdb.data(), -1., 0. );
      
      for ( int i=0; i<npts; ++i )
	{
	  double const t = ts[i];
	  double const w = ws[i];
	  ccdl::OrthonormalLegendreBasis( 0., 1., t, nbas, p.data() );
	  for ( int k=0; k<nbas; ++k )
	    {
	      //cx[k] += w * x[i]*p[k];
	      //cy[k] += w * y[i]*p[k];
	      dndX[i] += dndcx[k] * w * p[k];
	      dndY[i] += dndcy[k] * w * p[k];
	      for ( int kp=0; kp<nbas; ++kp )
		{
		  dndX[i] += dndcx[k] * dcdb[k+kp*nbas] * w * p[kp];
		  dndY[i] += dndcy[k] * dcdb[k+kp*nbas] * w * p[kp];
		};
	    }
	}
    }
  return fperp2;
}


double ccdl::ParametricLegendre::ParaForceSumSqPtDer
( double const * dfdx, double const * dfdy, 
  double * dndX, double * dndY )
{
  std::fill( dndX, dndX+npts, 0. );
  std::fill( dndY, dndY+npts, 0. );

  std::vector<double> dndcx(nbas,0.), dndcy(nbas,0.);
  double fperp2 = ParaForceSumSqCoefDer( dfdx, dfdy, dndcx.data(), dndcy.data() );

  std::vector<double> p(nbas,0.);

  if ( ! with_endpoint_values )
    {
      for ( int i=0; i<npts; ++i )
	{
	  double const t = ts[i];
	  double const w = ws[i];
	  ccdl::OrthonormalLegendreBasis( 0., 1., t, nbas, p.data() );
	  for ( int k=0; k<nbas; ++k )
	    {
	      //cx[k] += w * x[i]*p[k];
	      //cy[k] += w * y[i]*p[k];
	      dndX[i] += dndcx[k] * w * p[k];
	      dndY[i] += dndcy[k] * w * p[k];
	    }
	}
    }
  else
    {
      std::vector<double> D(2*nbas,0.);
      ccdl::OrthonormalLegendreBasis( 0., 1., 0., nbas, p.data() );
      for ( int i=0; i<nbas; ++i ) D[0+i*2] = p[i];
      ccdl::OrthonormalLegendreBasis( 0., 1., 1., nbas, p.data() );
      for ( int i=0; i<nbas; ++i ) D[1+i*2] = p[i];
      std::vector<double> DDt(4,0.);
      ccdl::ge_dot_gt( 2, nbas, D.data(), 2, nbas, D.data(), DDt.data() );
      ccdl::sdd_sym_power( -1., 2, 2, DDt.data() );
      std::vector<double> ID(2*nbas,0.);
      ccdl::ge_dot_ge( 2,2, DDt.data(), 2,nbas, D.data(), ID.data() );
      std::vector<double> dcdb(nbas*nbas,0.);
      ccdl::gt_dot_ge( 2,nbas, D.data(), 2,nbas, ID.data(), dcdb.data(), -1., 0. );
      
      for ( int i=0; i<npts; ++i )
	{
	  double const t = ts[i];
	  double const w = ws[i];
	  ccdl::OrthonormalLegendreBasis( 0., 1., t, nbas, p.data() );
	  for ( int k=0; k<nbas; ++k )
	    {
	      //cx[k] += w * x[i]*p[k];
	      //cy[k] += w * y[i]*p[k];
	      dndX[i] += dndcx[k] * w * p[k];
	      dndY[i] += dndcy[k] * w * p[k];
	      for ( int kp=0; kp<nbas; ++kp )
		{
		  dndX[i] += dndcx[k] * dcdb[k+kp*nbas] * w * p[kp];
		  dndY[i] += dndcy[k] * dcdb[k+kp*nbas] * w * p[kp];
		};
	    }
	}
    }
  return fperp2;
}


// double ccdl::ParametricLegendre::TotalForceSumSqPtDer
// ( double const * dfdx, double const * dfdy, 
//   double * dndX, double * dndY )
// {
//   std::vector<double> dx( npts, 0. ), dy( npts, 0. );
//   double f = PerpForceSumSqPtDer( dfdx, dfdy, dndX, dndY );
//   f += ParaForceSumSqPtDer( dfdx, dfdy, dx.data(), dy.data() );
//   for ( int i=0; i<npts; ++i )
//     {
//       dndX[i] += dx[i];
//       dndY[i] += dy[i];
//     };
//   return f;
// }
