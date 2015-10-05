#include <cstdio>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <tr1/memory>
#include <cstdlib>
#include <getopt.h>

#include <ccdl/spline.hpp>
#include <ccdl/constants.hpp>

#include <nlopt.hpp>



struct endpt
{
  endpt() : x(0), y(0), i(0), minimum(false) {}
  double x,y;
  int i;
  bool minimum;
  std::string GetName() const;
};

std::string endpt::GetName() const
{
  std::stringstream str;
  str << ( minimum ? "min" : "ts" ) << i;
  return str.str();
}

struct connection
{
  connection( endpt t0, endpt t1, int method = 1 ) : t0(t0), t1(t1), method(method), prev(-1) {};
  std::string GetName() const;
  endpt t0,t1;
  int method;
  int prev;
  std::tr1::shared_ptr< ccdl::ParametricLegendre > curve;
};

std::string connection::GetName() const
{
  std::stringstream str;
  str << t0.GetName() << "_" << t1.GetName();
  return str.str();
}











struct cli_options
{
  cli_options();
  std::vector< endpt > minima,maxima,saddle;
  std::vector< connection > connections;
  std::string meshfile;
  bool rezero;
  bool perx,pery;
};


cli_options::cli_options()
  : rezero(false),perx(false), pery(false)
{
  ////////////////////////////////////////
  {
    std::ifstream cin;
    cin.open( "minpts.dat" );
    if ( cin.good() )
      {
	std::string line;
	std::stringstream str;
	while ( std::getline(cin,line) )
	  {
	    if ( line.empty() ) continue; 
	    str.str( line );
	    str.clear();
	    endpt p;
	    p.minimum = true;
	    p.i = minima.size() + 1;
	    str >> p.x >> p.y;
	    minima.push_back( p );
	  }
      }
  }
  ////////////////////////////////////////
  {
    std::ifstream cin;
    cin.open( "tspts.dat" );
    if ( cin.good() )
      {
	std::string line;
	std::stringstream str;
	while ( std::getline(cin,line) )
	  {
	    if ( line.empty() ) continue; 
	    str.str( line );
	    str.clear();
	    endpt p;
	    p.minimum = false;
	    p.i = saddle.size() + 1;
	    str >> p.x >> p.y;
	    saddle.push_back( p );
	  }
      }
  }
  ////////////////////////////////////////
  {
    std::ifstream cin;
    cin.open( "maxpts.dat" );
    if ( cin.good() )
      {
	std::string line;
	std::stringstream str;
	while ( std::getline(cin,line) )
	  {
	    if ( line.empty() ) continue; 
	    str.str( line );
	    str.clear();
	    endpt p;
	    p.minimum = false;
	    p.i = maxima.size() + 1;
	    str >> p.x >> p.y;
	    maxima.push_back( p );
	  }
      }
  }
  ////////////////////////////////////////
  
}



void print_usage()
{
  std::printf("Usage: pmf2d [options] mesh\n\n");
  std::printf("Required arguments:\n");
  std::printf("  mesh\t\ta 2d mesh file of the form:\n");
  std::printf("\t\t   x(1) y(1) value(1,1)\n");
  std::printf("\t\t   x(2) y(1) value(2,1)\n");
  std::printf("\t\t       [...etc...]\n");
  std::printf("\t\t   x(4) y(1) value(N,1)\n");
  std::printf("\t\t      [a blank line]\n");
  std::printf("\t\t   x(1) y(2) value(1,2)\n");
  std::printf("\t\t   x(2) y(2) value(2,2)\n");
  std::printf("\t\t       [...etc...]\n");
  std::printf("\t\t   x(N) y(2) value(N,2)\n");
  std::printf("\t\t      [a blank line]\n");
  std::printf("\t\t   [...and so forth...]\n");
  std::printf("\n");
  std::printf("Options:\n");
  std::printf("  --zero\tif true, then write a new mesh file called mesh.0.dat\n");
  std::printf("          \tthat shifts all values such that the lowest minimum is 0\n");
  std::printf("  --con min=i\tendpoint set to the i'th value found in minpts.dat\n");
  std::printf("  --con ts=i\tendpoint set to the i'th value found in tspts.dat\n");
  std::printf("  --con method=i neb method to use (1 or 2) (default 1)\n");
  std::printf("           \tmethod=1 tends to be better for short paths\n");
  std::printf("         \tmethod=2 tends to be better for long paths\n");
  std::printf("\n");
  std::printf("\n\nA typical work-flow:\n");
  std::printf("  (1) pmf2d --zero mesh\n");
  std::printf("  (2) /bin/sh ./mesh.0.dat.sh\n");
  std::printf("  (3) evince ./mesh.0.dat.eps\n");
  std::printf("  (4) pmf2d --con min=1,ts=1,min=2 mesh.0.dat\n");
  std::printf("  (5) /bin/sh ./mesh.0.dat.sh\n");
  std::printf("  (6) evince ./mesh.0.dat.eps\n");
  std::printf("  (7) xmgrace min1_ts1.1d.dat ts1_min2.1d.dat\n");
  std::printf("Comments:\n");
  std::printf("  (1) pmf2d --zero mesh\n");
  std::printf("      Creates the files minpts.dat tspts.dat maxpts.dat\n");
  std::printf("      mesh.0.dat.sh mesh.0.dat\n");
  std::printf("      The minpts, tspts and maxpts files have the following format:\n");
  std::printf("         x y # f dfdx dfdy freq1 freq2\n");
  std::printf("      where\n");
  std::printf("         x y\tlocation\n");
  std::printf("         f\tfunction value @ x,y\n");
  std::printf("         dfdx dfdy\tfunction gradient @ x,y\n");
  std::printf("         freq1 freq2\tHessian eigenvalues @ x,y\n");
  std::printf("      The mesh.0.dat is a rewrite of mesh with the minimum shifted to zero\n");
  std::printf("  (2) /bin/sh ./mesh.0.dat.sh\n");
  std::printf("  (3) evince ./mesh.0.dat.eps\n");
  std::printf("      This will let you see the contour with red and blue numerals\n");
  std::printf("      The blue numerals are the minima in minpts.dat and the red\n");
  std::printf("      numerals are the points in tspts.dat\n");
  std::printf("  (4) pmf2d --con min=1,ts=1,min=2 mesh.0.dat\n");
  std::printf("      Overwrites mesh.0.dat.sh, such that it will create\n");
  std::printf("      a plot with a curve connecting the first and second minimum\n");
  std::printf("      through the first transition state\n");
  std::exit(0);
}



cli_options read_options( int argc, char ** argv )
{
  cli_options cli;

  char const * CON_opts[] = {
#define CON_MIN 0
    "min",
#define CON_TS 1
    "ts",
#define CON_METHOD 2
    "method",
    NULL };

  static struct option long_options[] =
    {
      { "help",      no_argument,       NULL, 'h'   },
      { "zero",      no_argument,       NULL, 'z'   },
      { "perx",      no_argument,       NULL, 'x'   },
      { "pery",      no_argument,       NULL, 'y'   },
#define CON     0100
      { "con",       required_argument, NULL, CON   },
      {NULL,0,NULL,0}
    };

  int opt = 0;
  int long_index = 0;
  char * subopts, * value;
  while ( (opt = getopt_long
          ( argc, argv, "hxy", 
            long_options, &long_index )) != -1 )
    {
      
      switch (opt)
        {
	  
        case 'h':     { print_usage(); break; }
	case 'z':     { cli.rezero = true; break; }
        case 'x':     { cli.perx = true; break; }
        case 'y':     { cli.pery = true; break; }
	case CON:
	  {
	    std::vector<endpt> ts;
	    int method = 1;
            subopts = optarg;
            while (*subopts != '\0')
              switch (getsubopt(&subopts, const_cast<char**>(CON_opts), &value))
                {
                case CON_MIN:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --con min= expects an int; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    int i = std::atoi( value );
		    if ( i < 1 or i > (int)(cli.minima.size()) )
		      std::printf("%s: --con min=%i out of bounds\n",argv[0],i);
		    ts.push_back( cli.minima[ i-1 ] );
		    break;
		  }
		case CON_TS:
 		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --con ts= expects an int; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    int i = std::atoi( value );
		    if ( i < 1 or i > (int)(cli.saddle.size()) )
		      std::printf("%s: --con ts=%i out of bounds\n",argv[0],i);
		    ts.push_back( cli.saddle[ i-1 ] );
		    break;
		  }
		case CON_METHOD:
 		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --con method= expects an int; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    int i = std::atoi( value );
		    if ( i < 1 or i > 2 )
		      std::printf("%s: --con method=%i out of bounds\n",argv[0],i);
		    method = i;
		    break;
		  }
		default:
		  {
		    std::printf("%s: Unknown subption '%s' given to --con\n",argv[0],value);
		    std::exit(EXIT_FAILURE);
		  }
		}
	    if ( ts.size() < 2 )
	      {
		std::printf("%s: --con expects two endpoints, but only 1 provided\n",argv[0]);
		std::exit(EXIT_FAILURE);
	      }
	    else
	      {
		for ( std::size_t i=1; i<ts.size(); ++i )
		  {
		    connection c( ts[i-1],ts[i] );
		    c.method = method;
		    c.prev = (int)(cli.connections.size())-1;
		    cli.connections.push_back( c );
		  };
	      }
	    break;
	  }
	case '?':
	  {
	    std::printf("%s: use -h for usage\n",argv[0]);
	    std::exit(EXIT_FAILURE);
	    break;
	  }
	default:
	  std::printf("%s: an command line input error occured\n",argv[0]);
	  std::exit(EXIT_FAILURE);
	};

    }

  if ( optind == argc )
    {
      std::printf("%s: missing non-optional argument -- mesh file\n",argv[0]);
      std::printf("%s: use -h for usage\n",argv[0]);
      std::exit(EXIT_FAILURE);
    }
  else if ( argc - optind > 1 )
    {
      std::printf("%s: expected two nonoptional args, but received %i args\n",argv[0],argc-optind); 
      std::printf("%s:",argv[0]);
      while (optind < argc)
	std::printf(" %s", argv[optind++]);
      std::printf("\n");
      std::exit(EXIT_FAILURE);
    }
  else
    {
      cli.meshfile = argv[optind++];
    }

  return cli;
}








struct NudgeElasticBandOpt 
{
  NudgeElasticBandOpt( ccdl::Mesh2d * pmesh, ccdl::ParametricLegendre * pcurve );
  double operator() ( double const * x );
  ccdl::Mesh2d * pmesh;
  ccdl::ParametricLegendre * pcurve;
  double SpringWt;
  double CurvatureWt;
  double EnergyWt;
  double MaximumWt;
  double Temp;
};

NudgeElasticBandOpt::NudgeElasticBandOpt
( ccdl::Mesh2d * mesh, 
  ccdl::ParametricLegendre * curve )
  : pmesh(mesh),
    pcurve(curve),
    SpringWt( 10. ),
    CurvatureWt( 0.1 ),
    EnergyWt( 1.0 ),
    MaximumWt( 100. ),
    Temp( 1. )
{
}

double NudgeElasticBandOpt::operator() ( double const * c )
{
  int npts = pcurve->GetNpts();
  pcurve->SetXY( c, c + npts );

  int naux = pcurve->GetNauxpts();
  std::vector<double> dfdx( naux, 0. ), dfdy( naux, 0. );
  double dtpen = 0.;
  double dxpen = 0.;
  double fsum = 0.;
  double arclen = 0.;

  std::vector< ccdl::Mesh2dHessian > vs( naux );
  ccdl::Mesh2dHessian vmin,vmax;
  for ( int i=0; i<naux; ++i )
    {
      double t = pcurve->GetAuxPt(i);
      double x,y;
      pcurve->GetXY( t, x, y );
      vs[i] = pmesh->GetHessian( x, y );
      if ( i == 0 or vs[i].f < vmin.f ) vmin = vs[i];
      if ( i == 0 or vs[i].f > vmax.f ) vmax = vs[i];
    }


  double beta = ccdl::AU_PER_KCAL_PER_MOL / ( Temp * ccdl::BOLTZMANN_CONSTANT_AU );
  std::vector<double> minw( naux, 0. ), maxw( naux, 0. );
  double minwsum = 0.;
  double maxwsum = 0.;
  for ( int i=0; i<naux; ++i )
    {
      minw[i] = std::exp( -beta * ( vs[i].f-vmin.f ) );
      maxw[i] = std::exp(  beta * ( vs[i].f-vmax.f ) );
      minwsum += minw[i];
      maxwsum += maxw[i];
    };
  double avgmaxw = 0.;
  for ( int i=0; i<naux; ++i )
    {
      minw[i] /= minwsum;
      maxw[i] /= maxwsum;
      avgmaxw += maxw[i];
    }
  avgmaxw /= naux;

  double fall=0.;
  for ( int i=0; i<naux; ++i )
    {
      double t = pcurve->GetAuxPt(i);
      double w = pcurve->xws[i];
      double x,y,dxdt,dydt,d2xdt2,d2ydt2;
      pcurve->GetXY( t, x, y, dxdt, dydt, d2xdt2, d2ydt2 );
      ccdl::Mesh2dHessian v =  vs[i];

      fall += w * v.f;
      fsum   += w * maxw[i] * v.f;
      arclen += w * std::sqrt( v.dfdx*v.dfdx + v.dfdy*v.dfdy );
      //double dt = 1.;
      //if ( i > 0 and i < naux-1 ) dt = (pcurve->GetAuxPt(i+1)-pcurve->GetAuxPt(i-1))/2.;
      //else if ( i > 0 ) dt = (pcurve->GetAuxPt(i)-pcurve->GetAuxPt(i-1));
      //else dt = (pcurve->GetAuxPt(i+1)-pcurve->GetAuxPt(i));

      //dtpen  += w * std::sqrt( ( d2xdt2*d2xdt2 + d2ydt2*d2ydt2 ) * std::pow(dt,2) );

      // perpendicular curvative
      double perp[2] = { dydt, -dxdt };
      double perpnrm = std::sqrt( perp[0]*perp[0] + perp[1]*perp[1] );
      perp[0] /= perpnrm;
      perp[1] /= perpnrm;
      
      double hp[2] = { v.h[0+0*2] * perp[0] + v.h[0+1*2] * perp[1],
		       v.h[1+0*2] * perp[0] + v.h[1+1*2] * perp[1] };
      double php = perp[0] * hp[0] + perp[1] * hp[1];
      if ( php < 0. ) dtpen += w * php*php;
      
      // double para[2] = { dxdt / perpnrm, dydt / perpnrm };
      // double hpp[2] = { v.h[0+0*2] * para[0] + v.h[0+1*2] * para[1],
      // 			v.h[1+0*2] * para[0] + v.h[1+1*2] * para[1] };
      // double php2 = para[0] * hpp[0] + para[1] * hpp[1];
      // if ( php2 > 0. ) dtpen += w * php2*php2;

      //double para[2] = { dxdt / perpnrm, dydt / perpnrm };
      //double uperp = v.evec[0+1*2] * perp[0] + v.evec[1+1*2] * perp[1];
      //double upara = v.evec[0+0*2] * para[0] + v.evec[1+0*2] * para[1];
      //dtpen += w * ( std::pow( uperp-upara, 2 ) - 1 );
      //std::printf("upp %20.10e %20.10e\n",uperp,upara);

      /*
      // parallel curvative
      double para[2] = { dxdt / perpnrm, dydt / perpnrm };
      para[0] /= perpnrm;
      para[1] /= perpnrm;
      double hpp[2] = { v.h[0+0*2] * para[0] + v.h[0+1*2] * para[1],
		       v.h[1+0*2] * para[0] + v.h[1+1*2] * para[1] };
      double php2 = para[0] * hpp[0] + para[1] * hpp[1];
      */
      // if ( v.eval[0] < 0. and v.eval[1] < 0. )
      // 	dtpen += w * (v.eval[0]*v.eval[0] + v.eval[1]*v.eval[1]);
      

      //dtpen += w * php2*php2;//std::exp( -0.01 * php );
      //std::printf("php2 %20.10e\n",php2);

      dfdx[i] = v.dfdx;
      dfdy[i] = v.dfdy;
    };
  // for ( int i=1; i<npts; ++i )
  //   {
  //     double x0,y0,x1,y1;
  //     double t1 = pcurve->GetQuadPt(i);
  //     double t0 = pcurve->GetQuadPt(i-1);
  //     pcurve->GetXY( t0, x0, y0 );
  //     pcurve->GetXY( t1, x1, y1 );
  //     double dx = x1-x0;
  //     double dy = y1-y0;
  //     double r2 = dx*dx+dy*dy;
  //     ccdl::Mesh2dValue v0 =  pmesh->GetValue( x0, y0 );
  //     ccdl::Mesh2dValue v1 =  pmesh->GetValue( x1, y1 );
  //     double pf1 = pcurve->ParaForceProj( t0, v0.dfdx, v0.dfdy );
  //     double pf2 = pcurve->ParaForceProj( t1, v1.dfdx, v1.dfdy );
  //     dxpen += r2 ;//* 100. * std::exp( -0.01 * ( pf1*pf1 + pf2*pf2 ) );
  //   };

  /*
  double X0,Y0,XN,YN;
  pcurve->GetXY( 0, X0, Y0 );
  pcurve->GetXY( 1, XN, YN );
  double DX = X0-XN;
  double DY = Y0-YN;
  double len2 = DX*DX+DY*DY;
  for ( int i=0; i<npts; ++i )
    for ( int j=0; j<i; ++j )
      {
	double x0,y0,x1,y1;
	double t1 = pcurve->GetQuadPt(i);
	double t0 = pcurve->GetQuadPt(j);
	pcurve->GetXY( t0, x0, y0 );
	pcurve->GetXY( t1, x1, y1 );
	double dx = x1-x0;
	double dy = y1-y0;
	double r2 = dx*dx+dy*dy;

	double r = std::sqrt(r2);
	double l = std::sqrt(len2);
	double percent = std::abs(i-j) / ( (double)(npts-1) );
	dxpen += 10.*std::pow( r/l - percent , 2 );
      };
  */

  {
    std::vector<double> gx(npts,0.), gy(npts,0.), cosa(npts,0.);
    double K = 1. / npts;
    for ( int i=0; i<npts; ++i )
      {
    	double x0,y0,xl,yl,xh,yh;
	double dx0,dy0,dx20,dy20;
    	double t0 = pcurve->GetQuadPt(i);
    	double th = 1.;
    	double tl = 0.;

	if ( i < npts-1 )
	  th = pcurve->GetQuadPt(i+1);
	if ( i > 0 )
	  tl = pcurve->GetQuadPt(i-1);

    	pcurve->GetXY( t0, x0, y0, dx0,dy0, dx20,dy20 );
    	pcurve->GetXY( tl, xl, yl );
    	pcurve->GetXY( th, xh, yh );
	//double h = std::sqrt( dx20*dx20 + dy20*dy20 );

	//double endpt = 1.;
	//if ( i < 1 or i > npts-2 ) endpt = 0.;
	//gx[i] = endpt * K * h * dx0 * std::pow( th-tl, 3 );
	//gy[i] = endpt * K * h * dy0 * std::pow( th-tl, 3 );

	double dxh = xh-x0;
	double dyh = yh-y0;
	double dxl = x0-xl;
	double dyl = y0-yl;

	gx[i] += K * dxh - K * dxl;
	gy[i] += K * dyh - K * dyl;

	double rh = std::sqrt( dxh*dxh+dyh*dyh );
	double rl = std::sqrt( dxl*dxl+dyl*dyl );

	cosa[i] = ( dxh*dxl + dyh*dyl ) / ( rh*rl );
	//if ( cosa[i] < 0. ) { dxpen += 10000. * std::pow(cosa[i],2); }
	if ( cosa[i] < 0. ) cosa[i] = 1.;
	else cosa[i] = 0.5*(1+std::cos(ccdl::PI*cosa[i]));
	//std::printf("%5i %20.10f %20.10f\n",i,cosa[i],0.5*(1+std::cos(ccdl::PI*cosa[i])));
      };
    for ( int i=0; i<npts; ++i )
      {
	double t = pcurve->GetQuadPt(i);
	double x,y;
	pcurve->GetXY( t, x, y );
	double pgx = gx[i];
	double pgy = gy[i];
	//double fpara = pcurve->ParaForceProj( t, pgx, pgy );
	// pgx and pgy are the parallel spring gradients 
	pgx = gx[i] * cosa[i] + pgx * ( 1.-cosa[i] );
	pgy = gy[i] * cosa[i] + pgy * ( 1.-cosa[i] );
	dxpen += pgx*pgx + pgy*pgy + cosa[i]*cosa[i];
      }

  }



  /*
  {
    std::vector<double> gx(npts,0.), gy(npts,0.);
    double K = 1.;
    double X0,Y0,XN,YN;
    pcurve->GetXY( 0, X0, Y0 );
    pcurve->GetXY( 1, XN, YN );
    double DX = X0-XN;
    double DY = Y0-YN;
    double len2 = DX*DX+DY*DY;
    for ( int i=0; i<npts; ++i )
      for ( int j=0; j<i; ++j )
	{
	  double x0,y0,x1,y1;
	  double t1 = pcurve->GetQuadPt(i);
	  double t0 = pcurve->GetQuadPt(j);
	  pcurve->GetXY( t0, x0, y0 );
	  pcurve->GetXY( t1, x1, y1 );
	  double dx = x1-x0;
	  double dy = y1-y0;
	  double r2 = dx*dx+dy*dy;
	  
	  double r = std::sqrt(r2);
	  double l = std::sqrt(len2);
	  double percent = std::abs(i-j) / ( (double)(npts-1) );
	  double du = r/l - percent;
	  // pen = du*du;
	  gx[i] += K*2*du * (dx/r) / l;
	  gy[i] += K*2*du * (dy/r) / l;
	  gx[j] -= K*2*du * (dx/r) / l;
	  gy[j] -= K*2*du * (dy/r) / l;
      };
  */

    /*

    }
    */

  double f = pcurve->PerpForceSumSq( dfdx.data(), dfdy.data() );

  // for ( int i=0; i<npts; ++i )
  //   {
  //     double t = pcurve->GetQuadPt( i );
  //     double x,y;
  //     pcurve->GetXY( t, x, y );
  //     std::printf("%20.10f %20.10f  (%20.10f + %20.10f) %20.10f %20.10f\n",c[i],c[i+npts],f,dtpen,x,y);
  //   }
  // std::printf("\n");
  //std::printf("f %20.8f  max %20.10f  spring %20.10f curve %20.10f\n",f,fsum,dxpen,dtpen);

  //std::printf("%20.10e\n",avgmaxw);
  //fsum = 0.;



  return f + EnergyWt * fall + MaximumWt * fsum + SpringWt * dxpen + (CurvatureWt) * dtpen;
  //return  (CurvatureWt) * dtpen;
  //return fall + SpringWt * dxpen;
}



double NLOPT_Chisq
( std::vector<double> const & x, 
std::vector<double> & grad, void * pin )
{

  NudgeElasticBandOpt * po = (NudgeElasticBandOpt *)( pin );
  NudgeElasticBandOpt & o = *po;

  double f = 0.;
  if ( grad.size() > 0 )
    {
      std::cerr << "NLOPT_Chisq does not support gradient-based minimizations\n";
      std::abort();
    }
  else
    f = o( x.data() );

  return f;
}





void NudgeElasticBand
( ccdl::Mesh2d & mesh, ccdl::ParametricLegendre & curve,
  double Kwt, double Cwt, double Ewt, double TSwt, double Temp )
{
  NudgeElasticBandOpt o( &mesh, &curve );
  o.Temp = Temp;
  o.MaximumWt = TSwt;
  o.EnergyWt = Ewt;
  o.SpringWt  = Kwt;
  o.CurvatureWt = Cwt;

  int npts = curve.GetNpts();
  int n = npts*2;

  nlopt::opt nlo( nlopt::LN_COBYLA, n );
  nlo.set_min_objective( &NLOPT_Chisq, &o );
  nlo.set_xtol_rel( 1.e-4 );
  nlo.set_ftol_abs( 1.e-4 );
  nlo.set_maxeval( 1000 );
  std::vector<double> s( n, 0.05 );
  nlo.set_initial_step( s );
  double const * x = curve.GetXpts();
  double const * y = curve.GetYpts();
  std::copy( x, x+npts, s.data() );
  std::copy( y, y+npts, s.data()+curve.GetNpts() );

  double chisq = 0.;
  //nlopt::result res = 
  nlo.optimize( s, chisq );
}







//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////




ccdl::ParametricLegendre
NudgeElasticBand_FitMethod1
( ccdl::Mesh2d & mesh, 
  double const x0, double const y0, 
  double const x1, double const y1 )
{
  double SpringWt          = 1.;
  double EnergyWt          = 0.01; 
  double CurvatureWt       = 0.01;
  double TransitionStateWt = 0.01;
  double Temperature       = 0.01;


  ccdl::ParametricLegendre curve( 3, 3, 100 );
  curve.SetEndptValues( x0,y0, x1,y1 );
  int cnt = 0;
  int orders[6] = { 3, 5, 8, 12, 15, 18 };
  for ( int iorder=0; iorder < 5; iorder++, ++cnt )
    {
      double TW = TransitionStateWt;
      double EW = EnergyWt;
      double SW = SpringWt;

      //if ( iorder < 3 ) { EW = 0.05; SW = 50.; TW = 250.; }
      //if ( iorder < 1 ) { EW = 0.1; SW = 100.; TW = 500.; }


      int order = orders[iorder];
      ccdl::ParametricLegendre new_curve( order, order, 100 );
      new_curve.SetEndptValues( x0,y0, x1,y1 );
      int npt = new_curve.GetNpts();
      std::vector<double> xt(npt,0.), yt(npt,0.); 
      if ( cnt == 0 )
	{
	  double mx = (x1-x0);
	  double bx = x0;
	  double my = (y1-y0);
	  double by = y0;
	  for ( int i=0; i<npt; ++i )
	    {
	      double t = new_curve.GetQuadPt(i);
	      xt[i] = mx * t + bx;
	      yt[i] = my * t + by;
	    }
	}
      else
	{
	  for ( int i=0; i<npt; ++i )
	    {
	      double t = new_curve.GetQuadPt(i);
	      curve.GetXY( t, xt[i], yt[i] );
	    }
	}
      new_curve.SetXY( xt.data(), yt.data() );
      NudgeElasticBand( mesh, new_curve, SW, CurvatureWt, EW, TW, Temperature );
      curve = new_curve;
    };


  return curve;
}



ccdl::ParametricLegendre
NudgeElasticBand_FitMethod2
( ccdl::Mesh2d & mesh, 
  double const x0, double const y0, 
  double const x1, double const y1 )
{
  double SpringWt          = 1.;
  double EnergyWt          = 0.05; 
  double CurvatureWt       = SpringWt * 0.01;
  double TransitionStateWt = EnergyWt * 10000.;
  double Temperature       = 1.;


  ccdl::ParametricLegendre curve( 3, 3, 100 );
  curve.SetEndptValues( x0,y0, x1,y1 );
  int cnt = 0;
  int orders[6] = { 3, 5, 8, 12, 15, 18 };
  for ( int iorder=0; iorder < 5; iorder++, ++cnt )
    {
      double TW = TransitionStateWt;
      double EW = EnergyWt;
      double SW = SpringWt;

      if ( iorder < 3 ) { EW = 0.05; SW = 50.; TW = 250.; }
      if ( iorder < 2 ) { EW = 0.1; SW = 100.; TW = 500.; }

      int order = orders[iorder];
      ccdl::ParametricLegendre new_curve( order, order, 100 );
      new_curve.SetEndptValues( x0,y0, x1,y1 );
      int npt = new_curve.GetNpts();
      std::vector<double> xt(npt,0.), yt(npt,0.); 
      if ( cnt == 0 )
	{
	  double mx = (x1-x0);
	  double bx = x0;
	  double my = (y1-y0);
	  double by = y0;
	  for ( int i=0; i<npt; ++i )
	    {
	      double t = new_curve.GetQuadPt(i);
	      xt[i] = mx * t + bx;
	      yt[i] = my * t + by;
	    }
	}
      else
	{
	  for ( int i=0; i<npt; ++i )
	    {
	      double t = new_curve.GetQuadPt(i);
	      curve.GetXY( t, xt[i], yt[i] );
	    }
	}
      new_curve.SetXY( xt.data(), yt.data() );
      NudgeElasticBand( mesh, new_curve, SW, CurvatureWt, EW, TW, Temperature );
      curve = new_curve;
    };


  return curve;
}



ccdl::ParametricLegendre
NudgeElasticBand_FitMethod1
( ccdl::Mesh2d & mesh, 
  ccdl::Mesh2dValue const & t0, ccdl::Mesh2dValue const & t1 )
{
  return NudgeElasticBand_FitMethod1( mesh, t0.x, t0.y, t1.x, t1.y );
}

ccdl::ParametricLegendre
NudgeElasticBand_FitMethod2
( ccdl::Mesh2d & mesh, 
  ccdl::Mesh2dValue const & t0, ccdl::Mesh2dValue const & t1 )
{
  return NudgeElasticBand_FitMethod2( mesh, t0.x, t0.y, t1.x, t1.y );
}








//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////




struct GradOpt 
{
  GradOpt( ccdl::Mesh2d * pmesh );
  double operator() ( double const * x );
  ccdl::Mesh2d * pmesh;
};

GradOpt::GradOpt
( ccdl::Mesh2d * mesh )
  : pmesh(mesh)
{
}


double GradOpt::operator() ( double const * x )
{
  ccdl::Mesh2dValue h( pmesh->GetValue( x[0], x[1] ) );
  return h.dfdx*h.dfdx + h.dfdy*h.dfdy;
}


double NLOPT_GradMinimize( std::vector<double> const & s, std::vector<double> & g, void * p )
{
  if ( g.size() > 0 )
    {
      std::cerr << "NLOPT_GradMinimize does not support gradients\n";
      std::abort();
    }
  return (*(GradOpt *)p)( s.data() );
}

// void GradMinimize
// ( ccdl::Mesh2d & mesh, double & x, double & y )
// {
//   GradOpt o( &mesh );

//   nlopt::opt nlo( nlopt::LN_COBYLA, 2 );
//   nlo.set_min_objective( &NLOPT_GradMinimize, &o );
//   nlo.set_xtol_rel( 1.e-5 );
//   nlo.set_ftol_abs( 1.e-5 );
//   nlo.set_maxeval( 1000 );
//   std::vector<double> s( 2, 0.005 );
//   nlo.set_initial_step( s );

//   s[0] = x;
//   s[1] = y;
//   double chisq = 0.;
//   nlopt::result res = nlo.optimize( s, chisq );
//   x = s[0];
//   y = s[1];
// }


bool sort_by_position( ccdl::Mesh2dHessian const & a, ccdl::Mesh2dHessian const & b )
{
  double dx = std::abs(a.x-b.x);
  bool xlt = dx < 1.e-4;
  return (!xlt) ? a.x < b.x : a.y < b.y;
}

/*
std::vector< ccdl::Mesh2dHessian >
GetTransitionStates( ccdl::Mesh2d & mesh, bool periodic )
{
  int nx = mesh.GetSizeX();
  int ny = mesh.GetSizeY();
  
  std::vector< ccdl::Mesh2dHessian > balls;
  std::vector< ccdl::Mesh2dHessian > vals(nx*ny);
  for ( int i=0; i<nx; ++i )
    for ( int j=0; j<ny; ++j )
      vals[j+i*ny] = mesh.GetHessian( mesh.GetX(i), mesh.GetY(j) );
        

  double xlow = mesh.GetLowX();
  double ylow = mesh.GetLowY();
  double lx = mesh.GetHighX() - xlow;
  double ly = mesh.GetHighY() - ylow;


  int const nxn = periodic ? nx : nx-1;
  int const nyn = periodic ? ny : ny-1;
  int const nx0 = periodic ? 0 : 1;
  int const ny0 = periodic ? 0 : 1;

  for ( int i=nx0; i<nxn; ++i )
    for ( int j=ny0; j<nyn; ++j )
      {
        int i1 = mesh.GetIndex( i+1, nx );
        int j1 = mesh.GetIndex( j+1, ny );

        ccdl::Mesh2dHessian & ll = vals[j +i *ny];
        ccdl::Mesh2dHessian & lr = vals[j +i1*ny];
        ccdl::Mesh2dHessian & ul = vals[j1+i *ny];
        ccdl::Mesh2dHessian & ur = vals[j1+i1*ny];

	
        double g1 = ll.dfdx * lr.dfdx ;
        double g2 = ul.dfdx * ur.dfdx ;
        double g3 = ll.dfdy * ul.dfdy ;
        double g4 = lr.dfdy * ur.dfdy ;
        int nsatisfied = (g1<=0) + (g2<=0) + (g3<=0) + (g4<=0);
	


        if ( nsatisfied > 0 )
          {
	    if ( ll.eval[0] < 0. and lr.eval[1] < 0. and ul.eval[0] < 0. and ur.eval[0] < 0. )
	      {
		std::printf("exr %8.3f %8.3f  %13.4e %13.4e %13.4e %13.4e\n",
			    ll.x,
			    ll.y,
			    ll.eval[1], lr.eval[1], ul.eval[1], ur.eval[1] );
	      };
            // there is an extremum in this quadrant

              }
          }


  std::sort( balls.begin(), balls.end(), sort_by_position );

  typedef std::vector< ccdl::Mesh2dHessian >::iterator iter;
  for ( iter p = balls.begin(); p != balls.end(); ++p )
    for ( iter q = p+1; q != balls.end(); )
      {
        double dx = std::abs(p->x - q->x);
        double dy = std::abs(p->y - q->y);
        if ( dx < 0.01 * lx and
             dy < 0.01 * ly )
          q = balls.erase( q );
        else
          ++q;
      }

  if ( ! periodic )
    for ( iter p = balls.begin(); p != balls.end(); )
      if ( std::abs( p->x - xlow ) / lx < 0.02 or
           std::abs( p->y - ylow ) / lx < 0.02 or
           std::abs( p->x - (xlow+lx) ) / lx < 0.02 or
           std::abs( p->y - (ylow+ly) ) / ly < 0.02 )
        p = balls.erase( p );
      else
        ++p;

  return balls;
}

*/










std::vector< ccdl::Mesh2dHessian > LocStationaryPts( ccdl::Mesh2d & mesh, bool periodic )
{

  GradOpt gopt( &mesh );

  std::vector< ccdl::Mesh2dHessian > balls;

  int nx = mesh.GetSizeX();
  int ny = mesh.GetSizeY();
  double xlow = mesh.GetLowX();
  double ylow = mesh.GetLowY();
  double lx = mesh.GetHighX() - xlow;
  double ly = mesh.GetHighY() - ylow;
  double delx = mesh.GetX(1)-mesh.GetX(0);
  double dely = mesh.GetY(1)-mesh.GetY(0);


  std::vector< ccdl::Mesh2dValue > vals(nx*ny);
  for ( int i=0; i<nx; ++i )
    for ( int j=0; j<ny; ++j )
      vals[j+i*ny] = mesh.GetValue( mesh.GetX(i), mesh.GetY(j) );
	
  int const nxn = periodic ? nx : nx-1;
  int const nyn = periodic ? ny : ny-1;
  int const nx0 = periodic ? 0 : 1;
  int const ny0 = periodic ? 0 : 1;


  for ( int i=nx0; i<nxn; ++i )
    {
      if ( i % 10 == 0 )
	std::printf("Optimizing %i / %i\n",i*(nyn-ny0),(nxn-nx0)*(nyn-ny0));
      for ( int j=ny0; j<nyn; ++j )
	{
	  int i1 = mesh.GetIndex( i+1, nx );
	  int j1 = mesh.GetIndex( j+1, ny );

	  ccdl::Mesh2dValue & ll = vals[j +i *ny];
	  ccdl::Mesh2dValue & lr = vals[j +i1*ny];
	  ccdl::Mesh2dValue & ul = vals[j1+i *ny];
	  ccdl::Mesh2dValue & ur = vals[j1+i1*ny];

	  double l2r[2] = { delx, 0. };
	  double u2l[2] = { 0., dely };
	  double ul2lr[2] = { delx, -dely };
	  double ll2ur[2] = { delx,  dely };

	  double g1 = ( l2r[0]*ll.dfdx + l2r[1]*ll.dfdy ) * ( l2r[0]*lr.dfdx + l2r[1]*lr.dfdy );
	  double g2 = ( u2l[0]*ul.dfdx + u2l[1]*ul.dfdy ) * ( u2l[0]*ll.dfdx + u2l[1]*ll.dfdy );
	  double g3 = ( ll2ur[0]*ll.dfdx + ll2ur[1]*ll.dfdy ) * ( ll2ur[0]*ur.dfdx + ll2ur[1]*ur.dfdy );
	  double g4 = ( ul2lr[0]*ul.dfdx + ul2lr[1]*ul.dfdy ) * ( ul2lr[0]*lr.dfdx + ul2lr[1]*lr.dfdy );
	  double g5 = ( l2r[0]*ul.dfdx + l2r[1]*ul.dfdy ) * ( l2r[0]*ur.dfdx + l2r[1]*ur.dfdy );
	  double g6 = ( u2l[0]*ur.dfdx + u2l[1]*ur.dfdy ) * ( u2l[0]*lr.dfdx + u2l[1]*lr.dfdy );


	  // if ( ll.x > -2.05 and ll.x < -1.9 and 
	  //      ll.y > -0.6 and ll.y < -0.4 )
	  //   {
	  //     std::printf("KK  %18.10f %18.10f  %11.2e %11.2e %11.2e %11.2e\n",
	  // 		ll.x,
	  // 		ll.y,
	  // 		g1,g2,g3,g4 );
	  //     std::printf("ll  %18.10f %18.10f  %11.2e %11.2e\n",
	  // 		ll.x, ll.y, ll.dfdx, ll.dfdy );
	  //     std::printf("lr  %18.10f %18.10f  %11.2e %11.2e\n",
	  // 		lr.x, lr.y, lr.dfdx, lr.dfdy );
	  //     std::printf("ul  %18.10f %18.10f  %11.2e %11.2e\n",
	  // 		ul.x, ul.y, ul.dfdx, ul.dfdy );
	  //     std::printf("ur  %18.10f %18.10f  %11.2e %11.2e\n",
	  // 		ur.x, ur.y, ur.dfdx, ur.dfdy );
	  //     std::printf("\n");
	  //   }

	  int nsatisfied = (g1<=0) + (g2<=0) + (g3<=0) + (g4<=0) + (g5<=0) + (g6<=0);

	  if ( nsatisfied > 2 )
	    {
	      // std::printf("exr %18.10f %18.10f  %13.4e %13.4e\n",
	      // 		ll.x + 0.5 * delx,
	      //   		ll.y + 0.5 * dely,
	      //   		ll.dfdx, ll.dfdy );
	      // there is an extremum in this quadrant
	      //if ( (ll.dfdx < 0 and ll.dfdy < 0) or (  )
	      {

		std::vector<double> o(2);
		o[0] = 0.01 * delx;
		o[1] = 0.01 * dely;

		nlopt::opt nlo( nlopt::LN_COBYLA, 2 );
		nlo.set_min_objective( &NLOPT_GradMinimize, &gopt );
		nlo.set_xtol_rel( 1.e-12 );
		nlo.set_ftol_abs( 1.e-15 );
		nlo.set_maxeval( 5000 );
		nlo.set_initial_step( o );

		o[0] = ll.x + 0. * delx;
		o[1] = ll.y + 0. * dely;

		double chisq = 0.;
		//nlopt::result res = 
		nlo.optimize( o, chisq );

		if ( chisq < 0.0001 )
		  balls.push_back( mesh.GetHessian( o[0], o[1] ) );
	      }
	    }
	}
    }


  std::sort( balls.begin(), balls.end(), ccdl::sort_hessian_by_position );

  typedef std::vector< ccdl::Mesh2dHessian >::iterator iter;


  for ( iter p = balls.begin(); p != balls.end(); ++p )
    for ( iter q = p+1; q != balls.end(); )
      {
	double dx = std::abs(p->x - q->x);
	double dy = std::abs(p->y - q->y);
	if ( (dx < 0.25 * delx and dy < 0.25 * dely) )
	  q = balls.erase( q );
	else
	  ++q;
      }

  
  if ( ! periodic )
    for ( iter p = balls.begin(); p != balls.end(); )
      if ( std::abs( p->x - xlow ) / lx < 0.03 or
	   std::abs( p->y - ylow ) / lx < 0.03 or
	   std::abs( p->x - (xlow+lx) ) / lx < 0.03 or
	   std::abs( p->y - (ylow+ly) ) / ly < 0.03 )
	p = balls.erase( p );
      else
	++p;
  

  return balls;
}

















bool sort_by_f( ccdl::Mesh2dHessian const & lhs, ccdl::Mesh2dHessian const & rhs )
{
  return lhs.f < rhs.f;
}

void GetStationaryPts
( ccdl::Mesh2d & mesh, bool periodic, bool rezero,
  std::vector<ccdl::Mesh2dHessian> & minima,
  std::vector<ccdl::Mesh2dHessian> & maxima,
  std::vector<ccdl::Mesh2dHessian> & saddle )
{
  minima.resize(0);
  maxima.resize(0);
  saddle.resize(0);
  std::vector<ccdl::Mesh2dHessian> all( LocStationaryPts( mesh, periodic ) );
  typedef std::vector<ccdl::Mesh2dHessian>::iterator iter;
  for ( iter p = all.begin(), pend=all.end(); p!=pend; ++p )
    if ( p->eval[0] > 0 and p->eval[1] > 0 )
      minima.push_back( *p );
    else if ( p->eval[0] < 0 and p->eval[1] < 0 )
      maxima.push_back( *p );
    else
      saddle.push_back( *p );
  std::sort( minima.begin(), minima.end(), sort_by_f );
  std::sort( maxima.begin(), maxima.end(), sort_by_f );
  std::sort( saddle.begin(), saddle.end(), sort_by_f );
  if ( rezero )
    {
      mesh.Add( - minima[0].f );
      for ( iter p=minima.begin(), pend=minima.end(); p!=pend; ++p )
	*p = mesh.GetHessian( p->x, p->y );
      for ( iter p=maxima.begin(), pend=maxima.end(); p!=pend; ++p )
	*p = mesh.GetHessian( p->x, p->y );
      for ( iter p=saddle.begin(), pend=saddle.end(); p!=pend; ++p )
	*p = mesh.GetHessian( p->x, p->y );
    }
}












#define FF(a) std::setw(20) << std::setprecision(10) << std::fixed << (a)



std::string main_points( char const * meshname, bool periodic, bool rezero )
{
  std::ifstream cin;
  cin.open( meshname );
  if ( ! cin.good() )
    {
      std::cerr << meshname << " file not found!\n";
      std::abort();
    };

  ccdl::Mesh2d mesh( cin, periodic );

  typedef std::vector< ccdl::Mesh2dHessian > vec;
  typedef vec::iterator iter;


  vec minima, maxima, saddle;
  GetStationaryPts( mesh, periodic, rezero, minima, maxima, saddle );
  std::string fname( meshname );
  if ( rezero )
    {
      fname += ".0.dat";
      std::ofstream cout;
      cout.open( fname.c_str() );
      mesh.Write( cout );
    }

  {
    std::ofstream cout;
    cout.open("minpts.dat");
    for ( iter p = minima.begin(), pend = minima.end(); p!=pend; ++p )
      cout << FF(p->x) << FF(p->y) 
	   << " # " << FF(p->f) << FF(p->dfdx) << FF(p->dfdy) 
	   << FF(p->eval[0]) << FF(p->eval[1]) << "\n";
  }
  {
    std::ofstream cout;
    cout.open("maxpts.dat");
    for ( iter p = maxima.begin(), pend = maxima.end(); p!=pend; ++p )
      cout << FF(p->x) << FF(p->y) 
	   << " # "  << FF(p->f) << FF(p->dfdx) << FF(p->dfdy) 
	   << FF(p->eval[0]) << FF(p->eval[1]) << "\n";
  }
  {
    std::ofstream cout;
    cout.open("tspts.dat");
    for ( iter p = saddle.begin(), pend = saddle.end(); p!=pend; ++p )
      cout << FF(p->x) << FF(p->y) 
	   << " # "  << FF(p->f) << FF(p->dfdx) << FF(p->dfdy) 
	   << FF(p->eval[0]) << FF(p->eval[1]) << "\n";
  }

  {
    std::string gnuplot( fname );
    gnuplot += ".sh";
    std::ofstream cout;
    cout.open( gnuplot.c_str() );

    cout << "#!/bin/sh" << "\n";
    cout << "cat <<EOF | gnuplot" << "\n";
    cout << "# make contours" << "\n";
    cout << "set contour base" << "\n";
    cout << "set cntrparam level incremental 0, 1, 50" << "\n";
    cout << "unset surface" << "\n";
    cout << "set table 'contours.dat'" << "\n";
    cout << "splot '" << fname << "' using 1:2:3" << "\n";
    cout << "unset table" << "\n";
    cout << "reset" << "\n";
    cout << "unset key" << "\n";
    cout << "#set palette rgbformulae 33,13,10" << "\n";
    cout << "set palette defined ( 0 '#0034a9', 1 '#1c61ff', 2 '#3371ff', 3 '#009cff', 4 '#00d8ff', 5 '#00ffba', 6 '#00ff6c', 7 '#00ff0c', 8 '#aeff00', 9 '#fff600', 10 '#ffc000', 11 '#ff9600', 12 '#ff6c00', 13 '#ff3c00', 14 '#ff0000', 15 '#7f0000', 16 '#530000')" << "\n";

    cout << "" << "\n";
    cout << "# minpts" << "\n";
    cout << "set style line 2 lc rgb 'black' pt 7 ps 1.2" << "\n";
    cout << "# maxpts" << "\n";
    cout << "set style line 3 lc rgb 'black' pt 6 ps 1.2 lw 1.5" << "\n";
    cout << "# tspts" << "\n";
    cout << "set style line 4 lc rgb 'black' pt 2 ps 1.  lw 1.5" << "\n";
    cout << "# paths" << "\n";
    cout << "set style line 5 lc rgb 'red' pt 7 ps 0.7 lw 1.5" << "\n";
    cout << "" << "\n\n\n";
    cout << "set style rect fillcolor lt -3 fillstyle solid 1.0 noborder front\n";

    cout << "# minima labels\n";
    int i=1;
    for ( iter p=minima.begin(),pend=minima.end(); p!=pend; ++p, ++i )
      {
	cout << "set label " << i << " '" << i << "' at "
	     << FF( p->x ) << "," << FF( p->y ) << " center tc rgb \"blue\" font \",10\" front" << "\n";
	cout << "set object " << i << " rect center " 
	     << FF(p->x ) << "," << FF(p->y) << " size char strlen('" << i << "'),char 1 front" << "\n";
      }

    cout << "# saddle pt labels\n";
    int j=1;
    for ( iter p=saddle.begin(),pend=saddle.end(); p!=pend; ++p, ++i, ++j )
      {
	cout << "set label " << i << " '" << j << "' at "
	     << FF( p->x ) << "," << FF( p->y ) << "  center tc rgb \"red\" font \",10\" front" << "\n";
	cout << "set object " << i << " rect center " 
	     << FF(p->x) << "," << FF(p->y) << " size char strlen('" << j << "'),char 1 front" << "\n";
      }


    cout << "" << "\n\n\n";
    cout << "set terminal postscript eps enhanced solid color 'Helvetica' 12" << "\n";
    cout << "set output '" << fname << ".eps'" << "\n";
    cout << "" << "\n";

    std::stringstream cblo;
    std::stringstream cbhi;
    if ( minima.size() ) cblo << FF( minima[0].f );
    if ( maxima.size() ) cbhi << FF( maxima.back().f );
    else if ( saddle.size() ) cbhi << FF( saddle.back().f );

    cout << "set cbrange[" << cblo.str() << ":" << cbhi.str() << "]" << "\n";
    cout << "set cblabel 'kcal/mol'\n";
    cout << "set ylabel 'Y-axis label'\n";
    cout << "set xlabel 'X-axis label'\n";
    cout << "p '" << fname << "' using 1:2:3 with image,\\" << "\n";
    cout << "  'contours.dat' w l lt -1 lw 0.5,\\" << "\n";
    cout << "  'minpts.dat'   w p ls 2,\\" << "\n";
    cout << "  'maxpts.dat'   w p ls 3,\\" << "\n";
    cout << "  'tspts.dat'    w p ls 4,\\" << "\n";
    cout << "#  'min1_ts1.dat' w l ls 5" << "\n";
    //cout << "pause -1\n";
    cout << "EOF" << "\n";
  }

  return fname;
}








int main( int argc, char ** argv)
{

  cli_options cli( read_options( argc, argv ) );

  if ( (! cli.connections.size()) or cli.rezero )
    cli.meshfile = main_points( cli.meshfile.c_str(), cli.perx or cli.pery, cli.rezero );
  
  if ( cli.connections.size() > 0 )
    {
      std::ifstream cin;
      cin.open( cli.meshfile.c_str() );
      if ( ! cin.good() )
	{
	  std::cerr << "could not open " << cli.meshfile << "\n";
	  std::abort();
	}

      ccdl::Mesh2d mesh( cin, cli.perx or cli.pery );

      int icon = 0;
      int ncon = cli.connections.size();
      double arcsum = 0.;
      for ( std::vector< connection >::iterator 
	      p=cli.connections.begin(), pend=cli.connections.end();
	    p != pend; ++p, ++icon )
	{
	  std::cout << "Performing NEB " 
		    << icon+1 << " / " << ncon
		    << " : " << p->GetName() << "\n";
	  
	  if ( p->method == 1 )
	    p->curve.reset( new ccdl::ParametricLegendre
			    (  NudgeElasticBand_FitMethod1
			       ( mesh, p->t0.x, p->t0.y, p->t1.x, p->t1.y ) ) );
	  else if ( p->method == 2 )
	    p->curve.reset( new ccdl::ParametricLegendre
			    (  NudgeElasticBand_FitMethod2
			       ( mesh, p->t0.x, p->t0.y, p->t1.x, p->t1.y ) ) );
	  

	  ccdl::ParametricLegendre & neb = *(p->curve);
	  std::ofstream xyfile,tfile;
	  std::stringstream str;
	  str << p->GetName() << ".2d.dat";
	  xyfile.open( str.str().c_str() );
	  for ( int i=0; i<501; ++i )
	    {
	      double t = i / 500.;
	      double x,y;
	      neb.GetXY( t, x, y );
	      ccdl::Mesh2dValue val( mesh.GetValue( x,y ) );
	      xyfile << FF(x) << FF(y) 
		     << " # " << FF(val.f) <<  "\n";
	    }
	  str.str("");
	  str.clear();
	  str << p->GetName() << ".1d.dat";

	  if ( p->prev < 0 ) arcsum = 0.;

	  tfile.open( str.str().c_str() );
	  for ( int i=0; i<501; ++i )
	    {
	      double t = i / 500.;
	      double x,y;
	      neb.GetXY( t, x, y );
	      double arclen = neb.GetArcLength( 0, t );
	      ccdl::Mesh2dValue val( mesh.GetValue( x,y ) );
	      tfile << FF(arclen+arcsum) << FF(val.f) 
		    << " # " << FF(val.x) << FF(val.y) << "\n";
	    }
	  arcsum += neb.GetArcLength( 0., 1. );
	};




      std::vector<std::string> colors;
      colors.push_back( std::string("black") );
      colors.push_back( std::string("red") );
      colors.push_back( std::string("green") );
      colors.push_back( std::string("blue") );
      colors.push_back( std::string("yellow") );
      colors.push_back( std::string("brown") );
      colors.push_back( std::string("grey") );
      colors.push_back( std::string("#4B0082") );
      colors.push_back( std::string("cyan") );
      colors.push_back( std::string("magenta") );
      colors.push_back( std::string("orange") );



      {
	std::string gnuplot( cli.meshfile );
	gnuplot += ".sh";
	std::ofstream cout;
	cout.open( gnuplot.c_str() );

	cout << "#!/bin/sh" << "\n";
	cout << "cat <<EOF | gnuplot" << "\n";
	cout << "# make contours" << "\n";
	cout << "set contour base" << "\n";
	cout << "set cntrparam level incremental 0, 1, 50" << "\n";
	cout << "unset surface" << "\n";
	cout << "set table 'contours.dat'" << "\n";
	cout << "splot '" << cli.meshfile << "' using 1:2:3" << "\n";
	cout << "unset table" << "\n";
	cout << "reset" << "\n";
	cout << "unset key" << "\n";
	cout << "#set palette rgbformulae 33,13,10" << "\n";
	cout << "set palette defined ( 0 '#0034a9', 1 '#1c61ff', 2 '#3371ff', 3 '#009cff', 4 '#00d8ff', 5 '#00ffba', 6 '#00ff6c', 7 '#00ff0c', 8 '#aeff00', 9 '#fff600', 10 '#ffc000', 11 '#ff9600', 12 '#ff6c00', 13 '#ff3c00', 14 '#ff0000', 15 '#7f0000', 16 '#530000')" << "\n";

	cout << "" << "\n";
	// cout << "# minpts" << "\n";
	// cout << "set style line 2 lc rgb 'black' pt 7 ps 1.2" << "\n";
	// cout << "# maxpts" << "\n";
	// cout << "set style line 3 lc rgb 'black' pt 6 ps 1.2 lw 3" << "\n";
	// cout << "# tspts" << "\n";
	// cout << "set style line 4 lc rgb 'black' pt 2 ps 1.  lw 3" << "\n";
	// cout << "# paths" << "\n";
	// cout << "set style line 5 lc rgb 'red' pt 7 ps 0.7 lw 3" << "\n";

	cout << "# minpts" << "\n";
	cout << "set style line 2 lc rgb 'black' pt 7 ps 1.2" << "\n";
	cout << "set style line 12 lc rgb 'white' pt 7 ps 0.7" << "\n";
	cout << "# maxpts" << "\n";
	cout << "set style line 3 lc rgb 'black' pt 12 ps 1. lw 3" << "\n";
	cout << "set style line 13 lc rgb 'white' pt 12 ps 0.75 lw 1" << "\n";
	cout << "# tspts" << "\n";
	cout << "set style line 4 lc rgb 'black' pt 2 ps 1.  lw 6" << "\n";
	cout << "set style line 14 lc rgb 'white' pt 2 ps 0.75  lw 1.5" << "\n";
	cout << "# paths" << "\n";
	cout << "set style line 100 lc rgb 'red' pt 7 ps 0.7 lw 2.5" << "\n";
	cout << "set style line 110 lc rgb 'white' pt 7 ps 0.7 lw 5.5" << "\n";

	cout << "" << "\n\n\n";



	cout << "" << "\n\n\n";
	cout << "set terminal postscript eps enhanced solid color 'Helvetica' 14" << "\n";
	cout << "set output '" << cli.meshfile << ".eps'" << "\n";
	cout << "" << "\n";

	std::stringstream cblo;
	std::stringstream cbhi;
	if ( cli.minima.size() ) cblo << FF( mesh.GetValue( cli.minima[0].x, cli.minima[0].y ).f );
	if ( cli.maxima.size() ) cbhi << FF( mesh.GetValue( cli.maxima.back().x, cli.maxima.back().y ).f );
	else if ( cli.saddle.size() ) cbhi << FF( mesh.GetValue( cli.saddle.back().x, cli.saddle.back().y ).f );

	cout << "set cbrange[" << cblo.str() << ":" << cbhi.str() << "]" << "\n";
	cout << "set cblabel 'kcal/mol'\n";
	cout << "set ylabel 'Y-axis label'\n";
	cout << "set xlabel 'X-axis label'\n";
	cout << "p '" << cli.meshfile << "' using 1:2:3 with image,\\" << "\n";
	cout << "  'contours.dat' w l lt -1 lw 0.5,\\" << "\n";
	int icolor = 0;
	std::stringstream xmgrace;
	for ( std::vector< connection >::iterator 
		p=cli.connections.begin(), pend=cli.connections.end();
	      p != pend; ++p, ++icon, ++icolor )
	  {
	    if ( p->prev < 0 ) 
	      { 
		icolor = 0;
		if ( xmgrace.str().size() > 0 )
		  std::cout << xmgrace.str() << "\n";
		xmgrace.str("");
		xmgrace.clear();
		xmgrace << "xmgrace " << p->GetName() << ".1d.dat";
	      }
	    else 
	      { 
		icolor = icolor % (int)(colors.size());
		xmgrace << " " << p->GetName() << ".1d.dat";
	      }
	    cout << "  '" << p->GetName() << ".2d.dat' w l ls 110,\\" << "\n";
	    cout << "  '" << p->GetName() << ".2d.dat' w l ls 100 lc rgb '" << colors[icolor] << "',\\" << "\n";
	  }
	std::cout << xmgrace.str() << "\n";
	cout << "  'minpts.dat'   w p ls 2,\\" << "\n";
	cout << "  'minpts.dat'   w p ls 12,\\" << "\n";
	cout << "  'maxpts.dat'   w p ls 3,\\" << "\n";
	cout << "  'maxpts.dat'   w p ls 13,\\" << "\n";
	cout << "  'tspts.dat'    w p ls 4,\\" << "\n";
	cout << "  'tspts.dat'    w p ls 14,\\" << "\n";
	cout << "#\n";
	cout << "EOF" << "\n";
      }


    }
  

  

  /*
  bool periodic = false;
  ccdl::Mesh2d mesh( std::cin, periodic );


  typedef std::vector< ccdl::Mesh2dHessian > vec;
  typedef vec::iterator iter;

  vec minima, maxima, saddle;
  GetStationaryPts( mesh, minima, maxima, saddle );
  {
    std::ofstream cout;
    cout.open("minpts.dat");
    for ( iter p = minima.begin(), pend = minima.end(); p!=pend; ++p )
      cout << FF(p->x) << FF(p->y) 
	   << " # " << FF(p->f) << FF(p->dfdx) << FF(p->dfdy) 
	   << FF(p->eval[0]) << FF(p->eval[1]) << "\n";
  }
  {
    std::ofstream cout;
    cout.open("maxpts.dat");
    for ( iter p = maxima.begin(), pend = maxima.end(); p!=pend; ++p )
      cout << FF(p->x) << FF(p->y) 
	   << " # "  << FF(p->f) << FF(p->dfdx) << FF(p->dfdy) 
	   << FF(p->eval[0]) << FF(p->eval[1]) << "\n";
  }
  {
    std::ofstream cout;
    cout.open("tspts.dat");
    for ( iter p = saddle.begin(), pend = saddle.end(); p!=pend; ++p )
      cout << FF(p->x) << FF(p->y) 
	   << " # "  << FF(p->f) << FF(p->dfdx) << FF(p->dfdy) 
	   << FF(p->eval[0]) << FF(p->eval[1]) << "\n";
  }

  */


  return 0;


  /*
  std::vector< ccdl::Mesh2dValue > minima = mesh.GetMinima();
  {
    std::ofstream cout;
    cout.open("minpts.dat");
    for ( std::size_t i=0; i<minima.size(); ++i )
    {
      ccdl::Mesh2dValue v = minima[i];
      cout << FF(v.x) << FF(v.y) << "\n";
    }
  }

  std::vector< ccdl::Mesh2dValue > maxima = mesh.GetMaxima();
  {
   std::ofstream cout;
    cout.open("maxpts.dat");
    for ( std::size_t i=0; i<maxima.size(); ++i )
    {
      ccdl::Mesh2dValue v = maxima[i];
      cout << FF(v.x) << FF(v.y) << "\n";
    }
  }
  */



  /*

  // double SpringWt          = 1.;
  // double EnergyWt          = 0.05; 
  // double CurvatureWt       = SpringWt * 0.01;
  // double TransitionStateWt = EnergyWt * 10000.;
  // double Temperature       = 1.;

  double SpringWt          = 1.;
  double EnergyWt          = 0.01; 
  double CurvatureWt       = 0.01;
  double TransitionStateWt = 0.01;
  double Temperature       = 0.01;




  double x0 = minima[4].x;
  double y0 = minima[4].y;
  double x1 = minima[7].x;
  double y1 = minima[7].y;

  ccdl::ParametricLegendre curve( 3, 3, 100 );
  curve.SetEndptValues( x0,y0, x1,y1 );

  int cnt = 0;
  //int orders[3] = { 5, 10, 15 };
  int orders[6] = { 3, 5, 8, 12, 15, 18 };
  for ( int iorder=0; iorder < 5; iorder++, ++cnt )
    {
      double TW = TransitionStateWt;
      double EW = EnergyWt;
      double SW = SpringWt;

      //if ( iorder < 3 ) { EW = 0.05; SW = 50.; TW = 250.; }
      //if ( iorder < 1 ) { EW = 0.1; SW = 100.; TW = 500.; }



      int order = orders[iorder];
      ccdl::ParametricLegendre new_curve( order, order, 100 );
      new_curve.SetEndptValues( x0,y0, x1,y1 );
      int npt = new_curve.GetNpts();
      std::vector<double> xt(npt,0.), yt(npt,0.); 
      if ( cnt == 0 )
	{
	  double mx = (x1-x0);
	  double bx = x0;
	  double my = (y1-y0);
	  double by = y0;
	  for ( int i=0; i<npt; ++i )
	    {
	      double t = new_curve.GetQuadPt(i);
	      xt[i] = mx * t + bx;
	      yt[i] = my * t + by;
	      //std::printf("%20.10f %20.10f\n",xt[i],yt[i]);
	    }
	}
      else
	{
	  for ( int i=0; i<npt; ++i )
	    {
	      double t = new_curve.GetQuadPt(i);
	      curve.GetXY( t, xt[i], yt[i] );
	    }
	}
      new_curve.SetXY( xt.data(), yt.data() );
      NudgeElasticBand( mesh, new_curve, SW, CurvatureWt, EW, TW, Temperature );
      curve = new_curve;
    };



  // for ( int i=0; i<curve.GetNpts(); ++i )
  //   {
  //     double t = curve.GetQuadPt(i);
  //     double x,y;
  //     curve.GetXY( t, x, y );
  //     std::printf("GREP %20.10f %20.10f\n",x,y);
  //   }
  
  */













  /*
  ccdl::ParametricLegendre neb1
    ( NudgeElasticBand_FitMethod1
      ( mesh,
	minima[0],minima[7] ) );
  
  {
    std::ofstream xyfile,tfile;
    xyfile.open( "neb1_0_7_xy.dat" );
    for ( int i=0; i<501; ++i )
      {
	double t = i / 500.;
	double x,y;
	neb1.GetXY( t, x, y );
	ccdl::Mesh2dValue val( mesh.GetValue( x,y ) );
	xyfile << FF(x) << FF(y) << " # " << FF(val.f) <<  "\n";
      }
    tfile.open( "neb1_0_7_t.dat" );
    for ( int i=0; i<501; ++i )
      {
	double t = i / 500.;
	double x,y;
	neb1.GetXY( t, x, y );
	ccdl::Mesh2dValue val( mesh.GetValue( x,y ) );
	tfile << FF(t) << FF(val.f) << " # " << FF(val.x) << FF(val.y) << "\n";
      }

  };
  */
  
  return 0;
}
