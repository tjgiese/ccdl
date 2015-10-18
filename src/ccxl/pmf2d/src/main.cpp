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
  endpt() : x(0), y(0), i(0), minimum(false), arclen(0.),value(0.) {}
  double x,y;
  int i;
  std::string label;
  bool minimum;
  double arclen;
  double value;
  std::string GetName() const;
  bool operator== ( endpt const rhs ) const { return i == rhs.i; }
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




struct RefinementPt
{
  RefinementPt() : x(0), y(0) {};
  RefinementPt( double x, double y ) : x(x), y(y) {};
  double x,y;
};

bool operator< ( RefinementPt const & lhs, RefinementPt const & rhs )
{
  return (lhs.x == rhs.x) ? lhs.y < rhs.y : lhs.x < rhs.x;
}

bool operator== ( RefinementPt const & lhs, RefinementPt const & rhs )
{
  return (lhs.x == rhs.x) and (lhs.y == rhs.y);
}


struct RefinementGrid
{
  RefinementGrid( double xlo, double ylo, double dx, double dy, int n ) 
    : xlo(xlo), ylo(ylo), dx(dx),dy(dy),n(n) {};

  double xlo,ylo;
  double dx,dy;
  int n;
};


struct BoxCrd
{
  BoxCrd() : xlo(0),ylo(0),xhi(0),yhi(0) {}
  BoxCrd( double xl, double yl, double xh, double yh )
    : xlo(xl),ylo(yl),xhi(xh),yhi(yh) {}

  bool IsInside( double x, double y ) const;

  double xlo,ylo,xhi,yhi;
};

bool BoxCrd::IsInside( double x, double y ) const
{
  return (x>=xlo and x<=xhi and y>=ylo and y<=yhi);
}

bool IsInside( double x, double y, std::vector<BoxCrd> const & blanks )
{
  bool inside = false;
  for ( std::vector<BoxCrd>::const_iterator 
	  p=blanks.begin(),pend=blanks.end(); 
	p!=pend; ++p )
    {
      inside = p->IsInside( x,y );
      if ( inside ) break;
    };
  return inside;
}


struct cli_options
{
  cli_options();
  std::vector< endpt > minima,maxima,saddle;
  std::vector< connection > connections;
  std::tr1::shared_ptr< RefinementGrid > refine;
  std::vector<BoxCrd> blanks;
  std::vector<RefinementPt> extra_searches;
  std::string meshfile;
  bool rezero;
  bool perx,pery;
  bool delunused;
  std::string range_xlo,range_xhi,range_ylo,range_yhi,range_clo,range_chi;
};


cli_options::cli_options()
  : rezero(false),perx(false), pery(false),delunused(false)
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
  std::printf("  --refine x=XLO,y=YLO,dx=DX,dy=DY,n=N\n");
  std::printf("           \tReturns a set of x,y values whose locations are within\n");
  std::printf("           \t|x(t)-(XLO+i*DX)| <= 2*N*DX and |y(t)-(YLO+j*DY)| <= 2*N*DY\n");
  std::printf("           \tof the paths define by the --con commands\n");
  std::printf("           \tFor example, if you have a coarse-grained guess at the PMF\n");
  std::printf("           \tand you want create new windows in a box from 1.4 to 4.0,\n");
  std::printf("           \tseparated by 0.1 Ang, then\n");
  std::printf("           \t--con min=1,ts=1,min=2 --refine x=1.4,y=1.4,dx=0.1,dy=0.1,n=1\n");
  std::printf("  --blank xlo=XLO,ylo=YLO,xhi=XHI,yhi=YHI\n");
  std::printf("           \tCover the selected area with a white box\n");
  std::printf("  --opt x=X,y=Y\n");
  std::printf("           \tPerform a stationary point search around X,Y\n");
  std::printf("  -d\n");
  std::printf("  --delete\n");
  std::printf("           \tExcludes unused stationary points when using --con\n");
  std::printf("  --range xlo=XLO,ylo=YLO,xhi=XHI,yhi=YHI,clo=CLO,chi=CHI\n");
  std::printf("           \tSets the contour axis range; unused values use gnuplot defaults\n");


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

  char const * REFINE_opts[] = {
#define REFINE_X 0
    "x",
#define REFINE_Y 1
    "y",
#define REFINE_DX 2
    "dx",
#define REFINE_DY 3
    "dy",
#define REFINE_N 4
    "n",
    NULL };

  char const * BLANK_opts[] = {
#define BLANK_XLO 0
    "xlo",
#define BLANK_YLO 1
    "ylo",
#define BLANK_XHI 2
    "xhi",
#define BLANK_YHI 3
    "yhi",
    NULL };

  char const * RANGE_opts[] = {
#define RANGE_XLO 0
    "xlo",
#define RANGE_YLO 1
    "ylo",
#define RANGE_XHI 2
    "xhi",
#define RANGE_YHI 3
    "yhi",
#define RANGE_CLO 4
    "clo",
#define RANGE_CHI 5
    "chi",
    NULL };


  char const * OPT_opts[] = {
#define OPT_X 0
    "x",
#define OPT_Y 1
    "y",
    NULL };

  static struct option long_options[] =
    {
      { "help",      no_argument,       NULL, 'h'   },
      { "zero",      no_argument,       NULL, 'z'   },
      { "perx",      no_argument,       NULL, 'x'   },
      { "pery",      no_argument,       NULL, 'y'   },
      { "delete",    no_argument,       NULL, 'd' },
#define CON     0100
      { "con",       required_argument, NULL, CON   },
#define REFINE  0200
      { "refine",    required_argument, NULL, REFINE },
#define BLANK   0300
      { "blank",     required_argument, NULL, BLANK },
#define OPT     0400
      { "opt",       required_argument, NULL, OPT },
#define RANGE   0500
      { "range",     required_argument, NULL, RANGE },

      {NULL,0,NULL,0}
    };

  int opt = 0;
  int long_index = 0;
  char * subopts, * value;
  while ( (opt = getopt_long
          ( argc, argv, "hxyd", 
            long_options, &long_index )) != -1 )
    {
      
      switch (opt)
        {
	  
        case 'h':     { print_usage(); break; }
	case 'z':     { cli.rezero = true; break; }
        case 'x':     { cli.perx = true; break; }
        case 'y':     { cli.pery = true; break; }
        case 'd':     { cli.delunused = true; break; }
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
		int iprev = -1;
		for ( std::size_t i=1; i<ts.size(); ++i )
		  {
		    connection c( ts[i-1],ts[i] );
		    c.method = method;
		    c.prev = iprev;
		    iprev = cli.connections.size();
		    cli.connections.push_back( c );
		  };
	      }
	    break;
	  }

	case REFINE:
	  {
	    if ( cli.refine )
	      {
		std::printf("%s: Only one instance of --refine is possible\n",argv[0]);
		std::exit(EXIT_FAILURE);
	      }
	    bool hasx = false;
	    bool hasy = false;
	    double x=0;
	    double y=0;
	    double dx=-1;
	    double dy=-1;
	    int n = 1;
            subopts = optarg;
            while (*subopts != '\0')
              switch (getsubopt(&subopts, const_cast<char**>(REFINE_opts), &value))
                {
                case REFINE_X:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --refine x= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    x = std::atof(value);
		    hasx=true;
		    break;
		  }
                case REFINE_Y:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --refine y= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    y = std::atof(value);
		    hasy=true;
		    break;
		  }
                case REFINE_DX:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --refine dx= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    dx = std::atof(value);
		    break;
		  }
                case REFINE_DY:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --refine dy= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    dy = std::atof(value);
		    break;
		  }
                case REFINE_N:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --refine n= expects an int; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    n = std::atoi(value);
		    break;
		  }
		default:
		  {
		    std::printf("%s: Unknown subption '%s' given to --refine\n",argv[0],value);
		    std::exit(EXIT_FAILURE);
		  }
		}
	    {
	      bool error = false;
	      if ( dx <= 0 )
		{
		  std::printf("%s: invalid value of --refine dx=%.3e\n",argv[0],dx);
		  error = true;
		}
	      if ( dy <= 0 )
		{
		  std::printf("%s: invalid value of --refine dy=%.3e\n",argv[0],dy);
		  error = true;
		}
	      if ( ! hasx )
		{
		  std::printf("%s: --refine used without specifying x\n",argv[0]);
		  error = true;
		}
	      if ( ! hasy )
		{
		  std::printf("%s: --refine used without specifying y\n",argv[0]);
		  error = true;
		}
	      if ( error ) std::exit(EXIT_FAILURE);

	      cli.refine.reset( new RefinementGrid( x,y,dx,dy,n ) );
	    }
	    break;
	  }


	case BLANK:
	  {
	    double xlo=-1.e+10;
	    double ylo=-1.e+10;
	    double xhi=1.e+10;
	    double yhi=1.e+10;
            subopts = optarg;
            while (*subopts != '\0')
              switch (getsubopt(&subopts, const_cast<char**>(BLANK_opts), &value))
                {
                case BLANK_XLO:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --blank xlo= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    xlo = std::atof(value);
		    break;
		  }
                case BLANK_YLO:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --blank ylo= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    ylo = std::atof(value);
		    break;
		  }
                case BLANK_XHI:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --blank xhi= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    xhi = std::atof(value);
		    break;
		  }
                case BLANK_YHI:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --blank yhi= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    yhi = std::atof(value);
		    break;
		  }
		default:
		  {
		    std::printf("%s: Unknown subption '%s' given to --blank\n",argv[0],value);
		    std::exit(EXIT_FAILURE);
		  }
		}

	    cli.blanks.push_back( BoxCrd( xlo,ylo,xhi,yhi ) );

	    break;
	  }






	case RANGE:
	  {
            subopts = optarg;
            while (*subopts != '\0')
              switch (getsubopt(&subopts, const_cast<char**>(RANGE_opts), &value))
                {
                case RANGE_XLO:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --range xlo= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    cli.range_xlo = value;
		    break;
		  }
                case RANGE_YLO:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --range ylo= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    cli.range_ylo = value;
		    break;
		  }
                case RANGE_XHI:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --range xhi= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    cli.range_xhi = value;
		    break;
		  }
                case RANGE_YHI:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --range yhi= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    cli.range_yhi = value;
		    break;
		  }
                case RANGE_CLO:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --range clo= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    cli.range_clo = value;
		    break;
		  }
                case RANGE_CHI:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --range chi= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    cli.range_chi = value;
		    break;
		  }
		default:
		  {
		    std::printf("%s: Unknown subption '%s' given to --range\n",argv[0],value);
		    std::exit(EXIT_FAILURE);
		  }
		}


	    break;
	  }










	case OPT:
	  {
	    bool hasx = false;
	    bool hasy = false;
	    double x=0;
	    double y=0;
            subopts = optarg;
            while (*subopts != '\0')
              switch (getsubopt(&subopts, const_cast<char**>(OPT_opts), &value))
                {
                case OPT_X:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --opt x= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    x = std::atof(value);
		    hasx=true;
		    break;
		  }
                case REFINE_Y:
		  {
		    if ( value == NULL ) 
		      { std::printf("%s: --opt y= expects a float; use -h for usage\n",argv[0]); 
			std::exit(EXIT_FAILURE); }
		    y = std::atof(value);
		    hasy=true;
		    break;
		  }
		default:
		  {
		    std::printf("%s: Unknown subption '%s' given to --opt\n",argv[0],value);
		    std::exit(EXIT_FAILURE);
		  }
		}

	    if ( !hasx or !hasy )
	      {
		std::printf("%s: --opt expects *both* x and y suboptions\n",argv[0]);
		std::exit(EXIT_FAILURE);
	      }

	    cli.extra_searches.push_back( RefinementPt( x, y ) );

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


  int imin = 0;
  for ( std::vector<connection>::iterator 
	  p = cli.connections.begin(), pend=cli.connections.end();
	p != pend; ++p )
    {
      if ( p->t0.minimum )
	{
	  int i = p->t0.i - 1;
	  //std::printf("%s %i %i\n",p->t0.GetName().c_str(),i,(int)(cli.minima.size()));
	  if ( cli.minima[i].label.size() == 0 )
	    {
	      ++imin;
	      std::stringstream str;
	      str << imin;
	      cli.minima[i].label = str.str();
	    }
	}
      if ( p->t1.minimum )
	{
	  int i = p->t1.i - 1;
	  if ( cli.minima[i].label.size() == 0 )
	    {
	      ++imin;
	      std::stringstream str;
	      str << imin;
	      cli.minima[i].label = str.str();
	    }
	}
    };
  for ( std::vector<connection>::iterator 
	  p = cli.connections.begin(), pend=cli.connections.end();
	p != pend; ++p )
    {
      if ( p->t0.minimum )
	p->t0.label = cli.minima[ p->t0.i - 1 ].label;
      if ( p->t1.minimum )
	p->t1.label = cli.minima[ p->t1.i - 1 ].label;
    }

  endpt prevmin;
  prevmin.label = "?";
  for ( std::vector<connection>::iterator 
	  p = cli.connections.begin(), pend=cli.connections.end();
	p != pend; ++p )
    {
      endpt nextmin;
      nextmin.label = "?";
      
      if ( p->t0.minimum and ! p->t1.minimum )
	prevmin = p->t0;

      if ( ! p->t0.minimum )
	{
	  if ( p->t1.minimum )
	    nextmin = p->t1;
	  int i = p->t0.i - 1;
	  if ( cli.saddle[i].label.size() == 0 )
	    {
	      std::stringstream str;
	      str << prevmin.label << "-" << nextmin.label;
	      cli.saddle[i].label = str.str();
	    }
	}
      else if ( ! p->t1.minimum )
	{
	  if ( p != pend )
	    {
	      std::vector<connection>::iterator q = p+1;
	      if ( q != pend )
		{
		  if ( q->t0.minimum )
		    nextmin = q->t0;
		  else if ( q->t0.i == p->t1.i and q->t1.minimum )
		    nextmin = q->t1;
		};
	    };
	  int i = p->t1.i - 1;
	  if ( cli.saddle[i].label.size() == 0 )
	    {
	      std::stringstream str;
	      str << prevmin.label << "-" << nextmin.label;
	      cli.saddle[i].label = str.str();
	    }
	}
      prevmin = nextmin;
    }

  for ( std::vector<connection>::iterator 
	  p = cli.connections.begin(), pend=cli.connections.end();
	p != pend; ++p )
    {
      if ( ! p->t0.minimum )
	p->t0.label = cli.saddle[ p->t0.i - 1 ].label;
      if ( ! p->t1.minimum )
	p->t1.label = cli.saddle[ p->t1.i - 1 ].label;
    }

  for ( std::vector<connection>::iterator 
	  p = cli.connections.begin(), pend=cli.connections.end();
	p != pend; ++p )
    {
      std::printf("%2i %8s %8s %8s %8s\n",
		  p->prev,
		  p->t0.GetName().c_str(),
		  p->t1.GetName().c_str(),
		  p->t0.label.c_str(),
		  p->t1.label.c_str() );

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










std::vector< ccdl::Mesh2dHessian > LocStationaryPts
( ccdl::Mesh2d & mesh, bool periodic, std::vector<BoxCrd> const & blanks, std::vector<RefinementPt> const & extra_searches )
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

	  if ( IsInside( ll.x, ll.y, blanks ) or IsInside( ur.x, ur.y, blanks ) )
	    continue;

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

	  if ( nsatisfied > 1 )
	    {
	      // std::printf("exr %18.10f %18.10f  %13.4e %13.4e\n",
	      // 		ll.x + 0.5 * delx,
	      //   		ll.y + 0.5 * dely,
	      //   		ll.dfdx, ll.dfdy );
	      // there is an extremum in this quadrant
	      //if ( (ll.dfdx < 0 and ll.dfdy < 0) or (  )



	      {

		std::vector<double> o(2);
		o[0] = 0.005 * delx;
		o[1] = 0.005 * dely;

		nlopt::opt nlo( nlopt::LN_COBYLA, 2 );
		nlo.set_min_objective( &NLOPT_GradMinimize, &gopt );
		nlo.set_xtol_rel( 1.e-12 );
		nlo.set_ftol_abs( 1.e-15 );
		nlo.set_maxeval( 5000 );
		nlo.set_initial_step( o );

		o[0] = ll.x + 0.25 * delx;
		o[1] = ll.y + 0.25 * dely;

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
  

  for ( std::vector<RefinementPt>::const_iterator
	  p = extra_searches.begin(), pend = extra_searches.end();
	p != pend; ++p )
    {

      std::vector<double> o(2);
      o[0] = 0.005;
      o[1] = 0.005;

      nlopt::opt nlo( nlopt::LN_COBYLA, 2 );
      nlo.set_min_objective( &NLOPT_GradMinimize, &gopt );
      nlo.set_xtol_rel( 1.e-12 );
      nlo.set_ftol_abs( 1.e-15 );
      nlo.set_maxeval( 5000 );
      nlo.set_initial_step( o );

      o[0] = p->x;
      o[1] = p->y;

      double chisq = 0.;
      nlo.optimize( o, chisq );

      if ( chisq < 0.0001 )
	{
	  balls.push_back( mesh.GetHessian( o[0], o[1] ) );
	  std::cout <<"(Extra search point located @ "
		    << std::fixed << std::setprecision(3) << std::setw(12) << o[0]
		    << std::fixed << std::setprecision(3) << std::setw(12) << o[1]
		    << " )\n";
	}
      else
	{
	  std::cout <<"(failed to location search point @ "
		    << std::fixed << std::setprecision(3) << std::setw(12) << p->x
		    << std::fixed << std::setprecision(3) << std::setw(12) << p->y
		    << " )\n";
	}
    }

  std::sort( balls.begin(), balls.end(), ccdl::sort_hessian_by_position );


  return balls;
}

















bool sort_by_f( ccdl::Mesh2dHessian const & lhs, ccdl::Mesh2dHessian const & rhs )
{
  return lhs.f < rhs.f;
}

void GetStationaryPts
( ccdl::Mesh2d & mesh, bool periodic, std::vector<BoxCrd> const & blanks,
  std::vector<RefinementPt> const & extra_searches,
  bool rezero,
  std::vector<ccdl::Mesh2dHessian> & minima,
  std::vector<ccdl::Mesh2dHessian> & maxima,
  std::vector<ccdl::Mesh2dHessian> & saddle )
{
  minima.resize(0);
  maxima.resize(0);
  saddle.resize(0);
  std::vector<ccdl::Mesh2dHessian> all( LocStationaryPts( mesh, periodic, blanks, extra_searches ) );
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



std::string main_points( char const * meshname, 
			 bool periodic, 
			 std::vector<BoxCrd> const & blanks,
			 std::vector<RefinementPt> const & extra_searches,
			 bool rezero )
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
  GetStationaryPts( mesh, periodic, blanks, extra_searches, rezero, minima, maxima, saddle );
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
    cout << "set cntrparam level incremental -50, 2.5, 50" << "\n";
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

    cout << "set cblabel 'E (kcal/mol)'\n";
    cout << "set ylabel 'Y-axis label'\n";
    cout << "set xlabel 'X-axis label'\n";
    cout << "p '" << fname << "' using 1:2:3 with image,\\" << "\n";
    cout << "  'contours.dat' w l lt -1 lw 0.5,\\" << "\n";
    cout << "  'minpts.dat'   w p ls 2,\\" << "\n";
    cout << "  'maxpts.dat'   w p ls 3,\\" << "\n";
    cout << "  'tspts.dat'    w p ls 4" << "\n";
    //cout << "#  'min1_ts1.dat' w l ls 5" << "\n";
    //cout << "pause -1\n";
    cout << "EOF" << "\n";

    std::cout << "Now  run: /bin/sh \"" << gnuplot << "\"; evince \"" << fname << ".eps\"" << " &\n";
    std::cout << "Then run: pmf2d --con min=1,ts=1,min=2 \"" << fname << "\"" << "\n";
    std::cout << "(or something analogous)\n";

  }


  return fname;
}




void GetValidPosition( ccdl::Mesh2d const & mesh, endpt & pt, std::vector<double> & offlimits )
{
  bool ok = false;
  double xlo = mesh.GetLowX();
  double ylo = mesh.GetLowY();
  double xhi = mesh.GetHighX();
  double yhi = mesh.GetHighY();

  double buf = 0.033;
  if ( ! pt.minimum ) buf = 0.038;
  double bufx = (xhi-xlo)*buf;
  double bufy = (yhi-ylo)*buf;

  int n = offlimits.size() / 2;

  int orders[18*2] = { 0, 1, 
		   1, 0,
		   0, -1,
		   -1, 0,
		   1, 1,
		   -1, 1,
		   1, -1,
		   -1, -1,
		   0, 2,
		   2, 0,
		   1, 2,
		   2, 1,
		   -1, 2,
		   2, -1,
		   1, -2,
		   -2, 1,
		   -1, -2,
		   -2, -1 };

  double scales[10] = { 0.66, 0.8, 1.0, 1.25, 
		       1.33, 1.50, 1.66, 1.75, 
		       1.9, 2.0 };
  
  for ( int iscale = 0; iscale < 11; ++iscale )
    for ( int ij=0; ij<19; ++ij )
      {
	int i = orders[0+ij*2];
	int j = orders[1+ij*2];
	
	if ( i == 0 and j == 0 ) continue;
	double x = pt.x + bufx * i * scales[iscale];
	double y = pt.y + bufy * j * scales[iscale];
	if ( x < xlo + bufx ) continue;
	if ( x > xhi - bufx ) continue;
	if ( y < ylo + bufy ) continue;
	if ( y > yhi - bufy ) continue;
	
	
	ok = true;
	for ( int k=0; k<n; ++k )
	  {
	    double dx = x - offlimits[0+k*2];
	    double dy = y - offlimits[1+k*2];
	    if ( std::abs(dx) < bufx and std::abs(dy) < bufy )
	      {
		//std::printf("too close %12.3f %12.3f %12.3f %12.3f (%12.3f %12.3f)\n",
		//	      x,y,offlimits[0+k*2],offlimits[1+k*2],bufx,bufy);
		ok = false;
		break;
	      }
	  }
	if ( ok )
	  {
	    //std::printf("ok\n");
	    offlimits.push_back( x );
	    offlimits.push_back( y );
	    pt.x = x;
	    pt.y = y;
	    return;
	  };
      }
}



bool sort_by_i( endpt const & lhs, endpt const & rhs )
{
  return lhs.i < rhs.i;
}




int main( int argc, char ** argv)
{

  cli_options cli( read_options( argc, argv ) );


  //std::printf("ncon %i\n",(int)cli.connections.size());

  if ( (! cli.connections.size()) or cli.rezero )
    cli.meshfile = main_points( cli.meshfile.c_str(), 
				cli.perx or cli.pery, 
				cli.blanks,
				cli.extra_searches,
				cli.rezero );
  
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

      for ( std::vector<BoxCrd>::iterator 
	      p=cli.blanks.begin(), pend=cli.blanks.end(); 
	    p!=pend; ++p )
	{
	  p->xlo = std::max( p->xlo, mesh.GetLowX() );
	  p->ylo = std::max( p->ylo, mesh.GetLowY() );
	  p->xhi = std::min( p->xhi, mesh.GetHighX() );
	  p->yhi = std::min( p->yhi, mesh.GetHighY() );
	}



      std::vector<double> offlimits;

      int icon = 0;
      int ncon = cli.connections.size();
      double arcsum = 0.;
      std::vector<double> arclens;
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
	  for ( int i=0; i<201; ++i )
	    {
	      double t = i / 200.;
	      double x,y;
	      neb.GetXY( t, x, y );
	      ccdl::Mesh2dValue val( mesh.GetValue( x,y ) );
	      xyfile << FF(x) << FF(y) 
		     << " # " << FF(val.f) <<  "\n";
	    }

	  for ( int i=0; i<51; ++i )
	    {
	      double t = i / 50.;
	      double x,y;
	      neb.GetXY( t, x, y );
	      offlimits.push_back( x );
	      offlimits.push_back( y );
	    }

	  str.str("");
	  str.clear();
	  str << p->GetName() << ".1d.dat";

	  if ( p->prev < 0 ) 
	    {
	      if ( p != cli.connections.begin() )
		arclens.push_back( arcsum );
	      arcsum = 0.;
	    }

	  tfile.open( str.str().c_str() );
	  for ( int i=0; i<201; ++i )
	    {
	      double t = i / 200.;
	      double x,y;
	      neb.GetXY( t, x, y );
	      double arclen = neb.GetArcLength( 0, t );
	      ccdl::Mesh2dValue val( mesh.GetValue( x,y ) );
	      if ( i == 0 )
		p->t0.value = val.f;
	      else if ( i == 200 )
		p->t1.value = val.f;
	      tfile << FF(arclen+arcsum) << FF(val.f) 
		    << " # " << FF(val.x) << FF(val.y) << "\n";
	    }
	  p->t0.arclen = arcsum;
	  arcsum += neb.GetArcLength( 0., 1. );
	  p->t1.arclen = arcsum;
	};
      arclens.push_back( arcsum );


      for ( std::vector<endpt>::iterator 
	      p=cli.minima.begin(), pend=cli.minima.end(); 
	    p!=pend; ++p )
	{
	  offlimits.push_back( p->x );
	  offlimits.push_back( p->y );
	}
      for ( std::vector<endpt>::iterator 
	      p=cli.saddle.begin(), pend=cli.saddle.end(); 
	    p!=pend; ++p )
	{
	  offlimits.push_back( p->x );
	  offlimits.push_back( p->y );
	}
      for ( std::vector<endpt>::iterator 
	      p=cli.maxima.begin(), pend=cli.maxima.end(); 
	    p!=pend; ++p )
	{
	  offlimits.push_back( p->x );
	  offlimits.push_back( p->y );
	}




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
      std::vector<int> outline(colors.size(),110);
      // outline[2]  = 120;
      // outline[4]  = 120;
      // outline[6]  = 120;
      // outline[8]  = 120;
      // outline[10] = 120;


      typedef std::vector< connection >::iterator conit;
      typedef std::pair< conit, conit > conpair;
      typedef std::vector< conpair > convec;
      typedef convec::iterator convecit;

      
      convec paths;
      conit begin = cli.connections.begin();
      for ( conit p=begin+1, pend=cli.connections.end(); p != pend; ++p )
	{
	  if ( p->prev < 0 )
	    {
	      paths.push_back( std::make_pair( begin, p ) );
	      begin = p;
	    }
	}
      paths.push_back( std::make_pair( begin, cli.connections.end() ) );


      std::vector< RefinementPt > refinept;
      if ( cli.refine )
	{
	  double x0 = cli.refine->xlo;
	  double y0 = cli.refine->ylo;
	  double dx = cli.refine->dx;
	  double dy = cli.refine->dy;
	  int n = cli.refine->n;
	  for ( convecit ppath = paths.begin(); ppath != paths.end(); ++ppath )
	    {
	      conit path_begin = ppath->first;
	      conit path_end   = ppath->second;
	      for ( conit p = path_begin; p != path_end; ++p )
		{
		  ccdl::ParametricLegendre & neb = *(p->curve);
		  int nt = 201;
		  for ( int i=0; i<nt; ++i )
		    {
		      double t = ( i / (double)(nt-1) );
		      double x,y;
		      neb.GetXY( t, x, y );
		      int ixlo = (int)((x - x0) / dx + 0.5) - (n);
		      int iylo = (int)((y - y0) / dy + 0.5) - (n);
		      for ( int ix=ixlo; ix < ixlo + 2*n+1; ++ix )
			for ( int iy=iylo; iy < iylo + 2*n+1; ++iy )
			  refinept.push_back( RefinementPt( x0 + ix*dx, y0 + iy*dy ) );
		    }
		  std::sort( refinept.begin(), refinept.end() );
		  refinept.erase( std::unique( refinept.begin(), refinept.end() ), refinept.end() );
		}
	    }
	  std::ofstream cout;
	  cout.open( "refinement.dat" );
	  for ( std::vector<RefinementPt>::iterator 
		  pt=refinept.begin(), ptend=refinept.end();
		pt != ptend; ++pt )
	    {
	      //std::printf("%10.3f %10.3f\n",pt->x,pt->y);
	      cout << std::fixed << std::setw(8) << std::setprecision(3) << pt->x
		   << std::fixed << std::setw(8) << std::setprecision(3) << pt->y
		   << "\n";
	      offlimits.push_back( pt->x );
	      offlimits.push_back( pt->y );
	    }
	}
    




      {
	std::string gnuplot( cli.meshfile );
	gnuplot += ".sh";
	std::ofstream cout;
	cout.open( gnuplot.c_str() );
	cout << "#!/bin/sh" << "\n";

	std::vector<double> offlimits_backup( offlimits );


	std::vector<endpt> usedMinima;
	std::vector<endpt> usedSaddle;
	{
	  for ( convecit ppath = paths.begin(); ppath != paths.end(); ++ppath )
	    {
	      conit path_begin = ppath->first;
	      conit path_end   = ppath->second;
	      
	      for ( conit p = path_begin; p != path_end; ++p )
		{
		  if ( p->t0.minimum )
		    usedMinima.push_back( p->t0 );
		  else
		    usedSaddle.push_back( p->t0 );
		  if ( p->t1.minimum )
		    usedMinima.push_back( p->t1 );
		  else
		    usedSaddle.push_back( p->t1 );
		};
	    };
	  std::sort( usedMinima.begin(), usedMinima.end(), sort_by_i );
	  std::sort( usedSaddle.begin(), usedSaddle.end(), sort_by_i );
	  usedMinima.erase( std::unique( usedMinima.begin(), usedMinima.end() ), usedMinima.end() );
	  usedSaddle.erase( std::unique( usedSaddle.begin(), usedSaddle.end() ), usedSaddle.end() );
	  
	  
	  std::ofstream fdel;
	  fdel.open( "minpts_used.dat" );
	  for ( std::vector<endpt>::iterator 
		  it = usedMinima.begin(), itend = usedMinima.end();
		it != itend; ++it )
	    fdel << FF(it->x) << ", " << FF(it->y) << "\n";
	  fdel.close();
	  fdel.open( "tspts_used.dat" );
	  for ( std::vector<endpt>::iterator 
		  it = usedSaddle.begin(), itend = usedSaddle.end();
		it != itend; ++it )
	    fdel << FF(it->x) << ", " << FF(it->y) << "\n";
	  fdel.close();
	  
	}

	// bounds common to all paths

	std::stringstream cblo;
	std::stringstream cbhi;
	double flo = 0.;
	double fhi = 0.;
	if ( cli.minima.size() )
	  {
	    if ( (! cli.delunused) or (! usedMinima.size()) )
	      flo = mesh.GetValue( cli.minima[0].x, cli.minima[0].y ).f;
	    else
	      flo = mesh.GetValue( usedMinima[0].x, usedMinima[0].y ).f;

	    cblo << FF( flo );
	  }
	if ( cli.maxima.size() ) 
	  {
	    for ( std::vector<endpt>::reverse_iterator 
		    pm=cli.maxima.rbegin(), pmend=cli.maxima.rend(); 
		  pm != pmend; ++pm )
	      {
		bool in_blank = false;
		for ( std::vector<BoxCrd>::iterator 
			pb=cli.blanks.begin(), pbend=cli.blanks.end();
		      pb != pbend; ++pb )
		  {
		    if ( pm->x >= pb->xlo and pm->x <= pb->xhi and
			 pm->y >= pb->ylo and pm->y <= pb->yhi )
		      {
			in_blank = true;
			break;
		      }
		  }
		if ( ! in_blank )
		  {
		    fhi = mesh.GetValue( pm->x, pm->y ).f;
		    cbhi.str("");
		    cbhi.clear();
		    cbhi << FF( fhi );
		    break;
		  }
	      }
	  }
	if ( cli.saddle.size() or (! usedSaddle.size() ) ) 
	  {
	    for ( std::vector<endpt>::reverse_iterator 
		    pm=cli.saddle.rbegin(), pmend=cli.saddle.rend(); 
		  pm != pmend; ++pm )
	      {
		bool in_blank = false;
		for ( std::vector<BoxCrd>::iterator 
			pb=cli.blanks.begin(), pbend=cli.blanks.end();
		      pb != pbend; ++pb )
		  {
		    if ( pm->x >= pb->xlo and pm->x <= pb->xhi and
			 pm->y >= pb->ylo and pm->y <= pb->yhi )
		      {
			in_blank = true;
			break;
		      }
		  }
		if ( ! in_blank )
		  {
		    double tmp = mesh.GetValue( pm->x, pm->y ).f;
		    if ( tmp > fhi )
		      {
			fhi = tmp;
			fhi += 0.275 * ( fhi-flo );
			cbhi.str("");
			cbhi.clear();
			cbhi << FF( fhi );
			break;
		      };
		  }
	      }
	    //fhi = mesh.GetValue( cli.saddle.back().x, cli.saddle.back().y ).f;
	    //cbhi << FF( fhi );
	  }
	if ( cli.delunused and usedSaddle.size() )
	  {
	    fhi = 0.;
	    for ( std::vector<endpt>::reverse_iterator 
		    pm=usedSaddle.rbegin(), pmend=usedSaddle.rend(); 
		  pm != pmend; ++pm )
	      {
		double tmp = mesh.GetValue( pm->x, pm->y ).f;
		if ( tmp > fhi )
		  {
		    fhi = tmp;
		    fhi += 0.275 * ( fhi-flo );
		    cbhi.str("");
		    cbhi.clear();
		    cbhi << FF( fhi );
		    break;
		  };
	      }
	  }
	

	if ( cli.range_clo.size() ) cblo.str( cli.range_clo );
	if ( cli.range_chi.size() ) cbhi.str( cli.range_chi );






	// collect a list of grace commands as we go -- printed at the end
	std::vector<std::string> xmgrace_cmds;



	
	int ipath = 1;
	for ( convecit ppath = paths.begin(); ppath != paths.end(); ++ppath, ++ipath )
	  {
	    conit path_begin = ppath->first;
	    conit path_end   = ppath->second;

	    offlimits = offlimits_backup;

	    cout << "#\n# BEGIN PATH " << ipath << "\n#\n\n";

	    cout << "cat <<EOF | gnuplot" << "\n";
	    cout << "# make contours" << "\n";
	    cout << "set contour base" << "\n";
	    cout << "set cntrparam level incremental -50, 2.5, 50" << "\n";
	    cout << "unset surface" << "\n";
	    cout << "set table 'contours.dat'" << "\n";
	    cout << "splot '" << cli.meshfile << "' using 1:2:3" << "\n";
	    cout << "unset table" << "\n";
	    cout << "reset" << "\n";
	    cout << "unset key" << "\n";
	    cout << "#set palette rgbformulae 33,13,10" << "\n";
	    cout << "#set palette defined ( 0 '#0034a9', 1 '#1c61ff', 2 '#3371ff', 3 '#009cff', 4 '#00d8ff', 5 '#00ffba', 6 '#00ff6c', 7 '#00ff0c', 8 '#aeff00', 9 '#fff600', 10 '#ffc000', 11 '#ff9600', 12 '#ff6c00', 13 '#ff3c00', 14 '#ff0000', 15 '#7f0000', 16 '#530000')" << "\n";
	    cout << "set palette defined ( 0 '#0034a9', 1 '#6695ff', 2 '#96e8ff', 3 '#96ffe0', 4 '#60ff73', 5 '#fcff24', 6 '#feffce', 7 '#ffdc99', 8 '#ffa07b', 9 '#ff5353', 10 '#bd0009' )" << "\n";
	    

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
	    cout << "set style line 2 lc rgb 'black' pt 7 ps 2.2" << "\n";
	    cout << "set style line 12 lc rgb 'white' pt 7 ps 1.4" << "\n";
	    cout << "# maxpts" << "\n";
	    cout << "set style line 3 lc rgb 'black' pt 12 ps 1.8 lw 9" << "\n";
	    cout << "set style line 13 lc rgb 'white' pt 12 ps 1.2 lw 4" << "\n";
	    cout << "# tspts" << "\n";
	    cout << "set style line 4 lc rgb 'black' pt 2 ps 2.  lw 12" << "\n";
	    cout << "set style line 14 lc rgb 'white' pt 2 ps 1.4  lw 4" << "\n";
	    cout << "# refinement pts" << "\n";
	    cout << "set style line 5 lc rgb 'black' pt 1 ps 0.75  lw 1.8" << "\n";
	    cout << "set style line 15 lc rgb 'white' pt 1 ps 0.5  lw 0.5" << "\n";
	    cout << "# paths" << "\n";
	    cout << "set style line 100 lc rgb 'red' lw 8" << "\n";
	    cout << "set style line 110 lc rgb 'white' lw 16" << "\n";
	    cout << "set style line 120 lc rgb 'black' lw 16" << "\n";

	    cout << "" << "\n\n\n";


	    int iblank=0;
	    for ( std::vector<BoxCrd>::iterator pb=cli.blanks.begin(), pbend=cli.blanks.end(); pb!=pbend; ++pb )
	      {
		++iblank;
		cout << "set obj " << iblank 
		     << " rect from "
		     << std::fixed << std::setprecision(10) << pb->xlo
		     << ","
		     << std::fixed << std::setprecision(10) << pb->ylo
		     << " to "
		     << std::fixed << std::setprecision(10) << pb->xhi
		     << ","
		     << std::fixed << std::setprecision(10) << pb->yhi
		     << "fs solid border lc rgb \"white\" lw 4 fc rgb \"white\" front\n";
	      }

	    //cout << "labelMacro(i,x,y,l) = sprintf('set obj %d rect at %f,%f size char strlen(\"%s\"), char 1 fs solid noborder 01 front fc rgb \"black\" ; set label %d at %f,%f \"%s\" front center tc rgb \"white\" font \"Helvetica-Bold,18\"', i, x, y, l, i, x, y, l)\n\n";
	    
	    cout << "labelMacro(i,x,y,l) = sprintf('set obj %d rect at %f,%f size char 1*(strlen(\"%s\")>2?2.425:strlen(\"%s\")), char 0.9 fs solid noborder 01 front fc rgb \"black\"; set label %d at %f,%f \"%s\" front center tc rgb \"white\" font \"Helvetica-Bold,24\"', i, x, y, l, l, i, x, y, l)\n\n";


	    std::vector<endpt> myMinima;
	    std::vector<endpt> mySaddle;
	    for ( conit p = path_begin; p != path_end; ++p )
	      {
		if ( p->t0.minimum )
		  myMinima.push_back( p->t0 );
		else
		  mySaddle.push_back( p->t0 );
		if ( p->t1.minimum )
		  myMinima.push_back( p->t1 );
		else
		  mySaddle.push_back( p->t1 );
	      };

	    std::sort( myMinima.begin(), myMinima.end(), sort_by_i );
	    std::sort( mySaddle.begin(), mySaddle.end(), sort_by_i );
	    myMinima.erase( std::unique( myMinima.begin(), myMinima.end() ), myMinima.end() );
	    mySaddle.erase( std::unique( mySaddle.begin(), mySaddle.end() ), mySaddle.end() );

		  

	    int ilabel = 10;
	    for ( std::vector<endpt>::iterator it = myMinima.begin(), itend = myMinima.end();
		  it != itend; ++it )
	      {
		//std::printf("%i %s\n",it->i,it->label.c_str());
		if ( it->label.size() == 0 ) continue;
		endpt pt( *it );
		//std::printf("pre %12.4f %12.4f\n",pt.x,pt.y);
		GetValidPosition( mesh, pt, offlimits );
		//std::printf("post %12.4f %12.4f\n",pt.x,pt.y);
		ilabel++;
		cout << "eval labelMacro(" << ilabel << ", "
		     << FF(pt.x) << ", " << FF(pt.y) << ", \"" 
		     << pt.label << "\")\n";
	      }

	    for ( std::vector<endpt>::iterator it = mySaddle.begin(), itend = mySaddle.end();
		  it != itend; ++it )
	      {
		if ( it->label.size() == 0 ) continue;
		endpt pt( *it );
		//std::printf("pre %12.4f %12.4f\n",pt.x,pt.y);
		GetValidPosition( mesh, pt, offlimits );
		//std::printf("post %12.4f %12.4f\n",pt.x,pt.y);
		ilabel++;
		cout << "eval labelMacro(" << ilabel << ", "
		     << FF(pt.x) << ", " << FF(pt.y) << ", \"" 
		     << pt.label << "\")\n";
	      }

	    
	    cout << "" << "\n\n\n";
	    cout << "set terminal postscript eps enhanced solid color 'Helvetica' 26 size 6,5.4" << "\n";
	    cout << "set output 'path2d." << ipath << ".eps'" << "\n";
	    cout << "" << "\n";


	    cout << "set cbrange[" << cblo.str() << ":" << cbhi.str() << "]" << "\n";
	    cout << "set xrange[ " << cli.range_xlo << " : " << cli.range_xhi << " ]\n";
	    cout << "set yrange[ " << cli.range_ylo << " : " << cli.range_yhi << " ]\n";

	    cout << "set cblabel 'E (kcal/mol)'\n";
	    cout << "set ylabel 'Y-axis label'\n";
	    cout << "set xlabel 'X-axis label'\n";
	    cout << "p '" << cli.meshfile << "' using 1:2:3 with image,\\" << "\n";
	    cout << "  'contours.dat' w l lt -1 lw 1.5,\\" << "\n";


	    int icolor = 0;
	    int icon = 0;
	    std::stringstream xmgrace;
	    xmgrace << "xmgrace";

	    for ( conit p = path_begin, pend=path_end; p != pend; ++p, ++icon, ++icolor )
	      {

		icolor = icolor % (int)(colors.size());
		xmgrace << " \"" << p->GetName() << ".1d.dat\" \\\n";
		
		int nep = 1;
		if ( p+1 == pend )
		  nep = 2;
		else
		  {
		    std::vector< connection >::iterator q = p+1;
		    if ( q->prev < 0 )
		      nep = 2;
		  }
		
		for ( int iep = 0; iep < nep; ++iep )
		  {
		    endpt pt = p->t0;
		    if ( iep == 1 ) pt = p->t1;

		    double fperc = (pt.value - flo) / ( fhi-flo );
		    double buf = 0.05;
		    std::stringstream xline;
		    if ( fperc < 0.5 and ! pt.minimum ) 
		      {
			if ( fperc > 0.3 )
			  buf += 0.2;
			else
			  buf += 0.25;
			xline << " -pexec 'with line; line on; line loctype world; line "
			      << FF(pt.arclen) << ","
			      << FF(pt.value) << ","
			      << FF(pt.arclen) << ","
			      << FF(pt.value + (buf-0.01)*(fhi-flo)) << ";"
			      << " line linewidth 2.0; line linestyle 2; line color 1;"
			      << " line arrow 0; line arrow type 0; line arrow length 1.000000;"
			      << " line arrow layout 1.000000, 1.000000; line def'";
		      }
		    if ( fperc > 0.5 and ! pt.minimum ) buf -= 0.45;
		    buf *= (fhi-flo);
		    
		    xmgrace << " -pexec \"g0.s" << icolor << " line linewidth 4\" \\\n"
			    << " -pexec \'with string; string on; string loctype world; string g0; string " 
			    << FF(pt.arclen)
			    << "," 
			    << FF(pt.value + buf)
			    << "; string color 1; string rot 90; string font 4; string just 6;"
			    << " string char size 1.120000; string def \"" 
			    << std::fixed << std::setprecision(2) << pt.x
			    << "," 
			    << std::fixed << std::setprecision(2) << pt.y
			    << " (\\6" << pt.label << "\\0)\"\' \\\n";
		    xmgrace << xline.str();
		  };

		cout << "  '" << p->GetName() << ".2d.dat' w l ls " << outline[icolor] << ",\\" << "\n";
		cout << "  '" << p->GetName() << ".2d.dat' w l ls 100 lc rgb '" << colors[icolor] << "',\\" << "\n";
	      }
	    //std::cout << xmgrace.str() << "\n";
	    xmgrace_cmds.push_back( xmgrace.str() );

	    if ( ! cli.delunused )
	      {
		cout << "  'minpts.dat'   w p ls 2,\\" << "\n";
		cout << "  'minpts.dat'   w p ls 12,\\" << "\n";
		cout << "  'maxpts.dat'   w p ls 3,\\" << "\n";
		cout << "  'maxpts.dat'   w p ls 13,\\" << "\n";
		cout << "  'tspts.dat'    w p ls 4,\\" << "\n";
		cout << "  'tspts.dat'    w p ls 14";
	      }
	    else
	      {
		cout << "  'minpts_used.dat'   w p ls 2,\\" << "\n";
		cout << "  'minpts_used.dat'   w p ls 12,\\" << "\n";
		cout << "  'tspts_used.dat'    w p ls 4,\\" << "\n";
		cout << "  'tspts_used.dat'    w p ls 14";
	      }
	    if ( ! cli.refine )
	      cout << "\n";
	    else
	      {
		cout << ",\\" << "\n";
		cout << "  'refinement.dat'    w p ls 5,\\" << "\n";
		cout << "  'refinement.dat'    w p ls 15" << "\n";
		//cout << "  'refinement.dat'    w p ls 5" << "\n";

	      }
	    cout << "EOF" << "\n\n";

	  }


	// now we can print all the xmgrace commands
	{
	  int icon = 0;
	  for ( std::vector<std::string>::iterator 
		  p=xmgrace_cmds.begin(), pend=xmgrace_cmds.end(); p!=pend; ++p, ++icon )
	    {
	      cout << *p 
		   << " -pexec \"page size 612, 612\" -pexec \"view 0.12, 0.11, 0.94, 0.96\" \\\n"
		   << " -world -1 " << cblo.str() << " " << arclens[icon]+1 << " " << cbhi.str() << " \\\n"
		   << " -autoscale x \\\n" 
		   << " -pexec \'yaxis tick major 5\' \\\n"
		   << " -pexec \'yaxis tick minor ticks 1\' \\\n"
		   << " -pexec \'xaxis label font 4\' \\\n"
		   << " -pexec \'xaxis ticklabel font 4\' \\\n"
		   << " -pexec \'xaxis label \"Arc Length\"\' \\\n"
		   << " -pexec \"xaxis label char size 1.20\" \\\n"
		   << " -pexec \"xaxis ticklabel char size 1.20\" \\\n"
		   << " -pexec \'yaxis label font 4\' \\\n"
		   << " -pexec \'yaxis ticklabel font 4\' \\\n"
		   << " -pexec \'yaxis label \"E \\(kcal/mol\\)\"\' \\\n"
		   << " -pexec \"yaxis label char size 1.20\" \\\n"
		   << " -pexec \"yaxis ticklabel char size 1.20\" \\\n"
		   << " -pexec \"yaxis label place 0.000000, 0.095000\" \\\n"
		   << " -hardcopy -noprint -saveall path1d." << icon+1 << ".agr\n\n";
	      cout << "xmgrace -hardcopy -printfile path1d." << icon+1 << ".eps path1d." << icon+1 << ".agr &> /dev/null\n\n";
	      
	      cout << "echo \"Now run: evince  path2d." << icon+1 << ".eps &\"\n\n";
	      cout << "echo \"Now run: xmgrace path1d." << icon+1 << ".agr &\"\n\n";

	    }
	}
	int ncon = xmgrace_cmds.size();

	std::cout << "Now run: /bin/sh \"" << gnuplot << "\"\n";


	// now we can print all the tex file



	cout << "\n\n\n";
	cout << "cat <<EOF > figure.tex\n";

	cout << "\\makeatletter\n";
	cout << "\\newcommand{\\dontusepackage}[2][]{%\n";
	cout << "\\@namedef{ver@#2.sty}{9999/12/31}%\n";
	cout << "\\@namedef{opt@#2.sty}{#1}}\n";
	cout << "\\makeatother\n";
	cout << "\\dontusepackage{mciteplus}\n";
	cout << "\\documentclass[journal=jctcce,manuscript=article,layout=twocolumn]{achemso}\n";
	cout << "\\makeatletter\n";
	cout << "\\renewcommand*\\acs@etal@firstonly{\\acs@etal@truncatetrue}\n";
	cout << "\\renewcommand*\\acs@maxauthors{0}\n";
	cout << "\\makeatother\n";
	cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
	cout << "\\usepackage{amsmath}\n";
	cout << "\\usepackage{graphicx}\n";
	cout << "\\usepackage{multirow}\n";
	cout << "\\usepackage{bm}\n";
	cout << "\\usepackage{wasysym}\n";
	cout << "\\usepackage{placeins}\n";
	cout << "\\usepackage{framed}\n";
	cout << "\\usepackage{comment}\n";
	cout << "\\usepackage{color}\n";
	cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
	cout << "\\title{Title}\n";
	cout << "\\author{Author}\n";
	cout << "\\email{Email}\n";
	cout << "\\affiliation{Center for Integrative Proteomics Research, \n";
	cout << "BioMaPS Institute for Quantitative Biology and Department of \n";
	cout << "Chemistry and Chemical Biology, \n";
	cout << "Rutgers University, Piscataway, NJ 08854-8087 USA}\n";
	cout << "\\begin{document}\n";
	cout << "\\begin{abstract}\n";
	cout << "\\end{abstract}\n";

	for ( int icon=0; icon < ncon; ++icon )
	  {
	    conit path_begin = paths[icon].first;
	    conit path_end = paths[icon].second;

	    double emin = 1.e+30;
	    std::vector< std::pair<double,double> > fwdrev;
	    std::vector< endpt > pts;
	    pts.push_back( path_begin->t0 );
	    for ( conit p=path_begin; p!=path_end; ++p )
	      {
		emin = std::min(emin, p->t0.value);
		emin = std::min(emin, p->t1.value);
		pts.push_back( p->t1 );
		if ( ! p->t0.minimum )
		  {
		    double eprev = 1.e+10;
		    double enext = 1.e+10;
		    if ( p != path_begin )
		      {
			conit q = p-1;
			if ( q->t0.minimum ) eprev = q->t0.value;
		      }
		    if ( p->t1.minimum ) enext = p->t1.value;
		    fwdrev.push_back( std::make_pair( p->t0.value - eprev,  p->t0.value - enext ) );
		  }
		else if ( p+1 == path_end )
		  {
		    if ( ! p->t1.minimum )
		      fwdrev.push_back( std::make_pair( p->t1.value - p->t0.value,  1.e+10 ) );
		  }
	      };

	    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
	    cout << "\\begin{figure*}[tb]\n";
	    cout << "\\includegraphics[clip,width=3.75in]{path2d." << icon+1 << ".eps}\n";
	    cout << "\\hspace{-0.375in}\n";
	    cout << "\\includegraphics[clip,width=3.35in]{path1d." << icon+1 << ".eps}\n";
	    cout << "\\caption{Path " << icon+1 << "}\n";
	    cout << "\\end{figure*}\n";
	    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";

	    cout << "\\begin{table*}" << "\n";
	    cout << "\\begin{tabular}{c r r r r r r}" << "\n";
	    cout << "\\hline" << "\n";
	    cout << "Step & x & y & \\$E$ & $\\Delta E_{\\text{path}}$ & $\\Delta E_{\\text{fwd}}$ & $\\Delta E_{\\text{rev}}$ \\\\\\\\" << "\n";
	    cout << "& Unit & Unit & kcal/mol & kcal/mol & kcal/mol & kcal/mol  \\\\\\\\" << "\n";
	    cout << "\\hline" << "\n";
	    int tscnt = 0;
	    for ( std::vector<endpt>::iterator p=pts.begin(), pend=pts.end(); p!=pend; ++p )
	      {
		cout << "({\\bf " << p->label << "}) & "
		     << std::fixed << std::setprecision(2) << p->x << " & " 
		     << std::fixed << std::setprecision(2) << p->y << " & " 
		     << std::fixed << std::setprecision(2) << p->value << " & " 
		     << std::fixed << std::setprecision(2) << p->value - emin << " & " ;
		if ( p->minimum )
		  cout << " $\\cdots$ & $\\cdots$";
		else
		  {
		    double fwd = 1.e+30;
		    double rev = 1.e+30;
		    if ( tscnt + 1 <= (int)( fwdrev.size() ) )
		      {
			fwd = fwdrev[tscnt].first;
			rev = fwdrev[tscnt].second;
			++tscnt;
		      }
		    if ( std::abs( fwd ) > 1.e+5 )
		      cout << " $\\cdots$ &";
		    else
		      cout << std::fixed << std::setprecision(2) << fwd << " & ";
		    if ( std::abs( rev ) > 1.e+5 )
		      cout << " $\\cdots$";
		    else
		      cout << std::fixed << std::setprecision(2) << rev ;
		  
		  }
		cout << " \\\\\\\\" << "\n";
	      }

	    cout << "\\hline" << "\n";
	    cout << "\\end{tabular}" << "\n";
	    cout << "\\caption{Path 1}" << "\n";
	    cout << "\\end{table*}" << "\n";
	    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
	    cout << "\\clearpage\n";
	    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";

	  }
	cout << "\\end{document}\n";

	cout << "EOF\n\n";

	cout << "latex figure.tex &>/dev/null; latex figure.tex &>/dev/null; \n";
	cout << "dvips -j0 -Ppdf -G0 -tletter -D 1200 -Z  figure.dvi &>/dev/null\n";
	cout << "ps2pdf -dCompatibilityLevel=1.3 -sPAPERSIZE=letter -dMAxSubsetPct=100 -dSubsetFonts=true -dEmbedAllFonts=true -dDetectBlends=true -dOptimize=true -dDownsampleColorImages=true -dColorImageResolution=1200 -dColorImageDownsampleType=/Average -dColorImageFilter=/FlateEncode -dAutoFilterColorImages=false -dAntiAliasColorImages=false -dColorImageDownsampleThreshold=1.50000 -dDownsampleGrayImages=true -dGrayImageResolution=1200 -dGrayImageDownsampleType=/Average -dGrayImageFilter=/FlateEncode -dAutoFilterGrayImages=false -dAntiAliasGrayImages=false -dGrayImageDownsampleThreshold=1.50000 -dDownsampleMonoImages=true -dMonoImageResolution=1200 -dMonoImageDownsampleType=/Average -dMonoImageFilter=/FlateEncode -dAutoFilterMonoImages=false -dAntiAliasMonoImages=false -dMonoImageDownsampleThreshold=1.50000 figure.ps &>/dev/null\n";


	cout << "echo \"Now run: evince figure.pdf &\"\n\n";


      }

    }
  

  
  return 0;
}
