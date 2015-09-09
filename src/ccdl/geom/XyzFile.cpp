#include "XyzFile.hpp"
#include "../constants.hpp"
#include <string>
#include <sstream>
#include <cstddef>
#include <iomanip>
#include <cstdlib>
#include <cstdio>

void ccdl::ReadXyz
( std::istream & cin, 
  std::vector<int> & z, 
  std::vector<double> & crd, 
  std::vector<double> & xtra )
{
  std::string elems("H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKr");
  int nat = 0;
  std::string line;
  std::stringstream sline;
  std::getline( cin, line );
  sline.str(line);
  sline.clear();
  sline >> nat;


  std::getline( cin, line );
  crd.resize(3*nat);
  xtra.resize(0);
  z.resize(nat);
  for ( int i=0; i<nat; ++i )
    {
      std::string e;
      std::getline( cin, line );
      sline.str(line);
      sline.clear();
      sline >> e >> crd[0+i*3] >> crd[1+i*3] >> crd[2+i*3];
      sline >> std::ws; // eat white space
      int c = sline.peek();
      if ( c == EOF )
        {
          if ( xtra.size() > 0 )
            {
              std::cerr << "Could not read charge field on line '" 
			<< line << "'\n";
              std::abort();
            };
        }
      else
        {
          double q = 0;
          sline >> q;
          xtra.push_back( q );
        };

      if ( e.size() == 1 )
        e += " ";
      std::size_t n = elems.find( e );
      if ( n == std::string::npos )
        {
          if ( e.compare("X") != 0 )
            std::cerr << "Could not determine element '" << e << "'\n";
          z[i] = 0;
        }
      else
        {
          z[i] = n/2 + 1;
        };
    }
  for ( int i=0; i<3*nat; ++i )
    crd[i] *= ccdl::AU_PER_ANGSTROM;
}



void ccdl::WriteXyz
( std::ostream & cout, 
  int const nat,
  int const * z,
  double const * crd,
  double const * xtra )
{
  std::string elems("H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKr");
  bool do_xtra = (xtra != NULL);
  cout << nat << "\n\n";
#define FMT std::fixed << std::setw(17) << std::setprecision(10) 
  for ( int a=0; a<nat; ++a )
    {
      std::string ele = "X ";
      if ( z[a] > 0 )
	ele = elems.substr( 2*(z[a]-1), 2 );

      cout << ele 
	   << FMT << crd[0+a*3] / ccdl::AU_PER_ANGSTROM
	   << FMT << crd[1+a*3] / ccdl::AU_PER_ANGSTROM
	   << FMT << crd[2+a*3] / ccdl::AU_PER_ANGSTROM;
      if ( do_xtra )
	cout << FMT << xtra[a];
      cout << "\n";
    };
#undef FMT
}

