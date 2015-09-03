#ifndef _XyzFile_hpp_
#define _XyzFile_hpp_

#include <iostream>
#include <vector>

namespace ccdl
{
  void ReadXyz( std::istream & cin, 
		std::vector<int> & z, 
		std::vector<double> & crd, 
		std::vector<double> & xtra );

  void WriteXyz( std::ostream & cout, 
		 int const nat,
		 int const * z,
		 double const * crd,
		 double const * xtra = NULL );
}

#endif
