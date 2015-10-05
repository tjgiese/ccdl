#!/usr/bin/env python
from atizer import *
from atizer.dbase import *

class pmf2d(autoprog):
    def __init__( self, srcdir=None ):
        super( pmf2d, self ).__init__( "pmf2d", srcdir )

        self.copyright_holder = "Timothy J. Giese"
        self.license = licenses.MIT

        ## @brief List of autolib objects representing library dependencies
        self.libs = [ ccdl() , nlopt() ]

        ## @brief List of filenames to be distributed, but not installed
        self.dist_noinst_SCRIPTS = []
        self.EXTRA_DIST = []

        ## @brief If True, then compile the target without installing it
        #         (default False)
        self.noinst = False

        ## @brief If True, "make doxygen-doc" create documentation html
        self.doxygen = False

        # self.enable_openmp()
        # self.enable_mpi()

package = autopackage(
    "pmf2d",
    targets=[ pmf2d( here() ) ],
    subdirs=[],
    version="0.1",
    apiversion="0:0:0")

if __name__ == "__main__":
    package.configure()

