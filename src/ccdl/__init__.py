#!/usr/bin/env python
from atizer import *
from atizer.dbase import *

package = autopackage(
    "ccdl",
    targets=[ ccdl( here() ) ],
    subdirs=[],
    version="0.1",
    apiversion="0:0:0")

package.license = licenses.Nonfree

if __name__ == "__main__":
    package.configure()
