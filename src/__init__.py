#!/usr/bin/env python
from atizer import *
from atizer.dbase import *
from ccdl import package as ccdl

package = autopackage(
    "src",
    targets=[],
    subdirs=[ccdl],
    version="0.1",
    apiversion="0:0:0")

if __name__ == "__main__":
    package.configure()
