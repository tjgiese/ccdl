#!/usr/bin/env python
from atizer import *
from atizer.dbase import *
from src import package as src

package = autopackage(
    "ccdl",
    targets=[],
    subdirs=[src],
    version="0.1",
    apiversion="0:0:0")

if __name__ == "__main__":
    package.configure()
