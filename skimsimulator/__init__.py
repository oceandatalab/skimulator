# #!/usr/bin/env python
# =======================================================================
#                        General Documentation

"""Utilities for SKIM Science Simulator for the ocean

   Some useful online help commands for the package:
   * help(skimsimulator):  Help for the package.  A list of all modules in
     this package is found in the "Package Contents" section of the
     help output.
   * help(skimsimulator.M):  Details of each module "M", where "M" is the
     module's name.

#-----------------------------------------------------------------------
#                       Additional Documentation
# Authors: Lucile Gaultier
#
# Modification History:
# - Mar 2017:  Original by Lucile Gaultier, ODL
#
# Notes:
# - Written for Python 2.7,  Python 3.5, tested with Python 2.7, Python 3.5
#
# Copyright (c)
#
#-----------------------------------------------------------------------
"""
# -----------------------------------------------------------------------

# ---------------- Module General Import and Declarations ---------------

# - Set module version to package version:

__version__ = '$1.00: Original version $'
__author__ = 'Lucile Gaultier <lucile.gaultier@oceandatalab.com>'
__date__ = '2017-03-01'
__email__ = 'lucile.gaultier@oceandatalab.com'
__url__ = ''
__description__ = ('SKIM Simulator')
__author_email__ = ('lucile.gaultier@oceandatalab.com')
__keywords__ = ()

# - If you're importing this module in testing mode, or you're running
#  pydoc on this module via the command line, import user-specific
#  settings to make sure any non-standard libraries are found:

import os
import sys
if (__name__ == "__main__") or \
   ("pydoc" in os.path.basename(sys.argv[0])):
    import user


# - Find python version number
__python_version__ = sys.version[:3]

# - Import numerical array formats
try:
    import numpy
except ImportError:
    print(''' Numpy is not available on this platform,
          ''')

# - Import scientific librairies
try:
    import scipy
except ImportError:
    print("""Scipy is not available on this platform,
          """)


# - Import netcdf reading librairies
try:
    import netCDF4
except ImportError:
    print(''' netCDF4 is not available on this machine,
          ''')
    # reading and writing netcdf functions in rw_data.py won't work'''
