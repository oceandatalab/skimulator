# =======================================================================
#                        General Documentation

"""Utilities for SKIM Science Simulator for the ocean

   Some useful online help commands for the package:
   * help(skimulator):  Help for the package.  A list of all modules in
     this package is found in the "Package Contents" section of the
     help output.
   * help(skimulator.M):  Details of each module "M", where "M" is the
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
#-----------------------------------------------------------------------
Copyright (C) 2017-2018 OceanDataLab
This file is part of skimulator.

skimulator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

skimulator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with skimulator.  If not, see <http://www.gnu.org/licenses/>.
"""
# -----------------------------------------------------------------------

# ---------------- Module General Import and Declarations ---------------

# - Set module version to package version:

__version__ = '1.0'
__author__ = 'Lucile Gaultier <lucile.gaultier@oceandatalab.com>'
__date__ = '2017-03-01'
__email__ = 'lucile.gaultier@oceandatalab.com'
__url__ = ''
__description__ = ('SKIM Simulator')
__author_email__ = ('lucile.gaultier@oceandatalab.com')
__keywords__ = ()

import os

# Try to improve accuracy of the version number by reading the build number in
# the share/VERSION.txt file if it is available
try:
    import pkg_resources
    _version_path = pkg_resources.resource_filename('skimulator',
                                                    'share/VERSION.txt')
    if os.path.exists(_version_path):
        with open(_version_path, 'rt') as f:
            _commits = f.readline().strip()
            __version__ = '{}.{}'.format(__version__, _commits)
            _ = f.readline()  # Commit hash, useful for debug
            __date__ = f.readline().strip()
except ImportError:
    pass

# - If you're importing this module in testing mode, or you're running
#  pydoc on this module via the command line, import user-specific
#  settings to make sure any non-standard libraries are found:

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
