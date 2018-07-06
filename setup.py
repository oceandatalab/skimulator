"""
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

"""Build and install the SKIM Simulator for Ocean Science package."""

import os
import sys
import shutil
import logging
from setuptools import setup, find_packages
logger = logging.getLogger()
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

# Check Python version
if not 3 == sys.version_info[0]:
    logger.error('This package is only available for Python 3.x')
    sys.exit(1)

__package_name__ = 'skimulator'
project_dir = os.path.dirname(__file__)
git_exe = shutil.which('git')
git_dir = os.path.join(project_dir, '.git')
has_git = (git_exe is not None and os.path.isdir(git_dir))
# version_file = os.path.join(project_dir, 'VERSION.txt')
readme_file = os.path.join(project_dir, 'README')
package_dir = os.path.join(project_dir, __package_name__)
init_file = os.path.join(package_dir, '__init__.py')

# - Read in the package version and author fields from the Python
#  main __init__.py file:
metadata = {}
with open(init_file, 'rt') as f:
    exec(f.read(), metadata)

requirements = []
with open('requirements.txt', 'r') as f:
    lines = [x.strip() for x in f if 0 < len(x.strip())]
    requirements = [x for x in lines if x[0].isalpha()]

with open(readme_file, 'rt') as f:
    long_description = f.read()

optional_dependencies = {'plot': ['matplotlib', ], 'carto': ['matplotlib',
                         'cartopy']}

cmds = ['skimulator = {}.cli:run_script'.format(__package_name__),
        ]

setup(name='skimulator',
      version=metadata['__version__'],
      description=metadata['__description__'],
      author=metadata['__author__'],
      author_email=metadata['__author_email__'],
      url=metadata['__url__'],
      license='COPYING',
      keywords=metadata['__keywords__'],
      long_description=long_description,
      packages=find_packages(),
      install_requires=requirements,
      setup_require=(),
      entry_points={'console_scripts': cmds},
      extras_require=optional_dependencies,
      )
