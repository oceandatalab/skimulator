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
import errno
import shutil
import logging
from setuptools import setup
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
readme_file = os.path.join(project_dir, 'README')
package_dir = os.path.join(project_dir, __package_name__)
init_file = os.path.join(package_dir, '__init__.py')
share_dir = os.path.join(package_dir, 'share')
version_file = os.path.join(share_dir, 'VERSION.txt')

# Regenerate a version file from git history
if has_git:
    import subprocess
    import datetime
    _githash = (git_exe, 'rev-parse', '--short', 'HEAD')
    githash = subprocess.check_output(_githash).decode('utf-8').strip()
    gitrev = (git_exe, 'rev-list', 'HEAD', '--count')
    commits = subprocess.check_output(gitrev).decode('utf-8').strip()

    # Note: this block could be replaced by
    #   os.makedirs(share_dir, exist_ok=True)
    # if support for Python < 3.5 is dropped
    if not os.path.isdir(share_dir):
        try:
            os.makedirs(share_dir)
        except OSError:
            _, e, _ = sys.exc_info()
            if e.errno == errno.EEXIST:
                pass

    with open(version_file, 'wt') as f:
        f.write('{}\n'.format(commits))
        f.write('{}\n'.format(githash))
        f.write('{}\n'.format(datetime.datetime.utcnow().strftime('%Y-%m-%d')))

# - Read in the package version and author fields from the Python
#  main __init__.py file:
#
# IMPORTANT: this must be done AFTER trying to generate the VERSION.txt file
# as the __version__ variable contained in the __init__.py will not only have
# the major and minor, but also the commit number if VERSION.txt exists.
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

cmds = ['skimul2b = {}.cli:run_script'.format(__package_name__),
        'skimul2c = {}.cli:run_l2c'.format(__package_name__),
        'skimul2d = {}.cli:run_l2d'.format(__package_name__),
        'interpl2d = {}.cli:interpolate_l2d'.format(__package_name__),
        ]

# VERSION.txt must be added to the package if the file has been generated
pkg_data = None
pkg_data = {__package_name__: ['share/coeff.npy',
                               'share/coeffr.npy',
                               'share/Spline_128_64_TED_TAS_6_DEG.npy',
                               'share/Spline_128_64_TED_TAS_12_DEG.npy',
                               'share/Spline_128_64_TED_CB_12_DEG.npy',
                               'share/Spline_128_64_TED_CB_6_DEG.npy',
                               ]}
if os.path.exists(version_file):
    pkg_data[__package_name__].append('share/VERSION.txt')

setup(name=__package_name__,
      version=metadata['__version__'],
      description=metadata['__description__'],
      author=metadata['__author__'],
      author_email=metadata['__author_email__'],
      url=metadata['__url__'],
      license='COPYING',
      keywords=metadata['__keywords__'],
      long_description=long_description,
      packages=(__package_name__,),
      install_requires=requirements,
      setup_require=(),
      entry_points={'console_scripts': cmds},
      extras_require=optional_dependencies,
      package_data=pkg_data,
      )
