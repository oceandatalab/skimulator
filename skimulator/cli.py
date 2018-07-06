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
import sys
import skimulator.mod_tools as mod_tools
# import argparse
import logging
logger = logging.getLogger()
handler = logging.StreamHandler()


def run_script():
    """Run SKIM Simulator"""
    import skimulator.run_simulator as run_simulator
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', action='store_true', default=False,
                        help='Display debug log messages')
    args = parser.parse_args()
    if args.debug is True:
        logger.setLevel(logging.DEBUG)
    '''
    if len(sys.argv) < 2:
        logger.error('Please specify a parameter file')
        sys.exit(1)
    else:
        file_param = str(sys.argv[1])

    p = mod_tools.load_python_file(file_param)
    run_simulator.run_simulator(p)
