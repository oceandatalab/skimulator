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
import argparse
import logging

logger = logging.getLogger(__name__)


def run_script():
    """Run SKIM Simulator"""
    import skimulator.run_simulator as run_simulator

    # Setup logging
    main_logger = logging.getLogger()
    main_logger.handlers = []
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    main_logger.addHandler(handler)
    main_logger.setLevel(logging.INFO)

    parser = argparse.ArgumentParser()
    parser.add_argument('params_file', nargs='?', type=str, default=None,
                        help='Path of the parameters file')
    parser.add_argument('--die-on-error', action='store_true', default=False,
                        help='Force simulation to quit on first error')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='Display debug log messages')

    args = parser.parse_args()

    if args.params_file is None:
        logger.error('Please specify a parameter file')
        sys.exit(1)

    if args.debug is True:
        main_logger.setLevel(logging.DEBUG)

    file_param = args.params_file

    p = mod_tools.load_python_file(file_param)
    try:
        run_simulator.run_simulator(p, args.die_on_error)
    except KeyboardInterrupt:
        logger.error('\nInterrupted by user (Ctrl+C)')
        sys.exit(1)
    sys.exit(0)


def run_l2c():
    """Run L2C reconstruction"""
    import skimulator.regridding as regridding

    # Setup logging
    main_logger = logging.getLogger()
    main_logger.handlers = []
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    main_logger.addHandler(handler)
    main_logger.setLevel(logging.INFO)

    parser = argparse.ArgumentParser()
    parser.add_argument('params_file', nargs='?', type=str, default=None,
                        help='Path of the parameters file')
    parser.add_argument('--die-on-error', action='store_true', default=False,
                        help='Force simulation to quit on first error')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='Display debug log messages')

    args = parser.parse_args()

    if args.params_file is None:
        logger.error('Please specify a parameter file')
        sys.exit(1)

    if args.debug is True:
        main_logger.setLevel(logging.DEBUG)

    file_param = args.params_file

    p = mod_tools.load_python_file(file_param)
    try:
        regridding.run_l2c(p, args.die_on_error)
    except KeyboardInterrupt:
        logger.error('\nInterrupted by user (Ctrl+C)')
        sys.exit(1)
    sys.exit(0)
