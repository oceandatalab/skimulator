import sys
import os
import shutil

#if  (len(sys.argv) < 2):
#    file_param = os.getcwd() + os.sep + 'example' + os.sep + 'params_example.txt'
#    print("no params file specified, default is " + file_param)
#else:
#    file_param = str(sys.argv[1])
#if os.path.isfile(file_param):
    #    basedir=os.path.dirname(swotsimulator.__file__)
#    shutil.copyfile(file_param, 'params.py')
    # os.path.join(basedir,'params.py'))
#else:
#    print("Error: No such file: '%s'" % file_param)
#    sys.exit()
import params as p
import logging
import argparse
logger = logging.getLogger()
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)
logger.setLevel(logging.INFO)
import skimsimulator.run_simulator as run_simulator

def parse_cli():
    """"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', action='store_true', default=False,
                        help='Display debug log messages')
    args = parser.parse_args()
    if args.debug is True:
        logger.setLevel(logging.DEBUG)


parse_cli()
run_simulator.run_simulator(p)
