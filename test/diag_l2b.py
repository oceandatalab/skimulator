import numpy
import mod_diag
import glob
import os
import matplotlib
import json
import sys

if len(sys.argv) < 1:
    print('Provide json file for diagnostics')
    sys.exit(1)
file_param = sys.argv[1]
with open(file_param, 'r') as f:
    params = json.load(f)

modelbox = params['l2b']['modelbox']
config = params['l2b']['config']
indatadir = params['l2b']['indatadir']
outdatadir = params['l2b']['outdatadir']
list_angle = params['l2b']['list_angle']
# indatadir = params['/tmp/key/data/skim_eq_output/{}'.format(config)
filesgrid = os.path.join(indatadir, '{}_grid'.format(config))

listfile = glob.glob('{}*.nc'.format(filesgrid))
files = os.path.join(indatadir, '{}_c'.format(config))
listfiles = glob.glob('{}*.nc'.format(files))
listfiles = sorted(listfiles)
output =  os.path.join(outdatadir, 'Grid_{}.png'.format(config))
output2 =  os.path.join(outdatadir, 'c_{}.png'.format(config))
print(config)
# listfiles.remove('/tmp/key/data/skim_eq_output/WW3_EQ_metop_2018_8a_c01_p022.nc')

mod_diag.diag_rms(listfiles[:], modelbox, output, list_angle)
