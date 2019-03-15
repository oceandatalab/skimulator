import numpy
import mod_diag
import glob
import os
import matplotlib

#import params as p
#indatadir = p.outdatadir
#config = p.config
#modelbox = p.modelbox

modelbox=None
config = 'WW3_EQ_metop_2018_8a'
indatadir = '/tmp/key/data/skim_eq_output/{}'.format(config)
filesgrid = os.path.join(indatadir, '{}_grid'.format(config))

listfile = glob.glob('{}*.nc'.format(filesgrid))
files = os.path.join(indatadir, '{}_c01'.format(config))
listfiles = glob.glob('{}*.nc'.format(files))
listfiles = sorted(listfiles)
output =  os.path.join(indatadir, 'Grid_{}.png'.format(config))
output2 =  os.path.join(indatadir, 'c_{}.png'.format(config))
print(config)
# listfiles.remove('/tmp/key/data/skim_eq_output/WW3_EQ_metop_2018_8a_c01_p022.nc')

mod_diag.diag_rms(listfiles[:], modelbox, output)
