import numpy
import mod_diag
import params as p
import glob
import os
import matplotlib

indatadir = p.outdatadir
config = p.config
modelbox = p.modelbox
filesgrid = os.path.join(indatadir, '{}_grid'.format(config))
listfile = glob.glob('{}*.nc'.format(filesgrid))
files = os.path.join(indatadir, '{}_c01'.format(config))
listfiles = glob.glob('{}*.nc'.format(files))
listfiles = sorted(listfiles)
output =  os.path.join(indatadir, 'Grid_{}.png'.format(config))
output2 =  os.path.join(indatadir, 'c_{}.png'.format(config))
print(p.config)
# listfiles.remove('/tmp/key/data/skim_eq_output/WW3_EQ_metop_2018_8a_c01_p022.nc')

mod_diag.diag_rms(listfiles[:], modelbox, output)
