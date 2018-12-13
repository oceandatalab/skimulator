import mod_plot
import params as p
import numpy
import glob
import os


# Initialize variables
indatadir = p.outdatadir
config = p.config
modelbox = p.modelbox
#modelbox = [-5, 5, 75, 85]
#modelbox = [-65, 55, 40, 45]
modelbox[0] = numpy.mod(modelbox[0] + 180.0, 360.0) - 180.0
modelbox[1] = numpy.mod(modelbox[1] + 180.0, 360.0) - 180.0
filesgrid = os.path.join(indatadir, '{}_grid'.format(config))
listfile = glob.glob('{}*.nc'.format(filesgrid))
output =  os.path.join(indatadir, 'Grid_{}.png'.format(config))

# Plot all Grids:
mod_plot.plot_all(listfile, modelbox, output)
print('Plot saved in {}'.format(output))

# Plot first 6 passes"
listfile = sorted(listfile)
output =  os.path.join(indatadir, 'Grid_6{}.png'.format(config))
mod_plot.plot_all(listfile[:6], modelbox, output)
print('Plot saved in {}'.format(output))

