'''
   FIG. 4: Radial currents and instrumental noise.
'''

import netCDF4
import numpy
from matplotlib import pyplot
import skimsimulator.rw_data as rw
import glob
import cartopy
import os

# Initialize color
listcolor = ['c', 'y', 'b', 'g', 'k', 'r', 'c', 'y']

# List files
indatadir = '/mnt/data/project/'
indatadir = '/tmp/key/project/'
indatadir = os.path.join(indatadir, 'skim', 'skim_output')
config="WW3_GS_6b108az"
config="WW3_GS_8b105az"
filesgrid = os.path.join(indatadir, '{}_'.format(config))
ipass = 59
indatapath = '{}c01_p{:03d}.nc'.format(filesgrid, ipass)
listfile = glob.glob(indatapath)
listfile = sorted(listfile)
outdatadir = '../'
modelbox = [-90, -70., 32., 40.]
scale = 8
#modelbox = [270, 290., 32., 40.]

# Prepare figure
pyplot.figure(figsize=(10, 15))
projection = cartopy.crs.Mercator()
transform = cartopy.crs.PlateCarree()
projection = transform
ax1 = pyplot.subplot(111, projection=projection)
#ax.coastlines()
ax1.add_feature(cartopy.feature.OCEAN, zorder=1)
ax1.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
#ax.set_extent([modelbox[0], modelbox[1],  modelbox[2], modelbox[3]],
#              projection)
ax1.set_extent([-74., -70., 34, 37], projection)
gl = ax1.gridlines(crs=transform, draw_labels=True, color='gray',
             linestyle='--', alpha=0.5)
gl.xlabels_top = False
gl.ylabels_left = False

data = netCDF4.Dataset(indatapath, 'r')
indatapath = '{}grid_p{:03d}.nc'.format(filesgrid, ipass)
datag = netCDF4.Dataset(indatapath, 'r')
lon = data['lon'][:]
lon[lon > 180] = lon[lon > 180] - 360
lat = data['lat'][:]
ur = data['ur_model'][:]
angle = datag['radial_angle'][:]
ur = data['ur_instr'][:]
urnoise = ur + urinstr
uur = ur * numpy.cos(corrangle)
vur = ur * numpy.sin(corrangle)
uurnoise = urnoise * numpy.cos(corrangle)
vurnoise = urnoise * numpy.sin(corrangle)
for i in range(numpy.shape(lon)[1]):
    style_color = '{}+'.format(listcolor[i])
    pyplot.plot(lon[:, i], lat[:, i], style_color,
                transform=cartopy.crs.PlateCarree())
    corrangle = angle[:, i]
    pyplot.quiver(lon[:, i], lat[:, i], uur[:, i], vur[:, i], color='red',
                  scale=scale, transform=transform)
    pyplot.quiver(lon[:, i], lat[:, i], uurnoise[:, i], vurnoise[:, i],
                  color='green', scale=scale, transform=transform)
pyplot.savefig(os.path.join(outdatadir, 'Fig3.png'))
