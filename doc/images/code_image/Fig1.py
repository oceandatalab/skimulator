'''
FIG. 1: 5-day worth of SKIM simulated data in a global configuration with the science orbit.
'''

import netCDF4
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import skimsimulator.rw_data as rw
import glob
import cartopy
import os

# List files
outdatadir = '/tmp/key'
outdatadir = '/mnt/data/project/'
outdatadir = os.path.join(outdatadir, 'skim', 'skim_output')
config="WW3_GLOB_6b108az"
filesgrid = os.path.join(outdatadir, '{}_grid'.format(config))
indatapath = '{}_*'.format(filesgrid)
listfile = glob.glob(indatapath)
listfile = sorted(listfile) #.sort()
outdatadir = '../'
modelbox = [0., 360., -90., 90.]

# Prepare figure
pyplot.figure(figsize=(10, 15))
ax = pyplot.axes(projection=cartopy.crs.Orthographic(0, 0))
#ax.set_extent([modelbox[0], modelbox[1],  modelbox[2], modelbox[3]],
#              crs=cartopy.crs.Mercator())
#ax.coastlines()
ax.add_feature(cartopy.feature.OCEAN, zorder=1)
ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
ax.set_global()
#ax.gridlines()
for ifile in listfile[:90]:
    data = netCDF4.Dataset(ifile, 'r')
    lon = data['lon'][:, 0]
    lon[lon > 180] = lon[lon > 180] - 360
    lat = data['lat'][:, 0]
    lon_nadir = data['lon_nadir'][:]
    lon_nadir[lon_nadir > 180] = lon_nadir[lon_nadir > 180] - 360
    lat_nadir = data['lat_nadir'][:]
    pyplot.plot(lon[:], lat[:], '.', color='#4EE2EC', markersize=0.3,
                transform=cartopy.crs.PlateCarree())
    pyplot.plot(lon_nadir, lat_nadir, '.', color='#565051', markersize=0.3,
                transform=cartopy.crs.PlateCarree())
pyplot.savefig(os.path.join(outdatadir, 'Fig1.png'))
