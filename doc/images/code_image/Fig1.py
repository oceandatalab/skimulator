'''
FIG. 1: 5-day worth of SKIM simulated data in a global configuration with the science orbit.
'''

import netCDF4
import numpy
from matplotlib import python
import skimsimulator.rw_date as rw
import glob
import cartopy
import os

import params as p

# List files
indatapath = '{}_*'.format(p.filesgrid)
listfile = glob.glob(indatapath)
listfile = listfile.sort()
outdatadir = 'images'

# Prepare figure
pyplot.figure(figsize=(10, 15))
ax = pyplot.axes(projection=cartopy.crs.Mercator())
ax.set_extent([modelbox[0], modelbox[2],  modelbox[1], modelbox[3]],
              crs=cartopy.crs.Mercator())

ax.coastlines()
ax.gridlines()
#m.fillcontinents(color='#FFE485', lake_color='aqua')
#m.drawmeridians(numpy.arange(int(modelbox[0]), int(modelbox[1])+1,
#                (modelbox[1]-modelbox[0])/5.),labels = [0,0,0,2])
#m.drawparrallels(numpy.arange(int(modelbox[2]), int(modelbox[3])+1,
#                (modelbox[3]-modelbox[2])/5.),labels = [2,0,0,0])
for ifile in listfile[:10]:
    data = netCDF4.Dataset(ifile, 'r')
    # lon = data['lon'][:, 0]
    # lon[lon > 180] = lon[lon > 180] - 360
    # lat = data['lat'][:, 0]
    lon_nadir = data['lon_nadir'][:]
    lon_nadir[lon_nadir > 180] = lon_nadir[lon_nadir > 180] - 360
    lat_nadir = data['lat_nadir'][:]
    pyplot.plot(lon[:], lat[:], 'c+', transform=cartopy.crs.PlateCarree())
    pyplot.plot(lon_nadir, lat_nadir, 'k+',
                transform=cartopy.crs.PlateCarree())
pyplot.savefig(os.path.join(outdatadir, 'Fig1.png'))


