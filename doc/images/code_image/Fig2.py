'''
   FIG. 2: scheme of the SKIM geometry with 4 beams at 12 degrees and 1 beam at 6 degree and 5 beams at 12 degrees and 2 beams at 6 degree.
'''

import netCDF4
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import skimulator.rw_data as rw
import glob
import cartopy
import os

# Initialize color
listcolor = ['c', 'y', 'b', 'g', 'k', 'r', 'c', 'y']

# List files
config="WW3_AT_metop_2018_8b"
indatadir = '/mnt/data/project/'
indatadir = '/tmp/key/data/skim_at_output/{}'.format(config)
#indatadir = os.path.join(indatadir, 'skim', 'skim_output')
filesgrid = os.path.join(indatadir, '{}_grid'.format(config))
ipass = 58
indatapath = '{}_p{:03d}.nc'.format(filesgrid, ipass)
outdatadir = '../'
modelbox = [-90, -70., 32., 40.]
#modelbox = [270, 290., 32., 40.]

# Prepare figure
pyplot.figure(figsize=(10, 15))
projection = cartopy.crs.Mercator()
transform = cartopy.crs.PlateCarree()
projection = transform
ax1 = pyplot.subplot(121, projection=projection)
#ax.coastlines()
ax1.add_feature(cartopy.feature.OCEAN, zorder=1)
ax1.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
#ax.set_extent([modelbox[0], modelbox[1],  modelbox[2], modelbox[3]],
#              projection)
ax1.set_extent([-74., -68., 34, 37], projection)
gl = ax1.gridlines(crs=transform, draw_labels=True, color='gray',
             linestyle='--', alpha=0.5)
gl.xlabels_top = False
gl.ylabels_left = False
ax1.set_title('(a) 6 beams configuration')

data = netCDF4.Dataset(indatapath, 'r')
lon = data['lon'][:]
lon[lon > 180] = lon[lon > 180] - 360
lat = data['lat'][:]
lon_nadir = data['lon_nadir'][:]
lon_nadir[lon_nadir > 180] = lon_nadir[lon_nadir > 180] - 360
lat_nadir = data['lat_nadir'][:]
pyplot.plot(lon_nadir, lat_nadir, 'k+', transform=transform)
for i in range(numpy.shape(lon)[1]):
    style_color = '{}+'.format(listcolor[i])
    pyplot.plot(lon[:, i], lat[:, i], style_color,
                transform=cartopy.crs.PlateCarree())
ax2 = pyplot.subplot(122,  projection=projection)
#ax.coastlines()
ax2.add_feature(cartopy.feature.OCEAN, zorder=1)
ax2.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
#ax.set_extent([modelbox[0], modelbox[1],  modelbox[2], modelbox[3]],
#              projection)
#config="WW3_GS_8b105az"
config="WW3_AT_metop_2018_6a"
indatadir = '/tmp/key/data/skim_at_output/{}'.format(config)
filesgrid = os.path.join(indatadir, '{}_grid'.format(config))
indatapath = '{}_p{:03d}.nc'.format(filesgrid, ipass)
print(indatapath)
ax2.set_extent([-74., -68., 34, 37], projection)
gl = ax2.gridlines(crs=transform, draw_labels=True, color='gray',
             linestyle='--', alpha=0.5)
gl.xlabels_top = False
gl.ylabels_left = False
ax2.set_title('(b) 8 beams configuration')
data = netCDF4.Dataset(indatapath, 'r')
lon = data['lon'][:]
lon[lon > 180] = lon[lon > 180] - 360
lat = data['lat'][:]
lon_nadir = data['lon_nadir'][:]
lon_nadir[lon_nadir > 180] = lon_nadir[lon_nadir > 180] - 360
lat_nadir = data['lat_nadir'][:]
pyplot.plot(lon_nadir, lat_nadir, 'k+', transform=transform)
for i in range(numpy.shape(lon)[1]):
    style_color = '{}+'.format(listcolor[i])
    pyplot.plot(lon[:, i], lat[:, i], style_color,
                transform=cartopy.crs.PlateCarree())
pyplot.savefig(os.path.join(outdatadir, 'Fig2.png'))
