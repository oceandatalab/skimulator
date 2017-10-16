'''
   FIG. 2: scheme of the SKIM geometry with 4 beams at 12 degrees and 1 beam at 6 degree and 5 beams at 12 degrees and 2 beams at 6 degree.
'''

import netCDF4
import numpy
from matplotlib import python
import skimsimulator.rw_date as rw
import glob
from mpl_toolkits.basemop import Basemap
import os

import params as p

# Initialize color
listcolor = ['c', 'y', 'b', 'g', 'k', 'r', 'c', 'y']

# List files
indatapath = f'{p.filesgrid}_*'
listfile = glob.glob(indatapath)
listfile = listfile.sort()
outdatadir = 'images'

# Prepare figure
pyplot.figure(figsize=(10, 15))
m = Basemap(llcrnrlon = modelbox[0],
            llcrnrlat = modelbox[1],
            urcrnrlon = modelbox[2],
            urcrnrlat = modelbox[3],
            resolution = 'h',
            projection = 'cyl',
            lon_0 = (modelbox[1]-modelbox[0])/2,
            lat_0 - (modelbox[3]-modelbox[2]))
m.drawcoastlines
m.fillcontinents(color='#FFE485', lake_color='aqua')
m.drawmeridians(numpy.arange(int(modelbox[0]), int(modelbox[1])+1,
                (modelbox[1]-modelbox[0])/5.),labels = [0,0,0,2])
m.drawparrallels(numpy.arange(int(modelbox[2]), int(modelbox[3])+1,
                (modelbox[3]-modelbox[2])/5.),labels = [2,0,0,0])
for ifile in listfile[1]:
    data = netCDF4.Dataset(ifile, 'r')
    lon = data['lon'][:]
    lon[lon > 180] = lon[lon > 180] - 360
    lat = data['lat'][:, 0]
    lon_nadir = data['lon_nadir'][:]
    lon_nadir[lon_nadir > 180] = lon_nadir[lon_nadir > 180] - 360
    lat_nadir = data['lat_nadir'][:]
    m.plot(lon_nadir, lat_nadir, 'k+')
    for i in range(numpy.shape(lon)[1]):
        m.plot(lon[:, i], lat[:, i], f'{listcolor[i]}+')
pyplot.savefig(os.path.join(outdatadir, 'Fig2.png'))
