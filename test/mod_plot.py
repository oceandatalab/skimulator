import numpy
import netCDF4
import os
import sys
import glob
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import cartopy
import skimulator.rw_data as rw
import params as p

# Initialize color
listcolor = ['c', 'y', 'b', 'g', 'k', 'r', 'c', 'y']

projection = cartopy.crs.PlateCarree()
transform = cartopy.crs.PlateCarree()


def plot_all(listfile, modelbox, output):
    # Initialize plot
    projection = cartopy.crs.PlateCarree()
    transform = cartopy.crs.PlateCarree()
    pyplot.figure(figsize=(12,12))
    if modelbox is None:
        projection = cartopy.crs.Orthographic(0, 0)
    ax = pyplot.axes(projection=projection)
    if modelbox is not None:
        ax.set_extent([modelbox[0], modelbox[1],  modelbox[2], modelbox[3]],
                      crs=transform)
        norder = 6
    else:
        ax.set_global()
        norder = 1
    #ax.add_feature(cartopy.feature.OCEAN, zorder=norder)
    #ax.add_feature(cartopy.feature.LAND, zorder=norder, edgecolor='black')
    gl = ax.gridlines(crs=transform, draw_labels=True, color='gray',
                       linestyle='--', alpha=0.5)
    gl.xlabels_top = False
    gl.ylabels_left = False
    # Loop on files
    for ifile in listfile:
        data = netCDF4.Dataset(ifile, 'r')
        lon = data['lon'][:]
        lat = data['lat'][:]
        lon = numpy.mod(lon + 180.0, 360.0) - 180.0
        lon_nadir = data['lon_nadir'][:]
        lon_nadir = numpy.mod(lon_nadir + 180.0, 360) - 180.0
        lat_nadir = data['lat_nadir'][:]
        pyplot.plot(lon[:, 0], lat[:, 0], '.', color=listcolor[0], markersize=1,
                    transform=transform)
        pyplot.plot(lon_nadir, lat_nadir, '.', color=listcolor[4], markersize=1,
                    transform=transform)
    # Save figure
    pyplot.savefig(output)
    return None