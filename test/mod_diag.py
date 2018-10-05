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


def diag_rms(listfile, modelbox, output):
    nanmean12 = 0
    nanmean06 = 0
    nanstd12 = 0
    nanstd06 = 0
    nanstd12norm = 0
    nanstd06norm = 0
    i06 = 0
    i12 = 0
    for ifile in listfile:
        data = netCDF4.Dataset(ifile, 'r')
        lon = data['lon'][:]
        lat = data['lat'][:]
        lon = numpy.mod(lon + 180.0, 360.0) - 180.0
        lon_nadir = data['lon_nadir'][:]
        lon_nadir = numpy.mod(lon_nadir + 180.0, 360) - 180.0
        lat_nadir = data['lat_nadir'][:]
        vartrue = data['ur_true'][:]
        varobs = data['ur_obs'][:]
        mask = ((abs(vartrue)>10) | (abs(varobs)>10))
        vartrue[mask] = numpy.nan
        varobs[mask] = numpy.nan
        beam_angle = p.list_angle #data['beam_angle'][:]
        data.close()
        for i in range(numpy.shape(lon)[1]):
            _tmp = numpy.nanmean(((varobs[:, i] - vartrue[:,i])/vartrue[:,i])**2)
            _tmp = numpy.nanmean(((varobs[:, i] - vartrue[:,i]))**2)
            _tmp2 =  numpy.nanmean(abs(vartrue[:, i]))
            if beam_angle[i] == 6:
                nanstd06norm += _tmp/_tmp2
                nanstd06 += _tmp
                i06 += 1
            else:
                nanstd12 += _tmp
                nanstd12norm += _tmp/_tmp2
                i12 += 1
    nanstd06 = numpy.sqrt(nanstd06/i06)
    nanstd12 = numpy.sqrt(nanstd12/i12)
    nanstd06 = numpy.sqrt(nanstd06/i06)
    nanstd12 = numpy.sqrt(nanstd12/i12)
    return nanstd06, nanstd12, nanstd06norm, nanstd12norm

