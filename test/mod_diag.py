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

# Initialize color
listcolor = ['c', 'y', 'b', 'g', 'k', 'r', 'c', 'y']

projection = cartopy.crs.PlateCarree()
transform = cartopy.crs.PlateCarree()


def diag_rms(listfile, modelbox, output, list_angle):
    nanmean12 = 0
    nanmean06 = 0
    nanmeaninstr12 = 0
    nanmeaninstr06 = 0
    nanmeanuwb12 = 0
    nanmeanuwb06 = 0
    nanmeanuwbc12 = 0
    nanmeanuwbc06 = 0
    nanstd12 = 0
    nanstd06 = 0
    nanstd12norm = 0
    nanstd06norm = 0
    i06 = 0
    i12 = 0
    for ifile in listfile:
        print(ifile)
        data = netCDF4.Dataset(ifile, 'r')
        lon = data.variables['lon'][:]
        lat = data.variables['lat'][:]
        lon = numpy.mod(lon + 180.0, 360.0) - 180.0
        lon_nadir = data.variables['lon_nadir'][:]
        lon_nadir = numpy.mod(lon_nadir + 180.0, 360) - 180.0
        lat_nadir = data.variables['lat_nadir'][:]
        vartrue = data.variables['ur_true'][:]
        varobs = data.variables['ur_obs'][:]
        if 'instr' in data.variables.keys():
            instr = True
            varinstr = data.variables['instr'][:]
        if 'uwb' in data.variables.keys():
            uwb = True
            varuwb = data.variables['uwb'][:]
        if 'uwb_corr' in data.variables.keys():
            uwbc = True
            varuwbc = data.variables['uwb_corr'][:]

        mask = ((abs(vartrue)>100) | (abs(varobs)>100))
        #vartrue[mask] = numpy.nan
        #varobs[mask] = numpy.nan
        beam_angle = list_angle #data['beam_angle'][:]
        data.close()
        if numpy.shape(lon)[0]<100:
            continue
        for i in range(numpy.shape(lon)[1]):
            if varobs[:, i].mask.all():
                continue
            if varuwbc[:, i].mask.all():
                continue
            try:
                _tmp = numpy.nanmean(((varobs[:, i] - vartrue[:,i])/vartrue[:,i])**2)
                _tmp = numpy.nanmean(((varobs[:, i] - vartrue[:,i]))**2)
                _tmp2 =  numpy.nanmean(abs(vartrue[:, i])**2)
            except:
                continue
            if instr is True:
                _tmpinstr = numpy.nanmean(abs(varinstr[:, i])**2)
            if uwb is True:
                try:
                    _tmpuwb = numpy.nanmean(abs(varuwb[:, i])**2)
                except:
                    continue
            if uwbc is True:
                _tmpuwbc = numpy.nanmean(abs(varuwbc[:, i])**2)
            if beam_angle[i] == 6:
                nanstd06norm += _tmp/_tmp2
                nanstd06 += _tmp
                if instr is True:
                    nanmeaninstr06 += _tmpinstr
                if uwb is True:
                    nanmeanuwb06 += _tmpuwb
                if uwbc is True:
                    nanmeanuwbc06 += _tmpuwbc
                i06 += 1
            else:
                nanstd12 += _tmp
                nanstd12norm += _tmp/_tmp2
                if instr is True:
                    nanmeaninstr12 += _tmpinstr
                if uwb is True:
                    nanmeanuwb12 += _tmpuwb
                if uwbc is True:
                    nanmeanuwbc12 += _tmpuwbc
                i12 += 1
    nanstd06 = numpy.sqrt(nanstd06/i06)
    nanstd12 = numpy.sqrt(nanstd12/i12)
    # nanstd06 = numpy.sqrt(nanstd06/i06)
    # nanstd12 = numpy.sqrt(nanstd12/i12)
    if instr is True:
        nanmeaninstr06 = numpy.sqrt(nanmeaninstr06/i06)
        nanmeaninstr12 = numpy.sqrt(nanmeaninstr12/i12)
    if uwb is True:
        nanmeanuwb06 = numpy.sqrt(nanmeanuwb06/i06)
        nanmeanuwb12 = numpy.sqrt(nanmeanuwb12/i12)
    if uwbc is True:
        nanmeanuwbc06 = numpy.sqrt(nanmeanuwbc06/i06)
        nanmeanuwbc12 = numpy.sqrt(nanmeanuwbc12/i12)
    print('RMSE 06: {}, RMSE 12: {}'.format(nanstd06, nanstd12))
    print('NRMSE 06: {}, NRMSE 12: {}'.format(nanstd06norm, nanstd12norm))
    if instr is True:
        print('instr rms 06: {}, instr rms 12: {}'.format(nanmeaninstr06, nanmeaninstr12))
    if uwb is True:
        print('uwb rms 06: {}, uwb rms 12: {}'.format(nanmeanuwb06, nanmeanuwb12))
    if uwbc is True:
        print('uwbc rms 06: {}, uwbc rms 12: {}'.format(nanmeanuwbc06, nanmeanuwbc12))
    return nanstd06, nanstd12, nanstd06norm, nanstd12norm


def diag_azimuth_rms(listfile, modelbox, output, list_angle):
    nanmean12 = 0
    nanmean06 = 0
    nanmeaninstr12 = 0
    nanmeaninstr06 = 0
    nanmeanuwb12 = 0
    nanmeanuwb06 = 0
    nanmeanuwbc12 = 0
    nanmeanuwbc06 = 0
    nanstd12 = 0
    nanstd06 = 0
    nanstd12norm = 0
    nanstd06norm = 0
    i06 = 0
    i12 = 0
    for ifile in listfile:
        print(ifile)
        data = netCDF4.Dataset(ifile, 'r')
        lon = data.variables['lon'][:]
        lat = data.variables['lat'][:]
        lon = numpy.mod(lon + 180.0, 360.0) - 180.0
        lon_nadir = data.variables['lon_nadir'][:]
        lon_nadir = numpy.mod(lon_nadir + 180.0, 360) - 180.0
        lat_nadir = data.variables['lat_nadir'][:]
        vartrue = data.variables['ur_true'][:]
        varobs = data.variables['ur_obs'][:]
        if 'instr' in data.variables.keys():
            instr = True
            varinstr = data.variables['instr'][:]
        if 'uwb' in data.variables.keys():
            uwb = True
            varuwb = data.variables['uwb'][:]
        if 'uwb_corr' in data.variables.keys():
            uwbc = True
            varuwbc = data.variables['uwb_corr'][:]

        mask = ((abs(vartrue)>100) | (abs(varobs)>100))
        #vartrue[mask] = numpy.nan
        #varobs[mask] = numpy.nan
        beam_angle = list_angle #data['beam_angle'][:]
        data.close()
        if numpy.shape(lon)[0]<100:
            continue
        for i in range(numpy.shape(lon)[1]):
            if varobs[:, i].mask.all():
                continue
            if varuwbc[:, i].mask.all():
                continue
            try:
                _tmp = numpy.nanmean(((varobs[:, i] - vartrue[:,i])/vartrue[:,i])**2)
                _tmp = numpy.nanmean(((varobs[:, i] - vartrue[:,i]))**2)
                _tmp2 =  numpy.nanmean(abs(vartrue[:, i])**2)
            except:
                continue
            if instr is True:
                _tmpinstr = numpy.nanmean(abs(varinstr[:, i])**2)
            if uwb is True:
                try:
                    _tmpuwb = numpy.nanmean(abs(varuwb[:, i])**2)
                except:
                    continue
            if uwbc is True:
                _tmpuwbc = numpy.nanmean(abs(varuwbc[:, i])**2)
            if beam_angle[i] == 6:
                nanstd06norm += _tmp/_tmp2
                nanstd06 += _tmp
                if instr is True:
                    nanmeaninstr06 += _tmpinstr
                if uwb is True:
                    nanmeanuwb06 += _tmpuwb
                if uwbc is True:
                    nanmeanuwbc06 += _tmpuwbc
                i06 += 1
            else:
                nanstd12 += _tmp
                nanstd12norm += _tmp/_tmp2
                if instr is True:
                    nanmeaninstr12 += _tmpinstr
                if uwb is True:
                    nanmeanuwb12 += _tmpuwb
                if uwbc is True:
                    nanmeanuwbc12 += _tmpuwbc
                i12 += 1
    nanstd06 = numpy.sqrt(nanstd06/i06)
    nanstd12 = numpy.sqrt(nanstd12/i12)
    # nanstd06 = numpy.sqrt(nanstd06/i06)
    # nanstd12 = numpy.sqrt(nanstd12/i12)
    if instr is True:
        nanmeaninstr06 = numpy.sqrt(nanmeaninstr06/i06)
        nanmeaninstr12 = numpy.sqrt(nanmeaninstr12/i12)
    if uwb is True:
        nanmeanuwb06 = numpy.sqrt(nanmeanuwb06/i06)
        nanmeanuwb12 = numpy.sqrt(nanmeanuwb12/i12)
    if uwbc is True:
        nanmeanuwbc06 = numpy.sqrt(nanmeanuwbc06/i06)
        nanmeanuwbc12 = numpy.sqrt(nanmeanuwbc12/i12)
    print('RMSE 06: {}, RMSE 12: {}'.format(nanstd06, nanstd12))
    print('NRMSE 06: {}, NRMSE 12: {}'.format(nanstd06norm, nanstd12norm))
    if instr is True:
        print('instr rms 06: {}, instr rms 12: {}'.format(nanmeaninstr06, nanmeaninstr12))
    if uwb is True:
        print('uwb rms 06: {}, uwb rms 12: {}'.format(nanmeanuwb06, nanmeanuwb12))
    if uwbc is True:
        print('uwbc rms 06: {}, uwbc rms 12: {}'.format(nanmeanuwbc06, nanmeanuwbc12))
    return nanstd06, nanstd12, nanstd06norm, nanstd12norm


def bin_rms(listfile, listvar, bin_in, modelbox):
    lonp = numpy.arange(modelbox[0], modelbox[1] + modelbox[2], modelbox[2])
    latp = numpy.arange(modelbox[3], modelbox[4] + modelbox[5], modelbox[5])
    resol = numpy.sqrt((modelbox[2] * numpy.cos(numpy.deg2rad(latp)) **2 + modelbox[5]**2)
    for ifile in listfile:
        data = netCDF4.Dataset(ifile, 'r')
        lon = data['lon'][:]
        lat = data['lat'][:]
        for j in len(lon):
            for i in len(lat):
                ind_key = 10000 * int(i) + int(j)
                dist = numpy.sqrt(((lonp[i, j] - lon)
                                   *numpy.cos(numpy.deg2rad(lat))**2
                                   + (latp[i, j] - lat)**2)
                iiobs = numpy.where(dist < resol[i])
                if iiobs.isempty():
                    continue
                if ind_key not in dic_v.keys()
                    dic_v[ind_j] = {}
                for ivar in listvar:
                    var = data[ivar][:]
                    if ivar not in dic_v[ind_key].keys():
                        dic_v[ind_key][ivar] = []
                    dic_v[ind_key][ivar].append(var[iiobs]
        data.close()
    with open(bin_in, 'wb') as f:
        pickle.dump(dic_v, f)


def read_bin(bin_in, list_var, list_diag, modelbox)
    with open(bin_in, 'rb') as f:
        dic_v = pickle.load(f)
    lonp = numpy.arange(modelbox[0], modelbox[1] + modelbox[2], modelbox[2])
    latp = numpy.arange(modelbox[3], modelbox[4] + modelbox[5], modelbox[5])
    for var in list_var:
        rms[var] = numpy.full((len{latp), len(lonp)), numpy.nan)
        std[var] = numpy.full((len{latp), len(lonp)), numpy.nan)
        snr[var] = numpy.full((len{latp), len(lonp)), numpy.nan)
    for j in len(lon):
        for i in len(lat):
            ind_key = 10000 * int(i) + int(j)
            for ivar in listvar:
                var = dic_v[ind_key][ivar]
                rms[var] = numpy.nanrms(var)
                std[var] = numpy.nanstd(var)
                std[var] = numpy.nanstd(var)

