import numpy
import netCDF4
import os
import sys
import glob
import pickle
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


def bin_variables(listfile, listvar, bin_in, modelbox):
    lonp = numpy.arange(modelbox[0], modelbox[1] + modelbox[2], modelbox[2])
    latp = numpy.arange(modelbox[3], modelbox[4] + modelbox[5], modelbox[5])
    resol = numpy.sqrt((modelbox[2] * numpy.cos(numpy.deg2rad(latp))) **2
                        + modelbox[5]**2)
    dic_v = {}
    for ifile in listfile:
        data = netCDF4.Dataset(ifile, 'r')
        lon = data['lon'][:]
        lat = data['lat'][:]
        for j in range(len(lonp)):
            for i in range(len(latp)):
                ind_key = 10000 * int(i) + int(j)
                lon = numpy.mod(lon + 180, 360) - 180
                dist = numpy.sqrt(((lonp[j] - lon)
                                   *numpy.cos(numpy.deg2rad(lat)))**2
                                   + (latp[i] - lat)**2)
                iiobs = numpy.where(dist < resol[i])
                if not iiobs:
                    continue
                if ind_key not in dic_v.keys():
                    dic_v[ind_key] = {}
                for ivar in listvar:
                    var = data[ivar][:]
                    mask = var.mask[iiobs]
                    var = numpy.array(var.data[iiobs])
                    if ivar not in dic_v[ind_key].keys():
                        dic_v[ind_key][ivar] = []
                    if var.any():
                        dic_v[ind_key][ivar].append(var[~mask])
        data.close()
    with open(bin_in, 'wb') as f:
        pickle.dump(dic_v, f)


def compute_rms(bin_in, bin_out, listvar, modelbox):
    with open(bin_in, 'rb') as f:
        dic_v = pickle.load(f)
    lonp = numpy.arange(modelbox[0], modelbox[1] + modelbox[2], modelbox[2])
    latp = numpy.arange(modelbox[3], modelbox[4] + modelbox[5], modelbox[5])
    rms = {}
    std = {}
    snr = {}
    rms['lon'] = lonp
    rms['lat'] = latp
    std['lon'] = lonp
    std['lat'] = latp
    snr['lon'] = lonp
    snr['lat'] = latp
    for ivar in listvar:
        rms[ivar] = numpy.full((len(latp), len(lonp)), numpy.nan)
        std[ivar] = numpy.full((len(latp), len(lonp)), numpy.nan)
        snr[ivar] = numpy.full((len(latp), len(lonp)), numpy.nan)
    for j in range(len(rms['lon'])):
        for i in range(len(rms['lat'])):
            ind_key = 10000 * int(i) + int(j)
            if ind_key  not in dic_v.keys():
                continue
            var1 = numpy.array(dic_v[ind_key]['ur_true']).flatten()
            print(var1)
            print(numpy.nanmean(var1))
            rms_true = numpy.sqrt(numpy.nanmean(var1**2))
            for ivar in listvar:
                var = numpy.array(dic_v[ind_key][ivar]).flatten()
                rms[ivar] = numpy.sqrt(numpy.nanmean(var**2))
                std[ivar] = numpy.nanstd(var)
                if ivar == 'ur_obs':
                    _std = numpy.nanstd(var - var1)
                    snr[ivar] = rms_true / _std
                elif ivar == 'instr' or ivar == 'uwb_corr':
                    snr[ivar]  = rms_true / std[ivar]
    with open('rms_{}'.format(bin_out), 'wb') as f:
        pickle.dump(rms, f)

    with open('std_{}'.format(bin_out), 'wb') as f:
        pickle.dump(std, f)

    with open('snr_{}'.format(bin_out), 'wb') as f:
        pickle.dump(snr, f)


def plot_rms(pfile, list_var, outfile, rms=True, std=True, snr=True):
    import mod_plot
    if rms is True:
        with open('rms_{}'.format(pfile), 'rb') as f:
            rms = pickle.load(f)
    if std is True:
        with open('std_{}'.format(pfile), 'rb') as f:
            std = pickle.load(f)
    if snr is True:
        with open('snr_{}'.format(pfile), 'rb') as f:
            snr = pickle.load(f)

    for ivar in list_var:
        if ivar in rms.keys() and rms is True:
            lon = rms['lon']
            lat = rms['lat']
            var = rms[ivar]
            _outfile = 'rms_{}_{}.png'.format(ivar, outfile)
            mod_plot.plot_diag(lon, lat, var, _outfile)
        if ivar in std.keys() and std is True:
            lon = std['lon']
            lat = std['lat']
            var = std[ivar]
            _outfile = 'std_{}_{}.png'.format(ivar, outfile)
            mod_plot.plot_diag(lon, lat, var, _outfile)
        if ivar in snr.keys():
            lon = snr['lon']
            lat = snr['lat']
            var = snr[ivar]
            _outfile = 'snr_{}_{}.png'.format(ivar, outfile)
            mod_plot.plot_diag(lon, lat, var, _outfile)
    return None


