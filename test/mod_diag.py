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
    nanmeandsig12 = 0
    nanmeandsig06 = 0
    nanmeangsig12 = 0
    nanmeangsig06 = 0
    nanstd12 = 0
    nanstd06 = 0
    nanstd12norm = 0
    nanstd06norm = 0
    i06 = 0
    i12 = 0
    for ifile in listfile:
        print(ifile)
        data = netCDF4.Dataset(ifile, 'r')
        ind0 = 100
        ind1 = -100
        sig = data.variables['sigma0'][ind0:ind1, :]
        uwnd = data.variables['uwnd'][ind0:ind1, :]
        vwnd = data.variables['vwnd'][ind0:ind1, :]
        wnd = numpy.sqrt(uwnd**2 + vwnd**2)
        lon = data.variables['lon'][ind0:ind1, :]
        lat = data.variables['lat'][ind0:ind1, :]
        lon = numpy.mod(lon + 180.0, 360.0) - 180.0
        lon_nadir = data.variables['lon_nadir'][ind0:ind1]
        lon_nadir = numpy.mod(lon_nadir + 180.0, 360) - 180.0
        lat_nadir = data.variables['lat_nadir'][ind0:ind1]
        vartrue = data.variables['ur_true'][ind0:ind1, :]
        varobs = data.variables['ur_obs'][ind0:ind1, :]
        if len(lon_nadir) < 700:
            continue
        instr = False
        if 'instr' in data.variables.keys():
            instr = True
            varinstr = data.variables['instr'][ind0:ind1, :]
            varinstr[numpy.where(abs(varinstr) > 1)] = numpy.nan
        uwb = False
        if 'uwd' in data.variables.keys():
            uwb = True
            varuwb = data.variables['uwd'][ind0:ind1, :]
        uwbc = False
        if 'uwd_est' in data.variables.keys():
            uwbc = True
            varuwbc = data.variables['uwd_est'][ind0:ind1, :]
            varuwbc = varuwbc - varuwb
        dsigma = False
        if 'dsigma' in data.variables.keys():
            dsigma = True
            vardsig = data.variables['dsigma'][ind0:ind1, :]
        gsig_atm = False

        if 'gsig_atm_err' in data.variables.keys():
            gsig_atm = False
            vargsig = data.variables['gsig_atm_err'][ind0:ind1, :]
        mask = ((abs(vartrue)>100) | (abs(varobs)>100) | numpy.isnan(varinstr))
        vartrue[mask] = numpy.nan
        varobs[mask] = numpy.nan
        beam_angle = list_angle #data['beam_angle'][:]
        data.close()
        # Remove borders
        if numpy.shape(lon)[0]<100:
            continue
        for i in range(numpy.shape(lon)[1]):
            if varobs[:, i].mask.all():
                continue
            if varuwbc[:, i].mask.all():
                continue
            _ind = numpy.where((wnd[:, i]>7) & (wnd[:, i]<10))
            #_ind = numpy.where((sig[:, i]>10**-5) & (wnd[:, i]>3))
            try:
                _tmp = numpy.nanmean(((varobs[_ind, i] - vartrue[_ind,i])/vartrue[_ind,i])**2)
                _tmp = numpy.nanmean(((varobs[_ind, i] - vartrue[_ind,i]))**2)
                _tmp2 =  numpy.nanmean(abs(vartrue[_ind, i])**2)
            except:
                continue
            if instr is True:
                _tmpinstr = numpy.nanmean(abs(varinstr[_ind, i])**2)
            if uwb is True:
                try:
                    _tmpuwb = numpy.nanmean(abs(varuwb[_ind, i])**2)
                except:
                    continue
            if uwbc is True:
                _tmpuwbc = numpy.nanmean(abs(varuwbc[_ind, i])**2)
            if dsigma is True:
                _tmpdsig = numpy.nanmean(abs(vardsig[_ind, i])**2)
            if gsig_atm is True:
                _tmpgsig = numpy.nanmean(abs(vargsig[_ind, i])**2)
            if beam_angle[i] == 6:
                nanstd06norm += _tmp2
                nanstd06 += _tmp
                if instr is True:
                    nanmeaninstr06 += _tmpinstr
                if uwb is True:
                    nanmeanuwb06 += _tmpuwb
                if uwbc is True:
                    nanmeanuwbc06 += _tmpuwbc
                if dsigma is True:
                    nanmeandsig06 += _tmpdsig
                if gsig_atm is True:
                    nanmeangsig06 += _tmpgsig
                i06 += 1
            else:
                nanstd12 += _tmp
                nanstd12norm += _tmp2
                if instr is True:
                    nanmeaninstr12 += _tmpinstr
                if uwb is True:
                    nanmeanuwb12 += _tmpuwb
                if uwbc is True:
                    nanmeanuwbc12 += _tmpuwbc
                if dsigma is True:
                    nanmeandsig12 += _tmpdsig
                if gsig_atm is True:
                    nanmeangsig12 += _tmpgsig
                i12 += 1
    nanstd06 = numpy.sqrt(nanstd06/i06)
    nanstd12 = numpy.sqrt(nanstd12/i12)
    # nanstd06 = numpy.sqrt(nanstd06/i06)
    # nanstd12 = numpy.sqrt(nanstd12/i12)
    nanstd06norm = numpy.sqrt(nanstd06norm/i06)
    nanstd12norm = numpy.sqrt(nanstd12norm/i12)
    if instr is True:
        nanmeaninstr06 = numpy.sqrt(nanmeaninstr06/i06)
        nanmeaninstr12 = numpy.sqrt(nanmeaninstr12/i12)
    if uwb is True:
        nanmeanuwb06 = numpy.sqrt(nanmeanuwb06/i06)
        nanmeanuwb12 = numpy.sqrt(nanmeanuwb12/i12)
    if uwbc is True:
        nanmeanuwbc06 = numpy.sqrt(nanmeanuwbc06/i06)
        nanmeanuwbc12 = numpy.sqrt(nanmeanuwbc12/i12)
    if dsigma is True:
        nanmeandsig06 = numpy.sqrt(nanmeandsig06/i06)
        nanmeandsig12 = numpy.sqrt(nanmeandsig12/i12)
    if gsig_atm is True:
        nanmeangsig06 = numpy.sqrt(nanmeangsig06/i06)
        nanmeangsig12 = numpy.sqrt(nanmeangsig12/i12)
    print('RMSE 06: {}, RMSE 12: {}'.format(nanstd06, nanstd12))
    print('NRMSE 06: {}, NRMSE 12: {}'.format(nanstd06norm, nanstd12norm))
    if instr is True:
        print('instr rms 06: {}, instr rms 12: {}'.format(nanmeaninstr06, nanmeaninstr12))
    if uwb is True:
        print('uwb rms 06: {}, uwb rms 12: {}'.format(nanmeanuwb06, nanmeanuwb12))
    if uwbc is True:
        print('uwbc rms 06: {}, uwbc rms 12: {}'.format(nanmeanuwbc06, nanmeanuwbc12))
    if dsigma is True:
        print('dsigma rms 06: {}, dsigma rms 12: {}'.format(nanmeandsig06, nanmeandsig12))
    if gsig_atm is True:

        print('gsig rms 06: {}, gsig rms 12: {}'.format(nanmeangsig06, nanmeangsig12))
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
        if 'uwd' in data.variables.keys():
            uwb = True
            varuwb = data.variables['uwd'][:]
        if 'uwd_est' in data.variables.keys():
            uwbc = True
            varuwbc = data.variables['uwd_est'][:]

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
                    try:
                        var = data[ivar][:]
                    except: print(ifile, ivar)
                    _mask = numpy.ma.getmaskarray(var)
                    mask = _mask[iiobs]
                    if ivar == 'ur_obs':
                        ivar2 = 'diff_ur'
                        var1 = data['ur_true'][:]
                        var2 = numpy.array(abs(var.data[iiobs] - var1.data[iiobs]))
                        if var2.any():
                            if ivar2 not in dic_v[ind_key].keys():
                                dic_v[ind_key][ivar2] = []
                            dic_v[ind_key][ivar2].append(var2[~mask])
                        ivar2 = 'cov_ur'
                        cov = numpy.cov(var.data[iiobs], var1.data[iiobs])
                        var2 = (200*cov[0, 1] / (numpy.trace(cov)))
                        if var2.any():
                            if ivar2 not in dic_v[ind_key].keys():
                                dic_v[ind_key][ivar2] = []
                            dic_v[ind_key][ivar2].append(var2)
                    if ivar == 'uwd_est':
                        ivar2 = 'diff_uwdr'
                        var1 = data['uwd'][:]
                        var2 = numpy.array(abs(var.data[iiobs] - var1.data[iiobs]))
                        if var2.any():
                            if ivar2 not in dic_v[ind_key].keys():
                                dic_v[ind_key][ivar2] = []
                            dic_v[ind_key][ivar2].append(var2[~mask])
                        ivar2 = 'cov_uwdr'
                        cov = numpy.cov(var.data[iiobs], var1.data[iiobs])
                        var2 = (200*cov[0, 1] / (numpy.trace(cov)))
                        if var2.any():
                            if ivar2 not in dic_v[ind_key].keys():
                                dic_v[ind_key][ivar2] = []
                            dic_v[ind_key][ivar2].append(var2)
                    var = numpy.array(var.data[iiobs])
                    if var.any():
                        if ivar not in dic_v[ind_key].keys():
                            dic_v[ind_key][ivar] = []
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
    cov = {}
    rms['lon'] = lonp
    rms['lat'] = latp
    std['lon'] = lonp
    std['lat'] = latp
    snr['lon'] = lonp
    snr['lat'] = latp
    if 'uwd_est' in listvar:
        listvar.append('diff_uwdr')
    listvar.append('diff_ur')
    listvar.append('cov_ur')
    listvar.append('cov_uwdr')
    for ivar in listvar:
        rms[ivar] = numpy.full((len(latp), len(lonp)), numpy.nan)
        std[ivar] = numpy.full((len(latp), len(lonp)), numpy.nan)
        snr[ivar] = numpy.full((len(latp), len(lonp)), numpy.nan)
        cov[ivar] = numpy.full((len(latp), len(lonp)), numpy.nan)
    for j in range(len(rms['lon'])):
        for i in range(len(rms['lat'])):
            ind_key = 10000 * int(i) + int(j)
            if ind_key  not in dic_v.keys():
                continue
            if 'ur_true' not in dic_v[ind_key].keys():
                continue
            if 'uwd' not in dic_v[ind_key].keys():
                continue
            var1 = numpy.concatenate(dic_v[ind_key]['ur_true']).ravel()
            if 'uwd' in dic_v[ind_key].keys():
                var2 = numpy.concatenate(dic_v[ind_key]['uwd']).ravel()
            if not var1.any() or not var2.any():
                continue
            rms_true = numpy.sqrt(numpy.nanmean(var1**2))
            # import pdb ; pdb.set_trace()
            for ivar in listvar:
                if ivar not in dic_v[ind_key].keys():
                    continue
                try:
                    var = numpy.concatenate(dic_v[ind_key][ivar]).ravel()
                except:
                    continue
                rms[ivar][i, j] = numpy.sqrt(numpy.nanmean(var**2))
                std[ivar][i, j] = numpy.nanstd(var)
                if ivar == 'diff_ur':
                    snr[ivar][i, j] = rms_true / std[ivar][i, j]
                elif ivar == 'instr' or ivar == 'diff_uwdr':
                    snr[ivar][i, j]  = rms_true / std[ivar][i, j]
                elif 'cov' in ivar:
                    cov[ivar][i, j] = numpy.nanmean(var)
    with open('rms_{}'.format(bin_out), 'wb') as f:
        pickle.dump(rms, f)

    with open('std_{}'.format(bin_out), 'wb') as f:
        pickle.dump(std, f)

    with open('snr_{}'.format(bin_out), 'wb') as f:
        pickle.dump(snr, f)

    with open('cov_{}'.format(bin_out), 'wb') as f:
        pickle.dump(cov, f)

def plot_rms(pfile, list_var, outfile, isrms=True, isstd=True, issnr=True, 
             isvar=True):
    import mod_plot
    if isrms is True:
        with open('rms_{}'.format(pfile), 'rb') as f:
            rms = pickle.load(f)
    if isstd is True:
        with open('std_{}'.format(pfile), 'rb') as f:
            std = pickle.load(f)
    if issnr is True:
        with open('snr_{}'.format(pfile), 'rb') as f:
            snr = pickle.load(f)
    if isvar is True:
        with open('cov_{}'.format(pfile), 'rb') as f:
            cov = pickle.load(f)

    if 'uwd_est' in list_var:
        list_var.append('diff_uwdr')
    list_var.append('diff_ur')
    for ivar in list_var:
        if ivar in rms.keys() and isrms is True:
            lon = rms['lon']
            lat = rms['lat']
            var = rms[ivar]
            _outfile = 'rms_{}_{}.png'.format(ivar, outfile)
            vmin = numpy.nanpercentile(var, 5)
            vmax = numpy.nanpercentile(var, 95)
            mod_plot.plot_diag(lon, lat, var, _outfile, vmin=vmin, vmax=vmax, cmap='jet')
        if ivar in std.keys() and isstd is True:
            lon = std['lon']
            lat = std['lat']
            var = std[ivar]
            vmin = numpy.nanpercentile(var, 5)
            vmax = numpy.nanpercentile(var, 95)
            _outfile = 'std_{}_{}.png'.format(ivar, outfile)
            mod_plot.plot_diag(lon, lat, var, _outfile, vmin=vmin, vmax=vmax, cmap='jet')
        if ivar in snr.keys():
            lon = snr['lon']
            lat = snr['lat']
            var = snr[ivar]
            vmin = numpy.nanpercentile(var, 5)
            vmax = numpy.nanpercentile(var, 95)
            _outfile = 'snr_{}_{}.png'.format(ivar, outfile)
            mod_plot.plot_diag(lon, lat, var, _outfile, vmin=vmin, vmax=vmax, cmap='jet')
    for ivar in ['cov_uwdr', 'cov_ur']:
        if ivar in cov.keys():
            lon = cov['lon']
            lat = cov['lat']
            var = cov[ivar]
            vmin = numpy.nanpercentile(var, 5)
            vmax = numpy.nanpercentile(var, 95)
            _outfile = 'snr_{}_{}.png'.format(ivar, outfile)
            mod_plot.plot_diag(lon, lat, var, _outfile, vmin=vmin, vmax=vmax, cmap='jet')
    return None


