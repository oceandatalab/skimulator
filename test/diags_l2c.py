import numpy
from matplotlib import pyplot
import netCDF4
import glob
import os
import sys
import json
import scipy.signal
from scipy.fftpack import fft


def cpsd1d(hh1=None, hh2=None, dx=1.,tap=0.05, detrend=True):


    hh1 = hh1 - numpy.mean(hh1)
    hh2 = hh2 - numpy.mean(hh2)
    nx = numpy.shape(hh1)[0]


    if detrend:
        hh1 = scipy.signal.detrend(hh1)
        hh2 = scipy.signal.detrend(hh2)

    if tap>0:
        ntaper = numpy.int(tap * nx + 0.5)
        taper = numpy.zeros(nx)+1.
        taper[:ntaper] = numpy.cos(numpy.arange(ntaper) / (ntaper-1.)
                                   *numpy.pi/2 + 3*numpy.pi/2)
        taper[-ntaper:] = numpy.cos(-numpy.arange(-ntaper+1, 1)/(ntaper-1.)
                                    *numpy.pi/2 + 3*numpy.pi/2)
        hh1 = hh1 * taper
        hh2 = hh2 * taper

    ss1 = fft(hh1)
    ss2 = fft(hh2)

    ff = numpy.arange(1, nx/2 - 1) / (nx*dx)


    C = numpy.cos(numpy.angle(ss1[1: int(nx/2)-1])
                  - numpy.angle(ss2[1: int(nx/2)-1]))

    return ff, C


def coherency_l2c(datadir_input, config, var, nal_min,
                  posting, outfile='coherency', fsize=1):
    list_files_all = []
    pyplot.figure()
    for indir, iconfig, ivar in zip(datadir_input, config, var):
        print(indir, iconfig, ivar)
        list_files = glob.glob(os.path.join(indir, '{}_l2c_c*p*.nc'.format(iconfig)))

        fid = netCDF4.Dataset(list_files[0], 'r')
        try:
            _tmp = numpy.array(fid.variables['u_ac_true'][:])
        except:
            continue
        _tmp[_tmp < -10] = numpy.nan
        fid.close()
        nac = _tmp.shape[1]
        nac_mid = int(nac/2)
        nal = nac


        dal = posting #km

        MCuac = []
        countuac = 0
        MCual = []
        countual = 0
        for ifile in list_files:
            ref = {}
            skim = {}
            fid = netCDF4.Dataset(ifile, 'r')
            try:
                ref['uac'] = numpy.array(fid.variables['u_ac_true'][:])
                ref['uac'] = scipy.ndimage.uniform_filter(ref['uac'], size=fisze)
            except:
                print(ifile)
                continue
            ref['uac'][ref['uac']<-10] = numpy.nan
            ref['uac'] = numpy.ma.masked_invalid(ref['uac'])
            ref['ual'] = numpy.array(fid.variables['u_al_true'][:])
            ref['ual'] = scipy.ndimage.uniform_filter(ref['ual'], size=fsize)
            ref['ual'][ref['ual']<-10] = numpy.nan
            ref['ual'] = numpy.ma.masked_invalid(ref['ual'])
            skim['uac'] = numpy.array(fid.variables['u_ac_{}'.format(ivar)][:])
            skim['uac'][numpy.abs(skim['uac']) > 10] = numpy.nan
            skim['uac'] = numpy.ma.masked_invalid(skim['uac'])
            skim['ual'] = numpy.array(fid.variables['u_al_{}'.format(ivar)][:])
            skim['ual'][numpy.abs(skim['ual']) > 10] = numpy.nan
            skim['ual'] = numpy.ma.masked_invalid(skim['ual'])
            fid.close()
            if numpy.shape(skim['ual'])[1] < nac:
                continue

            for i in range(3, nac-3):
                checknanref = + ref['uac'][:, i]
                checknanobs = + skim['uac'][:, i]
                merged_mask = (checknanref.mask | checknanobs.mask)
                indok = numpy.where(~merged_mask)[0]
                _ind = numpy.where(indok[1:] - indok[:-1] > 1)
                ensind = numpy.split(indok, numpy.cumsum(_ind) + 1)
                for chunk in range(len(ensind)):
                    if len(ensind[chunk]) > nal_min:
                        s1 = ref['uac'][ensind[chunk][:nal_min], i]
                        s2 = skim['uac'][ensind[chunk][:nal_min], i]
                        ffac, C = cpsd1d(hh1=s1, hh2=s2, dx=dal,
                                       tap=0.5, detrend=True)
                        countuac += 1
                        try:
                            MCuac += C
                        except:
                            MCuac = +C
                        s1 = ref['ual'][ensind[chunk][:nal_min], i]
                        s2 = skim['ual'][ensind[chunk][:nal_min], i]
                        try:
                            ffal, C = cpsd1d(hh1=s1, hh2=s2, dx=dal,
                                           tap=0.5, detrend=True)
                            countual += 1
                        except:
                            continue
                        try:
                            MCual += C
                        except:
                            MCual = +C

        MCuac = numpy.array(MCuac) / countuac
        MCual = numpy.array(MCual) / countual

        pyplot.semilogx(ffac, MCuac, label='{} {}'.format(ivar, iconfig)) # ,c='k')
        #pyplot.semilogx(ffal, MCual, label='along {} {}'.format(ivar, iconfig)) # ,c='k')
    pyplot.grid()
    pyplot.legend()
    pyplot.xlabel('/cy/km')
    pyplot.ylabel('coherency')
    pyplot.savefig('{}.png'.format(outfile))

def rms_l2c(datadir_input, config, output, threshold=0.1, fsize=1):
    datadir_output = output
    glob_files = os.path.join(datadir_input, '{}_l2c_c01*p*.nc'.format(config))
    list_files = glob.glob(glob_files)
    ref = {}
    skim = {}

    fid = netCDF4.Dataset(list_files[5], 'r')
    try:
        _tmp = numpy.array(fid.variables['u_ac_true'][:])
    except:
        print(list_files[0])
    #_tmp = numpy.array(fid.variables['u_ac_true'][:])
    _tmp[_tmp < -10] = numpy.nan
    fid.close()
    nac = _tmp.shape[1]
    nal = nac

    std_uac = numpy.zeros(nac)
    std_ual = numpy.zeros(nal)
    ntot_ac = numpy.zeros(nac)
    ntot_al = numpy.zeros(nal)
    std_uacm = numpy.zeros(nac)
    std_ualm = numpy.zeros(nal)
    ntot_acm = numpy.zeros(nac)
    ntot_alm = numpy.zeros(nal)
    print(nac)
    list_key = ['u_al_wd', 'u_ac_wd', 'u_al_wdrem', 'u_ac_wdrem',
                'u_ac_instr', 'u_al_instr', 'u_ac_dsigma', 'u_al_dsigma',
                'u_ac_true', 'u_al_true']
    list_fkey = ['uwnd', 'vwnd', 'rain', 'mssu', 'mssc']
    std_err = {}
    ntot_err = {}
    for ikey in list_key:
        std_err[ikey] = numpy.zeros(nac)
        ntot_err[ikey] = numpy.zeros(nac)

    for filev in list_files:
        fid = netCDF4.Dataset(filev, 'r')
        ipath = int(filev[-6:-3])
        #if ipath > 400:
        #    continue
       # if ipath%2==0:
        #   continue
        try:
            ref['uac'] = numpy.array(fid.variables['u_ac_true'][:])
            ref['uac'] = scipy.ndimage.uniform_filter(ref['uac'], size=fsize)
        except:
            print(filev)
            continue
        ref['uac'][numpy.abs(ref['uac']) > 10] = numpy.nan
        ref['ual'] = numpy.array(fid.variables['u_al_true'][:])
        ref['ual'] = scipy.ndimage.uniform_filter(ref['ual'], size=fsize)
        ref['ual'][numpy.abs(ref['ual']) > 10] = numpy.nan
        skim['uacm'] = numpy.array(fid.variables['u_ac_noerr'][:])
        skim['uacm'][numpy.abs(skim['uacm']) > 10] = numpy.nan
        skim['ualm'] = numpy.array(fid.variables['u_al_noerr'][:])
        skim['ualm'][numpy.abs(skim['ualm']) > 10] = numpy.nan
        skim['uac'] = numpy.array(fid.variables['u_ac_obs'][:])
        skim['uac'][numpy.abs(skim['uac']) > 10] = numpy.nan
        skim['ual'] = numpy.array(fid.variables['u_al_obs'][:])
        skim['ual'][numpy.abs(skim['ual']) > 10] = numpy.nan
        for ikey in list_fkey:
            skim[ikey] = numpy.array(fid.variables[ikey][:])
        wnd = numpy.sqrt(skim['uwnd']**2 + skim['vwnd']**2)
        for ikey in list_key:
            skim[ikey] = numpy.array(fid.variables[ikey][:])
            skim[ikey][numpy.abs(skim[ikey]) > 100] = numpy.nan
            skim['uac'][numpy.isnan(skim[ikey])] = numpy.nan
            skim['ual'][numpy.isnan(skim[ikey])] = numpy.nan
            skim[ikey][numpy.where((wnd > 5) & (wnd < 10))] = numpy.nan
            skim[ikey][numpy.where((skim['rain'] > 0.1))] = numpy.nan
        skim['uac'][numpy.where((wnd > 5) & (wnd < 10))] = numpy.nan
        skim['ual'][numpy.where((wnd > 5) & (wnd < 10))] = numpy.nan
        skim['uac'][numpy.where((skim['rain'] > 0.1))] = numpy.nan
        skim['ual'][numpy.where((skim['rain'] > 0.1))] = numpy.nan
        _indwdal = numpy.where(abs(skim['u_al_wdrem']) > 1)
        _indwdac = numpy.where(abs(skim['u_ac_wdrem']) > 1)
        skim['uac'][_indwdac] = numpy.nan
        skim['ual'][_indwdal] = numpy.nan
        #skim['u_ac_wdrem'][_indwdac] = numpy.nan
        #skim['u_al_wdrem'][_indwdal] = numpy.nan

        # skim['ual'][numpy.isnan(skim[ikey])] = numpy.nan
        nuac = ref['uac'].shape[1]
        fid.close()
        if nuac < 61:
            continue
        ind0 = 50
        ind1 = -50
       # skim['ual'] = skim['ualm'] + skim['u_al_instr'] +  skim['u_al_dsigma'] + skim['u_al_wdrem'] 
       # skim['uac'] = skim['uacm'] + skim['u_ac_instr'] +  skim['u_ac_dsigma'] + skim['u_ac_wdrem'] 
        for i in range(nuac):
            it_ac = len(numpy.where(numpy.isnan(skim['uac'][ind0:ind1, i]) == False)[0])
            if it_ac >= 100:
                _std = numpy.nanstd(skim['uac'][ind0:ind1, i] - ref['uac'][ind0:ind1, i])*it_ac
                if numpy.isfinite(_std):
                    std_uac[i] += _std
                    ntot_ac[i] += it_ac
                else:
                    std_uac[i] += 0
                    ntot_ac[i] += 0
            it_al = len(numpy.where(numpy.isnan(skim['ual'][ind0:ind1, i]) == False)[0])
            if it_al >= 100:
                _std = numpy.nanstd(skim['ual'][ind0:ind1, i]
                                    - ref['ual'][ind0:ind1, i])*it_al
                if numpy.isfinite(_std):
                    std_ual[i] += _std
                    ntot_al[i] += it_al
                else:
                    std_ual[i] += 0
                    ntot_al[i] += 0
            it_acm = len(numpy.where(numpy.isnan(skim['uacm'][ind0:ind1, i]) == False)[0])
            if it_acm >= 100:
                _std = numpy.nanstd(skim['uacm'][ind0:ind1, i]
                                    - ref['uac'][ind0:ind1, i])*it_acm
                if numpy.isfinite(_std):
                    std_uacm[i] += _std
                    ntot_acm[i] += it_acm
                else:
                    std_uacm[i] += 0
                    ntot_acm[i] += 0
            it_alm = len(numpy.where(numpy.isnan(skim['ualm'][ind0:ind1, i]) == False)[0])
            if it_alm >= 100:
                _std = numpy.nanstd(skim['ualm'][ind0:ind1, i]
                                    - ref['ual'][ind0:ind1, i])*it_alm
                if numpy.isfinite(_std):
                    std_ualm[i] += _std
                    ntot_alm[i] += it_alm
                else:
                    std_ualm[i] += 0
                    ntot_alm[i] += 0
            for ikey in list_key:
                it_alm = len(numpy.where(numpy.isnan(skim[ikey][ind0:ind1, i]) == False)[0])
                if it_alm >= 100:
                    _std = numpy.nanstd(skim[ikey][ind0:ind1, i]) * it_alm
                    if numpy.isfinite(_std):
                        std_err[ikey][i] += _std
                        ntot_err[ikey][i] += it_alm
                    else:
                        std_err[ikey][i] += 0
                        ntot_err[ikey][i] += 0

    std_uac = std_uac/ntot_ac
    std_ual = std_ual/ntot_al
    std_uacm = std_uacm/ntot_acm
    std_ualm = std_ualm/ntot_alm
    f, (ax1, ax2, ax3) = pyplot.subplots(1, 3, sharey= True, figsize=(20,5 ))
    xac = numpy.arange(-(nac - 1) * posting/2, (nac + 1)* posting/2, posting)
    indx = numpy.where(abs(xac)<= 40)
    #std_uac[indx] = numpy.nan
    #std_ual[indx] = numpy.nan
    print(xac[numpy.where(std_uac > threshold)])
    print(xac[numpy.where(std_ual > threshold)])
    _ind = numpy.where(numpy.abs(xac)>40)
    print(config, 'uac', numpy.nanmean(std_uac[_ind]))
    print(config, 'ual', numpy.nanmean(std_ual[_ind]))
    print(config, 'uacm', numpy.nanmean(std_uacm[_ind]))
    print(config, 'ualm', numpy.nanmean(std_ualm[_ind]))
    for ikey in list_key:
        std_err[ikey] = std_err[ikey] /  ntot_err[ikey]
        print(config, ikey, numpy.nanmean(std_err[ikey][_ind]))
    ax1.plot(xac, std_uac, 'r', label='across track')
    ax1.plot(xac, std_ual, 'b', label='along track')
    ax1.axhline(y=0.15, color="0.5")
    ax1.set_title('Observation {}'.format(config))
    ax1.set_ylim([0.00, 0.4])
    ax1.legend()
    ax2.plot(xac, std_uacm, 'r', label='across track')
    ax2.plot(xac, std_ualm, 'b', label='along track')
    ax2.axhline(y=0.15, color="0.5")
    ax2.set_title('Error-free {}'.format(config))
    ax2.legend()
    ax3.plot(xac, std_uacm, 'r', label='across track regridding')
    ax3.plot(xac, std_ualm, 'b', label='along track regridding')
    for ikey in list_key:
        std_err[ikey][indx] = numpy.nan
        ax3.plot(xac, std_err[ikey], label=ikey)
    ax3.axhline(y=0.15, color="0.5")
    ax3.set_title('Error decomposition {}'.format(config))
    ax3.legend()
    pyplot.savefig('{}.pdf'.format(datadir_output))

def bin_variables(listfile, listvar, bin_in, modelbox):
    lonp = numpy.arange(modelbox[0], modelbox[1] + modelbox[2], modelbox[2])
    latp = numpy.arange(modelbox[3], modelbox[4] + modelbox[5], modelbox[5])
    resol = numpy.sqrt((modelbox[2] * numpy.cos(numpy.deg2rad(latp))) **2
                        + modelbox[5]**2)
    dic_v = {}
    for ifile in listfile:
        print(ifile)
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
                if iiobs[0].shape == (0, ):
                    continue
                if ind_key not in dic_v.keys():
                    dic_v[ind_key] = {}
                for ivar in listvar:
                    ivaro = '{}_obs'.format(ivar)
                    ivart = '{}_noerr'.format(ivar)
                    ivarr = '{}_true'.format(ivar)
                    idvaro = '{}_diff_obs'.format(ivar)
                    idvart = '{}_diff_noerr'.format(ivar)
                    if ivar == 'norm':
                        varo = numpy.sqrt(data['u_ac_obs'][iiobs]**2
                                          + data['u_al_obs'][iiobs]**2)
                        vart = numpy.sqrt(data['u_ac_noerr'][iiobs]**2
                                          + data['u_al_noerr'][iiobs]**2)
                        varr = numpy.sqrt(data['u_ac_true'][iiobs]**2
                                          + data['u_al_true'][iiobs]**2)
                    else:
                        try:
                            varo = numpy.ma.array(data[ivaro][iiobs])
                        except: print(ifile, ivar, iiobs)
                        vart = data[ivart][iiobs]
                        varr = data[ivarr][iiobs]
                    mask = varo.mask
                    dvaro = numpy.array(abs(varo.data - varr.data))
                    dvart = numpy.array(abs(vart.data - varr.data))
                    if varo.any():
                        if ivaro not in dic_v[ind_key].keys():
                            dic_v[ind_key][ivaro] = []
                        if ivart not in dic_v[ind_key].keys():
                            dic_v[ind_key][ivart] = []
                        if ivarr not in dic_v[ind_key].keys():
                            dic_v[ind_key][ivarr] = []
                        if idvaro not in dic_v[ind_key].keys():
                            dic_v[ind_key][idvaro] = []
                        if idvart not in dic_v[ind_key].keys():
                            dic_v[ind_key][idvart] = []
                        dic_v[ind_key][ivaro].append(varo[~mask])
                        dic_v[ind_key][ivart].append(vart[~mask])
                        dic_v[ind_key][ivarr].append(varr[~mask])
                        dic_v[ind_key][idvart].append(dvart[~mask])
                        dic_v[ind_key][idvaro].append(dvaro[~mask])
                        dic_v[ind_key][ivart].append(varo[~mask])
                        dic_v[ind_key][ivart].append(varo[~mask])
                        dic_v[ind_key][ivart].append(varo[~mask])
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
    for ivar in listvar:
        _list = ['_obs', '_noerr', '_true', '_diff_obs', 'diff_noerr']
        for ivar2 in _list:
            ivar3 = '{}{}'.format(ivar, ivar2)
            rms[ivar3] = numpy.full((len(latp), len(lonp)), numpy.nan)
            std[ivar3] = numpy.full((len(latp), len(lonp)), numpy.nan)
            snr[ivar3] = numpy.full((len(latp), len(lonp)), numpy.nan)
    for j in range(len(rms['lon'])):
        for i in range(len(rms['lat'])):
            ind_key = 10000 * int(i) + int(j)
            if ind_key  not in dic_v.keys():
                continue
            # import pdb ; pdb.set_trace()
            for ivar in listvar:
                ivaro = '{}_obs'.format(ivar)
                ivart = '{}_noerr'.format(ivar)
                ivarr = '{}_true'.format(ivar)
                idvaro = '{}_diff_obs'.format(ivar)
                idvart = '{}_diff_noerr'.format(ivar)
                if ivaro not in dic_v[ind_key].keys():
                    continue
                varo = numpy.concatenate(dic_v[ind_key][ivar]).ravel()
                vart = numpy.concatenate(dic_v[ind_key][ivar]).ravel()
                varr = numpy.concatenate(dic_v[ind_key][ivar]).ravel()
                rms_true = numpy.sqrt(numpy.nanmean(varr**2))
                rms[ivaro][i, j] = numpy.sqrt(numpy.nanmean(varo**2))
                std[ivaro][i, j] = numpy.nanstd(varo)
                snr[ivaro][i, j] = rms_true / std[idvaro][i, j]
                rms[ivart][i, j] = numpy.sqrt(numpy.nanmean(vart**2))
                std[ivart][i, j] = numpy.nanstd(vart)
                snr[ivart][i, j] = rms_true / std[idvart][i, j]
                rms[ivarr][i, j] = numpy.sqrt(numpy.nanmean(varr**2))
                std[ivarr][i, j] = numpy.nanstd(varr)
    with open('rms_{}'.format(bin_out), 'wb') as f:
        pickle.dump(rms, f)

    with open('std_{}'.format(bin_out), 'wb') as f:
        pickle.dump(std, f)

    with open('snr_{}'.format(bin_out), 'wb') as f:
        pickle.dump(snr, f)


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
    list_var = rms.keys()
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
    return None


if '__main__' == __name__:
    if len(sys.argv) < 1:
        print('Provide json file for diagnostics')
        sys.exit(1)
    file_param = sys.argv[1]
    with open(file_param, 'r') as f:
        params = json.load(f)
 #   rms_l2c(p.outdatadir, p.config)
    pl2c = params['l2c']
    length_al = pl2c['alongtrack_length']
    list_config = pl2c['list_config']

    list_dir = pl2c['indatadir']
    posting = pl2c['posting']
    outdir = pl2c['outdatadir']
    list_filter_size = pl2c['filter_size']
    for i, iconfig in enumerate(list_config):
        filter_size = list_filter_size[i]
        indatadir = list_dir[i]
 #       print(indatadir, iconfig)
        outfile = os.path.join(outdir, 'std_{}'.format(iconfig))
        rms_l2c(indatadir, iconfig, outfile, threshold=0.15, fsize=filter_size)
        outfile = os.path.join(outdir, 'coherency_{}_obs_model'.format(iconfig))
        #coherency_l2c((indatadir, indatadir), (iconfig, iconfig),
        #             ('obs','noerr'), length_al,
        #             posting, outfile=outfile, fsize=filter_size)

    list_var = ('obs', 'obs', 'obs', 'obs') #, 'obs')
    nal_min = length_al
#    print(list_dir, list_config)
    outfile = os.path.join(outdir, 'coherency_obs_{}'.format(list_config[0][:-8]))
    #coherency_l2c(list_dir, list_config, list_var, nal_min,
    #              posting, outfile=outfile, fsize=filter_size)
    list_var = ('noerr', 'noerr', 'noerr', 'noerr') #, 'obs')
    print(list_dir, list_config)
    outfile = os.path.join(outdir, 'coherency_noerr_{}'.format(list_config[0][:-8]))
    #coherency_l2c(list_dir, list_config, list_var, nal_min,
    #              posting, outfile=outfile, fsize=filter_size)
    listvar = ['u_ac', 'u_al', 'ux', 'uy', 'norm']
    modelbox2 = params['l2b']['modelbox_bin']
    for i, iconfig in enumerate(list_config):
        filter_size = list_filter_size[i]
        indatadir = list_dir[i]
        print(indatadir, iconfig)
        files = os.path.join(indatadir, '{}_l2c_c01'.format(iconfig))
        listfiles = glob.glob('{}*.nc'.format(files))
        bin_file = os.path.join(outdir, '{}_l2c.pyo'.format(iconfig))
    #    bin_variables(listfiles[:10], listvar, bin_file, modelbox2)
    #    bin_file2 = '{}_l2c.pyo'.format(iconfig)
    #    compute_rms(bin_file, bin_file2, listvar, modelbox2)

    #    plot_rms(bin_file2, listvar, iconfig)

