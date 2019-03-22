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
                  posting, outfile='coherency'):
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
            except:
                print(ifile)
                continue
            ref['uac'][ref['uac']<-10] = numpy.nan
            ref['uac'] = numpy.ma.masked_invalid(ref['uac'])
            ref['ual'] = numpy.array(fid.variables['u_al_true'][:])
            ref['ual'][ref['ual']<-10] = numpy.nan
            ref['ual'] = numpy.ma.masked_invalid(ref['ual'])
            skim['uac'] = numpy.array(fid.variables['u_ac_{}'.format(ivar)][:])
            skim['uac'][skim['uac']<-10] = numpy.nan
            skim['uac'] = numpy.ma.masked_invalid(skim['uac'])
            skim['ual'] = numpy.array(fid.variables['u_al_{}'.format(ivar)][:])
            skim['ual'][skim['ual']<-10] = numpy.nan
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

def rms_l2c(datadir_input, config, output, threshold=0.1):
    datadir_output = output
    glob_files = os.path.join(datadir_input, '{}_l2c_c01*p*.nc'.format(config))
    list_files = glob.glob(glob_files)
    ref = {}
    skim = {}

    fid = netCDF4.Dataset(list_files[0], 'r')
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

    for filev in list_files:
        fid = netCDF4.Dataset(filev, 'r')
        ipath = int(filev[-6:-3])
        if ipath > 400:
            continue
       # if ipath%2==0:
        #   continue
        try:
            ref['uac'] = numpy.array(fid.variables['u_ac_true'][:])
        except:
            print(filev)
            continue
        ref['uac'][ref['uac'] < -10] = numpy.nan
        ref['ual'] = numpy.array(fid.variables['u_al_true'][:])
        ref['ual'][ref['ual'] < -10] = numpy.nan
        skim['uacm'] = numpy.array(fid.variables['u_ac_noerr'][:])
        skim['uacm'][skim['uacm'] < -10] = numpy.nan
        skim['ualm'] = numpy.array(fid.variables['u_al_noerr'][:])
        skim['ualm'][skim['ualm'] < -10] = numpy.nan
        skim['uac'] = numpy.array(fid.variables['u_ac_obs'][:])
        skim['uac'][skim['uac'] < -10] = numpy.nan
        skim['ual'] = numpy.array(fid.variables['u_al_obs'][:])
        skim['ual'][skim['ual'] < -10] = numpy.nan
        nuac = ref['uac'].shape[1]
        fid.close()
        if nuac < 61:
            continue
        for i in range(nuac):
            it_ac = len(numpy.where(numpy.isnan(skim['uac'][:, i]) == False)[0])
            if it_ac >= 62:
                std_uac[i] += numpy.nanstd(skim['uac'][:, i]
                                           - ref['uac'][:, i])*it_ac
                ntot_ac[i] += it_ac
            it_al = len(numpy.where(numpy.isnan(skim['ual'][:, i]) == False)[0])
            if it_al >= 62:
                std_ual[i] += numpy.nanstd(skim['ual'][:, i]
                                           - ref['ual'][:, i])*it_al
                ntot_al[i] += it_al
            it_acm = len(numpy.where(numpy.isnan(skim['uacm'][:, i]) == False)[0])
            if it_acm >= 62:
                std_uacm[i] += numpy.nanstd(skim['uacm'][:, i]
                                            - ref['uac'][:, i])*it_acm
                ntot_acm[i] += it_acm
            it_alm = len(numpy.where(numpy.isnan(skim['ualm'][:, i]) == False)[0])
            if it_alm >= 62:
                std_ualm[i] += numpy.nanstd(skim['ualm'][:, i]
                                            - ref['ual'][:, i])*it_alm
                ntot_alm[i] += it_alm

    std_uac = std_uac/ntot_ac
    std_ual = std_ual/ntot_al
    std_uacm = std_uacm/ntot_acm
    std_ualm = std_ualm/ntot_alm
    f, (ax1, ax2) = pyplot.subplots(1, 2, sharey= True, figsize=(12,5 ))
    xac = numpy.arange(-(nac - 1) * posting/2, (nac + 1)* posting/2, posting)
    print(xac[numpy.where(std_uac > threshold)])
    print(xac[numpy.where(std_ual > threshold)])
    _ind = numpy.where(numpy.abs(xac)>40)
    print(config, 'uac', numpy.nanmean(std_uac[_ind]))
    print(config, 'ual', numpy.nanmean(std_ual[_ind]))
    print(config, 'uacm', numpy.nanmean(std_uacm[_ind]))
    print(config, 'ualm', numpy.nanmean(std_ualm[_ind]))
    ax1.plot(xac, std_uac, 'r', label='across track')
    ax1.plot(xac, std_ual, 'b', label='along track')
    ax1.axhline(y=0.1, color="0.5")
    ax1.set_title('Observation {}'.format(config))
    ax1.set_ylim([0.00, 0.18])
    ax1.legend()
    ax2.plot(xac, std_uacm, 'r', label='across track')
    ax2.plot(xac, std_ualm, 'b', label='along track')
    ax2.set_title('Error-free {}'.format(config))
    ax2.legend()
    pyplot.savefig('{}.png'.format(datadir_output))


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
    for i, iconfig in enumerate(list_config):
        indatadir = list_dir[i]
        print(indatadir, iconfig)
        outfile = os.path.join(outdir, 'std_{}'.format(iconfig))
        rms_l2c(indatadir, iconfig, outfile)
        outfile = os.path.join(outdir, 'coherency_{}_obs_model'.format(iconfig))
        coherency_l2c((indatadir, indatadir), (iconfig, iconfig),
                     ('obs','noerr'), length_al,
                     posting, outfile=outfile)

    list_var = ('obs', 'obs', 'obs', 'obs') #, 'obs')
    nal_min = length_al
    print(list_dir, list_config)
    outfile = os.path.join(outdir, 'coherency_obs_{}'.format(list_config[0][:-8]))
    coherency_l2c(list_dir, list_config, list_var, nal_min,
                  posting, outfile=outfile)
    list_var = ('noerr', 'noerr', 'noerr', 'noerr') #, 'obs')
    print(list_dir, list_config)
    outfile = os.path.join(outdir, 'coherency_noerr_{}'.format(list_config[0][:-8]))
    coherency_l2c(list_dir, list_config, list_var, nal_min,
                  posting, outfile=outfile)
