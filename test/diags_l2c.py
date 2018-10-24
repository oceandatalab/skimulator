import numpy
from matplotlib import pyplot
from netCDF4 import Dataset
import params as p
import glob
import os

datadir_input = p.outdatadir
datadir_output = './'

config = p.config
list_files = glob.glob(os.path.join(datadir_input, '{}_L2C_*.nc'.format(config)))
nac = 63
nal = 63
ref = {}
skim = {}

fid = Dataset(list_files[0], 'r')
_tmp = numpy.array(fid.variables['u_ac_true'][:])
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
    fid = Dataset(filev, 'r')
    ipath = int(filev[-6:-3])
    if ipath%2==0:
       continue
    print(filev)
    ref['uac'] = numpy.array(fid.variables['u_ac_true'][:])
    ref['uac'][ref['uac'] < -10] = numpy.nan
    ref['ual'] = numpy.array(fid.variables['u_al_true'][:])
    ref['ual'][ref['ual'] < -10] = numpy.nan
    skim['uacm'] = numpy.array(fid.variables['u_ac_model'][:])
    skim['uacm'][skim['uacm'] < -10] = numpy.nan
    skim['ualm'] = numpy.array(fid.variables['u_al_model'][:])
    skim['ualm'][skim['ualm'] < -10] = numpy.nan
    skim['uac'] = numpy.array(fid.variables['u_ac_obs'][:])
    skim['uac'][skim['uac'] < -10] = numpy.nan
    skim['ual'] = numpy.array(fid.variables['u_al_obs'][:])
    skim['ual'][skim['ual'] < -10] = numpy.nan
    nuac = ref['uac'].shape[1]
    print(nuac, nac)
    fid.close()
    for i in range(nuac):
        it_ac = len(numpy.where(numpy.isnan(skim['uac'][:, i]) == False)[0])
        if it_ac >= 62:
            std_uac[i] += numpy.nanstd(skim['uac'][:, i] - ref['uac'][:, i])*it_ac
            ntot_ac[i] += it_ac
        it_al = len(numpy.where(numpy.isnan(skim['ual'][:, i]) == False)[0])
        if it_al >= 62:
            std_ual[i] += numpy.nanstd(skim['ual'][:, i] - ref['ual'][:, i])*it_al
            ntot_al[i] += it_al
        it_acm = len(numpy.where(numpy.isnan(skim['uacm'][:, i]) == False)[0])
        if it_acm >= 62:
            std_uacm[i] += numpy.nanstd(skim['uacm'][:, i] - ref['uac'][:, i])*it_acm
            ntot_acm[i] += it_acm
        it_alm = len(numpy.where(numpy.isnan(skim['ualm'][:, i]) == False)[0])
        if it_alm >= 62:
            std_ualm[i] += numpy.nanstd(skim['ualm'][:, i] - ref['ual'][:, i])*it_alm
            ntot_alm[i] += it_alm

std_uac = std_uac/ntot_ac
std_ual = std_ual/ntot_al
std_uacm = std_uacm/ntot_acm
std_ualm = std_ualm/ntot_alm

f, (ax1, ax2) = pyplot.subplots(1, 2, sharey= True, figsize=(12,5 ))
xac = numpy.arange(-(nac - 1) * p.posting/2, (nac + 1)* p.posting/2, p.posting)
ax1.plot(xac, std_uac, 'r', label='across track')
ax1.plot(xac, std_ual, 'b', label='along cross')
ax1.set_title('Observation {}'.format(config))
ax1.set_ylim([0.05, 0.24])
ax1.legend()
ax2.plot(xac, std_uacm, 'r', label='across track')
ax2.plot(xac, std_ualm, 'b', label='along cross')
ax2.set_title('Error-free {}'.format(config))
ax2.legend()
pyplot.savefig('std_{}.png'.format(config))
