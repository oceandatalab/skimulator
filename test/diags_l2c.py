import numpy
from matplotlib import pyplot
from netCDF4 import Dataset
import params as p
import glob
import os

datadir_input = p.outdatadir
datadir_output = './'

config = p.config
nac = 63
nal = 63
ref = {}
skim = {}

std_uac = numpy.zeros(nac)
std_ual = numpy.zeros(nal)
ntot_ac = numpy.zeros(nac)
ntot_al = numpy.zeros(nal)

list_files = glob.glob(os.path.join(datadir_input, '{}_L2C_*'.format(config)))
for filev in list_files:
    fid = Dataset(filev, 'r')
    ref['uac'] = numpy.array(fid.variables['u_ac_model'][:])
    ref['uac'][ref['uac'] < -10] = numpy.nan
    ref['ual'] = numpy.array(fid.variables['u_al_model'][:])
    ref['ual'][ref['ual'] < -10] = numpy.nan
    skim['uac'] = numpy.array(fid.variables['u_ac_obs'][:])
    skim['uac'][skim['uac'] < -10] = numpy.nan
    skim['ual'] = numpy.array(fid.variables['u_al_obs'][:])
    skim['ual'][skim['ual'] < -10] = numpy.nan

    for i in range(ref['uac'].shape[1]):
        it_ac = len(numpy.where(numpy.isnan(skim['uac'][:, i]) == False)[0])
        if it_ac >= 1:
            std_uac[i] += numpy.nanstd(skim['uac'][:, i] - ref['uac'][:, i])*it_ac
            ntot_ac[i] += it_ac
        it_al = len(numpy.where(numpy.isnan(skim['ual'][:, i]) == False)[0])
        if it_al >= 1:
          std_ual[i] += numpy.nanstd(skim['ual'][:, i] - ref['ual'][:, i])*it_al
          ntot_al[i] += it_al

std_uac = std_uac/ntot_ac
std_ual = std_ual/ntot_al

pyplot.figure()
xac = numpy.arange(-(nac - 1) * p.posting/2, (nac + 1)* p.posting/2, p.posting)
pyplot.plot(xac, std_uac, 'r', label='across track')
pyplot.plot(xac, std_ual, 'b', label='along cross')
pyplot.title(config)
pyplot.legend()
pyplot.savefig('std_{}.png'.format(config))
