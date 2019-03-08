'''Module to create one beam data:


\n
\n
#-----------------------------------------------------------------------
#                       Additional Documentation
# Authors: Lucile Gaultier
#
# Modification History:
# - Mar 2017:  Original by Lucile Gaultier, ODL
#
# Notes:
# - Written for Python 2.7,  Python 3.5, tested with Python 2.7, Python 3.5
#
#-----------------------------------------------------------------------
Copyright (C) 2017-2018 OceanDataLab
This file is part of skimulator.

skimulator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

skimulator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with skimulator.  If not, see <http://www.gnu.org/licenses/>.
'''
import os
from scipy import interpolate
import numpy
import math
import glob
import sys
import time
import pickle
import datetime
import logging
import skimulator.build_swath as build_swath
import skimulator.rw_data as rw_data
import skimulator.build_error as build_error
import skimulator.mod_tools as mod_tools
import skimulator.const as const
import multiprocessing
# Define logger level for debug purposes
logger = logging.getLogger(__name__)
#logger = multiprocessing.log_to_stderr()
#logger.setLevel(logging.DEBUG)


def load_error(p):
    '''Initialize random coefficients that are used to compute
    random errors following the specified spectrum. \n
    If a random coefficient file is specified, random coefficients
    are loaded from this file.
    '''
    err = build_error.error(p)
    err.init_error(p)
    errnad = build_error.errornadir(p)
    errnad.init_error(p)
    return err, errnad


def load_sgrid(sgridfile, p):
    '''Load SKIM swath and Nadir data for file sgridfile '''

    # Load SKIM swath file
    sgrid = rw_data.Sat_SKIM(ifile=sgridfile)
    cycle = 0
    x_al = []
    x_ac = []
    al_cycle = 0
    timeshift = 0
    sgrid.load_swath(p, cycle=cycle, x_al=x_al, x_al_nadir=[], x_ac=x_ac,
                     al_cycle=al_cycle, timeshift=timeshift)
    # sgrid.loncirc = []
    # for i in range(len(sgrid.lon)):
    #     sgrid.loncirc.append(numpy.rad2deg(numpy.unwrap(sgrid.lon[i])))
    # Extract the pass number from the file name
    ipass = int(sgridfile[-6: -3])
    sgrid.ipass = ipass
    return sgrid


def load_coordinate_model(p):
    model = p.model
    if p.file_input is not None:
        list_file = [line.strip() for line in open(p.file_input)]
    else:
        list_file = None
    # - Read model input coordinates '''
    # If a list of model files are specified, read model file coordinates
    if p.file_input is not None:
        model_data_ctor = getattr(rw_data, model)
        if p.file_grid_model is not None:
            _filename = list(p.file_grid_model)
        else:
            logger.info("WARNING: First file of list of files is used for"
                        "coordinates only")
            _filename = list_file[0].split(',')
        filename = []
        for ifile in _filename:
            filename.append(os.path.join(p.indatadir, ifile))

        #filename = os.path.join(p.indatadir, list_file[0])
        model_data = model_data_ctor(p, ifile=filename,
                                     lon=list(p.lon), lat=list(p.lat))
    return model_data, list_file


def interpolate_regular_1D(lon_in, lat_in, var, lon_out, lat_out, Teval=None):
    ''' Interpolation of data when grid is regular and coordinate in 1D. '''
    # To correct for IDL issues
    ind_sort = numpy.argsort(lon_in)
    lon_in = lon_in[ind_sort]
    var = var[:, ind_sort]
    if numpy.max(lon_in) > 359 and numpy.min(lon_in) < 1:
        ind_in1 = numpy.where(lon_in <= 180)
        ind_in2 = numpy.where(lon_in > 180)
        # lon_in[lon_in > 180] = lon_in[lon_in > 180] - 360
        # lon_in = np.mod(lon_in - (lref - 180), 360) + (lref - 180)
        # lon_in = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(lon_in)))
        ind_out1 = numpy.where(lon_out <= 180)
        ind_out2 = numpy.where(lon_out > 180)
        # lon_out[lon_out > 180] = lon_out[lon_out > 180] - 360
        # lon_out = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(lon_out)))
        interp = interpolate.RectBivariateSpline
        mask_teval = (numpy.isnan(var) | numpy.ma.getmaskarray(var))
        if Teval is None:
            Teval = numpy.zeros(numpy.shape(lon_out))
            if ind_out1[0].any():
                _tmp = interp(lat_in, lon_in[ind_in1],
                              mask_teval[:, ind_in1[0]], kx=1, ky=1,
                              s=0)
                Teval[ind_out1] = _tmp.ev(lat_out[ind_out1], lon_out[ind_out1])
            if ind_out2[0].any():
                _tmp = interp(lat_in, lon_in[ind_in2],
                              mask_teval[:, ind_in2[0]], kx=1, ky=1,
                              s=0)
                Teval[ind_out2] = _tmp.ev(lat_out[ind_out2], lon_out[ind_out2])
        # Trick to avoid nan in interpolation
        var_mask = + var
        var_mask[numpy.isnan(var_mask) | numpy.ma.getmaskarray(var_mask)] = 0.
        # Interpolate variable
        var_out = numpy.full(numpy.shape(lon_out), numpy.nan)
        if ind_out1[0].any():
            _tmp = interp(lat_in, lon_in[ind_in1], var_mask[:, ind_in1[0]],
                          kx=1, ky=1, s=0)
            var_out[ind_out1] = _tmp.ev(lat_out[ind_out1], lon_out[ind_out1])
        if ind_out2[0].any():
            _tmp = interp(lat_in, lon_in[ind_in2], var_mask[:, ind_in2[0]],
                          kx=1, ky=1, s=0)
            var_out[ind_out2] = _tmp.ev(lat_out[ind_out2], lon_out[ind_out2])

    else:
        mask_teval = (numpy.isnan(var) | numpy.ma.getmaskarray(var))
        # Interpolate mask if it has not been done (Teval is None)
        interp = interpolate.RectBivariateSpline
        if Teval is None:
            _Teval = interp(lat_in, lon_in, mask_teval, kx=1, ky=1, s=0)
            Teval = _Teval.ev(lat_out, lon_out)
        # Trick to avoid nan in interpolation
        var_mask = + var
        var_mask[numpy.isnan(var_mask) | numpy.ma.getmaskarray(var_mask)] = 0.
        # Interpolate variable
        _var_out = interp(lat_in, lon_in, var_mask, kx=1, ky=1, s=0)
        var_out = _var_out.ev(lat_out, lon_out)
    # Mask variable with Teval
    var_out[Teval > 0] = numpy.nan
    #var_out[Teval > 0] = numpy.nan
    #var_out[abs(var_out) > 1000] = numpy.nan
    return var_out, Teval


def interpolate_irregular_pyresample(swath_in, var, grid_out, radius,
                                     interp_type='nearest'):
    ''' Interpolation of data when grid is irregular and pyresample is
    installed.'''
    import pyresample as pr
    interp = pr.kd_tree.resample_nearest
    if interp_type == 'nearest':
        interp = pr.kd_tree.resample_nearest
        radius_n = radius * 10**3
        var_out = interp(swath_in, var, grid_out, radius_of_influence=radius_n,
                         epsilon=100)
    else:
        interp = pr.kd_tree.resample_gauss
        radius_g = radius * 3 * 10**3
        sigma_g = radius * 10**3
        var_out = interp(swath_in, var, grid_out, radius_of_influence=radius_g,
                         sigmas=sigma_g, fill_value=None)
    return var_out


def create_SKIMlikedata(cycle, list_file, modelbox,
                        sgrid, model_data, modeltime,
                        p):
    '''Create SKIM and nadir errors err and errnad, interpolate model velocity\
    model_data on swath and nadir track,
    compute SKIM-like and nadir-like data for cycle, SKIM grid sgrid and
    ngrid. '''
    #   Initialiaze errors and velocity
    output_var_i = {}
    shape_i = numpy.shape(sgrid.lon)[0]
    for key in p.list_output:
        output_var_i[key] = numpy.full(shape_i, numpy.nan)
    #for key in p.list_err:
    #    err_var_i[key] = numpy.full(shape_i, numpy.nan)

    date1 = cycle * sgrid.cycle
    # Definition of the time in the model
    time = sgrid.time / 86400. + date1  # in days
    lon = sgrid.lon
    lat = sgrid.lat
    timeshift = sgrid.timeshift / 86400.  # in days
    # Look for satellite data that are beween step-p.timestep/2 and
    # step+p.timestep/2
    if p.file_input is not None:
        lon2D = {}
        lat2D = {}
        # meshgrid in 2D for interpolation purposes
        for key in model_data.vlon.keys():
            if p.grid == 'irregular':
                lon2D[key], lat2D[key] = numpy.meshgrid(model_data.vlon[key],
                                                        model_data.vlat[key])
        index_filemodel = numpy.where(((time[-1] - timeshift) >=
                                      (modeltime-p.timestep/2.))
                                      & ((time[0] - timeshift) <
                                      (modeltime+p.timestep/2.)))
        # local variable to find time record in file for WW3
        nfile = 0
        time_offset = 0
        # At each step, look for the corresponding time in the satellite data
        for ifile in index_filemodel[0]:
            # If there are satellite data, Get true velcoity from model
            if numpy.shape(index_filemodel)[1] > 0:
                # Select part of the track that corresponds to the time of the
                # model (+-timestep/2)
                ind_time = numpy.where(((time-timeshift) >=
                                       (modeltime[ifile]-p.timestep/2.))
                                       & ((time-timeshift) <
                                       (modeltime[ifile]+p.timestep/2.)))
            else:
                logger.error('No model file is found around time')
                sys.exit(1)
            # Load data from this model file
            # if output from ww3, time dimension is >1 (hourly outputs,
            # one file per month): conversion of ifile into file number
            # and record number
            ## TODO : Clean
            filetime = (ifile - time_offset)%p.dim_time
            # read next file when the time dimension is reached
            if filetime >= (time_offset + p.dim_time):
                time_offset += p.dim_time
                nfile += 1
                filetime = (ifile - time_offset)%p.dim_time
            nfile = int(ifile /p.dim_time)
            filetime = ifile - nfile * p.dim_time

            _tmpfilename = list_file[nfile]
            filename = os.path.join(p.indatadir, _tmpfilename)

            model_step_ctor = getattr(rw_data, model_data.model)
            model_step = model_step_ctor(p, ifile=(filename, ),
                                         list_input_var=p.list_input_var,
                                         time=filetime)
            input_var_i = {}
            if p.grid == 'regular':
                model_step.read_var(p, ind_lon=model_data.ind_lon)
                for key in model_step.input_var.keys():
                    grid_key = model_step.numgrid[key]
                    _indlat = model_data.model_index_lat[grid_key]
                    _tmp = model_step.input_var[key][_indlat, :]
                    _indlon = model_data.model_index_lon[grid_key]
                    input_var_i[key] = +_tmp[:, _indlon]

            else:
                model_step.read_var(p, index=None)
                for key in model_step.input_var.keys():
                    _ind = model_data.model_index[model_step.numgrid[key]]
                    input_var_i[key] = + model_step.input_var[key][_ind]
            # - Interpolate Model data on a SKIM grid and/or along the
            #   nadir track
            # if grid is regular, use interpolate.RectBivariateSpline to
            # interpolate
            if p.grid == 'regular' and \
                    len(numpy.shape(model_data.vlon[0])) == 1:
                # Flatten satellite grid and select part of the track
                # corresponding to the model time
                for key in model_step.input_var.keys():
                    mlon = model_data.vlon[model_step.numgrid[key]]
                    mlat = model_data.vlat[model_step.numgrid[key]]
                    _tmp, Teval_u = interpolate_regular_1D(mlon, mlat,
                                                           input_var_i[key],
                                                           lon[ind_time[0]],
                                                           lat[ind_time[0]])
                    output_var_i[key][ind_time[0]] = + _tmp
            else:
                # Grid is irregular, interpolation can be done using
                # pyresample module if it is installed or griddata
                # function from scipy.
                # Note that griddata is slower than pyresample functions.
                try:
                    import pyresample as pr
                    lon_wrap = pr.utils.wrap_longitudes
                    geom = pr.geometry.SwathDefinition
                    lon = lon_wrap(lon)

                    grid_def = geom(lons=lon, lats=lat)
                    for key in model_step.input_var.keys():
                        grid_key = model_step.numgrid[key]
                        _ind =  model_data.model_index[grid_key]
                        _mlon = lon_wrap(model_data.vlon[grid_key])
                        _mlat = model_data.vlat[grid_key]
                        if len(_mlon[0]) <= 1:
                            _mlon = lon_wrap(lon2D[grid_key])
                            _mlat = lat2D[grid_key]
                        swath_def = geom(lons=_mlon[_ind], lats=_mlat[_ind])
                        _tmp = interpolate_irregular_pyresample(
                                          swath_def, input_var_i[key],
                                          grid_def,
                                          p.posting,
                                          interp_type=p.interpolation)
                        output_var_i[key][ind_time[0]] = + _tmp
                except ImportError:
                    for key in model_step.input_var.keys():
                        grid_key = model_step.numgrid[key]
                        _ind =  model_data.model_index[grid_key]
                        _mlon = model_data.vlon[grid_key]
                        _mlat = model_data.vlat[grid_key]
                        if len(_mlon) <= 1:
                            _mlon = lon2D[grid_key]
                            _mlat = lat2D[grid_key]
                        lonravel = + _mlon[_ind].ravel()
                        latravel = + _mlat[_ind].ravel()
                        _tmp = + input_var_i[key].ravel()
                        interp = interpolate.griddata((lonravel, latravel),
                                                      _tmp, (lon[ind_time[0]],
                                                      lat[ind_time[0]]),
                                                      method=p.interpolation)
                        output_var_i[key][ind_time[0]] = + interp
            # Force value outside modelbox at nan
            if modelbox[0] > modelbox[1]:
                for key in model_step.input_var.keys():
                    _ind = numpy.where(((lon > modelbox[0])
                                       & (lon < modelbox[1]))
                                       | (lat < modelbox[2])
                                       | (lat > modelbox[3]))
                    output_var_i[key][_ind] = numpy.nan
            else:
                for key in model_step.input_var.keys():
                    _ind = numpy.where((lon < modelbox[0])
                                       | (lon > modelbox[1])
                                       | (lat < modelbox[2])
                                       | (lat > modelbox[3]))
                    output_var_i[key][_ind] = numpy.nan
            output_var_i['vindice'][ind_time[0]] = ifile

    else:
        pass
    # TODO to proof: creation of empty mss if no mss and p.instr is True
    return output_var_i, time


def load_rain(rain_file):
    with open(rain_file, 'rb') as frain:
        dic = pickle.load(frain)
    size_dic = len(dic['xac'].keys())
    return dic, size_dic


def compute_rain(p, time, sgrid, dic, size_dic):
    hour = int((time - numpy.floor(time))*24)
    size_dic = len(dic['xal'][hour])
    rr_ind = int(numpy.random.random_sample() * size_dic)
    xal = dic['xal'][hour][rr_ind]
    var = dic['rr'][hour][rr_ind]
    xac = dic['xac'][hour][rr_ind]
    x_al_g_tot = + sgrid.x_al
    for i in range(numpy.shape(sgrid.x_al)[1]):
        x_al_g_tot[:, i] = sgrid.x_al[:, i] + sgrid.x_al_nadir
    xal_g = numpy.mod(x_al_g_tot - numpy.min(x_al_g_tot), numpy.max(xal))
    interp = interpolate.RectBivariateSpline
    _Teval = interp(xal, xac, numpy.isnan(var), kx=1, ky=1, s=0)
    Teval = _Teval.ev(xal_g, sgrid.x_ac)
        # Trick to avoid nan in interpolation
    var_mask = + var
    var_mask[numpy.isnan(var_mask)] = 0.
    # Interpolate variable
    _var_out = interp(xal, xac, var_mask, kx=1, ky=1, s=0)
    var_out = _var_out.ev(xal_g, sgrid.x_ac)
    # Mask variable with Teval
    var_out[Teval > 0] = numpy.nan
    xal_n = numpy.mod(sgrid.x_al_nadir - numpy.min(sgrid.x_al_nadir),
                      numpy.max(xal))
    var_nad = numpy.zeros(numpy.shape(xal_n))
    return var_out, var_nad


def compute_beam_noise_skim(p, output_var_i, radial_angle, beam_angle):
    output_var_i['ur_true'] = mod_tools.proj_radial(output_var_i['ucur'],
                                                    output_var_i['vcur'],
                                                    radial_angle)
    output_var_i['ur_obs'] = + output_var_i['ur_true']
    output_var_i['radial_angle'] = radial_angle
    if p.instr is True:
        # Compute sigma0:
        sigma0 = compute_sigma(output_var_i, beam_angle, radial_angle, p)
        output_var_i['sigma0'] = sigma0
        coeff_random = p.snr_coeff * output_var_i['sigma0']
        cshape = numpy.shape(coeff_random)
        center = numpy.zeros(cshape)
        output_var_i['instr'] = numpy.random.rand(cshape[0]) * coeff_random
        output_var_i['ur_obs'] += output_var_i['instr']

    # del output_var_i['mssx']
    # del output_var_i['mssy']
    # del output_var_i['mssxy']
    # Radial projection
    if p.uwb is True:
        output_var_i['ur_uss'] = mod_tools.proj_radial(output_var_i['uuss'],
                                                        output_var_i['vuss'],
                                                        radial_angle)
        nwr = numpy.sqrt((output_var_i['uwnd'] - output_var_i['ucur'])**2
                         + (output_var_i['vwnd'] - output_var_i['vcur'])**2)
        nwr[nwr==0] = numpy.nan
        _angle = numpy.deg2rad((beam_angle - 25) / 10)
        GR = 25 * (0.82 * numpy.log(0.2 + 7/nwr)) * (1 - numpy.tanh(_angle))
        GP = 0
        output_var_i['uwb'] = GR * output_var_i['ur_uss']
        #output_var_i['ur_obs'] +=  output_var_i['uwb']
    return None


def compute_nadir_noise_skim(p, output_var_i, sgrid, cycle):
    output_var_i['ssh_obs'] = + output_var_i['ssh']
    output_var_i['ssh_true'] = + output_var_i['ssh']
    if p.nadir is True:
        errnad = build_error.errornadir(p)
        errnad.init_error(p)
        errnad.make_error(sgrid, cycle, p)
        output_var_i['instr'] = errnad.nadir
    output_var_i['ssh_obs'] += output_var_i['instr']
    shape_0 = numpy.shape(sgrid.lon)
    for key in p.list_output:
        if key not in output_var_i.keys():
            output_var_i[key] = numpy.full(shape_0, numpy.nan)
    #if p.wet_tropo is True:
    #    output_var_i['wet_tropo'] = 
    return None


def compute_sigma_water(input_var, beam_angle, radial_angle):
    required = ('mssx', 'mssy', 'mssd', 'uwnd', 'vwnd', 'uwnd', 'ucur',
                'vcur')
    missing = [_ for _ in required if _ not in input_var.keys()]
    if 0 < len(missing):
        logger.info('Missing file to compute sigma, instrumental error not'
                    ' computed')
        logger.info('Missing parameters: {}'.format(', '.join(missing)))
        return None
    else:
        mssd = numpy.deg2rad(input_var['mssd'])
        mssu = + input_var['mssx']
        mssc = + input_var['mssy']
        uwnd = + input_var['uwnd']
        vwnd = + input_var['vwnd']
        ucur = + input_var['ucur']
        vcur = + input_var['vcur']
        mssxl = mssu * numpy.cos(mssd)**2 + mssc * numpy.sin(mssd)**2
        mssyl = mssu * numpy.sin(mssd)**2 + mssc * numpy.cos(mssd)**2
        mssxyl = (mssu - mssc) * numpy.sin(2 * mssd) / 2
        nwr = numpy.sqrt((uwnd - ucur)**2 + (vwnd - vcur)**2)
        mask = (nwr == 0)
        nwr[mask] = numpy.nan
        uwnd[mask] = numpy.nan
        vwnd[mask] = numpy.nan
        ucur[mask] = numpy.nan
        vcur[mask] = numpy.nan
        wrd = numpy.pi / 2 - numpy.arctan2(vwnd - vcur, uwnd - ucur)
        mssshort = numpy.log(nwr + 0.7) * 0.009
        # Replace nan values by 0 to avoid a runtime warning (nan values
        # will be restored afterwards)
        mssshort_nanmask_ind = numpy.where(numpy.isnan(mssshort))
        mssshort[mssshort_nanmask_ind] = 0
        mssshort[mssshort < 0] = 0
        mssshort[mssshort_nanmask_ind] = numpy.nan  # restore nan values
        # Directionality for short wave mss (if 0.5: isotrophic)
        facssdw = 0.6
        mssds = facssdw * mssshort
        msscs = mssshort - mssds
        mssxs = msscs * numpy.sin(wrd)**2 + mssds * numpy.cos(wrd)**2
        mssys = mssds * numpy.sin(wrd)**2 + msscs * numpy.cos(wrd)**2
        mssxys = abs(mssds - msscs) * numpy.sin(2* wrd)
        input_var['mssx'] = mssxs + mssxl
        input_var['mssy'] = mssys + mssyl
        input_var['mssxy'] = mssxys + mssxyl
        R2 = 0.55
        rbeam_angle = numpy.deg2rad(beam_angle)
        mssx = mssxs + mssxl
        mssy = mssys + mssyl
        mssxy = mssxys + mssxyl
        mask = ((mssx == 0) | (mssy == 0))
        mssx[mask] = numpy.nan
        mssy[mask] = numpy.nan
        expo = (-0.5 * numpy.tan(rbeam_angle)**2 * (numpy.cos(radial_angle)**2
                *mssy + numpy.sin(radial_angle)**2 * mssx
                - numpy.sin(2 * radial_angle) * mssxy) / (mssx * mssy))
        coeff = R2 / (2 * numpy.cos(rbeam_angle)**4 * numpy.sqrt(mssx * mssy))
        sigma_water = coeff * numpy.exp(expo)
        del input_var['mssd']
        return sigma_water


def compute_sigma(output_var_i, beam_angle, radial_angle, p):
    sigma_water = compute_sigma_water(output_var_i, beam_angle, radial_angle)
    if (p.ice is True) and ('ice' in output_var_i.keys()):
        if beam_angle == 6:
            sigma_ice = 2.5
        elif beam_angle == 12:
            sigma_ice = 1
        else:
            logger.info('beam angle is {} but should be either 6 or 12, '
                        'sigma_ice is set to 1'.format(beam_angle))
        mask = (numpy.isnan(output_var_i['ice']))
        c_ice = numpy.ma.MaskedArray(output_var_i['ice'], mask=mask)
        c_ice[c_ice.mask] = 0
        sigma0 = (1 - c_ice) * sigma_water + c_ice * sigma_ice
    else:
        sigma0 = sigma_water
    return sigma0

def save_SKIM(cycle, sgrid, time, outdata, p):
    file_output = '{}_c{:02d}_p{:03d}.nc'.format(p.file_output, cycle + 1,
                                                 sgrid.ipass)
    OutputSKIM = rw_data.Sat_SKIM(ifile=file_output, lon=sgrid.lon,
                                  lat=sgrid.lat, time=time,
                                  x_al=sgrid.x_al, cycle=sgrid.cycle)
    OutputSKIM.gridfile = sgrid.gridfile
    OutputSKIM.ipass = sgrid.ipass
    OutputSKIM.ncycle = sgrid.ncycle

    OutputSKIM.write_data(p, outdata) #ur_model=ur_model, index=vindice,
                       #   uss_err=err_uss, ur_uss=ur_uss, std_uss=std_uss,
                       #   nadir_err=[err.nadir, ], ur_obs=ur_obs,
                       #   instr=err_instr, u_model=u_model, v_model=v_model,
                       #   errdcos=errdcos)
    return None
