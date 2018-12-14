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

# - Define global variables for progress bars


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
    al_cycle = 0
    timeshift = 0
    sgrid.load_swath(p, cycle=cycle, x_al=x_al, x_al_nadir=[],
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
        if len(_filename) > 1:
            filename_u = os.path.join(p.indatadir, _filename[0])
            filename_v = os.path.join(p.indatadir, _filename[1])
        else:
            filename_u = os.path.join(p.indatadir, _filename[0])
            filename_v = os.path.join(p.indatadir, _filename[0])
        #filename = os.path.join(p.indatadir, list_file[0])
        model_data = model_data_ctor(p, ifile=(filename_u, filename_v),
                                     lonu=p.lonu, lonv=p.lonv, latu=p.latu,
                                     latv=p.latv)
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
        if Teval is None:
            Teval = numpy.zeros(numpy.shape(lon_out))
            if ind_out1[0].any():
                _tmp = interp(lat_in, lon_in[ind_in1],
                              numpy.isnan(var[:, ind_in1[0]]), kx=1, ky=1,
                              s=0)
                Teval[ind_out1] = _tmp.ev(lat_out[ind_out1], lon_out[ind_out1])
            if ind_out2[0].any():
                _tmp = interp(lat_in, lon_in[ind_in2],
                              numpy.isnan(var[:, ind_in2[0]]), kx=1, ky=1,
                              s=0)
                Teval[ind_out2] = _tmp.ev(lat_out[ind_out2], lon_out[ind_out2])
        # Trick to avoid nan in interpolation
        var_mask = + var
        var_mask[numpy.isnan(var_mask)] = 0.
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
        # Interpolate mask if it has not been done (Teval is None)
        interp = interpolate.RectBivariateSpline
        if Teval is None:
            _Teval = interp(lat_in, lon_in, numpy.isnan(var), kx=1, ky=1, s=0)
            Teval = _Teval.ev(lat_out, lon_out)
        # Trick to avoid nan in interpolation
        var_mask = + var
        var_mask[numpy.isnan(var_mask)] = 0.
        # Interpolate variable
        _var_out = interp(lat_in, lon_in, var_mask, kx=1, ky=1, s=0)
        var_out = _var_out.ev(lat_out, lon_out)
    # Mask variable with Teval
    var_out[Teval > 0] = numpy.nan
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
                        p, progress_bar=True):
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
        # meshgrid in 2D for interpolation purposes
        if p.grid == 'irregular':
            lon2Du, lat2Du = numpy.meshgrid(model_data.vlonu, model_data.vlatu)
            if (p.lonu == p.lonv) and (p.latu == p.latv):
                lon2Dv, lat2Dv = lon2Du, lat2Du
            else:
                lon2Dv, lat2Dv = numpy.meshgrid(model_data.vlonv,
                                                model_data.vlatv)
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
            filetime = (ifile - time_offset)%p.dim_time
            # read next file when the time dimension is reached
            if filetime >= (time_offset + p.dim_time):
                time_offset += p.dim_time
                nfile += 1
                filetime = (ifile - time_offset)%p.dim_time
            nfile = int(ifile /p.dim_time)
            filetime = ifile - nfile * p.dim_time

            _tmpfilename = list_file[nfile].split(',')
            if len(_tmpfilename) > 1:
                filename_u = os.path.join(p.indatadir, _tmpfilename[0])
                filename_v = os.path.join(p.indatadir, _tmpfilename[1])
            else:
                filename_u = os.path.join(p.indatadir, _tmpfilename[0])
                filename_v = os.path.join(p.indatadir, _tmpfilename[0])

            model_step_ctor = getattr(rw_data, model_data.model)
            model_step = model_step_ctor(p, ifile=(filename_u, filename_v),
                                         list_input_var=p.list_input_var,
                                         time=filetime)
            input_var_i = {}
            if p.grid == 'regular':
                model_step.read_var(p)
                if p.instr is True:
                    model_step.compute_mss()
                for key in model_step.input_var.keys():
                    _indlat = model_data.model_index_latu
                    _tmp = model_step.input_var[key][_indlat, :]
                    input_var_i[key] = +_tmp[:, model_data.model_index_lonu]

            else:
                model_step.read_var(p, index=None)
                if p.instr is True:
                    model_step.compute_mss()
                for key in model_step.input_var.keys():
                    _ind = model_data.model_indexu
                    input_var_i[key] = + model_step.input_var[key][_ind]
            # - Interpolate Model data on a SKIM grid and/or along the
            #   nadir track
            # if grid is regular, use interpolate.RectBivariateSpline to
            # interpolate
            if p.grid == 'regular' and \
                    len(numpy.shape(model_data.vlonu)) == 1:
                # Flatten satellite grid and select part of the track
                # corresponding to the model time
                for key in model_step.input_var.keys():
                    _tmp, Teval_u = interpolate_regular_1D(
                                                         model_data.vlonu,
                                                         model_data.vlatu,
                                                         input_var_i[key],
                                                         lon[ind_time[0]],
                                                         lat[ind_time[0]])
                    output_var_i[key][ind_time[0]] = + _tmp
            else:
                # Grid is irregular, interpolation can be done using
                # pyresample module if it is installed or griddata
                # function from scipy.
                # Note that griddata is slower than pyresample functions.
                lonuravel = + model_data.vlonu.ravel()
                laturavel = + model_data.vlatu.ravel()
                if p.lonu == p.lonv:
                    lonvravel = lonuravel
                    latvravel = laturavel
                else:
                    lonvravel = lon2Dv.ravel()
                    latvravel = lat2Dv.ravel()

                try:
                    import pyresample as pr
                    model_data.vlonu = pr.utils.wrap_longitudes(
                                                              model_data.vlonu)
                    lon = pr.utils.wrap_longitudes(lon)
                    if len(numpy.shape(model_data.vlonu)) <= 1:
                        model_data.vlonu = lon2Du
                        model_data.vlatu = lat2Du
                        if p.lonu == p.lonv:
                            model_data.vlonv = model_data.vlonu
                            model_data.vlatv = model_data.vlatu
                        else:
                            model_data.vlonv = lon2Dv
                            model_data.vlatv = lat2Dv

                    swath_defu = pr.geometry.SwathDefinition(
                                 lons=model_data.vlonu, lats=model_data.vlatu)
                    swath_defv = pr.geometry.SwathDefinition(
                                 lons=model_data.vlonv, lats=model_data.vlatv)
                    grid_def = pr.geometry.SwathDefinition(lons=lon,
                                                           lats=lat)
                    for key in model_step.input_var.keys():
                        _tmp = interpolate_irregular_pyresample(
                                          swath_defu, input_var_i[key],
                                          grid_def,
                                          p.posting,
                                          interp_type=p.interpolation)
                        output_var_i[key][ind_time[0]] = + _tmp
                except ImportError:
                    for key in model_step.input_var.keys():
                         _tmp = + input_var_i[key].ravel()
                    interp = interpolate.griddata((lonuravel, laturavel),
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


def compute_beam_noise_skim(p, output_var_i, radial_angle, beam_angle):
    output_var_i['ur_true'] = mod_tools.proj_radial(output_var_i['ucur'],
                                                    output_var_i['vcur'],
                                                    radial_angle)
    output_var_i['ur_obs'] = + output_var_i['ur_true']
    output_var_i['radial_angle'] = radial_angle
    if p.instr is True:
        # Compute sigma0:
        R2 = 0.5
        rbeam_angle = numpy.deg2rad(beam_angle)
        mssx = output_var_i['mssx']
        mssy = output_var_i['mssy']
        mssxy = output_var_i['mssxy']
        expo = (-0.5 * numpy.tan(rbeam_angle)**2 * (numpy.cos(radial_angle)**2
                *mssy + numpy.sin(radial_angle)**2 *mssx
                - numpy.sin(2 * radial_angle) * mssxy))
        coeff = R2 / (2 * numpy.cos(rbeam_angle)**4 * numpy.sqrt(mssx * mssy))
        sigma_water = coeff * numpy.exp(expo)
        if p.ice is True:
            if beam_angle == 6:
                sigma_ice = 2.5
            elif beam_angle == 12:
                sigma_ice = 1
            else:
                logger.info('beam angle is {} but should be either 6 or 12, '
                            'sigma_ice is set to 1'.format(beam_angle))
            c_ice = output_var_i['ice']
            sigma0 = (1 - c_ice) * sigma_water + c_ice * sigma_ice
        else:
            sigma0 = sigma_water
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
