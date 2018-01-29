'''Main program:
Usage: run_simulator(file_param)  \n
If no param file is specified, the default one is exemple/params_exemple.txt \n
In the first part of the program, model coordinates are read and the
SKIM swath is computing accordingly. \n
The SKIM grid parameters are saved in netcdf files, if you don't want to
recompute them, set maksgrid (in params file) to False.\n

In the second part of the program, errors are computed on SKIM grid for
each pass, for each cycle. The error free velocity is the velocity interpolated
from the model at each timestep. Note that there is no temporal interpolation
between model files and thus if several files are used in the velocity
interpolation, some discontinuities may be visible. \n

OUTPUTS are netcdf files containing the requested errors, the error free
radial velocity and the radial velocity with errors. There is one file every
pass and every cycle.

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
This file is part of skimsimulator.

skimsimulator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

skimsimulator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with skimsimulator.  If not, see <http://www.gnu.org/licenses/>.
'''
import os
from scipy import interpolate
import numpy
import math
import glob
import sys
import logging
import skimsimulator.build_swath as build_swath
import skimsimulator.rw_data as rw_data
import skimsimulator.build_error as build_error
import skimsimulator.mod_tools as mod_tools
import skimsimulator.const as const

# Define logger level for debug purposes
logger = logging.getLogger(__name__)

# - Define global variables for progress bars
istep = 0
ntot = 1
ifile = 0


def run_simulator(p):
    '''Main routine to run simulator, input is the imported parameter file,
    no outputs are returned but netcdf grid and data files are written as well
    as a skimsimulator.output file to store all used parameter.
    '''
    # - Initialize some parameters values
    p.shift_lon = getattr(p, 'shift_lon', None)
    p.shift_time = getattr(p, 'p.shift_time', None)
    model = getattr(p, 'model', 'NETCDF_MODEL')
    p.model = model
    p.model_nan = getattr(p, 'model_nan', 0)
    p.vel_factor = getattr(p, 'vel_factor', 1.)
    p.nadir = getattr(p, 'nadir', True)
    p.grid = getattr(p, 'grid', 'regular')
    p.rms_instr_factor = getattr(p, 'rms_instr_factor', 1)
    p.cycle = getattr(p, 'cycle', 0.0368)
    p.formula = getattr(p, 'formula', False)
    p.footprint_std = getattr(p, 'footprint_std', 0)
    p.timeshift = getattr(p, 'timeshift', 0)
    if p.shift_time is None:
        p.timeshift = 0

    # - Progress bar variables are global
    global istep
    global ntot

    # - Read list of user model files """
    if p.file_input is not None:
        list_file = [line.strip() for line in open(p.file_input)]
    else:
        list_file = None
    if p.uss is True and p.input_uss is not None:
        list_file_uss = [line.strip() for line in open(p.input_uss)]
    else:
        list_file_uss = None
    # - Read model input coordinates '''
    # If a list of model files are specified, read model file coordinates
    if p.file_input:
        # model_data = eval('rw_data.' + model
        #                  + '(p, file=p.indatadir+os.sep+list_file[0])')
        model_data_ctor = getattr(rw_data, model)
        filename = os.path.join(p.indatadir, list_file[0])
        model_data = model_data_ctor(p, file=filename)
    # if no modelbox is specified (modelbox=None), the domain of the input
    # data is taken as a modelbox
    # coordinates from the region defined by modelbox are selected
    if p.modelbox is not None:
        modelbox = numpy.array(p.modelbox, dtype='float')
        # Use convert to 360 data
        modelbox[0] = (modelbox[0]+360) % 360
        if modelbox[1] != 360:
            modelbox[1] = (modelbox[1]+360) % 360
    else:
        if p.file_input is not None:
            modelbox = model_data.calc_box()
        else:
            logger.error('modelbox should be provided if no model file is'
                         'provided')
            sys.exit(1)
    # - Extract data on modelbox
    # TODO: do only this step if modelbox is defined? Do it later?
    if p.file_input is not None:
        model_data.read_coordinates()
        # Select model data in the region modelbox
        if p.grid == 'regular':
            if modelbox[0] < modelbox[1]:
                _tmp = numpy.where(((modelbox[0]-1) <= model_data.vlonu)
                                   & (model_data.vlonu <= (modelbox[1]+1)))[0]
                model_data.model_index_lonu = + _tmp
            else:
                _tmp = numpy.where(((modelbox[0]-1) <= model_data.vlonu)
                                   | (model_data.vlonu <= (modelbox[1]+1)))[0]
                model_data.model_index_lonu = + _tmp
            _tmp = numpy.where(((modelbox[2]-1) <= model_data.vlatu)
                               & (model_data.vlatu <= (modelbox[3]+1)))[0]
            model_data.model_index_latu = + _tmp
            model_data.vlonu = model_data.vlonu[model_data.model_index_lonu]
            model_data.vlatu = model_data.vlatu[model_data.model_index_latu]
            if p.lonu != p.lonv:
                if modelbox[0] < modelbox[1]:
                    _tmp = numpy.where(((modelbox[0]-1) <= model_data.vlonv)
                                       & (model_data.vlonv
                                          <= (modelbox[1]+1)))[0]
                    model_data.model_index_lonv = + _tmp
                else:
                    _tmp = numpy.where(((modelbox[0]-1) <= model_data.vlonv)
                                       | (model_data.vlonv
                                          <= (modelbox[1]+1)))[0]
                    model_data.model_index_lonv = + _tmp
            else:
                model_data.model_index_lonv = model_data.model_index_lonu
            if p.latu != p.latv:
                _tmp = numpy.where(((modelbox[2]-1) <= model_data.vlatv)
                                   & (model_data.vlatv <= (modelbox[3]+1)))[0]
                model_data.model_index_latv = + _tmp
            else:
                model_data.model_index_latv = model_data.model_index_latu
            model_data.vlonv = model_data.vlonv[model_data.model_index_lonv]
            model_data.vlatv = model_data.vlatv[model_data.model_index_latv]

        else:
            if modelbox[0] < modelbox[1]:
                _tmp = numpy.where(((modelbox[0]-1) <= model_data.vlonu)
                                   & (model_data.vlonu <= (modelbox[1]+1))
                                   & ((modelbox[2]-1) <= model_data.vlatu)
                                   & (model_data.vlatu <= (modelbox[3]+1)))
                model_data.model_indexu = + _tmp
                model_data.model_indexv = model_data.model_indexu
                if p.lonu != p.lonv:
                    _tmp = numpy.where(((modelbox[0]-1) <= model_data.vlonv)
                                       & (model_data.vlonv <= (modelbox[1]+1))
                                       & ((modelbox[2]-1) <= model_data.vlatv)
                                       & (model_data.vlatv <= (modelbox[3]+1)))
                    model_data.model_indexv = + _tmp
            else:
                _tmp = numpy.where(((modelbox[0]-1) <= model_data.vlonu)
                                   | (model_data.vlonu <= (modelbox[1]+1))
                                   & ((modelbox[2]-1) <= model_data.vlatu)
                                   & (model_data.vlatu <= (modelbox[3]+1)))
                model_data.model_indexu = + _tmp
                model_data.model_indexv = model_data.model_indexu
                if p.lonu != p.lonv:
                    _tmp = numpy.where(((modelbox[0]-1) <= model_data.vlonv)
                                       | (model_data.vlonv <= (modelbox[1]+1))
                                       & ((modelbox[2]-1) <= model_data.vlatv)
                                       & (model_data.vlatv <= (modelbox[3]+1)))
                    model_data.model_indexv = + _tmp

        model_data.model = model
        # prevent IDL issues
        model_data.vloncircu = numpy.rad2deg(numpy.unwrap(model_data.vlonu))
        model_data.vloncircv = numpy.rad2deg(numpy.unwrap(model_data.vlonv))
        # If corrdinates are 1D and local std needs to be computed for uss bias
        # Grid coordinates in 2D
        if (p.uss is True and p.footprint_std is not None
             and p.footprint_std != 0):
            if len(numpy.shape(model_data.vlonu)) == 1:
                model_data.lon2D, model_data.lat2D = numpy.meshgrid(
                                                        model_data.vlonu,
                                                        model_data.vlatu)
            else:
                model_data.lon2D = model_data.vlonu
                model_data.lat2D = model_data.vlatu

    # avoid issue with 0=360 for global modelbox
    if modelbox[1] == 0:
        modelbox[1] = 359.99
    # - Make SKIM grid if necessary """
    if p.makesgrid is True:
        logger.info('\n Force creation of SKIM grid')
        # make nadir orbit
        orb = build_swath.makeorbit(modelbox, p, orbitfile=p.filesat)
        # build swath around this orbit
        build_swath.orbit2swath(modelbox, p, orb)
        logger.info("\n SKIM Grids and nadir tracks have been written in"
                    "{}".format(p.outdatadir))
        logger.info("-----------------------------------------------")

    # Initialize errors
    err, errnad = load_error(p)

    # - Compute interpolated velocity and errors for each pass, at each
    #   cycle
    logger.info('Compute interpolated velocity and errors:')
    #   load all SKIM grid files (one for each pass)
    listsgridfile = sorted(glob.glob(p.filesgrid + '_p*.nc'))
    if not listsgridfile:
        logger.error('\n There is no SKIM grid file in {}, run simulator with'
                     'option makesgrid set to true in your params'
                     'file'.format(p.outdatadir))
        sys.exit(1)
    # Build model time steps from parameter file
    modeltime = numpy.arange(0, p.nstep*p.timestep, p.timestep)
    #   Remove the grid from the list of model files
    if p.file_input:
        list_file.remove(list_file[0])
    #   Initialize progress bar variables
    istep = 0
    ntot = 1

    # - Loop on SKIM grid files
    for sgridfile in listsgridfile:
        #   Load SKIM grid files (Swath and nadir)
        sgrid = load_sgrid(sgridfile, p)
        # duplicate SKIM grids to assure that data are not modified and are
        # saved properly
        sgrid_tmp = load_sgrid(sgridfile, p)
        sgrid.gridfile = sgridfile
        # Convert cycle in days
        sgrid.cycle /= 86400.
        sgrid_tmp.cycle /= 86400.
        sgrid_tmp.indi = None
        # Select model data around the swath to reduce interpolation cost in
        # griddata ### TODO comment to be removed?
        # - Generate SKIM like data
        # Compute number of cycles needed to cover all nstep model timesteps
        rcycle = (p.timestep * p.nstep)/float(sgrid.cycle)
        ncycle = int(rcycle)
        #  Loop on all cycles
        for cycle in range(0, ncycle+1):
            if ifile > (p.nstep*p.timestep + 1):
                break
            # Create SKIM-like data
            if p.file_input is None:
                model_data = []
            # Initialize all list of variables (each beam is appended to the
            # list)
            # Initialize velocity, indices and mask to empty lists
            ur_true_all = []
            u_true_all = []
            v_true_all = []
            vindice_all = []
            ur_obs = []
            mask = []
            # Initialize noise to None if this noise is not computed
            # (so that the variable is not written in the netcdf)
            err_instr = None
            err_uss = None
            std_uss = None
            ur_uss = None
            errdcos_tot = None
            # Initialize noises to empty lists if the noise is set to True
            if p.instr is True:
                err_instr = []
            if p.uss is True:
                err_uss = []
                std_uss = []
                ur_uss = []
                errdcos_tot = []
            # Loop over the beams
            for i in range(len(p.list_pos) + 1):
                sgrid_tmp.lon = sgrid.lon[i]
                sgrid_tmp.lat = sgrid.lat[i]
                sgrid_tmp.time = sgrid.time[i]
                # If nadir, compute a different noise, not implemented so far,
                # Initialize at zero.
                if i == 0:
                    ur_true = numpy.full(numpy.shape(sgrid_tmp.lon), numpy.nan)
                    u_true = numpy.full(numpy.shape(sgrid_tmp.lon), numpy.nan)
                    v_true = numpy.full(numpy.shape(sgrid_tmp.lon), numpy.nan)
                    vindice = numpy.full(numpy.shape(sgrid_tmp.lon), numpy.nan)
                    err.ur_obs = numpy.full(numpy.shape(sgrid_tmp.lon),
                                            numpy.nan)
                    err.instr = numpy.full(numpy.shape(sgrid_tmp.lon),
                                           numpy.nan)
                    mask_tmp = numpy.full(numpy.shape(sgrid_tmp.lon),
                                          numpy.nan)
                    if p.uss is True:
                        err.ur_uss = numpy.full(numpy.shape(sgrid_tmp.lon),
                                                numpy.nan)
                        err.err_uss = numpy.full(numpy.shape(sgrid_tmp.lon),
                                                 numpy.nan)
                        err.std_uss = numpy.full(numpy.shape(sgrid_tmp.lon),
                                                 numpy.nan)
                    # If the footprint of the std is defined, find the cycles
                    # loacated in the area of the footprint to compute the std
                    # later. Only use when the approximation with errdcos is
                    # used.
                    if (p.uss is True and p.footprint_std is not None
                          and p.footprint_std != 0 and sgrid_tmp.indi is None
                          and p.formula is True):
                        sgrid_tmp.indi = numpy.zeros((numpy.shape(
                                                     sgrid_tmp.lon)[0], 2))
                        sgrid_tmp.indj = numpy.zeros((numpy.shape(
                                                     sgrid_tmp.lon)[0], 2))
                        for ib in range(len(sgrid_tmp.lon)):
                            ilon = sgrid_tmp.lon[ib]
                            ilat = sgrid_tmp.lat[ib]
                            _dist = (((model_data.lon2D-ilon)
                                     * numpy.cos(model_data.lat2D))**2
                                     + (model_data.lat2D - ilat)**2)
                            footprint_std = p.footprint_std / const.deg2km
                            indi, indj = numpy.where(_dist < footprint_std)
                            if indi.any() and indj.any():
                                sgrid_tmp.indi[ib, 0] = indi[0]
                                sgrid_tmp.indi[ib, 1] = indi[-1]
                                sgrid_tmp.indj[ib, 0] = indj[0]
                                sgrid_tmp.indj[ib, 1] = indj[-1]
                            else:
                                sgrid_tmp.indi[ib-1, :]

                                sgrid_tmp.indi[ib, 0] = sgrid_tmp.indi[ib-1, 0]
                                sgrid_tmp.indi[ib, 1] = sgrid_tmp.indi[ib-1, 1]
                                sgrid_tmp.indj[ib, 0] = sgrid_tmp.indj[ib-1, 0]
                                sgrid_tmp.indj[ib, 1] = sgrid_tmp.indj[ib-1, 1]
                # Interpolate the velocity and compute the noise for each beam
                else:
                    # Stoke drift bias impact factor
                    Gvar = p.G[i - 1]
                    # Instrument noise file
                    rms_instr = p.rms_instr[i - 1]
                    # Geometrical error errdcos
                    errdcos = 1
                    # If the formula is used and errdcos is constant for each
                    # beam. Rough temporary approximation, should not be used.
                    if p.errdcos is not None and p.uss is True:
                        errdcos = p.errdcos[i - 1]
                    # Read radial angle for projection on lon, lat reference
                    radial_angle = sgrid.radial_angle[:, i - 1]
                    ##############################
                    # Compute SKIM like data data
                    shape_all = (numpy.shape(listsgridfile)[0] * rcycle
                                 * (len(p.list_pos) + 1))
                    create = create_SKIMlikedata(cycle, shape_all, list_file,
                                                 list_file_uss, modelbox,
                                                 sgrid_tmp, model_data,
                                                 modeltime, err, Gvar,
                                                 rms_instr, errdcos,
                                                 radial_angle, p,
                                                 progress_bar=True)
                    ur_true, u_true, v_true, vindice, time, progress = create
                    mask_tmp = numpy.isnan(err.ur_uss)
                # Append variables for each beam
                if p.instr is True:
                    err_instr.append(err.instr)
                if p.uss is True:
                    ur_uss.append(err.ur_uss)
                    std_uss.append(err.std_uss)
                    err_uss.append(err.err_uss)
                ur_true_all.append(ur_true)
                u_true_all.append(u_true)
                v_true_all.append(v_true)
                vindice_all.append(vindice)
                mask.append(mask_tmp)
                # ur_obs.append(err.ur_obs)
            # Compute uss bias
            #   Compute errdcos if Formula is True
            if p.uss is True and p.errdcos is None and p.formula is True:
                errdcos_tot, err_uss2 = compute_errdcos(p, sgrid, mask,
                                                        err_uss)
            # Compute directly bias if formula is False
            if p.uss is True and p.formula is False:
                err_uss2 = compute_errussr(p, sgrid, mask, ur_uss, err_uss)
                # errdcos_tot.append(errdcos)
                # err_uss.append(err.err_uss)
                # , err_uss)
            for i in range(1, len(p.list_pos) + 1):
                make_err = build_error.make_vel_error
                ur_obs_i = make_err(ur_true_all[i], p, instr=err_instr[i],
                                    err_uss=err_uss2[i])
                ur_obs.append(ur_obs_i)
            #   Save outputs in a netcdf file
            if ((~numpy.isnan(numpy.array(vindice_all))).any()
                  or not p.file_input):
                sgrid.ncycle = cycle
                save_SKIM(cycle, sgrid, err, p, time=time, vindice=vindice_all,
                          ur_model=ur_true_all, ur_obs=ur_obs, std_uss=std_uss,
                          err_instr=err_instr, ur_uss=ur_uss, err_uss=err_uss2,
                          u_model=u_true_all, v_model=v_true_all,
                          errdcos=errdcos_tot)
            del time
            # if p.file_input: del index
        if p.file_input:
            model_data.vlonu = (model_data.vlonu + 360) % 360
            model_data.vlonv = (model_data.vlonv + 360) % 360

        modelbox[0] = (modelbox[0] + 360) % 360
        modelbox[1] = (modelbox[1] + 360) % 360
        del sgrid
    if progress != 1:
        progress = mod_tools.update_progress(1,
                                             'All passes have been processed',
                                             '')
    # - Write Selected parameters in a txt file
    rw_data.write_params(p, os.path.join(p.outdatadir,
                                         'skim_simulator.output'))
    logger.info("\n Simulated skim files have been written in "
                "{}".format(p.outdatadir))
    logger.info("----------------------------------------------------------")


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
    sgrid.load_swath(p, cycle=cycle, x_al=x_al, al_cycle=al_cycle,
                     timeshift=timeshift)
    # sgrid.loncirc = []
    # for i in range(len(sgrid.lon)):
    #     sgrid.loncirc.append(numpy.rad2deg(numpy.unwrap(sgrid.lon[i])))
    # Extract the pass number from the file name
    ipass = int(sgridfile[-6: -3])
    sgrid.ipass = ipass
    return sgrid


def interpolate_regular_1D(lon_in, lat_in, var, lon_out, lat_out, Teval=None):
    ''' Interpolation of data when grid is regular and coordinate in 1D. '''
    # To correct for IDL issues
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
        radius_n = radius * 10**3
        var_out = interp(swath_in, var, grid_out, radius_of_influence=radius_n,
                         epsilon=100)
    else:
        radius_g = radius * 3 * 10**3
        sigma_g = radius * 10**3
        var_out = interp(swath_in, var, grid_out, radius_of_influence=radius_g,
                         sigmas=sigma_g, fill_value=None)
    return var_out


def create_SKIMlikedata(cycle, ntotfile, list_file, list_file_uss, modelbox,
                        sgrid, model_data, modeltime, err, Gvar, rms_instr,
                        errdcos,
                        radial_angle, p, progress_bar=True):
    '''Create SKIM and nadir errors err and errnad, interpolate model velocity\
    model_data on swath and nadir track,
    compute SKIM-like and nadir-like data for cycle, SKIM grid sgrid and
    ngrid. '''
    # - Progress bar variables are global
    global istep
    global ntot
    #   Initialiaze errors and velocity
    progress = 0
    err.instr = numpy.full(numpy.shape(sgrid.lon)[0], numpy.nan)
    err.ur_uss = numpy.full(numpy.shape(sgrid.lon)[0], numpy.nan)
    std_uss = numpy.full(numpy.shape(sgrid.lon)[0], numpy.nan)
    err.wet_tropo1nadir = numpy.full(numpy.shape(sgrid.lon)[0], numpy.nan)
    err.wet_tropo2nadir = numpy.full(numpy.shape(sgrid.lon)[0], numpy.nan)
    err.wtnadir = numpy.full(numpy.shape(sgrid.lon)[0], numpy.nan)
    err.nadir = numpy.full(numpy.shape(sgrid.lon)[0], numpy.nan)
    date1 = cycle * sgrid.cycle
    u_true = numpy.full(numpy.shape(sgrid.lon)[0], numpy.nan)
    v_true = numpy.full(numpy.shape(sgrid.lon)[0], numpy.nan)
    u_uss = numpy.full(numpy.shape(sgrid.lon)[0], numpy.nan)
    v_uss = numpy.full(numpy.shape(sgrid.lon)[0], numpy.nan)
    ur_true = numpy.full(numpy.shape(sgrid.lon)[0], numpy.nan)
    vindice = numpy.full(numpy.shape(sgrid.lon)[0], numpy.nan)
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
            perc = float(istep)/float(ntot * ntotfile)
            strpass = 'pass: {}'.format(sgrid.ipass)
            strcycle = 'model file: {}, cycle: {}'.format(ifile, cycle + 1)
            progress = mod_tools.update_progress(perc, strpass, strcycle)
            # If there are satellite data, Get true velcoity from model
            if numpy.shape(index_filemodel)[1] > 0:
                # number of file to be processed used in the progress bar
                # ntot = ntot + numpy.shape(index_filemodel)[1]-1
                ntot = numpy.shape(index_filemodel)[1]
                # if numpy.shape(index)[1]>1:
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
            if model_data.model == 'WW3':
                filetime = ifile - time_offset
                # read next file when the time dimension is reached
                if filetime >= (time_offset + p.dim_time[nfile]):
                    time_offset += p.dim_time[nfile]
                    nfile += 1
                    filetime = ifile - time_offset
                filename = os.path.join(p.indatadir, list_file[nfile])
                model_step = rw_data.WW3(p, file=filename, varu=p.varu,
                                         varv=p.varv, time=filetime)
                if p.uss is True:
                    filename = os.path.join(p.indatadir, list_file_uss[nfile])
                    uss_step = rw_data.WW3(p, file=filename, varu='uuss',
                                           varv='vuss', time=filetime)
            else:
                model_step_ctor = getattr(rw_data, model_data.model)
                model_step = model_step_ctor(file=os.path.join(p.indatadir,
                                             list_file[ifile]), varu=p.varu,
                                             varv=p.varv)
                if p.uss is True:
                    model_step = model_step_ctor(file=os.path.join(p.indatadir,
                                                 list_file_uss[ifile]),
                                                 varu='uuss', varv='vuss')
            if p.grid == 'regular':
                model_step.read_var()
                u_model_tmp = model_step.vvaru[model_data.model_index_latu, :]
                u_model = u_model_tmp[:, model_data.model_index_lonu]
                v_model_tmp = model_step.vvarv[model_data.model_index_latv, :]
                v_model = v_model_tmp[:, model_data.model_index_lonv]
                if p.uss is True:
                    uss_step.read_var()
                    u_uss_mod_tmp = uss_step.vvaru[model_data.model_index_latu,
                                                   :]
                    u_uss_mod = u_uss_mod_tmp[:, model_data.model_index_lonu]
                    v_uss_mod_tmp = uss_step.vvaru[model_data.model_index_latv,
                                                   :]
                    v_uss_mod = v_uss_mod_tmp[:, model_data.model_index_lonv]
            else:
                model_step.read_var(index=None)
                u_model = model_step.vvaru[model_data.model_indexu]
                v_model = model_step.vvarv[model_data.model_indexv]
                if p.uss is True:
                    uss_step.read_var(index=None)
                    u_uss_mod = uss_step.vvaru[model_data.model_indexu]
                    v_uss_mod = uss_step.vvarv[model_data.model_indexv]
            # - Interpolate Model data on a SKIM grid and/or along the
            #   nadir track
            # if grid is regular, use interpolate.RectBivariateSpline to
            # interpolate
            if p.grid == 'regular' and \
                    len(numpy.shape(model_data.vlonu)) == 1:
                # Flatten satellite grid and select part of the track
                # corresponding to the model time
                u_true_ind_time, Teval_u = interpolate_regular_1D(
                                                         model_data.vlonu,
                                                         model_data.vlatu,
                                                         u_model,
                                                         lon[ind_time[0]],
                                                         lat[ind_time[0]])
                u_true[ind_time[0]] = + u_true_ind_time
                if p.uss is True:
                    u_uss_ind_time, Teval = interpolate_regular_1D(
                                                         model_data.vlonu,
                                                         model_data.vlatu,
                                                         u_uss_mod,
                                                         lon[ind_time[0]],
                                                         lat[ind_time[0]],
                                                         Teval=Teval_u)
                    u_uss[ind_time[0]] = + u_uss_ind_time
                v_true_ind_time, Teval_v = interpolate_regular_1D(
                                                         model_data.vlonv,
                                                         model_data.vlatv,
                                                         v_model,
                                                         lon[ind_time[0]],
                                                         lat[ind_time[0]])
                v_true[ind_time[0]] = + v_true_ind_time
                if p.uss is True:
                    v_uss_ind_time, Teval = interpolate_regular_1D(
                                                         model_data.vlonu,
                                                         model_data.vlatu,
                                                         v_uss_mod,
                                                         lon[ind_time[0]],
                                                         lat[ind_time[0]],
                                                         Teval=Teval_v)
                    v_uss[ind_time[0]] = + v_uss_ind_time
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
                    u_true[ind_time[0]] = interpolate_irregular_pyresample(
                                          swath_defu, u_model, grid_def,
                                          p.delta_al,
                                          interp_type=p.interpolation)
                    v_true[ind_time[0]] = interpolate_irregular_pyresample(
                                          swath_defv, v_model, grid_def,
                                          p.delta_al,
                                          interp_type=p.interpolation)
                    if p.uss is True:
                        u_uss[ind_time[0]] = interpolate_irregular_pyresample(
                                          swath_defu, u_uss_mod, grid_def,
                                          p.delta_al,
                                          interp_type=p.interpolation)
                        v_uss[ind_time[0]] = interpolate_irregular_pyresample(
                                          swath_defu, v_uss_mod, grid_def,
                                          p.delta_al,
                                          interp_type=p.interpolation)
                except ImportError:
                    uravel = + u_model.ravel()
                    interp = interpolate.griddata((lonuravel, laturavel),
                                                  uravel, (lon[ind_time[0]],
                                                  lat[ind_time[0]]),
                                                  method=p.interpolation)
                    u_true[ind_time[0]] = interp
                    vravel = + u_model.ravel()
                    interp = interpolate.griddata((lonuravel, laturavel),
                                                  vravel, (lon[ind_time[0]],
                                                  lat[ind_time[0]]),
                                                  method=p.interpolation)
                    v_true[ind_time[0]] = interp
                    if p.uss is True:
                        u_ussravel = + u_uss_mod.ravel()
                        interp = interpolate.griddata((lonuravel, laturavel),
                                                      u_ussravel,
                                                      (lon[ind_time[0]],
                                                      lat[ind_time[0]]),
                                                      method=p.interpolation)
                        u_uss[ind_time[0]] = interp
                        v_ussravel = + v_uss_mod.ravel()
                        interp = interpolate.griddata((lonvravel, latvravel),
                                                      v_ussravel,
                                                      (lon[ind_time[0]],
                                                      lat[ind_time[0]]),
                                                      method=p.interpolation)
                        v_uss[ind_time[0]] = interp
            # Force value outside modelbox at nan
            if modelbox[0] > modelbox[1]:
                    u_true[numpy.where(~((lon > modelbox[0])
                                       | (lon < modelbox[1]))
                                       | (lat < modelbox[2])
                                       | (lat > modelbox[3]))] = numpy.nan
                    v_true[numpy.where(((lon > modelbox[0])
                                       & (lon < modelbox[1]))
                                       | (lat < modelbox[2])
                                       | (lat > modelbox[3]))] = numpy.nan
            else:
                    u_true[numpy.where((lon < modelbox[0])
                                       | (lon > modelbox[1])
                                       | (lat < modelbox[2])
                                       | (lat > modelbox[3]))] = numpy.nan
                    v_true[numpy.where((lon < modelbox[0])
                                       | (lon > modelbox[1])
                                       | (lat < modelbox[2])
                                       | (lat > modelbox[3]))] = numpy.nan
            vindice[ind_time[0]] = ifile
            # del u_true, v_true, model_step
            # Compute std of uss at each beam to compute the corrected bias of
            # the stoke drift
            if p.uss is True:
                if (p.footprint_std is not None and p.footprint_std != 0
                       and p.formula is True):
                    std_uss_local = numpy.full(numpy.shape(ind_time[0]),
                                               numpy.nan)
                    for ib in range(len(ind_time[0])):
                        indi = sgrid.indi[ind_time[0][ib], :]
                        indj = sgrid.indj[ind_time[0][ib], :]
                        if indi[0] != -1:
                            slice_i = slice(int(indi[0]), int(indi[1] + 1))
                            slice_j = slice(int(indj[0]), int(indj[1] + 1))
                            # TODO Change if grid is not the same for u and v
                            u_uss_local = u_uss_mod[slice_i, slice_j]
                            v_uss_local = v_uss_mod[slice_i, slice_j]
                            std_uss_local[ib] = numpy.nanstd(numpy.sqrt(
                                                             u_uss_local**2
                                                             + v_uss_local**2))
                        else:
                            std_uss_local[ib] = numpy.nan
                else:
                    std_uss_local = (numpy.ones(numpy.shape(ind_time[0]))
                                     * numpy.nanstd(numpy.sqrt(
                                     u_uss_mod**2 + v_uss_mod**2)))
            std_uss[ind_time[0]] = + std_uss_local
        istep += 1
    else:
        istep += 1
        progress = mod_tools.update_progress(float(istep)/float(ntotfile*ntot),
                                             'pass: {}'.format(sgrid.ipass),
                                             'no model file provided, cycle: '
                                             '{}'.format(cycle+1))
    ur_true = mod_tools.proj_radial(u_true, v_true, radial_angle)
    if p.uss is not True:
        u_uss = None
        v_uss = None
    err.make_error(ur_true, p, radial_angle, Gvar, rms_instr,
                   uss=(u_uss, v_uss), std_local=std_uss, errdcos=errdcos)
    return ur_true, u_true, v_true, vindice, time, progress


def compute_errdcos(p, sgrid, mask, err_uss):
    '''Compute the weigthed distance for the uss reconstruction formula. '''
    # Compute number of azimuth per mega cycle
    naz = (60 / (p.rotation_speed * p.cycle * len(p.list_pos)))
    # Initialize errdcos
    errdcos_tot = []
    errdcos_tot.append(numpy.full(numpy.shape(mask[0][:]), numpy.nan))
    # Convert list into arrays
    lon_array = numpy.transpose(numpy.array(sgrid.lon[1:][:]))
    lat_array = numpy.transpose(numpy.array(sgrid.lat[1:][:]))
    mask_array = numpy.transpose(numpy.array(mask[1:][:]))
    for i in range(1, len(p.list_pos) + 1):
        mask_minicycle = mask[i][:]
        radial_angle = sgrid.radial_angle[:, i - 1]
        N = len(mask_minicycle)
        dtheta = numpy.mean(abs(radial_angle[1:] - radial_angle[:-1]))
        errdcos = numpy.full(numpy.shape(radial_angle), numpy.nan)
        for ib in range(N):
            if mask_minicycle[ib] is True:
                continue
            theta = radial_angle[ib] % (2 * math.pi)
            errdcos[ib] = 0
            ntheta2 = 0
            for theta2 in numpy.arange(theta - math.pi/2,
                                       theta + math.pi/2, dtheta):
                theta2 = theta2 % (2 * math.pi)
                start_az = int(max(0, (ib - naz*2)))
                end_az = int(min(N, (ib + naz*2)))
                slice_az = slice(start_az, end_az)

                lon = lon_array[slice_az, :]
                lat = lat_array[slice_az, :]
                mask_ind = mask_array[slice_az, :]
                lon[numpy.where(mask_ind is True)] = -1.36*10**9
                lat[numpy.where(mask_ind is True)] = -1.36*10**9
                angle = sgrid.radial_angle[slice_az, :]
                angle = numpy.mod(angle, 2 * math.pi)
                ind_angle = numpy.where((angle >= (theta2 - dtheta/2.))
                                        & (angle < (theta2 + dtheta/2.)))
                # To be tested: ind_angle can be empty near coast?
                if len(ind_angle) == 0:
                    logger.debug('WTF')
                    continue
                lon = lon[ind_angle]
                lat = lat[ind_angle]
                dlon_km = ((lon - sgrid.lon[i][ib])*111
                           * numpy.cos(sgrid.lat[i][ib] * math.pi / 180.))
                dlat_km = (lat - sgrid.lat[i][ib])*111
                dist = numpy.sqrt(dlon_km**2 + dlat_km**2)
                if len(dist) > 0:
                    ind_dist = numpy.argmin(dist)
                    errdcos[ib] += (dist[ind_dist] * numpy.cos(theta
                                    - theta2))**2
                else:
                    errdcos[ib] = numpy.nan

                ntheta2 += 1
            errdcos[ib] /= ntheta2
        errdcos = numpy.sqrt(errdcos)
        err_uss[i][:] *= errdcos / 20
        errdcos_tot.append(errdcos)
    return errdcos_tot, err_uss


def compute_errussr(p, sgrid, mask, uss_r, err_uss):
    ''' Compute the reconstructed error of ussr using information from the uss
    calculated with closed beams. '''
    # Compute number of azimuth per mega cycle
    naz = (60 / (p.rotation_speed * p.cycle * len(p.list_pos)))
    # Initialize errdcos
    # Convert list into arrays
    lon_array = numpy.transpose(numpy.array(sgrid.lon[1:][:]))
    lat_array = numpy.transpose(numpy.array(sgrid.lat[1:][:]))
    uss_r_array = numpy.transpose(numpy.array(uss_r[1:][:]))
    mask_array = numpy.transpose(numpy.array(mask[1:][:]))
    for i in range(1, len(p.list_pos) + 1):
        mask_minicycle = mask[i][:]
        radial_angle = sgrid.radial_angle[:, i - 1]
        N = len(mask_minicycle)
        dtheta = numpy.mean(abs(radial_angle[1:] - radial_angle[:-1]))
        # errdcos = numpy.zeros(numpy.shape(radial_angle)) * numpy.nan
        errussr = numpy.full(numpy.shape(radial_angle), numpy.nan)
        for ib in range(N):
            if mask_minicycle[ib] is True:
                continue
            theta = radial_angle[ib] % (2 * math.pi)
            errussr[ib] = 0
            ntheta2 = 0
            for theta2 in numpy.arange(theta - math.pi/2,
                                       theta + math.pi/2, dtheta):
                theta2rad = theta2 % (2 * math.pi)
                start_az = int(max(0, (ib - naz*3)))
                end_az = int(min(N, (ib + naz*3)))
                slice_az = slice(start_az, end_az)

                lon = lon_array[slice_az, :]
                lat = lat_array[slice_az, :]
                uss_r_loc = uss_r_array[slice_az, :]
                mask_ind = mask_array[slice_az, :]
                lon[numpy.where(mask_ind is True)] = -1.36*10**9
                lat[numpy.where(mask_ind is True)] = -1.36*10**9
                angle = sgrid.radial_angle[slice_az, :]
                angle = numpy.mod(angle, 2 * math.pi)
                ind_angle = numpy.where((angle >= (theta2rad - dtheta))
                                        & (angle < (theta2rad + dtheta)))
                # To be tested: ind_angle can be empty near coast?
                if len(ind_angle) == 0:
                    logger.debug('no ind_angle found, pass is to small?')
                    continue
                lon = lon[ind_angle]
                lat = lat[ind_angle]
                if lon.size == 0:
                    continue
                uss_r_loc2 = uss_r_loc[ind_angle]
                sgridlon = sgrid.lon[i][ib]
                if (numpy.max(lon)) > 359 and (numpy.min(lon) < 1):
                    lon[lon > 180] = lon[lon > 180] - 360
                    lon = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(lon)))
                    if sgridlon > 180:
                        sgridlon = sgridlon - 360
                        sgridlon = numpy.rad2deg(numpy.deg2rad(sgridlon))
                dlon_km = ((lon - sgridlon)*111
                           * numpy.cos(sgrid.lat[i][ib] * math.pi / 180.))
                dlat_km = (lat - sgrid.lat[i][ib]) * 111
                dist = numpy.sqrt(dlon_km**2 + dlat_km**2)
                if len(dist) > 0:
                    ind_dist = numpy.argmin(dist)
                    errussr[ib] += (numpy.cos(theta
                                    - angle[ind_angle[0][ind_dist],
                                            ind_angle[1][ind_dist]])
                                    * uss_r_loc2[ind_dist] * dtheta
                                    / (math.pi) * 2)
                else:
                    errussr[ib] = numpy.nan
                ntheta2 += 1
        err_uss[i][:] = (errussr - uss_r[i][:])
    return err_uss


def save_SKIM(cycle, sgrid, err, p, time=(), vindice=(), ur_model=(),
              ur_obs=(), err_instr=(), ur_uss=(), err_uss=(),
              u_model=(), v_model=(), std_uss=(), errdcos=()):
    file_output = (p.file_output + '_c' + str(cycle+1).zfill(2) + '_p'
                   + str(sgrid.ipass).zfill(3) + '.nc')
    OutputSKIM = rw_data.Sat_SKIM(ifile=file_output, lon=sgrid.lon,
                                  lat=sgrid.lat, time=sgrid.time,
                                  x_al=sgrid.x_al, cycle=sgrid.cycle)
    OutputSKIM.gridfile = sgrid.gridfile
    OutputSKIM.ipass = sgrid.ipass
    OutputSKIM.ncycle = sgrid.ncycle
    OutputSKIM.write_data(p, ur_model=ur_model, index=vindice,
                          uss_err=err_uss, ur_uss=ur_uss, std_uss=std_uss,
                          nadir_err=[err.nadir, ], ur_obs=ur_obs,
                          instr=err_instr, u_model=u_model, v_model=v_model,
                          errdcos=errdcos)
    return None
