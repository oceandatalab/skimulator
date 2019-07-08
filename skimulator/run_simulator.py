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
# - written for Python 3.5, tested with Python 3.5, 3.7
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
import traceback
import pkg_resources
import skimulator.build_swath as build_swath
import skimulator.rw_data as rw_data
import skimulator.build_error as build_error
import skimulator.mod_tools as mod_tools
import skimulator.const as const
import skimulator.mod_run as mod
import skimulator.mod_parallel
import skimulator.mod_uwb_corr as mod_uwb_corr
# Define logger level for debug purposes
logger = logging.getLogger(__name__)
#logger = multiprocessing.log_to_stderr()
#logger.setLevel(logging.DEBUG)

# - Define global variables for progress bars
istep = 0
ntot = 1
ifile = 0


def run_simulator(p, die_on_error=False):
    '''Main routine to run simulator, input is the imported parameter file,
    no outputs are returned but netcdf grid and data files are written as well
    as a skimulator.output file to store all used parameter.
    '''
    # - Initialize some parameters values
    timestart = datetime.datetime.now()
    mod_tools.initialize_parameters(p)
    mod_tools.check_path(p)
    model = p.model
    '''
    p.timeshift = getattr(p, 'timeshift', 0)
    if p.shift_time is None:
        p.timeshift = 0
    '''

    # - Read list of user model files """
    if p.file_input is not None:
        model_data, list_file = mod.load_coordinate_model(p)
    # if no modelbox is specified (modelbox=None), the domain of the input
    # data is taken as a modelbox
    # coordinates from the region defined by modelbox are selected
    if p.modelbox is not None:
        modelbox = numpy.array(p.modelbox, dtype='float')
        # Use convert to 360 data
        modelbox[0] = (modelbox[0] + 360) % 360
        if modelbox[1] != 360:
            modelbox[1] = (modelbox[1] + 360) % 360
    else:
        if p.file_input is not None:
            modelbox = model_data.calc_box(p)
        else:
            logger.error('modelbox should be provided if no model file is'
                         'provided')
            sys.exit(1)
    # - Extract data on modelbox
    # TODO: do only this step if modelbox is defined? Do it later?
    if p.file_input is not None:
        model_data.read_coordinates(p)
        # Select model data in the region modelbox
        model_data.model_index_lon = {}
        model_data.model_index_lat = {}
        model_data.model_index = {}
        model_data.vloncirc = {}
        for key in model_data.vlon.keys():
            _lon = + model_data.vlon[key]
            _lat = + model_data.vlat[key]
            if p.grid == 'regular':
                if modelbox[0] < modelbox[1]:
                    _tmp = numpy.where(((modelbox[0]-1) <= _lon)
                                       & (_lon <= (modelbox[1]+1)))[0]
                    model_data.model_index_lon[key] = + _tmp
                else:
                    _tmp = numpy.where(((modelbox[0]-1) <= _lon)
                                       | (_lon <= (modelbox[1]+1)))[0]
                    model_data.model_index_lon[key] = + _tmp
                _tmp = numpy.where(((modelbox[2]-1) <= _lat)
                                   & (_lat <= (modelbox[3]+1)))[0]
                model_data.model_index_lat[key] = + _tmp
                model_data.vlon[key] = + _lon[model_data.model_index_lon[key]]
                model_data.vlat[key] = + _lat[model_data.model_index_lat[key]]
            else:
                if modelbox[0] < modelbox[1]:
                    _tmp = numpy.where(((modelbox[0]-1) <= _lon)
                                       & (_lon <= (modelbox[1]+1))
                                       & ((modelbox[2]-1) <= _lat)
                                       & (_lat <= (modelbox[3]+1)))
                    model_data.model_index[key] = + _tmp
                else:
                    _tmp = numpy.where(((modelbox[0]-1) <= _lon)
                                       | (_lon <= (modelbox[1]+1))
                                       & ((modelbox[2]-1) <= _lat)
                                       & (_lat <= (modelbox[3]+1)))
                    model_data.model_index[key] = + _tmp
            # prevent IDL issues
            _wlon = model_data.vlon[key]
            model_data.vloncirc[key] = numpy.rad2deg(numpy.unwrap(_wlon))

        model_data.model = model
    else:
        model_data = []

    # avoid issue with 0=360 for global modelbox
    if modelbox[1] == 0:
        modelbox[1] = 359.99
    # - Make SKIM grid if necessary """
    if p.makesgrid is True:
        logger.info('\n Force creation of SKIM grid')
        # make nadir orbit
        orb = build_swath.makeorbit(modelbox, p, orbitfile=p.filesat)
        # build swath around this orbit
        ok = False
        try:
            ok = build_swath.orbit2swath(modelbox, p, orb, die_on_error)
        except skimulator.mod_parallel.MultiprocessingError:
            logger.error('An error occurred with the multiprocessing '
                         'framework')
            traceback.print_exception(*sys.exc_info())
            sys.exit(1)
        except skimulator.mod_parallel.DyingOnError:
            logger.error('An error occurred and all errors are fatal')
            sys.exit(1)

        if not ok:
            logger.error('Errors occurred while generating grid files')
            sys.exit(1)

        logger.info("\n SKIM Grids and nadir tracks have been written in"
                    "{}".format(p.outdatadir))
        logger.info("-----------------------------------------------")

    # Initialize errors
    # err, errnad = mod.load_error(p)

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
    if p.file_input and p.file_grid_model is None:
        logger.info("WARNING: the first file is not used to build data")
        list_file.remove(list_file[0])
    # Add compulsary key # TODO proof that
    list_compulsary_key = ['ur_obs', 'ur_true', 'radial_angle', 'ussr', 'uwd',
                           'dsigma', 'ussr_est']
    for key in list_compulsary_key:
        if key not in p.list_output:
            p.list_output.append(key)


    # - Loop on SKIM grid files
    jobs = []
    p2 = mod_tools.todict(p)
    time_yaw = None
    vac_yaw = None
    if p.attitude is True and os.path.isfile(p.yaw_file):
        time_yaw, vac_yaw = build_error.load_yaw_aocs(p.yaw_file)
        # time_yaw = time_yaw / 86400
    for sgridfile in listsgridfile:
        jobs.append([sgridfile, p2, listsgridfile, list_file,
                     modelbox, model_data, modeltime, time_yaw, vac_yaw])
    ok = False
    try:
        ok = make_skim_data(p.proc_count, jobs, die_on_error, p.progress_bar)
    except skimulator.mod_parallel.MultiprocessingError:
        logger.error('An error occurred with the multiprocessing framework')
        traceback.print_exception(*sys.exc_info())
        sys.exit(1)
    except skimulator.mod_parallel.DyingOnError:
        logger.error('An error occurred and all errors are fatal')
        sys.exit(1)

    # - Write Selected parameters in a txt file
    timestop = datetime.datetime.now()
    timestop = timestop.strftime('%Y%m%dT%H%M%SZ')
    timestart = timestart.strftime('%Y%m%dT%H%M%SZ')
    rw_data.write_params(p, os.path.join(p.outdatadir,
                                         'skim_simulator.output'))
    if ok is True:
        if p.progress_bar is True:
            __ = mod_tools.update_progress(1, 'All passes have been processed',
                                           '')
        else:
            __ = logger.info('All passes have been processed')
        logger.info("\n Simulated skim files have been written in "
                    "{}".format(p.outdatadir))
        logger.info(''.join(['-'] * 61))
        sys.exit(0)
    logger.error('\nERROR: At least one of the outputs was not saved.')
    sys.exit(1)


def make_skim_data(_proc_count, jobs, die_on_error, progress_bar):
    """ Compute SWOT-like data for all grids and all cycle, """
    # - Set up parallelisation parameters
    proc_count = min(len(jobs), _proc_count)

    status_updater = mod_tools.update_progress_multiproc
    jobs_manager = skimulator.mod_parallel.JobsManager(proc_count,
                                                       status_updater,
                                                       exc_formatter,
                                                       err_formatter)
    ok = jobs_manager.submit_jobs(worker_method_skim, jobs, die_on_error,
                                  progress_bar)

    if not ok:
        # Display errors once the processing is done
        jobs_manager.show_errors()

    return ok


def exc_formatter(exc):
    """Format exception returned by sys.exc_info() as a string so that it can
    be serialized by pickle and stored in the JobsManager."""
    error_msg = traceback.format_exception(exc[0], exc[1], exc[2])
    return error_msg


def err_formatter(pid, grid, cycle, exc):
    """Transform errors stored by the JobsManager into readable messages."""
    msg = None
    if cycle < 0:
        msg = '/!\ Error occurred while processing grid {}'.format(grid)
    else:
        _msg = '/!\ Error occurred while processing cycle {} on grid {}'
        msg = _msg.format(cycle, grid)
    return msg


def worker_method_skim(*args, **kwargs):
    msg_queue, sgridfile, p2, listsgridfile = args[:4]
    list_file, modelbox, model_data, modeltime, time_yaw, vac_yaw = args[4:]
    p = mod_tools.fromdict(p2)
    #   Load SKIM grid files (Swath and nadir)
    sgrid = mod.load_sgrid(sgridfile, p)
    # duplicate SKIM grids to assure that data are not modified and are
    # saved properly
    sgrid_tmp = mod.load_sgrid(sgridfile, p)

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
    if p.rain is True and p.rain_file is not None:
        rain_dic, rain_size = build_error.load_rain(p.rain_file)
    else:
        p.rain = False
    mss_path =  pkg_resources.resource_filename('skimulator',
                                                'share/noise_pdf_mss1d.npy')
    mssr_noise = numpy.load(mss_path)[()]
    #  Loop on all cycles
    for cycle in range(0, ncycle+1):
        # if ifile > (p.nstep*p.timestep + 1):
        #    break
        # Create SKIM-like data
        msg_queue.put((os.getpid(), sgridfile, cycle + 1, None))
        if p.file_input is None:
            model_data = []
        # Initialize all list of variables (each beam is appended to the
        # list)
        # Initialize velocity, indices and mask to empty lists
        output_var = {}
        for key in p.list_output:
            output_var[key] = []
        for key in p.list_input_var.keys():
            output_var[key] = []
        # Initialize noise to None if this noise is not computed
        # (so that the variable is not written in the netcdf)
        # err_var = {}
        # for key in p.list_err:
        #     err_var[key] = None
        # Initialize noises to empty lists if the noise is set to True
        if p.instr is True:
            output_var['instr'] = []
            output_var['dsigma'] = []
        if p.uwb is True:
            output_var['uwd'] = []
        if p.rain is True:
            output_var['rain'] = []
            output_var['gsig_atm_err'] = []
        if p.attitude is True:
            output_var['yaw'] = []
            output_var['yaw_aocs'] = []
            output_var['yaw_ted'] = []
        #if 'radial_angle' in p.list_output:
        #    output_var['radial_angle'] = sgrid.radial_angle
        # Loop over the beams
        time_all = []
        for i in range(len(p.list_pos) + 1):
            err_var_i = {}
            sgrid_tmp.lon = sgrid.lon[i]
            sgrid_tmp.lat = sgrid.lat[i]
            sgrid_tmp.time = sgrid.time[i]
            # If nadir, compute a different noise, not implemented so far,
            # Initialize at zero.
            if i == 0:
                output_var_i = {}
                # TODO remove the next line
                # for key in p.list_err:
                #     output_var_i[key ] = numpy.full(shape_0, numpy.nan)
                time = sgrid_tmp.time / 86400. + sgrid.cycle * cycle

              # mask_tmp = numpy.full(numpy.shape(sgrid_tmp.lon),
                #                      numpy.nan)
                #ssh_i, vindice = create_nadir_data()numpy create a mask from an a
                create = mod.create_SKIMlikedata(cycle, list_file,
                                                 modelbox, sgrid_tmp,
                                                 model_data, modeltime, p,
                                                 )
                output_var_i, time = create
                mod.compute_nadir_noise_skim(p, output_var_i, sgrid, cycle)

                if p.attitude is True:
                    yaw_aocs = build_error.make_yaw_aocs(time_yaw, vac_yaw, time)
                    # first_time = datetime.datetime.strptime(p.first_time,
                    #                                        '%Y-%m-%dT%H:%M:%SZ')
                    yaw_ted = 0 * yaw_aocs
                    yaw_total = 0 * yaw_aocs

            # Interpolate the velocity and compute the noise for each beam
            else:
                #Beam angle value to correct for attenuation in radial velocity
                beam_angle = p.list_angle[i - 1]
                # Read radial angle for projection on lon, lat reference
                radial_angle = sgrid.radial_angle[:, i - 1]


                ac_angle =  sgrid.angle[:, i - 1]
                ##############################
                # Compute SKIM like data data
                # TODO remove the next line
                # shape_all = (numpy.shape(listsgridfile)[0] * rcycle
                #              * (len(p.list_pos) + 1))
                create = mod.create_SKIMlikedata(cycle, list_file,
                                                 modelbox, sgrid_tmp,
                                                 model_data, modeltime, p,
                                                 )
                output_var_i, time = create
                build_error.compute_beam_noise_skim(p, output_var_i,
                                                    radial_angle, beam_angle,
                                                    ac_angle)
                if p.attitude is True:
                    yaw_aocs = build_error.make_yaw_aocs(time_yaw, vac_yaw, time)
                    first_time = datetime.datetime.strptime(p.first_time,
                                                            '%Y-%m-%dT%H:%M:%SZ')
                    yaw_ted = + build_error.make_yaw_ted(time, sgrid_tmp.cycle,
                                                         ac_angle,
                                                         first_time, beam_angle, p.instr_configuration)
                    # Conversion from microrad to m/s
                    yaw_total = ((yaw_aocs + yaw_ted) * const.vsat * 10**(-6)
                                 * numpy.cos(ac_angle))

            # Append variables for each beam
            time_all.append(time)
            if p.attitude is True:
                output_var['yaw'].append(yaw_total)
                output_var['yaw_aocs'].append(yaw_aocs)
                output_var['yaw_ted'].append(yaw_ted)
                #output_var['yaw_corr'].append(err_yaw)
            for key in output_var_i.keys():
                #if key in output_var_i.keys():
                output_var[key].append(output_var_i[key])
        # Compute correction with errdcos formulae

        if p.uwb is True:
            #ouput_var['uwb_corr']
            lon = numpy.transpose(numpy.array(sgrid.lon))[:, 1:]
            lat = numpy.transpose(numpy.array(sgrid.lat))[:, 1:]
            lon_nadir = numpy.array(sgrid.lon[0])
            lat_nadir = numpy.array(sgrid.lat[0])
            uwnd = numpy.transpose(numpy.array(output_var['uwnd']))
            vwnd = numpy.transpose(numpy.array(output_var['vwnd']))
            wnd_dir = numpy.mod(numpy.arctan2(vwnd, uwnd)[:, 1:], 2*numpy.pi)
            mss = (numpy.transpose(numpy.array(output_var['mssu']))
                   + numpy.transpose(numpy.array(output_var['mssc'])))
            mssx = numpy.transpose(numpy.array(output_var['mssx']))[:, 1:]
            mssy = numpy.transpose(numpy.array(output_var['mssy']))[:, 1:]
            mssxy = numpy.transpose(numpy.array(output_var['mssxy']))[:, 1:]
            hs = output_var['hs'][0]
            usr = numpy.transpose(numpy.array(output_var['ussr']))[:, 1:]
            ice = numpy.transpose(numpy.array(output_var['ice']))
            p.delta_azim = 20
            incl = sgrid.incl
            _angle = +  sgrid.angle
            if (sgrid.ipass %2) != 0:
                _angle = sgrid.angle + numpy.pi
            _combine = mod_uwb_corr.combine_usr
            usr_comb, usp_comb, mssr_comb = _combine(lon, lat, usr, mssx, mssy,
                                                     mssxy, mssr_noise,
                                                     p.delta_azim,
                                                     _angle, incl, wnd_dir)
            _closest = mod_uwb_corr.find_closest
            mssclose, hsclose = _closest(lon, lat, lon_nadir, lat_nadir, mss,
                                         mssr_comb, ice, hs, p.list_angle)
            #
            #hsclose = numpy.transpose(numpy.array(output_var['hs']))[:, 1:]
            #mssclose = + mss[:, 1:]
            # Temporary trick to compensate for bad usr correction
            #usr_comb = usr_comb / 3 + 2 * usr / 3
            uwd_est, usr_est2 = mod_uwb_corr.estimate_uwd(usr_comb, output_var, hsclose,
                                                mssclose, sgrid.radial_angle,
                                                p.list_angle)
            output_var['uwd_est'] = uwd_est
            output_var['usr_est2'] = usr_est2
            output_var['ussr_est'] = []
            output_var['mssr_est'] = []
            for i in range(len(output_var['ur_obs'])):
                if i == 0:
                    output_var['ussr_est'].append(usr_comb[:, 0])
                    output_var['mssr_est'].append(mssclose[:, 0])
                else:
                    output_var['ussr_est'].append(usr_comb[:, i-1])
                    output_var['mssr_est'].append(mssclose[:, i-1])
                output_var['uwd_est'][i] = output_var['uwd_est'][i] # / 3 + output_var['uwd'][i] * 2 /3

                corr = output_var['uwd'][i] - output_var['uwd_est'][i]
                output_var['ur_obs'][i][:]  = (output_var['ur_obs'][i][:]
                                               + corr)
        if p.rain is True:
            if p.rain_file is None:
                if 'rain' in output_var.keys():
                    for i in range(len(output_var['ur_obs'])):
                        _rain = output_var['rain'][i]
                        output_var['ur_obs'][i][numpy.where(_rain > p.rain_threshold)] = numpy.nan
            else:
                mean_time = numpy.mean(time)
                rain, rain_nad, gpia, gpia_nad = build_error.compute_rain(p, mean_time, sgrid,
                                                          rain_dic, rain_size)
                for i in range(len(output_var['ur_obs'])):
                    if i == 0:
                        _rain = rain_nad[:]
                        _gpia_err = gpia_nad[:]
                    else:
                        #Beam angle value to correct for attenuation in radial velocity
                        beam_angle = p.list_angle[i - 1]
                        ac_angle =  sgrid.angle[:, i - 1]
                        _rain = rain[:, i - 1]
                        _gpia_err = gpia[:, i - 1]
                    _gpia_err = mod_tools.convert_dbkm2ms(_gpia_err, ac_angle, beam_angle)
                    output_var['rain'].append(_rain)
                    output_var['gsig_atm_err'].append(_gpia_err)
                    output_var['ur_obs'][i][_rain > p.rain_threshold] = numpy.nan
                    output_var['ur_obs'][i][abs(_gpia_err) > 1] = numpy.nan
                    output_var['ur_obs'][i] += _gpia_err

        #   Compute errdcos if Formula is True

        # Compute directly bias if formula is False
        # if p.uss is True and p.formula is False:
        #    err_uss2 = compute_errussr(p, sgrid, mask, ur_uss, err_uss)
        #    # errdcos_tot.append(errdcos)
        #    # err_uss.append(err.err_uss)
        #    # , err_uss)
        # err_uss2 should be deleted, spectrum will provide stoke error
        #   Save outputs in a netcdf file
        if ((~numpy.isnan(numpy.array(output_var['vindice']))).any()
              or not p.file_input):
            sgrid.ncycle = cycle
            mod.save_SKIM(cycle, sgrid, time_all, output_var, p)
                      #    time=time, vindice=vindice_all,
                  #ur_model=ur_true_all, ur_obs=ur_obs, std_uss=std_uss,
                  #err_instr=err_instr, ur_uss=ur_uss, err_uss=err_uss2,
                  #u_model=u_true_all, v_model=v_true_all,
                  #errdcos=errdcos_tot)
        del time
        # if p.file_input: del index

    if p.file_input is not None:
        for key in model_data.vlon.keys():
            model_data.vlon[key] = (model_data.vlon[key] + 360) % 360

    modelbox[0] = (modelbox[0] + 360) % 360
    modelbox[1] = (modelbox[1] + 360) % 360
    del sgrid
    msg_queue.put((os.getpid(), sgridfile, None, None))

