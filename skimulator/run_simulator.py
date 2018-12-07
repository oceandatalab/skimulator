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
import skimulator.build_swath as build_swath
import skimulator.rw_data as rw_data
import skimulator.build_error as build_error
import skimulator.mod_tools as mod_tools
import skimulator.const as const
import skimulator.mod_run as mod
import skimulator.mod_uwb_corr as mod_uwb_corr
import multiprocessing
# Define logger level for debug purposes
logger = logging.getLogger(__name__)
#logger = multiprocessing.log_to_stderr()
#logger.setLevel(logging.DEBUG)

# - Define global variables for progress bars
istep = 0
ntot = 1
ifile = 0


def run_simulator(p):
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
        # if (p.uss is True and p.footprint_std is not None
        #     and p.footprint_std != 0):
        #    if len(numpy.shape(model_data.vlonu)) == 1:
        #        model_data.lon2D, model_data.lat2D = numpy.meshgrid(
        #                                                model_data.vlonu,
        #                                                model_data.vlatu)
        #    else:
        #        model_data.lon2D = model_data.vlonu
        #        model_data.lat2D = model_data.vlatu
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
        build_swath.orbit2swath(modelbox, p, orb)
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

    # - Loop on SKIM grid files
    jobs = []
    p2 = mod_tools.todict(p)
    for sgridfile in listsgridfile:
        jobs.append([sgridfile, p2, listsgridfile, list_file,
                     modelbox, model_data, modeltime])
    ok = make_skim_data(p.proc_count, jobs)

    # - Write Selected parameters in a txt file
    timestop = datetime.datetime.now()
    timestop = timestop.strftime('%Y%m%dT%H%M%SZ')
    timestart = timestart.strftime('%Y%m%dT%H%M%SZ')
    rw_data.write_params(p, os.path.join(p.outdatadir,
                                         'skim_simulator.output'))
    if ok is True:
        __ = mod_tools.update_progress(1, 'All passes have been processed', '')
        logger.info("\n Simulated skim files have been written in "
                    "{}".format(p.outdatadir))
        logger.info(''.join(['-'] * 61))
        sys.exit(0)
    logger.error('\nERROR: At least one of the outputs was not saved.')
    sys.exit(1)


def make_skim_data(_proc_count, jobs):
    """ Compute SWOT-like data for all grids and all cycle, """
    # - Set up parallelisation parameters
    proc_count = min(len(jobs), _proc_count)

    manager = multiprocessing.Manager()
    msg_queue = manager.Queue()
    pool = multiprocessing.Pool(proc_count)
    # Add the message queue to the list of arguments for each job
    # (it will be removed later)
    [j.append(msg_queue) for j in jobs]
    chunk_size = int(math.ceil(len(jobs) / proc_count))
    status = {}
    for n, w in enumerate(pool._pool):
        status[w.pid] = {'done': 0, 'total': 0, 'grids': None, 'extra': ''}
        proc_jobs = jobs[n::proc_count]
        status[w.pid]['grids'] = [j[0] for j in proc_jobs]
        status[w.pid]['total'] = len(proc_jobs)
    sys.stdout.write('\n' * proc_count)
    tasks = pool.map_async(worker_method_skim, jobs, chunksize=chunk_size)
    sys.stdout.flush()
    ok = True
    while not tasks.ready():
        if not msg_queue.empty():
            msg = msg_queue.get()
            _ok = mod_tools.update_progress_multiproc(status, msg)
            ok = ok and _ok

        time.sleep(0.1)

    while not msg_queue.empty():
        msg = msg_queue.get()
        mod_tools.update_progress_multiproc(status, msg)

    sys.stdout.flush()
    pool.close()
    pool.join()
    return ok


def worker_method_skim(*args, **kwargs):
    _args = list(args)[0]
    msg_queue = _args.pop()
    sgridfile = _args[0]
    p2, listsgridfile, list_file, modelbox, model_data, modeltime = _args[1:]
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


    #  Loop on all cycles
    for cycle in range(0, ncycle+1):
        # if ifile > (p.nstep*p.timestep + 1):
        #    break
        # Create SKIM-like data
        msg_queue.put((os.getpid(), sgridfile, cycle + 1))
        if p.file_input is None:
            model_data = []
        # Initialize all list of variables (each beam is appended to the
        # list)
        # Initialize velocity, indices and mask to empty lists
        output_var = {}
        for key in p.list_output:
            output_var[key] = []
        # Initialize noise to None if this noise is not computed
        # (so that the variable is not written in the netcdf)
        # err_var = {}
        # for key in p.list_err:
        #     err_var[key] = None
        # Initialize noises to empty lists if the noise is set to True
        if p.instr is True:
            output_var['instr'] = []
        if p.uwb is True:
            output_var['uwb'] = []
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
                shape_0 = numpy.shape(sgrid_tmp.lon)
                for key in p.list_output:
                    if ('ssh' not in key) or ('indice' not in key):
                        output_var_i[key ] = numpy.full(shape_0, numpy.nan)
                # for key in p.list_err:
                #     output_var_i[key ] = numpy.full(shape_0, numpy.nan)
                time = sgrid_tmp.time / 86400. + sgrid.cycle * cycle

                # mask_tmp = numpy.full(numpy.shape(sgrid_tmp.lon),
                #                      numpy.nan)
                #ssh_i, vindice = create_nadir_data()numpy create a mask from an a
            # Interpolate the velocity and compute the noise for each beam
            else:
                #Beam angle value to correct for attenuation in radial velocity
                beam_angle = p.list_angle[i - 1]
                # Read radial angle for projection on lon, lat reference
                radial_angle = sgrid.radial_angle[:, i - 1]


                ac_angle =  sgrid.angle[:, i - 1]
                ##############################
                # Compute SKIM like data data
                try:
                    shape_all = (numpy.shape(listsgridfile)[0] * rcycle
                                 * (len(p.list_pos) + 1))
                    create = mod.create_SKIMlikedata(cycle, shape_all, list_file,
                                                 modelbox, sgrid_tmp,
                                                 model_data, modeltime,
                                                 radial_angle, ac_angle,
                                                 beam_angle, p,
                                                 progress_bar=True)
                    output_var_i, time = create
                except:
                    import sys
                    e = sys.exc_info()
                    logger.error('bouh', exc_info=e)
                    msg_queue.put((os.getpid(), sgridfile, -1))
            # Append variables for each beam
            time_all.append(time)
            for key in p.list_output:
                if key in output_var_i.keys():
                    output_var[key].append(output_var_i[key])
        # Compute correction with errdcos formulae

        if p.uwb is True:
            #ouput_var['uwb_corr']
            _tmp = mod_uwb_corr.compute_erruwb(p, sgrid, output_var)
            output_var['uwb_corr'] = _tmp
            for i in range(len(output_var['ur_obs'])):
                output_var['ur_obs'][i][:]  = (output_var['ur_obs'][i][:]
                                               + output_var['uwb_corr'][i][:])

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
            try:
                mod.save_SKIM(cycle, sgrid, time_all, output_var, p)
                          #    time=time, vindice=vindice_all,
                      #ur_model=ur_true_all, ur_obs=ur_obs, std_uss=std_uss,
                      #err_instr=err_instr, ur_uss=ur_uss, err_uss=err_uss2,
                      #u_model=u_true_all, v_model=v_true_all,
                      #errdcos=errdcos_tot)
            except:
                msg_queue.put((os.getpid(), sgridfile, -1))
        del time
        # if p.file_input: del index

    if p.file_input is not None:
        model_data.vlonu = (model_data.vlonu + 360) % 360
        model_data.vlonv = (model_data.vlonv + 360) % 360

    modelbox[0] = (modelbox[0] + 360) % 360
    modelbox[1] = (modelbox[1] + 360) % 360
    del sgrid
    msg_queue.put((os.getpid(), sgridfile, None))

