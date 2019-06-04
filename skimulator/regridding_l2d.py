import numpy
import os
import sys
import glob
import ctypes
import datetime
import time
import math
import netCDF4
import pickle
import skimulator.const as const
import skimulator.rw_data as rw
import skimulator.mod_tools as mod_tools
import skimulator.mod_run as mod
import skimulator.mod_parallel
import logging
import traceback
import multiprocessing

logger = logging.getLogger(__name__)

# Using global variables is discouraged but this one will only be used as a
# simple way to share the observations vector with worker processes (the main
# process writes to this variable, the workers should only read from it)
obs_vector = {}


def run_l2d(p, die_on_error=False):
    config = p.config
    mod_tools.initialize_parameters(p)

    resols = p.resol_spatial_l2d # km
    resolt = p.resol_temporal_l2d  # days

    # Domain (spatial grid)
    if p.spatial_domain is not None:
        modelbox = numpy.array(p.spatial_domain, dtype='float')
        # Use convert to 360 data
        #modelbox[0] = (modelbox[0] + 360) % 360
        #if modelbox[1] != 360:
        #    modelbox[1] = (modelbox[1] + 360) % 360
    else:
        logger.error('Please provide modelbox_l2d for L2d reconstruction')
        sys.exit(1)
    lon0, lon1, lat0, lat1 = modelbox
    global_domain = False
    if (lon1 - lon0) > 200:
        global_domain = True
    # Use convert to 360 data
    lon0o = (lon0 + 360) % 360
    if lon1 != 360:
        lon1o = (lon1 + 360) % 360
            #modelbox[0] = (modelbox[0] + 360) % 360
                    #if modelbox[1] != 360:
                            #    modelbox[1] = (modelbox[1] + 360) % 360

    dlon, dlat = p.posting_l2d
    # Time domain
    time0, time1, dt = p.time_domain

    # #####################################
    # Spatial grid construction
    # #####################################
    grd = make_grid(lon0, lon1, dlon, lat0, lat1, dlat)
    dlatlr = 1 # dlat
    dlonlr = 1 # dlon
    grd['dlon'] = dlon
    grd['dlonlr'] = dlonlr
    grd['dlat'] = dlat
    grd['dlatlr'] = dlatlr
    grd['nxlr'] = int((lon1-lon0)/dlonlr)
    grd['nylr'] = int((lat1-lat0)/dlatlr)
    # Filtering as a function of latitude (10d at Equator, 14 d at 60º)
    resolt = (-4 * numpy.cos(numpy.deg2rad(grd['lat'])) + 9) * resolt

    # #####################################
    # READ L2B OBS
    # #####################################
    tcycle = p.satcycle
    # TODO change
    resoltmax = numpy.max(resolt)
    c0 = int((max(time0 - resoltmax, 0)) / tcycle) + 1
    c1 = int((time1 + resoltmax) / tcycle) + 1
    obs = {}
    # c0 = 1 ; c1 =  2
    p2 = mod_tools.todict(p)
    jobs = []
    print(datetime.datetime.now())
    for cycle in range(c0, c1 + 1, 1):
        _pat = '{}_c{:02d}_p*'.format(p.config, cycle)
        pat = os.path.join(p.outdatadir, _pat)
        listfile = glob.glob(pat)
        # Create specific funtion to retrieve time_unit
        filev = os.path.join(p.outdatadir, listfile[0])
        _, _, time_unit = read_l2b(filev)
        for pattern in listfile:
            jobs.append([p2, pattern, lon0, lon1, lat0, lat1, time0, time1,
                         resoltmax, dlonlr, dlatlr])
    ok = False
    task = 0
    results = []
    try:
        ok, results = build_obs(p.proc_count, jobs, die_on_error,
                                p.progress_bar)
    except skimulator.mod_parallel.MultiprocessingError:
        logger.error('An error occurred with the multiprocessing framework')
        traceback.print_exception(*sys.exc_info())
        sys.exit(1)
    except skimulator.mod_parallel.DyingOnError:
        logger.error('An error occurred and all errors are fatal')
        sys.exit(1)

    print()
    print('Aggregating results...', datetime.datetime.now())
    for i in range(len(results)):
        obs_rjobs = results[i]
        #key, ind_key, obs_r = results[i]
        for j in range(len(obs_rjobs)):
            try:
                obs_r = obs_rjobs[j]
            except:
                obs_r = obs_rjobs
            obs_r_keys = list(obs_r.keys())
            for key in obs_r_keys:
                if key not in obs.keys():
                    obs[key] = {}
                for ind_key in obs_r[key].keys():
                    if ind_key not in obs[key].keys():
                        obs[key][ind_key] = []
                    obs[key][ind_key] = numpy.concatenate((obs[key][ind_key],
                                                          obs_r[key][ind_key]))
                del(obs_r[key])
    del(results)

    print('Building shared obs vector...', datetime.datetime.now())
    # Populate the shared observations vector
    global obs_vector
    obs_vector = {}
    for key in obs:
        obs_vector[key] = {}
        obs_ind_keys = list(obs[key].keys())
        for ind_key in obs_ind_keys:
            _buffer = multiprocessing.RawArray(ctypes.c_float,
                                               numpy.size(obs[key][ind_key]))
            narray = numpy.frombuffer(_buffer, dtype='float32')
            numpy.copyto(narray, obs[key][ind_key])
            obs_vector[key][ind_key] = _buffer

            # Release memory from original array to keep the total amount of
            # RAM as low as possible
            if 'time' != key:
                del(obs[key][ind_key])  # release memory as soon as possible

    print(datetime.datetime.now())

    print('Memory housekeeping...', datetime.datetime.now())
    keys = list(obs.keys())
    for key in keys:
        if 'time' != key:
            del(obs[key])
    print(datetime.datetime.now())

    #import pdb ; pdb.set_trace()
    # #####################################
    # Independant time loop for each analysis
    # #####################################
    enstime = numpy.arange(time0, time1 + dt, dt)
    time_model = numpy.arange(time0, time1, p.timestep)
    # Center at noon
    window = dt / 2.
    enstime = enstime + window
    indice_mask = make_mask(p, 'ucur', grd)
    list_key = {'ucur':'ux_true', 'vcur':'uy_true'} #, 'uwnd':'uwnd',
                #'vwnd': 'vwnd'}
    list_input = {}
    for key in list_key:
        list_input[key] = p.list_input_var_l2d[key]

    for it, timeref in enumerate(enstime):
        print(it, timeref)
        iobs = {}
        for ind_key in obs['time'].keys():
            iobs[ind_key] = numpy.where((obs['time'][ind_key] > (timeref
                                         - numpy.max(resoltmax)))
                                         & (obs['time'][ind_key] < (timeref
                                         + numpy.max(resoltmax))))[0]
        make_oi(p, p2, grd, iobs, resols, timeref, resolt, indice_mask,
                die_on_error=die_on_error)
        # Interpolate model data to have a true reference
        #indice_model = it
        #print('read model')
        #model_data, out_var, list_file = read_model(p, it, p.dim_time,
        #                                            list_input)
        #print('interpolate model')
        #if global_domain is False:
        #    interpolate_model(p, model_data, out_var, grd, list_key)

        pattern = os.path.join(p.outdatadir, '{}_{}_l2d_'.format(config,
                                                                 p.config_l2d))
        save_l2d(pattern, timeref, window, time_unit, grd)
        print(datetime.datetime.now())


def offline_interpolation(p):
    pattern = os.path.join(p.outdatadir, '{}_{}_l2d_'.format(p.config,
                                                             p.config_l2d))
    listfile = glob.glob(os.path.join(p.outdatadir, '{}*.nc'.format(pattern)))
    str_format = '%Y-%m-%dT%H:%M:%SZ'
    list_key = {'ucur':'ux_true', 'vcur':'uy_true'}
    list_input = {}
    dimtime = 'time'
    dimlon = 'lon'
    dimlat = 'lat'
    longname = { "ux_noerr": "Error-free zonal velocity",
                 "uy_noerr": "Error-free meridional velocity",
                "ux_obs": "Observed zonal velocity",
                "uy_obs": "Observed meridional velocity",
                "ux_model": "Error-free zonal velocity",
                "uy_model": "Error-free meridional velocity",
                "ux_true": "True zonal velocity",
                "uy_true": "True meridional velocity",
                }
    unit = {"ux_noerr": "m/s", "ux_obs": "m/s",
            "uy_noerr": "m/s", "uy_obs": "m/s",
            "ux_true": "m/s", "uy_true": "m/s"
            }
    for key in list_key:
        list_input[key] = p.list_input_var_l2d[key]

    for ifile in listfile:
        fid = netCDF4.Dataset(ifile, 'a')
        start_date = datetime.datetime.strptime(fid.time_coverage_start,
                                                str_format)
        end_date = datetime.datetime.strptime(fid.time_coverage_end,
                                              str_format)
        model_start =  datetime.datetime.strptime(p.first_time, str_format)
        it = -(model_start - start_date).days
        print('read l2d coordinates')
        grd = {}
        grd['lon'] = fid.variables['lon'][:]
        grd['lat'] = fid.variables['lat'][:]
        grd['lon2'], grd['lat2'] = numpy.meshgrid(grd['lon'], grd['lat'])
        print('read model')
        print(start_date, model_start)
        model_data, out_var, list_file = read_model(p, it, p.dim_time,
                                                    list_input)
        print('interpolate model')
        interpolate_model(p, model_data, out_var, grd, list_key)
        list_var = ('ux_true', 'uy_true')
        for key in list_var:
            nvar = '{}'.format(key)
            var = fid.createVariable(nvar, 'f4', (dimtime, dimlat, dimlon),
                                     fill_value=-1.36e9)
            value = grd[key]
            try:
                var.units = unit[str(key)]
            except:
                var.units = ''
            try:
                var.long_name = longname[str(key)]
            except:
                var.long_name = str(key)
            if value.any():
                mask = numpy.isnan(value)
                value[numpy.where(mask)] = -1.36e9
                mask_ind = numpy.where(value < -1e7)
                value[mask_ind] = -1.36e9
                mask_ind = numpy.where(value > 1e7)
                value[mask_ind] = -1.36e9
                mask_ind = numpy.where(value == numpy.PINF)
                value[mask_ind] = -1.36e9
                var[0, :, :] = value




def worker_build_obs(*args, **kwargs):
    msg_queue, p2, pattern = args[:3]
    p = mod_tools.fromdict(p2)
    lon0, lon1, lat0, lat1, time0, time1, resolt, dlonlr, dlatlr = args[3:12]
    res_queue = args[-1]

    filev = os.path.join(p.outdatadir, pattern)
    obs = {}
    if os.path.isfile(filev):
        obs_i, mask_data, time_unit = read_l2b(filev)
        _lon = numpy.mod(obs_i['lon'] + 180, 360) - 180
#                if lon1o < lon0o:
#                    mask_lon = ((obs_i['lon'] < lon0o) & (obs_i['lon'] > lon1o))
#                else:
#                    mask_lon = ((obs_i['lon'] > lon0o) & (obs_i['lon'] < lon1o))
        if lon1 < lon0:
            mask_lon = ((_lon < lon0) & (_lon > lon1))
        else:
            mask_lon = ((_lon > lon0) & (_lon < lon1))
        mask_lat = ((obs_i['lat'] > lat0) & (obs_i['lat'] < lat1))
        mask_time = ((obs_i['time'] > time0 - resolt)
                     & (obs_i['time'] < time1 + resolt))
        mask = (mask_time & mask_lat & mask_lon & mask_data)
        #if not numpy.any(mask):
        #    msg_queue.put((os.getpid(), filev, None, None))
        #    return None
        for key in obs_i.keys():
            obs_i[key] = obs_i[key][mask]
        _lon = _lon[mask]
        ind_i = numpy.floor((obs_i['lat'] - lat0) / dlatlr)
        ind_j = numpy.floor(numpy.mod(_lon - lon0, 360) / dlonlr)
        unique_i = list(set(ind_i))
        unique_j = list(set(ind_j))
        for key in obs_i.keys():
             if key not in obs.keys():
                obs[key] = {}
        for i in unique_i:
            for j in unique_j:
                mask_ind = numpy.where((ind_i==i) & (ind_j == j))
                if mask_ind[0].any():
                    ind_key = 10000 * int(i) + int(j)  # '{}_{}'.format(int(i), int(j))
                    for key in obs_i.keys():
                        if ind_key not in obs[key].keys():
                            obs[key][ind_key] = []
                        obs[key][ind_key] = numpy.concatenate((obs[key][ind_key],
                                           obs_i[key][mask_ind]))

    res_queue.put(obs)  # Pass result to parent method
    msg_queue.put((os.getpid(), filev, None, None))


def build_obs(_proc_count, jobs, die_on_error, progress_bar):
    """ Compute SWOT-like data for all grids and all cycle, """
    # - Set up parallelisation parameters
    proc_count = min(len(jobs), _proc_count)

    status_updater = mod_tools.update_progress_multiproc
    jobs_manager = skimulator.mod_parallel.JobsManager(proc_count,
                                                       status_updater,
                                                       exc_formatter,
                                                       err_formatter)

    # Include the results queue as last argument for each job
    for j in jobs:
        j.append(jobs_manager.res_queue)

    res = []
    ok = jobs_manager.submit_jobs(worker_build_obs, jobs, die_on_error,
                                  progress_bar, results=res.append)

    if not ok:
        # Display errors once the processing is done
        jobs_manager.show_errors()

    return ok, res



def save_l2d(filenc, timeref, window, time_unit, grd):
    dateformat = '%Y-%m-%dT%H:%M:%SZ'
    datewformat = '%Y%m%dT%H%M%S'
    time_std = netCDF4.num2date(timeref - window, time_unit)
    time_start = time_std.strftime(format=dateformat)
    time_std = netCDF4.num2date(timeref + window, time_unit)
    time_end = time_std.strftime(format=dateformat)
    time_std = netCDF4.num2date(timeref, time_unit)
    time_middle = time_std.strftime(format=datewformat)
    pattern = '{}{}.nc'.format(filenc, time_middle)
    grd['time'] = timeref
    metadata = {}
    metadata['first_time'] = time_unit
    metadata['file'] = pattern
    metadata['time_coverage_start'] = time_start
    metadata['time_coverage_end'] = time_end
    metadata['file'] = pattern
    if os.path.exists(pattern):
        os.remove(pattern)
    if 'uy_true' in grd.keys():
        rw.write_l2d(metadata, grd, ux_noerr=grd['ux_noerr'],
                     uy_noerr=grd['uy_noerr'], ux_obs=grd['ux_obs'],
                     uy_obs=grd['uy_obs'], ux_true=grd['ux_true'],
                     uy_true=grd['uy_true'])
    else:
        rw.write_l2d(metadata, grd, ux_noerr=grd['ux_noerr'],
                     uy_noerr=grd['uy_noerr'], ux_obs=grd['ux_obs'],
                     uy_obs=grd['uy_obs']) #, ux_true=grd['ux_true'],
                     #uy_true=grd['uy_true'])


def read_l2b(nfile, model_nan=0):
    fid = netCDF4.Dataset(nfile, 'r')
    obs_i = {}
    obs_i['ux'] = numpy.ma.array(fid.variables['ucur'][:]).flatten()
    obs_i['uy'] = numpy.ma.array(fid.variables['vcur'][:]).flatten()
    if 'ur_true' in fid.variables.keys():
        obs_i['ur_true'] = numpy.ma.array(fid.variables['ur_true'][:]).flatten()
    if 'ur_obs' in fid.variables.keys():
        obs_i['ur_obs'] = numpy.ma.array(fid.variables['ur_obs'][:]).flatten()
    obs_i['lon'] = numpy.ma.array(fid.variables['lon'][:]).flatten()
    #obs_i['lon'] = numpy.mod(obs_i['lon'] -180, 360) + 180
    obs_i['lon'] = numpy.mod(obs_i['lon'] + 360, 360)
    obs_i['lat'] = numpy.ma.array(fid.variables['lat'][:]).flatten()
    obs_i['time'] = numpy.ma.array(fid.variables['time'][:]).flatten()
    angle = numpy.ma.array(fid.variables['radial_angle'][:]).flatten()
    obs_i['angle'] = numpy.mod(angle, 2 * numpy.pi)
    mask_invalid = (numpy.ma.getmaskarray(obs_i['ux'])
                    | numpy.ma.getmaskarray(obs_i['uy'])
                    | numpy.ma.getmaskarray(obs_i['ur_true'])
                    | numpy.ma.getmaskarray(obs_i['ur_obs']))
    mask_invalid = (mask_invalid | (obs_i['ux']==model_nan)
                    | (obs_i['uy']==model_nan)
                    | (obs_i['ur_true']==model_nan)
                    | (obs_i['ur_obs']==model_nan))
    mask_data = ~(mask_invalid)
    time_unit = fid.variables['time'].units
    fid.close()
    return obs_i, mask_data, time_unit


def make_grid(lon0, lon1, dlon, lat0, lat1, dlat):
    grd = {}
    grd['lon'] = numpy.arange(lon0, lon1 + dlon, dlon)
    grd['lon'] = numpy.mod(grd['lon'] + 360, 360) - 360
    grd['lat'] = numpy.arange(lat0, lat1 + dlat, dlat)
    grd['lon'] = numpy.mod(grd['lon'] +360, 360)
    grd['lon2'], grd['lat2'] = numpy.meshgrid(grd['lon'], grd['lat'])
    grd['nx'] = len(grd['lon'])
    grd['ny'] = len(grd['lat'])
    return grd


def make_oi(p, p2, grd, iobs, resols, timeref, resolt, index,
            die_on_error=False):
    indices = [(j, i) for j, i in zip(index[0], index[1])]

    # Split the list of (j, i) couples into chunks that the workers will
    # process in parallel. The method to decide how many chunks are required is
    # purely empirical
    chunk_size = int(len(indices) / (3 * p.proc_count))  # 3 is a magic number
    chunk_size = max(chunk_size, 1)  # avoid chunk_size = 0 due to rounding
    ind_count = math.ceil(len(indices) / chunk_size)
    jobs = []
    for k in range(0, ind_count):
        jobs.append([p2, grd, iobs, resols, timeref, resolt, index,
                     indices[k * chunk_size :(k + 1) * chunk_size]])

    ok = False
    try:
        ok = par_make_oi(grd, p.proc_count, jobs, die_on_error, p.progress_bar)
    except skimulator.mod_parallel.MultiprocessingError:
        logger.error('An error occurred with the multiprocessing framework')
        traceback.print_exception(*sys.exc_info())
        sys.exit(1)
    except skimulator.mod_parallel.DyingOnError:
        logger.error('An error occurred and all errors are fatal')
        sys.exit(1)


def init_oi_worker(obs_dict):
    """Initialize the obs_vector global variable with the value passed by the
    parent process.
    obs_vector is not modified by the worker, so sharing this object between
    the main process and the workers should not consume additional memory due
    to Copy On Write fork/spawn behavior (at least on Linux)."""
    global obs_vector
    obs_vector = obs_dict


def par_make_oi(grd, _proc_count, jobs, die_on_error, progress_bar):
    """ Compute SWOT-like data for all grids and all cycle, """
    global obs_vector

    # - Set up parallelisation parameters
    proc_count = min(len(jobs), _proc_count)

    status_updater = mod_tools.update_progress_multiproc
    jobs_manager = skimulator.mod_parallel.JobsManager(proc_count,
                                                       status_updater,
                                                       exc_formatter,
                                                       err_formatter,
                                                       init_oi_worker,
                                                       (obs_vector,))

    for j in jobs:
        j.append(jobs_manager.res_queue)

    jobs_res = []
    ok = jobs_manager.submit_jobs(worker_make_oi, jobs, die_on_error,
                                  progress_bar, results=jobs_res.append)

    print()
    print('par_make_oi complete', datetime.datetime.now())
    if not ok:
        # Display errors once the processing is done
        jobs_manager.show_errors()

    # Reconstruct results grids
    grd['ux_noerr'] = numpy.full((grd['ny'], grd['nx']), numpy.nan)
    grd['uy_noerr'] = numpy.full((grd['ny'], grd['nx']), numpy.nan)
    grd['ux_obs'] = numpy.full((grd['ny'], grd['nx']), numpy.nan)
    grd['uy_obs'] = numpy.full((grd['ny'], grd['nx']), numpy.nan)
    for k in range(0, len(jobs_res)):
        for res in jobs_res[k]:
            j, i, ux_noerr, uy_noerr, ux_obs, uy_obs = res
            grd['ux_noerr'][j, i] = ux_noerr
            grd['uy_noerr'][j, i] = uy_noerr
            grd['ux_obs'][j, i] = ux_obs
            grd['uy_obs'][j, i] = uy_obs

    del(jobs_res)

    return ok


def worker_make_oi(*args, **kwargs):
    msg_queue = args[0]
    p2, grd, iobs, resols, timeref, resolt, index, indices = args[1:9]
    res_queue = args[-1]
    p = mod_tools.fromdict(p2)

    global obs_vector
    jobs_res = []

    for j, i in indices:
        obs = {}
        ni = 1
        nj = 1
        dlat = grd['dlat']
        dlatlr = grd['dlatlr']
        dlon = grd['dlon']
        dlonlr = grd['dlonlr']
        resoltij = resolt[j]
        # Filtering as a function of latitude (90km at Eqautor, 40 km at 60º)
        resolsij = (100 * numpy.cos(numpy.deg2rad(grd['lat'][j])) - 10) * resols
        jlr = int(numpy.floor(j*dlat/dlatlr))
        ilr = int(numpy.floor(i*dlon/dlonlr))
        for jx in range(max(0, jlr-nj), min(jlr + nj + 1, grd['nylr'])):
            for jy in range(max(0, ilr-ni), min(ilr + ni + 1, grd['nxlr'])):
                # ind_key = '{}_{}'.format(jx, jy)
                ind_key = 10000 * int(jx) + int(jy)  # '{}_{}'.format(int(i), int(j))
                if ind_key not in iobs.keys():
                    continue
                if ~iobs[ind_key].any():
                    continue
                for key in obs_vector.keys():
                    if key not in obs.keys():
                        obs[key] = []

                    # Rebuild numpy array from shared memory (no copy)
                    _buffer = obs_vector[key][ind_key]
                    narray = numpy.frombuffer(_buffer, dtype='float32')

                    _obs = narray[numpy.array(iobs[ind_key], dtype='int')]
                    obs[key] = numpy.concatenate(([obs[key], _obs]))

        # TODO: to be optimized, especially for global, ...
        flat_ind = j + i*grd['nylr']
        if 'lon' not in obs.keys():
            continue

        # Handle IDL and Greenwich line
        _lonobs = numpy.mod(obs['lon'] + 360, 360)
        _longrd = numpy.mod(grd['lon2'][j, i] + 360, 360)
        dlon = 180 - abs(abs(_lonobs - _longrd) - 180)

        # dist = 110. * (numpy.cos(numpy.deg2rad(grd['lat'][j]))**2 * (dlon)**2
        #               + (obs['lat'] - grd['lat2'][j, i])**2)**0.5
        fac = 1
        if abs(grd['lat'][j]) < 10:
             fac = fac + (10 - abs(grd['lat'][j])) / 8
        dist = 110. * (1 / fac*numpy.cos(numpy.deg2rad(grd['lat'][j]))**2 * (dlon)**2
                       + fac * (obs['lat'] - grd['lat2'][j, i])**2)**0.5
        iiobs=numpy.where((dist < resolsij))[0]
        if len(iiobs)>=2:
            H = numpy.zeros((len(iiobs), 2))
            H[:, 0] = numpy.cos(obs['angle'][iiobs])
            H[:, 1] = numpy.sin(obs['angle'][iiobs])
            win_s = numpy.exp(-dist[iiobs]**2/(0.5*resolsij)**2)
            time_cen = obs['time'][iiobs] - timeref
            win_t = numpy.exp(-time_cen**2/(0.5 * resoltij)**2) # exp window
            Ri = win_s * win_t
            RiH = numpy.tile(Ri, (2, 1)).T*H
            M = numpy.dot(H.T, RiH)
            if numpy.linalg.cond(M) < 1e3:
                Mi = numpy.linalg.inv(M)
                eta_true = numpy.dot(numpy.dot(Mi, RiH.T),
                                     obs['ur_true'][iiobs])
                eta_obs = numpy.dot(numpy.dot(Mi, RiH.T),
                                    obs['ur_obs'][iiobs])

                jobs_res.append((j, i, eta_true[0], eta_true[1], eta_obs[0],
                                 eta_obs[1]))

    # Pass results for all the processed (j, i) to the parent process
    res_queue.put(jobs_res)

    # Notify parent process about job completion
    msg_queue.put((os.getpid(), j, i, None))


def read_model(p, ifile, dim_time, list_input):
    model_data, list_file = mod.load_coordinate_model(p)
    model_data.read_coordinates(p)
    nfile = int(ifile)
    print(nfile)
    filename = os.path.join(p.indatadir, list_file[nfile])
    print(filename)
    model_step_ctor = getattr(rw, p.model)
    out_var = {}
    for i in range(dim_time):
        model_step = model_step_ctor(p, ifile=(filename, ),
                                     list_input_var=list_input,
                                     time=i)
        model_step.read_var(p)
        for key in list_input.keys():
            if not key in out_var.keys():
                out_var[key] = model_step.input_var[key] / dim_time
            else:
                out_var[key] = (out_var[key]
                                + model_step.input_var[key] / dim_time)
    return model_data, out_var, list_file


def make_mask(p, key, grid):
    """ Return indices of points on the ocean (non masked value in the model)
    """
    list_input = {key: p.list_input_var_l2d[key]}
    model_data, out_var, list_file = read_model(p, 0, 1, list_input)
    mask_ucur = numpy.ma.getmaskarray(out_var['ucur'])
    mask_ucur = (mask_ucur | numpy.isnan(out_var['ucur']))
    out_var['ucur'] = numpy.ma.array(out_var['ucur'], mask=mask_ucur)
    list_key = {'ucur':'mask'}
    interpolate_model(p, model_data, out_var, grid, list_key)
    mask_index = numpy.where(~numpy.ma.getmaskarray(grid['mask']))
                            # & (grid['mask'] != 0)
                            # & (~numpy.isnan(grid['mask'])))
    # TODO TO proof if mask_index is empty
    return mask_index


def interpolate_model(p, model_data, model_var, grd, list_key):
    import pyresample as pr
    wrap_lon = pr.utils.wrap_longitudes
    geom = pr.geometry.SwathDefinition
    interp = mod.interpolate_irregular_pyresample
    lon = wrap_lon(grd['lon2'])

    for ikey, okey in list_key.items():
        if len(p.list_input_var[ikey]) > 2:
            grid_number = p.list_input_var[ikey][2]
        else:
            grid_number = 0
        _lon = wrap_lon(model_data.vlon[grid_number])
        _lat = model_data.vlat[grid_number]
        if len(numpy.shape(_lon)) <= 1:
            _lon, _lat = numpy.meshgrid(_lon, _lat)
        swath_def = geom(lons=_lon, lats=_lat)
        grid_def = geom(lons=lon, lats=grd['lat2'])

        var = model_var[ikey]
        grd[okey] = interp(swath_def, var, grid_def, p.resol,
                           interp_type=p.interpolation)
        grd[okey][grd[okey] == p.model_nan] = numpy.nan
    return grd


def exc_formatter(exc):
    """Format exception returned by sys.exc_info() as a string
so that it can
    be serialized by pickle and stored in the JobsManager."""
    error_msg = traceback.format_exception(exc[0], exc[1], exc[2])
    return error_msg


def err_formatter(pid, it, cycle, exc):
    """Transform errors stored by the JobsManager into readable messages."""
    msg = '/!\ Error occurred while processing it {}'.format(it)
    #if cycle < 0:
    #    msg = '/!\ Error occurred while processing grid'
    #else:
    #    _msg = '/!\ Error occurred while processing cycle {}'
    #    msg = _msg.format(cycle)
    return msg
