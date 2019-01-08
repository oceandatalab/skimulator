import numpy
import os
import sys
import glob
import datetime
import time
import netCDF4
import skimulator.const as const
import skimulator.rw_data as rw
import skimulator.mod_tools as mod_tools
import skimulator.mod_run as mod
import skimulator.mod_parallel
import logging
import traceback
logger = logging.getLogger(__name__)


def run_l2d(p, die_on_error=False):

    config = p.config
    mod_tools.initialize_parameters(p)

    L2Berror = False
    resols = p.resol_spatial_l2d # km
    resolt = p.resol_temporal_l2d  # days

    # Domain (spatial grid)
    if p.modelbox_l2d is not None:
        modelbox = numpy.array(p.modelbox_l2d, dtype='float')
        # Use convert to 360 data
        modelbox[0] = (modelbox[0] + 360) % 360
        if modelbox[1] != 360:
            modelbox[1] = (modelbox[1] + 360) % 360
    else:
        logger.error('Please provide modelbox_l2d for L2d reconstruction')
        sys.exit(1)
    lon0, lon1, lat0, lat1 = modelbox
    dlon, dlat = p.posting_l2d
    # Time domain
    time0, time1, dt = p.time_domain

    # #####################################
    # READ L2B OBS
    # #####################################
    tcycle=29.
    norb=412
    c0 = int((max(time0 - resolt, 0)) / tcycle) + 1
    p0 = int(numpy.mod((max(time0 - resolt,0)), tcycle) * 2*norb/tcycle) + 1
    c1 = int((time1+resolt)/tcycle) + 1
    p1 = int(numpy.mod((time1+resolt),tcycle) * 2*norb/tcycle) + 1
    obs = {}
    for cycle in range(c0, c1 + 1, 1):
        p0a = 1
        p1a = 2*norb
        if cycle == c0:
            p0a=p0
        if cycle == c1:
            p1a=p1
        for pa in range(p0a, p1a + 1, 1):
            pattern = '{}_c{:02d}_p{:03d}.nc'.format(p.config, cycle, pa)
            filev = os.path.join(p.outdatadir, pattern)
            if os.path.isfile(filev):
                obs_i, mask_data, time_unit = read_l2b(filev)
                if lon1 < lon0:
                    mask_lon = ((obs_i['lon'] < lon0) & (obs_i['lon'] > lon1))
                else:
                    mask_lon = ((obs_i['lon'] > lon0) & (obs_i['lon'] < lon1))
                mask_lat = ((obs_i['lat'] > lat0) & (obs_i['lat'] < lat1))
                mask_time = ((obs_i['time'] > time0 - resolt)
                             & (obs_i['time'] < time1 + resolt))
                mask = (mask_time & mask_lat & mask_lon & mask_data)
               # print(obs_i['time'], time0, time1)
                for key in obs_i.keys():
                    if key not in obs.keys():
                        obs[key] = []
                    obs[key] = numpy.concatenate(([obs[key],
                                                   obs_i[key][mask]]))

    # #####################################
    # Spatial grid construction
    # #####################################
    grd = make_grid(lon0, lon1, dlon, lat0, lat1, dlat)

    # #####################################
    # Independant time loop for each analysis
    # #####################################
    enstime = numpy.arange(time0, time1 + dt, dt)
    time_model = numpy.arange(time0, time1, p.timestep)
    # Center at noon
    window = dt / 2.
    enstime = enstime + window

    for it, timeref in enumerate(enstime):
        iobs = numpy.where((obs['time'] > timeref - resolt)
                         & (obs['time'] < timeref + resolt))[0]
        make_oi(grd, obs, iobs, resols, timeref, resolt)
        # Interpolate model data to have a true reference
        indice_model = it
        model_data, out_var, list_file = read_model(p, it)
        interpolate_model(p, model_data, out_var, grd)

        pattern = os.path.join(p.outdatadir, '{}_l2d_'.format(config))
        save_l2d(pattern, timeref, window, time_unit, grd)


def save_l2d(filenc, timeref, window, time_unit, grd):
    dateformat = '%Y-%m-%dT%H:%M:%SZ'
    time_std = netCDF4.num2date(timeref - window, time_unit)
    time_start = time_std.strftime(format=dateformat)
    time_std = netCDF4.num2date(timeref + window, time_unit)
    time_end = time_std.strftime(format=dateformat)
    time_std = netCDF4.num2date(timeref, time_unit)
    time_middle = time_std.strftime(format=dateformat)
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
    rw.write_l2d(metadata, grd, ux_noerr=grd['ux_noerr'],
                 uy_noerr=grd['uy_noerr'], ux_obs=grd['ux_obs'],
                 uy_obs=grd['uy_obs'], ux_true=grd['ux_true'],
                 uy_true=grd['uy_true'])



def read_l2b(nfile):
    fid = netCDF4.Dataset(nfile)
    obs_i = {}
    obs_i['ux'] = numpy.array(fid.variables['ucur'][:]).flatten()
    obs_i['uy'] = numpy.array(fid.variables['vcur'][:]).flatten()
    if 'ur_true' in fid.variables.keys():
        obs_i['ur_true'] = numpy.array(fid.variables['ur_true'][:]).flatten()
    if 'ur_obs' in fid.variables.keys():
        obs_i['ur_obs'] = numpy.array(fid.variables['ur_obs'][:]).flatten()
    obs_i['lon'] = numpy.array(fid.variables['lon'][:]).flatten()
    #obs_i['lon'] = numpy.mod(obs_i['lon'] -180, 360) + 180
    obs_i['lon'] = numpy.mod(obs_i['lon'] + 360, 360)
    obs_i['lat'] = numpy.array(fid.variables['lat'][:]).flatten()
    obs_i['time'] = numpy.array(fid.variables['time'][:]).flatten()
    angle = numpy.array(fid.variables['radial_angle'][:]).flatten()
    obs_i['angle'] = numpy.mod(angle, 2 * numpy.pi)
    mask_data = ((obs_i['ux'] > -1000) & (obs_i['uy'] > -1000)
                 & (obs_i['ur_true'] > -1000) & (obs_i['ur_obs'] > -1000))
    time_unit = fid.variables['time'].units
    return obs_i, mask_data, time_unit


def make_grid(lon0, lon1, dlon, lat0, lat1, dlat):
    grd = {}
    grd['lon'] = numpy.arange(lon0, lon1 + dlon, dlon)
    grd['lat'] = numpy.arange(lat0, lat1 + dlat, dlat)
    grd['lon2'], grd['lat2'] = numpy.meshgrid(grd['lon'], grd['lat'])
    grd['nx'] = len(grd['lon'])
    grd['ny'] = len(grd['lat'])
    grd['ux_noerr'] = numpy.full((grd['ny'], grd['nx']), numpy.nan)
    grd['uy_noerr'] = numpy.full((grd['ny'], grd['nx']), numpy.nan)
    grd['ux_obs'] = numpy.full((grd['ny'], grd['nx']), numpy.nan)
    grd['uy_obs'] = numpy.full((grd['ny'], grd['nx']), numpy.nan)
    return grd


def make_oi(grd, obs, iobs, resols, timeref, resolt):
    # print(grd['nx'], grd['ny'], numpy.shape(grd['ux_noerr']))
    for j in range(grd['ny']):
        for i in range(grd['nx']):
          # TODO: to be optimized, especially for global, ...
          dist = 110. * (numpy.cos(numpy.deg2rad(grd['lat'][j]))**2
                         * (obs['lon'][iobs] - grd['lon2'][j, i])**2
                         + (obs['lat'][iobs] - grd['lat2'][j, i])**2)**0.5
          iiobs=numpy.where((dist < resols))[0]
          if len(iiobs)>=2:
              H = numpy.zeros((len(iiobs), 2))
              H[:, 0] = numpy.cos(obs['angle'][iobs][iiobs])
              H[:, 1] = numpy.sin(obs['angle'][iobs][iiobs])
              win_s = numpy.exp(-dist[iiobs]**2/(0.5*resols)**2)
              time_cen = obs['time'][iobs][iiobs] - timeref
              win_t = numpy.exp(-time_cen**2/(0.5 * resolt)**2) # exp window
              Ri = win_s * win_t
              RiH = numpy.tile(Ri, (2, 1)).T*H
              M = numpy.dot(H.T, RiH)
              if numpy.linalg.cond(M) < 1e3:
                  Mi = numpy.linalg.inv(M)
                  eta_true = numpy.dot(numpy.dot(Mi, RiH.T),
                                       obs['ur_true'][iobs][iiobs])
                  eta_obs = numpy.dot(numpy.dot(Mi, RiH.T),
                                      obs['ur_obs'][iobs][iiobs])
                  grd['uy_noerr'][j, i] = eta_true[1]
                  grd['ux_noerr'][j, i] = eta_true[0]
                  grd['uy_obs'][j, i] = eta_obs[1]
                  grd['ux_obs'][j, i] = eta_obs[0]


def read_model(p, ifile):
    model_data, list_file = mod.load_coordinate_model(p)
    model_data.read_coordinates(p)
    nfile = int(ifile)
    filename = os.path.join(p.indatadir, list_file[nfile])
    model_step_ctor = getattr(rw, p.model)
    out_var = {}
    for i in range(p.dim_time):
        model_step = model_step_ctor(p, ifile=(filename, ),
                                     list_input_var=p.list_input_var_l2d,
                                     time=i)
        model_step.read_var(p)
        for key in p.list_input_var_l2d.keys():
            if not key in out_var.keys():
                out_var[key] = model_step.input_var[key]/p.dim_time
            else:
                out_var[key] = (out_var[key]
                                + model_step.input_var[key] / p.dim_time)
    return model_data, out_var, list_file


def interpolate_model(p, model_data, model_var, grd):
    import pyresample as pr
    wrap_lon = pr.utils.wrap_longitudes
    geom = pr.geometry.SwathDefinition
    interp = mod.interpolate_irregular_pyresample
    lon = wrap_lon(grd['lon2'])
    list_key = {'ucur':'ux_true', 'vcur':'uy_true'}

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
    return grd
