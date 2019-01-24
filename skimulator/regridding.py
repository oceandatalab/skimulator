import numpy
import os
import sys
import glob
import scipy.interpolate
import datetime
import time
import skimulator.const as const
import skimulator.rw_data as rw
import skimulator.mod_tools as mod_tools
import skimulator.mod_run as mod
import skimulator.mod_parallel
import logging
import traceback
logger = logging.getLogger(__name__)


def make_obs(data, grid, obs, ind):
    #obs['vrt'] = obs['vrt'][ind]
    obs['vxt'] = numpy.array(data.ucur).flatten()[ind]
    obs['vyt'] = numpy.array(data.vcur).flatten()[ind]
    obs['vmodr'] = numpy.array(data.ur_true).flatten()[ind]
    obs['lon'] = numpy.array(data.lon).flatten()[ind]
    obs['lat'] = numpy.array(data.lat).flatten()[ind]
    obs['lon_nadir'] = numpy.array(data.lon_nadir[:])
    obs['lat_nadir'] = numpy.array(data.lat_nadir[:])
    obs['time'] = numpy.array(data.time).flatten()[ind] #+ (cycle-1)*self.tcycle 
    obs['time_nadir'] = numpy.array(data.time_nadir[:])
    obs['vindice'] = numpy.array(data.vindice).flatten()[ind]

    obs['dir'] = numpy.mod(grid.angle + numpy.pi/2, 2*numpy.pi).flatten()[ind]
    obs['angle'] = numpy.mod(grid.radial_angle, 2*numpy.pi).flatten()[ind]
    obs['al_nadir'] = grid.x_al_nadir
    obs['al'] = (numpy.array(grid.x_al)
              + numpy.tile(obs['al_nadir'], (obs['nbeam'], 1)).transpose())
    obs['al'] = obs['al'].flatten()[ind]
    obs['ac'] = numpy.array(grid.x_ac).flatten()[ind]
    return obs


def across_track(lon, lat, posting, max_ac, desc=False):
    npoints = 1
    nind = len(lon)
    SatDir = numpy.zeros((int(nind/npoints), 3))
    SatLoc = numpy.zeros((int((nind)/npoints), 3))

    s2cart = mod_tools.spher2cart(lon[::npoints], lat[::npoints])
    SatLoc[:, 0], SatLoc[:, 1], SatLoc[:, 2] = s2cart
    # Compute satellite direction (SatLoc is periodic)
    SatDir[1: -1, 0] = ((SatLoc[2:, 0] - SatLoc[: -2, 0])
                        / numpy.sqrt(SatLoc[1: -1, 0]**2
                        + SatLoc[1: -1, 1]**2 + SatLoc[1: -1, 2]**2))
    SatDir[1: -1, 1] = ((SatLoc[2:, 1] - SatLoc[: -2, 1])
                        / numpy.sqrt(SatLoc[1: -1, 0]**2
                        + SatLoc[1: -1, 1]**2 + SatLoc[1: -1, 2]**2))
    SatDir[1: -1, 2] = ((SatLoc[2:, 2] - SatLoc[: -2, 2])
                        / numpy.sqrt(SatLoc[1: -1, 0]**2
                        + SatLoc[1: -1, 1]**2 + SatLoc[1: -1, 2]**2))
    SatDir[-1, :] = SatDir[-2, :]
    SatDir[0, :] = SatDir[1, :]
    # Rotate from earth center around satellite direction to compute
    # swath points of angles between the borders of the swath in left
    # and right swath
    nhalfswath = int(max_ac/posting)
    new_lon = numpy.zeros((nind, 2*nhalfswath+1))
    new_lat = numpy.zeros((nind, 2*nhalfswath+1))
    ac = numpy.zeros((2*nhalfswath+1))
    for i in range(0, nind, npoints):
        for j in range(0, int(nhalfswath)+1):
            #j = int(nhalfswath) - j
            ac[nhalfswath + j] = j*posting
            ac[nhalfswath - j] = -j*posting
            R = mod_tools.rotationmat3D(float((j*posting)
                                        / (const.Rearth*10**-3)),
                                        SatDir[int(i/npoints), :])
            ObsLoc = numpy.dot(R, SatLoc[int(i/npoints)])
            cs = mod_tools.cart2spher(ObsLoc[0], ObsLoc[1], ObsLoc[2])
            new_lon[i, nhalfswath+j], new_lat[i, nhalfswath+j] = cs
            ObsLoc = numpy.dot(numpy.transpose(R), SatLoc[int(i/npoints)])
            cs = mod_tools.cart2spher(ObsLoc[0], ObsLoc[1], ObsLoc[2])
            new_lon[i, nhalfswath-j], new_lat[i, nhalfswath-j] = cs
    new_lon = new_lon[:, ::-1]
    new_lat = new_lat[:, ::-1]
    return new_lon, new_lat, ac


def make_grid(grid, obs, posting, desc=False):
    grd = {}

    if   obs['ac'].max() < 0:
        return None
    max_ac = abs(obs['ac']).max()
    # OI grid 
    grd['dal'] = posting
    # if desc is True:
    #    ind = numpy.where(obs['al_nadir'] > obs['al'].min())
    #    obs['al_nadir'] = obs['al_nadir'][ind]
    grd['al'] = numpy.arange(obs['al_nadir'].min(), obs['al_nadir'].max(),
                             grd['dal'])

    #grd['lon'] = scipy.interpolate.griddata(_alac, obs['lon'], _inalac,
    #                                        method='linear')
    #grd['lat'] = scipy.interpolate.griddata(_alac, obs['lat'], _inalac,
    #                                     method='linear')
    lon360 = numpy.mod(obs['lon_nadir'] + 180, 360) - 180
    lon = scipy.interpolate.griddata(obs['al_nadir'], obs['lon_nadir'],
                                      (grd['al']), method='linear')
    lat = scipy.interpolate.griddata(obs['al_nadir'], obs['lat_nadir'],
                                      (grd['al']), method='linear')
    grd['lon'], grd['lat'], grd['ac'] = across_track(lon, lat, posting,
                                                     max_ac + grd['dal'],
                                                     desc=desc)
    grd['ac2'], grd['al2'] = numpy.meshgrid(grd['ac'], grd['al'])
    grd['nal'], grd['nac'] = numpy.shape(grd['ac2'])
    grd['time'] = scipy.interpolate.griddata(obs['al_nadir'], obs['time_nadir'],
                                             (grd['al']), method='linear')
    grd['angle'] = numpy.full((grd['nal'], grd['nac']), numpy.nan)
    grd['angle'][:, :-1]  = numpy.angle((grd['lon'][:, 1:] - grd['lon'][:, :-1])
                             + 1j * (grd['lat'][:, 1:] - grd['lat'][:, :-1]))
    grd['angle'][:, -1]= grd['angle'][:, -2]


    grd['vobsal'] = numpy.full((grd['nal'], grd['nac']), numpy.nan)
    grd['vobsac'] = numpy.full((grd['nal'], grd['nac']), numpy.nan)
    grd['vmodal'] = numpy.full((grd['nal'], grd['nac']), numpy.nan)
    grd['vmodac'] = numpy.full((grd['nal'], grd['nac']), numpy.nan)
    return grd


def perform_oi_1(grd, obs, resol, desc=False):
    for j in range(grd['nac']):
        for i in range(grd['nal']):
            dist = numpy.sqrt((obs['ac'] - grd['ac2'][i, j])**2
                              + (obs['al'] - grd['al2'][i, j])**2)
            ios = numpy.where((dist < resol))[0]
            if len(ios) >= 4:
                H = numpy.zeros((len(ios), 2))
                H[:, 0] = numpy.cos(obs['dir'][ios])
                H[:, 1] = numpy.sin(obs['dir'][ios])
                std_err = numpy.ones((len(ios)))
                # rectangular filtering window for now...
                # Ri = std_err**-2
                Ri=numpy.exp(-dist[ios]**2/(0.5*resol)**2) # exp window
                RiH = numpy.tile(Ri, (2, 1)).T * H
                M = numpy.dot(H.T, RiH)
                Mi = numpy.linalg.inv(M)
                eta_obs = numpy.dot(numpy.dot(Mi, RiH.T), obs['vobsr'][ios])
                eta_mod = numpy.dot(numpy.dot(Mi, RiH.T), obs['vmodr'][ios])
                j2 = j
                #if desc is True:
                #    j2 = grd['nac'] - 1 - j
                grd['vobsal'][i, j2]=eta_obs[0]
                grd['vobsac'][i, j2]=eta_obs[1]
                grd['vmodal'][i, j2]=eta_mod[0]
                grd['vmodac'][i, j2]=eta_mod[1]
    return grd

"""
def perform_oi_2(grd, obs, resol):
    # - In parameter file ## TODO -
    # Number of pixel (resolution for healpix)
    nside = 256
    # Number of diamonds for healpix
    ndiam = 12
    ntotpixel = nside * nside * ndiam
    # Conditionning threshold
    thresh_cond = 10
    ph = 2 * numpy.pi - numpy.deg2rad(lon)
    th = numpy.pi / 2 - numpy.deg2rad(lat)
    pidx = heal.ang2pix(nside, th, ph)

    for i in range(nbeam):
        for j in range(ndata):
            if ur[j, i] > -1E9:
                ip = pidx[j,i]
                # compute imulated model
                im[ip, 1] += u[j, i]
                im[ip, 2] += v[j,i]
                nim[ip] += 1
                # compute covariance(s) model
                co = numpy.cos(rangle[j,i])
                si = numpy.sin(rangle[j,i])
                w = ww[j,i]
                cov[ip, 0, 0] += co * co
                cov[ip, 1, 0] += si * co
                cov[ip, 0, 1] += si * co
                cov[ip, 1, 1] += si * si

                cov2[ip, 0, 0] += w * co * co
                cov2[ip, 1, 0] += w * si * co
                cov2[ip, 0, 1] += w * si * co
                cov2[ip, 1, 1] += w * si * si

                # compute data vector model
                vec[ip, 0] += co * ur[j,i]
                vec[ip, 1] += si * ur[j,i]

                # compute data noise vector model
                vec2[ip, 0] += w* co * uro[j,i]
                vec2[ip, 1] += w * si * uro[j,i]

                # compute doppler projection
                for k in range(3):
                    vecdop[k, ip, 0] += w * co * tdop[j,i,k]
                    vecdop[k, ip, 1] += w * si * tdop[j,i,k]
"""


def read_model(p, indice):
    model_data, list_file = mod.load_coordinate_model(p)
    model_data.read_coordinates(p)
    list_model_step = []
    for ifile in indice:
        nfile = int(ifile /p.dim_time)
        filetime = ifile - nfile * p.dim_time
        _tmpfilename = list_file[nfile]
        filename = os.path.join(p.indatadir, _tmpfilename)
        model_step_ctor = getattr(rw, p.model)
        model_step = model_step_ctor(p, ifile=(filename, ),
                                     list_input_var=p.list_input_var,
                                     time=filetime)
        model_step.read_var(p)
        list_model_step.append(model_step)
    return model_data, list_model_step, list_file


def interpolate_model(p, model_data, list_model_step, grd, list_obs,
                      desc=False):
    import pyresample as pr
    wrap_lon = pr.utils.wrap_longitudes
    geom = pr.geometry.SwathDefinition
    interp = mod.interpolate_irregular_pyresample
    lon = wrap_lon(grd['lon'])
    list_key = {'ucur':'u_model', 'vcur':'v_model'}

    for ikey, okey in list_key.items():
        grd[okey] = numpy.full(numpy.shape(grd['lat']), numpy.nan)
        if len(p.list_input_var[ikey]) > 2:
            grid_number = p.list_input_var[ikey][2]
        else:
            grid_number = 0
        _lon = wrap_lon(model_data.vlon[grid_number])
        _lat = model_data.vlat[grid_number]
        if len(numpy.shape(_lon)) <= 1:
            _lon, _lat = numpy.meshgrid(_lon, _lat)
        swath_def = geom(lons=_lon, lats=_lat)
        for i in range(len(list_model_step)):
            model_step = list_model_step[i]
            if desc is False:
                ind_lat = numpy.where(grd['lat']> list_obs[i])
            else:
                ind_lat = numpy.where(grd['lat']< list_obs[i])
            model_step.read_var(p)
            grid_def = geom(lons=lon[ind_lat], lats=grd['lat'][ind_lat])
            #grid_def = pr.geometry.SwathDefinition(lons=lon, #[ind_lat],
            #                                     lats=grd['lat']) #[ind_lat])
            var = model_step.input_var[ikey]
            _tmp = interp(swath_def, var, grid_def, p.resol,
                          interp_type=p.interpolation)
            grd[okey][ind_lat] = _tmp
    return grd



def write_l2(outfile, grd, obs, cycle, passn, firsttime):
    if os.path.exists(outfile):
        os.remove(outfile)
    metadata = {}
    metadata['file'] = outfile
    dateformat = '%Y-%m-%dT%H:%M:%SZ'
    time_model = datetime.datetime.strptime(firsttime, '%Y-%m-%dT%H:%M:%SZ')
    grdtime0 = numpy.nanmin(obs['time_nadir'])
    grdtime1 = numpy.nanmax(obs['time_nadir'])
    if numpy.isnan(grdtime0):
        grdtime0 = 0
    if numpy.isnan(grdtime1):
        grdtime1 = 0
    day = numpy.floor(grdtime0)
    seconds = (grdtime0 - day) * 86400
    time0 = time_model + datetime.timedelta(day, seconds)
    day = numpy.floor(grdtime1)
    seconds = (grdtime1 - day) * 86400
    time1 = time_model + datetime.timedelta(day, seconds)

    metadata['time_coverage_start'] = time0.strftime(format=dateformat)
    metadata['time_coverage_end'] = time1.strftime(format=dateformat)
    metadata['cycle'] = cycle
    metadata['pass'] = passn
    metadata['first_time'] = firsttime
    #geolocation = {}
    #geolocation['lon']
    rw.write_l2c(metadata, grd, u_ac_obs=grd['vobsac'], u_al_obs=grd['vobsal'],
                 u_ac_model=grd['vmodac'], u_al_model=grd['vmodal'],
                 angle=grd['angle'], u_obs=grd['vobsx'], v_obs=grd['vobsy'],
                 u_model=grd['vmodx'], v_model=grd['vmody'],
                 u_true=grd['u_model'], v_true=grd['v_model'],
                 u_ac_true=grd['vtrueac'], u_al_true=grd['vtrueal'],
                 )


def run_l2c(p, die_on_error=False):
    timestart = datetime.datetime.now()
    pattern = os.path.join(p.outdatadir, '{}_c*'.format(p.config))
    list_file = glob.glob(pattern)
    gpath = os.path.join(p.outdatadir, '{}_grid'.format(p.config))
    mod_tools.initialize_parameters(p)
    iterate = 0

    # - Loop on SKIM grid files
    jobs = []
    p2 = mod_tools.todict(p)
    for ifile in list_file:
        jobs.append([ifile, p2, gpath])
    ok = False
    try:
        ok = make_skim_l2c(p.proc_count, jobs, die_on_error, p.progress_bar)
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
    rw.write_params(p, os.path.join(p.outdatadir,
                                         'skim_l2c.output'))
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


def make_skim_l2c(_proc_count, jobs, die_on_error, progress_bar):
    """ Compute SWOT-like data for all grids and all cycle, """
    # - Set up parallelisation parameters
    proc_count = min(len(jobs), _proc_count)

    status_updater = mod_tools.update_progress_multiproc
    jobs_manager = skimulator.mod_parallel.JobsManager(proc_count,
                                                       status_updater,
                                                       exc_formatter,
                                                       err_formatter)
    ok = jobs_manager.submit_jobs(worker_method_l2c, jobs, die_on_error,
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


def worker_method_l2c(*args, **kwargs):
    msg_queue, ifile, p2, gpath = args[:]
    p = mod_tools.fromdict(p2)

    passn = int(ifile[-6:-3])
    if passn % 2 == 0:
        desc = True
    else:
        desc = False
    cycle = int(ifile[-10:-8])
    fileg = '{}_p{:03d}.nc'.format(gpath, passn)
    data = rw.Sat_SKIM(ifile=ifile)
    grid = rw.Sat_SKIM(ifile=fileg)
    data.load_data(p, ur_true=[], ur_obs=[], ucur=[],
                   vcur=[], time=[], lon_nadir=[], lat_nadir=[],
                   lon=[], lat=[], time_nadir=[], vindice=[])
    grid.load_swath(p, radial_angle=[], angle=[], x_al=[], x_al_nadir=[],
                    x_ac=[])

    obs = {}
    test = numpy.array(data.ucur)
    nil, nbeams = numpy.shape(test)
    sbeam_incid = numpy.zeros((nil, nbeams))
    ### TODO Change this
    obs['vobsr'] = numpy.array(data.ur_obs)
    obs['nsamp'], obs['nbeam'] = numpy.shape(obs['vobsr'])
    obs['vobsr'] = obs['vobsr'].flatten()
    ind = numpy.where((obs['vobsr'] > -1000))[0]
    obs['vobsr'] = obs['vobsr'][ind]
    if len(ind) > 2 and len(data.lon_nadir) >2:
        try:
            obs = make_obs(data, grid, obs, ind)
            grd = make_grid(grid, obs, p.posting, desc=desc)
            ## TODO proof error
            if grd is None:
                return

            # OI
            grd = perform_oi_1(grd, obs, p.resol, desc=desc)
        except:
            print(passn)
            return
        vindice = obs['vindice']
        diff_indice = vindice[1:] - vindice[:-1]
        ind = numpy.where(diff_indice != 0)[0]
        first_lat = numpy.min(grd['lat'])
        sign_uv = 1
        if desc is True:
            first_lat = numpy.max(grd['lat'])
            sign_uv = -1

        if ind.any():
            vindice = [obs['vindice'][0], obs['vindice'][ind[0] + 1]]
            ind_lat = [first_lat, obs['lat'][ind[0] + 1]]
        else:
            vindice = [obs['vindice'][0],]
            ind_lat = [first_lat,]
        #vindice = [obs['vindice'][0],]
        #ind_lat = [first_lat,]
        model_data, model_step, list_file2 = read_model(p, vindice)
        grd = interpolate_model(p, model_data, model_step, grd, ind_lat,
                                desc=desc)
        ac_thresh = p.ac_threshold
        grd['vobsac'][numpy.abs(grd['ac2']) < ac_thresh] = numpy.nan
        grd['vobsx'] = sign_uv * (grd['vobsac'] * numpy.cos(grd['angle'])
                        + grd['vobsal'] * numpy.cos(grd['angle'] + numpy.pi/2))
        grd['vobsy'] = sign_uv * (grd['vobsac'] * numpy.sin(grd['angle'])
                        + grd['vobsal'] * numpy.sin(grd['angle'] + numpy.pi/2))
        grd['vmodac'][numpy.abs(grd['ac2']) < ac_thresh] = numpy.nan
        grd['vmodx'] = sign_uv * (grd['vmodac'] * numpy.cos(grd['angle'])
                        + grd['vmodal'] * numpy.cos(grd['angle'] + numpy.pi/2))
        grd['vmody'] = sign_uv * (grd['vmodac'] * numpy.sin(grd['angle'])
                        + grd['vmodal'] * numpy.sin(grd['angle'] + numpy.pi/2))
        grd['vtrueac'] = sign_uv *(grd['u_model']*numpy.cos(grd['angle'])
                                   + grd['v_model'] * numpy.sin(grd['angle']))
        grd['vtrueal'] = sign_uv * (-grd['u_model']*numpy.sin(grd['angle'])
                                    + grd['v_model']*numpy.cos(grd['angle']))
        mask = ((grd['u_model'] == 0) | (grd['v_model'] == 0)
                | (abs(grd['u_model']) > 10) | (abs(grd['v_model']) > 10))
        grd['u_model'][mask] = numpy.nan
        grd['v_model'][mask] = numpy.nan
        grd['vtrueac'][mask] = numpy.nan
        grd['vtrueal'][mask] = numpy.nan
        grd['vobsac'][mask] = numpy.nan
        grd['vobsal'][mask] = numpy.nan
        grd['vobsx'][mask] = numpy.nan
        grd['vobsy'][mask] = numpy.nan
        grd['vmodac'][mask] = numpy.nan
        grd['vmodal'][mask] = numpy.nan
        grd['vmodx'][mask] = numpy.nan
        grd['vmody'][mask] = numpy.nan
        pattern_out = '{}_L2C_c{:02d}_p{:03d}.nc'.format(p.config, cycle, passn)
        outfile = os.path.join(p.outdatadir, pattern_out)
        write_l2(outfile, grd, obs, cycle, passn, p.first_time)
    msg_queue.put((os.getpid(), ifile, None, None))
