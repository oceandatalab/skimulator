import numpy
import os
import glob
import scipy.interpolate
import datetime
import skimulator.const as const
import skimulator.rw_data as rw
import skimulator.mod_tools as mod_tools
import skimulator.mod_run as mod
import logging
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
            ios = numpy.where((dist < 0.5 * resol))[0]
            if len(ios) >= 4:
                H = numpy.zeros((len(ios), 2))
                H[:, 0] = numpy.cos(obs['dir'][ios])
                H[:, 1] = numpy.sin(obs['dir'][ios])
                std_err = numpy.ones((len(ios)))
                # rectangular filtering window for now...
                Ri = std_err**-2
                RiH = numpy.tile(Ri, (2, 1)).T * H
                M = numpy.dot(H.T, RiH)
                Mi = numpy.linalg.inv(M)
                eta_obs = numpy.dot(numpy.dot(Mi, RiH.T), obs['vobsr'][ios])
                eta_mod = numpy.dot(numpy.dot(Mi, RiH.T), obs['vmodr'][ios])
                j2 = j
                if desc is True:
                    j2 = grd['nac'] - 1 - j
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
        _tmpfilename = list_file[nfile].split(',')
        if len(_tmpfilename) > 1:
            filename_u = os.path.join(p.indatadir, _tmpfilename[0])
            filename_v = os.path.join(p.indatadir, _tmpfilename[1])
        else:
            filename_u = os.path.join(p.indatadir, _tmpfilename[0])
            filename_v = os.path.join(p.indatadir, _tmpfilename[0])

        model_step_ctor = getattr(rw, p.model)
        model_step = model_step_ctor(p, ifile=(filename_u, filename_v),
                                     time=filetime)
        model_step.read_var(p)
        list_model_step.append(model_data)
    return model_data, list_model_step, list_file


def interpolate_model(p, model_data, list_model_step, grd, list_obs,
                      desc=False):
    import pyresample as pr
    model_data.vlonu = pr.utils.wrap_longitudes(model_data.vlonu)
    lon = pr.utils.wrap_longitudes(grd['lon'])
    if len(numpy.shape(model_data.vlonu)) <= 1:
        model_data.vlonu, model_data.vlatu = numpy.meshgrid(model_data.vlonu,
                                                            model_data.vlatu)
        if p.lonu == p.lonv:
            model_data.vlonv = model_data.vlonu
            model_data.vlatv = model_data.vlatu
        else:
            model_data.vlonv, model_data.vlatv = numpy.meshgrid(model_data.vlonu,
                                                                model_data.vlatu)
    swath_defu = pr.geometry.SwathDefinition(lons=model_data.vlonu,
                                             lats=model_data.vlatu)
    swath_defv = pr.geometry.SwathDefinition(lons=model_data.vlonv,
                                             lats=model_data.vlatv)
    grd['u_model'] = numpy.full(numpy.shape(grd['lat']), numpy.nan)
    grd['v_model'] = numpy.full(numpy.shape(grd['lat']), numpy.nan)
    for i in range(len(list_model_step)):
        model_step = list_model_step[i]
        if desc is False:
            ind_lat = numpy.where(grd['lat']> list_obs[i])
        else:
            ind_lat = numpy.where(grd['lat']< list_obs[i])
        model_step.read_var(p)
        grid_def = pr.geometry.SwathDefinition(lons=lon[ind_lat],
                                               lats=grd['lat'][ind_lat])
        var = model_step.input_var['ucur']
        _tmp = mod.interpolate_irregular_pyresample(swath_defu, var, grid_def,
                                                    p.posting,
                                                    interp_type=p.interpolation)
        grd['u_model'][ind_lat] = _tmp
        var = model_step.input_var['vcur']
        _tmp = mod.interpolate_irregular_pyresample(swath_defu, var,
                                                    grid_def, p.posting,
                                                    interp_type=p.interpolation)
        grd['v_model'][ind_lat] = _tmp
    return grd



def write_l2(outfile, grd, cycle, passn, firsttime):
    if os.path.exists(outfile):
        os.remove(outfile)
    metadata = {}
    metadata['file'] = outfile
    dateformat = '%Y-%m-%dT%H:%M:%SZ'
    time_model = datetime.datetime.strptime(firsttime, '%Y-%m-%dT%H:%M:%SZ')
    time0 = time_model + datetime.timedelta(0, grd['time'][0])
    time1 = time_model + datetime.timedelta(0, grd['time'][-1])

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
                 u_ac_true=grd['vtrueac'], u_al_true=grd['vtrueal'])


def run_l2c(p):
    pattern = os.path.join(p.outdatadir, '{}_c*'.format(p.config))
    list_file = glob.glob(pattern)
    gpath = os.path.join(p.outdatadir, '{}_grid'.format(p.config))
    mod_tools.initialize_parameters(p)
    iterate = 0
    for ifile in list_file:
        progress = iterate / numpy.float(len(list_file))
        arg1 = os.path.basename(ifile)
        arg2 = '{} {} {}'.format(progress, iterate, len(list_file))
        arg2 = ''
        mod_tools.update_progress(progress, arg1, arg2)
        iterate += 1
        passn = int(ifile[-6:-3])
        if passn % 2 == 0:
            desc = True
        else:
            desc = False
        cycle = int(ifile[-10:-8])
        fileg = '{}_p{:03d}.nc'.format(gpath, passn)
        data = rw.Sat_SKIM(ifile=ifile)
        grid = rw.Sat_SKIM(ifile=fileg)
        data.load_data(p, ur_true=[], ur_obs=[], instr=[], ucur=[],
                       vcur=[], time=[], lon_nadir=[], lat_nadir=[],
                       lon=[], lat=[], time_nadir=[], vindice=[])
        grid.load_swath(p, radial_angle=[], angle=[], x_al=[], x_al_nadir=[],
                        x_ac=[])

        obs = {}
        test = numpy.array(data.ucur)
        nil, nbeams=numpy.shape(test)
        sbeam_incid=numpy.zeros((nil, nbeams))
        ### TODO Change this
        obs['vobsr'] = numpy.array(data.ur_obs)
        obs['nsamp'], obs['nbeam'] = numpy.shape(obs['vobsr'])
        obs['vobsr'] = obs['vobsr'].flatten()
        ind = numpy.where((obs['vobsr'] > -1000))[0]
        obs['vobsr'] = obs['vobsr'][ind]
        if len(ind) > 1 and len(data.lon_nadir) >2:
            obs = make_obs(data, grid, obs, ind)
            grd = make_grid(grid, obs, p.posting, desc=desc)
            if grd is None:
                continue

            # OI
            grd = perform_oi_1(grd, obs, p.resol, desc=desc)
            vindice = obs['vindice']
            diff_indice = vindice[1:] - vindice[:-1]
            ind = numpy.where(diff_indice != 0)[0]
            first_lat = numpy.min(grd['lat'])
            if desc is True:
                first_lat = numpy.max(grd['lat'])

            if ind.any():
                vindice = [obs['vindice'][0], obs['vindice'][ind[0] + 1]]
                ind_lat = [first_lat, obs['lat'][ind[0] + 1]]
                print(ind, ind[0])
            else:
                vindice = [obs['vindice'][0],]
                ind_lat = [first_lat,]
            model_data, model_step, list_file2 = read_model(p, vindice)
            grd = interpolate_model(p, model_data, model_step, grd, ind_lat,
                                    desc=desc)
            ac_thresh = p.ac_threshold
            grd['vobsac'][numpy.abs(grd['ac2']) < ac_thresh] = numpy.nan
            grd['vobsx'] = (grd['vobsac'] * numpy.cos(grd['angle'])
                         + grd['vobsal'] * numpy.cos(grd['angle'] + numpy.pi/2))
            grd['vobsy'] = (grd['vobsac'] * numpy.sin(grd['angle'])
                         + grd['vobsal'] * numpy.sin(grd['angle'] + numpy.pi/2))
            grd['vmodac'][numpy.abs(grd['ac2']) < ac_thresh] = numpy.nan
            grd['vmodx'] = (grd['vmodac'] * numpy.cos(grd['angle'])
                         + grd['vmodal'] * numpy.cos(grd['angle'] + numpy.pi/2))
            grd['vmody'] = (grd['vmodac'] * numpy.sin(grd['angle'])
                         + grd['vmodal'] * numpy.sin(grd['angle'] + numpy.pi/2))
            grd['vtrueac'] = (grd['u_model']*numpy.cos(grd['angle'])
                              + grd['v_model'] * numpy.sin(grd['angle']))
            grd['vtrueal'] = (-grd['u_model']*numpy.sin(grd['angle'])
                              + grd['v_model']*numpy.cos(grd['angle']))
            pattern_out = '{}_L2C_c{:02d}_p{:03d}.nc'.format(p.config, cycle, passn)
            outfile = os.path.join(p.outdatadir, pattern_out)
            write_l2(outfile, grd, cycle, passn, p.first_time)
    __ = mod_tools.update_progress(1, 'All passes have been processed', '')
    logger.info("\n Simulated skim files have been written in "
                "{}".format(p.outdatadir))
