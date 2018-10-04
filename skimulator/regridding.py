import numpy
import os
import glob
import scipy.interpolate
import datetime
import skimulator.const as const
import skimulator.rw_data as rw


def make_obs(data, grid, obs, ind):
    #obs['vrt'] = obs['vrt'][ind]
    obs['vxt'] = numpy.array(data.u_model).flatten()[ind]
    obs['vyt'] = numpy.array(data.v_model).flatten()[ind]
    obs['vmodr'] = numpy.array(data.ur_model).flatten()[ind]
    obs['lon'] = numpy.array(data.lon).flatten()[ind]
    obs['lat'] = numpy.array(data.lat).flatten()[ind]
    obs['time'] = numpy.array(data.time).flatten()[ind] #+ (cycle-1)*self.tcycle 
    obs['time_nadir'] = numpy.array(data.time_nadir[:])

    obs['dir'] = numpy.mod(grid.angle + numpy.pi/2, 2*numpy.pi).flatten()[ind]
    obs['angle'] = numpy.mod(grid.radial_angle, 2*numpy.pi).flatten()[ind]
    obs['al_nadir'] = grid.x_al_nadir
    obs['al'] = (numpy.array(grid.x_al)
              + numpy.tile(obs['al_nadir'], (obs['nbeam'], 1)).transpose())
    obs['al'] = obs['al'].flatten()[ind]
    obs['ac'] = numpy.array(grid.x_ac).flatten()[ind]
    return obs


def make_grid(grid, obs, posting):
    grd = {}

    # OI grid 
    grd['dal'] = posting
    grd['ac'] = numpy.arange(0, obs['ac'].max(), grd['dal'])
    grd['ac']= numpy.concatenate((-grd['ac'][1:][::-1], grd['ac']))
    grd['al'] = numpy.arange(obs['al'].min(),obs['al'].max(),grd['dal'])
    grd['ac2'], grd['al2'] = numpy.meshgrid(grd['ac'],grd['al'])
    grd['nal'], grd['nac'] = numpy.shape(grd['ac2'])

    # Lucile: à faire directement avec la projection lon,lat en xal,xac
    _alac = numpy.array([obs['ac'], obs['al']]).transpose()
    _inalac = (grd['ac2'], grd['al2'])
    grd['lon'] = scipy.interpolate.griddata(_alac, obs['lon'], _inalac,
                                            method='linear')
    grd['lat'] = scipy.interpolate.griddata(_alac, obs['lat'], _inalac,
                                         method='linear')
    grd['time'] = scipy.interpolate.griddata(obs['al_nadir'], obs['time_nadir'],
                                             (grd['al']), method='linear')
    grd['angle'] = numpy.full((grd['nal'], grd['nac']), numpy.nan)
    for j in range(grd['nac'] - 1):
      for i in range(grd['nal']):
        _complex = (grd['lon'][i, j + 1] - grd['lon'][i, j]
                    + 1j * (grd['lat'][i, j + 1] - grd['lat'][i, j]))
        grd['angle'][i, j] = numpy.angle(_complex)
    grd['angle'][:, -1]= grd['angle'][:, -2]



    grd['vobsal'] = numpy.full((grd['nal'],grd['nac']), numpy.nan)
    grd['vobsac'] = numpy.full((grd['nal'],grd['nac']), numpy.nan)
    grd['vmodal'] = numpy.full((grd['nal'],grd['nac']), numpy.nan)
    grd['vmodac'] = numpy.full((grd['nal'],grd['nac']), numpy.nan)
    return grd


def perform_oi_1(grd, obs, resol):
    for j in range(grd['nac']):
        for i in range(grd['nal']):
            dist = numpy.sqrt((obs['ac'] - grd['ac2'][i, j])**2
                              + (obs['al'] - grd['al2'][i, j])**2)
            ios = numpy.where((dist < 0.5 * resol))[0]
            if len(ios) >= 3:
                H = numpy.zeros((len(ios), 2))
                H[:, 0] = numpy.cos(obs['dir'][ios])
                H[:, 1] = numpy.sin(obs['dir'][ios])
                std_err = numpy.ones((len(ios)))
                # rectangular filtering window for now...
                Ri = std_err**-2
                RiH = numpy.tile(Ri, (2, 1)).T * H
                M = numpy.dot(H.T, RiH)
                #try:
                Mi = numpy.linalg.inv(M)
                #except:
                #    grd['vobsal'][i, j]=eta_obs[0]
                #    grd['vobsac'][i, j]=eta_obs[1]
                #    grd['vmodal'][i, j]=eta_mod[0]
                #    grd['vmodac'][i, j]=eta_mod[1]
                #    print(i, j, M)
                    #continue
                    #print(M)
                #    import sys
                #    sys.exit(1)
                eta_obs = numpy.dot(numpy.dot(Mi, RiH.T), obs['vobsr'][ios])
                eta_mod = numpy.dot(numpy.dot(Mi, RiH.T), obs['vmodr'][ios])
                grd['vobsal'][i, j]=eta_obs[0]
                grd['vobsac'][i, j]=eta_obs[1]
                grd['vmodal'][i, j]=eta_mod[0]
                grd['vmodac'][i, j]=eta_mod[1]
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



def write_l2(outfile, grd, cycle, passn):
                grd['vobsal'][i, j]=eta_obs[0]
                grd['vobsal'][i, j]=eta_obs[0]
                grd['vobsac'][i, j]=eta_obs[1]
                grd['vmodal'][i, j]=eta_mod[0]
                grd['vmodac'][i, j]=eta_mod[1]
                grd['vobsac'][i, j]=eta_obs[1]
                grd['vmodal'][i, j]=eta_mod[0]
                grd['vmodac'][i, j]=eta_mod[1]
    if os.path.exists(outfile):
        os.remove(outfile)
    metadata = {}
    metadata['file'] = outfile
    metadata['time_coverage_start'] = '20170101'
    metadata['time_coverage_end'] = '20170101'
    metadata['cycle'] = cycle
    metadata['pass'] = passn
    #geolocation = {}
    #geolocation['lon']
    rw.write_l2c(metadata, grd, u_ac=grd['vobsac'], u_al=grd['vobsal'],
                 angle=grd['angle'], u_obs=grd['vobsx'], v_obs=grd['vobsy'],
                 u_mod=grd['vmodx'], v_mod=grd['vmody'])


def run_l2c(p):
    pattern = os.path.join(p.outdatadir, '{}_c*'.format(p.config))
    list_file = glob.glob(pattern)
    gpath = os.path.join(p.outdatadir, '{}_grid'.format(p.config))

    p.resol=40. # km
    p.posting=5. # km

    for ifile in list_file:
        passn = int(ifile[-6:-3])
        cycle = int(ifile[-10:-8])
        print(ifile, passn, cycle)
        fileg = '{}_p{:03d}.nc'.format(gpath, passn)
        data = rw.Sat_SKIM(ifile=ifile)
        grid = rw.Sat_SKIM(ifile=fileg)
        data.load_data(p, ur_model=[], ur_obs=[], instr=[], u_model=[],
                       v_model=[], time=[],
                       lon=[], lat=[], time_nadir=[])
        grid.load_swath(p, radial_angle=[], angle=[], x_al=[], x_al_nadir=[],
                        x_ac=[])

        obs = {}
        test = numpy.array(data.u_model)
        nil, nbeams=numpy.shape(test)
        sbeam_incid=numpy.zeros((nil, nbeams))
        ### TODO Change this
        obs['vobsr'] = numpy.array(data.ur_model)
        obs['nsamp'], obs['nbeam'] = numpy.shape(obs['vobsr'])
        obs['vobsr'] = obs['vobsr'].flatten()
        ind = numpy.where((obs['vobsr'] > -1000))[0]
        if len(ind) > 1:
            obs = make_obs(data, grid, obs, ind)
            grd = make_grid(grid, obs, p.posting)

            # OI
            grd = perform_oi_1(grd, obs, p.resol)
            grd['vobsac'][numpy.abs(grd['ac2']) < 20] = numpy.nan

            grd['vobsx'] = (grd['vobsac'] * numpy.cos(grd['angle'])
                         + grd['vobsal'] * numpy.cos(grd['angle'] + numpy.pi/2))
            grd['vobsy'] = (grd['vobsac'] * numpy.sin(grd['angle'])
                         + grd['vobsal'] * numpy.sin(grd['angle'] + numpy.pi/2))
            grd['vmodac'][numpy.abs(grd['ac2']) < 20] = numpy.nan

            grd['vmodx'] = (grd['vmodac'] * numpy.cos(grd['angle'])
                         + grd['vmodal'] * numpy.cos(grd['angle'] + numpy.pi/2))
            grd['vmody'] = (grd['vmodac'] * numpy.sin(grd['angle'])
                         + grd['vmodal'] * numpy.sin(grd['angle'] + numpy.pi/2))
            pattern_out = '{}_L2C_c{:02d}_p{:03d}.nc'.format(p.config, cycle, passn)
            outfile = os.path.join(p.outdatadir, pattern_out)
            write_l2(outfile, grd, cycle, passn)
