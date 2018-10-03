import numpy
import os
import glob
import netCDF4 as nc
from scipy.interpolate import griddata
from netCDF4 import Dataset
import datetime
import skimulator.const as const
import skimulator.rw_data as rw


def make_obs(data, grid, obs, ind):
    obs['vrt'] = obs['vrt'][ind]
    obs['vxt'] = numpy.array(data.u_model).flatten()[ind]
    obs['vyt'] = numpy.array(data.v_model).flatten()[ind]
    obs['vr'] = numpy.array(data.ur_obs).flatten()[ind]
    obs['lon'] = numpy.array(data.lon).flatten()[ind]
    obs['lat'] = numpy.array(data.lat).flatten()[ind]
    obs['time'] = numpy.array(data.time).flatten()[ind] #+ (cycle-1)*self.tcycle 
    time_nadir = numpy.array(data.time_nadir][:])

    obs['dir'] = numpy.mod(grid.angle + numpy.pi/2, 2*numpy.pi).flatten()[ind]
    obs['angle'] = numpy.mod(grid.radial_angle, 2*numpy.pi).flatten()[ind]
    al_nadir = grid.x_al_nadir
    obs['al'] = (numpy.array(grid.x_al)
              + numpy.tile(al_nadir, (obs['nbeam'], 1)).transpose())
    obs['al'] = obs['al'].flatten()[ind]
    obs['ac'] = numpy.array(grid.x_ac).flatten()[ind]
    return obs

def make_grid(grid, obs):
    grd = {}

    # OI grid 
    grd['dal'] = posting
    grd['ac'] = numpy.arange(0, obs['ac'].max(), grd['dal'])
    grd['ac']= numpy.concatenate((-grd['ac'][1:][::-1], grd]'ac']))
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
    grd['time'] = scipy.interpolate.griddata(al_nadir, time_nadir, (grd['al']),
                                          method='linear')
    grd['angle'] = numpy.full((grd['nal'], grd['nac']), numpy.nan)
    for j in range(grd['nac'] - 1):
      for i in range(grd['nal']):
        _complex = (grd['lon'][i, j + 1] - grd['lon'][i, j]
                    + 1j * (grd['lat'][i, j + 1] - grd['lat'][i, j]))
        grd['angle'][i, j] = numpy.angle(_complex)
    grd['angle'][:, -1]= grd['angle'][:, -2]



    grd['val'] = numpy.full((grd.nal,grd.nac), numpy.nan)
    grd['vac'] = numpy.full((grd.nal,grd.nac), numpy.nan)
    return grd


def perform_oi(grd, obs):
    for j in range(grd['nac']):
        for i in range(grd['nal']):
            dist = numpy.sqrt((obs['ac'] - grd['ac2'][i, j])**2
                              + (obs['al'] - grd['al2'][i, j])**2)
            ios = numpy.where((dist < 0.5 * resol))[0]
            if len(ios) >= 2:
                H = numpy.zeros((len(ios), 2))
                H[:, 0] = numpy.cos(obs['dir'][ios])
                H[:, 1] = numpy.sin(obs['dir'][ios])
                std_err = numpy.ones((len(ios)))
                # rectangular filtering window for now...
                Ri = std_err**-2
                RiH = numpy.tile(Ri, (2, 1)).T * H
                M = numpy.dot(H.T, RiH)
                Mi = numpy.linalg.inv(M)
                eta = numpy.dot(numpy.dot(Mi, RiH.T), obs['vr'][ios])
                grd['val'][i, j]=eta[0]
                grd['vac'][i, j]=eta[1]
    return grd


def write_l2(outfile, grd):
    if os.path.exists(outfile):
        os.remove(outfile)
    metadata = {}
    metadata['file'] = outfile
    geolocation = {}
    geolocation['lon']
    rw_data.write_l2c(metadata, grd, u_ac=grd['vac'], u_al=grd['val'],
                      angle=grd['angle'], u_obs=grd['vx'], v_obs=grd['vy'])

def run_l2c(p):
    pattern = os.path.join(p.outdatadir, '{}_*'.format(p.config))
    list_file = glob.glob(pattern)
    gpath = os.path.join(p.outdatadir, '{}_grid'.format(p.config))

    resol=40. # km
    posting=5. # km

    for ifile in list_file:
        passn = ifile[-6:-3]
        cycle = ifile[-10:-8]
        fileg = '{}_p{:03d}.nc'.config(gpath, passn)
        data = rw.Sat_SKIM(ifile=ifile)
        grid = rw.Sat_SKIM(ifile=fileg)
        data.load_data(p, ur_model=[], ur_obs=[], instr=[], u_model=[], v_model=[],
                       lon=[], lat=[], time_nadir=time_nadir)
        grid.load_swath(p, radial_angle=[], angle=[], x_al=[], x_al_nadir=[],
                        x_ac=[])

        obs = type('', (), {})()
        obs = {}
        test = numpy.array(data.u_model)
        nil, nbeams=numpy.shape(test)
        sbeam_incid=numpy.zeros((nil, nbeams))
        obs['vrt'] = numpy.array(data.ur_model)
        obs['nsamp'], obs['nbeam'] = numpy.shape(obs['vrt'])
        obs['vrt'] = obs['vrt'].flatten()
        ind = numpy.where((obs['vrt'] > -1000))[0]
        if len(ind) > 1:
            obs = make_obs(data, grid, obs, ind)
            grd = make_grid(grid, obs)

        # OI
        grd = perform_oi(grd, obs)
        grd['vac'][numpy.abs(grd['ac2']) < 20] = numpy.nan

        grd['vx'] = (grd['vac'] * numpy.cos(grd['angle'])
                     + grd['val'] * numpy.cos(grd['angle'] + numpy.pi/2))
        grd['vy'] = (grd['vac'] * numpy.sin(grd['angle'])
                     + grd['val'] * numpy.sin(grd['angle'] + numpy.pi/2)
        outfile = os.path.join(datadir_output,
                               '{}_L2C_c{:02d}_p{:03d}.nc'.format(p.config, cycle,
                                                                  passn)
        write_l2(outfile, grd)
