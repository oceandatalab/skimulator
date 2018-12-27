"""
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
"""
import os
import sys
import numpy
import math
from scipy import interpolate
import skimulator.mod_tools as mod_tools
import skimulator.const as const
import skimulator.rw_data as rw_data
import skimulator.mod_parallel
import time
import logging
import traceback

# Define logger level for debug purposes
logger = logging.getLogger(__name__)


def makeorbit(modelbox, p, orbitfile='orbit_292.txt', filealtimeter=None):
    '''Computes the swath of SKIM satellites on a subdomain.
    The path of the satellite is given by the orbit file and the subdomain
    corresponds to the one in the model. Note that a subdomain can be manually
    added in the parameters file. \n
    Inputs are satellite orbit (p.filesat), subdomain (modelbox), List of
    postion of six degree beams, list of position of twelve degree beams,
    rotation speed. \n
    Outputs are netcdf files containing SKIM position (lon, lat number of days
    in a cycle, distance crossed in a cycle, time)
    '''
    # - Load SKIM orbit ground track
    logger.info('Load data from orbit file')
    if p.order_orbit_col is None:
        volon, volat, votime = numpy.loadtxt(orbitfile, usecols=(0, 1, 2),
                                             comments='#', unpack=True)
    else:
        ncols = p.order_orbit_col
        volon, volat, votime = numpy.loadtxt(orbitfile, usecols=ncols,
                                             comments='#', unpack=True)
    dic_sat = {}
    with open(orbitfile, 'r') as fh:
        for i, line in enumerate(fh):
            if line.strip().startswith('#'):
                key, value = line.strip().split('=')
                dic_sat[key[1:].strip()] = float(value.strip())
            else:
                break
    if 'cycle' in dic_sat.keys() and 'elevation' in dic_sat.keys():
        p.satcycle = dic_sat['cycle']
        p.sat_elev = dic_sat['elevation']

    # - If orbit is at low resolution, interpolate at cycle (s) resolution
    cycle = p.cycle
    votime = votime * const.secinday
    # macrocycle is fixed or depend on number of beam
    logger.info('Interpolate orbit at {} seconds'.format(cycle))
    # votime2 = votime + 0
    # - If orbit is at low resolution, interpolate at 0.5 s resolution
    hr_step = 0.5
    if numpy.mean(votime[1:] - votime[:-1]) > hr_step:
        x, y, z = mod_tools.spher2cart(volon, volat)
        time_hr = numpy.arange(0., votime[-1], hr_step)
        # Parallelisation of this loop?
        # jobs = []
        # jobs.append([votime, time_hr, x])
        # jobs.append([votime, time_hr, y])
        # jobs.append([votime, time_hr, z])
        # worker
        f = interpolate.interp1d(votime, x)
        x_hr = f(time_hr)
        f = interpolate.interp1d(votime, y)
        y_hr = f(time_hr)
        f = interpolate.interp1d(votime, z)
        z_hr = f(time_hr)
        lon_hr = numpy.zeros(len(x_hr)) + numpy.nan
        lat_hr = numpy.zeros(len(x_hr)) + numpy.nan
        lon_hr, lat_hr = mod_tools.cart2sphervect(x_hr, y_hr, z_hr)
        # Cut orbit if more than an orbit cycle
        if p.satcycle is None:
            p.satcycle = const.satcycle
        ind = numpy.where((time_hr < p.satcycle * const.secinday))
        volon = lon_hr[ind]
        volat = lat_hr[ind]
        votime = time_hr[ind]

    # - Get number of points in orbit
    nop = numpy.shape(votime)[0]

    # - Get cycle period.
    tcycle = votime[nop-1] + votime[1] - votime[0]

    # shift time if the user needs to shift the time of the orbit
    if p.shift_time is not None:
        shift_index = numpy.where(votime >= p.shift_time)[0]
        volon = numpy.hstack([volon[shift_index[0][0]:],
                             volon[:shift_index[0][0]]])
        volat = numpy.hstack([volat[shift_index[0][0]:],
                             volat[:shift_index[0][0]]])
    # shift lon if the user needs to shift the localisation of the orbit
    if p.shift_lon is not None:
        volon = volon + p.shift_lon
    volon = (volon + 360) % 360

    # - Rearrange orbit starting from pass 1
    # Detect the beginning of pass 1 in orbit txt file. By definition, it is
    # the first passage at southernmost latitude.
    dlat = numpy.roll(volat, 1) - volat
    ind = numpy.where((dlat < 0) & (numpy.roll(dlat, 1) >= 0))
    # Shift coordinates, so that the first point of the orbit is the beginning
    # of pass 1
    decal = ind[0][-1]
    # timeshift = votime[-1] - votime[decal]
    volon = numpy.hstack([volon[decal:], volon[:decal]])
    volat = numpy.hstack([volat[decal:], volat[:decal]])
    votime = numpy.hstack([votime[decal:], votime[:decal]])
    votime = (votime - votime[0]) % tcycle
    if votime[numpy.where(votime < 0)]:
        logger.warn('WARNING: there are negative times in your orbit')
    del ind
    # Compute the initial time of each pass
    dlat = numpy.roll(volat, 1) - volat
    ind = numpy.where(((dlat < 0) & (numpy.roll(dlat, 1) >= 0)) | ((dlat > 0)
                      & (numpy.roll(dlat, 1) <= 0)))
    # index=numpy.hstack([0,ind[0]-1])
    index = ind[0]
    passtime = votime[index]  # Times of pass

    # - Extract region
    matnpbox = numpy.zeros((nop))
    if modelbox[0] > modelbox[1]:
        matnpbox[numpy.where((((modelbox[0] - 1) <= volon)
                 | (volon <= (modelbox[1] + 1)))
                 & ((modelbox[2] - 1) <= volat)
                 & ((modelbox[3] + 1) >= volat))] = 1
    else:
        matnpbox[numpy.where(((modelbox[0] - 1) <= volon)
                 & (volon <= (modelbox[1] + 1))
                 & ((modelbox[2] - 1) <= volat)
                 & ((modelbox[3] + 1) >= volat))] = 1
    norp = int(numpy.sum(matnpbox))
    # Initialize total distance travelled by the satellite since the first
    # point of the cycle in the subdomain at low (orbital file) resolution
    x_al_lr = numpy.zeros((norp))
    lon_lr = numpy.zeros((norp))
    lat_lr = numpy.zeros((norp))
    stime_lr = numpy.zeros((norp))

    # Initialize vector with accumulated distance travelled by the satellite
    indp = 0
    distance = numpy.zeros((nop))
    # Compute and store distances and coordinates that are in the defined
    # subdomain
    logger.info('Compute SKIM nadir coordinate in the new domain')
    for i in range(0, nop - 1):
        if p.progress_bar is True:
            mod_tools.update_progress(float(i) / float(nop-1), None, None)
        if abs(volon[i + 1] - volon[i]) > 1:
            if volon[i + 1] > 180.:
                volon[i + 1] = volon[i + 1] - 360
            if volon[i] > 180.:
                volon[i] = volon[i] - 360
        distance[i+1] = (distance[i] + numpy.sqrt(((volon[i+1]-volon[i])
                         * const.deg2km*numpy.cos(volat[i+1]
                         * 2*math.pi/360.))**2 + ((volat[i+1] - volat[i])
                         * const.deg2km)**2))  # numpy.sum(dl[:i])
        volon[i + 1] = (volon[i + 1] + 360) % 360
        if matnpbox[i]:
            x_al_lr[indp] = distance[i]
            lon_lr[indp] = (volon[i] + 360) % 360
            lat_lr[indp] = volat[i]
            stime_lr[indp] = votime[i]
            indp += 1

    # - Interpolate orbit at cycle resolution (default is cycle=0.0096)
    # Detect gap in time in stime (to detect step in x_al, lon and lat)
    dstime = stime_lr[:] - numpy.roll(stime_lr[:], 1)
    ind = numpy.where(dstime > 3*(votime[1] - votime[0]))
    index = numpy.hstack([0, ind[0]])
    nindex = numpy.shape(index)[0]
    # Initialize along track distance, time and coordinates at cycle
    # resolution
    if nindex > 1:
        dgap = numpy.zeros((nindex))
        for i in range(1, nindex):
            dgap[i] = stime_lr[index[i]] - stime_lr[max(index[i] - 1, 0)]
        Ninterp = int(((stime_lr[-1] - stime_lr[0] - sum(dgap))
                      / float(p.cycle)) + 1)
        x_al = numpy.zeros((Ninterp))
        stime = numpy.zeros((Ninterp))
        lon = numpy.zeros((Ninterp))
        lat = numpy.zeros((Ninterp))
        imin = 0
        imax = 0
        substract_point = 0
        for i in range(0, nindex - 1):
            imax = imin + int((stime_lr[index[i+1]-1] - stime_lr[index[i]])
                              / float(p.cycle)) + 1
            if imax <= (imin + 1):
                x_al[imin] = x_al_lr[index[i]]
                stime[imin] = stime_lr[index[i]]
                lon[imin] = lon_lr[index[i]]
                lat[imin] = lat_lr[index[i]]
            else:
                slicei = slice(index[i], index[i + 1])
                nofp = (stime_lr[index[i+1] - 1]-stime_lr[index[i]]) / p.cycle
                if (nofp).is_integer():
                    imax = imax - 1
                    substract_point += 1
                stime[imin: imax] = numpy.arange(stime_lr[index[i]],
                                                 stime_lr[index[i+1] - 1],
                                                 p.cycle)
                x_al[imin: imax] = numpy.interp(stime[imin: imax],
                                                stime_lr[slicei],
                                                x_al_lr[slicei])
                loncirc = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(
                                        lon_lr[slicei])))
                lon[imin: imax] = numpy.interp(stime[imin: imax],
                                               stime_lr[slicei],
                                               loncirc)
                lat[imin: imax] = numpy.interp(stime[imin: imax],
                                               stime_lr[slicei],
                                               lat_lr[slicei])
            imin = imax
        if substract_point > 0:
            substract_point = 1
        _tmp = numpy.arange(stime_lr[index[-1]], stime_lr[index[-1]]
                            + (Ninterp - imin)*p.cycle
                            - p.cycle * substract_point, p.cycle)
        if len(_tmp) > len(stime[imin:]):
            _tmp2 = _tmp[:-1]
        elif len(_tmp) < len(stime[imin:]):
            print(len(_tmp), len(stime[imin:]))
            _tmp2 = numpy.arange(stime_lr[index[-1]], stime_lr[index[-1]]
                                + (Ninterp - imin)*p.cycle, p.cycle)
            #logger.error('damnit, there is a bug, contact me')
        else:
            _tmp2 = _tmp

        stime[imin:] = _tmp2 #numpy.arange(stime_lr[index[-1]], stime_lr[index[-1]]
                              #      + (Ninterp - imin)*p.cycle
                               #     - p.cycle * substract_point, p.cycle)
        x_al[imin:] = numpy.interp(stime[imin:], stime_lr[index[-1]:],
                                   x_al_lr[index[-1]:])
        loncirc = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(
                                                          lon_lr[index[-1]:])))
        lon[imin:] = numpy.interp(stime[imin:], stime_lr[index[-1]:], loncirc)
        lat[imin:] = numpy.interp(stime[imin:], stime_lr[index[-1]:],
                                  lat_lr[index[-1]:])
    else:
        Ninterp = int(((stime_lr[-2] - stime_lr[0]) / float(p.cycle)) + 1)
        stime = numpy.arange(stime_lr[0], stime_lr[-2], p.cycle)
        x_al = numpy.interp(stime, stime_lr[:-1], x_al_lr[:-1])
        loncirc = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(lon_lr[:-1])))
        lon = numpy.interp(stime, stime_lr[:-1], loncirc)
        lat = numpy.interp(stime, stime_lr[:-1], lat_lr[:-1])
        #x_al = numpy.arange(x_al_lr[0], x_al_lr[-2], p.delta_al)
        #stime = numpy.interp(x_al, x_al_lr[:-1], stime_lr[:-1])
        #lon = numpy.interp(x_al, x_al_lr[:-1], loncirc)
        #lat = numpy.interp(x_al, x_al_lr[:-1], lat_lr[:-1])
    lon = lon % 360
    # Save orbit data in Sat_SKIM object
    orb = rw_data.Sat_SKIM(ifile='{}.nc'.format(os.path.basename(orbitfile)))
    orb.x_al = x_al
    orb.time = stime
    orb.lon = lon
    orb.lat = lat
    orb.cycle = tcycle
    orb.al_cycle = distance[-1]
    orb.passtime = numpy.sort(passtime)
    orb.timeshift = p.timeshift
    return orb


def orbit2swath(modelbox, p, orb, die_on_error):
    '''Computes the swath of SKIM satellites on a subdomain from an orbit.
    The path of the satellite is given by the orbit file and the subdomain
    corresponds to the one in the model. Note that a subdomain can be manually
    added in the parameters file. \n
    Inputs are satellite orbit (p.filesat), subdomain (modelbox), Swath
    parameters (half gap distance p.halfgap, half swath distance p.halfswath,
    along track
    resolution p.delta_al, across track resolution p.delta_ac). \n
    Outputs are netcdf files containing SKIM grid (along track distance x_al,
    radial angle, longitude lon and latitude lat,
    number of days in a cycle cycle, distance crossed in a cycle cycle_al,
    time'''
    ''' Compute orbit from Swath '''
    # - Load altimeter orbit
    x_al = orb.x_al
    stime = orb.time
    lon = orb.lon
    lat = orb.lat
    tcycle = orb.cycle
    al_cycle = orb.al_cycle
    passtime = orb.passtime
    # - Computation of SKIM grid and storage by passes
    logger.info('\n Compute SKIM grid')
    # Detect first pass that is in the subdomain
    ipass0 = 0
    # Check that list of position, shift and angle have the same
    # dimension
    if len(p.list_pos) != len(p.list_shift) or \
          len(p.list_pos) != len(p.list_angle) or \
          len(p.list_angle) != len(p.list_shift):
        logger.error('Wrong length in list_pos, list_shift'
                     'or list_angle')
        sys.exit(1)
    # Loop on all passes after the first pass detected (note that ipass is
    # actually stored as ipass + 1 to have the first pass at 1 and ascending
    jobs = []
    p2 = mod_tools.todict(p)
    for ipass in range(ipass0, numpy.shape(passtime)[0]):
        jobs.append([ipass, p2, passtime, stime, x_al, tcycle, al_cycle, lon,
                     lat, orb.timeshift])
    ok = make_skim_grid(p.proc_count, jobs, die_on_error, p.progress_bar)
    if p.progress_bar is True:
        mod_tools.update_progress(1,  'All swaths have been processed', ' ')
    else:
        logger.info('All swaths have been processed')
    return ok


def make_skim_grid(_proc_count, jobs, die_on_error, progress_bar):
    """ Compute SWOT grids for every pass in the domain"""
    # - Set up parallelisation parameters
    proc_count = min(len(jobs), _proc_count)

    status_updater = mod_tools.update_progress_multiproc
    jobs_manager = skimulator.mod_parallel.JobsManager(proc_count,
                                                       status_updater,
                                                       exc_formatter,
                                                       err_formatter)
    ok = jobs_manager.submit_jobs(worker_method_grid, jobs, die_on_error,
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


def err_formatter(pid, ipass, cycle, exc):
    """Transform errors stored by the JobsManager into readable messages."""
    msg = None
    if cycle < 0:
        msg = '/!\ Error occurred while processing pass {}'.format(ipass)
    else:
        _msg = '/!\ Error occurred while processing cycle {} on pass {}'
        msg = _msg.format(cycle, ipass)
    return msg


def worker_method_grid(*args, **kwargs):
    msg_queue, ipass, p2 = args[:3]
    passtime, stime, x_al, tcycle, al_cycle, lon, lat, timeshift = args[3:]

    p = mod_tools.fromdict(p2)
    # Compute rotating beams
    omega = p.rotation_speed * 2 * math.pi / 60.
    # Number of beam to lighten:
    nbeam = len(p.list_shift) + 1
    sat_elev= p.sat_elev
    if sat_elev is None:
        sat_elev = const.sat_elev
    # Detect indices corresponding to the pass
    if ipass == numpy.shape(passtime)[0]-1:
        ind = numpy.where((stime >= passtime[ipass]))[0]
    else:
        ind = numpy.where((stime >= passtime[ipass])
                          & (stime < passtime[ipass+1]))[0]
    nind = numpy.shape(ind)[0]
    # Compute swath grid if pass is in the subdomain
    if nind > 5:
        # Initialize SKIM grid, grid variables
        filesgrid = '{}_p{:03d}.nc'.format(p.filesgrid, ipass + 1)
        sgrid = rw_data.Sat_SKIM(ifile=filesgrid)
        sgrid.x_al = x_al[ind]
        x_al_nad = x_al[ind]
        x_al_nad = x_al_nad[0::nbeam]
        sgrid.cycle = tcycle
        sgrid.al_cycle = al_cycle
        sgrid.ipass = ipass + 1

        # Compute nadir coordinate and initialize angles
        lonnad = (lon[ind] + 360) % 360
        latnad = + lat[ind]
        timenad = + stime[ind]
        lon_beam = [lonnad[0::nbeam]]
        lat_beam = [latnad[0::nbeam]]
        time_beam = [timenad[0::nbeam]]
        # n_nad = numpy.shape(lon_beam[0])[0]
        angle_beam = [numpy.zeros(numpy.shape(lon_beam[0]))]
        radial_angle_tot = [numpy.zeros(numpy.shape(lon_beam[0]))]
        xal_beam = [x_al_nad]
        xac_beam = [x_al_nad * 0]
        xal_beam_tot = [x_al_nad]
        # x_al_tmp = x_al_nad - x_al_nad[0]
        inclination_angle = numpy.zeros(numpy.shape(lonnad))
        inclination_angle[1:] = numpy.arctan((latnad[1:] - latnad[:-1])
                                             / numpy.cos(latnad[1:]
                                             * math.pi / 180.)
                                             / (lonnad[1:] - lonnad[:-1]))
        inclination_angle[0] = inclination_angle[1]
        # Loop on beam to construct cycloid
        for angle, shift, beam in zip(p.list_pos, p.list_shift,
                                      p.list_angle):
            # Angle projected on the earth
            rc = (const.Rearth * (beam * math.pi/180
                  - numpy.arcsin(const.Rearth * numpy.sin(math.pi - beam
                  * math.pi/180) / (const.Rearth + sat_elev)))
                  * 10**(-3))
            timebeamshift = timenad[shift::nbeam]
            beam_angle = omega * timebeamshift + angle
            xal = -(rc * numpy.sin(beam_angle)) / const.deg2km
            xac = (rc * numpy.cos(beam_angle)) / const.deg2km
            sign_ac = 1
            # Even pass: descending
            if ((ipass + 1) % 2 == 0):
                # inclination = -inclination_angle[shift::nbeam] + math.pi
                inclination = inclination_angle[shift::nbeam]
                inclination_save = inclination_angle[0::nbeam]
                radial_angle = -beam_angle + inclination - math.pi/2.
                sign_ac = -1
            # Odd pass: ascending
            else:
                inclination = math.pi + inclination_angle[shift::nbeam]
                # inclination = + inclination_angle[shift::nbeam]
                inclination_save = math.pi + inclination_angle[0::nbeam]
                radial_angle = -beam_angle + inclination - math.pi / 2.
            lon_tmp = (lonnad[shift::nbeam] + (xal * numpy.cos(inclination)
                       + xac * numpy.sin(inclination))
                       / numpy.cos(latnad[shift::nbeam] * math.pi / 180.))
            lat_tmp = (latnad[shift::nbeam] + (xal * numpy.sin(inclination)
                       - xac * numpy.cos(inclination)))
            lon_tmp = (lon_tmp + 360) % 360
            # Concatenate list for each beam angle
            lon_beam.append(lon_tmp)
            lat_beam.append(lat_tmp)
            xal_beam.append(xal * const.deg2km * sign_ac)
            xac_beam.append(xac * const.deg2km * sign_ac)
            xal_beam_tot.append(sgrid.x_al[shift::nbeam])
            time_beam.append(timebeamshift)
            angle_beam.append(beam_angle)
            radial_angle_tot.append(radial_angle)
        # Save Sgrid object
        sgrid.timeshift = timeshift
        sgrid.lon = lon_beam
        sgrid.lat = lat_beam
        sgrid.time = time_beam
        sgrid.x_al = xal_beam
        sgrid.x_al_tot = xal_beam_tot
        sgrid.x_ac = xac_beam
        sgrid.list_angle = p.list_angle
        sgrid.list_pos = p.list_pos
        sgrid.beam_angle = angle_beam
        sgrid.radial_angle = radial_angle_tot
        sgrid.angle = inclination_save
        # Remove grid file if it exists and save it
        if os.path.exists(filesgrid):
            os.remove(filesgrid)
        sgrid.write_swath(p)
    msg_queue.put((os.getpid(), ipass, None, None))
    return None
