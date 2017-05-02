import numpy
import math
from scipy import interpolate
import skimsimulator.mod_tools as mod_tools
import skimsimulator.const as const
import skimsimulator.rw_data as rw_data
import os
import logging

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
    npoints = 1
    # - Load SKIM orbit ground track
    logger.info('Load data from orbit file')
    bnorbit = os.path.basename(orbitfile)
    if ((bnorbit == 'orbs1a.txt') | (bnorbit == 'orbs1a_test.txt')
         | (bnorbit == 'orbits1_ifremer')
         | (bnorbit == 'orbits1_ifremer_test')):
        volon, volat, votime = numpy.loadtxt(orbitfile, usecols=(0, 1, 2),
                                         unpack=True)
    else:
        volon, volat, votime = numpy.loadtxt(orbitfile, usecols=(1, 2, 0),
                                         unpack=True)
    # - If orbit is at low resolution, interpolate at cycle (s) resolution
    cycle = const.cycle
    votime = votime * const.secinday
    ## macrocycle is fixed or depend on number of beam
    if numpy.mean(votime[1:] - votime[:-1]) != cycle:
        x, y, z = mod_tools.spher2cart(volon, volat)
        time_hr = numpy.arange(0., votime[-1], cycle)
        f = interpolate.interp1d(votime, x)
        x_hr = f(time_hr)
        f = interpolate.interp1d(votime, y)
        y_hr = f(time_hr)
        f = interpolate.interp1d(votime, z)
        z_hr = f(time_hr)
        lon_hr = numpy.zeros(len(x_hr)) + numpy.nan
        lat_hr = numpy.zeros(len(x_hr)) + numpy.nan
        for ii in range(len(x_hr)):
            lon_hr[ii], lat_hr[ii] = mod_tools.cart2spher(x_hr[ii], y_hr[ii],
                                                          z_hr[ii])
        # Cut orbit if too long
        ind = numpy.where((time_hr < const.satcycle * const.secinday))
        volon = lon_hr[ind]
        volat = lat_hr[ind]
        votime = time_hr[ind]

    # - Get number of points in orbit
    nop = numpy.shape(votime)[0]

    # - Get cycle period.
    tcycle = votime[nop-1] + votime[1] - votime[0]

    # shift time if the user needs to shift the time of the orbit
    try:
        pshift_time = p.shift_time
        if pshift_time:
            shift_index = numpy.where(votime >= pshift_time)[0]
            volon = numpy.hstack([volon[shift_index[0][0]:],
                                 volon[:shift_index[0][0]]])
            volat = numpy.hstack([volat[shift_index[0][0]:],
                                 volat[:shift_index[0][0]]])
    except:
        p.shift_time = None
    # shift lon if the user needs to shift the localisation of the orbit
    try:
        pshift_lon = p.shift_lon
        if pshift_lon is not None:
            volon = volon + pshift_lon
    except:
        p.shift_lon = None
    volon = (volon + 360) % 360

    # - Rearrange orbit starting from pass 1
    # Detect the beginning of pass 1 in orbit txt file. By definition, it is
    # the first passage at southernmost latitude.
    dlat = numpy.roll(volat, 1) - volat
    ind = numpy.where((dlat < 0) & (numpy.roll(dlat, 1) >= 0))
    # Shift coordinates, so that the first point of the orbit is the beginning
    # of pass 1
    decal = ind[0][-1]
    timeshift = votime[-1] - votime[decal]
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


    ## - Extract region
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

    # - Interpolate orbit at delta_al km resolution (default is delta_al=1)
    # Detect gap in time in stime (to detect step in x_al, lon and lat)
    dstime = stime_lr[:] - numpy.roll(stime_lr[:], 1)
    ind = numpy.where(dstime > 3*(votime[1] - votime[0]))
    index = numpy.hstack([0, ind[0]])
    nindex = numpy.shape(index)[0]
    lon = lon_lr
    lat = lat_lr
    stime = stime_lr
    x_al = x_al_lr
    lon = lon % 360
    # Save orbit data in Sat_SKIM object
    orb = rw_data.Sat_SKIM(file='{}.nc'.format(os.path.basename(orbitfile)))
    orb.x_al = x_al
    orb.time = stime
    orb.lon = lon
    orb.lat = lat
    orb.cycle = tcycle
    orb.al_cycle = distance[-1]
    orb.passtime = numpy.sort(passtime)
    orb.timeshift = timeshift
    return orb


def orbit2swath(modelbox, p, orb):
    # - Load altimeter orbit
    npoints = 1
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
    strpass = []
    # Compute rotating beams
    omega = p.rotation_speed * 2 * math.pi  / 60.
    # Number of beam to lighten:
    nbeam = len(p.list_shift) + 1
    # Loop on all passes after the first pass detected (note that ipass is
    # actually stored as ipass + 1 to have the first pass at 1 and ascending
    for ipass in range(ipass0, numpy.shape(passtime)[0]):
        # Detect indices corresponding to the pass
        if ipass == numpy.shape(passtime)[0]-1:
            ind = numpy.where((stime >= passtime[ipass]))[0]
        else:
            ind = numpy.where((stime >= passtime[ipass])
                              & (stime < passtime[ipass+1]))[0]
        nind = numpy.shape(ind)[0]
        # Compute swath grid if pass is in the subdomain
        if nind > 5:
            mod_tools.update_progress(float(ipass+1)
                                      / float(numpy.shape(passtime)[0]),
                                      'selected pass: ' + str(ipass+1), None)
            # Initialize SKIM grid, grid variables
            filesgrid = p.filesgrid + '_p' + str(ipass+1).zfill(3) + '.nc'
            sgrid = rw_data.Sat_SKIM(file=filesgrid)
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
            n_nad = numpy.shape(lon_beam[0])[0]
            angle_beam = [numpy.zeros(numpy.shape(lon_beam[0]))]
            radial_angle_tot = [numpy.zeros(numpy.shape(lon_beam[0]))]
            xal_beam = [x_al_nad]
            xac_beam = [x_al_nad * 0]
            x_al_tmp = x_al_nad - x_al_nad[0]
            inclination_angle = numpy.zeros(numpy.shape(lonnad))
            inclination_angle[1:] = numpy.arctan((latnad[1:] - latnad[:-1])
                                           / numpy.cos(latnad[1:] * math.pi / 180.)
                                           /(lonnad[1:] - lonnad[:-1]))
            inclination_angle[0] = inclination_angle[1]
            # Check that list of position, shift and angle have the same
            # dimension
            if len(p.list_pos) != len(p.list_shift) or \
                  len(p.list_pos) != len(p.list_angle) or \
                  len(p.list_angle) != len(p.list_shift):
                logger.error('Wrong length in list_pos, list_shift'\
                             'or list_angle')
                sys.exit(1)
            # Loop on beam to construct cycloid
            for angle, shift, beam in zip(p.list_pos, p.list_shift,
                                          p.list_angle):
                # Angle projected on the earth
                rc = const.Rearth* (beam * math.pi / 180 - numpy.arcsin(const.Rearth
                           *numpy.sin(math.pi- beam * math.pi / 180)
                           /(const.Rearth + const.sat_elev))) * 10**(-3)
                timebeamshift = timenad[shift::nbeam]
                beam_angle = omega * timebeamshift + angle
                xal = -( rc * numpy.sin(beam_angle))  /const.deg2km
                xac = (rc * numpy.cos(beam_angle))  / const.deg2km
                # Even pass: descending
                if ((ipass + 1) % 2 ==0):
                    #inclination = -inclination_angle[shift::nbeam] + math.pi
                    inclination = inclination_angle[shift::nbeam]
                    inclination_save = inclination_angle[0::nbeam]
                    radial_angle = -beam_angle + inclination - math.pi/2.
                # Odd pass: ascending
                else:
                    inclination = math.pi + inclination_angle[shift::nbeam]
                    #inclination = + inclination_angle[shift::nbeam]
                   # lon_tmp = (lonnad[shift::nbeam] + (xal * numpy.cos(inclination)
                   #                 + xac * numpy.sin(inclination))
                   #                 / numpy.cos(latnad[shift::nbeam]
                   #                 * math.pi / 180.))
                   # lat_tmp = (latnad[shift::nbeam] + (xal * numpy.sin(inclination)
                   #                 - xac * numpy.cos(inclination)))
                    inclination_save = math.pi + inclination_angle[0::nbeam]
                    radial_angle = -beam_angle + inclination - math.pi / 2.
                lon_tmp = (lonnad[shift::nbeam]
                                + (xal * numpy.cos(inclination)
                                + xac * numpy.sin(inclination))
                                / numpy.cos(latnad[shift::nbeam]
                                * math.pi / 180.))
                lat_tmp = (latnad[shift::nbeam]
                                + (xal * numpy.sin(inclination)
                                - xac * numpy.cos(inclination)))
                lon_tmp = (lon_tmp + 360) % 360
                # Concatenate list for each beam angle
                lon_beam.append(lon_tmp)
                lat_beam.append(lat_tmp)
                xal_beam.append(xal * const.deg2km)
                xac_beam.append(xac * const.deg2km)
                time_beam.append(timebeamshift)
                angle_beam.append(beam_angle)
                radial_angle_tot.append(radial_angle)
            # Save Sgrid object
            sgrid.timeshift = orb.timeshift
            sgrid.lon = lon_beam
            sgrid.lat = lat_beam
            sgrid.time = time_beam
            sgrid.x_al = xal_beam
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
    mod_tools.update_progress(1,  'All swaths have been processed', ' ')
    return None
