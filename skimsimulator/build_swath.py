import numpy
import math
from scipy import interpolate
import skimsimulator.mod_tools as mod_tools
import skimsimulator.const as const
import skimsimulator.rw_data as rw_data
import os


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
    print('Load data from orbit file')
    bnorbit = os.path.basename(orbitfile)
    if ((bnorbit == 'orbs1a.txt') | (bnorbit == 'orbs1a_test.txt')
         | (bnorbit == 'orbits1_ifremer')):
        volon, volat, votime = numpy.loadtxt(orbitfile, usecols=(0, 1, 2),
                                         unpack=True)
    else:
        volon, volat, votime = numpy.loadtxt(orbitfile, usecols=(1, 2, 0),
                                         unpack=True)
    # - If orbit is at low resolution, interpolate at 0.5 s resolution
    cycle = 0.0368
    votime = votime * 86400
    ## macrocycle is fixed or depend on number of beam
    macrocycle = cycle * 6 #(len(p.list_pos_12) + len(p.list_pos_6) + 1)
    if numpy.mean(votime[1:] - votime[:-1]) != macrocycle:
        x, y, z = mod_tools.spher2cart(volon, volat)
        time_hr = numpy.arange(0., votime[-1], macrocycle)
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
#        time_hr = time_hr / 86400.
        ind = numpy.where((time_hr < 12.*86400))
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
        print('WARNING: there are negative times in your orbit')
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
    print('Compute SKIM nadir coordinate in the new domain')
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
    # Initialize along track distance, time and coordinates at delta_al
    # resolution
    '''
    if nindex > 1:
        dgap = numpy.zeros((nindex))
        for i in range(1, nindex):
            dgap[i] = x_al_lr[index[i]] - x_al_lr[max(index[i] - 1, 0)]
        Ninterp = (int((x_al_lr[-1] - x_al_lr[0] - sum(dgap))
                   / float(p.delta_al)) + 1)
        x_al = numpy.zeros((Ninterp))
        stime = numpy.zeros((Ninterp))
        lon = numpy.zeros((Ninterp))
        lat = numpy.zeros((Ninterp))
        imin = 0
        imax = 0
        for i in range(0, nindex - 1):
            imax = imin + int((x_al_lr[index[i+1]-1] - x_al_lr[index[i]])
                              / float(p.delta_al)) + 1
            if imax <= (imin + 1):
                x_al[imin] = x_al_lr[index[i]]
                stime[imin] = stime_lr[index[i]]
                lon[imin] = lon_lr[index[i]]
                lat[imin] = lat_lr[index[i]]
            else:
                x_al[imin: imax] = x_al_lr[index[i]], index[i+1]]
                stime[imin: imax] = stime_lr[index[i]: index[i+1]])
                loncirc = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(
                                        lon_lr[index[i]: index[i+1]])))
                # if numpy.min(lon_lr[index[i]:index[i+1]])<1.
                # and numpy.max(lon_lr[index[i]:index[i+1]])>359.:
                # lontmp=lon_lr[index[i]:index[i+1]]
                # lontmp[numpy.where(lontmp>180.)]=lontmp[numpy.where(
                # lontmp>180.)]-360.
                # lon[imin:imax]=numpy.interp(x_al[imin:imax],
                # x_al_lr[index[i]:index[i+1]], lontmp)
                #    lon[imin:imax]=(lon[imin:imax]+360)%360
                # else:
                #    lon[imin:imax]=numpy.interp(x_al[imin:imax],
                # x_al_lr[index[i]:index[i+1]], lon_lr[index[i]:index[i+1]])
                lon = lon_lr[index[i]: index[i+1]]
                lat[imin: imax] = lat_lr[index[i]: index[i+1]])
            imin = imax
        x_al[imin:] = x_al_lr[index[-1]:]
        stime[imin:] = stime_lr[index[-1]:]
        loncirc = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(
                                lon_lr[index[-1]:])))
        # if numpy.min(lon_lr[index[-1]:])<1.
        # and numpy.max(lon_lr[index[-1]:])>539:
        # lontmp=lon_lr[index[-1]:]
        # lontmp[numpy.where(lontmp>180.)]=lontmp[numpy.where(lontmp>180.)]-360
        # lon[imin:]=numpy.interp(x_al[imin:], x_al_lr[index[-1]:], lontmp)
        #    lon[imin:]=(lon[imin:]+360)%360
        # else:
        #    lon[imin:]=numpy.interp(x_al[imin:], x_al_lr[index[-1]:], lon_lr[index[-1]:])
        lon[imin:] = lon_lr[index[-1]:]
        lat[imin:] = lat_lr[index[-1]:]
    else:
        x_al = x_al_lr[:-1]
        stime = stime_lr[:-1]
        loncirc = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(lon_lr[:-1])))
        lon = lon_lr[:-1]
        lat = lat_lr[:-1]
    '''
    lon = lon_lr
    lat = lat_lr
    stime = stime_lr
    x_al = x_al_lr
    lon = lon % 360
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
    print('\n Compute SKIM grid')
    # Detect first pass that is in the subdomain
    ipass0 = 57
    strpass = []
    # Compute rotating beams
    omega = p.rotation_speed * 2 * math.pi  / 60.
    #omega = p.rotation_speed / 60.
    rc_12 = 140
    rc_6 = 70
    #list_pos_12 = [0, math.pi / 2., math.pi, math.pi * 3. / 2.)
    #list_angle_6 = [math.pi / 4.]
    #list_pos_12 = p.list_pos_12 # [0, math.pi / 2., math.pi, math.pi * 3. / 2.)
    #list_pos_6 = p.list_pos_6 # [math.pi / 4.]
    # Loop on all passes after the first pass detected
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
            sgrid.cycle = tcycle
            sgrid.al_cycle = al_cycle

            # Project in cartesian coordinates satellite ground location
            lonnad = (lon[ind] + 360) % 360
            latnad = lat[ind]
            lon_beam = [lonnad]
            lat_beam = [latnad]
            time_beam = [stime[ind]]
            x_al_tmp = sgrid.x_al - sgrid.x_al[0]
            sgrid.angle = lonnad * 0
            sgrid.angle[1:] = numpy.arctan((lonnad[1:] - lonnad[:-1])
                                        /(latnad[1:] - latnad[:-1]))
            sgrid.angle[0] = sgrid.angle[1]
            inclination = p.inclination
            # inclination = sgrid.angle + math.pi/2
            x_nad, y_nad, z_nad = mod_tools.spher2cart(lonnad, latnad)
            #import pdb ; pdb.set_trace()
            for angle_12, shift_12 in zip(p.list_pos_12, p.list_shift_12):
                timebeamshift = stime[ind] + shift_12
                #lon_tmp = (lonnad + (rc_12 * numpy.cos(omega * timebeamshift
                #            + angle_12))/const.deg2km/numpy.cos(lat[ind] * math.pi / 180.))
                #lat_tmp = (latnad + (rc_12 * numpy.sin(omega * timebeamshift
                #           + angle_12))/const.deg2km)
                if (ipass % 2 ==0):
                    xal = ( rc_12 * numpy.cos(omega * timebeamshift
                          + angle_12))  /const.deg2km
                    xac = (rc_12 * numpy.sin(omega * timebeamshift
                           + angle_12))  / const.deg2km
                    lon_tmp = (lonnad - (xal * numpy.cos(inclination)
                                    + xac * numpy.sin(inclination))
                                    / numpy.cos(latnad * math.pi / 180.))
                    lat_tmp = (latnad + (xal * numpy.sin(inclination)
                                    - xac * numpy.cos(inclination)))
                else:
                    xal = ( rc_12 * numpy.sin(omega * timebeamshift
                          + angle_12))  /const.deg2km
                    xac = (rc_12 * numpy.cos(omega * timebeamshift
                           + angle_12))  / const.deg2km
                    lon_tmp = (lonnad - (xal * numpy.cos(inclination)
                                    + xac * numpy.sin(inclination))
                                    / numpy.cos(latnad * math.pi / 180.))
                    lat_tmp = (latnad - (xal * numpy.sin(inclination)
                                    - xac * numpy.cos(inclination)))
                #lon_tmp, lat_tmp = mod_tools.cart2sphervect(x, y, z_nad)
                lon_tmp = (lon_tmp + 360) % 360
                lon_beam.append(lon_tmp)
                lat_beam.append(lat_tmp)
                time_beam.append(timebeamshift)
            for angle_6, shift_6 in zip(p.list_pos_6, p.list_shift_6):
                timebeamshift = stime[ind] + shift_6
                # x = (x_nad + rc_6 * numpy.cos(omega * timebeamshift
                    #+ angle_6))
                # y = (y_nad + rc_6 * numpy.sin(omega * timebeamshift
                    #+ angle_6))
                if (ipass % 2 ==0):
                    xal = ( rc_6 * numpy.cos(omega * timebeamshift
                          + angle_6)) /const.deg2km
                    xac = (rc_6 * numpy.sin(omega * timebeamshift
                           + angle_6)) / const.deg2km
                    lon_tmp = (lonnad - (xal * numpy.cos(inclination)
                                    + xac * numpy.sin(inclination))
                                    / numpy.cos(latnad * math.pi / 180.))
                    lat_tmp = (latnad + (xal * numpy.sin(inclination)
                                    - xac * numpy.cos(inclination)))
                else:
                    xal = ( rc_6 * numpy.sin(omega * timebeamshift
                          + angle_6)) /const.deg2km
                    xac = (rc_6 * numpy.cos(omega * timebeamshift
                           + angle_6)) / const.deg2km
                    lon_tmp = (lonnad - (xal * numpy.cos(inclination)
                                    + xac * numpy.sin(inclination))
                                    / numpy.cos(latnad * math.pi / 180.))
                    lat_tmp = (latnad - (xal * numpy.sin(inclination)
                                    - xac * numpy.cos(inclination)))
                # lon_tmp, lat_tmp = mod_tools.cart2sphervect(x, y, z_nad)
                lon_tmp = (lon_tmp + 360) % 360
                lon_beam.append(lon_tmp)
                lat_beam.append(lat_tmp)
                time_beam.append(timebeamshift)
            sgrid.timeshift = orb.timeshift
            sgrid.lon = lon_beam
            sgrid.lat = lat_beam
            sgrid.time = time_beam
            if os.path.exists(filesgrid):
                os.remove(filesgrid)
            sgrid.write_swath(p)
    mod_tools.update_progress(1,  'All swaths have been processed', ' ')
    return None
