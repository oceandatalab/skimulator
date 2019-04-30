import numpy
import math
import logging
import skimulator.mod_tools as mod_tools


# Define logger level for debug purposes
logger = logging.getLogger(__name__)


def compute_erruwb(p, sgrid, ouput_var):
    ''' Compute the reconstructed error of uwb using information from the uwb
    calculated with closed beams. '''
    # Compute number of azimuth per mega cycle
    naz = (60 / (p.rotation_speed * p.cycle * len(p.list_pos)))
    # Initialize errdcos
    # Convert list into arrays
    lon_array = numpy.transpose(numpy.array(sgrid.lon[1:][:]))
    lat_array = numpy.transpose(numpy.array(sgrid.lat[1:][:]))
    uwb_r_array = numpy.transpose(numpy.array(ouput_var['uwb'][1:][:]))
    mask_array = (numpy.isnan(uwb_r_array))
    err_uwb = []
    err_uwb.append(numpy.full(numpy.shape(sgrid.lat[0][:]), numpy.nan))
    # mask_array = numpy.transpose(uwb_r_array.mask)
    for i in range(1, len(p.list_pos) + 1):
        mask_minicycle = mask_array[:, i - 1]
        radial_angle = sgrid.radial_angle[:, i - 1]
        N = len(mask_minicycle)
        dtheta = numpy.mean(abs(radial_angle[1:] - radial_angle[:-1]))
        # errdcos = numpy.zeros(numpy.shape(radial_angle)) * numpy.nan
        erruwb = numpy.full(numpy.shape(radial_angle), numpy.nan)
        for ib in range(N):
            if mask_minicycle[ib] is True:
                continue
            theta = radial_angle[ib] % (2 * math.pi)
            erruwb[ib] = 0
            ntheta2 = 0
            for theta2 in numpy.arange(theta - math.pi/2,
                                       theta + math.pi/2, dtheta):
                theta2rad = theta2 % (2 * math.pi)
                start_az = int(max(0, (ib - naz*3)))
                end_az = int(min(N, (ib + naz*3)))
                slice_az = slice(start_az, end_az)

                lon = lon_array[slice_az, :]
                lat = lat_array[slice_az, :]
                uwb_r_loc = uwb_r_array[slice_az, :]
                mask_ind = mask_array[slice_az, :]
                lon[numpy.where(mask_ind is True)] = -1.36*10**9
                lat[numpy.where(mask_ind is True)] = -1.36*10**9
                angle = sgrid.radial_angle[slice_az, :]
                angle = numpy.mod(angle, 2 * math.pi)
                ind_angle = numpy.where((angle >= (theta2rad - dtheta))
                                        & (angle < (theta2rad + dtheta)))
                # To be tested: ind_angle can be empty near coast?
                if len(ind_angle) == 0:
                    logger.debug('no ind_angle found, pass is to small?')
                    continue
                lon = lon[ind_angle]
                lat = lat[ind_angle]
                if lon.size == 0:
                    continue
                uwb_r_loc2 = uwb_r_loc[ind_angle]
                sgridlon = sgrid.lon[i][ib]
                if (numpy.max(lon)) > 359 and (numpy.min(lon) < 1):
                    lon[lon > 180] = lon[lon > 180] - 360
                    lon = numpy.rad2deg(numpy.unwrap(numpy.deg2rad(lon)))
                    if sgridlon > 180:
                        sgridlon = sgridlon - 360
                        sgridlon = numpy.rad2deg(numpy.deg2rad(sgridlon))
                dlon_km = ((lon - sgridlon)*111
                           * numpy.cos(sgrid.lat[i][ib] * math.pi / 180.))
                dlat_km = (lat - sgrid.lat[i][ib]) * 111
                dist = numpy.sqrt(dlon_km**2 + dlat_km**2)
                if len(dist) > 0:
                    ind_dist = numpy.argmin(dist)
                    erruwb[ib] += (numpy.cos(theta
                                    - angle[ind_angle[0][ind_dist],
                                            ind_angle[1][ind_dist]])
                                    * uwb_r_loc2[ind_dist] * dtheta
                                    / (math.pi) * 2)
                else:
                    erruwb[ib] = numpy.nan
                ntheta2 += 1
        err_uwb.append(erruwb - uwb_r_array[:, i-1])
    return err_uwb


def combine_usr(lon, lat, usr, dazi, angle, wnd_dir):
    '''
    PURPOSE : recombine measurements of various azimuths and ambiguity
    removal using wind direction
    Provides a single estimate (not applicable to a set of points)
    INPUT : - usrazi : radial Stokes drift for each azimuth in azibin [m/s], len = nazi
    dazi : azimuth increment [degrees]
    azimuthat : along track azimuth [radians]
    thwnd : direction of the wind [radians], without ambiguity. Must be in the [0;2*pi[ interval
    OUTPUT : usrhat : estimate of the radial stokes drift [m/s]
    '''

    dazi = 33
    azibin = numpy.arange(0, 180, dazi)
    # Ambiguity from SKIM measurment
    usrabs = numpy.abs(usr)
    # From rad to deg and considering pi ambiguity
    degangle = numpy.rad2deg(numpy.mod(angle, numpy.pi))
    nsample, nbeam = numpy.shape(usrabs)
    usr_comb = numpy.full((nsample, nbeam), numpy.nan)
    dmax = 300
    for ibeam in range(nbeam):
        for isample in range(nsample):
            if not numpy.isfinite(usrabs[isample, ibeam]):
                continue
            plon = lon[isample, ibeam]
            plat = lat[isample, ibeam]
            pangle = degangle[isample, ibeam]
            rangle = angle[isample, ibeam]
            dist = mod_tools.dist_sphere(plon, lon, plat, lat)
            Idist0 = numpy.where(dist < dmax)
            #azimuthr = azimuth[Idist0].ravel()
            lonr = lon[Idist0].ravel()
            latr = lat[Idist0].ravel()
            usrr = usrabs[Idist0].ravel()
            distr = dist[Idist0].ravel()
            ang_azi = numpy.digitize(degangle[Idist0].ravel(), azibin) - 1
            Ifinite = numpy.isfinite(usrr)
            # for each bin, find the closest point and append radial Stokes drift with ambiguity on sign
            usrazi = numpy.zeros(len(azibin))
            for iazi in range(len(azibin)):
                I1 = numpy.where(numpy.logical_and(ang_azi==iazi, Ifinite))
                distazi = mod_tools.dist_sphere(plon, lonr[I1], plat, latr[I1])
                try:
                    Idist = numpy.argmin(distazi)
                    usrazi[iazi] = usrr[I1][Idist]
                except:
                    usrazi[iazi] = numpy.nan
            iaziat = int(numpy.where(numpy.logical_and(90 > (azibin),
                                     90 < (azibin + dazi)))[0])
            Iaziatnot = numpy.logical_not(numpy.arange(len(azibin))==iaziat)
            # 2. Interpolate input data in the along track direction (will be to avoid as the noise cnnot be removed in that direction)
            usrazi[iaziat] = numpy.interp(azibin[iaziat], azibin[Iaziatnot],
                                          usrazi[Iaziatnot], period=180)
            # 3. Perform an estimate of the Stokes drift direction and magnitude
            c = numpy.fft.fft(usrazi)/len(usrazi)
            rd2 = -0.5 * numpy.angle(c[1]) # estimate of Us direction
            thrbeam = numpy.mod(rd2, 2 * numpy.pi)
            if numpy.cos(wnd_dir[isample, ibeam] - thrbeam) < 0.:
                thrbeam = numpy.arctan2(-numpy.sin(thrbeam),
                                        -numpy.cos(thrbeam))
            # 4. Ambiguity removal using closeness to estimated direction
            radazibin = numpy.deg2rad(azibin[iazi])
            cond = numpy.where(numpy.cos(thrbeam - radazibin) > 0)[0]
            usfull = - usrazi
            if cond.any():
                usfull[(cond)] = usrazi
            usr_comb[isample, ibeam] = numpy.nansum(usfull * numpy.cos(rangle)) * 2 * dazi / 180.
    return usr_comb


def find_closest(lon, lat, lon_nadir, lat_nadir, mss, hs, beam_angle):
    ''' Find closest mss at 6degree or nadir.
    '''
    nsample, nbeam = numpy.shape(mss[:, 1:])
    mssclose = numpy.full((nsample, nbeam), numpy.nan)
    hsclose = numpy.full((nsample, nbeam), numpy.nan)
    ind_6 = numpy.where(numpy.array(beam_angle) == 6)[0]
    mss6 = mss[:, (ind_6 + 1)].ravel()
    mss_nadir = mss[:, 0]
    for isample in range(nsample):
        for ibeam in range(nbeam):
            if not numpy.isfinite(mss[isample, ibeam]):
                continue
            plon = lon[isample, ibeam]
            plat = lat[isample, ibeam]
            dist_nadir = mod_tools.dist_sphere(plon, lon_nadir, plat, lat_nadir)
            #dist_nadir[numpy.where(numpy.isnan(mss_nadir))] = numpy.nan
            # TODO remove ugly try except
            try:
                inadir = numpy.nanargmin(dist_nadir)
            except:
                continue
            dnadir = dist_nadir[inadir]
            dist_6 = mod_tools.dist_sphere(plon, lon[:, (ind_6)].ravel(), plat,
                                           lat[:, (ind_6)].ravel())
            try:
                i6 = numpy.nanargmin(dist_6)
            except:
                continue
            d6 = dist_6[i6]
            hsclose[isample, ibeam] = hs[inadir]
            if d6 < dnadir:
                mssclose[isample, ibeam] = mss6[i6]
            else:
                mssclose[isample, ibeam] = mss_nadir[inadir]
    return mssclose, hsclose


def estimate_uwd(usr, output_var, hs_nadir, mssclose, radial_angle,
                 beam_angle):
    import skimulator.build_error as build_error
    nsample, nbeam = numpy.shape(usr)
    uwdr_est = []
    uwdr_est.append(numpy.full((nsample), numpy.nan))
    for ibeam in range(nbeam):
        from matplotlib import pyplot
        '''
        pyplot.figure()
        pyplot.plot(usr[:, ibeam], label='usr')
        pyplot.plot(output_var['uwnd'][ibeam + 1], label='wnd')
        pyplot.plot(hs_nadir[:, ibeam], label='hs')
        pyplot.plot(mssclose[:, ibeam], label='mss')
        pyplot.legend()
        pyplot.savefig('test.png')
        '''
        output_var2 ={}
        output_var2['ussr'] = usr[:, ibeam]
        output_var2['uwnd'] = output_var['uwnd'][ibeam + 1]
        output_var2['vwnd'] = output_var['vwnd'][ibeam + 1]
        output_var2['hs'] = hs_nadir[:, ibeam]
        output_var2['mssclose'] = mssclose[:, ibeam]
        est_wd = build_error.compute_wd_ai_par
        uwdr_est.append(est_wd(output_var2, radial_angle[:, ibeam],
                                    beam_angle[ibeam]))
    return uwdr_est
