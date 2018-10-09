import numpy
import math
import logging

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
