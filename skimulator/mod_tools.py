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
''' Spectral  and algebra tools  for SKIM simulator. \n
Contains the following functions:
- load_python_file: load and parse parameter file \n
- gen_signal1d: compute 1d random error from spectrum \n
- gen_signal2d: compute 2d random error from spectrum \n
- rotationmat3d: rotate data  \n
- spher2cart: convert spherical to cartesian coordinates \n
- cart2spher: convert cartesian to spherical coordinates \n
- cart2spher: convert cartesian to spherical coordinates for vectors \n
- proj_radial: projection on radial direction \n
- update_progress: Progress bar'''

import numpy
import math
import logging
import sys
import os
import types
import datetime
import skimulator.const as const

# Define logger level for debug purposes
logger = logging.getLogger(__name__)


def load_python_file(file_path):
    """Load a file and parse it as a Python module."""
    if not os.path.exists(file_path):
        raise IOError('File not found: {}'.format(file_path))

    full_path = os.path.abspath(file_path)
    python_filename = os.path.basename(full_path)
    module_name, _ = os.path.splitext(python_filename)
    module_dir = os.path.dirname(full_path)
    if module_dir not in sys.path:
        sys.path.append(module_dir)

    module = __import__(module_name, globals(), locals(), [], 0)
    return module


def initialize_parameters(p):
    p.shift_lon = getattr(p, 'shift_lon', None)
    p.shift_time = getattr(p, 'shift_time', None)
    if p.shift_time is None:
        p.timeshift = 0
    p.timeshift = getattr(p, 'timeshift', 0)
    if p.shift_time is None:
        p.timeshift = 0
    model = getattr(p, 'model', 'NETCDF_MODEL')
    p.model = model
    p.model_nan = getattr(p, 'model_nan', 0)
    p.vel_factor = getattr(p, 'vel_factor', 1.)
    p.first_time = getattr(p, 'first_time', '2018-01-01T00:00:00Z')
    p.list_input_var = getattr(p, 'list_input_var', None)
    p.grid = getattr(p, 'grid', 'regular')
    p.cycle = getattr(p, 'cycle', 0.0368)
    p.order_orbit_col = getattr(p, 'order_orbit_col', None)
    p.satcycle = getattr(p, 'satcycle', None)
    p.sat_elev = getattr(p, 'sat_elev', None)
    p.ice_mask = getattr(p, 'ice_mask', True)
    listo = ['wlv', 'ssh_obs', 'ur_true', 'ucur', 'vcur', 'uuss', 'vuss',
             'radial_angle', 'vwnd', 'mssx', 'mssy', 'mssxy', 'uwb',
             'vindice', 'ur_obs', 'mask', 'uwnd', 'sigma0', 'ice']
    p.list_output = getattr(p, 'list_output', listo)
    p.uwb =  getattr(p, 'uwb', True)  #TODO: remove this
    p.snr_coeff = getattr(p, 'snr_coeff', 3e-3)
    p.nadir = getattr(p, 'nadir', False)
    p.ice = getattr(p, 'ice', False)
    p.proc_count = getattr(p, 'proc_number', 1)
    p.progress_bar = getattr(p, 'progress_bar', True)
    p.ac_threshold = getattr(p, 'ac_threshold', 20)  # in km
    #p.resol = getattr(p, 'resol_spatial_l2c', 40) # in km
    p.resol = getattr(p, 'resol', 40) # in km
    p.posting = getattr(p, 'posting_l2c', 5) # in km
    p.resol_spatial_l2d = getattr(p, 'resol_spatial_l2d', 50) # in km
    p.resol_temporal_l2d = getattr(p, 'resol_temporal_l2d', 8) # in days
    p.posting_l2d = getattr(p, 'posting_l2d', (0.1, 0.1)) # in degrees
    # Time domain: (start_time, end_time, dtime) in days
    p.time_domain = getattr(p, 'time_domain', (5, 25, 1)) # in days
    p.list_input_var_l2c = getattr(p, 'list_input_var_l2c', p.list_input_var)
    p.list_input_var_l2d = getattr(p, 'list_input_var_l2d', p.list_input_var)
    p.config_l2d = getattr(p, 'config_l2d', '')
    p.config_l2c = getattr(p, 'config_l2c', '')
    p.config = getattr(p, 'config', 'SKIM')
    p.rain = getattr(p, 'rain', False)
    p.rain_file = getattr(p, 'rain_file', None)
    p.rain_threshold = getattr(p, 'rain_threshold', 0.1)
    p.attitude = getattr(p, 'attitude', False)
    p.yaw_file = getattr(p, 'yaw_file', None)
    p.instr_configuration = getattr(p, 'instr_configuration', 'A')
    make_list_output(p)
    return None


def make_list_output(p):
    #if p.rain is True:
    compulsory = ('ur_true', 'ur_obs', 'radial_angle', 'vindice', 'ur_obs')
    for key in compulsory:
        if key not in p.list_output:
            p.list_output.append(key)
    if p.nadir is True:
        compulsory = ('ssh_obs',)
        for key in compulsory:
            if key not in p.list_output:
                p.list_output.append(key)
    if p.instr is True:
        compulsory = ('instr', 'mssxy', 'sigma0')
        for key in compulsory:
            if key not in p.list_output:
                p.list_output.append(key)
    if p.uwb is True:
        compulsory = ('mssc', 'mssu', 'ur_uss', 'uwb_noerr', 'uwb')
        for key in compulsory:
            if key not in p.list_output:
                p.list_output.append(key)
    return None


def check_path(p):
    if os.path.isdir(p.dir_setup) is False:
        logger.error('Data directory {} not found'.format(p.dir_setup))
        sys.exit(1)
    if os.path.isdir(p.indatadir) is False:
        logger.error('Input directory {} not found'.format(p.indatadir))
        sys.exit(1)
    if os.path.isdir(p.outdatadir) is False:
        logger.warn('Output directory {} did not exist and was '
                    'created'.format(p.outdatadir))
        os.makedirs(p.outdatadir)
    filesat_path = os.path.join(p.dir_setup, p.filesat)
    if os.path.isfile(filesat_path) is False:
        logger.error('Orbit file {} not found'.format(filesat_path))
        sys.exit(1)
    if p.file_input is not None:
        if os.path.isfile(p.file_input) is False:
            logger.error('Model file list {} not found'.format(p.file_input))
            sys.exit(1)
    return None


def gen_coeff_signal1d(f, PS, nc):
    '''Generate nc random coefficient from a spectrum PS
    with frequencies f. \n
    Return Amplitude, phase and frequency of nc realisations'''

    f1 = min(f)
    f2 = max(f)
    logf = numpy.log10(f)
    logf1 = numpy.log10(f1)
    logf2 = numpy.log10(f2)
    # ''' Compute nc random vectors in [logf1, logf2] '''
    logfr = (logf2 - logf1)*numpy.random.random(nc) + logf1
    fr = 10**(logfr)
    # ''' Compute amplitude for random vectors '''
    logPS = numpy.log10(PS)
    logPSr = numpy.interp(logfr, logf, logPS)
    PSr = 10**(logPSr)
    A = numpy.sqrt(0.5 * PSr * fr * ((f2/f1)**(1./nc)-1))
    # ''' Compute nc phases in (0,2*pi)'''
    phi = 2 * math.pi * numpy.random.random(nc)
    return A, phi, fr


def gen_coeff_signal2d(f, PS, nc):
    '''Generate nc random coefficient from a spectrum PS
    with frequencies f. \n
    Inputs are: frequency [f], spectrum [PS], number of realisation [nc]
    Return Amplitude, phase and frequency in 2D (frx, fry) of nc
    realisations'''

    # ''' Compute nc random vectors in an annular
    # (radius is between (min(f), max(f)) '''
    f1 = min(f)
    f2 = max(f)
    logf = numpy.log10(f)
    fr = (f2 - f1)*numpy.random.random(nc) + f1
    logfr = numpy.log10(fr)
    dir = 2. * math.pi*numpy.random.random(nc)
    frx = fr * numpy.cos(dir)
    fry = fr * numpy.sin(dir)

    # ''' Compute coefficients corresponding to random vectors '''
    logPS = numpy.log10(PS)
    logPSr = numpy.interp(logfr, logf, logPS)
    PSr = 10.**logPSr
    A = numpy.sqrt(0.5*PSr*2*math.pi*(f2 - f1)/(nc))

    # ''' Compute nc phases in (0,2*pi)'''
    phi = 2*math.pi*numpy.random.random(nc)

    return A, phi, frx, fry


def rotationmat3D(theta, axis):
    ''' Creates a rotation matrix: Slow method. \n
    Inputs are rotation angle theta and rotation axis axis.
    The rotation matrix correspond to a rotation of angle theta
    with respect to axis axis. \n
    Return the rotation matrix.'''
    # mat = numpy.eye(3, 3)
    axis = axis / math.sqrt(numpy.dot(axis, axis))
    a = math.cos(theta/2.)
    b, c, d = -axis*math.sin(theta/2.)

    return numpy.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
                       [2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
                       [2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c]])


def sign(a):
    return (a > 0) - (a < 0)


def spher2cart(lon, lat):
    ''' Convert spherical coordinates to cartesian coordinates.\n
    Inputs are longitude, latitude. \n
    Return x, y, z'''
    x = numpy.cos(lon*math.pi/180.) * numpy.cos(lat*math.pi/180.)
    y = numpy.sin(lon*math.pi/180.) * numpy.cos(lat*math.pi/180.)
    z = numpy.sin(lat*math.pi/180.)
    return x, y, z


def cart2sphervect(x, y, z):
    ''' Convert cartesian coordinates to spherical coordinates. \n
    Inputs are cartiesian coordinates x, y, z. \n
    Return lon, lat. '''
    norm = numpy.sqrt(x*x + y*y + z*z)
    lat = numpy.arcsin(z/norm) * 180./math.pi
    lon = numpy.arctan(y/x) % (2*math.pi)
    if (x < 0).any():
        lon[x < 0] = (numpy.arctan(y[x < 0] / x[x < 0]) % (2*math.pi)
                      + x[x < 0] / x[x < 0]*math.pi)
    lon = lon * 180/math.pi
    return lon % 360, lat


def cart2spher(x, y, z):
    ''' Convert cartesian coordinates to spherical coordinates. \n
    Inputs are cartiesian coordinates x, y, z. \n
    Return lon, lat. '''
    norm = numpy.sqrt(x*x + y*y + z*z)
    lat = numpy.arcsin(z/norm) * 180./math.pi
    if (x < 0):
        lon = (numpy.arctan(y/x) % (2*math.pi)
               + max(-numpy.sign(x), 0)*math.pi)
    else:
        lon = numpy.arctan(y/x) % (2*math.pi)

    lon = lon * 180/math.pi
    return lon % 360, lat

def dist_sphere(lon1,lon2,lat1,lat2):
    R = const.Rearth
    radlat1 = numpy.deg2rad(lat1)
    radlat2 = numpy.deg2rad(lat2)
    raddifflon = numpy.deg2rad(lon2 - lon1)
    d = numpy.arccos(numpy.sin(radlat2) * numpy.sin(radlat1)
                     + numpy.cos(radlat2) * numpy.cos(radlat1)
                     * numpy.cos(raddifflon))
    return d * R / 1000. # in km


def proj_radial(u, v, radial_angle):
    ur = u * numpy.cos(radial_angle) + v * numpy.sin(radial_angle)
    return ur

def cross_product(mat, ncoeff, cshape):
    ncross = sum(range(ncoeff + 1)) + ncoeff
    cross = numpy.full((cshape, ncross), numpy.nan)
    cross[:, 0:ncoeff] = mat[:, :ncoeff]
    niter = ncoeff
    for k in range(ncoeff):
        for l in range(ncoeff - k):
            cross[:, niter] = mat[:, k]*mat[:, k + l]
            niter += 1
    return cross


#def reconstruct_var(ncoeff, cshape, f1, b1, cross, xlabel):
#    shape_f1 = numpy.shape(f1)[0]
#    shape_b1 = numpy.shape(b1)[0]
#    #f1 = f1.reshape(int(shape_f1//float(shape_b1)), shape_f1)
#    proba = numpy.exp(numpy.dot(cross, f1) + b1)
#    proba_tot = numpy.tile(numpy.nansum(proba, axis=1), (shape_b1, 1))
#    proba = proba / numpy.transpose(proba_tot)
#    var = numpy.sum(proba * numpy.tile(xlabel, (cshape, 1)),axis=1)
 #   return var


def reconstruct_var(ncoeff, cshape, f1, b1, cross, xlabel):
    shape_f1 = numpy.shape(f1)[0]
    shape_b1 = numpy.shape(b1)[0]
    #f1 = f1.reshape(int(shape_f1//float(shape_b1)), shape_f1)
    proba = numpy.exp(numpy.dot(cross, f1) + b1)
    proba_tot = numpy.tile(numpy.nansum(proba, axis=1), (shape_b1, 1))
    proba = proba / numpy.transpose(proba_tot)
    var = numpy.sum(proba * numpy.tile(xlabel, (cshape, 1)),axis=1)
    b1 = numpy.tile(b1, (1, 1))
    train_weight = numpy.concatenate((f1,b1),axis=0)
    l_train = numpy.zeros((cshape, ncoeff + 1)).astype('float64')
    l_train[:, : ncoeff] = cross
    l_train[:, ncoeff] = 1.
    res = numpy.dot(l_train, train_weight)
    xx2D = numpy.tile(xlabel, (cshape, 1))
    var = (numpy.sum(xx2D * numpy.exp(res), axis=1)
           / numpy.sum(numpy.exp(res),axis=1))
    return var


def todict(p):
    result = {}
    for attr in dir(p):
        value = getattr(p, attr)
        if (isinstance(value, types.ModuleType)
             or isinstance(value, types.MethodType)
             or isinstance(value, types.FunctionType)
             or attr.startswith('__')):
            continue
        result[attr] = value
    return result


def fromdict(result):
    p = type('skimparam', (object,), result)
    return p


def convert_dbkm2ms(var_out, ac_angle, beam_angle):
    # Computation of lobe antenna
    th3 = numpy.deg2rad(0.65)
    an = th3 / numpy.sqrt(4 * numpy.log(2))
    # Computation of mispointing
    delta = an**2 / numpy.sin(numpy.deg2rad(beam_angle))**2 * var_out
    out = delta *const.vsat * numpy.cos(ac_angle)
    return out



def update_progress_multiproc(status, info):
    """Creation of progress bar: print on screen progress of run, optimized
    for parrallelised tasks"""
    if info[3] is not None:
        # There has been an error
        return False
    pid = info[0]
    grid_name = info[1]
    if isinstance(grid_name, str):
        ipass = grid_name[-6:-3]
    else:
        ipass = '{:03d}'.format(grid_name)

    cycle = info[2]

    count = len(status.keys())
    sys.stdout.write(_term_move_up() * count)
    now = datetime.datetime.now().strftime('%H:%M:%S')
    for pid in status:
        if grid_name in status[pid]['jobs']:
            if cycle is None:
                # Grid has been completely processed
                status[pid]['done'] += 1
                status[pid]['extra'] = '{}|> pass: {} ....... DONE'.format(
                                                                    now, ipass)
            else:
                # Just update extra info
                status[pid]['extra'] = '{}|> pass: {}, cycle: {:04d}'.format(
                                                             now, ipass, cycle)

    bar_size = 20
    for pid, proc_state in status.items():
        done = math.floor(bar_size * proc_state['done'] / proc_state['total'])
        todo = bar_size - done
        proc_elems = ['\n[']
        proc_elems.extend(['#'] * int(math.ceil(done)))
        proc_elems.extend([' '] * int(math.ceil(todo)))
        proc_elems.extend(['] {}'.format(proc_state['extra'])])
        sys.stdout.write(''.join(proc_elems))
        sys.stdout.flush()
    return True


def update_progress(progress, arg1, arg2):
    '''Creation of a progress bar: print on screen the progress of the run'''
    barLength = 20  # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    if arg1 and arg2:
        text = "\r[{0}] {1}%, {2}, {3}".format("#"*block + "-"*(barLength-block), "%.2f" % (progress*100), arg1 + ', ' + arg2, status)
    elif arg1:
        text = "\r[{0}] {1}%, {2}, {3}".format("#"*block + "-"*(barLength-block), "%.2f" % (progress*100), arg1, status)
    else:
        text = "\r[{0}] {1}%, {2} ".format("#"*block + "-"*(barLength-block), "%.2f" % (progress*100), status)
    sys.stdout.write(text)
    sys.stdout.flush()
    return progress


def _term_move_up():  # pragma: no cover
    """Borrowed from https://github.com/tqdm/tqdm/blob/master/tqdm/_tqdm.py
    MIT 2016 (c) [PR #96] on behalf of Google Inc.
    MIT 2013 (c) Noam Yorav-Raphael, original author."""
    colorama = None
    return '' if (os.name == 'nt') and (colorama is None) else '\x1b[A'
