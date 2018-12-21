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
'''
Module to read and write data
Contains functions to read variables and coordinates from a netcdf files. \n
Contains model classes: \n
                       -ROMS \n
                       -NEMO \n
                       -NETCDF_MODEL \n
                       -WW3 \n
Contains satellite class: Sat_SKIM \n
Contains file instrumentation class: file_instr \n
'''
import skimulator
import skimulator.grid_check
from netCDF4 import Dataset
import numpy
import sys
import time as ti
import logging
import datetime
import os
version = skimulator.__version__
logger = logging.getLogger(__name__)


class IncompatibleGridError(Exception):
    """Raised"""
    def __init__(self, path, grid_hash, params_hash, *args, **kwargs):
        """"""
        self.path = path
        self.grid = skimulator.grid_check.revert_b64_gzipped_hash(grid_hash)
        self.p = skimulator.grid_check.revert_b64_gzipped_hash(params_hash)


def write_params(params, pfile):
    """ Write parameters that have been selected to run swot_simulator. """
    with open(pfile, 'w') as f:
        for key in dir(params):
            if not key[0:2] == '__':
                f.write('{} = {}\n'.format(key, params.__dict__[key]))


def read_coordinates(nfile,  nlon, nlat, twoD=True):
    ''' General routine to read coordinates in a netcdf file. \n
    Inputs are file name, longitude name, latitude name. \n
    Outputs are longitudes and latitudes (2D arrays).'''

    # - Open Netcdf file
    try:
        fid = Dataset(nfile, 'r')
    except IOError:
        logger.error('There was an error opening the file {}'.format(nfile))
        sys.exit(1)

    # - Read 1d or 2d coordinates
    try:
        vartmp = fid.variables[nlat]
    except:
        logger.error('Coordinates {} not found in file {}'.format(nlat, nfile))
        sys.exit(1)
    try:
        vartmp = fid.variables[nlon]
    except:
        logger.error('Coordinates {} not found in file {}'.format(nlon, nfile))
        sys.exit(1)
    if len(vartmp.shape) == 1:
        lon_tmp = numpy.array(fid.variables[nlon][:]).squeeze()
        lat_tmp = numpy.array(fid.variables[nlat][:]).squeeze()
        if twoD:
            lon, lat = numpy.meshgrid(lon_tmp, lat_tmp)
        else:
            lon = lon_tmp
            lat = lat_tmp
    elif len(vartmp.shape) == 2:
        lon = numpy.array(fid.variables[nlon][:, :]).squeeze()
        lat = numpy.array(fid.variables[nlat][:, :]).squeeze()
        if not twoD:
            lon = lon[0, :]
            lat = lat[:, 0]
    else:
        logger.warn('unknown dimension for lon and lat')
    fid.close()
    return lon, lat


def read_var(nfile, var, index=None, time=0, depth=0, model_nan=None):
    ''' General routine to read variables in a netcdf file. \n
    Inputs are file name, variable name, index=index to read part
    of the variables, time=time to read a specific time, depth=depth to read a
    specific depth, model_nan=nan value '''

    # - Open Netcdf file
    try:
        fid = Dataset(nfile, 'r')
    except IOError:
        logger.error('There was an error opening the file {}'.format(nfile))
        sys.exit(1)

    # - Check dimension of variable
    try:
        vartmp = fid.variables[var]
    except:
        logger.error('Variable {} not found in file {}'.format(var, nfile))
        sys.exit(1)
    # - Read variable
    if index is None:
        if len(vartmp.shape) == 1:
            T = numpy.array(fid.variables[var][:]).squeeze()
        elif len(vartmp.shape) == 2:
            T = numpy.array(fid.variables[var][:, :]).squeeze()
        elif len(vartmp.shape) == 3:
            if time is None:
                T = numpy.array(fid.variables[var][:, :, :]).squeeze()
            else:
                T = numpy.array(fid.variables[var][time, :, :]).squeeze()
        elif len(vartmp.shape) == 4:
            if time is None:
                if depth is None:
                    T = numpy.array(fid.variables[var][:, :, :, :]).squeeze()
                else:
                    T = numpy.array(fid.variables[var][:, depth,
                                    :, :]).squeeze()
            elif depth is None:
                T = numpy.array(fid.variables[var][time, :, :, :]).squeeze()
            else:
                T = numpy.array(fid.variables[var][time,depth, :, :]).squeeze()
        else:
            logger.error('wrong dimension in variables {}'.format(var))
            sys.exit(1)
    else:
            if len(vartmp.shape) == 1:
                Ttmp = numpy.array(fid.variables[var][:]).squeeze()
                T = Ttmp[index]
            elif len(vartmp.shape) == 2:
                Ttmp = numpy.array(fid.variables[var][:, :]).squeeze()
                T = Ttmp[index]
            elif len(vartmp.shape) == 3:
                if time is None:
                    U = numpy.array(fid.variables[var][:, :, :]).squeeze()
                    T = U[:, index]
                else:
                    U = numpy.array(fid.variables[var][time, :, :]).squeeze()
                    T = U[index]
            elif len(vartmp.shape) == 4:
                if time is None:
                    if depth is None:
                        U = numpy.array(fid.variables[var][:, :, :,
                                                           :]).squeeze()
                        T = U[:, :, index]
                    else:
                        U = numpy.array(fid.variables[var][:, depth, :,
                                                           :]).squeeze()
                        T = U[:, index]
                elif depth is None:
                    U = numpy.array(fid.variables[var][time, :, :,
                                                       :]).squeeze()
                    T = U[:, index]
                else:
                    U = numpy.array(fid.variables[var][time, depth, :,
                                                       :]).squeeze()
                    T = U[index]
            else:
                logger.error('Wrong dimension')
                sys.exit(1)

    fid.close()
    # - Mask value that are NaN
    if model_nan is not None:
        T[numpy.where(T == model_nan)] = numpy.nan
    return T


def write_l2c(metadata, geolocation, **kwargs):
    ti = datetime.datetime.now()
    lat = geolocation['lat']
    lon = geolocation['lon']
    tim = geolocation['time']
    # - Open Netcdf file in write mode
    fid = Dataset(metadata['file'], 'w', format='NETCDF4_CLASSIC')
    # - Create Global attribute
    fid.title = 'SKIM L2C simulated by SKIM simulator'
    fid.keywords = 'SKIM, Doppler'  # Check keywords
    fid.Conventions = "CF-1.6"
    fid.summary = 'SKIM data produced'
    fid.description = "SKIM fixed swath"
    fid.Metadata_Conventions = "Unidata Dataset Discovery v1.0"
    fid.history = 'L2C File created by skimulator version ' + version
    fid.processing_level = 'L2C'
    fid.standard_name_vocabulary = "CF-1.6"
    fid.creator_name = "Lucile Gaultier"
    fid.creator_email = "lucile.gaultier@gmail.com"
    fid.publisher_url = ""
    fid.time_coverage_start = metadata['time_coverage_start']
    # p.date0+"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
    fid.time_coverage_end = metadata['time_coverage_end']
    # p.date0 +"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
    fid.geospatial_lat_min = "{:.2f}".format(numpy.min(lat))
    fid.geospatial_lat_max = "{:.2f}".format(numpy.max(lat))
    fid.geospatial_lat_units = "degrees_north"
    fid.geospatial_lon_max = "{:.2f}".format(numpy.max(lon))
    fid.geospatial_lon_min = "{:.2f}".format(numpy.min(lon))
    fid.geospatial_lon_units = "degrees_east"
    fid.project = "SKIM"
    fid.date_created = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
    fid.date_modified = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
    fid.keywords_vocabulary = ""
    fid.references = ""
    fid.cycle = "{0:d}".format(int(metadata['cycle']))
    fid.track = "{} th pass".format(metadata['pass'])
    # - Create dimensions
    # if (not os.path.isfile(self.file)):
    dimlon = 'al'
    dimlat = 'ac'
    dimtime = 'time'
    nlon = numpy.shape(lon)[0]
    nlat = numpy.shape(lat)[1]
    ntime = None
    fid.createDimension(dimlat, nlat)
    fid.createDimension(dimtime, ntime)
    # - Create and write Variables
    vtime = fid.createVariable('time', 'f8', (dimtime))
    vtime.axis = "T"
    vtime.units = "days since {}".format(metadata['first_time'])
    vtime.long_name = "Time"
    vtime.standard_name = "time"
    vtime.calendar = "gregorian"
    vtime[:] = tim
    vlon = fid.createVariable('lon', 'f4', (dimtime, dimlat))
    vlon.axis = "X"
    vlon.long_name = "Longitude"
    vlon.standard_name = "longitude"
    vlon.units = "degrees_east"
    vlon[:, :] = lon
    val = fid.createVariable('x_al', 'f4', (dimtime))
    val.axis = "Y"
    val.long_name = "Along track distance from beginning of cycle"
    val.standard_name = "Along track distance"
    val.units = "km"
    val[:] = geolocation['al']
    vac = fid.createVariable('x_ac', 'f4', (dimlat))
    vac.axis = "Y"
    vac.long_name = "Along track distance from nadir"
    vac.standard_name = "Across track distance"
    vac.units = "km"
    vac[:] = geolocation['ac']
    vlat = fid.createVariable('lat', 'f4', (dimtime, dimlat))
    vlat.axis = "Y"
    vlat.long_name = "Latitude"
    vlat.standard_name = "latitude"
    vlat.units = "degrees_north"
    vlat[:, :] = lat
    longname = { "ux_noerr": "Error-free zonal velocity",
                 "uy_noerr": "Error-free meridional velocity",
                "ux_obs": "Observed zonal velocity",
                "uy_obs": "Observed meridional velocity",
                "uac_obs": "Observed across track velocity",
                "ual_obs": "Observed along track velocity",
                "ux_model": "Error-free zonal velocity",
                "uy_model": "Error-free meridional velocity",
                "uac_model": "Error-free across track velocity",
                "ual_model": "Error-free along track velocity",
                "angle": "angle of xac with eastward vector",
                "ux_true": "True zonal velocity",
                "uy_true": "True meridional velocity",
                "x_al": "Along track distance from beginning of cycle",
                "x_ac": "Across track distance from nadir"
                }
    unit = {"ux_noerr": "m/s", "ux_obs": "m/s", "uac_obs": "m/s",
            "uy_noerr": "m/s", "uy_obs": "m/s", "ual_obs": "m/s",
            "angle": "rad", "ux_model": "m/s", "uy_model": "m/s",
            "uac_model": "m/s", "ual_model": "m/s",
            "ux_true": "m/s", "uy_true": "m/s", "x_al": "km", "x_ac": "km"
            }
    for key, value in kwargs.items():
        if value is not None:
            nvar = '{}'.format(key)
            var = fid.createVariable(nvar, 'f4', (dimtime, dimlat),
                                     fill_value=-1.36e9)
            try:
                var.units = unit[str(key)]
            except:
                var.units = ''
            try:
                var.long_name = longname[str(key)]
            except:
                var.long_name = str(key)
            if value.any():
                mask = numpy.isnan(value)
                value[numpy.where(mask)] = -1.36e9
                mask_ind = numpy.where(value < -1e7)
                value[mask_ind] = -1.36e9
                mask_ind = numpy.where(value > 1e7)
                value[mask_ind] = -1.36e9
                mask_ind = numpy.where(value == numpy.PINF)
                value[mask_ind] = -1.36e9
                var[:, :] = value
    fid.close()


class Sat_SKIM():
    ''' Sat_SKIM class: to read and write data that has been
    created by SKIM simulator '''
    def __init__(self, ifile=None, lon=None, lat=None, lon_nadir=None,
                 lat_nadir=None, time=None, cycle=None, al_cycle=None,
                 x_al=None, timeshift=None):
        self.file = ifile
        self.lon = lon
        self.lat = lat
        self.lon_nadir = lon_nadir
        self.lat_nadir = lat_nadir
        self.time = time
        self.cycle = cycle
        self.x_al = x_al
        self.al_cycle = al_cycle
        self.timeshift = timeshift

    def write_swath(self, p, **kwargs):
        '''Write swath location in Satellite grid file sgridfile.\n
        Dimensions are  time (i.e. along track), x_ac (across
        track) and cycle (1). \n
        Variables are longitude, latitude, number of days in a cycle,
        distance crossed in a cycle, time, along track and across track
        distances are stored.'''
        grid_params_hash = skimulator.grid_check.get_b64_gzipped_hash(p)
        # - Open Netcdf file in write mode
        fid = Dataset(self.file, 'w', format='NETCDF4_CLASSIC')
        # - Create Global attribute
        fid.title = 'SKIM swath simulated by SKIM simulator'
        fid.keywords = 'check keywords'  # Check keywords
        fid.Conventions = "CF-1.6"
        fid.summary = 'SKIM data produced'
        fid.description = "SKIM fixed swath"
        fid.Metadata_Conventions = "Unidata Dataset Discovery v1.0"
        fid.history = 'Grid File created by skimulator version ' + version
        fid.processing_level = 'L2'
        fid.standard_name_vocabulary = "CF-1.6"
        fid.creator_name = "Lucile Gaultier"
        fid.creator_email = "lucile.gaultier@gmail.com"
        fid.publisher_url = ""
        fid.time_coverage_start = self.time[0][0]
        # p.date0+"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.time_coverage_end = self.time[-1][-1]
        # p.date0 +"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.geospatial_lat_min = "{:.2f}".format(numpy.min(self.lat[0]))
        fid.geospatial_lat_max = "{:.2f}".format(numpy.max(self.lat[0]))
        fid.geospatial_lat_units = "degrees_north"
        fid.geospatial_lon_max = "{:.2f}".format(numpy.max(self.lon[0]))
        fid.geospatial_lon_min = "{:.2f}".format(numpy.min(self.lon[0]))
        fid.geospatial_lon_units = "degrees_east"
        fid.project = "SKIM"
        fid.date_created = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.date_modified = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.keywords_vocabulary = ""
        fid.references = ""
        fid.cycle = "{0:d}".format(int(self.al_cycle))
        fid.track = "{} th pass".format(self.ipass)
        fid.grid_params_hash = grid_params_hash
        # - Create dimensions
        # if (not os.path.isfile(self.file)):
        dimsample = 'sample'
        maxpos = numpy.argmax(numpy.array(p.list_shift))
        nsample = numpy.shape(self.lon[maxpos + 1])[0]
        fid.createDimension(dimsample, nsample)
        # fid.createDimension('time_nadir', numpy.shape(self.lon)[0])
        dimcycle = 'cycle'
        fid.createDimension(dimcycle, 1)
        nbeam = len(self.list_pos)
        dimnbeam = 'nbeam'
        fid.createDimension(dimnbeam, nbeam)
        # - Create and write Variables
        vtime = fid.createVariable('time', 'f8', (dimsample, dimnbeam))
        vtime.axis = "T"
        vtime.units = "seconds since {}".format(p.first_time)
        vtime.long_name = "Time"
        vtime.standard_name = "time"
        vtime.calendar = "gregorian"
        vtime_nadir = fid.createVariable('time_nadir', 'f8', (dimsample,))
        vtime_nadir.axis = "T"
        vtime_nadir.units = "seconds since {}".format(p.first_time)
        vtime_nadir.long_name = "Time at nadir"
        vtime_nadir.standard_name = "time"
        vtime_nadir.calendar = "gregorian"
        vlon = fid.createVariable('lon', 'f4', (dimsample, dimnbeam))
        vlon.axis = "X"
        vlon.long_name = "Longitude"
        vlon.standard_name = "longitude"
        vlon.units = "degrees_east"
        vlon_nadir = fid.createVariable('lon_nadir', 'f4', (dimsample,))
        vlon_nadir.axis = "X"
        vlon_nadir.long_name = "Longitude at nadir"
        vlon_nadir.standard_name = "longitude"
        vlon_nadir.units = "degrees_east"
        vlat = fid.createVariable('lat', 'f4', (dimsample, dimnbeam))
        vlat.axis = "Y"
        vlat.long_name = "Latitude"
        vlat.standard_name = "latitude"
        vlat.units = "degrees_north"
        vlat_nadir = fid.createVariable('lat_nadir', 'f4', (dimsample,))
        vlat_nadir.axis = "Y"
        vlat_nadir.long_name = "Latitude at nadir"
        vlat_nadir.standard_name = "latitude"
        vlat_nadir.units = "degrees_north"
        vx_al = fid.createVariable('x_al', 'f4', (dimsample, dimnbeam))
        vx_al.units = "km"
        vx_al.long_name = "Along track distance from the nadir"
        # vx_al_tot = fid.createVariable('x_al_total', 'f4', (dimsample,
        #                                dimnbeam))
        # vx_al_tot.units = "km"
        # vx_al_tot.long_name = "Along track distance from the beginning of "\
        #                       "the pass projected on nadir"
        vx_ac = fid.createVariable('x_ac', 'f4', (dimsample, dimnbeam))
        vx_ac.units = "km"
        vx_ac.long_name = "Across track distance from the nadir"
        vx_al_nadir = fid.createVariable('x_al_nadir', 'f4', (dimsample,))
        vx_al_nadir.units = "km"
        vx_al_nadir.long_name = "Nadir along track distance from the"\
                                "beginning of the cycle"
        vangle = fid.createVariable('angle', 'f4', (dimsample, dimnbeam))
        vangle.units = "rad"
        vangle.long_name = "Angle of the beam refered to the across track"\
                           " direction"
        vrangle = fid.createVariable('radial_angle', 'f4',
                                     (dimsample, dimnbeam))
        vrangle.units = "rad"
        vrangle.long_name = "Radial angle refered to (longitude towards east,"\
                            " latitude toward north)"
        for i in range(nbeam + 1):
            if i == 0:
                vtime_nadir[:] = self.time[i][:nsample]
                vlon_nadir[:] = self.lon[i][:nsample]
                vlat_nadir[:] = self.lat[i][:nsample]
                vx_al_nadir[:] = self.x_al[i][:nsample]

            else:
                vtime[:, i - 1] = self.time[i][:nsample]
                vlon[:, i - 1] = self.lon[i][:nsample]
                vlat[:, i - 1] = self.lat[i][:nsample]
                vx_al[:, i - 1] = self.x_al[i][:nsample]
                vx_ac[:, i - 1] = self.x_ac[i][:nsample]
                vangle[:, i - 1] = self.beam_angle[i][:nsample]
                vrangle[:, i - 1] = self.radial_angle[i][:nsample]
        vcycle = fid.createVariable('cycle', 'f4', (dimcycle, ))
        valcycle = fid.createVariable('al_cycle', 'f4', (dimcycle, ))
        vtimeshift = fid.createVariable('timeshift', 'f4', (dimcycle, ))
        vcycle[:] = self.cycle
        vcycle.units = "second"
        vcycle.long_name = "seconds during a cycle"
        valcycle[:] = self.al_cycle
        valcycle.units = "km"
        valcycle.long_name = " Distance travelled during the pass"
        vtimeshift[:] = self.timeshift
        vtimeshift.units = "day"
        vtimeshift.long_name = "Shift time to match model time"
        vlistpos = fid.createVariable('beam_position', 'f4', (dimnbeam, ))
        vlistpos[:] = self.list_pos
        vlistpos.units = ""
        vlistpos.long_name = "Beam position"
        vlistangle = fid.createVariable('beam_angle', 'f4', (dimnbeam, ))
        vlistangle[:] = self.list_angle
        vlistangle.units = ""
        vlistangle.long_name = "Beam angle"
        vincl = fid.createVariable('inclination', 'f4', (dimsample, ))
        vincl.units = "rad"
        vincl.long_name = "Track inclination at nadir"
        # vincl[:] = numpy.array(self.angle)[0,:nsample]
        vincl[:] = self.angle[:nsample]
        fid.close()
        return None

    def write_data(self, p, outdata): # **kwargs):
        '''Write SKIM data in output file file_output
        Dimensions are x_al (along track distance), x_ac (across
        track distance). \n
        Variables are longitude, latitude, index (file number),
        error-free radial velocity (velocity interpolated from the model and
        projected with the radial angle), selected
        errors (instrument, uss bias, radial uss) and velocity with errors. \n
        '''
        # - Open netcdf file in write mode
        fid = Dataset(self.file, 'w', format='NETCDF4_CLASSIC')
        fid.description = "Ouptut from SKIM simulator"
        try:
            fid.corresponding_grid = self.gridfile
        except:
            pass
        fid.title = 'SKIM-like data simulated by SKIM simulator'
        fid.keywords = 'SKIM, altimetry, satellite, remote sensing'
        fid.Conventions = "CF-1.6"
        fid.summary = 'SKIM grid data produced'
        fid.description = "SKIM fixed grid"
        fid.Metadata_Conventions = "Unidata Dataset Discovery v1.0"
        fid.history = 'Grid File created by skimulator version ' + version
        fid.processing_level = 'L2'
        fid.standard_name_vocabulary = "CF-1.6"
        fid.creator_name = "Lucile Gaultier"
        fid.creator_email = "lucile.gaultier@gmail.com"
        fid.publisher_url = ""
        dateformat = '%Y-%m-%dT%H:%M:%SZ'
        time_model = datetime.datetime.strptime(p.first_time, '%Y-%m-%dT%H:%M:%SZ')
        mintime = numpy.nanmin(self.time)
        maxtime = numpy.nanmax(self.time)
        day = numpy.floor(mintime)
        seconds = (mintime - day) * 86400
        time0 = time_model + datetime.timedelta(day, seconds)
        fid.time_coverage_start = time0.strftime(format=dateformat)
        # p.date0+"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        day = numpy.floor(maxtime)
        seconds = (maxtime - day) * 86400
        time0 = time_model + datetime.timedelta(day, seconds)
        fid.time_coverage_end = time0.strftime(format=dateformat)
        # p.date0 +"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.geospatial_lat_min = "{:.2f}".format(numpy.min(self.lat))
        fid.geospatial_lat_max = "{:.2f}".format(numpy.max(self.lat))
        fid.geospatial_lat_units = "degrees_north"
        fid.geospatial_lon_max = "{:.2f}".format(numpy.max(self.lon))
        fid.geospatial_lon_min = "{:.2f}".format(numpy.min(self.lon))
        fid.geospatial_lon_units = "degrees_east"
        fid.project = "SKIM"
        fid.date_created = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.date_modified = ti.strftime("%Y-%m-%dT%H:%M:%SZ")
        fid.keywords_vocabulary = ""
        fid.references = ""
        fid.cycle = "{} th cycle".format(self.ncycle)
        fid.track = "{} th pass".format(self.ipass)
        dimsample = 'sample'
        fid.createDimension(dimsample, numpy.shape(self.lon[0])[0])
        # fid.createDimension('time_nadir', numpy.shape(self.lon)[0])
        dimcycle = 'cycle'
        fid.createDimension(dimcycle, 1)
        nbeam = len(p.list_pos)
        dimnbeam = 'nbeam'
        fid.createDimension(dimnbeam, nbeam)
        # - Create and write Variables
        vtime = fid.createVariable('time', 'f8', (dimsample, dimnbeam))
        vtime.axis = "T"
        vtime.units = "days since {}".format(p.first_time)
        vtime.long_name = "Time"
        vtime.standard_name = "time"
        vtime.calendar = "gregorian"
        vtime_nadir = fid.createVariable('time_nadir', 'f8', (dimsample,))
        vtime_nadir.axis = "T"
        vtime_nadir.units = "days since {}".format(p.first_time)
        vtime_nadir.long_name = "Time at nadir"
        vtime_nadir.standard_name = "time"
        vtime_nadir.calendar = "gregorian"
        vlon = fid.createVariable('lon', 'f4', (dimsample, dimnbeam))
        vlon.axis = "X"
        vlon.long_name = "Longitude"
        vlon.standard_name = "longitude"
        vlon.units = "degrees_east"
        vlon_nadir = fid.createVariable('lon_nadir', 'f4', (dimsample,))
        vlon_nadir.axis = "X"
        vlon_nadir.long_name = "Longitude at nadir"
        vlon_nadir.standard_name = "longitude"
        vlon_nadir.units = "degrees_east"
        vlat = fid.createVariable('lat', 'f4', (dimsample, dimnbeam))
        vlat.axis = "Y"
        vlat.long_name = "Latitude"
        vlat.standard_name = "latitude"
        vlat.units = "degrees_north"
        vlat_nadir = fid.createVariable('lat_nadir', 'f4', (dimsample,))
        vlat_nadir.axis = "Y"
        vlat_nadir.long_name = "Latitude at nadir"
        vlat_nadir.standard_name = "latitude"
        vlat_nadir.units = "degrees_north"
        for i in range(nbeam + 1):
            if i == 0:
                vtime_nadir[:] = self.time[i][:]
                vlon_nadir[:] = self.lon[i][:]
                vlat_nadir[:] = self.lat[i][:]

            else:
                vtime[:, i - 1] = self.time[i][:]
                vlon[:, i - 1] = self.lon[i][:]
                vlat[:, i - 1] = self.lat[i][:]
        longname = {"sigma0": "sigma0",
                    "ur_true": "Radial velocity interpolated from model",
                    "ucur": "Zonal velocity interpolated from model",
                    "vcur": "Meridional velocity interpolated from model",
                    "ur_obs": "Observed radial velocity (Ur_model+errors)",
                    "index": "Equivalent model output number in list of file",
                    "ur_uss": "Stokes drift radial velocity bias",
                    "uwb": "Wave bias componant",
                    "uwb_corr": "Remaining wave bias current after inversion",
                    "uwnd": "Eastward wind at 10m ",
                    "vwnd": "Northward wind at 10m ",
                    "nadir_err": "Nadir error",
                    "wlv": "SSH interpolated from model",
                    "ssh_obs": "Observed SSH",
                    "ice": "Sea ice concentration",
                    "uwb": "Current Wave bias"
                    }
        unit = {"sigma0": "", "ur_true": "m/s", "ur_obs": "m/s",
                "index": " ", "ur_uss": "m/s", "uwnd": "m/s",
                "vwnd": "m/s", "uwb": "m/s", "ucur": "m/s",
                "vcur": "m/s", "ssh_obs": "m", "wlv": "m",
                "nadir_err": "m", "ssh_obs":"m", "uwb_corr": "m/s"
                }
        list_nadir = ("nadir_err", "ssh_true", "ssh_obs")
        for key, value in outdata.items():
            if value is not None:
                if key in list_nadir or key == "vindice" or key == "instr":
                    nvar_nadir = '{}_nadir'.format(key)
                    var_nadir = fid.createVariable(nvar_nadir, 'f4',
                                                   (dimsample, ),
                                                   fill_value=-1.36e9)
                    try:
                        var_nadir.units = unit[str(key)]
                    except:
                        var_nadir.units = ''
                    try:
                        var_nadir.long_name = longname[str(key)]
                    except:
                        var_nadir.long_name = str(key)
                if key not in list_nadir:
                    nvar = '{}'.format(key)
                    var = fid.createVariable(nvar, 'f4', (dimsample, dimnbeam),
                                             fill_value=-1.36e9)
                    try:
                        var.units = unit[str(key)]
                    except:
                        var.units = ''
                    try:
                        var.long_name = longname[str(key)]
                    except:
                        var.long_name = str(key)
                for i in range(len(value)):
                    if value[i].any():
                        value_tmp = value[i][:]
                        # vmin = numpy.nanmin(value_tmp)
                        # vmax = numpy.nanmax(value_tmp)
                        mask = numpy.isnan(value_tmp)
                        value_tmp[numpy.where(mask)] = -1.36e9
                        mask_ind = numpy.where(value_tmp < -1e7)
                        value_tmp[mask_ind] = -1.36e9
                        mask_ind = numpy.where(value_tmp > 1e7)
                        value_tmp[mask_ind] = -1.36e9
                        mask_ind = numpy.where(value_tmp == numpy.PINF)
                        value_tmp[mask_ind] = -1.36e9
                        if i == 0:
                            if key in list_nadir or key == "vindex" or key == 'instr':
                                var_nadir[:] = value_tmp
                            else:
                                continue
                        else:
                            if key not in list_nadir:
                                var[:, i - 1] = value_tmp

                # try:    var.missing_value = p.model_nan
                # except: var.missing_value = 0.
                # fid.setncattr('missing_value','-9999.f')
        # TODO set range values
        fid.close()
        return None

    def load_swath(self, p, **kwargs):
        '''Load swath variables stored in Satellite grid file sgridfile. \n
        (longitude, latitude, number of days in a cycle, crossed distance
        during a cycle, time, along track and across track position).'''

        # - Open Netcdf file
        try:
            fid = Dataset(self.file, 'r')
        except IOError:
            logger.error('There was an error opening the file '
                         '{}'.format(self.file))
            sys.exit(1)

        if 'grid_params_hash' in fid.ncattrs():
            grid_params_hash = skimulator.grid_check.get_b64_gzipped_hash(p)
            if fid.grid_params_hash != grid_params_hash:
                raise IncompatibleGridError(self.file, fid.grid_params_hash,
                                            grid_params_hash)

        # fid = Dataset(self.file, 'r')
        time = []
        lon = []
        lat = []
        # cycle = []
        # x_al = []
        listvar = {'time': time, 'lon': lon, 'lat': lat, }
        self.lon = []
        self.lat = []
        self.time = []
        # - Read variables in listvar and return them
        for stringvar in listvar:
            var = fid.variables['{}{}'.format(stringvar, '_nadir')]
            listvar[stringvar].append(numpy.array(var[:]).squeeze())
            var = fid.variables[stringvar]
            for i in range(len(p.list_pos)):
                listvar[stringvar].append(numpy.array(var[:, i]).squeeze())
            setattr(self, stringvar, listvar[stringvar])
        # - Read variables in arguments
        for key, value in kwargs.items():
            var = fid.variables[key]
            value = numpy.array(fid.variables[key][:]).squeeze()
            # value[value == var.fill_value] = numpy.nan
            setattr(self, key, value)
        self.pos = numpy.array(fid.variables['beam_position'][:])
        _dpos = numpy.abs(self.pos - numpy.array(p.list_pos,
                                                 dtype=numpy.float32))
        if (numpy.all(_dpos) > 0.0001):
            logger.error('List of beam positions has changed,'
                         ' reprocess the grids')
            sys.exit(1)
        self.angle = numpy.array(fid.variables['beam_angle'][:])
        _dangle = numpy.abs(self.angle - numpy.array(p.list_angle,
                                                     dtype=numpy.float32))
        if (numpy.all(_dangle) > 0.0001):
            logger.error('List of beam angles has changed,'
                         ' reprocess the grids')
            sys.exit(1)
        self.radial_angle = numpy.array(fid.variables['radial_angle'][:])
        self.angle = numpy.array(fid.variables['angle'][:])
        try:
            self.corresponding_grid = fid.corresponding_grid
        except AttributeError:
            pass
        try:
            self.incl = numpy.array(fid.variables['inclination'][:]).squeeze()
        except:
            logger.info('inclination variable not found')
        fid.close()
        return None


    def load_data(self, p, **kwargs):
        '''Load swath variables stored in Satellite grid file sgridfile. \n
        (longitude, latitude, number of days in a cycle, crossed distance
        during a cycle, time, along track and across track position).'''

        # - Open Netcdf file
        try:
            fid = Dataset(self.file, 'r')
        except IOError:
            logger.error('There was an error opening the file '
                         '{}'.format(self.file))
            sys.exit(1)
        # fid = Dataset(self.file, 'r')
        time = []
        lon = []
        lat = []
        # cycle = []
        # x_al = []
        #listvar = {'time': time, 'lon': lon, 'lat': lat, }
        #self.lon = []
        #self.lat = []
        #self.time = []
        # - Read variables in listvar and return them
        #for stringvar in listvar:
        #    var = fid.variables['{}{}'.format(stringvar, '_nadir')]
        #    listvar[stringvar].append(numpy.array(var[:]).squeeze())
        #    var = fid.variables[stringvar]
        #    for i in range(len(p.list_pos)):
        #        listvar[stringvar].append(numpy.array(var[:, i]).squeeze())
        #    setattr(self, stringvar, listvar[stringvar])
        # - Read variables in arguments
        for key, value in kwargs.items():
            var = fid.variables[key]
            value = numpy.array(fid.variables[key][:]).squeeze()
            # value[value == var.fill_value] = numpy.nan
            setattr(self, key, value)
        try:
            self.corresponding_grid = fid.corresponding_grid
        except AttributeError:
            pass
        fid.close()
        return None


class NETCDF_MODEL():
    '''Class to read any netcdf data.\n
    USAGE is NETCDF_MODEL(file=name of file ,var= variable name,
    lon=variable longitude, lat=variable latitude, units=).\n
    Argument file is mandatory, arguments var, lon, lat
    are specified in params file. \n
    '''
    def __init__(self, p, ifile=None, list_input_var=None, lonu='longitude',
                 latu='latitude', lonv='longitude', latv='latitude', depth=0,
                 time=0):
        if p.list_input_var is None:
            logger.error('Specify your list_input_var in parameter file')
            sys.exit(1)
            self.input_var_list = {'ucur': ['uo', ''],
                                   'vcur': ['vo', ''],
                                   }
        else:
            self.input_var_list = p.list_input_var
        self.input_var = {}
        self.nlonu = lonu
        self.nlatu = latu
        self.nlonv = lonv
        self.nlatv = latv
        self.nfile = ifile
        self.depth = depth
        self.time = time
        self.model_nan = getattr(p, 'model_nan', 0)
        p.model_nan = self.model_nan
        logger.debug('Nan Values {}, {}'.format(p.model_nan, self.model_nan))

    def read_var(self, p, index=None):
        '''Read variables from netcdf file \n
        Argument is index=index to load part of the variable.'''
        for key, value in self.input_var_list.items():
            nfile0 = self.nfile[0]
            _nfile = '{}{}.nc'.format(nfile0, value[1])
            if os.path.exists(_nfile):
                self.input_var[key] = read_var(_nfile, value[0], index=index,
                                               time=self.time,
                                               depth=self.depth,
                                               model_nan=self.model_nan)
            else:
                logger.info('{} not found'.format(_nfile))
        # self.vvar[numpy.where(numpy.isnan(self.vvar))]=0
        return None

    def read_coordinates(self, p, index=None):
        '''Read coordinates from netcdf file \n
        Argument is index=index to load part of the variable.'''
        if p.grid == 'regular':
            lonu, latu = read_coordinates(self.nfile[0], self.nlonu, self.nlatu,
                                          twoD=False)
            lonv, latv = read_coordinates(self.nfile[1], self.nlonv, self.nlatv,
                                          twoD=False)
        else:
            lonu, latu = read_coordinates(self.nfile[0], self.nlonu, self.nlatu)
            lonv, latv = read_coordinates(self.nfile[1], self.nlonv, self.nlatv)
        self.vlatu = latu
        self.vlonu = (lonu + 360) % 360
        self.vlatv = latv
        self.vlonv = (lonv + 360) % 360
        return None

    def calc_box(self, p):
        '''Calculate subdomain coordinates from netcdf file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates(p)
        if (numpy.min(self.vlonu) < 1.) and (numpy.max(self.vlonu) > 359.):
            _ind = numpy.where(self.vlonu > 180.)
            self.vlonu[_ind] = self.vlonu[_ind] - 360
            lon1 = (numpy.min(self.vlonu) + 360) % 360
            lon2 = (numpy.max(self.vlonu) + 360) % 360
        else:
            lon1 = numpy.min(self.vlonu)
            lon2 = numpy.max(self.vlonu)
        return [lon1, lon2, numpy.min(self.vlatu), numpy.max(self.vlatu)]


class WW3():
    '''Class to read ww3 netcdf data.\n
    USAGE is NETCDF_MODEL(file=name of file ,var= variable name,
    lon=variable longitude, lat=variable latitude, units=).\n
    Argument file is mandatory, arguments var, lon, lat
    are specified in params file. \n
    '''
    def __init__(self, p, ifile=None, list_input_var=None, lonu='longitude',
                 latu='latitude', lonv='longitude',
                 latv='latitude', depth=0, time=0):
        self.nlonu = 'longitude'
        self.nlatu = 'latitude'
        self.nlonv = self.nlonu
        self.nlatv = self.nlatu
        self.nfile = ifile
        self.depth = depth
        self.time = time
        if list_input_var is None:
            self.input_var_list = {'ucur': ['ucur', 'cur'],
                                   'vcur': ['vcur', 'cur'],
                                   'uuss': ['uuss', 'uss'],
                                   'vuss': ['vuss', 'uss'],
                                   'ice': ['ice', 'ice'],
                                   'mssd': ['mssd', 'msd'],
                                   'mssx': ['mssx', 'mss'],
                                   'mssy':['mssy', 'mss'],
                                   'ssh': ['wlv', 'wlv'],
                                   'uwnd': ['uwnd', 'wnd'],
                                   'vwnd': ['vwnd', 'wnd']}
        else:
            self.input_var_list = p.list_input_var
        self.input_var = {}
        self.model_nan = getattr(p, 'model_nan', 0.)
        p.model_nan = self.model_nan
        logger.debug('Nan Values {}, {}'.format(p.model_nan, self.model_nan))


    def read_var(self, p, index=None):
        '''Read variables from netcdf file \n
        Argument is index=index to load part of the variable.'''
        for key, value in self.input_var_list.items():
            nfile0 = self.nfile[0]
            _nfile = '{}{}.nc'.format(nfile0[:-6], value[1])
            if os.path.exists(_nfile):
                self.input_var[key] = read_var(_nfile, value[0], index=index,
                                               time=self.time,
                                               depth=self.depth,
                                               model_nan=self.model_nan)
        return None


    def read_coordinates(self, p, index=None):
        '''Read coordinates from netcdf file \n
        Argument is index=index to load part of the variable.'''
        if p.grid == 'regular':
            lonu, latu = read_coordinates(self.nfile[0], self.nlonu, self.nlatu,
                                          twoD=False)
            lonv, latv = read_coordinates(self.nfile[1], self.nlonv, self.nlatv,
                                          twoD=False)
        else:
            lonu, latu = read_coordinates(self.nfile, self.nlonu, self.nlatu)
            lonv, latv = read_coordinates(self.nfile, self.nlonv, self.nlatv)
        self.vlatu = latu
        self.vlonu = (lonu + 360) % 360
        self.vlatv = latv
        self.vlonv = (lonv + 360) % 360
        return None


    def calc_box(self, p):
        '''Calculate subdomain coordinates from netcdf file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates(p)
        if (numpy.min(self.vlonu) < 1.) and (numpy.max(self.vlonu) > 359.):
            _ind = numpy.where(self.vlonu > 180.)
            self.vlonu[_ind] = self.vlonu[_ind] - 360
            lon1 = (numpy.min(self.vlonu) + 360) % 360
            lon2 = (numpy.max(self.vlonu) + 360) % 360
        else:
            lon1 = numpy.min(self.vlonu)
            lon2 = numpy.max(self.vlonu)
        return [lon1, lon2, numpy.min(self.vlatu), numpy.max(self.vlatu)]


    def compute_mss(self):
        required = ('mssx', 'mssy', 'mssd', 'uwnd', 'vwnd', 'uwnd', 'ucur',
                    'vcur')
        missing = [_ for _ in required if _ not in self.input_var.keys()]
        if 0 < len(missing):
            logger.info('Missing file to compute sigma, instrumental error not'
                        ' computed')
            logger.info('Missing parameters: {}'.format(', '.join(missing)))
            return None
        else:
            mssd = numpy.deg2rad(self.input_var['mssd'])
            mssu = self.input_var['mssx']
            mssc = self.input_var['mssy']
            uwnd = self.input_var['uwnd']
            vwnd = self.input_var['vwnd']
            ucur = self.input_var['ucur']
            vcur = self.input_var['vcur']
            mssxl = mssu * numpy.cos(mssd)**2 + mssc * numpy.sin(mssd)**2
            mssyl = mssu * numpy.sin(mssd)**2 + mssc * numpy.cos(mssd)**2
            mssxyl = (mssu - mssc) * numpy.sin(2 * mssd) / 2
            nwr = numpy.sqrt((uwnd - ucur)**2 + (vwnd - vcur)**2)
            wrd = numpy.pi / 2 - numpy.arctan2(vwnd - vcur, uwnd - ucur)
            mssshort = numpy.log(nwr + 0.7) * 0.009
            # Replace nan values by 0 to avoid a runtime warning (nan values
            # will be restored afterwards)
            mssshort_nanmask_ind = numpy.where(numpy.isnan(mssshort))
            mssshort[mssshort_nanmask_ind] = 0
            mssshort[mssshort < 0] = 0
            mssshort[mssshort_nanmask_ind] = numpy.nan  # restore nan values
            #Directionality for short wave mss (if 0.5: isotrophic)
            facssdw = 0.6
            mssds = facssdw * mssshort
            msscs = mssshort - mssds
            mssxs = msscs * numpy.sin(wrd)**2 + mssds * numpy.cos(wrd)**2
            mssys = mssds * numpy.sin(wrd)**2 + msscs * numpy.cos(wrd)**2
            mssxys = abs(mssds - msscs) * numpy.sin(2* wrd)
            self.input_var['mssx'] = mssxs + mssxl
            self.input_var['mssy'] = mssys + mssyl
            self.input_var['mssxy'] = mssxys + mssxyl
            del self.input_var['mssd']


