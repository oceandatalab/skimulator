'''
Module to read and write data
Contains functions to read variables and coordinates from a netcdf files. \n
Contains model classes: \n
                       -ROMS \n
                       -NEMO \n
                       -NETCDF_MODEL \n
Contains satellite class: Sat_SKIM \n
Contains file instrumentation class: file_instr \n
'''
netcdf4 = True
version = '2.21'
try:
    from netCDF4 import Dataset
except ImportError:
    print('WARNING: package netCDF4 missing, scipy.io.netcdf is used instead'\
          'of netCDF4')
    from scipy.io.netcdf import netcdf_file as Dataset
    netcdf4 = False
#from scipy.io.netcdf import netcdf_file as Dataset
import skimsimulator.const as const
import numpy
import sys, os
import time as ti
import logging
logger = logging.getLogger(__name__)

def read_params(params_file):
    """ Read parameters from parameters file and store it in p.\n
    This program is not used in skim_simulator."""
    import imp
    f = open(params_file)
    #global p
    p = imp.load_source('p', '',f)
    f.close()
    return p

def write_params(params, pfile):
    """ Write parameters that have been selected to run skim_simulator. """
    f = open(pfile, 'w')
    for key in dir(params):
        if not key[0: 2] == '__':
            f.write('{} = {}\n'.format(key, params.__dict__[key]))
    f.close()

def read_coordinates(file,  nlon, nlat, twoD=True ):
    ''' General routine to read coordinates in a netcdf file. \n
    Inputs are file name, longitude name, latitude name. \n
    Outputs are longitudes and latitudes (2D arrays).'''

## - Open Netcdf file
    try:
        fid = Dataset(file, 'r')
    except IOError:
        print('There was an error opening the file {}'.format(file))
        sys.exit()
    #fid = Dataset(file, 'r')

## - Read 1d or 2d coordinates
    try:
        vartmp = fid.variables[nlat]
    except:
        sys.exit('Coordinates {} not found in file {}'.format(nlat, file))

    try:
        vartmp = fid.variables[nlon]
    except:
        sys.exit('Coordinates {} not found in file {}'.format(nlon, file))
        #sys.exit()
    if  len(vartmp.shape) == 1:
        lon_tmp = numpy.array(fid.variables[nlon][:]).squeeze()
        lat_tmp = numpy.array(fid.variables[nlat][:]).squeeze()
        if twoD:
            lon, lat = numpy.meshgrid(lon_tmp, lat_tmp)
        else:
            lon = lon_tmp
            lat = lat_tmp
    elif len(vartmp.shape) == 2:
        lon = numpy.array(fid.variables[nlon][:,:]).squeeze()
        lat = numpy.array(fid.variables[nlat][:,:]).squeeze()
        if not twoD:
            lon = lon[0,:]
            lat = lat[:,0]
    else:
        print('unknown dimension for lon and lat')
    fid.close()
    return lon, lat




def read_var(file, var, index=None, time=0, depth=0, model_nan=None):
    ''' General routine to read variables in a netcdf file. \n
    Inputs are file name, variable name, index=index to read part
    of the variables, time=time to read a specific time, depth=depth to read a
    specific depth, model_nan=nan value '''

## - Open Netcdf file
    try:
    #fid = nc.netcdf_file(file, 'r')
        fid = Dataset(file, 'r')
    except IOError:
        print('There was an error opening the file {}'.format(file))
        sys.exit()
    #fid = Dataset(file, 'r')

## - Check dimension of variable
    try : vartmp = fid.variables[var]
    except:
        sys.exit('Variable {} not found in file {}'.format(var, file))
## - Read variable
    if index is None:
        if len(vartmp.shape) == 1:
            T = numpy.array(fid.variables[var][:]).squeeze()
        elif len(vartmp.shape) == 2 :
            T = numpy.array(fid.variables[var][:, :]).squeeze()
        elif len(vartmp.shape) == 3 :
            if time is None:
                T = numpy.array(fid.variables[var][:, :, :]).squeeze()
            else:
                T = numpy.array(fid.variables[var][time, :, :]).squeeze()
        elif len(vartmp.shape) == 4 :
            if time is None:
                if depth is None:
                    T = numpy.array(fid.variables[var][:, :, :, :]).squeeze()
                else:
                    T = numpy.array(fid.variables[var][:, depth,
                                    :, :]).squeeze()
            elif depth is None:
                T = numpy.array(fid.variables[var][time, :, :, :]).squeeze()
            else:
                T = numpy.array(fid.variables[var][time, depth, :, :]).squeeze()
    else :
            if len(vartmp.shape) == 1:
                Ttmp = numpy.array(fid.variables[var][:]).squeeze()
                T = Ttmp[index]
            elif len(vartmp.shape) == 2 :
                Ttmp = numpy.array(fid.variables[var][:, :]).squeeze()
                T = Ttmp[index]
            elif len(vartmp.shape) == 3 :
                if time is None:
                    U = numpy.array(fid.variables[var][:, :, :]).squeeze()
                    T = U[:, index]
                else:
                    U = numpy.array(fid.variables[var][time, :, :]).squeeze()
                    T = U[index]
            elif len(vartmp.shape) == 4 :
                if time is None:
                    if depth is None:
                        U = numpy.array(fid.variables[var][:, :, :, :]).squeeze()
                        T = U[:, :, index]
                    else:
                        U = numpy.array(fid.variables[var][:, depth, :, :]).squeeze()
                        T = U[:, index]
                elif depth is None:
                    U = numpy.array(fid.variables[var][time, :, :, :]).squeeze()
                    T = U[:, index]
                else:
                    U = numpy.array(fid.variables[var][time,depth,:,:]).squeeze()
                    T = U[index]

            #print 'Unsupported number of dimensions'
            #sys.exit()
    fid.close()
## - Mask value that are NaN
    if not model_nan is None:
        T[numpy.where(T == model_nan)] = numpy.nan
    return T

class Sat_SKIM():
    ''' Sat_SKIM class: to read and write data that has been
    created by SKIM simulator '''
    def __init__(self,
                file=None,
                lon=None,
                lat=None,
                lon_nadir=None,
                lat_nadir=None,
                time=None,
                cycle=None,
                al_cycle=None,
                x_al=None,
                timeshift=None):
        self.file = file
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
## - Open Netcdf file in write mode
        if netcdf4:
          fid = Dataset(self.file, 'w', format='NETCDF4_CLASSIC')
        else:
          fid = Dataset(self.file, 'w')
        ## - Create Global attribute
        fid.title = 'SKIM swath simulated by SKIM simulator'
        fid.keywords = 'check keywords'  # Check keywords
        fid.Conventions = "CF-1.6"
        fid.summary = 'SKIM data produced'
        fid.description = "SKIM fixed swath"
        fid.Metadata_Conventions = "Unidata Dataset Discovery v1.0"
        fid.history = 'Grid File created by skimsimulator version ' + version
        fid.processing_level = 'L2'
        fid.standard_name_vocabulary = "CF-1.6"
        fid.creator_name = "Lucile Gaultier and Clement Ubelmann"
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
        ## - Create dimensions
        #if (not os.path.isfile(self.file)):
        dimsample = 'sample'
        nsample = numpy.shape(self.lon[-1])[0]
        fid.createDimension(dimsample, nsample)
        #fid.createDimension('time_nadir', numpy.shape(self.lon)[0])
        dimcycle = 'cycle'
        fid.createDimension(dimcycle, 1)
        nbeam = len(self.list_pos)
        dimnbeam = 'nbeam'
        fid.createDimension(dimnbeam, nbeam)
## - Create and write Variables
        vtime = fid.createVariable('time', 'f', (dimsample, dimnbeam))
        vtime.axis = "T"
        vtime.units = "seconds since the beginning of the sampling"
        vtime.long_name = "Time"
        vtime.standard_name = "time"
        vtime.calendar = "gregorian"
        vtime_nadir = fid.createVariable('time_nadir', 'f', (dimsample,))
        vtime_nadir.axis = "T"
        vtime_nadir.units = "seconds since the beginning of the sampling"
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
        vlat = fid.createVariable('lat', 'f4', (dimsample,dimnbeam))
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
        vx_ac = fid.createVariable('x_ac', 'f4', (dimsample,dimnbeam))
        vx_ac.units = "km"
        vx_ac.long_name = "Across track distance from the nadir"
        vx_al_nadir = fid.createVariable('x_al_nadir', 'f4', (dimsample,))
        vx_al_nadir.units = "km"
        vx_al_nadir.long_name = "Nadir along track distance from the"\
                                "beginning of the pass"
        vangle = fid.createVariable('angle', 'f4', (dimsample, dimnbeam))
        vangle.units = "rad"
        vangle.long_name = "Angle of the beam refered to the across track"\
                           " direction"
        vrangle = fid.createVariable('radial_angle', 'f4', (dimsample, dimnbeam))
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
        vcycle = fid.createVariable('cycle', 'f4', (dimcycle,))
        valcycle = fid.createVariable('al_cycle', 'f4', (dimcycle,))
        vtimeshift = fid.createVariable('timeshift', 'f4', (dimcycle,))
        vcycle[:] = self.cycle
        vcycle.units = "second"
        vcycle.long_name = "seconds during a cycle"
        valcycle[:] = self.al_cycle
        valcycle.units = "km"
        valcycle.long_name = " Distance travelled during the pass"
        vtimeshift[:] = self.timeshift
        vtimeshift.units = "day"
        vtimeshift.long_name = "Shift time to match model time"
        vlistpos = fid.createVariable('beam_position', 'f4', (dimnbeam,))
        vlistpos[:] = self.list_pos
        vlistpos.units = ""
        vlistpos.long_name = "Beam position"
        vlistangle = fid.createVariable('beam_angle', 'f4', (dimnbeam,))
        vlistangle[:] = self.list_angle
        vlistangle.units = ""
        vlistangle.long_name = "Beam angle"
        vincl = fid.createVariable('inclination', 'f4', (dimsample,))
        vincl.units = "rad"
        vincl.long_name = "Track inclination at nadir"
        #vincl[:] = numpy.array(self.angle)[0,:nsample]
        vincl[:] = self.angle[:nsample]
        fid.close()
        return None

    def write_data(self, p,  **kwargs):
        '''Write SKIM data in output file file_output
        Dimensions are x_al (along track distance), x_ac (across
        track distance). \n
        Variables are longitude, latitude, index (file number),
        error-free SSH (SSH interpolated from the model), selected
        errors (karin, wet_tropo, roll, baseline_dilation, phase,
        timing) and SSH with errors. \n
        '''
## - Open netcdf file in write mode
        if netcdf4:
          fid = Dataset(self.file, 'w', format='NETCDF4_CLASSIC')
        else:
          fid = Dataset(self.file, 'w')
        fid.description = "Ouptut from SKIM simulator"
        try:
            fid.corresponding_grid=self.gridfile
        except:
            pass
        fid.title = 'SKIM-like data simulated by SKIM simulator'
        fid.keywords = 'SKIM, altimetry, SSH, satellite, remote sensing'
        fid.Conventions = "CF-1.6"
        fid.summary = 'SKIM grid data produced'
        fid.description = "SKIM fixed grid"
        fid.Metadata_Conventions = "Unidata Dataset Discovery v1.0"
        fid.history = 'Grid File created by skimsimulator version ' + version
        fid.processing_level = 'L2'
        fid.standard_name_vocabulary = "CF-1.6"
        fid.creator_name = "Lucile Gaultier"
        fid.creator_email = "lucile.gaultier@gmail.com"
        fid.publisher_url = ""
        fid.time_coverage_start = self.time[0][0]
        # p.date0+"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.time_coverage_end = self.time[-1][-1]
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
        #fid.createDimension('time_nadir', numpy.shape(self.lon)[0])
        dimcycle = 'cycle'
        fid.createDimension(dimcycle, 1)
        nbeam = len(p.list_pos)
        dimnbeam = 'nbeam'
        fid.createDimension(dimnbeam, nbeam)
## - Create and write Variables
        vtime = fid.createVariable('time', 'f', (dimsample, dimnbeam))
        vtime.axis = "T"
        vtime.units = "seconds since the beginning of the sampling"
        vtime.long_name = "Time"
        vtime.standard_name = "time"
        vtime.calendar = "gregorian"
        vtime_nadir = fid.createVariable('time_nadir', 'f', (dimsample,))
        vtime_nadir.axis = "T"
        vtime_nadir.units = "seconds since the beginning of the sampling"
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
        vlat = fid.createVariable('lat', 'f4', (dimsample,dimnbeam))
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

        longname = {"instr": "Instrumental error",
                  "ur_model": "Radial velocity interpolated from model",
                  "u_model": "Zonal velocity interpolated from model",
                  "v_model": "Meridional velocity interpolated from model",
                  "ur_obs": "Observed radial velocity (Ur_model+errors)",
                  "index": "Equivalent model output number in list of file",
                  "ur_uss": "Stokes drift radial velocity bias",
                  "uss_err": "Stokes drift radial velocity bias_corrected",
                  "nadir_err": "Nadir error",
                  "std_uss": "Standard deviation of uss on a {} km"\
                             "radius".format(p.footprint_std),
                  "errdcos": "Weight for uss_err computation"}
        unit = {"instr": "m/s", "ur_model": "m/s", "ur_obs": "m/s",
                "index": " ", "ur_uss": "m/s", "uss_err": "m/s",
                "uss_err": "m/s", "nadir_err": "m/s", "u_model": "m/s",
                "v_model": "m/s", "std_uss": "m/s", "errdcos": "km",
                }
        list_nadir = {"instr", "nadir_err"}
        for key, value in kwargs.items():
            if value is not None:
                nvar_nadir = '{}_nadir'.format(key)
                var_nadir = fid.createVariable(nvar_nadir, 'f4', (dimsample,),
                                               fill_value=-1.36e9)
                nvar = '{}'.format(key)
                var = fid.createVariable(nvar, 'f4',
                                         (dimsample,dimnbeam),
                                         fill_value=-1.36e9)
                try:
                    var.units = unit[str(key)]
                    var_nadir.units = unit[str(key)]
                except:
                    var.units=''
                    var_nadir.units=''
                try:
                    var.long_name = longname[str(key)]
                    var_nadir.long_name = longname[str(key)]
                except:
                    var.long_name=str(key)
                    var_nadir.long_name=str(key)
                for i in range(len(value)):
                    if value[i].any():
                        value_tmp = value[i][:]
                        vmin = numpy.nanmin(value_tmp)
                        vmax = numpy.nanmax(value_tmp)
                        mask = numpy.isnan(value_tmp)
                        value_tmp[numpy.where(mask)] = -1.36e9
                        mask_ind = numpy.where(value_tmp == 0)
                        value_tmp[mask_ind] = -1.36e9
                        if i == 0:
                            if key in list_nadir:
                                var_nadir[:] = value_tmp
                            else:
                                continue
                        else:
                            try:
                                var[:, i - 1] = value_tmp
                            except:
                                import pdb ; pdb.set_trace()

                # try:    var.missing_value = p.model_nan
                # except: var.missing_value = 0.
                # fid.setncattr('missing_value','-9999.f')
        ###############TODO set range values
        fid.close()
        return None

    def load_swath(self, p, **kwargs):
        '''Load swath variables stored in Satellite grid file sgridfile. \n
        (longitude, latitude, number of days in a cycle, crossed distance
        during a cycle, time, along track and across track position).'''

## - Open Netcdf file
        try :
            fid = Dataset(self.file, 'r')
        except IOError:
            print('There was an error opening the file '+self.file)
            sys.exit(1)
        #fid = Dataset(self.file, 'r')
        time = []; lon = []; lat = []; cycle = []; x_al = []
        listvar={'time':time, 'lon': lon, 'lat': lat,}
        self.lon = []
        self.lat = []
        self.time = []
## - Read variables in listvar and return them
        for stringvar in listvar:
            var = fid.variables['{}{}'.format(stringvar, '_nadir')]
            listvar[stringvar].append(numpy.array(var[:]).squeeze())
            var = fid.variables[stringvar]
            for i in range(len(p.list_pos)):
                listvar[stringvar].append(numpy.array(var[:, i]).squeeze())
            #exec('self.'+stringvar+' = listvar[stringvar]')
            setattr(self, stringvar, listvar[stringvar])
## - Read variables in arguments
        for key, value in kwargs.items():
            var = fid.variables[key]
            value = numpy.array(fid.variables[key][:]).squeeze()
            #value[value == var.fill_value] = numpy.nan
            setattr(self, key, value)
        self.pos = numpy.array(fid.variables['beam_position'][:])
        if numpy.all(numpy.abs(self.pos - numpy.array(p.list_pos, dtype = numpy.float32)) > 0.0001):
            logger.error('List of beam positions has changed, reprocess the grids')
            sys.exit(1)
        self.angle = numpy.array(fid.variables['beam_angle'][:])
        if numpy.all(numpy.abs(self.angle - numpy.array(p.list_angle, dtype = numpy.float32)) > 0.0001):
            logger.error('List of beam angles has changed, reprocess the grids')
            sys.exit(1)
        self.radial_angle = numpy.array(fid.variables['radial_angle'][:])
        try:
            self.corresponding_grid=fid.corresponding_grid
        except:
            pass
        try:
            self.incl=numpy.array(fid.variables['inclination'][:]).squeeze()
        except:
            print('inclination variable not found')
        fid.close()
        return None


class NEMO():
    '''Class to read NEMO data \n
    USAGE is NEMO(file=name of file ,var= variable name,
    lon=longitude name, lat=latitude name, depth= depth name,
    time=time name).\n
    Argument file is mandatory, other arguments have default
    values var='sossheig', lon='nav_lon', lat='nav_lat', depth='depth',
    time='time. \n'''
    def __init__(self,
                file=None,
                varu='vomecrtx',
                lonu='nav_lonu',
                latu='nav_latu',
                varv='vozocrty',
                lonv='nav_lonv',
                latv='nav_latv',
                time='time',
                depth='depth',
                ):
        self.nvar = varu
        self.nlon = lonu
        self.nlat = latu
        self.nvar = varv
        self.nlon = lonv
        self.nlat = latv
        self.ntime = time
        self.nfile = file
        self.ndepth=depth
        try:
            self.model_nan = p.model_nan
        except:
            self.model_nan = 0.
            p.model_nan = 0.

    def read_var(self, index=None):
        '''Read variables from NEMO file \n
        Argument is index=index to load part of the variable.'''
        try:
            vel_factor = p.vel_factor
        except:
            vel_factor = 1.
            p.vel_factor = vel_factor
        self.vvaru = read_var(self.nfile, self.nvaru, index=index, time=0,
                              depth=0, model_nan=self.model_nan)*vel_factor
        self.vvarv = read_var(self.nfile, self.nvarv, index=index, time=0,
                              depth=0, model_nan=self.model_nan)*vel_factor
        return None

    def read_coordinates(self, index=None):
        '''Read coordinates from NEMO file \n
        Argument is index=index to load part of the variable.'''
        if p.grid == 'regular':
            lonu, latu = read_coordinates(self.nfile, self.nlonu, self.nlatu,
                                        twoD=False)
            lonv, latv = read_coordinates(self.nfile, self.nlonv, self.nlatv,
                                        twoD=False)
        else:
            lonu, latu = read_coordinates(self.nfile, self.nlonu, self.nlatu)
            lonv, latv = read_coordinates(self.nfile, self.nlonv, self.nlatv)

        self.vlatu = latu
        self.vlonu = (lonu + 360) % 360
        self.vlatv = latv
        self.vlonv = (lonv + 360) % 360
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from NEMO file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if (numpy.min(self.vlonu) < 1.) and (numpy.max(self.vlonu) > 359.):
            self.vlonu[numpy.where(self.vlonu > 180.)] = self.vlonu[numpy.where(self.vlonu > 180.)] - 360
            lon1 = (numpy.min(self.vlon) + 360) % 360
            lon2 = (numpy.max(self.vlon) + 360) % 360
            if lon1 == lon2:
              lon1 = 0
              lon2 = 360
        else:
            lon1 = numpy.min(self.vlon)
            lon2 = numpy.max(self.vlon)
        return [lon1, lon2, numpy.min(self.vlat), numpy.max(self.vlat)]

class ROMS():
    '''Class to read ROMS data \n
    USAGE is ROMS(file=name of file ,var= variable name,
    lon=longitude name, lat=latitude name, depth= depth name,
    time=time name).\n
    Argument file is mandatory, other arguments have default
    values var='rho', lon='x_rho', lat='y_rho', depth='depth',
    time='time. \n
    Variable units is specified in params file and default value
    is True (coordinates in degree). \n
    If units is False (coordinates in km), specify left
    low corner of the domain (lon0, lat0) in params file.'''
    def __init__(self,
                file=None,
                varu = 'u',
                varv = 'v',
                depth = 'depth',
                time = 'time',
                lonu = 'lon_u',
                latu = 'lat_u',
                lonv = 'lon_v',
                latv = 'lat_v',
                ):
        self.nvaru = varu
        self.nlonu = lonu
        self.nlatu = latu
        self.nvarv = varv
        self.nlonv = lonv
        self.nlatv = latv
        self.ntime = time
        self.nfile = file
        self.ndepth = depth
        try:
            self.model_nan = p.model_nan
        except:
            self.model_nan = 0.
            p.model_nan = 0.

    def read_var(self, index=None):
        '''Read variables from ROMS file\n
        Argument is index=index to load part of the variable.'''
        try:
            SSH_factor = p.SSH_factor
        except:
            vel_factor = 1.
            p.vel_factor = 1.
        self.vvaru = read_var(self.nfile, self.nvaru, index=index, time=0,
                              depth=0, model_nan=self.model_nan) *vel_factor
        self.vvarv = read_var(self.nfile, self.nvarv, index=index, time=0,
                              depth=0, model_nan=self.model_nan) * vel_factor
        return None


    def read_coordinates(self, index=None):
        '''Read coordinates from ROMS file \n
        Argument is index=index to load part of the variable.'''
        if p.grid == 'regular':
            lonu, latu = read_coordinates(self.nfile, self.nlonu, self.nlatu,
                                          twoD=False)
            lonv, latv = read_coordinates(self.nfile, self.nlonv, self.nlatv,
                                          twoD=False)
        else:
            lonu, latu = read_coordinates(self.nfile, self.nlonu, self.nlatu)
            lonv, latv = read_coordinates(self.nfile, self.nlonv, self.nlatv)
        self.vlatu = latu
        self.vlonu = (lonu + 360) % 360
        self.vlatv = latv
        self.vlonv = (lonv + 360) % 360
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from ROMS file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if (numpy.min(self.vlonu) < 1.) and (numpy.max(self.vlon) > 359.):
            self.vlonu[numpy.where(self.vlonu > 180.)] = self.vlonu[numpy.where(self.vlonu > 180.)] - 360
            lon1=(numpy.min(self.vlonu) + 360) % 360
            lon2=(numpy.max(self.vlonu) + 360) % 360
        else:
            lon1=numpy.min(self.vlonu)
            lon2=numpy.max(self.vlonu)
        return [lon1, lon2, numpy.min(self.vlatu), numpy.max(self.vlatu)]



class NETCDF_MODEL():
    '''Class to read any netcdf data.\n
    USAGE is NETCDF_MODEL(file=name of file ,var= variable name,
    lon=variable longitude, lat=variable latitude, units=).\n
    Argument file is mandatory, arguments var, lon, lat
    are specified in params file. \n
    '''
    def __init__(self,
                file=None,
                varu='u',
                lonu='lon_u',
                latu='lat_u',
                varv='v',
                lonv='lon_v',
                latv='lat_v',
                depth=0,
                time=0,
                ):
        self.nvaru = varu
        self.nlonu = lonu
        self.nlatu = latu
        self.nvarv = varv
        self.nlonv = lonv
        self.nlatv = latv
        self.nfile = file
        self.depth = depth
        self.time = time
        try:
            self.model_nan = p.model_nan
        except:
            self.model_nan = 0.
            p.model_nan = 0.

    def read_var(self, index=None):
        '''Read variables from netcdf file \n
        Argument is index=index to load part of the variable.'''
        try:
            vel_factor = p.vel_factor
        except:
            vel_factor=1.
            p.vel_factor=1.
        self.vvarv = read_var(self.nfile, self.nvarv, index=index,
                              time=self.time, depth=self.depth,
                              model_nan=self.model_nan) * vel_factor
        self.vvaru = read_var(self.nfile, self.nvaru, index=index,
                              time=self.time, depth=self.depth,
                              model_nan=self.model_nan) * vel_factor
        #self.vvar[numpy.where(numpy.isnan(self.vvar))]=0
        return None

    def read_coordinates(self, index=None):
        '''Read coordinates from netcdf file \n
        Argument is index=index to load part of the variable.'''
        if p.grid=='regular':
            lonu, latu = read_coordinates(self.nfile, self.nlonu, self.nlatu,
                                          twoD=False)
            lonv, latv = read_coordinates(self.nfile, self.nlonv, self.nlatv,
                                          twoD=False)
        else:
            lonu, latu = read_coordinates(self.nfile, self.nlonu, self.nlatu)
            lonv, latv = read_coordinates(self.nfile, self.nlonv, self.nlatv)
        self.vlatu = latu
        self.vlonu = (lonu + 360) % 360
        self.vlatv = latv
        self.vlonv = (lonv + 360) % 360
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from netcdf file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if (numpy.min(self.vlonu) < 1.) and (numpy.max(self.vlonu) > 359.):
            self.vlonu[numpy.where(self.vlonu > 180.)]=self.vlonu[numpy.where(self.vlonu > 180.)] - 360
            lon1=(numpy.min(self.vlonu) + 360) % 360
            lon2=(numpy.max(self.vlonu) + 360) % 360
        else:
            lon1=numpy.min(self.vlonu)
            lon2=numpy.max(self.vlonu)
        return [lon1, lon2, numpy.min(self.vlatu), numpy.max(self.vlatu)]

class WW3():
    '''Class to read ww3 netcdf data.\n
    USAGE is NETCDF_MODEL(file=name of file ,var= variable name,
    lon=variable longitude, lat=variable latitude, units=).\n
    Argument file is mandatory, arguments var, lon, lat
    are specified in params file. \n
    '''
    def __init__(self, p,
                file=None,
                varu='ucur',
                lonu='longitude',
                latu='latitude',
                varv='vcur',
                lonv='longitude',
                latv='latitude',
                depth=0,
                time=0,
                ):
        self.nvaru = varu
        self.nlonu = lonu
        self.nlatu = latu
        self.nvarv = varv
        self.nlonv = lonv
        self.nlatv = latv
        self.nfile = file
        self.depth = depth
        self.time = time
        self.p = p
        try:
            self.model_nan = p.model_nan
        except:
            self.model_nan = 0.
            p.model_nan = 0.

    def read_var(self, index=None):
        '''Read variables from netcdf file \n
        Argument is index=index to load part of the variable.'''
        try:
            vel_factor = self.p.vel_factor
        except:
            vel_factor=1.
            self.p.vel_factor=1.
        self.vvarv = read_var(self.nfile, self.nvarv, index=index,
                              time=self.time, depth=self.depth,
                              model_nan=self.model_nan) * vel_factor
        self.vvaru = read_var(self.nfile, self.nvaru, index=index,
                              time=self.time, depth=self.depth,
                              model_nan=self.model_nan) * vel_factor
        #self.vvar[numpy.where(numpy.isnan(self.vvar))]=0
        return None

    def read_coordinates(self, index=None):
        '''Read coordinates from netcdf file \n
        Argument is index=index to load part of the variable.'''
        if self.p.grid=='regular':
            lonu, latu = read_coordinates(self.nfile, self.nlonu, self.nlatu,
                                          twoD=False)
            lonv, latv = read_coordinates(self.nfile, self.nlonv, self.nlatv,
                                          twoD=False)
        else:
            lonu, latu = read_coordinates(self.nfile, self.nlonu, self.nlatu)
            lonv, latv = read_coordinates(self.nfile, self.nlonv, self.nlatv)
        self.vlatu = latu
        self.vlonu = (lonu + 360) % 360
        self.vlatv = latv
        self.vlonv = (lonv + 360) % 360
        return None

    def calc_box(self):
        '''Calculate subdomain coordinates from netcdf file
        Return minimum, maximum longitude and minimum, maximum latitude'''
        self.read_coordinates()
        if (numpy.min(self.vlonu) < 1.) and (numpy.max(self.vlonu) > 359.):
            self.vlonu[numpy.where(self.vlonu > 180.)]=self.vlonu[numpy.where(self.vlonu > 180.)] - 360
            lon1=(numpy.min(self.vlonu) + 360) % 360
            lon2=(numpy.max(self.vlonu) + 360) % 360
        else:
            lon1=numpy.min(self.vlonu)
            lon2=numpy.max(self.vlonu)
        return [lon1, lon2, numpy.min(self.vlatu), numpy.max(self.vlatu)]

