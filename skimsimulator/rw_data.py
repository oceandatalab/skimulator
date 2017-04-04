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
try: import params as p
except:
    print('params.py not found')
    sys.exit()
#from scipy.io import netcdf as nc
#from scipy.io.netcdf import netcdf as nc

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
        fid.time_coverage_start = self.time[0]
        # p.date0+"YYYY-MM-DDThh:mmZ"  #tim0 converted to format
        fid.time_coverage_end = self.time[-1]
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
        ## - Create dimensions
        #if (not os.path.isfile(self.file)):
        fid.createDimension('time', numpy.shape(self.lon[0])[0])
        #fid.createDimension('time_nadir', numpy.shape(self.lon)[0])
        fid.createDimension('cycle', 1)

## - Create and write Variables
        n12beam = len(p.list_pos_12)
        n6beam = len(p.list_pos_6)
        for i in range(n12beam + 1 + n6beam):
            if i == 0:
                ntime = 'time_nadir'
                nlon = 'lon_nadir'
                nlat = 'lat_nadir'

            elif i <= n12beam:
                ntime = 'time_12_{}'.format(i - 1)
                nlon = 'lon_12_{}'.format(i - 1)
                nlat = 'lat_12_{}'.format(i - 1)
            else:
                ntime = 'time_6_{}'.format(i - n12beam - 1)
                nlon = 'lon_6_{}'.format(i - n12beam - 1)
                nlat = 'lat_6_{}'.format(i - n12beam - 1)
            vtime = fid.createVariable(ntime, 'f', ('time',))
            vtime[:] = self.time[i][:]
            vtime.axis = "T"
            vtime.units = "days since the beginning of the sampling"
            vtime.long_name = "Time"
            vtime.standard_name = "time"
            vtime.calendar = "gregorian"
            vlon = fid.createVariable(nlon, 'f4', ('time',))
            vlon[:] = self.lon[i][:]
            vlon.axis = "X"
            vlon.long_name = "Longitude"
            vlon.standard_name = "longitude"
            vlon.units = "degrees_east"
            vlat = fid.createVariable(nlat, 'f4', ('time',))
            vlat[:] = self.lat[i][:]
            vlat.axis = "Y"
            vlat.long_name = "Latitude"
            vlat.standard_name = "latitude"
            vlat.units = "degrees_north"
        vcycle = fid.createVariable('cycle', 'f4', ('cycle',))
        valcycle = fid.createVariable('al_cycle', 'f4', ('cycle',))
        vtimeshift = fid.createVariable('timeshift', 'f4', ('cycle',))
        vx_al = fid.createVariable('x_al', 'f4', ('time',))
        vcycle[:] = self.cycle
        vcycle.units = "days during a cycle"
        vcycle.long_name = "Cycle"
        valcycle[:] = self.al_cycle
        valcycle.units = "km"
        valcycle.long_name = " Distance travelled during the pass"
        vtimeshift[:] = self.timeshift
        vtimeshift.units = "day"
        vtimeshift.long_name = "Shift time to match model time"
        vx_al[:] = self.x_al
        vx_al.units = "km"
        vx_al.long_name = "Along track distance from the beginning of the pass"
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
        fid.keywords_vocabulary = "NASA"
        fid.references = ""
        # fid.cycle = "{0:d}".format(int(self.al_cycle))
        fid.createDimension('time', numpy.shape(self.lon[0])[0])
        #fid.createDimension('time_nadir', numpy.shape(self.lon)[0])
        fid.createDimension('cycle', 1)

## - Create and write Variables
        n12beam = len(p.list_pos_12)
        n6beam = len(p.list_pos_6)
        for i in range(n12beam + 1 + n6beam):
            if i == 0:
                ntime = 'time_nadir'
                nlon = 'lon_nadir'
                nlat = 'lat_nadir'

            elif i <= n12beam:
                ntime = 'time_12_{}'.format(i - 1)
                nlon = 'lon_12_{}'.format(i - 1)
                nlat = 'lat_12_{}'.format(i - 1)
            else:
                ntime = 'time_6_{}'.format(i - n12beam - 1)
                nlon = 'lon_6_{}'.format(i - n12beam - 1)
                nlat = 'lat_6_{}'.format(i - n12beam - 1)
            vtime = fid.createVariable(ntime, 'f', ('time',))
            vtime[:] = self.time[i][:]
            vtime.axis = "T"
            vtime.units = "days since the beginning of the sampling"
            vtime.long_name = "Time"
            vtime.standard_name = "time"
            vtime.calendar = "gregorian"
            vlon = fid.createVariable(nlon, 'f4', ('time',))
            vlon[:] = self.lon[i][:]
            vlon.axis = "X"
            vlon.long_name = "Longitude"
            vlon.standard_name = "longitude"
            vlon.units = "degrees_east"
            vlat = fid.createVariable(nlat, 'f4', ('time',))
            vlat[:] = self.lat[i][:]
            vlat.axis = "Y"
            vlat.long_name = "Latitude"
            vlat.standard_name = "latitude"
            vlat.units = "degrees_north"
        longname = {"instr": "Instrumental error",
                  "ur_model":"Radial velocity interpolated from model",
                  "ur_obs":"Observed radial velocity (Ur_model+errors)",
                  "index":"Equivalent model output number in list of file",
                  "swh_err":"", "nadir_err":"Nadir error", "stoke_err": ""}
        unit = {"instr":"m/s", "ur_model":"m/s", "ur_obs":"m/s", "index":" ",
                "swh_err":"m/s", "nadir_err":"m/s", "stoke_err":"m/s",
                }
        for key, value in kwargs.items():
            #if not value is None:
            if value is not None:
                #for i in range(n12beam + 1 + n6beam):
                for i in range(len(value)):
                    if value[i].any():
                        if i == 0:
                            nvar = '{}_nadir'.format(key)
                        elif i <= n12beam:
                            nvar = '{}_12_{}'.format(key, i - 1)
                        else:
                            nvar = '{}_6_{}'.format(key, i - n12beam - 1)
                        var = fid.createVariable(nvar, 'f4', ('time',),
                                                 fill_value=-1.36e9)

                        value_tmp = value[i][:]
                        vmin=numpy.nanmin(value_tmp)
                        vmax=numpy.nanmax(value_tmp)
                        value_tmp[numpy.isnan(value_tmp)] = -1.36e9
                        value_tmp[value_tmp==0] = -1.36e9
                        var[:] = value_tmp
                        try:
                            var.units = unit[str(key)]
                        except:
                            var.units=''
                        try:
                            var.long_name = longname[str(key)]
                        except:
                            var.long_name=str(key)
                # try:    var.missing_value = p.model_nan
                # except: var.missing_value = 0.
                # fid.setncattr('missing_value','-9999.f')

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
            sys.exit()
        #fid = Dataset(self.file, 'r')
        time = []; lon = []; lat = []; cycle = []; x_al = []
        listvar={'time':time, 'lon': lon, 'lat': lat}
        self.lon = []
        self.lat = []
        self.time = []
        suffixe = ['_nadir']
        for i in range(len(p.list_pos_12)):
           suffixe.append('_12_{}'.format(i))
        for i in range(len(p.list_pos_6)):
           suffixe.append('_6_{}'.format(i))
## - Read variables in listvar and return them
        for stringvar in listvar:
            for suf in suffixe:
                var = fid.variables['{}{}'.format(stringvar, suf)]
                listvar[stringvar].append(numpy.array(var[:]).squeeze())
            exec('self.'+stringvar+' = listvar[stringvar]')

## - Read variables in arguments
        for key, value in kwargs.items():
            var = fid.variables[key]
            value = numpy.array(fid.variables[key][:]).squeeze()
            #value[value == var.fill_value] = numpy.nan
            exec('self.'+key+' = value')
        try:
            self.corresponding_grid=fid.corresponding_grid
        except:
            pass
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
    def __init__(self,
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

