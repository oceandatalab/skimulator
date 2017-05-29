# -----------------------#
# Files and directories
# -----------------------#
## -- Get the user home directory
from os.path import expanduser
import os
import math
home = expanduser("~") + '/src/'
# ------ Directory that contains orbit file:
dir_setup = os.path.join(home, 'skimsimulator', 'data')
# ------ Directory that contains your own inputs:
indatadir = os.path.join(home, 'skimsimulator', 'example', 'input_fields')
indatadir = '/mnt/data/model/ww3_fram/' #netcdf3/'
# ------ Directory that contains your outputs:
outdatadir = os.path.join(home, 'skimsimulator', 'example', 'skim_output')
# ------ Orbit file:
#filesat = os.path.join(dir_setup,'orbs1a.txt')
filesat = os.path.join(dir_setup,'orbits1_ifremer')
# , dir_setup+os.sep+'orbjason.txt', dir_setup+os.sep+'orbaltika.txt' ]
# ------ Name of the configuration (to build output files names) 
config="WW3_FRAM_8b60az"

# -----------------------# 
# SKIM swath parameters 
# -----------------------# 
# ------ Satellite grid file root name:
# 	 (Final file name is root_name_[numberofpass].nc)
filesgrid = os.path.join(outdatadir, '{}_grid'.format(config))
# ------ Force the computation of the satellite grid:
makesgrid = False
# ------ Give a subdomain if only part of the model is needed:
#	 (modelbox=[lon_min, lon_max, lat_min, lat_max])
# 	 (If modelbox is None, the whole domain of the model is considered)
modelbox =  [334,24,72,86]
#------- Rotation speed of the antenna (in tr/min)
#rotation_speed = 3.774  # * 180
rotation_speed = 3.396739
#------- List of position of beams:
list_pos = (0, 72*math.pi/180., 144*math.pi/180., 216*math.pi / 180.,
            288*math.pi/180., 0, math.pi)
#------- List of angle of beams in degrees:
list_angle = (12, 12, 12, 12, 12, 6, 6)
#------- List of timeshift as regard to nadir for 12 degree beams:
list_shift = (1, 2, 4, 5, 7, 3, 6)
# ------ Shift longitude of the orbit file if no pass is in the domain 
#        (in degree): Default value is None (no shift)
shift_lon = 0
# ------ Shift time of the satellite pass (in day):
#        Default value is None (no shift)
shift_time = None

# -----------------------#
# Model input parameters
# -----------------------#
# ------ List of model files:
#	 (The first file contains the grid and is not considered as model data)
#        To generate the noise alone, file_input=None and specify region 
#        in modelbox
file_input = os.path.join(indatadir, 'list_of_file.txt')
# ------ Type of model data: 
#	 (Optional, default is NETCDF_MODEL and reads netcdf3 and netcdf4 files)
#	 (Other options are ROMS, NEMO and CLS to read Nemo, roms or CLS)
model = 'WW3'
# ------ Type of grid: 
#        'regular' or 'irregular', if 'regular' only 1d coordinates 
#        are extracted from model       
grid = 'regular'
# ------ Specify velocities variable:
varu = 'ucur'
varv = 'vcur'
# ------ Specify factor to convert velocity values in m/s:
vel_factor = 1.
# ------ Specify longitude variable:
lonu = 'longitude'
lonv = 'longitude'
# ------ Specify latitude variable:
latu = 'latitude'
latv = 'latitude'
# ------ Specify number of time in file:
dim_time = (71,)
# ------ Time step between two model outputs (in days):
timestep = 1/24.
# ------ Number of outputs to consider:
#        (timestep*nstep=total number of days)
nstep = 70.
# ------ Not a number value:
model_nan = -32767.

# -----------------------# 
# SKIM output files  
# -----------------------# 
# ------ Output file root name:
#	 (Final file name is root_name_c[cycle].nc
file_output = os.path.join(outdatadir, config)
# ------ Interpolation of the SSH from the model (if grid is irregular and 
#         pyresample is not installed:
#        (either 'linear' or 'nearest', use 'nearest' for large region
#        as it is faster and use less memory.)
interpolation = 'linear'
# -----------------------# 
# SKIM error parameters 
# -----------------------# 
# ------ File containing random coefficients to compute and save 
#	 random error coefficients so that runs are reproducible:
#        If file_coeff is specified and does not exist, file is created
#	 If you don't want runs to be reproducible, file_coeff is set to None
file_coeff = None  # outdatadir+os.sep+'Random_coeff.nc'
# ------ Number of random realisations for instrumental and geophysical error 
#        (recommended ncomp=2000), ncomp1d is used for 1D spectrum, and ncomp2d
#        is used for 2D spectrum (wet troposphere computation):
ncomp1d = 3000
ncomp2d = 2000
# ------- Instrument white noise error
instr = True
# ------- Instrument white noise rms 
# rms_instr = [10 * 10**(-2), 10 * 10**(-2), 10 * 10**(-2), 10 * 10**(-2),
# 20 * 10 ** (-2)]
rms_instr = [os.path.join(dir_setup, 'instrumentnoise_12.dat'),
             os.path.join(dir_setup, 'instrumentnoise_12.dat'),
             os.path.join(dir_setup, 'instrumentnoise_12.dat'),
             os.path.join(dir_setup, 'instrumentnoise_12.dat'),
             os.path.join(dir_setup, 'instrumentnoise_12.dat'),
             os.path.join(dir_setup, 'instrumentnoise_06.dat'),
             os.path.join(dir_setup, 'instrumentnoise_06.dat')]
# ------- Stoke drift velocity [beam 12, beam 6]
uss = True
input_uss = os.path.join(indatadir, 'list_file_uss.txt')
G = [50, 50, 50, 50, 50, 50, 50]
bias_std = 0.09
errdcos = None
#[25.2006/20, 25.2747/20, 25.4763/20, 25.4271/20, 19.9728/20]
footprint_std = 0 #400
formula = False

## -- Geophysical error
## ----------------------
# ------ Wet tropo error (True to compute it):
wet_tropo = True
# ------ Beam print size (in km):
#        Gaussian footprint of sigma km
sigma = 8.
