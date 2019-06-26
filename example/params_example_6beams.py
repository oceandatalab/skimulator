# -----------------------#
# Files and directories
# -----------------------#
## -- Get the user home directory
from os.path import expanduser
import os
import math
home = expanduser("~")

# ------ Name of the configuration (to build output files names) 
#config="WW3_EQ_metop_2018_6a"
# 6 beams, 60 azimuths, 512 pulses and cycle length of 37/2 ms
#config="WW3_EQ_metop_2019_6a"
config="WW3_EQ_metop_2019_6b"
config="WW3_EQ_metop_2019_6c"
# 6 beams, ?? azimuths, 1024 pulses and cycle length of 37 ms
#config = "WW3_EQ_metop_2019_6b"
# ------ Directory that contains orbit file:
dir_setup = os.path.join(home, 'skimulator', 'data')
# ------ Directory that contains your own inputs:
indatadir = os.path.join(home, 'skimulator', 'example', 'input_fields')
# ------ Directory that contains your outputs:
outdatadir = os.path.join(home, 'skimulator', 'example', 'skim_output')
# ------ Orbit file:
#filesat = os.path.join(dir_setup,'orbs1a.txt')
filesat = os.path.join(dir_setup,'orbmetop_skim.txt')
# ------ Number of days in orbit (optional if specified in orbit file)
satcycle = 29
# ------ Satellite elevation (optional if specified in orbit file)
sat_elev = 817 * 10**3
# ------ Order of columns (lon, lat, time) in orbit file
# (default is (0, 1, 2) with order_orbit_col = None)
order_orbit_col = None
# , dir_setup+os.sep+'orbjason.txt', dir_setup+os.sep+'orbaltika.txt' ]
# ------ Number of processor for parallelisation
proc_number = 1
# ------ Deactivate printing of progress bar to avoid huge log
progress_bar = True

# -----------------------#
# SKIM swath parameters
# -----------------------#
# ------ Satellite grid file root name:
# 	 (Final file name is root_name_[numberofpass].nc)
filesgrid = os.path.join(outdatadir, '{}_grid'.format(config))
# ------ Force the computation of the satellite grid:
makesgrid = True
# ------ Give a subdomain if only part of the model is needed:
#	 (modelbox=[lon_min, lon_max, lat_min, lat_max])
# 	 (If modelbox is None, the whole domain of the model is considered)
modelbox = None # [329.,347., -8.,8.]
#------- Rotation speed of the antenna (in tr/min)
if '2019_6a' in config:
    rotation_speed = 3.3
    #rotation_speed = 6.3
elif '2019_6b' in config:
    rotation_speed = 6
elif '2019_6c' in config:
    rotation_speed = 6
else:
    rotation_speed = 9.26
# ------ Cycle duration
cycle = 0.0368
#------- List of position of beams:
if '2019_6a' in config:
    list_pos = (0, 90 * math.pi/180., 180*math.pi/180, 270*math.pi/180, 0)
elif '2019_6b' in config:
    list_pos = (0, 120 * math.pi/180, 240 * math.pi/180, 0, 180 * math.pi / 180)
elif '2019_6c' in config:
    list_pos = (0, 120 * math.pi/180, 240 * math.pi/180, 0, 180 * math.pi / 180)
else:
    list_pos = (0, 120*math.pi/180., 240*math.pi/180.,
            0, math.pi)
#------- List of angle of beams in degrees:
if '2019_6a' in config:
    list_angle = (12, 12, 12, 12, 6)
else:
    list_angle = (12, 12, 12, 6, 6)
#------- List of timeshift as regard to nadir for 12 degree beams:
if '2019_6a' in config:
    list_shift = (1, 2, 5, 3, 4)
elif '2019_6b' in config:
    list_shift = (1, 5, 2, 4, 3)
elif '2019_6c' in config:
    list_shift = (1, 4, 2, 5, 3)
else:
    list_shift = (5, 2, 3, 1, 4)
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
#	 (Other option is WW3)
model = 'WW3'
# ------ First time of the model
first_time = '2011-11-15T00:00:00Z'
# ------ Grid file name
file_grid_model = (os.path.join(indatadir, 'ww3.20111115_cur.nc'),)
# ------ Specify if there is a ice mask for high latitudes
#        (if true, mask is recomputed at each cycle)
ice_mask = False
# ------ Type of grid: 
#        'regular' or 'irregular', if 'regular' only 1d coordinates 
#        are extracted from model       
grid = 'regular'
# ------ Specify list of variable:
list_input_var = {'ucur': ['ucur', 'cur', 0], 'vcur': ['vcur', 'cur', 0],
                  'uuss': ['uuss', 'uss', 0], 'vuss': ['vuss', 'uss', 0],
                  'ice': ['ice', 'ice', 0], 'mssd': ['mssd', 'msd', 0],
                  'mssx': ['mssu', 'mss', 0], 'mssy':['mssc', 'mss', 0],
                  'ssh': ['wlv', 'wlv', 0], 'hs':['hs', 'hs', 0],
                  'uwnd': ['uwnd', 'wnd', 0], 'vwnd': ['vwnd', 'wnd', 0]}
# ------ Specify longitude variable:
lon = ('longitude',)
# ------ Specify latitude variable:
lat = ('latitude',)
# ------ Specify number of time in file:
dim_time = 24
# ------ Time step between two model outputs (in days):
timestep = 1/24.
# ------ Number of outputs to consider:
#        (timestep*nstep=total number of days)
nstep = 35*24
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
# ------ List of output variables:
list_output = ['ssh_obs', 'ur_true', 'instr',
               'radial_angle', 'uwb',
               'ssh_true', 'vindice', 'ur_obs', 'sigma0']

# -----------------------# 
# SKIM error parameters 
# -----------------------# 
# ------ File containing random coefficients to compute and save 
#	 random error coefficients so that runs are reproducible:
#        If file_coeff is specified and does not exist, file is created
#	 If you don't want runs to be reproducible, file_coeff is set to None
file_coeff = None  # outdatadir+os.sep+'Random_coeff.nc'
# Compute instrumental nadir noise:
nadir = True
# ------ Number of random realisations for instrumental and geophysical error 
#        (recommended ncomp=2000), ncomp1d is used for 1D spectrum, and ncomp2d
#        is used for 2D spectrum (wet troposphere computation):
ncomp1d = 3000
ncomp2d = 2000
# ------- Instrument white noise error
instr = True
# ------- Choice of instrument configuration
instr_configuration = 'A'
# ------- Coefficient SNR to retrieve instrumental noise from sigma, 
#         Recommanded value for 1024 pulses: 3e-2, for 512 pulses: 3sqrt(2)e-3
snr_coeff = 1.4142*6e-3
# ------- Attitude error
attitude = True
# ------- File which provide the AOCS error:
yaw_file = os.path.join(dir_setup, 'sample_req1.nc')

# ------- Wave bias
uwb = True


## -- Geophysical error
## ----------------------
# ------ Consider ice in sigma0 computation
ice = True
# ------ Rain error (True to compute it):
rain = True
# ------ Rain file containing scenarii (python file):
rain_file = os.path.join(dir_setup, 'rain_eq_atl.pyo')
# ------ Threshold to flag data:
rain_threshold = 0.15
wet_tropo = False

# -----------------------#
# L2C computation
# -----------------------#
# config name for L2d:
config_l2c = ''
# Length resolution to select neighbors (in km):
resol = 40
# Grid resolution for l2c (alongtrack, acrosstrack) grid (in km):
posting = 5
# Remove noisy data around nadir (in km):
ac_threshold = 20
# List of variables to be interpolated on the swath:
list_input_var_l2c = {'ucur': ['ucur', 'cur', 0], 'vcur': ['vcur', 'cur', 0]}

# -----------------------#
# L2D computation
# -----------------------#
# config name for L2d:
config_l2d = ''
# Length resolution to select neighbors (multiplication factor):
resol_spatial_l2d = 1
# Temporal resolution to select neighbors (multiplication factor):
resol_temporal_l2d = 1
# Grid resolution for l2d (lat, lon) grid (in degrees):
posting_l2d = (0.1, 0.1)
# Time domain: (start_time, end_time, dtime) in days:
time_domain = (7, 23, 1)
# Spatial domain (lon_min, lon_max, lat_min, lat_max):
spatial_domain = [0, 360, -90, 90]
# List of variables to be interpolated on the grid:
list_input_var_l2d = {'ucur': ['ucur', 'cur', 0], 'vcur': ['vcur', 'cur', 0]}

