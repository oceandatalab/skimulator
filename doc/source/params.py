# -----------------------#
# Files and directories
# -----------------------#
## -- Get the user home directory
from os.path import expanduser
import os
import math
home = expanduser("~")

# ------ Name of the configuration (to build output files names) 
config = [yourconfig]
# ------ Directory that contains orbit file:
dir_setup = os.path.join([yourpath], 'skimulator', 'data')
# ------ Directory that contains your own inputs:
indatadir = '[yourpath_to_yourdata]/'
# ------ Directory that contains your outputs:
outdatadir = '[yourpath_to_outputs]/'
# ------ Orbit file:
satname = [chosenorbit]
filesat = os.path.join(dir_setup, [chosenorbit])
# ------ Number of days in orbit (optional if specified in orbit file)
satcycle = 29
# ------ Satellite elevation (optional if specified in orbit file)
sat_elev = 817 * 10**3
# ------ Order of columns (lon, lat, time) in orbit file
# (default is (0, 1, 2) with order_orbit_col = None)
order_orbit_col = None
# , dir_setup+os.sep+'orbjason.txt', dir_setup+os.sep+'orbaltika.txt' ]
# ------ Number of processors for parallelisation purposes
proc_number = [number of processor (integer)]
# ------ Deactivate printing of progress bar to avoid huge log
progress_bar = True or False

# -----------------------#
# SKIM swath parameters
# -----------------------#
# ------ Satellite grid file root name:
# 	 (Final file name is root_name_[numberofpass].nc)
filesgrid = os.path.join(outdatadir, '{}_grid'.format(config))
or filesgrid = os.path.join(outdatadir, '[your_grid_root_name]')
# ------ Force the computation of the satellite grid:
makesgrid = True or False
# ------ Give a subdomain if only part of the model is needed:
#	 (modelbox=[lon_min, lon_max, lat_min, lat_max])
# 	 (If modelbox is None, the whole domain of the model is considered)
modelbox =  None or [yourlon_min, yourlon_max, yourlat_min, yourlat_max]
#------- Rotation speed of the antenna (in tr/min)
rotation_speed = rotation depends on the chosen config
#------- List of position of beams:
list_pos = (0, [angle_in_rad], [angle_in_rad] ...)
#------- List of angle of beams in degrees:
list_angle = ([incidence], [incidence], [incidence] ...)
#------- List of timeshift as regard to nadir for 12 degree beams:
list_shift = (1, 3, 2 ...)
#------- Cycle duration
cycle = 0.0096
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
file_input = os.path.join(indatadir, [your_list_of_file_name.txt]' or None
# ------ Type of model data:
#	 (Optional, default is NETCDF_MODEL and reads netcdf3 and netcdf4 files)
#	 (Other option is WW3)
model = 'WW3' or 'NETCDF_MODEL'
# ------ First time of the model
first_time = 'yyyy-mm-ddTHH:MM:SSZ'
# ------ Grid file name
file_grid_model = (os.path.join(indatadir, [yourgridfileu]),
                   os.path.join(indatadir, [yourgridfilev]),)
# ------ Specify if there is a ice mask for high latitudes
#        (if true, mask is recomputed at each cycle)
ice_mask = False or True
# ------ Type of grid: 
#        'regular' or 'irregular', if 'regular' only 1d coordinates 
#        are extracted from model       
grid = 'regular'
# ------ Specify list of variable:
list_input_var = {'ucur': [[u_var], [vel_ext], 0], 'vcur': [[v_var], [v_ext], 1],
                  'uuss': [[uuss_var], [uss_ext], 0], 'vuss': [[vuss_var], [uss_ext], 1],
                  'ice': [[ice_var], [ice_ext], 0], 'mssd': [[mssd_var], [msd_ext], 0],
                  'mssx': [[mssx_var], [mss_ext], 0], 'mssy':[[mssy_var], [mss_ext]],
                  'ssh': [[ssh_var], [ssh_ext], 0],
                  'uwnd': [[uwnd_var], [wnd_ext], 0], 'vwnd': [[vwnd_var], [wnd_ext], 1]}
# ------ Specify longitude variable:
lon = ('longitude', 'longitude')
# ------ Specify latitude variable:
lat = ('latitude', 'latitude')
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
interpolation = 'nearest' or 'linear'
list_output = ['ssh_obs', 'ur_true', 'ucur', 'vcur', 'uuss', 'vuss', 'instr',
               'radial_angle', 'vwnd', 'mssx', 'mssy', 'mssxy', 'uwb',
               'ssh_true', 'ssh', 'ice', 'mssd',
               'vindice', 'ur_obs', 'uwnd', 'sigma0']

# -----------------------# 
# SKIM error parameters 
# -----------------------# 
# ------ File containing random coefficients to compute and save 
#	 random error coefficients so that runs are reproducible:
#        If file_coeff is specified and does not exist, file is created
#	 If you don't want runs to be reproducible, file_coeff is set to None
file_coeff = None  or os.path.join(outdatadir, 'Random_coeff.nc')
# Compute instrumental nadir noise:
nadir = True
# ------ Number of random realisations for instrumental and geophysical error 
#        (recommended ncomp=2000), ncomp1d is used for 1D spectrum, and ncomp2d
#        is used for 2D spectrum (wet troposphere computation):
ncomp1d = 3000
ncomp2d = 2000
# ------- Instrument white noise error
instr = True or False
# ------- Coefficient SNR to retrieve instrumental noise from sigma, 
#         Recommanded value for 1024 pulses: 3e-2, for 512 pulses: 3sqrt(2)e-3
snr_coeff = 6e-3
if '2018_8b' in config:
    snr_coeff = 1.4142*6e-3

# ------- Wave bias
uwb = True or False


## -- Geophysical error
## ----------------------
# ------ Consider ice in sigma0 computation
ice = True or False
#### Not implemented yet
# ------ Rain error (True to compute it):
wet_tropo = False or True

# -----------------------#
# L2C computation
# -----------------------#
# Length resolution to select neighbors (in km):
resol = 40
# Grid resolution for l2c (alongtrack, acrosstrack) grid (in km):
posting = 5
# Remove noisy data around nadir (in km):
ac_threshold = 20
# List of variables to be interpolated on the swath:
list_input_var_l2c = {'ucur': ['ucur', 'cur', 0], 'vcur': ['vcur', 'cur', 1]}

# -----------------------#
# L2D computation
# -----------------------#
# Length resolution to select neighbors (in km):
resol_spatial_l2d = 50
# Temporal resolution to select neighbors (in days):
resol_temporal_l2d = 8
# Grid resolution for l2d (lat, lon) grid (in degrees):
posting_l2d = (0.1, 0.1)
# Time domain: (start_time, end_time, dtime) in days:
time_domain = (5, 25, 1)
# Spatial domain (lon_min, lon_max, lat_min, lat_max):
spatial_domain = [0, 360, -90, 90]
# List of variables to be interpolated on the grid:
list_input_var_l2d = {'ucur': ['ucur', 'cur', 0], 'vcur': ['vcur', 'cur', 1]}

