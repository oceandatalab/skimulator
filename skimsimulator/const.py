'''Constants that are defined for SKIM simulator. \n
It contains Earth parameters as well as SKIM instrument
and satellite caracteristics. '''

# ################################
# # EARTH CONSTANTS             ##
# ################################
# - Earth radius (m)
Rearth = 6378. * 10**3
# - Convert degree to km
deg2km = 111.11
# - Seconds in a day
secinday = 86400.


# ###################################
# # SKIM INSTRUMENT CARACTERISTICS ##
# ###################################
# - Satellite elevation (m)
sat_elev = 699*10**3
# - Baseline (m)
B = 10
# (in Hz)
Fka = 35.75 * 10**9
# Satellite cycle (S1) in days
satcycle = 12.

# - Light speed (m/s)
C = 2.998*10**8

# ###################################
# # OTHER PARAMETERS               ##
# ###################################
# - Radius to interpolate locally model data on the swath (in km)
#  data are selected every xal_step points and on a radius of radius_interp
radius_interp = 100.
# - Sampling to interpolate locally model data on the swath (in km)
#  Data are selected every xal_step points and on a radius of radius_interp
xal_step = 20.
