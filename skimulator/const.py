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
# - Light speed (m/s)
C = 2.998*10**8


# ###################################
# # SKIM INSTRUMENT CARACTERISTICS ##
# ###################################
# - Satellite elevation (m)
sat_elev = 817*10**3
# - Baseline (m)
B = 10
# (in Hz)
Fka = 35.75 * 10**9
# Satellite cycle (S1) in days
satcycle = 29.
# Plateform error model for reconstruction in 3D
theta1 = 0.7E-4
theta0 = 0.1
gamma0 = 1E-4


# ###################################
# # OTHER PARAMETERS               ##
# ###################################
# - Radius to interpolate locally model data on the swath (in km)
#  data are selected every xal_step points and on a radius of radius_interp
radius_interp = 100.
# - Sampling to interpolate locally model data on the swath (in km)
#  Data are selected every xal_step points and on a radius of radius_interp
xal_step = 20.
