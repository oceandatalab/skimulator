.. _science:
################################
SKIM Simulator
################################
Lucile Gaultier

OceanDataLab

.. role:: red
.. toctree::
   :maxdepth: 2
   :numbered:

Abstract:
=========
This software simulates sea surface radial current (Level-2) synthetic
observations of the proposed SKIM mission that can be applied to an ocean
general circulation model (OGCM), allowing the exploration of ideas and
methods to optimize information retrieval from the SKIM Mission in the future.
From OGCM currents and Stoke drift inputs, the software generates SKIM-like
outputs on a rotating geometry around the orbit ground track, as well as
outputs from a nadir altimeter. Some measurement error and noise are generated
according to technical characteristics published by the SKIM project team. Not
designed to directly simulate the payload instrument performance, this SKIM
simulator aims at providing statistically realistic outputs for the science
community with a simple software package released as an open source in Python.
The software is scalable and designed to support future evolution of orbital
parameters, error budget estimates from the project team and suggestions from
the science community.

Simulation of the SKIM sampling over synthetic Sea Surface current
==================================================================
From a global or regional OGCM configuration, the software generates radial
velocities, including instrumental and geophysical noise, for each beam.
Note that for an accurate instrumental ang geophysical noise, various forcings
 are needed such as mean square slope, Stockes drift, wind, ice ...

.. _Fig1:

.. figure:: ../images/Fig1.png 
   :alt: Science SKIM orbit

   FIG. 1: 5-day worth of SKIM simulated data in a global configuration with
             the science orbit.

.. _ProposedSKIMorbits:

Proposed SKIM orbits
---------------------
The software uses as an input the ground-tracks of the satellite orbit.

+-------------+--------------+--------------+------------+-------------+-----------+
|             | Repeat Cycle | Repeat Cycle | Sub-cycles | Inclination | Elevation |
|             | (days)       | (Orbits)     | (days)     |             | (km)      |
+=============+==============+==============+============+=============+===========+
| sentinel 1  |       12     |     175      |   6        |    90.18    |  698      |
| **metop**   |      **29**  |   **412**    |   **5**    |  **98.63**  | **817**   |
| fast sampl  |        3     |      43      |   0        |    98.5     |  775      |
| fast sampl  |        8     |     113      |   0        |    98.8     |  845      |
| scanning    |              |              |            |             |           |
+-------------+--------------+--------------+------------+-------------+-----------+

The ground-track coordinates corresponding to these orbits are given as input
ASCII files of 3 columns (longitude, latitude, time) for one complete cycle
sampled at every  ~5 km. The first ascending node has been arbitrarily set to
270 degree of longitude, but the user can shift the orbit by any value in
longitude. The default orbit is metop.

Other orbit files of the same format (time, longitude, latitude) can also be
used as an input. To avoid distortions in the grid, we recommend a minimum of
10km sampling between the ground-track points of the orbit.

Note that the first two commented lines of the files concerns the satellite
 cycle (in days) and elevation (in km).
.. code-block:: python
# cycle = 29
# elevation = 817000

If these lines does not exist, the skimulator will look for these values in the
 parameter file or take default value (cycle = 29 days and elevation = 817000)

The SKIM geometry
-----------------

From the orbit nadir ground track the software generates a grid covering the
swath over one satellite cycle. The longitude and latitude coordinates as well as
the time are referenced for each grid point. A scheme of the SKIM geometry is
presented on :ref:`Fig. 2 <Fig2>`.
The SKIM grid is stored by pass (e.g. 412 ascending passes and 412 descending
passes for the Metop orbit). A pass is defined by an orbit starting at
the lowest latitude for ascending track and at the highest latitude for
descending track. The first pass starts at the first lowest latitude crossing
in the input file, meaning that ascending passes are odd numbers and descending
passes are even numbers.

.. _Fig2:

.. figure:: ../images/Fig2.png
   :alt: SKIM geometry

   FIG. 2: scheme of the SKIM geometry with 4 beams at 12 degrees and 1 beam at
           6 degree for figure a and 5 beams at 12 degrees and 2 beams at 6
           degree for figure b.

Interpolation of model variables on the SKIM grid and nadir track
--------------------------------------------------------------------------
A list of model variables can be given to the skimulator. All the variables
 should have the same coordinates. The naming of the netcdf files should be
 :math:`[pattern_model]_[pattern_variable].nc`, where :math:`pattern_variable`
 is a 3-character string.
All input variables must be given at the same regular time step.

The absolute time of the first time step is zero
and corresponds to the beginning of pass 1. A first date can be provided in
order to have a consistent timestamps in the netcdf file.  All provided
 variables are
interpolated on the SKIM grid and nadir track for each pass and successive
cycles if the input data exceeds 1 cycle. Current and Stokes drift are then
 projected along the radial component as only measurement along the radial axis
 can be made.

No interpolation is made in time (the model file with the closest time step is
chosen). This avoids contaminations of the rapid signals (e.g. internal waves)
if they are under-sampled in the model outputs. However, note that locally,
sharp transitions of the variable along the swath may occur if the satellite
happens to be over the domain at the time of transition between two time steps.
By default a linear 2D spatial interpolation is performed to compute the variable
data on the SKIM grid.


:ref:`Fig. 3a <Fig3>` shows an input current as an example.
:ref:`Fig 3b <Fig3>` is the interpolated and radial component of the current.

.. _Fig3:

.. figure:: ../images/Fig3.png
   :alt: Model current and model current interpolated and projeced on SKIM grid

   FIG. 3: Model interpolated currents and the corresponding radial currents.

Simulation of errors
====================

Instrumental errors
````````````
The instrumental error corresponds to the geometric doppler.
This componant is proportional to sigma0 with a SNR specified in the parameter
 file. 
The following variables are needed to compute long range and short range mss:
 mssu, mssc, mssd, uwnd, vwnd, ucur, vcur.

Computation of long range MSS:
.. math::
    mssxl = mssu * \cos(mssd)^2 + mssc * \sin(mssd)^2
    mssyl = mssu * \sin(mssd)^2 + mssc * \cos(mssd)^2
    mssxyl = (mssu - mssc) * sin(2 * mssd) / 2


Computation of short range MSS:
.. math::
   nwr = \sqrt{(uwnd - ucur)^2 + (vwnd - vcur)^2}
   wrd = \pi / 2 - arctan2(vwnd - vcur, uwnd - ucur)
   mssshort = \log(nwr + 0.7) * 0.009
   mssshort[mssshort < 0] = 0

   #Directionality for short wave mss (if 0.5: isotrophic)
   facssdw = 0.6
   mssds = facssdw * mssshort
   msscs = mssshort - mssds
   mssxs = msscs * \sin(wrd)^2 + mssds * \cos(wrd)^2
   mssys = mssds * \sin(wrd)^2 + msscs * \cos(wrd)^2
   mssxys = \abs(mssds - msscs) * \sin(2* wrd)


Computation of total MSS:
.. math::
   mssx = mssxs + mssxl
   mssy = mssys + mssyl
   mssxy = mssxys + mssxyl

Computation of sigma0:

To compute this noise, the short and long range MSS (mean square slopes) needs
to be computed
 For now a random
noise has been implemented. The amplitude of the noise depends on the beam
angle and the azimuth.
.. _Fig4:

.. figure:: ../images/Fig4.png
   :alt: Model current Intrumental noise 

   FIG. 4: Model interpolated currents and the corresponding intrumental error.


Stoke drift remaining bias
````````````````````````
The geophysical noise contains only the bias on the doppler due to the Stoke
drift. The stoke drift is provided to the simulator, we retrieve the
corresponding radial Stoke componant measured by SKIM and correct it using
nearest negihbours in all azimuth. The corrected variable is called the Stoke 
drift remaing bias. Near the coast, not all azimuth are
available and thus the drift remaining bias is higher than in the open ocean.
.. _Fig5:

.. figure:: ../images/Fig5.png
   :alt: Model current geophysical noise

   FIG. 5: Model interpolated currents and the corresponding geophysical noise
           (for now, it contains only stoke drift remaining bias).

Total error
```````````
All previous errors are added to compute the total error. 

.. _Fig6:

.. figure:: ../images/Fig6.png
   :alt: Model current with noise

   FIG. 5: Model interpolated currents and the corresponding total noise
           (instrumental + geophysical).

Simulation of errors for the nadir altimeter
============================================
Two main components of the nadir altimetry error budget are simulated: the
 altimeter noise and residual wet-tropospheric errors. For the altimeter noise,
the noise follow a spectrum of error consistent with global estimates from the
Jason-2 altimeter. The wet-tropospheric residual errors are generated using the
simulated wet-tropospheric signal and radiometer beam convolution described in
SWOT Simulator documentation.

.. raw:: latex

    \newpage

The software
=============
The software is written in Python, and uses Numpy, Scipy amd netCDF4 python
libraries. All the parameters that can be modified by the user are read in a
params file (e.g. params.py) specified by the user. These parameters are
written in :ref:`yellow <params>` later on and are linked to their location in
the params file example.

The software is divided in 6 modules:

* :mod:`run_simulator.py` is the main program that runs the simulator.

* :mod:`build_swath.py` generates the SKIM geometry and save several
coordinates and angular variables in a netcdf file.

* :mod:`build_error.py` generates all the errors on the swath.

* :mod:`rw_data.py` contains all the classes to read and write model and SKIM
data (in netcdf).

* :mod:`mod_tools.py` contains miscellaneous functions (algebraic functions and
generation of random coefficients).


Inputs
-------
The inputs are current and Stoke drift model outputs in netcdf. Two lists of
file (in .txt format) are read by the software (one for the current and one for
the Stoke drift). It contains the grid file and all model outputs. The first
file in this list is the grid, therefore, if the grid is contained in the data
files, the first data file should be repeated. The directory that contains
input (:ref:`indatadir <params-file>`) and the names of the list of files
(:ref:`file_input <params-file>` for current and
:ref:`input_uss <params-file>`) are specified in the params file.

.. code-block:: python

   Grid_model.nc  
   model_0001.nc
   model_0002.nc
   model_0003.nc 

FIG 19: Example of a list of files, a real example is located in the example directory.

It is possible to generate the noise alone, without using any model as an
input. To generate the noise alone, the name of the list of files
(:ref:`file_input <params-file>`) should be set to `None`. Note that if you set
the `input_uss <params-error>` list to `None`, no Stoke drift bias will be
computed.


The module :mod:`rw_data.py` is used to read model data. For any standard
netcdf format, the data can be read using
:ref:`model <params-model>` = MODEL_NETCDF, which is the
:ref:`model<params-model>` default value. The user needs to specify the
latitude (:ref:`latu <params-model>` and :ref:`latv <params-model>`), longitude
(:ref:`lonu <params-model>` and :ref:`lonv <params-model>`), and velocity
(:ref:`varu <params-model>` and :ref:`varv <params-model>`) variable names for
both component. Netcdf data that follow WW3 format can automatically be read
using :ref:`model <params-model>` = WW3 and there is no need to specify the
longitude, latitude or current variables name. The coordinates are supposed to
be in degrees and current variables in m/s in the program, if the current is
not in m/s, specify the conversion factor in :ref:`vel_factor <params-model>`
(so that u*vel_factor is in m/s). If there is more than one time step in a
single file, a list of the time dimension for each file can be provided in
:ref:`dim_time <params-model>`. The time step between two inputs
(:ref:`timestep <params-model>`) and the number of steps that have to be
processed (nstep) can be modified in the params file. The value corresponding
to not a number can be specified in :ref:`model_nan <params-model>`.

Generation of the SKIM gemoetry
-------------------------------
The SKIM grid is generated in the :mod:`build_swath.py` module. The orbit file
(:ref:`filesat <params-file>`) is located in :ref:`dir_setup <params-file>` and
contains longitude, latitude and the corresponding time for each point of the
orbit (see section :ref:`ProposedSKIMorbits` for more details on the provided
orbit). The orbit is interpolated at the cycle duration time resolution. The
rotation speed of the antenna is specified in
:ref:`rotation_speed <params-skimswath>` in tr/min. The provided value has been
computed to keep an integer number of illumination in the macro-cycle. The
geometry of beams is provided with lists and give for each beam a position in
radian (:ref:`list_pos` <params-skimswath>), an angle on the sensor in degree
(:ref:`list_angle <params-skimswath>`), an order for the illumination
(:ref:`list_shift <params-skimswath>`) with the macro-cycle starting with the
nadir beam.
The generation of the SKIM grid can be made on the whole region of the model
(:ref:`modelbox <params-skimswath>` = `None`) or on a subdomain
(:ref:`modelbox <params-skimswath>` = [lon_min, lon_max, lat_min, lat_max]). To
compute the noise alone (:ref:`file_input = None <params-file>`), a
:ref:`modelbox <params-skimswath>` has to be provided. If there is no pass on
the required domain, the user can specify a shift in longitude
(:ref:`shift_lon <params-skimswath>`).

A netcdf file containing SKIM grid information is stored for each pass in the
output directory (:ref:`outdatadir <params-file>`) under the name
:ref:`filesgrid <params-skimswath>` _p[pass].nc. It contains nadir variables
(_nadir) and other beams variables (one beam per column in the same order as in
the parameter file) for the following variables: along track and across track
distance from the nadir (in km), longitudes (lon) and latitudes (lat), time,
inclination at nadir, the number of days in a cycle (cycle) and the distance
crossed by the satellite in a cycle (al_cycle). Once the SKIM grid has been
created, it is stored in :ref:`outdatadir <params-file>`. As long as the domain
(:ref:`modelbox <params-skimswath>` parameter) and the geometry do not change,
the grids do not have to be recomputed and :ref:`makesgrid <params-skimswath>`
can be set to False.
For a more convenient use, in the example the name of the output files are a
concatenation of a :ref:`config <params-skimswath>` name and
:ref:`satname <params-file>` orbit name.

Radial current and error fields
--------------------------------
At each pass, for each cycle, an output netcdf file containing the currents
interpolated from the model as well as the interpolated current projected on
the radial direction (if :ref:`file_input <params-file>` is not set to `None`)
and the different errors are created. The output file names are
:ref:`file_output <params-output>` _c[cycle]_p[pass].nc for the swath and
:ref:`file_output <params-output>` _c[cycle]_p[pass].nc for the nadir. The SSH
is interpolated on the SKIM grid. If the model grid is regular, option
:ref:`grid <param-file>` can bet set to `regular` and RectBivariateSpline
interpolation from scipy is used. In all cases, :ref:`grid <param-file>` option
can be set to `irregular` and pyresample is used for the interpolation if the
module is installed. If the grid is irregular and pyresample is not installed,
mgriddata from scipy interpolates data is used with either the 'linear'
(:ref:`interpolation <params-output>` ='linear') or 'nearest' neighbor
(:ref:`interpolation <params-output>` ='nearest') option. In case of large
domain, this last solution for the mgriddata, interpolation can be slow and
even trigger memory error. The use of the ‘nearest’ interpolation is necessary
to avoid memory error though the derivative of the current can be significantly
altered using this interpolation method.

To compute each error, set the corresponding parameter to True
(:ref:`instr <params-error>`, :ref:`uss <params-error>`.
The rms for the instrumental noise depends on the beam angle and on the
position related to the along track direction. .dat files (one file per beam
angle) contain the rms as a function of the position. The file names
(including the path) are provided as a list of file path in
:ref:`rms_instr <params-error>`.
There are two methods to compute the Stoke drift remaining bias. For both
methods, a stoke parameter :ref:`G <params-error>` is defined. It depends
highly on the wind intensity.
To use the parametrisation as a function of the distance and the uss of
neighboring beams, use the options :ref:`formula = True <params-error>`, and
define the :ref:`uss_std <params-error>` and `factor_errdcos <params-error>`
parameters corresponding to the region of interest. The bias is computed using
the following formula:
`errwb = G * uss_std * std(uss) *errdcos[beam] / factor_errdcos`
If :ref:`formula <params-error>` is set to False,
`errwb = G * sum(uss[beam] / errdcos[beam])`


All the computed errors are saved in the output netcdf file. The observed SSH
(SSH_obs) is also computed by adding all the computed errors to the SSH from
the model (SSH_model) when model data are provided. Note that if
:ref:`nbeam  <params-error>` ='both' (residual error due to path delay using
one beam and two beams are both computed), only the residual error due to path
delay using one beam is considered in the observed SSH.

Two errors are considered in the nadir. The first one is the instrument error,
which follows the 1d spectrum computed from the current altimeters. You have
to set (:ref:`nadir <params-error>`) to True to compute this error. The second
error is the path delay due to the wet troposphere and this error is computed
with the residual path delay error in the swath. The observed SSH (SSH_obs) is
computing by adding these two errors to the SSH interpolated from the model (SSH_model).

Getting started 
----------------
All information regarding the installation and running of the software are in
the :ref:`README <README>` file. An example of a :ref:`params.txt <params>` file is given
below.

Once you have installed skimulator, you can print help by typing in a python
or ipython window:

.. code-block:: python

   >>>import skimulator.M
   >>>help(skimulator.M)

with M the name of the module.

To run the example, type in any terminal:

.. code-block:: python

   >> skimulator ./example/params_example.py

for the SKIM simulator.


.. _params:

Example of Params.txt for SKIM-like data
``````````````````````````````````````````

.. _params-file:

.. literalinclude:: params.py
   :lines: 1-19

.. _params-swotswath:

.. literalinclude:: params.py
   :lines: 20-51

.. _params-model:

.. literalinclude:: params.py
   :lines: 53-90

.. _params-output:

.. literalinclude:: params.py
   :lines: 92-104

.. _params-error:

.. literalinclude:: params.py
   :lines: 106-


References:
===========
.. _SKIM_proposal_2017:
The SKIM team Sea surface KInematics Multiscale monitoring, full proposal for
ESA EE9, 2017, The SKIM team

