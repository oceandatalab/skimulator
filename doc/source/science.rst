.. _science:

################################
SKIM Simulator
################################
Lucile Gaultier (OceanDataLab)

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

   FIG. 1: 5-day worth of SKIM simulated data in a global configuration with the science orbit.

.. _ProposedSKIMorbits:

Proposed SKIM orbits
---------------------
The software uses as an input the ground-tracks of the satellite orbit.

+-------------+--------------+--------------+------------+-------------+-----------+
|             | Repeat Cycle | Repeat Cycle | Sub-cycles | Inclination | Elevation |
|             | (days)       | (Orbits)     | (days)     |             | (km)      |
+=============+==============+==============+============+=============+===========+
| sentinel 1  |       12     |     175      |   6        |    90.18    |  698      |
+-------------+--------------+--------------+------------+-------------+-----------+
| **metop**   |      **29**  |   **412**    |   **5**    |  **98.63**  | **817**   |
+-------------+--------------+--------------+------------+-------------+-----------+
| fast sampl  |        3     |      43      |   0        |    98.5     |  775      |
+-------------+--------------+--------------+------------+-------------+-----------+
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


::

    cycle = 29
    elevation = 817000


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

   FIG. 2: scheme of the SKIM geometry with 4 beams at 12 degrees and 1 beam at 6 degree for figure a and 5 beams at 12 degrees and 2 beams at 6 degree for figure b.

Interpolation of model variables on the SKIM grid and nadir track
--------------------------------------------------------------------------
A list of model variables should be given to the skimulator, as well as a list
grids if the coordinates differ from one variable to another.
 The naming of the netcdf files should be
:math:`[pattern\_model]\_[pattern\_variable].nc`, where :math:`pattern\_variable`
is a string.
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
    mssxl = mssu * \cos(mssd)^2 + mssc * \sin(mssd)^2 \\
    mssyl = mssu * \sin(mssd)^2 + mssc * \cos(mssd)^2 \\
    mssxyl = (mssu - mssc) * \frac{\sin(2 * mssd)}{2}


Computation of short range MSS:

.. math::
   nwr = \sqrt{(uwnd - ucur)^2 + (vwnd - vcur)^2} \\ 
   wrd = \pi / 2 - arctan2(vwnd - vcur,\ uwnd - ucur) \\
   mssshort = \log(nwr + 0.7) * 0.009 \\
   mssshort[mssshort < 0] = 0
    
Directionality for short wave mss (if 0.5: isotrophic)

.. math::
   facssdw = 0.6 \\
   mssds = facssdw * mssshort \\
   msscs = mssshort - mssds \\
   mssxs = msscs * \sin(wrd)^2 + mssds * \cos(wrd)^2 \\
   mssys = mssds * \sin(wrd)^2 + msscs * \cos(wrd)^2 \\
   mssxys = |mssds - msscs| * \sin(2* wrd) \\

Computation of total MSS:

.. math::
   mssx = mssxs + mssxl \\
   mssy = mssys + mssyl \\
   mssxy = mssxys + mssxyl \\

:math:`\sigma^0` on water is computed from the total MSS:

.. math::
    B = -0.5 * \tan(beam)^2 * \frac{(\cos(azimuth)^2 * mssy + \sin(azimuth)^2 *mssx -\sin(2*azimuth)*mssxy)}{mssx * mssy} \\

.. math::
    A = \frac{R^2}{(2 * \cos(beam)^4 * \sqrt{mssx * mssy}} \\
    \sigma^0_{water} =  A \exp(B)

with :math:`R^2=0.55` which is a typical value for the tropics in Ka band.
Note that R depends on the radar frequency, water temperature and salinity
(eg :math:`R^2=0.50` for 3ºC water).


In the presence of ice, we use the concentration of sea ice :math:`C_{ice}`
and assume that :math:`\sigma^0_{ice}` is constant (:math:`\sigma^0_{ice}=2.5`
for 6º beam and :math:`\sigma^0_{ice}=1` for 12º beam).

.. math::
   \sigma^0 = (1 - C_{ice}) * \sigma^0_{water} +  C_{ice} * `\sigma^0_{ice}

Finally, the instrumental error is a random number proportional to
:math:`\sigma^0`


.. _Fig4:

.. figure:: ../images/Fig4.png
   :alt: Model current Intrumental noise 

   FIG. 4: Model interpolated currents and the corresponding intrumental error.


Wave bias
````````````````````````
The geophysical doppler includes also part of the currents due to the Stokes
drift. This componant is later refered as the current wave bias :math:`Uwb`.
To compute the wave bias, the Stoke drift and the wind are necessary:
The relation between the wave bias and the Stoke drift is parametrized using a
radial and perpendicular G parameter.

Compute the wind stress on the surface of the ocean:

.. math::
   nwr = \sqrt((u_{wind} - u_{cur})^2 + (v_{wind} - v_{cur})^2)

Compute radial (:math:`G_r`) and perpendicular (:math:`G_p`) parameter

.. math::
    G_r = 25 * (0.82 * \log(0.2 + \frac{7}{nwr})) * (1-tanh((beam - 25)/10))\\
    G_p = 0

Compute the wave bias

.. math::
    Uwb = G_r * ur_{uss} + G_p * up_{uss}

This wave bias can be corrected assuming that we can compute it using
neighbors in different azimuthal direction.This corrected componant is called
the remaining wave bias.
Near the coast, not all azimuth are
available and thus the drift remaining bias is higher than in the open ocean.

.. _Fig5:

.. figure:: ../images/Fig5.png
   :alt: Model current geophysical noise

   FIG. 5: Model interpolated currents and the corresponding wave remaining bias.

Total error
```````````
All previous errors are added to compute the total error. 

.. _Fig6:

.. figure:: ../images/Fig6.png
   :alt: Model current with noise

   FIG. 6: Model interpolated currents and the corresponding total noise (instrumental + geophysical).


Simulation of errors for the nadir altimeter
============================================
Two main components of the nadir altimetry error budget are simulated: the
altimeter noise and residual wet-tropospheric errors. For the altimeter noise,
the noise follow a spectrum of error consistent with global estimates from the
Jason-2 altimeter.
The wet-tropospheric residual errors (not implemented yet) are generated using
the simulated wet-tropospheric signal and radiometer beam convolution described
in SWOT Simulator documentation.

.. raw:: latex

    \newpage

The software
=============
The software is written in Python, and uses Numpy, Scipy amd netCDF4 python
libraries. Pyresample is also required for L2C computation and faster
interpolation.
All the parameters that can be modified by the user are read in a
params file (e.g. params.py) specified by the user. These parameters are
written in :ref:`yellow <params>` later on and are linked to their location in
the params file example.

The software is divided in 6 modules:

* :mod:`run_simulator.py` is the main program that runs the simulator.

* :mod:`mod_run.py` contains interpolations and data construction functions.

* :mod:`build_swath.py` generates the SKIM geometry and save several coordinates and angular variables in a netcdf file.

* :mod:`build_error.py` generates all the errors on the swath.

* :mod:`rw_data.py` contains all the classes to read and write model and SKIM data (in netcdf).

* :mod:`mod_tools.py` contains miscellaneous functions (algebraic functions and generation of random coefficients).

* :mod:`regridding.py` contains reconstruction for L2c products functions

* :mod:`mod_uwb_corr.py` contains function to correct the wave bias using neighbors


Inputs
-------
You can provide to the simulator a list of model outputs in netcdf. You need
to have at least the meridional and zonal currents to compute error-free radial
L2B velocities and SSH if you want to simulate nadir data. Wind and MSS
are necessary to compute instrumental noise (proportional to :math:`sigma^0`),
Stokes drift and Wind are needed to compute wave bias. Ice concentration should
also be provided to improve the computation of sigma0 in polar areas.
Any other variables provided to the skimulator will be interpolated on the
SKIM points.

A list of files (in .txt format) is provided using :ref:`file_input <params-file>`
in the parameter file.

The extension should not be provided in the list_of_files:

::

   model_0001_
   model_0002_
   model_0003_

The corresponding model file for a variable `var` should be  model_0001_var.nc
For example, if all the variables are in the same file
:math:`myfile\_[date].nc`, the list of file will be:

::

   myfile_date1
   myfile_date2
   myfile_date3

FIG 19: Example of a list of files, a list is provided in the example directory.

The grid files are provided as a list in the parameter file, using the key
:ref:`file_grid_model  <params-model>`. Make a list of all grid files that are necessary for
your variables, the correspondance between the variable and the grid is given
in a number in the :ref:`list_input_var  <params-model>`.
If no file_grid_model is provided, The skimulator is going to use the first
file of your list and data in this file will be ignored.


It is possible to generate the SKIM sampling alone, without using any model as an
input. If the name of the list of files
(:ref:`file_input <params-file>`) is set to `None`, then only SKIM grid files
will be generated.


The module :mod:`rw_data.py` is used to read model data. For any standard
netcdf format, the data can be read using
:ref:`model <params-model>` = MODEL_NETCDF, which is the
:ref:`model<params-model>` default value. The user needs to specify the list of
latitude (:ref:`lat <params-model>`) and longitude
(:ref:`lon <params-model>`) variables names corresponding to the list of grid
files provided in :ref:`file_grid_model  <params-model>`.
All other variables that are to be read, are added to the dictionnary
:ref:`list_input_var <params-model>`:

::

  list_input_var = {'key': [[variable\_name], [variable\_extension_file], [number\_corresponding\_to_grid_file]]}

The following table summarizes the key that are required to compute
instrumental noise and wave bias:

+-------+---------------------------+-----------------------------------+
| Key   | corresponding variable    | Necessary to compute ...          |
+=======+===========================+===================================+
| ucur  | Zonal total current       | Wave bias, radial current         |
+-------+---------------------------+-----------------------------------+
| vcur  | Meridional total current  | Wave bias, radial current         |
+-------+---------------------------+-----------------------------------+
| uuss  | Zonal Stokes drift        | Wave bias                         |
+-------+---------------------------+-----------------------------------+
| vuss  | Meridional Stokes drift   | Wave bias                         |
+-------+---------------------------+-----------------------------------+
| ice   | Ice concentration         | Instrumental noise if there is ice|
+-------+---------------------------+-----------------------------------+
| mssd  | Direction long wave mss   | Instrumental noise                |
+-------+---------------------------+-----------------------------------+
| mssx  | Zonal MSS                 | Instrumental noise                |
+-------+---------------------------+-----------------------------------+
| mssy  | Meridional MSS            | Instrumental noise                |
+-------+---------------------------+-----------------------------------+
| ssh   | Sea Surface Height        | Nadir SSH                         |
+-------+---------------------------+-----------------------------------+
| uwnd  | Zonal wind                | Wave bias, Instrumental noise     |
+-------+---------------------------+-----------------------------------+
| vwnd  | Meridional wind           | Wave bias, Instrumental noise     |
+-------+---------------------------+-----------------------------------+



Netcdf data that follow WW3 format can automatically be read
using :ref:`model <params-model>` = WW3 and there is no need to specify the
longitude or latitude variables name.
Below is an example of :ref:`list_input_var  <params-model>` for WW3 model
(all variables are on the same grid):

::

  file_grid_model = ('grid.nc', )
  lon = ('longitude', )
  lat = ('latitude', )
  list_input_var = {'ucur': ['ucur', 'cur', 0],
                    'vcur': ['vcur', 'cur', 0],
                    'uuss': ['uuss', 'uss', 0],
                    'vuss': ['vuss', 'uss', 0],
                    'ice': ['ice', 'ice', 0],
                    'mssd': ['mssd', 'msd', 0],
                    'mssx': ['mssx', 'mss', 0],
                    'mssy':['mssy', 'mss', 0],
                    'ssh': ['wlv', 'wlv', 0],
                    'uwnd': ['uwnd', 'wnd', 0],
                    'vwnd': ['vwnd', 'wnd', 0]}

Below is an example of :ref:`list_input_var  <params-model>` for a model with
an Arakawa grid (type C):

::

  file_grid_model = ('grid_u.nc', 'grid_v.nc', 'grid_T.nc')
  lon = ('lon_u', 'lon_v', 'lon_t')
  lat = ('lat_u', 'lat_v', 'lat_t')
  list_input_var = {'ucur': ['u', 'cur', 0],
                    'vcur': ['v', 'cur', 1],
                    'uuss': ['uuss', 'uss', 0],
                    'vuss': ['vuss', 'uss', 1],
                    'ice': ['ice', 'ice', 2],
                    'ssh': ['wlv', 'wlv', 2],
                    'uwnd': ['u10', 'wnd', 0],
                    'vwnd': ['v10', 'wnd', 1]}

The coordinates are supposed to
be in degrees and current variables in m/s in the program.

To refer timestamp properly in netcdf files, fill in the
:ref:`first_time  <params-model>` key
following  :ref:`first_time`='yyyy-mm-ddTHH:MM:SSZ'  <params-model>`.  By default,
:ref:`first_time`='2011-11-15T00:00:00Z'  <params-model>`.

If there is a ice_mask varying in time, set :ref:`ice_mask  <params-model>`
to False to
recompute the mask at every cycle.

The number of time in each file should be constant for all the files
and specified in the
:ref:`dim_time  <params-model>` parameter. The time step between two inputs
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
radian (:ref:`list_pos <params-skimswath>`), an angle on the sensor in degree
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
:ref:`grid <params-file>` can bet set to `regular` and RectBivariateSpline
interpolation from scipy is used. In all cases, :ref:`grid <params-file>` option
can be set to `irregular` and pyresample is used for the interpolation if the
module is installed. If the grid is irregular and pyresample is not installed,
mgriddata from scipy interpolates data is used with either the 'linear'
(:ref:`interpolation <params-output>` ='linear') or 'nearest' neighbor
(:ref:`interpolation <params-output>` ='nearest') option. In case of large
domain, this last solution for the mgriddata, interpolation can be slow and
even trigger memory error. The use of the ‘nearest’ interpolation is necessary
to avoid memory error though the derivative of the current can be significantly
altered using this interpolation method.
The :ref:`list_output  <params-output>` key  concerns the list of variables that you want to be
store in the netcdf files. The most common key are summarized in the table
below, you can add any other key that you want to interpolate on the SKIM grid:

+------------+--------------------------------+-------------------------------+
|Key         | Long name                      | Required for ...              |
+============+================================+===============================+
| ur_true    | Radial error free velocity     |                               |
+------------+--------------------------------+-------------------------------+
| ur_obs     | Radial velocity with errors    |                               |
+------------+--------------------------------+-------------------------------+
| ucur       | Meridional true current        | Radial velocity, Wave bias    |
+------------+--------------------------------+-------------------------------+
| vcur       | Zonal true current             | Radial velocity, Wave bias    |
+------------+--------------------------------+-------------------------------+
| uuss       | Meridional Stokes drift        | ur_obs, Wave bias             |
+------------+--------------------------------+-------------------------------+
| vuss       | Zonal Stokes drift             | ur_obs, Wave bias             |
+------------+--------------------------------+-------------------------------+
| instr      | Instrumental error             | instrumental error, ur_obs    |
+------------+--------------------------------+-------------------------------+
|radial_angle| Azimuthal angle                | all                           |
+------------+--------------------------------+-------------------------------+
| uwnd       | Meridional wind                | wave bias, ur_obs, instr      |
+------------+--------------------------------+-------------------------------+
| vwnd       | Zonal wind                     | wave bias, ur_obs, instr      |
+------------+--------------------------------+-------------------------------+
| mssx       | Meridional MSS                 | instr, ur_obs                 |
+------------+--------------------------------+-------------------------------+
| mssy       | Zonal MSS                      | instr, ur_obs                 |
+------------+--------------------------------+-------------------------------+
| mssxy      | Mixed MSS                      | instr, ur_obs                 |
+------------+--------------------------------+-------------------------------+
| uwb        | Wave bias                      | wave bias, ur_obs, uwb_corr   |
+------------+--------------------------------+-------------------------------+
| uwb_corr   | Remaining wave bias            | wave bias, ur_obs             |
+------------+--------------------------------+-------------------------------+
| sigma0     | NRCS                           | instr, ur_obs                 |
+------------+--------------------------------+-------------------------------+
| ssh_obs    | Sea Surface Height with errors | ssh_obs, nadir                |
+------------+--------------------------------+-------------------------------+
| ssh_true   | Error free Sea Surface Height  | ssh_obs,                      |
+------------+--------------------------------+-------------------------------+


To compute each error, set the corresponding parameter to True
(:ref:`instr <params-error>`, :ref:`ice <params-error>`,
:ref:`uwb <params-error>`, :ref:`nadir <params-error>`).
If :ref:`nadir <params-error>` is True, then :ref:`ncomp1d  <params-error>`
should be specified to compute random error realisations. 
If :ref:`instr <params-error>` is True, the :ref:`snr_coeff <params-error>`
coefficient (to retrieve instrumental noise from sigma) should be specified.


All the computed errors are saved in the output netcdf file. The observed SSH
(SSH_obs) is also computed by adding all the computed errors to the SSH from
the model (SSH_model) when model data are provided.

Two errors are considered in the nadir. The first one is the instrument error,
which follows the 1d spectrum computed from the current altimeters. You have
to set (:ref:`nadir <params-error>`) to True to compute this error.


L2C 2d currents
---------------
L2b radial current can be projected on a along swath and across swath grid
using the skiml2c command with the same parameter file as the one used for
the L2b current production. It uses the L2b produced previously as an input.
The along track and across track grid resolution is specified in
:ref:`posting <params-l2c>`. By default, the spatial resolution of the grid is
5 km.
The L2c reconstrunction uses neighbors to project and interpolate the current. 
The length resolution to select neighbors can be set in
:ref:`resol <params-l2c>`. Coefficient for the OI are deacreasing exponentially
with the distance to the pixel.
As data around nadir are particularly noisy (all radial velocity are
parrallels), one can mask them by specifiying the distance in km from nadir
where data are to be masked :ref:`ac_threshold <params-l2c>`.

The L2c outputs contains along track, across track, meridional and zonal
current reconstructed from the error-free and radial velocity with errors.
True along track, across track meridional and zonal velocity interpolated
from the model inputs are also stored for diagnosis purposes.


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

To compute L2b products:

.. code-block:: python

   >>>skimulator [your_parameter_file]


You can compute L2c products after L2b files have been produced, keep the same
parameter file and run:

.. code-block:: python

   >>>skiml2c [your_parameter_file]



.. _params:

Example of Params.txt for SKIM-like data
``````````````````````````````````````````

.. _params-file:

.. literalinclude:: params.py
   :lines: 1-32

.. _params-skimswath:

.. literalinclude:: params.py
   :lines: 34-62

.. _params-model:

.. literalinclude:: params.py
   :lines: 64-107

.. _params-output:

.. literalinclude:: params.py
   :lines: 109-124

.. _params-error:

.. literalinclude:: params.py
   :lines: 125-160

.. _params-l2c:

.. literalinclude:: params.py
   :lines: 159-


References:
===========
.. _SKIM_proposal_2017:
The SKIM team Sea surface KInematics Multiscale monitoring, full proposal for
ESA EE9, 2017, The SKIM team

