'''Main program:
Usage: run_simulator(file_param)  \n
If no param file is specified, the default one is exemple/params_exemple.txt \n
In the first part of the program, model coordinates are read and the
SKIM swath is computing accordingly. \n
The SKIM grid parameters are saved in netcdf files, if you don't want to
recompute them, set maksgrid (in params file) to False.\n

In the second part of the program, errors are computed on SKIM grid for
each pass, for each cycle. The error free SSH is the velocity interpolated
from the model at each timestep. Note that there is no temporal interpolation
between model files and thus if several files are used in the velocity
interpolation, some discontinuities may be visible. \n

OUTPUTS are netcdf files containing the requested errors, the error free
radial velocity and the radial velocity with errors. There is one file every
pass and every cycle.

\n
\n
#-----------------------------------------------------------------------
#                       Additional Documentation
# Authors: Lucile Gaultier and Clement Ubelmann
#
# Modification History:
# - Mar 2017:  Original by Clement Ubelmann and Lucile Gaultier
#
# Notes:
# - Written for Python 2.7,  Python 3.5, tested with Python 2.7, Python 3.5
#
# Copyright (c)
#
#-----------------------------------------------------------------------
'''
import os
import shutil
from scipy import interpolate
import numpy
import glob
import sys
import logging
try:
    import params as p
except:
    if os.path.isfile('params.py'):
        print("There is a wrong entry in your params file")
        import params
    else:
        print("Error: No params.py module found")
        sys.exit()
import skimsimulator
import skimsimulator.build_swath as build_swath
import skimsimulator.rw_data as rw_data
import skimsimulator.build_error as build_error
import skimsimulator.mod_tools as mod_tools
import skimsimulator.const as const

# Define logger level for debug purposes
logger = logging.getLogger(__name__)

# - Define global variables for progress bars
istep = 0
ntot = 1
ifile = 0


def run_simulator():

    # - Initialize some parameters values
    try:
        p.shift_lon = p.shift_lon
    except:
        p.shift_lon = None
    try:
        p.shift_time = p.shift_time
    except:
        p.shift_time = None
    try:
        model = p.model
    except:
        model = 'NETCDF_MODEL'
        p.model = model
    try:
        p.model_nan = p.model_nan
    except:
        p.model_nan = 0.
    try:
        p.vel_factor = p.vel_factor
    except:
        p.vel_factor = 1.
    try:
        p.nadir = p.nadir
    except:
        p.nadir = True
    try:
        p.grid = p.grid
    except:
        p.grid = 'regular'
    # - Progress bar variables are global
    global istep
    global ntot

    # - Read list of user model files """
    if p.file_input is not None:
        list_file = [line.strip() for line in open(p.file_input)]
    if p.uss is True:
        list_file_uss = [line.strip() for line in open(p.input_uss)]
    else:
        list_file = None
    # - Read model input coordinates '''
    # if no modelbox is specified (modelbox=None), the domain of the input
    # data is taken as a modelbox
    # coordinates from the region defined by modelbox are selected
    if p.file_input:
        model_data = eval('rw_data.' + model
                          + '(file=p.indatadir+os.sep+list_file[0])')
    if p.modelbox:
        modelbox = numpy.array(p.modelbox, dtype='float')
        # Use convert to 360 data
        modelbox[0] = (modelbox[0]+360) % 360
        if modelbox[1] != 360:
            modelbox[1] = (modelbox[1]+360) % 360
    else:
        if p.file_input:
            modelbox = model_data.calc_box()
        else:
            print('modelbox should be provided if no model file is provided')
            sys.exit()
    if p.file_input:
        model_data.read_coordinates()
        # Select model data in the region modelbox
        if p.grid == 'regular':
            if modelbox[0] < modelbox[1]:
                model_data.model_index_lonu = numpy.where(((modelbox[0]-1) <= model_data.vlonu) & (model_data.vlonu <= (modelbox[1]+1)))[0]
            else:
                model_data.model_index_lonu = numpy.where(((modelbox[0]-1) <= model_data.vlonu) | (model_data.vlonu <= (modelbox[1]+1)))[0]
            model_data.model_index_latu = numpy.where(((modelbox[2]-1) <= model_data.vlatu) & (model_data.vlatu <= (modelbox[3]+1)))[0]
            model_data.vlonu = model_data.vlonu[model_data.model_index_lonu]
            model_data.vlatu = model_data.vlatu[model_data.model_index_latu]
            if p.lonu != p.lonv:
                if modelbox[0] < modelbox[1]:
                    model_data.model_index_lonv = numpy.where(((modelbox[0]-1) <= model_data.vlonv) & (model_data.vlonv <= (modelbox[1]+1)))[0]
                else:
                    model_data.model_index_lonv = numpy.where(((modelbox[0]-1) <= model_data.vlonv) | (model_data.vlonv <= (modelbox[1]+1)))[0]
            else:
                 model_data.model_index_lonv = model_data.model_index_lonu
            if p.latu != p.latv:
                model_data.model_index_latv = numpy.where(((modelbox[2]-1) <= model_data.vlatv) & (model_data.vlatv <= (modelbox[3]+1)))[0]
            else:
                model_data.model_index_latv = model_data.model_index_latu
            model_data.vlonv = model_data.vlonv[model_data.model_index_lonv]
            model_data.vlatv = model_data.vlatv[model_data.model_index_latv]

        else:
            if modelbox[0] < modelbox[1]:
                model_data.model_indexu = numpy.where(((modelbox[0]-1) <= model_data.vlonu) & (model_data.vlonu <= (modelbox[1]+1)) & ((modelbox[2]-1) <= model_data.vlatu) & (model_data.vlatu <= (modelbox[3]+1)))
                if p.lonu != p.lonv:
                    model_data.model_indexv = numpy.where(((modelbox[0]-1) <= model_data.vlonv) & (model_data.vlonv <= (modelbox[1]+1)) & ((modelbox[2]-1) <= model_data.vlatv) & (model_data.vlatv <= (modelbox[3]+1)))
            else:
                model_data.model_indexu = numpy.where(((modelbox[0]-1) <= model_data.vlonu) | (model_data.vlonu <= (modelbox[1]+1)) & ((modelbox[2]-1) <= model_data.vlatu) & (model_data.vlatu <= (modelbox[3]+1)))
                if p.lonu != p.lonv:
                    model_data.model_indexv = numpy.where(((modelbox[0]-1) <= model_data.vlonv) | (model_data.vlonv <= (modelbox[1]+1)) & ((modelbox[2]-1) <= model_data.vlatv) & (model_data.vlatv <= (modelbox[3]+1)))

        model_data.model = model
        model_data.vloncircu = numpy.rad2deg(numpy.unwrap(model_data.vlonu))
        if p.lonu != p.lonv:
            model_data.vloncircv = numpy.rad2deg(numpy.unwrap(model_data.vlonv))
    if modelbox[1] == 0:
        modelbox[1] = 359.99
    # - Make SKIM grid if necessary """
    if p.makesgrid:
        print('\n Force creation of SKIM grid')
        orb = build_swath.makeorbit(modelbox, p, orbitfile=p.filesat)
        build_swath.orbit2swath(modelbox, p, orb)
        print("\n SKIM Grids and nadir tracks have been written in"\
              "{}".format(p.outdatadir))
        print("-----------------------------------------------")

    ## Initialize errors
    err, errnad = load_error(p)

    # - Compute interpolated SSH and errors for each pass, at each
    #   cycle
    print('Compute interpolated SSH and errors:')
    #   load all SKIM grid files (one for each pass)
    listsgridfile = sorted(glob.glob(p.filesgrid + '_p*.nc'))
    if not listsgridfile:
        print('\n There is no SKIM grid file in ' + p.outdatadir
              + ', run simulator with option makesgrid set to true in your params file')
        sys.exit()
    #   Model time step
    modeltime = numpy.arange(0, p.nstep*p.timestep, p.timestep)
    #   Remove the grid from the list of model files
    if p.file_input:
        list_file.remove(list_file[0])
    #   Initialize progress bar variables
    istep = 0
    ntot = 1
    # - Loop on SKIM grid files
    for sgridfile in listsgridfile:
        #   Load SKIM grid files (Swath and nadir)
        sgrid = load_sgrid(sgridfile, p)
        sgrid_tmp = load_sgrid(sgridfile, p)
        sgrid.gridfile = sgridfile
    #   Select model data around the swath to reduce interpolation cost in
    #   griddata
    # - Generate SKIM like
    #   Compute number of cycles needed to cover all nstep model timesteps
        rcycle = (p.timestep * p.nstep)/float(sgrid.cycle)
        ncycle = int(rcycle)
    #   Loop on all cycles
        for cycle in range(0, ncycle+1):
            if ifile > (p.nstep/p.timestep + 1):
                break
            #   Create SKIM-like
            if p.file_input is None:
                model_data = []
            ur_true_all = []
            u_true_all = []
            v_true_all = []
            err_instr = []
            err_uss = []
            ur_obs = []
            for i in range(len(p.list_pos_12) + len(p.list_pos_6) + 1):
                sgrid_tmp.lon = sgrid.lon[i]
                sgrid_tmp.lat = sgrid.lat[i]
                sgrid_tmp.time = sgrid.time[i]
                if (i > 0) and i < len(p.list_pos_12):
                    pos = p.list_pos_12[i - 1]
                elif i >= len(p.list_pos_12):
                    pos = p.list_pos_6[i - 1 - len(p.list_pos_12)]
                if i == 0:
                    ur_true = numpy.zeros((numpy.shape(sgrid_tmp.lon)))
                    u_true = numpy.zeros((numpy.shape(sgrid_tmp.lon)))
                    v_true = numpy.zeros((numpy.shape(sgrid_tmp.lon)))
                    err.ur_obs = numpy.zeros((numpy.shape(sgrid_tmp.lon)))
                    err.instr = numpy.zeros((numpy.shape(sgrid_tmp.lon)))
                    err.ur_uss = numpy.zeros((numpy.shape(sgrid_tmp.lon)))
                else:
                    ur_true, u_true, v_true, vindice, time, progress = \
                    create_SKIMlikedata(cycle,
                                        numpy.shape(listsgridfile)[0]*rcycle,
                                        list_file, list_file_uss, modelbox,
                                        sgrid_tmp, model_data, modeltime, err,
                                        pos, p, progress_bar=True)
                err_instr.append(err.instr)
                err_uss.append(err.ur_uss)
                ur_true_all.append(ur_true)
                u_true_all.append(u_true)
                v_true_all.append(v_true)
                ### Make error here
                ur_obs.append(err.ur_obs)
            #   Save outputs in a netcdf file
            if (~numpy.isnan(vindice)).any() or not p.file_input:
                save_SKIM(cycle, sgrid, err, p, time=time, vindice=vindice,
                          ur_model=ur_true_all, ur_obs=ur_obs,
                          err_instr=err_instr, err_uss=err_uss,
                          u_model=u_true_all, v_model=v_true_all)
            del time
            # if p.file_input: del index
        if p.file_input:
            model_data.vlonu = (model_data.vlonu + 360) % 360
            if p.lonu != p.lonv:
                model_data.vlonv = (model_data.vlonv + 360) % 360

        modelbox[0] = (modelbox[0] + 360) % 360
        modelbox[1] = (modelbox[1] + 360) % 360
        del sgrid
    if progress != 1:
        progress = mod_tools.update_progress(1,
                                             'All passes have been processed',
                                             '')
    # - Write Selected parameters in a txt file
    rw_data.write_params(p, p.outdatadir + os.sep + 'skim_simulator.output')
    print("\n Simulated skim files have been written in " + p.outdatadir)
    print("----------------------------------------------------------")


def load_error(p):
    '''Initialize random coefficients that are used to compute
    random errors following the specified spectrum. \n
    If a random coefficient file is specified, random coefficients
    are loaded from this file.
    '''
    # import skimsimulator.build_error as build_error
    err = build_error.error(p)
    err.init_error(p)
    errnad = build_error.errornadir(p)
    errnad.init_error(p)
    return err, errnad


def load_sgrid(sgridfile, p):
    '''Load SKIM swath and Nadir data for file sgridfile '''
    # import skimsimulator.rw_data as rw_data

    # Load SKIM swath file
    sgrid = rw_data.Sat_SKIM(file=sgridfile)
    cycle = 0
    x_al = []
    al_cycle = 0
    timeshift = 0
    sgrid.load_swath(p, cycle=cycle ,x_al=x_al, al_cycle=al_cycle,
                     timeshift=timeshift)
    sgrid.loncirc = []
    for i in range(len(sgrid.lon)):
        sgrid.loncirc.append(numpy.rad2deg(numpy.unwrap(sgrid.lon[i])))
    # Extract the pass number from the file name
    ipass = int(sgridfile[-6: -3])
    sgrid.ipass=ipass
    return sgrid


def create_SKIMlikedata(cycle, ntotfile, list_file, list_file_uss, modelbox, sgrid,
                        model_data, modeltime, err, pos, p,
                        progress_bar=True):
    '''Create SKIM and nadir errors err and errnad, interpolate model SSH model_data on swath and nadir track,
    compute SKIM-like and nadir-like data for cycle, SKIM grid sgrid and ngrid. '''
    # - Progress bar variables are global
    global istep
    global ntot
    #   Initialiaze errors and SSH
    progress = 0
    #inclination = p.inclination #* math.pi / 180
    #omega = p.rotation_speed * numpy.pi * 2. / 60.
    err.instr = numpy.zeros((numpy.shape(sgrid.lon)[0]))
    err.ur_uss = numpy.zeros((numpy.shape(sgrid.lon)[0]))
    err.wet_tropo1nadir = numpy.zeros((numpy.shape(sgrid.lon)[0]))
    err.wet_tropo2nadir = numpy.zeros((numpy.shape(sgrid.lon)[0]))
    err.wtnadir = numpy.zeros((numpy.shape(sgrid.lon)[0]))
    err.nadir = numpy.zeros((numpy.shape(sgrid.lon)[0]))
    date1 = cycle * sgrid.cycle
    u_true = (numpy.zeros((numpy.shape(sgrid.lon)[0])))
    v_true = (numpy.zeros((numpy.shape(sgrid.lon)[0])))
    u_uss = (numpy.zeros((numpy.shape(sgrid.lon)[0])))
    v_uss = (numpy.zeros((numpy.shape(sgrid.lon)[0])))
    ur_true = (numpy.zeros((numpy.shape(sgrid.lon)[0])))
    vindice = numpy.zeros((numpy.shape(sgrid.lon)[0])) * numpy.nan
    # Definition of the time in the model
    time = sgrid.time / 86400. + date1
    lon = sgrid.lon
    lat = sgrid.lat
    sgrid.timeshift /= 86400
    # Look for satellite data that are beween step-p.timesetp/2 end setp+p.step/2
    if p.file_input is not None:
        index_filemodel = numpy.where(((time[-1]-sgrid.timeshift) >= (modeltime-p.timestep/2.))
                                      & ((time[0]-sgrid.timeshift) < (modeltime+p.timestep/2.)))  # [0]
        nfile=0
        # At each step, look for the corresponding time in the satellite data
        for ifile in index_filemodel[0]:
            progress = mod_tools.update_progress(float(istep)/float(ntot
                                                 * ntotfile), 'pass: '
                                                 + str(sgrid.ipass),
                                                 'model file: '+ str(ifile)
                                                 + ', cycle:'+str(cycle+1))
                # If there are satellite data, Get true SSH from model
            if numpy.shape(index_filemodel)[1] > 0:
                # number of file to be processed used in the progress bar
                ntot = ntot + numpy.shape(index_filemodel)[1]-1
                # if numpy.shape(index)[1]>1:
                # Select part of the track that corresponds to the time of the model (+-timestep/2)
                ind_time = numpy.where(((time-sgrid.timeshift) >= (modeltime[ifile]-p.timestep/2.)) & ((time-sgrid.timeshift) < (modeltime[ifile]+p.timestep/2.)) )
            # Load data from this model file
            #import pdb ; pdb.set_trace()
            if model_data.model == 'WW3':
                filetime = ifile - p.dim_time * nfile
                if filetime >= (p.dim_time * (nfile + 1)):
                    nfile += 1
                model_step = rw_data.WW3(file=p.indatadir+os.sep+list_file[0], varu=p.varu, varv=p.varv, time=filetime)
                if p.uss is True:
                    uss_step =  rw_data.WW3(file=p.indatadir+os.sep+list_file_uss[0], varu='uuss', varv='vuss', time=filetime)
            else:
                model_step = eval('rw_data.' + model_data.model
                              + '(file=p.indatadir+os.sep+list_file[ifile], varu=p.varu, varv=p.varv)')
            if p.grid == 'regular':
                model_step.read_var()
                u_model = model_step.vvaru[model_data.model_index_latu, :]
                u_model = u_model[:, model_data.model_index_lonu]
                v_model = model_step.vvarv[model_data.model_index_latv, :]
                v_model = v_model[:, model_data.model_index_lonv]
                if p.uss is True:
                    uss_step.read_var()
                    u_uss_mod = uss_step.vvaru[model_data.model_index_latu, :]
                    u_uss_mod = u_uss_mod[:, model_data.model_index_lonu]
                    v_uss_mod = uss_step.vvaru[model_data.model_index_latv, :]
                    v_uss_mod = v_uss_mod[:, model_data.model_index_lonv]
            else:
                model_step.read_var(index=None)
                u_model = model_step.vvaru[model_data.model_indexu]
                v_model = model_step.vvarv[model_data.model_indexv]
                if p.uss is True:
                    uss_step.read_var(index=None)
                    u_uss_mod = uss_step.vvaru[model_data.model_indexu]
                    v_uss_mod = uss_step.vvarv[model_data.model_indexv]
            # - Interpolate Model data on a SKIM grid and/or along the
            #   nadir track
            # if grid is regular, use interpolate.RectBivariateSpline to
            # interpolate
            if p.grid == 'regular' and \
                    len(numpy.shape(model_data.vlonu)) == 1:
                # ########################TODO
                # Flatten satellite grid and select part of the track
                # corresponding to the model time
                Teval = interpolate.RectBivariateSpline(model_data.vlatu, model_data.vlonu, numpy.isnan(u_model), kx=1, ky=1, s=0).ev(lat[ind_time[0]], lon[ind_time[0]])
                u_model_mask = + u_model
                u_model_mask[numpy.isnan(u_model_mask)] = 0.
                u_true_ind_time = interpolate.RectBivariateSpline(model_data.vlatu, model_data.vlonu, u_model_mask, kx=1, ky=1, s=0).ev(lat[ind_time[0]], lon[ind_time[0]])
                u_true_ind_time[Teval > 0] = numpy.nan
                u_true[ind_time[0]] = u_true_ind_time #[ind_time[0]]
                if p.uss is True:
                    u_uss_mod_mask = + u_uss_mod
                    u_uss_mod_mask[numpy.isnan(u_uss_mod_mask)] = 0.
                    u_uss_ind_time = interpolate.RectBivariateSpline(model_data.vlatu, model_data.vlonu, u_uss_mod_mask, kx=1, ky=1, s=0).ev(lat[ind_time[0]], lon[ind_time[0]])
                    u_uss_ind_time[Teval > 0] = numpy.nan
                    u_uss[ind_time[0]] = u_uss_ind_time #[ind_time[0]]
                Teval = interpolate.RectBivariateSpline(model_data.vlatv, model_data.vlonv, numpy.isnan(v_model), kx=1, ky=1, s=0).ev(lat[ind_time[0]], lon[ind_time[0]])
                v_model_mask = + v_model
                v_model_mask[numpy.isnan(v_model_mask)] = 0.
                v_true_ind_time = interpolate.RectBivariateSpline(model_data.vlatv, model_data.vlonv, v_model_mask, kx=1, ky=1, s=0).ev(lat[ind_time[0]], lon[ind_time[0]])
                v_true_ind_time[Teval > 0] = numpy.nan
                v_true[ind_time[0]] = v_true_ind_time #[ind_time[0]]
                if p.uss is True:
                    v_uss_mod_mask = + v_uss_mod
                    v_uss_mod_mask[numpy.isnan(v_uss_mod_mask)] = 0.
                    v_uss_ind_time = interpolate.RectBivariateSpline(model_data.vlatv, model_data.vlonv, v_uss_mod_mask, kx=1, ky=1, s=0).ev(lat[ind_time[0]], lon[ind_time[0]])
                    v_uss_ind_time[Teval > 0] = numpy.nan
                    v_uss[ind_time[0]] = v_uss_ind_time #[ind_time[0]]
            else:
                # Grid is irregular, interpolation can be done using
                # pyresample module if it is installed or griddata
                # function from scipy.
                # Note that griddata is slower than pyresample functions.

                try:
                    import pyresample as pr
                    model_data.vlon = pr.utils.wrap_longitudes(model_data.vlon)
                    lon = pr.utils.wrap_longitudes(lon)
                    if len(numpy.shape(model_data.vlon)) <= 1:
                        model_data.vlon, model_data.vlat = numpy.meshgrid(model_data.vlon, model_data.vlat)
                    swath_defu = pr.geometry.SwathDefinition(lons=model_data.vlonu, lats=model_data.vlatu)
                    swath_defv = pr.geometry.SwathDefinition(lons=model_data.vlonv, lats=model_data.vlatv)
                    grid_def = pr.geometry.SwathDefinition(lons=lon,
                                                           lats=lat)
                    if p.interpolation=='nearest':
                        u_true[ind_time[0]] = pr.kd_tree.resample_nearest(swath_defu, u_model, grid_def, radius_of_influence=max(p.delta_al)*10**3, epsilon=100)
                        v_true[ind_time[0]] = pr.kd_tree.resample_nearest(swath_defv, v_model, grid_def, radius_of_influence=max(p.delta_al)*10**3, epsilon=100)
                        if p.uss is True:
                            u_uss[ind_time[0]] = pr.kd_tree.resample_nearest(swath_defu, u_uss_mod, grid_def, radius_of_influence=(p.delta_al)*10**3, epsilon=100)
                            v_uss[ind_time[0]] = pr.kd_tree.resample_nearest(swath_defv, v_uss_mod, grid_def, radius_of_influence=max(p.delta_al)*10**3, epsilon=100)
                    else:
                        u_true[ind_time[0]] = pr.kd_tree.resample_gauss(swath_defu, u_model, grid_def, radius_of_influence=3*max(p.delta_al)*10**3, sigmas=max(p.delta_al, p.delta_ac)*10**3, fill_value=None)
                        v_true[ind_time[0]] = pr.kd_tree.resample_gauss(swath_defv, v_model, grid_def, radius_of_influence=3*max(p.delta_al)*10**3, sigmas=max(p.delta_al, p.delta_ac)*10**3, fill_value=None)
                        if p.uss is True:
                            u_uss[ind_time[0]] = pr.kd_tree.resample_gauss(swath_defu, u_uss_mod, grid_def, radius_of_influence=3*max(p.delta_al)*10**3, sigmas=max(p.delta_al, p.delta_ac)*10**3, fill_value=None)
                            v_uss[ind_time[0]] = pr.kd_tree.resample_gauss(swath_defv, v_uss_mod, grid_def, radius_of_influence=3*max(p.delta_al)*10**3, sigmas=max(p.delta_al, p.delta_ac)*10**3, fill_value=None)
                except:
                    u_true[ind_time[0]] = interpolate.griddata((model_data.vlonu.ravel(), model_data.vlatu.ravel()), u_model.ravel(), (lon[ind_time[0]], lat[ind_time[0]]), method=p.interpolation)
                    v_true[ind_time[0]] = interpolate.griddata((model_data.vlonv.ravel(), model_data.vlatv.ravel()), v_model.ravel(), (lon[ind_time[0]], lat[ind_time[0]]), method=p.interpolation)
                    if p.uss is True:
                        u_uss[ind_time[0]] = interpolate.griddata((model_data.vlonu.ravel(), model_data.vlatu.ravel()), u_uss_mod.ravel(), (lon[ind_time[0]], lat[ind_time[0]]), method=p.interpolation)
                        v_uss[ind_time[0]] = interpolate.griddata((model_data.vlonv.ravel(), model_data.vlatv.ravel()), v_uss_mod.ravel(), (lon[ind_time[0]], lat[ind_time[0]]), method=p.interpolation)
                    if p.interpolation == 'nearest':
                        if modelbox[0] > modelbox[1]:
                            u_true[numpy.where(((lon < modelbox[0])
                                     & (lon > modelbox[1]))
                                     | (lat < modelbox[2])
                                     | (lat > modelbox[3]))] = numpy.nan
                            v_true[numpy.where(((lon < modelbox[0])
                                     & (lon > modelbox[1]))
                                     | (lat < modelbox[2])
                                     | (lat > modelbox[3]))] = numpy.nan
                        else:
                            u_true[numpy.where((lon < modelbox[0])
                                    | (lon > modelbox[1])
                                    | (lat < modelbox[2])
                                    | (lat > modelbox[3]))] = numpy.nan
                            v_true[numpy.where((lon < modelbox[0])
                                    | (lon > modelbox[1])
                                    | (lat < modelbox[2])
                                    | (lat > modelbox[3]))] = numpy.nan

            vindice[ind_time[0]] = ifile
            # del u_true, v_true, model_step
        istep += 1
    else:
        istep += 1
        progress = mod_tools.update_progress(float(istep)/float(ntotfile*ntot),
                                             'pass: ' + str(sgrid.ipass),
                                             'no model file provided'
                                             + ', cycle:' + str(cycle+1))
    try:
       ur_true = mod_tools.proj_radial(u_true, v_true, time, pos,
                                       sgrid.incl, p)
    except: import pdb; pdb.set_trace()
    if p.uss is not True:
        u_uss = None
        v_uss = None
    err.make_error(ur_true, time, pos, p, sgrid.incl,
                   uss=(u_uss, v_uss))
    err.make_vel_error(ur_true, p)
    # if p.file_input: del ind_time, SSH_model, model_step
    return ur_true, u_true, v_true, vindice, time, progress


def save_SKIM(cycle, sgrid, err, p, time=(), vindice=(), ur_model=(),
              ur_obs=(), err_instr=(), err_uss=(), u_model=(), v_model=()):
    file_output = (p.file_output + '_c' + str(cycle+1).zfill(2) + '_p'
                   + str(sgrid.ipass).zfill(3) + '.nc')
    OutputSKIM = rw_data.Sat_SKIM(file=file_output, lon=sgrid.lon,
                                  lat=sgrid.lat, time=sgrid.time,
                                  x_al=sgrid.x_al, cycle=sgrid.cycle)
    OutputSKIM.gridfile = sgrid.gridfile
    OutputSKIM.write_data(p, ur_model=ur_model, index=[vindice,],
                          uss_err=err_uss,
                          nadir_err=[err.nadir, ], ur_obs=ur_obs,
                          instr=err_instr, u_model=u_model, v_model=v_model)
    return None

