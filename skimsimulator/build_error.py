# import params as p
import numpy
import scipy
import skimsimulator.rw_data as rw_data
import skimsimulator.const as const
import skimsimulator.mod_tools as mod_tools
from math import pi, sqrt
from scipy.io import netcdf as nc
import sys
from scipy.ndimage.filters import gaussian_filter
import skimsimulator.mod_tools as mod_tools
import skimsimulator.rw_data as rw_data


class error():
    '''Class error define all the possible errors that can be computed using
    SKIM simulator.
    Random realisation of errors can be initialized using init_error.
    If the random realisations have already been computed and stored in
    file file_coeff, the random realisations are read directly using
    load_coeff.  The corresponding errors on a swath can be computed using
    make_error. '''
    def __init__(self, p,
                 instr=None,
                 uss=None,
                 ):
        self.instr = instr
        self.ur_uss = uss

    def init_error(self, p):
        '''Initialization of errors: Random realisation of errors are
        computed using a known power spectrum.
        The outputs are the amplitude, the phase and the frequency of
        each random realisation.
        By default, there are ncomp1d=2000 random realisations for
        the instrumental errors (1d spectrum) and ncomp2d=2000 random
        realisations for the geophysical errors (2d spectrum)
        and nrandkarin*x_ac km of random number for KaRIN noise.'''
        # # Run reprodictible: generate or load nrand random numbers:
        # - If one instrumental error needs to be computed, load frequencies
        #   and spectrum from the instrumental netcdf file:
        pass


    def make_error(self, u_true, time, pos, p, incl, uss=(None, None)):
        ''' Build errors corresponding to each selected noise
        among the effect of the wet_tropo, the phase between the two signals,
        the timing error, the roll of the satellite, the sea surface bias,
        the distorsion of the mast,
        the karin noise due to the sensor itself. '''
        # - Errors array are the same size as the swath size
        nal = numpy.shape(u_true)
        # ind_al=numpy.arange(0,nal)
        if p.instr is True:
            self.instr = numpy.random.normal(0.0, p.rms_instr, nal)
        if p.uss is True and uss[0] is not None and uss[1] is not None:
            self.ur_uss = mod_tools.proj_radial(uss[0], uss[1],
                                                time, pos, incl, p)
            self.ur_uss *= p.G
        return None

    def make_vel_error(self, ur_true, p):
        '''Compute observed SSH adding all the computed error to the model SSH.
        '''
        numpy.seterr(invalid='ignore')
        self.ur_obs = ur_true
        if p.instr is True:
            self.ur_obs = self.ur_obs + self.instr
        if p.uss is True:
            self.ur_obs = self.ur_obs + self.ur_uss
        if p.file_input is not None:
            self.ur_obs[numpy.where(ur_true == p.model_nan)] = p.model_nan


class errornadir():
    '''Class errornadir defines the error on the nadir.
    Random realisation of errors can be initialized using init_error.
    The correspondg errors on a swath can be computed using make_error. '''
    def __init__(self,p,
                 nadir=None,
                 wet_tropo1=None,
                 wt=None):
        self.nadir = nadir
        self.wet_tropo1 = wet_tropo1
        self.wt = wt
        try:
            self.nrand = p.nrandkarin
        except:
            self.nrand = 1000
            p.nrandkarin = 1000
        try:
            self.ncomp2d = p.ncomp2d
        except:
            self.ncomp2d = 2000
            p.ncomp2d = 2000
        try:
            self.ncomp1d = p.ncomp1d
        except:
            self.ncomp1d = 2000
            p.ncomp1d = 2000

    def init_error(self, p):
        '''Initialization of errors: Random realisation of errors are computed
        using a known power spectrum.
        The outputs are the amplitude, the phase and the frequency of each
        random realisation.
        By default, there are ncomp2d=2000 random realisations for the
        wet tropo and ncomp1d=2000 random realisations for the nadir 1d
        spectrum error.'''
        # Run reprodictible: generate or load nrand random numbers:
        # - Compute random coefficients in 1D for the nadir error
        try:
            wnoise = p.wnoise
        except:
            wnoise = 100
            p.wnoise = 100
        # - Define the sepctrum of the nadir instrument error
        # self.A=numpy.random.normal(0.0,sqrt(p.wnoise)
        # /numpy.float64(sqrt(2*p.delta_al)), (self.nrand))*0.01
        p.delta_al = 1
        f = numpy.arange(1./3000., 1./float(2.*p.delta_al), 1./3000.)
        PSD=8 + 1.05 * 10**(-4) * f**(-2.2)
        indf = numpy.where(f < 0.00023627939582672978)
        PSD[indf] = 10**4
        # Convert spectrum in m2/cy
        PSD = PSD * 10**(-4)
        self.A, self.phi, self.f = mod_tools.gen_coeff_signal1d(f, PSD,
                                                                self.ncomp1d)
        if p.wet_tropo:
            # - Define power spectrum of error in path delay due to wet tropo
            f = numpy.arange(1./3000., 1./float(2.*p.delta_al), 1./3000.)
            # - Global mean wet tropo power spectrum in cm**2/(cycle/km)
            #   for L >= 100 km
            PSwt = 3.156 * 10**-5 * f**(-8./3.)  # *10**-4
            # - Wet tropo power spectrum in cm**2/(cycle/km) for L < 100 km
            indf = numpy.where(f > 10**-2)
            PSwt[indf] = 1.4875 * 10**-4 * f[indf]**(-2.33)  # *10**-4
            # - Compute random coefficients in 2D using the previously defined
            #   power spectrum
            self.A_wt, self.phi_wt, self.frx_wt, self.fry_wt = mod_tools.gen_coeff_signal2d(f, PSwt , self.ncomp2d)
            # - Define radiometer error power spectrum for a beam
            #   High frequencies are cut to filter the associated error during
            #   the reconstruction of the wet trop signal
            # f=numpy.arange(1./3000.,1./float(20.),1./3000.)
            PSradio = 9.5 * 10**-5 * f**-1.79
            PSradio[numpy.where((f < 1./1000.))] = 9.5 * 10**-5*(10**-3)**-1.79
            indf = numpy.where((f > 0.0023) & (f <= 0.0683))
            PSradio[indf] = 0.036 * f[indf]**-0.814
            PSradio[numpy.where(f > 0.0683)] = 0.32
            # - Compute random coefficients (1D) for the radiometer error power
            #   spectrum for right and left beams
            self.A_radio, self.phi_radio, self.fr_radio = mod_tools.gen_coeff_signal1d(f, PSradio , self.ncomp2d)
        return None

    def load_coeff(self, p):
        '''Load existing random realisations that has been stored in
        nadir+file_coeff. The outputs are the amplitude,
        the phase and the frequency of each random realisation.
        There are ncomp random realisations.'''
        try:
            fid = nc.netcdf_file(p.file_coeff[:-3] + '_nadir.nc', 'r')
        except:
            print('There was an error opening the file nadir'
                  + p.file_coeff[:-3] + '_nadir.nc')
            sys.exit()
        self.A = numpy.array(fid.variables['A'][:]).squeeze()
        self.f = numpy.array(fid.variables['f'][:]).squeeze()
        self.phi = numpy.array(fid.variables['phi'][:]).squeeze()
        if p.wet_tropo:
            self.A_wt = numpy.array(fid.variables['A_wt'][:]).squeeze()
            if numpy.shape(self.A_wt)[0] != self.ncomp2d:
              sys.exit(p.file_coeff + ' dimensions are different from ncomp2d='
                       + str(self.ncomp2d) + '\n remove ' + p.file_coeff
                       + ' or adjust ncomp2d number in the parameter file')
            self.phi_wt = numpy.array(fid.variables['phi_wt'][:]).squeeze()
            self.frx_wt = numpy.array(fid.variables['frx_wt'][:]).squeeze()
            self.fry_wt = numpy.array(fid.variables['fry_wt'][:]).squeeze()
            self.A_radio = numpy.array(fid.variables['A_radio'][:]).squeeze()
            self.phi_radio = numpy.array(fid.variables['phi_radio'][:]).squeeze()
            self.fr_radio = numpy.array(fid.variables['fr_radio'][:]).squeeze()
        fid.close()
        return None

    def make_error(self, orb, cycle, SSH_true,p):
        nal = numpy.shape(SSH_true)[0]
        errnadir = numpy.zeros((nal))
        # - Compute random noise of 10**2 cm**2/(km/cycle)
        # - Compute the correspond error on the nadir in m
        for comp in range(0, self.ncomp1d):
            phase_x_al = (2. * pi * float(self.f[comp])
                          *(numpy.float64(orb.x_al[:])
                          + float(cycle*orb.al_cycle))) % (2.*pi)
            errnadir[:] = (errnadir[:] + 2*self.A[comp]
                           *numpy.cos(phase_x_al[:]+self.phi[comp]))
        # - Compute the correspond timing error on the swath in m
        self.nadir = errnadir[:]
        if p.wet_tropo:
            # - Initialization of radiometer error in right and left beam
            err_radio = numpy.zeros((nal))
            # - Initialization of swath matrices and large swath matrices
            #   (which include wet tropo data around the nadir and outside the
            #   swath)
            #   x_ac_large and wt_large are necessary to compute the gaussian
            #   footprint of a beam on the nadir or near the edge of the swath
            x_ac_large = numpy.arange(-2. * p.sigma/float(p.delta_ac),
                                      2.*p.sigma/float(p.delta_ac)+p.delta_ac,
                                      p.delta_ac)
            wt_large = numpy.zeros((numpy.shape(orb.x_al[:])[0],
                                   numpy.shape(x_ac_large)[0]))
            x_large,y_large = numpy.meshgrid(x_ac_large,
                                             numpy.float64(orb.x_al[:])
                                             + float(cycle*orb.al_cycle))
            # - Compute path delay error due to wet tropo and radiometer error
            #   using random coefficient initialized with power spectrums
            for comp in range(0, self.ncomp2d):
                phase_x_al_large = (2. * pi * (float(self.frx_wt[comp])
                                    * (numpy.float64(x_large))
                                    + float(self.fry_wt[comp])
                                    * numpy.float64(y_large))) % (2.*pi)
                wt_large = (wt_large + self.A_wt[comp]
                            * numpy.cos(phase_x_al_large
                            + self.phi_wt[comp])*10**-2)
                phase_x_al = (2. * pi * float(self.fr_radio[comp])
                              * (numpy.float64(orb.x_al[:])
                              + float(cycle*orb.al_cycle))) % (2.*pi)
                err_radio = (err_radio + 2*self.A_radio[comp]
                             *numpy.cos(phase_x_al+self.phi_radio[comp])*10**-2)
            # - Compute Residual path delay error after a 1-beam radiometer
            #   correction
            beam = numpy.zeros((nal))
            diff_h1 = numpy.zeros((nal))
            # indac=numpy.where((sgrid.x_ac<2.*p.sigma)
            # & (sgrid.x_ac>-2.*p.sigma))[0]
            # - Find across track indices in the gaussian footprint of
            #   2.*p.sigma
            indac = numpy.where((x_ac_large < 2.*p.sigma)
                                & (x_ac_large > -2.*p.sigma))[0]
            for i in range(0, nal):
            # - Find along track indices in the gaussian footprint of
            #   2.*p.sigma
                indal = numpy.where(((orb.x_al[:]-orb.x_al[i]) <= (2*p.sigma))
                                    & ((orb.x_al[:]-orb.x_al[i]) > -2*p.sigma))[0]
                # indal=numpy.where(((sgrid.x_al[:]-sgrid.x_al[i])<=
                # (2*p.sigma)) & ((sgrid.x_al[:]-sgrid.x_al[i])>(-2*p.sigma)))[0]
                x, y = numpy.meshgrid(x_ac_large[min(indac): max(indac)+1],
                                      (orb.x_al[(min(indal)): (max(indal)+1)]-orb.x_al[i]))
                # - Compute path delay on gaussian footprint
                G = 1. / (2.*pi*p.sigma**2) * numpy.exp(-(x**2.+y**2.)
                                                        / (2.*p.sigma**2))
                beam[i] = sum(sum(G*wt_large[min(indal): max(indal)+1,
                              min(indac):max(indac)+1]))/sum(sum(G))+err_radio[i]
            # - Filtering beam signal to cut frequencies higher than 125 km
            beam = gaussian_filter(beam, 30./p.delta_al)
            # - Compute residual path delay
            diff_h1 = wt_large[:, int(numpy.shape(wt_large)[1]/2.)] - beam
            self.wet_tropo1 = diff_h1
            self.wt = wt_large[:, int(numpy.shape(wt_large)[1]/2.)]
        return None

    def save_coeff(self, p):
        '''Save random realisations to enable runs to be reproducible.
        The ncomp1d random phase phi, amplitude A and frequency fr for
        1D spectrum and ncomp2d random phase phi, amplitude A and frequencies
        frx and fry for 2D spectrum are saved in nadirfile_coeff for each error
        and can be loaded using load_coeff.
        '''
        # - Open Netcdf file in write mode
        fid = nc.netcdf_file(p.file_coeff[:-3] + '_nadir.nc', 'w')
        fid.description = "Random coefficients from orbit simulator"

## - Create dimensions
        fid.createDimension('nrand1d', self.ncomp1d)
        fid.createDimension('nrand2d', self.ncomp2d)
## - Create and write Variables
        var = fid.createVariable('A', 'f4', ('nrand1d', ))
        var[:] = self.A
        var = fid.createVariable('f', 'f4', ('nrand1d', ))
        var[:] = self.f
        var = fid.createVariable('phi', 'f4', ('nrand1d', ))
        var[:] = self.phi

        # var = fid.createVariable('phi', 'f4', ('ninstr',))
        # var[:] = self.phi
        # var = fid.createVariable('fr', 'f4', ('ninstr',))
        # var[:] = self.fr
        if p.wet_tropo:
            var = fid.createVariable('A_wt', 'f4', ('nrand2d', ))
            var[:] = self.A_wt
            var = fid.createVariable('phi_wt', 'f4', ('nrand2d', ))
            var[:] = self.phi_wt
            var = fid.createVariable('frx_wt', 'f4', ('nrand2d', ))
            var[:] = self.frx_wt
            var = fid.createVariable('fry_wt', 'f4', ('nrand2d', ))
            var[:] = self.fry_wt
            var = fid.createVariable('A_radio', 'f4', ('nrand2d', ))
            var[:] = self.A_radio
            var = fid.createVariable('phi_radio', 'f4', ('nrand2d', ))
            var[:] = self.phi_radio
            var = fid.createVariable('fr_radio', 'f4', ('nrand2d', ))
            var[:] = self.fr_radio
        fid.close()
        return None
