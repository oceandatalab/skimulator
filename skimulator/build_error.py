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
import sys
import math
import numpy
from scipy.io import netcdf as nc
from scipy.ndimage.filters import gaussian_filter
import skimulator.mod_tools as mod_tools
import logging
logger = logging.getLogger(__name__)


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
                 err_uss=None,
                 ):
        self.instr = instr
        self.ur_uss = uss
        self.err_uss = err_uss

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

    def make_error(self, u_true, p, radial_angle, Gvar, file_rms_instr,
                   uss=(None, None), std_local=None, errdcos=None):
        ''' Build errors corresponding to each selected noise
        among the effect of the wet_tropo, the phase between the two signals,
        the timing error, the roll of the satellite, the sea surface bias,
        the distorsion of the mast,
        the karin noise due to the sensor itself. '''
        # - Load instrumental rms file into dictionnary
        theta, rms = numpy.loadtxt(file_rms_instr, usecols=(0, 1), unpack=True)
        dic_rms = {}
        with open(file_rms_instr) as f:
            for line in f:
                (key, val) = line.split()
                dic_rms[int(key)] = float(val)

        # - Compute rms corresponding to each radial angle
        # Refer radial angle to along track direction
        radial_angle_along = numpy.degrees((radial_angle + math.pi / 2)
                                           % math.pi)
        # Interpolate radial_angle with theta values
        step = theta[1]-theta[0]
        if((theta[1:]-theta[:-1]) != step).any():
            logger.error('instrumental rms file noise has not a constant step')
            sys.exit(1)
        thetaprev = numpy.floor(radial_angle_along/step)
        fac = radial_angle_along - thetaprev
        rms_theta = radial_angle * 0
        i = 0
        for itheta, ifac in zip(thetaprev, fac):
            keytheta = list(dic_rms.keys())[int(itheta)]
            keythetanext = list(dic_rms.keys())[int(itheta + step) % 180]
            rms_theta[i] = (dic_rms[keytheta] * ifac / step
                            + (step - ifac) / step * dic_rms[keythetanext])
            i += 1
        # - Errors array are the same size as the swath size
        # nal = numpy.shape(u_true)
        # ind_al=numpy.arange(0, nal)
        if p.instr is True:
            self.instr = numpy.random.normal(0.0 * rms_theta, rms_theta
                                             * p.rms_instr_factor * 10**(-2))
        if p.uss is True and uss[0] is not None and uss[1] is not None:
            self.ur_uss = mod_tools.proj_radial(uss[0], uss[1],
                                                radial_angle)

            self.ur_uss *= Gvar
        if p.uss is True and std_local is not None:
            if p.formula is True:
                if errdcos is None:
                    errdcos = 1.
                err_uss_tmp = p.bias_std * std_local * Gvar * errdcos
                self.err_uss = numpy.random.normal(0.0 * std_local,
                                                   err_uss_tmp)
            self.std_uss = std_local
        return None

    def make_vel_error(self, ur_true, p):
        '''Compute observed velocity adding all the computed error to the model
        velocity.
        '''
        numpy.seterr(invalid='ignore')
        self.ur_obs = + ur_true
        if p.instr is True:
            self.ur_obs = self.ur_obs + self.instr
        if p.uss is True:
            self.ur_obs = self.ur_obs + self.err_uss
        if p.file_input is not None:
            self.ur_obs[numpy.where(ur_true == p.model_nan)] = p.model_nan


class errornadir():
    '''Class errornadir defines the error on the nadir.
    Random realisation of errors can be initialized using init_error.
    The correspondg errors on a swath can be computed using make_error. '''
    def __init__(self, p, nadir=None, wet_tropo1=None, wt=None):
        self.nadir = nadir
        self.wet_tropo1 = wet_tropo1
        self.wt = wt
        self.nrand = getattr(p, 'nrandkarin', 1000)
        p.nrandkarin = self.nrand
        self.ncomp2d = getattr(p, 'ncomp2d', 2000)
        p.ncomp2d = self.ncomp2d
        self.ncomp1d = getattr(p, 'ncomp1d', 2000)
        p.ncomp1d = self.ncomp1d

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
        wnoise = getattr(p, 'wnoise', 100)
        p.wnoise = wnoise
        # - Define the sepctrum of the nadir instrument error
        # self.A=numpy.random.normal(0.0,sqrt(p.wnoise)
        # /numpy.float64(sqrt(2*p.delta_al)), (self.nrand))*0.01
        p.delta_al = 1
        f = numpy.arange(1./3000., 1./float(2.*p.delta_al), 1./3000.)
        PSD = 8 + 1.05 * 10**(-4) * f**(-2.2)
        indf = numpy.where(f < 0.00023627939582672978)
        PSD[indf] = 10**4
        # Convert spectrum in m2/cy
        PSD = PSD * 10**(-4)
        self.A, self.phi, self.f = mod_tools.gen_coeff_signal1d(f, PSD,
                                                                self.ncomp1d)
        if p.wet_tropo is True:
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
            gen_coeff = mod_tools.gen_coeff_signal2d(f, PSwt, self.ncomp2d)
            self.A_wt, self.phi_wt, self.frx_wt, self.fry_wt = gen_coeff
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
            gen_coeff = mod_tools.gen_coeff_signal1d(f, PSradio, self.ncomp2d)
            self.A_radio, self.phi_radio, self.fr_radio = gen_coeff
        return None

    def load_coeff(self, p):
        '''Load existing random realisations that has been stored in
        nadir+file_coeff. The outputs are the amplitude,
        the phase and the frequency of each random realisation.
        There are ncomp random realisations.'''
        try:
            fid = nc.netcdf_file(p.file_coeff[:-3] + '_nadir.nc', 'r')
        except (FileNotFoundError, IOError):
            logger.error('There was an error opening the file nadir'
                         '{}_nadir.nc'.format(p.file_coeff[:-3]))
            sys.exit(1)
        self.A = numpy.array(fid.variables['A'][:]).squeeze()
        self.f = numpy.array(fid.variables['f'][:]).squeeze()
        self.phi = numpy.array(fid.variables['phi'][:]).squeeze()
        if p.wet_tropo is True:
            self.A_wt = numpy.array(fid.variables['A_wt'][:]).squeeze()
            if numpy.shape(self.A_wt)[0] != self.ncomp2d:
                logger.error('{} dimensions are different from ncomp2d = {} \n'
                             'remove {} or adjustncomp2d number in the '
                             'parameter file'.format(p.file_coeff,
                                                     self.ncomp2d,
                                                     p.file_coeff))
                sys.exit(1)
            self.phi_wt = numpy.array(fid.variables['phi_wt'][:]).squeeze()
            self.frx_wt = numpy.array(fid.variables['frx_wt'][:]).squeeze()
            self.fry_wt = numpy.array(fid.variables['fry_wt'][:]).squeeze()
            self.A_radio = numpy.array(fid.variables['A_radio'][:]).squeeze()
            self.phi_radio = numpy.array(fid.variables['phi_radio']
                                         [:]).squeeze()
            self.fr_radio = numpy.array(fid.variables['fr_radio'][:]).squeeze()
        fid.close()
        return None

    def make_error(self, orb, cycle, SSH_true, p):
        nal = numpy.shape(SSH_true)[0]
        errnadir = numpy.zeros((nal))
        # - Compute random noise of 10**2 cm**2/(km/cycle)
        # - Compute the correspond error on the nadir in m
        for comp in range(0, self.ncomp1d):
            phase_x_al = (2. * math.pi * float(self.f[comp])
                          * (numpy.float64(orb.x_al[:])
                          + float(cycle*orb.al_cycle))) % (2.*math.pi)
            errnadir[:] = (errnadir[:] + 2*self.A[comp]
                           * numpy.cos(phase_x_al[:] + self.phi[comp]))
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
            x_large, y_large = numpy.meshgrid(x_ac_large,
                                              numpy.float64(orb.x_al[:])
                                              + float(cycle*orb.al_cycle))
            # - Compute path delay error due to wet tropo and radiometer error
            #   using random coefficient initialized with power spectrums
            for comp in range(0, self.ncomp2d):
                phase_x_al_large = (2. * math.pi * (float(self.frx_wt[comp])
                                    * (numpy.float64(x_large))
                                    + float(self.fry_wt[comp])
                                    * numpy.float64(y_large))) % (2.*math.pi)
                wt_large = (wt_large + self.A_wt[comp]
                            * numpy.cos(phase_x_al_large
                            + self.phi_wt[comp])*10**-2)
                phase_x_al = (2. * math.pi * float(self.fr_radio[comp])
                              * (numpy.float64(orb.x_al[:])
                              + float(cycle*orb.al_cycle))) % (2.*math.pi)
                err_radio = (err_radio + 2*self.A_radio[comp]
                             * numpy.cos(phase_x_al + self.phi_radio[comp])
                             * 10**-2)
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
                x, y = numpy.meshgrid(x_ac_large[min(indac): max(indac)+1],
                                      (orb.x_al[(min(indal)):
                                                (max(indal)+1)]-orb.x_al[i]))
                # - Compute path delay on gaussian footprint
                G = 1. / (2.*math.pi*p.sigma**2) * numpy.exp(-(x**2.+y**2.)
                                                             / (2.*p.sigma**2))
                beam[i] = (sum(sum(G*wt_large[min(indal): max(indal)+1,
                                              min(indac):max(indac)+1]))
                           / sum(sum(G))+err_radio[i])
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

        # - Create dimensions
        fid.createDimension('nrand1d', self.ncomp1d)
        fid.createDimension('nrand2d', self.ncomp2d)
        # - Create and write Variables
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
        if p.wet_tropo is True:
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


def make_vel_error(ur_true, p, instr=None, err_uss=None):
        '''Compute observed velocity adding all the computed error to the model
        velocity.
        '''
        numpy.seterr(invalid='ignore')
        ur_obs = + ur_true
        if p.instr is True and instr is not None:
            ur_obs = ur_obs + instr
        if p.uss is True and err_uss is not None:
            ur_obs = ur_obs + err_uss
        if p.file_input is not None:
            ur_obs[numpy.where(ur_true == p.model_nan)] = p.model_nan
        return ur_obs