import os
import numpy
import healpy as heal
from numpy import linalg as LA
from netCDF4 import Dataset
import params as p
import skimulator.const as const
import skimulator.rw_data as rw

import matplotlib.pyplot as plt

theta1 = const.theta1
theta0 = const.theta0
gamma0 = const.gamma0

# - In parameter file ## TODO -
# Number of pixel (resolution for healpix)
nside = 256
# Number of diamonds for healpix
ndiam = 12
ntotpixel = nside * nside * ndiam
# Conditionning threshold
thresh_cond = 10
# List of pass
swathtab = [40] #,'p025','p042','p053','p055']
# 
cycle = 1
plot_res = 10


gamma = numpy.array(p.list_angle)
ind_12 = numpy.where(gamma == 12)
ind_6 = numpy.where(gamma == 6)
gamma = numpy.deg2rad(gamma)
n_ind_12 = len(ind_12)
n_ind_6 = len(ind_6)

path = os.path.join(p.outdatadir, '{}_c{:02d}'.format(p.config, cycle))

gpath = os.path.join(p.outdatadir, '{}_grid'.format(p.config))


#heal.mollview(np.arange(12*64*64))
#plt.show()

# rebin variance to beam
sig1 = 1E-2*numpy.loadtxt(p.rms_instr[ind_12[0][0]], usecols=(0, 1), unpack=True)
sig2 = 1E-2*numpy.loadtxt(p.rms_instr[ind_6[0][0]], usecols=(0, 1), unpack=True)
sig = numpy.zeros([len(gamma), 360])

for i in ind_12[0]:
    sig[i, : 180] = sig2[1, : 180]
    sig[i, 180 + numpy.arange(180)] = sig2[1, 180 - numpy.arange(180)]
for i in ind_6[0]:
    sig[i, 0: 180] = sig1[1, 0: 180]
    sig[i, 180 + numpy.arange(180)] = sig1[1, 180 - numpy.arange(180)]

#===== table to store image inversion

im = numpy.zeros((ntotpixel, 3))
nim = numpy.zeros((ntotpixel))
cov = numpy.zeros((ntotpixel, 2, 2))
cov2 = numpy.zeros([ntotpixel, 2, 2])
vec = numpy.zeros([ntotpixel, 2])
vec2 = numpy.zeros([ntotpixel, 2])
vec3 = numpy.zeros([ntotpixel, 2])
vecdop = numpy.zeros([3, ntotpixel, 2])

for swath in swathtab:
    data = rw.Sat_SKIM(ifile='{}_p{:03d}.nc'.format(path, swath))
    grid = rw.Sat_SKIM(ifile='{}_p{:03d}.nc'.format(gpath, swath))

    data.load_data(p, ur_model=[], ur_obs=[], instr=[], u_model=[], v_model=[],
                    lon_nadir=[], lat_nadir=[], time_nadir=[], lon=[], lat=[])
    grid.load_swath(p, radial_angle=[])
    ur = data.ur_model
    uro = data.ur_obs
    noise = data.instr
    u = data.u_model
    v = data.v_model

    rangle = numpy.mod(grid.radial_angle, numpy.pi *2)

    lon = data.lon
    lat = data.lat

    ndata, nbeam = ur.shape

    # COMPUTE SIGMA AT EACH SAMPLE

    sigma = numpy.zeros(ur.shape)
    arctan_rangle = numpy.arctan(numpy.sin(rangle + numpy.pi/2),
                                 numpy.cos(rangle + numpy.pi/2))
    _tmp = numpy.rad2deg(numpy.mod(arctan_rangle + 2 * numpy.pi,
                                       2 * numpy.pi))
    for i in range(len(gamma)):
        sigma[:, i] = sig[i, (_tmp[:, i]).astype(int)]


    #COMPUTE ORBITAL SPEED FROM NADIR INFORMATION
    # TODO MANAGE BUG OF -180 to +180 TODOTODOTODO
    lat_nadir = data.lat_nadir
    lon_nadir = data.lon_nadir
    time_nadir = data.time_nadir

    time_nadir = time_nadir - time_nadir[0]

    latpar = numpy.polyfit(time_nadir, lat_nadir, 3)
    lonpar = numpy.polyfit(time_nadir, lon_nadir, 3)

    dlatorb = latpar[2] + 2*time_nadir*latpar[1] + 3*time_nadir**2*latpar[0]
    dlonorb = lonpar[2] + 2*time_nadir*lonpar[1] + 3*time_nadir**2*lonpar[0]

    #altiorb=793.0
    alti_total = const.Rearth + const.sat_elev
    usat = (numpy.deg2rad(dlonorb) * numpy.cos(numpy.deg2rad(lat_nadir))
            * alti_total)
    vsat = numpy.deg2rad(dlatorb) * alti_total

    dtheta = theta1 - gamma0 * numpy.sin(rangle - theta0)
    dgamma = gamma0 * numpy.cos(rangle - theta0)

    # nongeophysical doppler to add to measure
    gamma_mat = numpy.array([gamma,] * ndata)
    ddop = numpy.zeros(ur.shape)
    usat = numpy.array([usat, ]*len(gamma)).transpose()
    vsat = numpy.array([vsat, ]*len(gamma)).transpose()
    ddop1 = dtheta * (usat * numpy.sin(rangle) + vsat * numpy.cos(rangle))
    ddop2 = dgamma * numpy.cos(gamma_mat) * (usat * numpy.cos(rangle)
                                             + vsat * numpy.sin(rangle))
    ddop = ddop1 + ddop2

    tdop=numpy.zeros((ndata, nbeam, 3))
    tdop[:, :, 0] = usat * numpy.sin(rangle) + vsat * numpy.cos(rangle)
    tdop[:, :, 1] = (-numpy.sin(rangle) * (usat*numpy.sin(rangle)
                     + vsat * numpy.cos(rangle)) + numpy.cos(rangle)
                     * numpy.cos(gamma_mat) * (usat*numpy.cos(rangle)
                     + vsat * numpy.sin(rangle)))
    tdop[:, :, 2] = (numpy.cos(rangle) * (usat * numpy.sin(rangle) + vsat
                     * numpy.cos(rangle)) + numpy.sin(rangle)
                     * numpy.cos(gamma_mat) * (usat * numpy.cos(rangle)
                     + vsat * numpy.sin(rangle)))

    # DEBUG PLOT
    #plt.plot(ddop-(tdop[:,:,0]*theta1+tdop[:,:,1]*gamma0*np.cos(theta0)+tdop[:,:,2]*gamma0*np.sin(theta0)))
    #plt.show()

    uro = ddop + uro

    ww = 1 / sigma**2

    ph = 2 * numpy.pi - numpy.deg2rad(lon)
    th = numpy.pi / 2 - numpy.deg2rad(lat)
    pidx = heal.ang2pix(nside, th, ph)

    for i in range(nbeam):
        for j in range(ndata):
            if ur[j, i] > -1E9:
                ip = pidx[j,i]
                # compute imulated model
                im[ip, 1] += u[j, i]
                im[ip, 2] += v[j,i]
                nim[ip] += 1
                # compute covariance(s) model
                co = numpy.cos(rangle[j,i])
                si = numpy.sin(rangle[j,i])
                w = ww[j,i]
                cov[ip, 0, 0] += co * co
                cov[ip, 1, 0] += si * co
                cov[ip, 0, 1] += si * co
                cov[ip, 1, 1] += si * si

                cov2[ip, 0, 0] += w * co * co
                cov2[ip, 1, 0] += w * si * co
                cov2[ip, 0, 1] += w * si * co
                cov2[ip, 1, 1] += w * si * si

                # compute data vector model
                vec[ip, 0] += co * ur[j,i]
                vec[ip, 1] += si * ur[j,i]

                # compute data noise vector model
                vec2[ip, 0] += w* co * uro[j,i]
                vec2[ip, 1] += w * si * uro[j,i]

                # compute doppler projection
                for k in range(3):
                    vecdop[k, ip, 0] += w * co * tdop[j,i,k]
                    vecdop[k, ip, 1] += w * si * tdop[j,i,k]

rim = 0 * im
rim2 = 0 * im
mask = numpy.zeros([ntotpixel])
rimdop = numpy.zeros([3, ntotpixel, 2])

for i in range(ntotpixel):
    if cov2[i, 0, 0] > 0:
        mask[i] = LA.cond(cov2[i, :, :])
        if mask[i] < thresh_cond:
            rim[i, 1: 3] = LA.solve(cov[i, :, :],vec[i, :])
            rim2[i, 1: 3] = LA.solve(cov2[i, :, :],vec2[i, :])
            for k in range(3):
                rimdop[k, i, :] = LA.solve(cov2[i, :, :], vecdop[k, i, :])

            im[i, 0] = numpy.sqrt((im[i, 1]/nim[i])**2 + (im[i, 2]/nim[i])**2)
            rim[i, 0] = numpy.sqrt(rim[i, 1]**2 + rim[i, 2]**2)
            rim2[i, 0] = numpy.sqrt(rim2[i, 1]**2 + rim2[i, 2]**2)

mat=numpy.zeros([3,3])
vres=numpy.zeros([3])
for i in range(nbeam):
    for j in range(ndata):
        ip = pidx[j,i]
        if mask[ip] < thresh_cond and ur[j, i] > -1E9:
            co = numpy.cos(rangle[j, i])
            si = numpy.sin(rangle[j, i])
            w = ww[j, i]

            vv = uro[j, i] - rim2[ip, 1] * co - rim2[ip, 2] * si
            tvv = numpy.zeros((3,))
            for k in range(3):
                tvv[k] = (tdop[j, i, k] - rimdop[k, ip, 0] * co
                          - rimdop[k, ip, 1] * si)

            for k in range(3):
                for l in range(3):
                    mat[k, l] += w * tvv[k] * tvv[l]

            for k in range(3):
                vres[k] += w * tvv[k] * vv

rr = LA.solve(mat, vres)
print(rr)
#print(rr[0], numpy.arctan(rr[1], rr[2]), numpy.sqrt(rr[1]**2 + rr[2]**2))
rot=[numpy.rad2deg(numpy.mean(ph)), numpy.rad2deg(numpy.pi/2 - numpy.mean(th))]

heal.gnomview(mask,rot=rot,reso=plot_res * 3,title='CONDITIONING',max=thresh_cond)
#plt.pcolormesh()
plt.savefig('{}_conditioning.png'.format(p.config))
heal.gnomview(im[:,0],rot=rot, reso=plot_res, title='MODEL',xsize=600)
plt.savefig('{}_model.png'.format(p.config))
heal.gnomview(rim[:,0],rot=rot,reso=plot_res,title='INVERSE',xsize=600)
plt.savefig('{}_inverse.png'.format(p.config))
heal.gnomview(rim2[:,0],rot=rot,reso=plot_res,title='INVERSE WITH NOISE+DOPPLER',xsize=600)
plt.savefig('{}_inverse_noise.png'.format(p.config))
