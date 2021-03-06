'''
   FIG. 3: Model interpolated currents and the corresponding radial currents.
'''

import netCDF4
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import skimulator.rw_data as rw
import glob
import cartopy
import os

# Initialize color
listcolor = ['c', 'y', 'b', 'g', 'k', 'r', 'c', 'y']

# List files
config="WW3_AT_metop_2018_8b"
indatadir = '/mnt/data/project/'
indatadir = '/tmp/key/data/skim_at_output/{}'.format(config)

modeldatapath = '/tmp/key/data/model/ww3_gs/ww3.201109_cur.nc'
filesgrid = os.path.join(indatadir, '{}_'.format(config))
ipass = 58
indatapath = '{}c01_p{:03d}.nc'.format(filesgrid, ipass)
listfile = glob.glob(indatapath)
listfile = sorted(listfile)
outdatadir = '../'
modelbox = [-71, -69.0, 34.0, 36.0]
scale = 11
is_cartopy = True
# Prepare figure
pyplot.figure(figsize=(10, 15))
projection = cartopy.crs.PlateCarree()
transform = cartopy.crs.PlateCarree()
if is_cartopy is True:
    #projection = transform
    ax1 = pyplot.subplot(121, projection=projection)
else:
    ax1 = pyplot.subplot(121)
if is_cartopy is True:
    ax1.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
    ax1.set_extent(modelbox, projection)
    gl = ax1.gridlines(crs=transform, draw_labels=True, color='gray',
                       linestyle='--', alpha=0.5)
    gl.xlabels_top = False
    gl.ylabels_left = False
ax1.set_title('(a) Model, "true" velocity')

data = netCDF4.Dataset(indatapath, 'r')
datam = netCDF4.Dataset(modeldatapath, 'r')
iifile = data['index'][:, 0]
lon = datam['longitude'][:]
lon[lon > 180] = lon[lon > 180] - 360
lat = datam['latitude'][:]
mlon, mlat = numpy.meshgrid(lon, lat)
mu = datam['ucur'][int(iifile[0]),:, :]
mv = datam['vcur'][int(iifile[0]),:, :]
mu=numpy.ma.masked_invalid(mu)
numpy.ma.masked_where(mu == datam['ucur']._FillValue, mu, copy=False)
mu[mu.mask] = 0
numpy.ma.masked_where(mv == datam['vcur']._FillValue, mv, copy=False)
mv[mv.mask] = 0
s0 = 10
if is_cartopy is True:
    pyplot.pcolor(mlon, mlat, numpy.sqrt(mu**2 + mv**2), vmax=2.0, cmap='jet',
                  transform=transform)
    pyplot.quiver(mlon[::s0, ::s0], mlat[::s0, ::s0], mu[::s0, ::s0],
                  mv[::s0, ::s0], scale=scale, transform=transform)
else:
    pyplot.pcolor(mlon, mlat, numpy.sqrt(mu**2 + mv**2), vmax=2.0, cmap='jet')
    pyplot.quiver(mlon[::s0, ::s0], mlat[::s0, ::s0], mu[::s0, ::s0],
                  mv[::s0, ::s0])



indatapath = '{}grid_p{:03d}.nc'.format(filesgrid, ipass)
datag = netCDF4.Dataset(indatapath, 'r')
lon = data['lon'][:]
lon[lon > 180] = lon[lon > 180] - 360
lat = data['lat'][:]
ax2 = pyplot.subplot(122,  projection=projection)
ax2.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
print(indatapath)
ax2.set_extent(modelbox, projection)
gl = ax2.gridlines(crs=transform, draw_labels=True, color='gray',
             linestyle='--', alpha=0.5)
gl.xlabels_top = False
gl.ylabels_left = False
ax2.set_title('(b) Interpolated (green) and projected (red) velocity')
lon = data['lon'][:]
lon[lon > 180] = lon[lon > 180] - 360
lat = data['lat'][:]
ur = data['ur_true'][:]
corrangle = datag['radial_angle'][:]
datau = data['u_model'][:]
datav = data['v_model'][:]
uur = ur * numpy.cos(corrangle)
uur[uur.mask] = 0
vur = ur * numpy.sin(corrangle)
vur[vur.mask] = 0
for i in range(numpy.shape(lon)[1]):
    style_color = '{}+'.format(listcolor[i])
    # pyplot.plot(lon[:, i], lat[:, i], style_color,
    #             transform=transform)
    pyplot.quiver(lon[:, i], lat[:, i], datau[:, i], datav[:, i],
                  color='green', scale=scale, transform=transform)
    pyplot.quiver(lon[:, i], lat[:, i], uur[:, i], vur[:, i], color='red',
                  scale=scale, transform=transform)
pyplot.savefig(os.path.join(outdatadir, 'Fig3.png'))
