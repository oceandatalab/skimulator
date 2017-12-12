'''
   FIG. 3: Model interpolated currents and the corresponding radial currents.
'''

import netCDF4
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import skimsimulator.rw_data as rw
import glob
import cartopy
import os

# Initialize color
listcolor = ['c', 'y', 'b', 'g', 'k', 'r', 'c', 'y']

# List files
indatadir = '/mnt/data/project/'
indatadir = '/tmp/key/project/'
indatadir = os.path.join(indatadir, 'skim', 'skim_output')
config="WW3_GS_6b108az"
config="WW3_GS_8b105az"
modeldatapath = '/tmp/key/model/ww3_gs/ww3.201109_cur.nc'
filesgrid = os.path.join(indatadir, '{}_'.format(config))
ipass = 59
indatapath = '{}c01_p{:03d}.nc'.format(filesgrid, ipass)
listfile = glob.glob(indatapath)
listfile = sorted(listfile)
outdatadir = '../'
modelbox = [-72, -70., 34., 35.5]
scale = 8
#modelbox = [270, 290., 32., 40.]
is_cartopy = False
# Prepare figure
pyplot.figure(figsize=(10, 15))
projection = cartopy.crs.Mercator()
transform = cartopy.crs.PlateCarree()
if is_cartopy is True:
    projection = cartopy.crs.Mercator()
    transform = cartopy.crs.PlateCarree()
    #projection = transform
    ax1 = pyplot.subplot(121, projection=projection)
else:
    ax1 = pyplot.subplot(121)
if is_cartopy is True:
    ax1.add_feature(cartopy.feature.OCEAN, zorder=1)
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
    pyplot.pcolor(mlon, mlat, numpy.sqrt(mu**2 + mv**2), cmap='jet',
                  transform=transform)
    pyplot.quiver(mlon[::s0, ::s0], mlat[::s0, ::s0], mu[::s0, ::s0],
                  mv[::s0, ::s0], transform=transform)
else:
    pyplot.pcolor(mlon, mlat, numpy.sqrt(mu**2 + mv**2), cmap='jet')
    pyplot.quiver(mlon[::s0, ::s0], mlat[::s0, ::s0], mu[::s0, ::s0],
                  mv[::s0, ::s0])



indatapath = '{}grid_p{:03d}.nc'.format(filesgrid, ipass)
datag = netCDF4.Dataset(indatapath, 'r')
lon = data['lon'][:]
lon[lon > 180] = lon[lon > 180] - 360
lat = data['lat'][:]
ax2 = pyplot.subplot(122,  projection=projection)
ax2.add_feature(cartopy.feature.OCEAN, zorder=1)
ax2.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='black')
print(indatapath)
#ax2.set_extent(modelbox, projection)
gl = ax2.gridlines(crs=transform, draw_labels=True, color='gray',
             linestyle='--', alpha=0.5)
gl.xlabels_top = False
gl.ylabels_left = False
ax2.set_title('(b) Interpolated velocity (in green) and projected velocity '
              'along the radial angle (in red)')
lon = data['lon'][:]
lon[lon > 180] = lon[lon > 180] - 360
lat = data['lat'][:]
ur = data['ur_model'][:]
corrangle = datag['radial_angle'][:]
datau = data['u_model'][:]
datav = data['v_model'][:]
uur = ur * numpy.cos(corrangle)
uur[uur.mask] = 0
vur = ur * numpy.sin(corrangle)
vur[vur.mask] = 0
for i in range(numpy.shape(lon)[1]):
    style_color = '{}+'.format(listcolor[i])
    pyplot.plot(lon[:, i], lat[:, i], style_color,
                transform=transform)
    pyplot.quiver(lon[:, i], lat[:, i], datau[:, i], datav[:, i],
                  color='green', scale=scale, transform=transform)
    pyplot.quiver(lon[:, i], lat[:, i], uur[:, i], vur[:, i], color='red',
                  scale=scale, transform=transform)
pyplot.savefig(os.path.join(outdatadir, 'Fig3.png'))
