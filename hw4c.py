import numpy             as np
import matplotlib.pyplot as plt
import netCDF4           as nc
import sys

#save input limit
n = int(sys.argv[1])
g = 9.8

#loading data from .nc file.
rootgrp = nc.Dataset('/home/B12/b12209017/hw4/hw4/ERA5_Easia.nc')
IVT  = np.zeros((31, 201, 321))
q_d  = rootgrp.variables['q']
t_d  = rootgrp.variables['t']
u_d  = rootgrp.variables['u']
v_d  = rootgrp.variables['v']
time = rootgrp.variables['time'][:]
lat  = rootgrp.variables['lat' ][:]
lon  = rootgrp.variables['lon' ][:]
lev  = rootgrp.variables['lev' ][:]*100

#save pressure data from 1d array to 3d array,
#where dlev(i,:,:) save each layer of pressure.
dlev = np.zeros((len(lev), len(lat), len(lon)))
for i in range(len(lev)):
    dlev[i, :, :][:] = np.full((201, 321), lev[i])

#calculate each day's IVT value, save in a 3d array with shape(31, 201, 321)
for j in range(len(time)):
    IVTu = (((q_d[j,:-1]*u_d[j,:-1]+q_d[j,1:]*u_d[j,1:])*(dlev[:-1]-dlev[1:]))*(-1/g)/2).sum(0)
    IVTv = (((q_d[j,:-1]*v_d[j,:-1]+q_d[j,1:]*v_d[j,1:])*(dlev[:-1]-dlev[1:]))*(-1/g)/2).sum(0)
    IVT[j] = np.sqrt(IVTu**2 + IVTv**2)

#if IVT>n modify it to 1, and others to 0
boolIVT    = np.where(IVT > n, IVT, 0)
boolIVT    = np.where(boolIVT == 0., boolIVT, 1)
#sum of 31 days' bool value, division by 31.
BoolIVT    = (boolIVT[:,:,:].sum(0))/31

#loading data of East Asia map
mlon, mlat = np.loadtxt('/home/B12/b12209017/hw4/hw4/Easia_coastline.txt', dtype=float, comments=None, delimiter=',', skiprows=1, unpack=True)
mlon       = np.where(mlon > -123456, mlon, np.nan)
mlat       = np.where(mlat > -123456, mlat, np.nan)

#figure set
plt.subplots(1, 1, figsize=(8,6))
plt.title('Heat Map of IVT (> %.1f kg/m/s)' %(n))
plt.xlim([80, 160])
plt.ylim([-10, 40])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.plot(mlon, mlat, 'k')
CS=plt.pcolormesh(lon, lat, BoolIVT, vmax=1.0, vmin=0, cmap = plt.cm.hot_r)
plt.colorbar(CS,orientation='vertical')
plt.savefig('heatmap_%05d.png' %(n), bbox_inches='tight', dpi=300)
