import numpy             as np
import matplotlib.pyplot as plt
import netCDF4           as nc
import sys

#save input date
tin        = str(sys.argv[1])
tins       = tin.split('.')
if(int(tins[1]) == 7):
    t = int(tins[2]) + 15
else:
    t = int(tins[2]) - 15

#loading data from .nc file.
rootgrp = nc.Dataset('/home/B12/b12209017/hw4/hw4/ERA5_Easia.nc')
q_d  = rootgrp.variables['q']
u_d  = rootgrp.variables['u']
v_d  = rootgrp.variables['v']
time = rootgrp.variables['time'][:]
lat  = rootgrp.variables['lat' ][:]
lon  = rootgrp.variables['lon' ][:]
lev  = rootgrp.variables['lev' ][:]*100
g    = 9.8

#save pressure data from 1d array to 3d array,
#where dlev(i,:,:) save each layer of pressure.
dlev = np.zeros((len(lev), len(lat), len(lon)))
for i in range(len(lev)):
    dlev[i, :, :][:] = np.full((201, 321), lev[i])

#calculate sum of IVTu/IVTv in 31 layers of pressure,
#and use them to calculate IVT(2d array(len(lat), len(lon))).
IVTu = (((q_d[t,:-1]*u_d[t,:-1]+q_d[t,1:]*u_d[t,1:])*(dlev[:-1]-dlev[1:]))*(-1/g)/2).sum(0)
IVTv = (((q_d[t,:-1]*v_d[t,:-1]+q_d[t,1:]*v_d[t,1:])*(dlev[:-1]-dlev[1:]))*(-1/g)/2).sum(0)
IVT  = np.sqrt(IVTu**2 + IVTv**2)

#loading data of East Asia map
mlon, mlat = np.loadtxt('/home/B12/b12209017/hw4/hw4/Easia_coastline.txt', dtype=float, comments=None, delimiter=',', skiprows=1, unpack=True)
mlon       = np.where(mlon > -123456, mlon, np.nan)
mlat       = np.where(mlat > -123456, mlat, np.nan)

#figure set
plt.subplots(1, 1, figsize=(8,6))
plt.title('IVT [kg/m/s] 2016.%02d.%02d' %(int(tins[1]), int(tins[2])))
plt.xlim([80, 160])
plt.ylim([-10, 40])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.plot(mlon, mlat, 'k')
r = np.linspace(100, 1000, 10)
CS=plt.contourf(lon, lat, IVT[:,:], levels = r, cmap = plt.cm.Blues, extend = 'both')
CS.set_cticks = r
plt.colorbar(CS,orientation='vertical')
plt.savefig('IVT_2016.%02d.%02d.png' %(int(tins[1]), int(tins[2])), bbox_inches='tight', dpi=300)
