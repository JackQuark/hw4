import numpy             as np
import matplotlib.pyplot as plt
import netCDF4           as nc
import sys

#tin        = str(sys.argv[1])
tin        = '2016.07.05'
tins       = tin.split('.')
if(int(tins[1]) == 7):
    t = int(tins[2]) + 15
else:
    t = int(tins[2]) - 15

IVT = IVTu = IVTv = IVTU = IVTV = np.zeros((201, 321))
mlon, mlat = np.loadtxt('Easia_coastline.txt', dtype=float, comments=None, delimiter=',', skiprows=1, unpack=True)
mlon       = np.where(mlon > -10000, mlon, np.nan)
mlat       = np.where(mlat > -10000, mlat, np.nan)
g          = 9.8

rootgrp = nc.Dataset('ERA5_Easia.nc')

q_d  = rootgrp.variables['q']
u_d  = rootgrp.variables['u']
v_d  = rootgrp.variables['v']
time = rootgrp.variables['time'][:]
lat  = rootgrp.variables['lat' ][:]
lon  = rootgrp.variables['lon' ][:]
lev  = rootgrp.variables['lev' ][:]*100

dlev = np.zeros((len(lev), len(lat), len(lon)))
for i in range(len(lev)):
    dlev[i, :, :][:] = np.full((201, 321), lev[i])

IVTu = (((q_d[t,:-1]*u_d[t,:-1]+q_d[t,1:]*u_d[t,1:])*(dlev[:-1]-dlev[1:]))*(-1/g)/2)
IVTv = (((q_d[t,:-1]*v_d[t,:-1]+q_d[t,1:]*v_d[t,1:])*(dlev[:-1]-dlev[1:]))*(-1/g)/2)

#sum of IVTu/IVTv in 31 layers of pressure, saving in a 2d array.
IVTU = IVTu[:,:,:].sum(0)
IVTV = IVTv[:,:,:].sum(0)

IVT  = np.sqrt(IVTU**2 + IVTV**2)

plt.subplots(1, 1, figsize=(8,6))
plt.title('IVT [kg/m/s] 2016.%02d.%02d' %(int(tins[1]), int(tins[2])))
plt.xlim([80, 160])
plt.ylim([-10, 40])
plt.plot(mlon, mlat, 'k')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
r = np.linspace(100, 1000, 10)
CS=plt.contourf(lon, lat, IVT[:,:], levels = r, cmap = plt.cm.Blues, extend = 'both')
CS.set_cticks = r
plt.colorbar(CS,orientation='vertical')
plt.savefig('IVT_2016.%02d.%02d.png' %(int(tins[1]), int(tins[2])), bbox_inches='tight', dpi=300)