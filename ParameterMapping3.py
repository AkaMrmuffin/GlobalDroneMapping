import numpy as np
import netCDF4 as nc
from os import listdir
from os.path import isfile, join
from scipy import stats


import numpy as np
import netCDF4 as nc
from os import listdir
from os.path import isfile, join
from scipy import stats



mypath1 = r'C:\Users\mozhou\Desktop\UAV Article Figures\ParameterMap\ColorExo'
files = [f for f in listdir(mypath1) if isfile(join(mypath1, f))]
l = len(files)

temp = np.zeros(l*721*1440).reshape(l,721,1440)
index =0
for file in files:

    dfile = nc.Dataset(r'C:\Users\mozhou\Desktop\UAV Article Figures\ParameterMap\ColorExo\%s' % file,'r')

    d = dfile.variables['color'][:]
    temp[index,:,:] = d

    index = index +1
    dfile.close()


CM = np.zeros(721*1440).reshape(721,1440)
i = 0
while i<721:
    j = 0
    while j<1440:

        cell = temp[:,i,j]

        m = stats.mode(cell)

        CM[i,j] = m[0]

        j += 1

    i += 1

data =nc.Dataset(r"C:\Users\mozhou\Desktop\UAV Article Figures\ParameterMap\Color_0.nc", 'w', format='NETCDF4_CLASSIC')
lat = data.createDimension('lat', 721)
lon = data.createDimension('lon', 1440)
dt = data.createVariable('color', np.int32, ( 'lat', 'lon'))
dt[:] = CM
data.close()
