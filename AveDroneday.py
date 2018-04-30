import netCDF4 as nc
import numpy as np
import datetime

T1 = datetime.datetime.now()
newd = np.zeros(721*1440).reshape(721,1440)
name = [2007,2008,2009,2010, 2011, 2012, 2013, 2014, 2015, 2016]
for n in name:
    f = nc.Dataset(r"C:\Users\mozhou.gao\Desktop\TyphoneH_1\Droneday%s" % n, 'r')
    i = 0
    while i<721:

        d = f.variables['daytime'][:,i,:]
        s = np.sum(d,axis =0)

        newd[i,:] += s
        i += 1
        print i

    f.close()
    print "Year: " + str(n)

newd = newd/10
print np.max(newd)

dronedata = nc.Dataset("Ave_TyphoneH_Drone_Day_p1.nc", 'w', format='NETCDF4_CLASSIC')
time = dronedata.createDimension('time', 1)
lat = dronedata.createDimension('lat', 721)
lon = dronedata.createDimension('lon', 1440)
dt = dronedata.createVariable('droneday', np.int32, ('time', 'lat', 'lon'))
dt[:] = newd

print "Average droneday file created!"
dronedata.close()
T2 = datetime.datetime.now()
print T2 - T1
