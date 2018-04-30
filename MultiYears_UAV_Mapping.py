import netCDF4 as nc
import numpy as np
import datetime
from netCDF4 import Dataset
T1 = datetime.datetime.now()

# Set Parameter Thresholds
P = 1
W = 25
T = [-30, 50]
# Create name list
name = [2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016]


for n in name:
    pfile = nc.Dataset(r"Y:\Mo_GlobalUAVMapping\High_res_data\Globe_%s_3hr_0.25_precip.nc" % n, 'r')
    wtfile = nc.Dataset(r"Y:\Mo_GlobalUAVMapping\High_res_data\Globe_%s_3hr_0.25_windtemp.nc" % n, 'r')
    # extract latitude, longitude, and time
    time = pfile.variables['time'][:]
    lat = pfile.variables['latitude'][:]
    lon = pfile.variables['longitude'][:]
    # Check the Dimension
    t = len(time)
    la = len(lat)
    lo = len(lon)
    # Leap Year checking
    if t == 2928:
        y = 366
        Daylightfile = r"Y:\Mo_GlobalUAVMapping\High_res_data\daylightleap.nc"
    else:
        y = 365
        Daylightfile = r"Y:\Mo_GlobalUAVMapping\High_res_data\daylight.nc"

    Dayfile = nc.Dataset(Daylightfile)
    print "All Files loaded"

    index = np.arange(t)
    data = np.zeros(t * la * lo).reshape(t, la, lo)
    for i in index:
        u = wtfile.variables['u10'][i, :, :]
        v = wtfile.variables['v10'][i, :, :]
        temp = wtfile.variables['t2m'][i, :, :]
        prep = pfile.variables['tp'][i, :, :]
        day = Dayfile.variables['daytime'][i, :, :]
        # PreProcessing
        wind = (u ** 2 + v ** 2) ** 0.5
        temp = temp - 273.15
        # Array Filter
        # Wind
        wa = wind <= W
        wb = wind > W
        wind[wa] = 1
        wind[wb] = 0
        print "wind filter finished"
        # Precipitation
        pa = prep < P
        pb = prep >= P
        prep[pa] = 1
        prep[pb] = 0
        print "precip filter finished"
        # Temperature
        ta = (temp <= T[1]) & (temp >= T[0])
        tb = (temp > T[1]) | (temp < T[0])
        temp[ta] = 1
        temp[tb] = 0
        print "temperature filter finished"
        # DroneDay
        Overall = wind + temp + prep + day
        oa = Overall == 4
        ob = Overall != 4
        Overall[oa] = 1
        Overall[ob] = 0
        data[i, :, :] = Overall
        print i

    # Covert Drone 3hour to Dronday

    ngrid = np.zeros(y * la * lo).reshape(y, la, lo)
    index = np.linspace(8, t, y)
    index = index.astype(int)

    k = 0
    for ti in index:

        ui = ti - 1
        li = ti - 8
        cell = data[li:ui,:,:]
        s = np.sum(cell, axis=0)
        s[s > 0] = 1
        ngrid[k, :, :] = s
        print k
        k = k + 1

    print "Converted to daily Scale!!!"

    dronedata = Dataset(r"Y:\Mo_GlobalUAVMapping\High_res_data\SkyRanger_1\Droneday%s" % n, 'w', format='NETCDF4_CLASSIC')
    time = dronedata.createDimension('time', y)
    lat = dronedata.createDimension('lat', la)
    lon = dronedata.createDimension('lon', lo)
    dt = dronedata.createVariable('daytime', np.int32, ('time', 'lat', 'lon'))
    dt[:] = ngrid
    print "droneday files created"
    # Close File
    dronedata.close()
    pfile.close()
    wtfile.close()
    Dayfile.close()
    print "Year is " + str(n)

print "All File Created"
T2 = datetime.datetime.now()
print T2 - T1