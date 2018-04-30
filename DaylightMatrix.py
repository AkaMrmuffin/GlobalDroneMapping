# import side packages
import netCDF4 as nc
from netCDF4 import Dataset
import numpy as np

NCfile = r"D:\Mo_GlobalUAVMapping\High_res_data\Globe_2016_3hr_0.25_precip.nc"
# 0.25 x 0.25 degree & 3 hour .nc file path
grid = nc.Dataset(NCfile)
# extract latitude, longitude, and time
lat = grid.variables['latitude'][:]
lon = grid.variables['longitude'][:]
time = grid.variables['time'][:]
# assign dimensions
t = len(time)
la = len(lat)
lo = len(lon)
print t,la,lo
# assign time based hours/days index
timed_re = np.zeros(t)
timeh_re = np.zeros(t)
inh = 0
ind = 0
i = 0
while i < t:
    timed_re[i] = ind + 1
    timeh_re[i] = inh
    inh = inh + 3
    if inh == 24:
        ind = ind + 1
        inh = 0
    i = i + 1
print "time index created"
# calculate the solar declination for each time step
dec = 0.04093*np.sin(0.0172*(timed_re-82.2))

# Create time and lat based array for daylight length
D = np.zeros(t*la).reshape(t,la)
i = 0
while i < t:
    j = 0
    while j < la:
        # inside the square bracket
        num = (-1*np.sin(lat[j]*np.pi/180)*np.sin(dec[i])-0.1047)/(np.cos(lat[j]*np.pi/180)*np.cos(dec[i]))
        if num > -0.87:
            D[i,j] = 7.639*np.arccos(num)
        else:
            D[i,j] = 7.639*np.arccos(-0.87)
        j = j+1

    i = i+1
print "Daylength Calculation Done"
# calculate UTC zone based on longitude
lon = lon.astype(int)
UTC = []
for l in lon:
    if l < 0:
        s = -1
    else:
        s = 1
    re = np.abs(l)%15
    if re <=7.5:
        z = int(np.abs(l))/15
    elif re> 7.5:
        z = (int(np.abs(l))/15) + 1

    UTC.append(s*z)

UTC = np.array(UTC)
UTC = UTC.astype(int)
print "UTC zone Created!!"

# Create local time array
localtime = np.zeros(t*lo).reshape(t,lo)
k=0 # time index
while k < t:
    UTC00_time = timeh_re[k]
    j = 0 # longitude index
    while j<lo:

        cal_time_zone = UTC[j]

        cal_time = UTC00_time - (0- cal_time_zone)

        if 0<= cal_time<24:

            localtime[k,j] = cal_time
        elif cal_time < 0:
            localtime[k,j] = cal_time + 24
        else:
            localtime[k,j] = cal_time - 24
        j = j+1

    k= k+ 1
print "local time created!!!"

lt = localtime.astype(int)
dl = D.astype(int)


# Create a daylight 3d-matrix
daylight = np.zeros(t*lo*la).reshape(t,la,lo)
k =0
while k<t:
    i = 0
    while i<la:
        j = 0
        d = 8+dl[k,i]
        while j<lo:
            local = lt[k,j]
            if 8<local<d:
                daylight[k,i,j] = 1
            else:
                daylight[k,i,j] = 0
            j = j+1

        i = i+1

    k = k+1

print "daylight matrix created"
# For normal year
# daylight365= daylight[0:2920,:,:]

# Save the daylight file
data = Dataset("daylightleap.nc",'w',format = 'NETCDF4_CLASSIC')
time = data.createDimension('time',t)
lat = data.createDimension('lat',la)
lon = data.createDimension('lon',lo)
dt = data.createVariable('daytime',np.int32,('time','lat','lon'))
data.variables['daytime'][:] = daylight
data.close()

# close the nc file
grid.close()
print "close all files"