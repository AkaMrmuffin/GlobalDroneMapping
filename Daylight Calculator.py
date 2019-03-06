# Name:     Daylight Calculator
# Purpose:  The purpose of this tool is to Calculate the determine wether the inputted location and time is daytime or not 
# Input:    ERA-interm grid    
# Output:   Daylight matrix based on ERA-interm grid  
# Author:   Mozhou Gao
# Project:  Global Dronedays Mapping 
# Created:  05/03/2019
# Copyright:(c) mozhou.gao 2019


# import the side packages 
# datetime ephem: calculates the solar angle 
import ephem
import datetime
# math, numpy, netCDF4 were used to manipulate ERA-interm grid
import math
import numpy as np
import netCDF4 as nc
# os was used to manipulate the operation system
from os import listdir
from os.path import isfile, join


def doyandhod (t,t_base):
    # All time are based on gregorian calender
    # t_base: The first time step of the gregorian calender time
    # t: A list or array that include all time step
    # calculate the timestep difference
    dif = t/24 - t_base/24
    # calculate the day of the year
    doy = int(dif)+1
    # calculate the hour of the day
    hod = t%24
    return doy,hod

def calsunlight(lat, lon, tobj):
    lat = str(lat)
    lon = str(lon)
    observation = tobj
    # create the location of the sun
    sun = ephem.Sun()
    # create the observer object that relates to the current lat, long and time
    observer = ephem.Observer()
    observer.lat, observer.lon, observer.elevation = lat, lon, 0
    # UTC time
    # date time object is a string
    observer.date = observation
    # calculate the solar angle 
    sun.compute(observer)
    current_sun_alt = sun.alt
    # convert the solar altitude to the degree
    return current_sun_alt*180/math.pi
	

pfile = nc.Dataset(r"ERA-interm grid file", 'r')

time = pfile.variables['time'][:]
lat = pfile.variables['latitude'][:]
lon = pfile.variables['longitude'][:]
t = len(time)
la = len(lat)
lo = len(lon)
print (t,la,lo)

daylight = np.empty((len(time),len(lat),len(lon)),dtype = np.float)
 
i= 0 
year = 2008
while i<t:
    dy, hd = doyandhod(time[i], time[0])
    # calculate the date time object
    dtobj = datetime.datetime(year, 1, 1) + datetime.timedelta(dy - 1)
    # extract date object
    dobj = dtobj.date()
    # get the time object
    tobj = datetime.time(hd)
    # combine to the new date time object
    nobj = datetime.datetime.combine(dobj, tobj)
    j = 0
    while j<la:
        latitude = lat[j]
        k = 0
        while k<lo:
            longitude = lon[k]
            ang = calsunlight(latitude, longitude, nobj)
            # we use the urban twilight as the boundary (-6 degree), if greater than -6 degree then it is a day time, 
            # otherwise it is a night time
            if ang > -6:
                daylight[i,j,k]=1
            else:
                daylight[i,j,k]=0
            k += 1
        j += 1
    print (i)
    i += 1
# Export data 
dronedata =nc.Dataset(r"D:\leapyeardaylightmatrix.nc", 'w', format='NETCDF4_CLASSIC')
# Create dimension variable
lat = dronedata.createDimension('lat', la)
lon = dronedata.createDimension('lon', lo)
time = dronedata.createDimension('time', t)
# Create variable 
dt = dronedata.createVariable('daytime', np.int32, ('time','lat', 'lon'))
# R is our result
dt[:] = daylight

# Close File
dronedata.close()
pfile.close()