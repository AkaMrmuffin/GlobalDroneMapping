# Name:     Dronedays Calculator
# Purpose:  The purpose of this tool is to Calculate the average dronedays by using ERA-Intntrim Grid weather data  
# Input:    ERA-interm grid, weather: temperature, wind speed, precipitation, Drone's specifications     
# Output:   Average dronedays 
# Author:   Mozhou Gao
# Project:  Global Dronedays Mapping
# Created:  05/03/2019
# Copyright:(c) mozhou.gao 2019

# Import side packages 
from numba import double 
from numba.decorators import jit, autojit
import ephem
import datetime
import math
import numpy as np
import netCDF4 as nc
from os import listdir
from os.path import isfile, join
import pandas as pd

# load Drone spec table
df = pd.read_excel(r'load your drone SPECs table')

print (df.head())
# create drone types index 
drone_index = np.arange(0, 12, 1)

# Read weather data
mypath1 = r'your wind and temperature data files'
Windfiles = [f for f in listdir(mypath1) if isfile(join(mypath1, f))]
mypath2 = r'your precipitation data files'
Precipfiles = [f for f in listdir(mypath2) if isfile(join(mypath2, f))]
name = [2007,2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016]
files = zip(name,Windfiles,Precipfiles)
files = list(files)
print(files)

def UAV_Mapping(files,W,P,T,D_type):
    '''
	UAV_Mapping can be used to batch calculate the dronedays 
	'''
    for n in files:
        print ("year: " + str (n[0]))
		# record the function start time 
        T1 = datetime.datetime.now()
        pfile = nc.Dataset(r"precipitation file folder"%n[2], 'r')
        wtfile = nc.Dataset(r"wind and temperature file folder"%n[1], 'r')
        # extract latitude, longitude, and time
        time = pfile.variables['time'][:]
        lat = pfile.variables['latitude'][:]
        lon = pfile.variables['longitude'][:]
        # Check the Dimension
        t = len(time)
        la = len(lat)
        lo = len(lon)
        # Check leap Year (if you dont want to add the daylight constraint, then you can ignore the loading of dayfile)
        if t == 2928:
            y = 366
            dayfile = nc.Dataset(r"leapyear file","r")
        else:
            y = 365
            dayfile = nc.Dataset(r"normal year","r")
       
        index = np.arange(t)
		# Create empty array for storing the averge drone days 
        data = np.empty((len(time),len(lat),len(lon)),dtype = np.float)
        for i in index:
			# assign all the variables 
            u = wtfile.variables['u10'][i, :, :]
            v = wtfile.variables['v10'][i, :, :]
            temp = wtfile.variables['t2m'][i, :, :]
            prep = pfile.variables['tp'][i, :, :]
			# if you dont want to add the daylight constraint, you can comment loading of the daytime file 
            day = dayfile.variables['daytime'][i,:,:]
            # PreProcessing
			# calculate the resultant wind speed 
            wind = (u ** 2 + v ** 2) ** 0.5
			# convert unit of temperature to celsius degree
            temp = temp - 273.15
            # Array Filter
            # Wind
            wa = wind <= W
            wb = wind > W
            wind[wa] = 1
            wind[wb] = 0

            # Precipitation
            pa = prep < P
            pb = prep >= P
            prep[pa] = 1
            prep[pb] = 0

            # Temperature
            ta = (temp <= T[1]) & (temp >= T[0])
            tb = (temp > T[1]) | (temp < T[0])
            temp[ta] = 1
            temp[tb] = 0
            
			
            # DroneDay (if you dont want to add daylight constraint, u can remove 'day')
            Overall = wind + temp + prep + day
            oa = Overall == 4
            ob = Overall != 4
            Overall[oa] = 1
            Overall[ob] = 0
            data[i, :, :] = Overall
            
            print (i)

        # Covert Drone 3hour to Dronday

        ngrid = np.zeros(shape=(y,la,lo))
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
            print (k)
            k = k + 1
         
        dday = np.sum(ngrid,axis= 0)

        dronedata =nc.Dataset(r"output folder" % D_type+str(n[0]), 'w', format='NETCDF4_CLASSIC')
        lat = dronedata.createDimension('lat', la)
        lon = dronedata.createDimension('lon', lo)
        dt = dronedata.createVariable('daytime', np.int32, ('lat', 'lon'))
        dt[:] = dday
    
        # Close File
        dronedata.close()
        pfile.close()
        wtfile.close()
        dayfile.close()
        T2 = datetime.datetime.now()
        
        print (T2-T1)

# Convert the normal function to multicore processing function 
uav_numba = autojit(UAV_Mapping)	
# run the numba function for each group of weather parameters 
for ind in drone_index:
    T = []
    T.append(df.iloc[ind, :]["TempL"])
    T.append(df.iloc[ind, :]["TempU"])
    W = df.iloc[ind, :]["Wind "]
    P = df.iloc[ind, :]["Precip"]
    typ = df.iloc[ind, :]["Drones"]
    uav_numba(files,W,P,T,typ)
		