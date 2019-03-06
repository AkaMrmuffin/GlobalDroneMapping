# Name:     calculate station based dronedays (Validation script) 
# Purpose:  The purpose of this tool is to validate the ERA-interm dronedays calculated by Dronedays Calculator with the NOAA weather station data
# Input:    NOAA 5 mins weather station data    
# Output:   weather station dronedays  
# Author:   Mozhou Gao
# Project:  Global Dronedays Mapping
# Created:  05/03/2019
# Copyright:(c) mozhou.gao 2019

import numpy as np 
import pandas as pd 
from os import listdir
from os.path import isfile, join

#Simplelize The Data Table 
def simptable(df): 
    UTCday = df.iloc[:,1]
    UTCday.pop(0)
    DOY =[]
    for d in UTCday:
        DOY.append(d.timetuple().tm_yday)
    
    UCTtime = df.iloc[:,2]
    UCTtime.pop(0)
    Long = df.iloc[:,6]
    Long.pop(0)
    Lat = df.iloc[:,7]
    Lat.pop(0)
    Temp = df.iloc[:,8]
    Precip = df.iloc[:,9]
    Wind = df.iloc[:,21]
    WindFlag = df.iloc[:,22]
    Temp.pop(0)
    Precip.pop(0)
    Wind.pop(0)
    WindFlag.pop(0)
    
    da = {'UCT_Day':DOY , 'UCT_Time':UCTtime, 'Long': Long ,'Lat':Lat,
       "Temp": Temp,"Precip":Precip,"Wind": Wind,"WindFlag": WindFlag}
    
    ndf = pd.DataFrame(data=da,columns =['UCT_Day','UCT_Time','Long',
                                        'Lat','Temp','Precip','Wind','WindFlag'])
    
    return ndf
 
#Check the missing data 
def fixmissing(p):
    p = p.where(p != -99.000)
    p = p.where(p != -9999.0)
    p = p.where(p != 9999)
    p = p.interpolate()
    return p

### Calculate Flyability ### 
def weatherfilter(param):
    temp,wind,precip = param
    temp = float(temp)
    wind = float(wind)
    precip = float(precip)
    if (0<=temp<=40) and (wind<=10) and (precip<=0): 
        return 1 
    else: 
        return 0 
		
# merge the time steps 
def consthreeh (lis):
    longest = 0
    current = 0
    for num in lis:
        if num == 1:
            current += 1
        else:
            longest = max(longest, current)
            current = 0

    return max(longest, current)

# convert string to float 
def filtstr(param):
    col = param 
    if type(col) == str:        
        return np.float(col[:-2])
    else: 
        return col 

		
		
# load the weather station data folder 		
mypath = r'weather station data folder'
files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
# create empty list for storing the corresponding info from dataframe 
droneday =[]
daystore =[]
city = [] 
Lat =[] 
Long = [] 
index = 0
for f in files: 
    # read the data 
    df = pd.read_excel(r'weather station data folder\%s' %f)
    city.append(f[:-5])
    # simplize the data table 
    ndf = simptable(df)
    # 
    ndf.to_csv(r'simplized weather station data folder\%s' %(f[:-5]+'.csv') ,sep=',')
    Lat.append(ndf['Lat'][1])
    Long.append(ndf['Long'][1])
    
    # Deal With '-' error 
    ndf['NewT'] = ndf['Temp'].apply(filtstr)
    ndf['NewP'] = ndf['Precip'].apply(filtstr)
    ndf['NewW'] = ndf['Wind'].apply(filtstr)
    
    # remove the missing value 
    ndf['NewP'] = fixmissing(ndf['NewP'])
    ndf['NewW'] = fixmissing(ndf['NewW'])
    ndf['NewT'] = fixmissing(ndf['NewT'])
    ndf['Flyornot'] = ndf[['NewT','NewW',"NewP"]].apply(weatherfilter,axis=1)
    
    i = 0 
    C= ndf['Flyornot']
    day = []
    templist = []
	# Merge the 5-mins time steps into 3-hours  
    for c in C:
        templist.append(c)
        if len(templist) == 36:
            dronehour = consthreeh(templist)
            # 3 hour = 36 5-mins
            if dronehour == 36:
                day.append(1)
            else:
                day.append(0)
            templist = []
    daystore.append(day)
	
    # determine the consecutive dronehour 
    d =[]
    templist =[]
    for hour in day:
        templist.append(hour)

        if len(templist) == 8:

            dronehour = consthreeh(templist)
            # 3 hour = 36 5-mins

            if dronehour >= 1:
                d.append(1)
            else:
                d.append(0)

            templist = []
    
    droneday.append(sum(d))
    print (index)
    index += 1 
    

