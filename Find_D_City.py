import pandas as pd
import numpy as np
import netCDF4 as nc
from os import listdir
from os.path import isfile, join

# read the city excel
city = pd.read_excel('NewCity.xls')
print city.head()

c_lat = city['lat']
c_lon = city['lng']

# Dronedays file directory
mypath = r'C:\Users\mozhou\Desktop\UAV_Global_Mapping\Tables\Table_1'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

temp = np.zeros(11*721*1440).reshape(11,721,1440)
index = 0
for file in onlyfiles:

    dfile = nc.Dataset(r'C:\Users\mozhou\Desktop\UAV_Global_Mapping\Tables\Table_1\%s' % file,'r')

    d = dfile.variables['droneday'][:]
    temp[index,:,:] = d
    index = index + 1

# calculate the average dronedays for all types of drone
Ave = np.mean(temp,axis = 0)

# create coordinator list
c_coor = zip(c_lat,c_lon)

index = 0
c_freq = []
c_coor = list(c_coor)
for coor in c_coor:
    lat = coor[0]
    lon = coor[1]

    lai =np.round((90-lat)/0.25)
    i = int(lai)
    loi = np.round((lon-(-180))/0.25)
    j = int(loi)

    freq= Ave[i-1,j-1]
    c_freq.append(freq)
    print index
    index+=1

city['Frequency_P_1'] = c_freq
print city.head()
city.to_csv(r'C:\Users\mozhou\Desktop\UAV_Global_Mapping\Tables\Freq2.csv',sep =',')
