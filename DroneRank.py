import pandas as pd
import numpy as np
import netCDF4 as nc
from os import listdir
from os.path import isfile, join

mypath = r'C:\Users\mozhou\Desktop\UAV_Global_Mapping\Tables\Table_0'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

x = 721
y = 1440

# concatenate files
temp = np.zeros(12*x*y).reshape(12,x,y)
index = 0
i = 0
for file in onlyfiles:

    dfile = nc.Dataset(r'C:\Users\mozhou\Desktop\UAV_Global_Mapping\Tables\Table_0\%s' % file,'r')

    d = dfile.variables['droneday'][:]
    temp[index,:,:] = d

    index = index +1
    dfile.close()

# create name list
name = ['Albris','Bramor','Phantom4','Ebee','H920','M100','M600','pixhawk2','Q500','SkyRanger','TyphoneH','X4p']
z = 12
k = 0
rank = []
while k<z:
    a = temp[k,:,:]
    mean = np.mean(a)

    rank.append(mean)

    k += 1

print rank,name

d = {'Products': name, 'Global_Ave_Freq': np.round(rank)}
df = pd.DataFrame(data=d)
ndf = df.sort_values(by=['Global_Ave_Freq'],ascending = False)

print ndf

ndf.to_csv(r'C:\Users\mozhou\Desktop\UAV_Global_Mapping\Tables\Freq2.csv',sep =',')
