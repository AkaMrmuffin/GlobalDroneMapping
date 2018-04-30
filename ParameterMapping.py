import pandas as pd
import numpy as np
import netCDF4 as nc
from os import listdir
from os.path import isfile, join
from scipy import stats

# Read weather data
mypath1 = r'Y:\Mo_GlobalUAVMapping\High_res_data\Wind_Temp'
Windfiles = [f for f in listdir(mypath1) if isfile(join(mypath1, f))]
mypath2 = r'Y:\Mo_GlobalUAVMapping\High_res_data\Precip'
Precipfiles = [f for f in listdir(mypath2) if isfile(join(mypath2, f))]

files = zip(Windfiles,Precipfiles)
files = list(files)

print files

# Read Drone Specs table
df = pd.read_excel(r'C:\Users\mozhou.gao\Desktop\NewDroneTable.xlsx')
print df.head()

for index,row in df.iterrows():
    # Loading the weather thresholds
    name = row['Drones']
    t1 = row['TempL']
    t2 = row['TempU']
    ws = row['Wind ']
    p = row['Precip'] ## for precipitation > 1, change Precip to Precip.1

    print t1,t2,ws,p

    for z in files:
        wtname = z[0]
        year  = wtname[6:10]
        print year
        pname = z[1]
        pfile = nc.Dataset(r"Y:\Mo_GlobalUAVMapping\High_res_data\Wind_Temp\%s" % wtname, 'r')
        wtfile = nc.Dataset(r"Y:\Mo_GlobalUAVMapping\High_res_data\Precip\%s" % pname, 'r')
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

        else:
            y = 365

        print "File loaded"

        # Crate empty Boolean Matrix
        boolT = np.zeros(t*la*lo).reshape(t,la,lo)
        boolV = np.zeros(t*la*lo).reshape(t,la,lo)
        boolP = np.zeros(t*la*lo).reshape(t,la,lo)

        # Create Time index
        index = np.arange(t)

        for i in index:
            u = wtfile.variables['u10'][i, :, :]
            v = wtfile.variables['v10'][i, :, :]
            temp = wtfile.variables['t2m'][i, :, :]
            prep = pfile.variables['tp'][i, :, :]
            wind = (u ** 2 + v ** 2) ** 0.5
            temp = temp - 273.15

            # Wind
            wa = wind <= ws
            wb = wind > ws
            wind[wa] = 0
            wind[wb] = 1

            # Precipitation
            pa = prep <= p
            pb = prep > p
            prep[pa] = 0
            prep[pb] = 1

            # Temperature
            ta = (temp <= t2) & (temp >= t1)
            tb = (temp > t2) | (temp < t1)
            temp[ta] = 0
            temp[tb] = 1

            #boolean
            boolT[i,:,:] = temp
            boolV[i,:,:] = wind
            boolP[i,:,:] = prep

        # Calculate the sum
        Tc = np.sum(boolT,axis = 0)
        Wc = np.sum(boolV,axis = 0)
        Pc = np.sum(boolP,axis = 0)

        # Create empty Color array
        Color = np.zeros(la*lo).reshape(la,lo)
        i = 0
        while i<la:
            j = 0
            while j<lo:

                tt = Tc[i,j]
                ww = Wc[i,j]
                pp = Pc[i,j]

                L = [tt,ww,pp]

                c = max(L)

                if c == ww:
                    Color[i,j] =1
                if c == tt:
                    Color[i,j] = 2
                if c == pp:
                    Color[i,j] = 3

                j += 1

            i += 1

        wtfile.close()
        pfile.close()
        N = name + year
        dronedata =nc.Dataset(r"C:\Users\mozhou.gao\Desktop\%s.nc" % N, 'w', format='NETCDF4_CLASSIC')
        lat = dronedata.createDimension('lat', la)
        lon = dronedata.createDimension('lon', lo)
        dt = dronedata.createVariable('color', np.int32, ( 'lat', 'lon'))
        dt[:] = Color
        dronedata.close()







