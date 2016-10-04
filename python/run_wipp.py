import numpy as np
# import pandas as pd
# import pickle
from scipy import interpolate
# import matplotlib.pyplot as plt
import os
# import itertools
# import random
import time
import datetime as dt
from index_helpers import load_TS_params
from index_helpers import load_Dst
from index_helpers import load_Kp
from index_helpers import load_ae
from spacepy import coordinates as coord
from spacepy.time import Ticktock

import bisect

project_root = '/shared/users/asousa/WIPP/3dWIPP/'


R_E = 6378.0    # km

# Run the WIPP code!
# Simulation time
ray_datenum = dt.datetime(2010, 3, 7, 11, 50, 00);

# yearday = '%d%03d'%(ray_datenum.year, ray_datenum.timetuple().tm_yday)
# milliseconds_day = (ray_datenum.second + ray_datenum.minute*60 + ray_datenum.hour*60*60)*1e3

iyr = ray_datenum.year
idoy= ray_datenum.timetuple().tm_yday
isec = (ray_datenum.second + ray_datenum.minute*60 + ray_datenum.hour*60*60)

inp_rayfile = '/shared/users/asousa/WIPP/3dWIPP/outputs/8_adjacent_model3_damped.ray'
# Get closest Kp value (Or should we interpolate?)
tvec, kvec = load_Kp()

tt = bisect.bisect_left(tvec, ray_datenum)
Kp = kvec[tt]
# ii = interpolate.interp1d(tvec, kvec)
# Kp = ii(ray_datenum)

# Get closest AE value
tvec, avec = load_ae()
tt = bisect.bisect_left(tvec, ray_datenum)
AE = np.log10(avec[tt])


print "year: ", iyr
print "day: ", idoy
print "sec: ", isec

inp_lats = [17.]
inp_lons = [32.]
launch_alt = R_E*1e3;

freqs    = np.array([1000]) 

inp_w = 2.0*np.pi*freqs


lats, lons, ws = np.meshgrid(inp_lats, inp_lons, inp_w)
lats = lats.flatten()
lons = lons.flatten()
ws   = ws.flatten()
alts = launch_alt*np.ones_like(lats)

# Create spacepy coordinate structures
inp_coords = coord.Coords(zip(alts, lons, lats), 'GEO', 'sph', units=['m','deg','deg'])
inp_coords.ticks = Ticktock(np.tile(ray_datenum, len(inp_coords))) # add ticks

print "geo spherical:"
print inp_coords.data
# Rotate to SM cartesian coordinates
inp_coords = inp_coords.convert('GEO','car')
print "geo cartesian:"
print inp_coords.data
inp_coords = inp_coords.convert('SM','car')
print "SM cartesian:"
print inp_coords.data
# Build wipp code
os.chdir(project_root)

print '-------- Building WIPP ------'
os.system('make')

# run it
wipp_cmd = 'bin/wipp -i %s -t %s -u %d -v %d -a %e -b %e -c %e'%(inp_rayfile, iyr, idoy, isec,
                inp_coords.data[0,0], inp_coords.data[0,1], inp_coords.data[0,2])
print wipp_cmd

os.system(wipp_cmd)
