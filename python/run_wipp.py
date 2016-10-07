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


R_E = 6371.0    # km

# Simulation time
ray_datenum = dt.datetime(2000, 01, 1, 0, 0, 00);

iyr = ray_datenum.year
idoy= ray_datenum.timetuple().tm_yday 
isec = (ray_datenum.second + (ray_datenum.minute)*60 + ray_datenum.hour*60*60)

inp_rayfile = '/shared/users/asousa/WIPP/3dWIPP/outputs/rayout_damped.ray'
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

inp_lats = [45]
inp_lons = [17]
launch_alt = 1;

freqs    = np.array([23]) 

inp_w = 2.0*np.pi*freqs


lats, lons, ws = np.meshgrid(inp_lats, inp_lons, inp_w)
lats = lats.flatten()
lons = lons.flatten()
ws   = ws.flatten()
alts = launch_alt*np.ones_like(lats)

# Create spacepy coordinate structures
inp_coords = coord.Coords(zip(alts, lats, lons), 'GEO', 'sph', units=['Re','deg','deg'])
inp_coords.ticks = Ticktock(np.tile(ray_datenum.isoformat(), len(inp_coords)),'ISO') # add ticks

os.chdir(project_root)

print '-------- Building WIPP ------'
os.system('make')

# run it
wipp_cmd = 'bin/wipp -i %s -t %s -u %d -v %d -a %g -b %g -c %g'%(inp_rayfile, iyr, idoy, isec,
                inp_coords.radi[0], inp_coords.lati[0], inp_coords.long[0])
print wipp_cmd


print "geo spherical (py):"
print inp_coords.data
# inp_coords = inp_coords.convert('MAG','car')
# Rotate to SM cartesian coordinates
inp_coords = inp_coords.convert('GEO','car')
print "geo cartesian (py):"
print inp_coords.data
inp_coords = inp_coords.convert('SM','car')
print "SM cartesian (py):"
print inp_coords.data
# Build wipp code

os.system(wipp_cmd)
