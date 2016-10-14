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
import sys
from index_helpers import load_TS_params
from index_helpers import load_Dst
from index_helpers import load_Kp
from index_helpers import load_ae
from spacepy import coordinates as coord
from spacepy.time import Ticktock

import xflib  # Fortran xform-double library (coordinate transforms)

import bisect

project_root = '/shared/users/asousa/WIPP/3dWIPP/'


R_E = 6371.0    # km

# Simulation time
ray_datenum = dt.datetime(2010, 06, 04, 03, 17, 00);

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
inp_lons = [180]
launch_alt = (R_E + 5)*1e3;

# freqs    = np.array([23]) 

# inp_w = 2.0*np.pi*freqs


lats, lons = np.meshgrid(inp_lats, inp_lons)
lats = lats.flatten()
lons = lons.flatten()
alts = launch_alt*np.ones_like(lats)

# # Create spacepy coordinate structures
# inp_coords = coord.Coords(zip(alts, lats, lons), 'GEO', 'sph', units=['Re','deg','deg'])
# inp_coords.ticks = Ticktock(np.tile(ray_datenum.isoformat(), len(inp_coords)),'ISO') # add ticks


# Create coordinates
inp_coords = zip(alts, lats, lons)  # Geographic

# Coordinate transformation library
xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')

print "Inputs (geomagnetic RLL)"
for r in inp_coords: print r



os.chdir(project_root)

print '-------- Building WIPP ------'
buildstatus = os.system('make')

if buildstatus != 0:
    sys.exit("Build failed!")

# run it
wipp_cmd = 'bin/wipp -i %s -t %s -u %d -v %d -a %g -b %g -c %g'%(inp_rayfile, iyr, idoy, isec,
                inp_coords[0][0], inp_coords[0][1], inp_coords[0][2])
print wipp_cmd


print "geo spherical (py):"
print inp_coords

inp_coords = [xf.s2c(r) for r in inp_coords]
print "geo cartesian (py):"
for r in inp_coords: print r


inp_coords = [xf.geo2sm(r, ray_datenum) for r in inp_coords]
print "SM (py):"
for r in inp_coords: print r
# Build wipp code

os.system(wipp_cmd)
