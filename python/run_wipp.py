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

# ------------------ Simulation params ---------------------

# Simulation time
ray_datenum = dt.datetime(2010, 06, 04, 07, 00, 00);

# Flash location
inp_lat = 45
inp_lon = 0
launch_alt = ((R_E + 5)*1e3)/R_E;
flash_I0 = -10e3

freqs = [1010, 1180];

# Output coordinates (geomagnetic)
out_lat = 50
out_lon = 0

model_number = 1        # b-field model (0 = dipole, 1 = IGRF)

ray_input_directory = os.path.join(project_root, "outputs", "rays2")
output_directory    = os.path.join(project_root, "outputs", "test_WIPP_outs")

# ----------------------------------------------------------





iyr = ray_datenum.year
idoy= ray_datenum.timetuple().tm_yday 
isec = (ray_datenum.second + (ray_datenum.minute)*60 + ray_datenum.hour*60*60)


if not os.path.exists(output_directory):
    os.mkdir(output_directory)

print "Clearing data from previous runs..."
os.system('rm %s/*'%(output_directory))

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

# Load solar wind parameters (for Tsykadenko corrections)
Pdyn, ByIMF, BzIMF, W = load_TS_params(ray_datenum)

# Load Dst
Dst = load_Dst(ray_datenum)

TS_params = [Pdyn, Dst, ByIMF, BzIMF, W[0], W[1], W[2], W[3], W[4], W[5]]
# print TS_params

print "Pdyn: ", Pdyn
print "ByIMF: ", ByIMF
print "BzIMF: ", BzIMF
print "W: ", W

print "year: ", iyr
print "day: ", idoy
print "sec: ", isec

# Create coordinates
inp_coords = [launch_alt, inp_lat, inp_lon]

# Coordinate transformation library
xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')

print "input coords (geomagnetic RLL):"
print inp_coords

os.chdir(project_root)

print '-------- Building WIPP ------'
buildstatus = os.system('make')

if buildstatus != 0:
    sys.exit("Build failed!")


# run it
# wipp_cmd = ['bin/wipp -i %s -o %s'%(ray_input_directory, output_directory) +
#             ' -t %s -u %d -v %d'%(iyr, idoy, isec) +
#             ' -a %g -b %g -c %g'%(inp_coords[0], inp_coords[1], inp_coords[2]) +
#             ' -e %g -f %g -d 0'%(out_lat, out_lon) +
#             ' -g %s -h %s'%(rayfile1, rayfile2)][0]

wipp_cmd = ['bin/wipp --out_dir %s'%(output_directory) +
            ' --iyr %s --idoy %d --isec %d --I0 %d'%(iyr, idoy, isec, flash_I0) +
            ' --f_alt %g --f_lat %g --f_lon %g'%(inp_coords[0], inp_coords[1], inp_coords[2]) +
            ' --out_lat %g --out_lon %g --b_model %d'%(out_lat, out_lon, model_number) +
            ' --f1 %g --f2 %g --ray_dir %s'%(freqs[0], freqs[1], ray_input_directory)][0]
            # ' --low_file %s --hi_file %s'%(rayfile1, rayfile2)][0]



print wipp_cmd

os.system(wipp_cmd)
