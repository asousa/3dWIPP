import numpy as np
from scipy import interpolate
# import matplotlib.pyplot as plt
import os
import itertools
from partition import partition
import time
import datetime as dt
import sys
from index_helpers import load_TS_params
from index_helpers import load_Dst
from index_helpers import load_Kp
from index_helpers import load_ae
# from spacepy import coordinates as coord
# from spacepy.time import Ticktock

import xflib  # Fortran xform-double library (coordinate transforms)
import bisect

import commands
import subprocess
# from mpi4py import MPI



project_root = '/shared/users/asousa/WIPP/3dWIPP/'
flux_file = os.path.join(project_root,'data','AE8MaxFlux_expanded.dat')
# flux_file = os.path.join(project_root,'data','EQFLUXMA.dat')

# ---------------------- Constants -------------------------
R_E = 6371.0    # km
R2D = 180./np.pi
D2R = np.pi/180.


# ------------------ Simulation params ---------------------

# Simulation time
ray_datenum = dt.datetime(2010, 06, 04, 07, 00, 00);

# Flash location
inp_lat = 30
inp_lon = 0
launch_alt = ((R_E + 5)*1e3)/R_E;
flash_I0 = -100e3

# Frequencies
f1 = 200; f2 = 30000;
num_freqs = 33
flogs = np.linspace(np.log10(f1), np.log10(f2), num_freqs)
freqs = np.round(pow(10, flogs)/10.)*10
freq_pairs = zip(freqs[0:], freqs[1:])

# Output coordinates (geomagnetic)
# out_lat = [40, 50]
# Select latitudes for a uniform spread across L-shells
# out_Lsh = np.arange(1.2, 8.0, 0.1)
# out_lat = np.round(10.0*np.arccos(1.0/out_Lsh)*R2D)/10.0
out_lat = [66.4]

out_lon = [0]

model_number = 0        # b-field model (0 = dipole, 1 = IGRF)
num_freq_steps = 20     # number of interpolating steps between 
                        # each guide frequency.

flux_dist = 0
alpha_dist = 0

vec_ind = 0     # Which set of default params to use for the gcpm model

ray_input_directory = '/shared/users/asousa/WIPP/rays/2d/nightside/mode6/kp0/'
output_directory    = os.path.join(project_root, "outputs", "main_2d_test")
log_directory       = os.path.join(output_directory, "logs")

# ----------------------------------------------------------

iyr = ray_datenum.year
idoy= ray_datenum.timetuple().tm_yday 
isec = (ray_datenum.second + (ray_datenum.minute)*60 + ray_datenum.hour*60*60)

# Flash input coordinates:
inp_coords = [launch_alt, inp_lat, inp_lon]

for olat in out_lat:
    olon = out_lon[0]
    flux_cmd = '%sbin/flux --inp_dir %s'%(project_root, output_directory) +\
               ' --out_dir %s'%(output_directory) +\
               ' --flux_file %s --lat %g --lon %g'%(flux_file, olat, olon) +\
               ' --flux_dist %d --alpha_dist %d'%(flux_dist, alpha_dist)

    for freq in freqs[:-1]:
        flux_cmd += ' --f %g'%freq

    print flux_cmd
    os.system(flux_cmd)




