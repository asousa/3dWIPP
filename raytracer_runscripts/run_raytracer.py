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

import bisect   # fast searching of sorted lists


project_root = '/shared/users/asousa/WIPP/3dWIPP/'
raytracer_root = '/shared/users/asousa/software/foust_raytracer/'
damping_root = '/shared/users/asousa/WIPP/3dWIPP/damping/'
ray_bin_dir    = os.path.join(raytracer_root, 'bin')

R_E = 6371.0    # km

# ----------- Simulation params ----------------
t_max = 10.     # Maximum duration in seconds

dt0 = 0.01      # Initial timestep in seconds
dtmax = 0.01     # Maximum allowable timestep in seconds
root = 2        # Which root of the Appleton-Hartree equation
                # (2=whistler in magnetosphere)
fixedstep = 0   # Don't use fixed step sizes, that's a bad idea.
maxerr = 1e-4   # Error bound for adaptive timestepping
maxsteps = 1e5  # Max number of timesteps (abort if reached)
modelnum = 3    # Which model to use (1 = ngo, 2=GCPM, 3=GCPM interp, 4=GCPM rand interp)
use_IGRF = 1    # Magnetic field model (1 for IGRF, 0 for dipole)
use_tsyg = 1    # Use the Tsyganenko magnetic field model corrections

minalt   = (R_E + 100)*1e3 # cutoff threshold in meters

# GCPM grid to use (plasmasphere model)
if modelnum==1:
    configfile = '/shared/users/asousa/WIPP/3dWIPP/outputs/newray.in'
# interpfile = os.path.join(project_root,'raytracer_runscripts','gcpm_models','gcpm_kp40_20010101_0000_MLD01.txt')
if modelnum==3:
    interpfile = '/shared/users/asousa/software/foust_raytracer/bin/gcpm_kp4_2001001_L10_80x80x80_noderiv.txt'
if modelnum==4:
    interpfile = '/shared/users/asousa/software/foust_raytracer/bin/gcpm_kp4_2001001_L10_random_5000_20000_0_200000_600000.txt'
    scattered_interp_window_scale = 1.5
    scattered_interp_order = 2
    scattered_interp_exact = 0
    scattered_interp_local_window_scale = 5


# Simulation time
ray_datenum = dt.datetime(2010, 06, 04, 03, 17, 00);

# Damping parameters:
damp_mode = 1  # 0 for old 2d damping code, 1 for modern code

# Change this once you set up parallel stuff!
ray_inpfile = os.path.join(ray_bin_dir,'ray_inputs.txt')
ray_outfile = os.path.join(project_root,'outputs','rayout.ray')
damp_outfile = os.path.join(project_root,'outputs','rayout_damped.ray')
dumpfile    = os.path.join(project_root,'output','dumpout.txt')

# Clean up previous files:
print "Cleaning previous runs..."
if os.path.exists(ray_inpfile):
    os.remove(ray_inpfile)
if os.path.exists(ray_outfile):
    os.remove(ray_outfile)
if os.path.exists(dumpfile):
    os.remove(dumpfile)

# Load solar wind parameters (for Tsykadenko corrections)
Pdyn, ByIMF, BzIMF, W = load_TS_params(ray_datenum)

# Load Dst
Dst = load_Dst(ray_datenum)
# print Dst

# ---------- Ray inputs -----------------

# inp_lats = [45, 46]
# inp_lons = [180, 181] #np.arange(0,360, step=45)
# launch_alt = (R_E + 2000.)*1e3 # meters
inp_lats = [45, 46]
inp_lons = [179, 180]
launch_alt = (R_E + 2000.)*1e3

freqs    = np.array([1000, 1100]) 

inp_w = 2.0*np.pi*freqs


lats, lons, ws = np.meshgrid(inp_lats, inp_lons, inp_w)
lats = lats.flatten()
lons = lons.flatten()
ws   = ws.flatten()
alts = launch_alt*np.ones_like(lats)

# Create spacepy coordinate structures
inp_coords = coord.Coords(zip(alts, lats, lons), 'GEO', 'sph', units=['m','deg','deg'])
inp_coords.ticks = Ticktock(np.tile(ray_datenum, len(inp_coords))) # add ticks

# Rotate to SM cartesian coordinates
inp_coords = inp_coords.convert('SM','car')
N = len(inp_coords)

# Write rays to the input file
f = open(ray_inpfile,'w+')
for pos0, w0 in zip(inp_coords.data, ws):
    dir0 = pos0/np.linalg.norm(pos0)    # radial outward
    f.write('%g %g %g %g %g %g %g\n'%(pos0[0],pos0[1],pos0[2],dir0[0], dir0[1],dir0[2], w0))
f.close()

yearday = '%d%03d'%(ray_datenum.year, ray_datenum.timetuple().tm_yday)
milliseconds_day = (ray_datenum.second + ray_datenum.minute*60 + ray_datenum.hour*60*60)*1e3  \
                +  ray_datenum.microsecond*1e-3

# --------- Build the run command --------

cmd= '%s/raytracer --outputper=%d --dt0=%g --dtmax=%g'%(ray_bin_dir, 1, dt0, dtmax) + \
     ' --tmax=%g --root=%d --fixedstep=%d --maxerr=%g'%(t_max, root, fixedstep, maxerr) + \
     ' --maxsteps=%d --minalt=%d --inputraysfile=%s --outputfile=%s'%( maxsteps, minalt, ray_inpfile, ray_outfile) + \
     ' --modelnum=%d --yearday=%s --milliseconds_day=%d'%(modelnum, yearday, milliseconds_day) + \
     ' --use_tsyganenko=%d --use_igrf=%d --tsyganenko_Pdyn=%g'%(use_tsyg, use_IGRF, Pdyn) + \
     ' --tsyganenko_Dst=%g --tsyganenko_ByIMF=%g --tsyganenko_BzIMF=%g'%( Dst, ByIMF, BzIMF ) + \
     ' --tsyganenko_W1=%g --tsyganenko_W2=%g --tsyganenko_W3=%g'%(W[0], W[1], W[2]) + \
     ' --tsyganenko_W4=%g --tsyganenko_W5=%g --tsyganenko_W6=%g'%(W[3], W[4], W[5])

# Append model-specific parameters to the command line
if  modelnum == 1:
    cmd += ' --ngo_configfile=%s'%configfile
elif modelnum == 2:
    cmd += ' --gcpm_kp=%g'%kp
elif modelnum == 3: 
    cmd += ' --interp_interpfile=%s'%interpfile
elif modelnum == 4:
    cmd += ' --interp_interpfile=%s'%interpfile
    cmd += ' --scattered_interp_window_scale=%g'%scattered_interp_window_scale
    cmd += ' --scattered_interp_order=%d'%scattered_interp_order
    cmd += ' --scattered_interp_exact=%d'%scattered_interp_exact
    cmd += ' --scattered_interp_local_window_scale=%g'%scattered_interp_local_window_scale

print cmd
# Start the raytracer
os.system(cmd)


# ------- Run Damping Code ------------

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

# print Kp, AE

damp_cmd = '%sbin/damping -i %s -o %s -k %g -a %g -m %d'%(damping_root, ray_outfile, damp_outfile, Kp, AE, damp_mode)
print damp_cmd

# Start the damping code
os.system(damp_cmd)

