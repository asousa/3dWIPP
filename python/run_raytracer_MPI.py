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

import xflib  # Fortran xform-double library (coordinate transforms)

# from spacepy import coordinates as coord
# from spacepy.time import Ticktock

import bisect   # fast searching of sorted lists

import commands
import subprocess
from mpi4py import MPI


project_root = '/shared/users/asousa/WIPP/3dWIPP/'
raytracer_root = '/shared/users/asousa/software/foust_raytracer/'
damping_root = '/shared/users/asousa/WIPP/3dWIPP/damping/'
ray_bin_dir    = os.path.join(raytracer_root, 'bin')
ray_out_dir = '/shared/users/asousa/WIPP/3dWIPP/outputs/1lon_ngo_10_65'

R_E = 6371.0    # km

# ----------- Simulation params ----------------
t_max = 10.     # Maximum duration in seconds

dt0 = 1e-4      # Initial timestep in seconds
dtmax = 1e-1    # Maximum allowable timestep in seconds
root = 2        # Which root of the Appleton-Hartree equation
                # (1 = negative, 2 = positive)
                # (2=whistler in magnetosphere)
fixedstep = 0   # Don't use fixed step sizes, that's a bad idea.
maxerr = 1e-5   # Error bound for adaptive timestepping
maxsteps = 1e5  # Max number of timesteps (abort if reached)
modelnum = 1    # Which model to use (1 = ngo, 2=GCPM, 3=GCPM interp, 4=GCPM rand interp)
use_IGRF = 1    # Magnetic field model (1 for IGRF, 0 for dipole)
use_tsyg = 0    # Use the Tsyganenko magnetic field model corrections

minalt   = (R_E + 100)*1e3 # cutoff threshold in meters

# ---------- Ray inputs -----------------

# Geomagnetic please.
inp_lats = np.arange(10, 65, 1) #[40, 41, 42, 43]
inp_lons = [0, 1]
launch_alt = (R_E + 1000.)*1e3

# freqs    = np.array([1000, 2000]) 

# freqs = [200,240,289,347,418,502,603,725,872,1048,1259,1514,1819,2187,2629,3160,3798,4565,5487,6596,7928,9530,11455,13769,16550,19893,23912,28742,34549,41528,49916,60000]

f1 = 200; f2 = 30000;
num_freqs = 32
flogs = np.linspace(np.log10(f1), np.log10(f2), num_freqs)
freqs = np.round(pow(10, flogs)/10.)*10


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

# -------------- set up MPI -----------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")


# Split frequency vector into smaller chunks, pass each chunk to a process
nProcs = 1.0*comm.Get_size()
nFreqs = 1.0*np.shape(freqs)[0]
nSteps = np.ceil(nFreqs/nProcs).astype(int)


# Shuffle the frequency vector (adjacent frequencies take about as long to run)
# np.random.shuffle(freqs)

chunks = [freqs[i:i+nSteps] for i in range(0, len(freqs), nSteps)]

if rank==0:
    print "We have %d processes available"%(nProcs)
    # Clean up previous files:
    if (not os.path.exists(ray_out_dir)):
        os.mkdir(ray_out_dir);
    if (not os.path.exists(os.path.join(ray_out_dir, "logs"))):
        os.mkdir(os.path.join(ray_out_dir,"logs"))


# Load solar wind parameters (for Tsykadenko corrections)
Pdyn, ByIMF, BzIMF, W = load_TS_params(ray_datenum)

# Load Dst
Dst = load_Dst(ray_datenum)

yearday = '%d%03d'%(ray_datenum.year, ray_datenum.timetuple().tm_yday)
milliseconds_day = (ray_datenum.second + ray_datenum.minute*60 + ray_datenum.hour*60*60)*1e3 + ray_datenum.microsecond*1e-3

# Coordinate transformation library
xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')

# # Each subprocess does a subset of frequencies
if (rank < len(chunks)):
    for freq in chunks[rank]:
        print "Subprocess %s on %s: doing frequency %g"%(rank, host, freq) 



        inp_w = 2.0*np.pi*freq


        lats, lons, ws = np.meshgrid(inp_lats, inp_lons, inp_w)
        lats = lats.flatten()
        lons = lons.flatten()
        ws   = ws.flatten()
        alts = launch_alt*np.ones_like(lats)

        # VERY UNSCIENTIFIC BAND-AID:
        # Launch rays below 600 hz at 4000 km instead of 1000.
        alts[ws < 600*2*np.pi] += 3000e3



        # Create coordinates
        inp_coords = zip(alts, lats, lons)  # Geomagnetic pls.

        print "Frequency = ", ws/(2.*np.pi)
        print "Inputs (geomagnetic RLL)"
        for r in inp_coords: print r




        working_path = os.path.join(os.path.expanduser("~"),"rayTmp_%d"%(freq))
        ray_inpfile = os.path.join(working_path,'ray_inputs.txt')
        ray_outfile = os.path.join(ray_out_dir, 'ray_%d.ray'%(freq))
        damp_outfile = os.path.join(ray_out_dir,'damp_%d.ray'%(freq))
        # dumpfile    = os.path.join(project_root,'output','dumpout.txt')
        


        if (not os.path.exists(working_path)):
           os.mkdir(working_path)


        print "Cleaning previous runs..."
        if os.path.exists(ray_inpfile):
            os.remove(ray_inpfile)
        if os.path.exists(ray_outfile):
            os.remove(ray_outfile)
        if os.path.exists(damp_outfile):
            os.remove(damp_outfile)

        # Rotate from geomagnetic to SM cartesian coordinates
        inp_coords = [xf.rllmag2sm(r, ray_datenum) for r in inp_coords]

        # inp_coords = inp_coords.convert('SM','car')
        N = len(inp_coords)

        # Write rays to the input file (used by the raytracer):
        f = open(ray_inpfile,'w+')
        for pos0, w0 in zip(inp_coords, ws):        
            dir0 = pos0/np.linalg.norm(pos0)    # radial outward
            f.write('%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n'%(pos0[0], pos0[1], pos0[2], dir0[0], dir0[1], dir0[2], w0))
        f.close()



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
        # os.system(cmd)
        runlog = subprocess.check_output(cmd, shell=True)
        file = open(os.path.join(ray_out_dir, 'logs', "ray%g.log"%(freq)),'w+')
        file.write(runlog)
        file.close()


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

        damp_cmd = '%sbin/damping -i %s -o %s -k %g -a %g -m %d -t %s -u %d'%(damping_root, ray_outfile, damp_outfile, Kp, AE, damp_mode, yearday, milliseconds_day)
        print damp_cmd

        # Start the damping code
        # os.system(damp_cmd)
        damplog = subprocess.check_output(damp_cmd, shell=True)
        file = open(os.path.join(ray_out_dir, "logs", "damp%g.log"%(freq)),'w+')
        file.write(damplog)
        file.close()

        # # remove undamped file
        # os.remove(ray_outfile)
