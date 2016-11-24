import numpy as np
# import pandas as pd
# import pickle
from scipy import interpolate
# import matplotlib.pyplot as plt
import os
import shutil
# import itertools
# import random
import time
import datetime as dt
from index_helpers import load_TS_params
from index_helpers import load_Dst
from index_helpers import load_Kp
from index_helpers import load_ae

from distutils.dir_util import mkpath

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
ray_out_dir = '/shared/users/asousa/WIPP/3dWIPP/outputs/rays_dipole_florida'

R_E = 6371.0    # km

# ----------- Simulation params ----------------
t_max = 3.     # Maximum duration in seconds

dt0 = 1e-4      # Initial timestep in seconds
dtmax = 1e-1    # Maximum allowable timestep in seconds
root = 2        # Which root of the Appleton-Hartree equation
                # (1 = negative, 2 = positive)
                # (2=whistler in magnetosphere)
fixedstep = 0   # Don't use fixed step sizes, that's a bad idea.
maxerr = 1e-5   # Error bound for adaptive timestepping
maxsteps = 1e5  # Max number of timesteps (abort if reached)
modelnum = 1    # Which model to use (1 = ngo, 2=GCPM, 3=GCPM interp, 4=GCPM rand interp)
use_IGRF = 0    # Magnetic field model (1 for IGRF, 0 for dipole)
use_tsyg = 0    # Use the Tsyganenko magnetic field model corrections

minalt   = (R_E + 100)*1e3 # cutoff threshold in meters

# ---------- Ray inputs -----------------

#  7am UTC at 0 magnetic longitude should be about midnight in the Florida gulf

# Geomagnetic please.
inp_lats = np.arange(30, 60, 5) #[40, 41, 42, 43]
inp_lons = np.arange(0, 1, 1)
# freqs    = [200, 2000]
launch_alt = (R_E + 1000)*1e3
# freqs    = np.array([1000, 2000]) 

# freqs = [200,240,289,347,418,502,603,725,872,1048,1259,1514,1819,2187,2629,3160,3798,4565,5487,6596,7928,9530,11455,13769,16550,19893,23912,28742,34549,41528,49916,60000]

f1 = 200; f2 = 30000;
num_freqs = 32
flogs = np.linspace(np.log10(f1), np.log10(f2), num_freqs)
freqs = np.round(pow(10, flogs)/10.)*10

lats, lons, fs = np.meshgrid(inp_lats, inp_lons, freqs)
lats = lats.flatten()
lons = lons.flatten()
fs   = fs.flatten()

alts = launch_alt*np.ones_like(lats)
alts[fs < 600] += 3000e3



# Simulation time
ray_datenum = dt.datetime(2010, 06, 04, 07, 00, 00);

# Damping parameters:
damp_mode = 1  # 0 for old 2d damping code, 1 for modern code


# -------------------------------------------

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


# -------------- set up MPI -----------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")

# print "hi from %s %d"%(host, rank)


inps = zip(alts, lats, lons, fs)



# Split frequency vector into smaller chunks, pass each chunk to a process
nProcs = 1.0*comm.Get_size()
nInps  = 1.0*np.shape(inps)[0]
nSteps = np.ceil(nInps/nProcs).astype(int)

chunks = [inps[i:i+nSteps] for i in range(0, int(nInps), nSteps)]


if rank==0:
    print "We have %d processes available"%(nProcs)
    print "we have %d inputs to do"%(nInps)
    print "Each node gets %d jobs"%nSteps
    # Clean up previous files:
    os.system("rm %s -r"%ray_out_dir)
    os.mkdir(ray_out_dir)
    os.mkdir(os.path.join(ray_out_dir,"logs"))

    # build output file tree
    for f in freqs:
        for lo in inp_lons:
            ray_out_subdir = os.path.join(ray_out_dir, "f_%d"%f, "lon_%d"%lo)
            if (not os.path.exists(ray_out_subdir)):
                mkpath(ray_out_subdir)



# Load solar wind parameters (for Tsykadenko corrections)
Pdyn, ByIMF, BzIMF, W = load_TS_params(ray_datenum)

# Load Dst
Dst = load_Dst(ray_datenum)

yearday = '%d%03d'%(ray_datenum.year, ray_datenum.timetuple().tm_yday)
milliseconds_day = (ray_datenum.second + ray_datenum.minute*60 + ray_datenum.hour*60*60)*1e3 + ray_datenum.microsecond*1e-3

# Coordinate transformation library
xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')


# # Sync up:
comm.Barrier();

# # Each subprocess does a subset of frequencies
if (rank < len(chunks)):
    working_path = os.path.join(os.path.expanduser("~"),"rayTmp_%d"%(rank))

    if (not os.path.exists(working_path)):
        mkpath(working_path)

    print "Subprocess %s on %s: doing %d rays"%(rank, host, len(chunks[rank]))
    for inp in chunks[rank]:

        print "Subprocess %s on %s: doing ray from (%d, %d) at %d hz"%(rank, host, inp[1], inp[2], inp[3]) 

        ray_out_subdir = os.path.join(ray_out_dir, "f_%d"%inp[3], "lon_%d"%(inp[2]))


        ray_inpfile = os.path.join(working_path,'ray_inputs_%d_%d_%d.txt'%(inp[3], inp[1], inp[2]))
        ray_outfile = os.path.join(ray_out_subdir, 'ray_%d_%d_%d.ray'%(inp[3], inp[1], inp[2]))
        damp_outfile = os.path.join(ray_out_subdir,'damp_%d_%d_%d.ray'%(inp[3], inp[1], inp[2]))
        # dumpfile    = os.path.join(project_root,'output','dumpout.txt')
        


        # # print "Cleaning previous runs..."
        # if os.path.exists(ray_inpfile):
        #     os.remove(ray_inpfile)

        # Rotate from geomagnetic to SM cartesian coordinates
        # print "inp: ", inp
        inp_coords = xf.rllmag2sm(inp, ray_datenum)
        # print "inp_coords (sm): ", inp_coords

        # Write ray to the input file (used by the raytracer):
        f = open(ray_inpfile,'w')
        # for pos0, w0 in zip(inp_coords, ws):        
        pos0 = inp_coords
        # print "pos0: ", pos0
        dir0 = pos0/np.linalg.norm(pos0)    # radial outward
        w0   = inp[3]*2.0*np.pi
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

        # print cmd
        # Start the raytracer
        # os.system(cmd)
        runlog = subprocess.check_output(cmd, shell=True)
        file = open(os.path.join(ray_out_dir, 'logs', "ray_%g_%g_%g.log"%(inp[3], inp[1], inp[2])),'w+')
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
        # print damp_cmd

        # Start the damping code
        # os.system(damp_cmd)
        # damplog = subprocess.check_output(damp_cmd, shell=True)
        # file = open(os.path.join(ray_out_dir, "logs", "damp_%g_%g_%g.log"%(inp[3], inp[1], inp[2])),'w+')
        # file.write(damplog)
        # file.close()
 
    # shutil.rmtree(working_path)

