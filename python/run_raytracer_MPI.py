from mpi4py import MPI
import numpy as np
from scipy import interpolate
import os
from partition import partition
import itertools
import random
import time
import datetime as dt
import commands
import subprocess
import shutil

from index_helpers import load_TS_params
from index_helpers import load_Dst
from index_helpers import load_Kp
from index_helpers import load_ae

import xflib  # Fortran xform-double library (coordinate transforms)
import bisect   # fast searching of sorted lists

# ------------ Start MPI -------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")
nProcs = 1.0*comm.Get_size()


# ------------------Constants --------------------------
R_E = 6371


# -------------- Simulation params ---------------------
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

# ---------------- Input parameters --------------------

inp_lats = np.arange(10, 60, 1) #[40, 41, 42, 43]
inp_lons = np.arange(-5, 6, 1)

launch_alt = (R_E + 1000)*1e3

f1 = 200; f2 = 30000;
num_freqs = 32
flogs = np.linspace(np.log10(f1), np.log10(f2), num_freqs)
freqs = np.round(pow(10, flogs)/10.)*10

# Simulation time
ray_datenum = dt.datetime(2010, 06, 04, 07, 00, 00);

# Damping parameters:
damp_mode = 1  # 0 for old 2d damping code, 1 for modern code

project_root = '/shared/users/asousa/WIPP/3dWIPP/'
raytracer_root = '/shared/users/asousa/software/foust_raytracer/'
damping_root = '/shared/users/asousa/WIPP/3dWIPP/damping/'
ray_bin_dir    = os.path.join(raytracer_root, 'bin')
ray_out_dir = '/shared/users/asousa/WIPP/3dWIPP/outputs/rays_IGRF_florida'

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



# -------- Partition tasks for MPI --------------------
if rank == 0:
    lats, lons, fs = np.meshgrid(inp_lats, inp_lons, freqs)
    lats = lats.flatten()
    lons = lons.flatten()
    fs   = fs.flatten()

    alts = launch_alt*np.ones_like(lats)
    alts[fs < 600] += 3000e3

    tasklist = zip(alts, lats, lons, fs)

    # tasklist = [(launch_alt, w,x,y) for w,x,y in itertools.product(inp_lats, inp_lons, freqs)]
    chunks = partition(tasklist, nProcs)

else:
    tasklist = None
    chunks = None

tasklist = comm.bcast(tasklist, root=0)
chunks   = comm.bcast(chunks, root=0)

nTasks  = 1.0*len(tasklist)
nSteps = np.ceil(nTasks/nProcs).astype(int)



# print "Subprocess %s on %s thinks tasklist has length %d"%(rank, host, len(tasklist))
# print "Subprocess %s on %s thinks chunks has length %d"%(rank, host, len(chunks))

# ---------- Set up output directory tree -------------
if rank==0:
    if (not os.path.exists(ray_out_dir)):
        os.mkdir(ray_out_dir)

    if (not os.path.exists(os.path.join(ray_out_dir, "logs"))):
        os.mkdir(os.path.join(ray_out_dir,"logs"))

    for f in freqs:
        ray_out_subdir = os.path.join(ray_out_dir,'f_%d'%f)
        if (not os.path.exists(ray_out_subdir)):
            os.mkdir(ray_out_subdir)

    for f in freqs:
        for lo in inp_lons:
            ray_out_subdir = os.path.join(ray_out_dir, 'f_%d'%f, 'lon_%d'%lo)
            if (not os.path.exists(ray_out_subdir)):
                os.mkdir(ray_out_subdir)



# ------------ Load Kp, Dst, etc ---------------------
# Load solar wind parameters (for Tsykadenko corrections)
Pdyn, ByIMF, BzIMF, W = load_TS_params(ray_datenum)

# Load Dst
Dst = load_Dst(ray_datenum)

# Get closest Kp value (Or should we interpolate?)
tvec, kvec = load_Kp()

tt = bisect.bisect_left(tvec, ray_datenum)
Kp = kvec[tt]


# Get closest AE value
tvec, avec = load_ae()
tt = bisect.bisect_left(tvec, ray_datenum)
AE = np.log10(avec[tt])



yearday = '%d%03d'%(ray_datenum.year, ray_datenum.timetuple().tm_yday)
milliseconds_day = (ray_datenum.second + ray_datenum.minute*60 + ray_datenum.hour*60*60)*1e3 + ray_datenum.microsecond*1e-3

# Coordinate transformation library
xf = xflib.xflib(lib_path='/shared/users/asousa/WIPP/3dWIPP/python/libxformd.so')


if rank == 0:
    print "Dst: ", Dst
    print "Kp:  ", Kp
    print "Ae:  ", AE


# Sync up:
time.sleep(30)
comm.Barrier()

working_path = "/tmp";

# Each subprocess does a subset of frequencies
if (rank < len(chunks)):
    # working_path = os.path.join(os.path.expanduser("~"),"rayTmp_%d"%(rank))
    print "Subprocess %s on %s: doing %d rays"%(rank, host, len(chunks[rank]))
    for inp in chunks[rank]:

        # print inp
        ray_out_subdir = os.path.join(ray_out_dir, "f_%d"%inp[3], "lon_%d"%(inp[2]))

        ray_inpfile   = os.path.join(working_path,'ray_inputs_%d_%d_%d.txt'%(inp[3], inp[1], inp[2]))
        ray_outfile   = os.path.join(ray_out_subdir, 'ray_%d_%d_%d.ray'%(inp[3], inp[1], inp[2]))
        ray_tempfile  = os.path.join(working_path,   'ray_%d_%d_%d.ray'%(inp[3], inp[1], inp[2]))
        damp_outfile  = os.path.join(ray_out_subdir,'damp_%d_%d_%d.ray'%(inp[3], inp[1], inp[2]))
        damp_tempfile = os.path.join(working_path,  'damp_%d_%d_%d.ray'%(inp[3], inp[1], inp[2]))

        ray_templog   = os.path.join(working_path,        "ray_%g_%g_%g.log"%( inp[3], inp[1], inp[2])) 
        ray_logfile   = os.path.join(ray_out_dir, "logs", "ray_%g_%g_%g.log"%( inp[3], inp[1], inp[2]))         
        damp_templog  = os.path.join(working_path,        "damp_%g_%g_%g.log"%(inp[3], inp[1], inp[2])) 
        damp_logfile  = os.path.join(ray_out_dir, "logs", "damp_%g_%g_%g.log"%(inp[3], inp[1], inp[2]))

        # print "Cleaning previous runs..."
        if os.path.exists(ray_inpfile):
            os.remove(ray_inpfile)
        if os.path.exists(ray_outfile):
            os.remove(ray_outfile)

        # Rotate from geomagnetic to SM cartesian coordinates
        inp_coords = xf.rllmag2sm(inp, ray_datenum)

        # Write ray to the input file (used by the raytracer):
        f = open(ray_inpfile,'w')
        pos0 = inp_coords
        dir0 = pos0/np.linalg.norm(pos0)    # radial outward
        w0   = inp[3]*2.0*np.pi
        f.write('%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e\n'%(pos0[0], pos0[1], pos0[2], dir0[0], dir0[1], dir0[2], w0))
        f.close()

        # --------- Run raytracer --------

        cmd= '%s/raytracer --outputper=%d --dt0=%g --dtmax=%g'%(ray_bin_dir, 1, dt0, dtmax) + \
             ' --tmax=%g --root=%d --fixedstep=%d --maxerr=%g'%(t_max, root, fixedstep, maxerr) + \
             ' --maxsteps=%d --minalt=%d --inputraysfile=%s --outputfile=%s'%( maxsteps, minalt, ray_inpfile, ray_tempfile) + \
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


        # Run the raytracer
        runlog = subprocess.check_output(cmd, shell=True)


        with open(ray_templog,"w") as file:
            file.write(runlog)
            file.close()

        # ------- Run Damping Code ------------

        damp_cmd =  '%sbin/damping --inp_file %s --out_file %s '%(damping_root, ray_tempfile, damp_tempfile) + \
                    '--Kp %g --AE %g --mode %d'%(Kp, AE, damp_mode) + \
                    ' --yearday %s --msec %d'%(yearday, milliseconds_day)

        damplog = subprocess.check_output(damp_cmd, shell=True)
 
        with open(damp_templog,"w") as file:
            file.write(damplog)
            file.close()


        # Move completed files to the output directory
        if (not os.path.exists(ray_out_subdir)):
            print "Process %d on %s can't find the path at %s"%(rank, host, ray_out_subdir)
        else:
            shutil.move(ray_tempfile,  ray_outfile)     # ray
            shutil.move(damp_tempfile, damp_outfile)    # damp
            shutil.move(ray_templog,   ray_logfile)     # ray log
            shutil.move(damp_templog,  damp_logfile)    # damp log


        os.remove(ray_inpfile)

