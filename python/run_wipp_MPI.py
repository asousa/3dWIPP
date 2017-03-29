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
from mpi4py import MPI


project_root = '/shared/users/asousa/WIPP/3dWIPP/'

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
out_Lsh = np.arange(1.2, 8.0, 0.1)
out_lat = np.round(10.0*np.arccos(1.0/out_Lsh)*R2D)/10.0


out_lon = [0]

model_number = 0        # b-field model (0 = dipole, 1 = IGRF)
num_freq_steps = 20     # number of interpolating steps between 
                        # each guide frequency.

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

# -------------- set up MPI -----------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")
nProcs = 1.0*comm.Get_size()


if rank==0:
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    if not os.path.exists(log_directory):
        os.mkdir(log_directory)

    # print "Clearing data from previous runs..."
    # os.system('rm %s/*'%(output_directory))

    # print '-------- Building WIPP ------'
    # os.chdir(project_root)
    # buildstatus = os.system('make')

    # if buildstatus != 0:
    #     sys.exit("Build failed!")
    # os.chdir( os.path.expanduser('~'))
else:
    pass


comm.Barrier()

# --------------- Prep input vector and make chunks -----------


if rank==0:

    tasklist = [(w,x,y) for w,x,y in itertools.product(out_lat, out_lon, freq_pairs)]
    np.random.shuffle(tasklist)
    chunks = partition(tasklist, nProcs)
else:
    tasklist = None
    chunks   = None

comm.Barrier()

tasklist = comm.bcast(tasklist, root=0)
chunks   = comm.bcast(chunks, root=0)

# Split frequency vector into smaller chunks, pass each chunk to a process
nTasks  = 1.0*len(tasklist)
nSteps = np.ceil(nTasks/nProcs).astype(int)


# ------------ Load Kp, Dst, etc ---------------------
# Load solar wind parameters (for Tsykadenko corrections)
# Pdyn, ByIMF, BzIMF, W = load_TS_params(ray_datenum)

# Load Dst
# Dst = load_Dst(ray_datenum)

# # Load Kp
# tvec, kvec = load_Kp()
# tt = bisect.bisect_left(tvec, ray_datenum)
# Kp = kvec[tt]

# # Get closest AE value
# tvec, avec = load_ae()
# tt = bisect.bisect_left(tvec, ray_datenum)
# AE = np.log10(avec[tt])


# Mean parameter vals for set Kp:
Kpvec = [0, 2, 4, 6, 8]
Aevec = [1.6, 2.2, 2.7, 2.9, 3.0]
Dstvec= [-3, -15, -38, -96, -215]
Pdynvec=[1.4, 2.3, 3.4, 5.8, 7.7]
ByIMFvec=[-0.1, -0.1, 0.1, 0.5, -0.2]
BzIMFvec=[1.0, 0.6, -0.5, -2.3, -9.2]


Kp   = Kpvec[vec_ind]
AE   = Aevec[vec_ind]
Pdyn = Pdynvec[vec_ind]
Dst  = Dstvec[vec_ind]
ByIMF= ByIMFvec[vec_ind]
BzIMF= BzIMFvec[vec_ind]
W = np.zeros(6)   # Doesn't matter if we're not using Tsyg

# Tsyganenko input vector
TS_params = [Pdyn, Dst, ByIMF, BzIMF, W[0], W[1], W[2], W[3], W[4], W[5]]


# Run each set of jobs on current node:
if (rank < len(chunks)):
    print "Process %d on host %s, doing %g jobs"%(rank, host, len(chunks[rank]))

    for job in chunks[rank]:
        print job

        olat = job[0]
        olon = job[1]
        flo  = job[2][0]
        fhi  = job[2][1]

        wipp_cmd = '%sbin/wipp --out_dir %s'%(project_root, output_directory) + \
                    ' --iyr %s --idoy %d --isec %d --I0 %d'%(iyr, idoy, isec, flash_I0) + \
                    ' --f_alt %g --f_lat %g --f_lon %g'%(inp_coords[0], inp_coords[1], inp_coords[2]) + \
                    ' --out_lat %g --out_lon %g --b_model %d'%(olat, olon, model_number) + \
                    ' --f1 %g --f2 %g --ray_dir %s'%(flo, fhi, ray_input_directory) +\
                    ' --num_freq_steps %d'%num_freq_steps

        print wipp_cmd

        logfile = os.path.join(log_directory, "wipp_%g_%g_%g.log"%(olat, olon, flo));
        # # os.system(wipp_cmd)
        # file = open(logfile, "w+")
        # subprocess.call(wipp_cmd, shell=True, stdout=file)
        # file.close()

        wipp_log = subprocess.check_output(wipp_cmd, shell=True)
        file = open(logfile,'w')
        file.write(wipp_log)
        file.close()

        # Zip output files to keep everything from getting too huge
        n_filename = os.path.join(output_directory,'pN_%2.1f_%d_%d.dat'%(olat, olon, flo))
        s_filename = os.path.join(output_directory,'pS_%2.1f_%d_%d.dat'%(olat, olon, flo))   

        os.system("gzip %s"%n_filename)
        os.system("gzip %s"%s_filename)

comm.Barrier()

if rank == 0:
    print "-------- Finished WIPP --------"