from mpi4py import MPI
import numpy as np
import os
import datetime as dt
import time
from partition import partition
import commands


''' build_gcpm_grid.py
    A port of gcpm_dens_model_buildgrid_parallel.m, originally by Dan Golden and Forrest Foust.
    This version started 8.2016 by Austin Sousa.

    Calculates plasma parameters based on the GCPM model on a regular grid.
    Useful for performing multiple raytracing in the same plasma environment, since GCPM model
    is S L O W.

    This thing is parallelized by divvying up the Z-axis. Don't request any more nodes
    than you have Z entries!

    Key settings:
    L_max, mlt, kp.
'''



# Initialize MPI:
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")


raytracer_root = '/shared/users/asousa/software/foust_raytracer/'
output_dir = '/shared/users/asousa/WIPP/3dWIPP/raytracer_runscripts/foust_output/'
# tmp_dir = '/tmp'

# print "using raytracer at %s"%raytracer_root
R_E = 6370e3;

# Input params
L_max = 10.0
mlt = 1.0
kp = 4.0

# mlt_theta = np.mod(np.pi + (mlt/24.)*2.*np.pi, 2*np.pi);

# xl = sorted(np.append(L_max*R_E*np.cos(mlt_theta + np.array([-1,0,1])*2/24*2*np.pi), 0))
# yl = sorted(np.append(L_max*R_E*np.sin(mlt_theta + np.array([-1,0,1])*2/24*2*np.pi), 0))

minx = -1.0*L_max*R_E;
maxx =  L_max*R_E;
miny = -1.0*L_max*R_E;
maxy =  L_max*R_E;
minz = -1.0*L_max*R_E;
maxz =  L_max*R_E;

# minx = min(xl);
# maxx = max(xl);
# miny = min(yl);
# maxy = max(yl);
# minz = -L_max*R_E/2;
# maxz = L_max*R_E/2;

nx = 80
ny = 80
nz = 80
# nx = int(round((maxx - minx)/R_E)*16)
# ny = int(round((maxy - miny)/R_E)*16)
# nz = int(round((maxz - minz)/R_E)*16)

# print nz
compder = False
day_datevec = dt.datetime(2001, 01, 1, 0, 0, 0);

yearday = '%s%03d'%(day_datevec.year, day_datevec.timetuple().tm_yday)
milliseconds_day = (day_datevec.hour*3600 + day_datevec.minute*60 + day_datevec.second)*1e3 + day_datevec.microsecond*1e-3

    
final_output_filename = os.path.join(output_dir,'gcpm_kp%02.0f_%s_MLD%02d'%(kp*10., day_datevec.strftime("%Y%m%d_%H%M"),mlt))

# Move to running directory
cwd = os.curdir
os.chdir(os.path.join(raytracer_root,'bin'))

# Log the start time
t_start = time.time()



# ------------ Prep stuff on root node: ----------
if rank==0:
    # Clean parts from previous run:
    rmpart_str = 'rm -f ' + os.path.join(output_dir, '*PART*')
    os.system(rmpart_str)

    # Generate z-axis
    z = np.linspace(minz, maxz, nz)

    tasklist = zip(z, np.arange(0, len(z)))

    #np.random.shuffle(tasklist)
    nTasks = 1.0*len(tasklist)
    nProcs = comm.Get_size()
    nSteps = np.ceil(nTasks/nProcs).astype(int)

    chunks = partition(tasklist, nProcs)

else:
    nProcs = None
    chunks = None
    tasklist = None

# # Sync up
# comm.Barrier()
chunks = comm.bcast(chunks,root=0)
nProcs = comm.bcast(nProcs, root=0)
chunk = comm.scatter(chunks)
# print chunk


# # Set up temp directory:

# proc_tmp_path = '%s_%d'%(tmp_dir, rank+1)
# os.system('mkdir %s_%d'%proc_tmp_path)


# ------------ Do stuff on each subnode: -------------
print 'subprocess %d running z(%d:%d) on %s'%(rank, chunk[0][1], chunk[-1][1], host)

this_minz = chunk[ 0][0]
this_maxz = chunk[-1][0]
nz_part = len(chunk)
print " is between %2.2e and %2.2e"%(this_minz, this_maxz)


final_output_filename_part = final_output_filename + '_PART_%02d_of_%02d'%(rank + 1, nProcs)
progress_filename = final_output_filename_part + '_PROGRESS.txt'
# grid_script = (os.path.join(raytracer_root, 'bin/gcpm_dens_model_buildgrid'))
try:
    cmd = ['./gcpm_dens_model_buildgrid --minx=%e --maxx=%e --miny=%e'%(minx, maxx, miny) +
           ' --maxy=%e --minz=%e --maxz=%e' %(maxy, this_minz, this_maxz) +
           ' --nx=%d --ny=%d --nz=%d'%(nx, ny, nz_part) +
           ' --compder=%d --filename=%s --gcpm_kp=%0.1f'%(compder, final_output_filename_part, kp) +
           ' --yearday=%s --milliseconds_day=%d > %s'%(yearday, milliseconds_day, progress_filename)][0]

    print cmd   
    # Run it
    os.system(cmd);

    print "Completed job %d/%d on %s"%(rank +1, nProcs, host)
except:
    print "Something went wrong on %d/%d on %s"%(rank +1, nProcs, host)
# Sync up again
comm.Barrier()

# time.sleep(1)

# Combine work into single file
if rank==0:
        
    print "Merging individual jobs..."
    # Write header
    fid_out = open(final_output_filename + '.txt', 'w');
    nspecies = 4;
    fid_out.write('%10d %10d %10d %10d %10d\n'%(compder, nspecies, nx, ny, nz))
    fid_out.write('%E %E %E %E %E %E\n'%(minx, maxx, miny, maxy, minz, maxz))
    fid_out.close()


    for kk in range(nProcs):
        final_output_filename_part = final_output_filename + '_PART_%02d_of_%02d'%(kk + 1, nProcs)
        
        if kk==0:
            # % If this is the first file, include lines 3 and 4 of the header,
            # % which are fundamental constants
            cat_cmd = 'cat %s | sed -n \'3,$ p\' >> %s'%(final_output_filename_part, final_output_filename + ".txt")
        else:
            cat_cmd = 'cat %s | sed -n \'5,$ p\' >> %s'%(final_output_filename_part, final_output_filename + ".txt")

        os.system(cat_cmd)


    print "Completed!"

else:
    pass
    # print dir(comm)
    # MPI.Finalize()
    # comm.Disconnect()
# Move back home
# os.chdir(cwd)
