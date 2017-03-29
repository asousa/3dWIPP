# Load_phi_files.py
# 11.12.2015
# Austin Sousa
#
#   Loads all Phi files from a directory
#   Returns N, S, L
#   N and S are numpy ndarrays (num_T, num_E, num_L).
#   Data are in units of [counts / (cm^2 keV s)]
#   L is a vector of L-shells found in the directory.
#
#   V1.0. Seems to match the Matlab verson!
#   V1.1  Started 6.1.2016, adapting to work with constants object
#   v1.2  12.2016, added support for gzipped files (since it adds up fast)

import os
# import glob
import re
#import struct
import numpy as np
import gzip
#from matplotlib import pyplot as plt

def load_phi_files_L(dataDir, sc):
    # This version for files with L-shell in the filename:

    # allfiles = os.listdir(os.getcwd() + '/' + dataDir)
    allfiles = os.listdir(dataDir)

    n_files = sorted([f for f in allfiles if 'phi_' in f and '_N' in f])
    s_files = sorted([f for f in allfiles if 'phi_' in f and '_S' in f])

    if not n_files and not s_files:
        print "No phi files found!"



    L_n = sorted([float(f[4:(len(f) - 2)]) for f in n_files])
    L_s = sorted([float(f[4:(len(f) - 2)]) for f in s_files])
    #L_n = sorted([float(re.findall("\d+\.\d+",f)[0]) for f in n_files])
    #L_s = sorted([float(re.findall("\d+\.\d+",f)[0]) for f in s_files])
    
    assert (L_n == L_s), "North / South mismatch!"
    

    L_vec = np.array(L_n)
    N_arr = np.zeros([sc.NUM_E, sc.NUM_STEPS, len(L_vec)])
    S_arr = np.zeros([sc.NUM_E, sc.NUM_STEPS, len(L_vec)])

    #for ind, filename in enumerate(n_files):
    for L_ind, L in enumerate(L_vec):
        # print "L: ",L
        
        # Binary files -- little-endian, 4-byte floats 
        dt = np.dtype('<f4')
        phi_N = np.fromfile(os.path.join(dataDir,"phi_%g_N"%(L)),dtype=dt)
        phi_S = np.fromfile(os.path.join(dataDir,"phi_%g_S"%(L)),dtype=dt)

        N_arr[:,:,L_ind] = phi_N.reshape(sc.NUM_E, sc.NUM_STEPS, order='c')
        S_arr[:,:,L_ind] = phi_S.reshape(sc.NUM_E, sc.NUM_STEPS, order='c')

    print "N NaNs: %d" % np.sum(np.isnan(N_arr))
    print "S NaNs: %d" % np.sum(np.isnan(S_arr))
    return N_arr, S_arr, L_vec

def load_phi_files_latlon(dataDir, lon, sc, zipped=False):
    # This version for files with output latitude and longitude in the filename:
    
    # allfiles = os.listdir(os.getcwd() + '/' + dataDir)
    allfiles = os.listdir(dataDir)
    if zipped:
        n_files = sorted([f for f in allfiles if f.startswith('phi_') and f.endswith('_N.dat.gz')])
        s_files = sorted([f for f in allfiles if f.startswith('phi_') and f.endswith('_S.dat.gz')])
    else:
        n_files = sorted([f for f in allfiles if f.startswith('phi_') and f.endswith('_N.dat')])
        s_files = sorted([f for f in allfiles if f.startswith('phi_') and f.endswith('_S.dat')])

    if not n_files and not s_files:
        print "No phi files found!"

    # p = re.compile("\d+")
    
    # lat_n = sorted([float(p.findall(f)[0]) for f in n_files])
    # lat_s = sorted([float(p.findall(f)[0]) for f in s_files])
    # lon_n = sorted([float(p.findall(f)[1]) for f in n_files])
    # lon_s = sorted([float(p.findall(f)[1]) for f in s_files])

    lat_n = sorted([float(f.split('_')[1]) for f in n_files])
    lat_s = sorted([float(f.split('_')[1]) for f in s_files])
    lon_n = sorted([float(f.split('_')[2]) for f in n_files])
    lon_s = sorted([float(f.split('_')[2]) for f in s_files])
    print lat_n
    print lat_s
    # assert (lat_n == lat_s), "North / South mismatch!"
    # assert (lon_n == lon_s), "North / South mismatch!"

    

    latvec = np.unique(np.array(lat_s))
    lonvec = np.unique(np.array(lon_s))

    if lon in lonvec:
        N_arr = np.zeros([sc['NUM_E'], sc['NUM_TIMES'], len(latvec)])
        S_arr = np.zeros([sc['NUM_E'], sc['NUM_TIMES'], len(latvec)])

        #for ind, filename in enumerate(n_files):
        for lat_ind, lat in enumerate(latvec):
            
            # Binary files -- little-endian, 4-byte floats 
            dt = np.dtype('<f8')
            if zipped:
                gzf = gzip.open(os.path.join(dataDir,"phi_%g_%g_N.dat.gz"%(lat, lon)),'rb').read()
                phi_N = np.fromstring(gzf, dtype=dt)
                gzf = gzip.open(os.path.join(dataDir,"phi_%g_%g_S.dat.gz"%(lat, lon)),'rb').read()
                phi_S = np.fromstring(gzf, dtype=dt)
            else:
                phi_N = np.fromfile(os.path.join(dataDir,"phi_%g_%g_N.dat"%(lat, lon)),dtype=dt)
                phi_S = np.fromfile(os.path.join(dataDir,"phi_%g_%g_S.dat"%(lat, lon)),dtype=dt)

            N_arr[:,:,lat_ind] = phi_N.reshape(sc['NUM_E'], sc['NUM_TIMES'], order='c')
            S_arr[:,:,lat_ind] = phi_S.reshape(sc['NUM_E'], sc['NUM_TIMES'], order='c')

        print "N NaNs: %d" % np.sum(np.isnan(N_arr))
        print "S NaNs: %d" % np.sum(np.isnan(S_arr))
        return N_arr, S_arr, latvec


if __name__ == '__main__':
    # ----------------------------------
    # Main Program
    # ----------------------------------
    N, S, L = load_phi_files("outputs/run_GIOEEV/in_45", num_E = 128, num_T = 600)
