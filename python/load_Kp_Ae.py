import scipy.io
import numpy as np
import os.path 
from scipy.interpolate import interp1d
import datetime



def load_Kp(file = '/shared/users/asousa/WIPP/3dWIPP/data/Kp_1999_2016.dat'):

    # Kp files are ascii with a per-character format that doesn't split up easily
    if os.path.exists(file):

        tvec = [];
        kvec = [];
        avec = [];
        
        with open(file) as f:
            for line in f:
                yr = int(line[0:2])
                mo = int(line[2:4])
                dy = int(line[4:6])
                
                if 17 <= yr <= 99:
                    yr += 1900
                else:
                    yr += 2000
                    
    #             print yr, mo, dy
                # Each row contains entries for 3 hour increments.
                # Timestamp is in the middle of each bin.
                
                tvec.append(dt.datetime(yr,mo,dy,1,30,0))
                tvec.append(dt.datetime(yr,mo,dy,4,30,0))
                tvec.append(dt.datetime(yr,mo,dy,7,30,0))
                tvec.append(dt.datetime(yr,mo,dy,10,30,0))
                tvec.append(dt.datetime(yr,mo,dy,13,30,0))
                tvec.append(dt.datetime(yr,mo,dy,16,30,0))
                tvec.append(dt.datetime(yr,mo,dy,19,30,0))
                tvec.append(dt.datetime(yr,mo,dy,22,30,0))
                
                kvec.append(float(line[12:14])/10.)
                kvec.append(float(line[14:16])/10.)
                kvec.append(float(line[16:18])/10.)
                kvec.append(float(line[18:20])/10.)
                kvec.append(float(line[20:22])/10.)
                kvec.append(float(line[22:24])/10.)
                kvec.append(float(line[24:26])/10.)
                kvec.append(float(line[26:28])/10.)
                
    #             kvec.append(float(line[12:14]))
    #             kvec.append(float(line[14:16]))
    #             kvec.append(float(line[16:18]))
    #             kvec.append(float(line[18:20]))
    #             kvec.append(float(line[20:22]))
    #             kvec.append(float(line[22:24]))
    #             kvec.append(float(line[24:26]))
    #             kvec.append(float(line[26:28]))
                
                avec.append(float(line[31:34]))
                avec.append(float(line[34:37]))
                avec.append(float(line[37:40]))
                avec.append(float(line[40:43]))
                avec.append(float(line[43:46]))
                avec.append(float(line[46:49]))
                avec.append(float(line[49:52]))
                avec.append(float(line[52:55]))
                
    return tvec, kvec


def load_ae(file = '/shared/users/asousa/WIPP/3dWIPP/data/DST_AE.dat', daily_avg = False):
    # Load AE index:

    # Kp files are ascii with a per-character format that doesn't split up easily

    tvec = []
    avec = []

    if os.path.exists(file):
        with open(file) as f:
            for line in f:
                if line[0:2] == 'AE':
    #                 print line
                    yr = int(line[14:16] + line[3:5])
                    mo = int(line[5:7])
                    dy = int(line[8:10])
                
                    if daily_avg:
                        avec.append(int(line[116:120]))
                        tvec.append(dt.datetime(yr,mo,dy,12,0,0))
                    else:
                        for x in range(0,24):
                        
                            ae = int(line[20 + 4*x:20 + 4*(x+1)])
                            avec.append(ae)
                            tvec.append(dt.datetime(yr,mo,dy,x,0,0))

    return tvec, avec
