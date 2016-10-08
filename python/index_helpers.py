import scipy.io
import numpy as np
import os.path 
from scipy.interpolate import interp1d
import datetime


def load_TS_params(t, datapath = '/shared/users/asousa/software/vlf_matlabwork/indices/'):
    ''' Load values from the TS04 database files, and interpolate parameters'''

    filename = os.path.join(datapath,'%04d_OMNI_5m_with_TS05_variables.mat'%t.year)

    if os.path.exists(filename):
        dd = scipy.io.loadmat(filename)
        TS04_date = dd['TS04_date']
        Pdyn_vec  = dd['Pdyn']
        W_vec = dd['W']
        
        # Requested time as a datenum
        intime_dn = datetime2matlabdn(t)

        Pdyn   = interp1d(dd['TS04_date'].squeeze(), dd['Pdyn'].squeeze() ).__call__(intime_dn)
        ByIMF  = interp1d(dd['TS04_date'].squeeze(), dd['ByIMF'].squeeze()).__call__(intime_dn)
        BzIMF  = interp1d(dd['TS04_date'].squeeze(), dd['BzIMF'].squeeze()).__call__(intime_dn)
        TS04_W = np.array([interp1d(dd['TS04_date'][:,0], dd['W'][:,x]).__call__(intime_dn) for x in range(np.shape(dd['W'])[1])])

        return Pdyn, ByIMF, BzIMF, TS04_W

    else:
        print 'Can''t load TS04 file! Using defaults'
        #      Pdyn, ByIMF, BzIMF, W[0:5]
        return 4,    0,    -5,    [0.132,    0.303,    0.083,    0.070,    0.211,    0.308 ]


def load_Dst(t, file = '/shared/users/asousa/software/vlf_matlabwork/indices/dst.mat'):

    if os.path.exists(file):
        dd = scipy.io.loadmat(file)
        intime_dn = datetime2matlabdn(t)
        
        return interp1d(dd['dst_date'].squeeze(), dd['dst'].squeeze()).__call__(intime_dn)
    else:
        print 'can''t find Dst file! Returning 0'
        return 0

def datetime2matlabdn(dt):
    '''Convert a Python datetime object to a Matlab datenum'''
    mdn = dt + datetime.timedelta(days = 366)
    frac_seconds = (dt-datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
    frac_microseconds = dt.microsecond / (24.0 * 60.0 * 60.0 * 1000000.0)
    return mdn.toordinal() + frac_seconds + frac_microseconds

def matlabdn2datetime(datenum):
    '''Convert a matlab datenum to a Python datetime object'''
    return datetime.datetime.fromordinal(int(datenum)) + datetime.timedelta(days=datenum%1) - datetime.timedelta(days=366)



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
                
                tvec.append(datetime.datetime(yr,mo,dy,1,30,0))
                tvec.append(datetime.datetime(yr,mo,dy,4,30,0))
                tvec.append(datetime.datetime(yr,mo,dy,7,30,0))
                tvec.append(datetime.datetime(yr,mo,dy,10,30,0))
                tvec.append(datetime.datetime(yr,mo,dy,13,30,0))
                tvec.append(datetime.datetime(yr,mo,dy,16,30,0))
                tvec.append(datetime.datetime(yr,mo,dy,19,30,0))
                tvec.append(datetime.datetime(yr,mo,dy,22,30,0))
                
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
                        tvec.append(datetime.datetime(yr,mo,dy,12,0,0))
                    else:
                        for x in range(0,24):
                        
                            ae = int(line[20 + 4*x:20 + 4*(x+1)])
                            avec.append(ae)
                            tvec.append(datetime.datetime(yr,mo,dy,x,0,0))

    return tvec, avec