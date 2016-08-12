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