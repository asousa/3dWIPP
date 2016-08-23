# load_rayfile.py
# Loads a ray file from newray1.8, in python
# Started 3/2/2016, Austin Sousa

import numpy as np
import pandas as pd
import pickle
#from build_database import flux_obj
from scipy import interpolate
import matplotlib.pyplot as plt
import os
import itertools
import random

# def load_rayfile(directory, frequency):

#     return load_rayfile(directory, filename=['newray', frequency,'.dat'], damping_name = ['d',frequency,'.dat'])

def load_rayfile(directory, frequency):
    ''' 
    Loads a rayfile and damping file, as returned by the Fortran raytracer and C Landau damping calculation.
    Inputs: either frequency, or {filename, damping_name}
    '''

    filename = 'newray%g.dat'%frequency
    damping_name = 'd%g.dat'%frequency


    full_filename = os.path.join(directory,filename)
    print "loading ", full_filename

    # Load into a Pandas dataframe. Magic!
    df = pd.read_csv(full_filename,header=None,delim_whitespace=True)

    # Rename columns (keeping consistency with previous scripts, sorry for shouting)
    # This version as used in the raytracer
    #df.columns=['TG','DRE','LAT','DELTA','TP','ELE','PSIG','PSIRAY','PSIRG','AN','ANE','ANH','ANHE','ANO','GF','AP','AR','AL']    
    # This version as used in Jacob's code
    #df.columns=['tg','distre','lat','delta','tp','l_sh','psi','psiray','psires','mu','dens','anH','anHe','anO','fH','stixP','stixR','stixL']
    df.columns=['tg','distre','lat','delta','tp','l_sh','psi','psiray','psires','mu','dens','anH','anHe','anO','fH','stixP','stixR','stixL']    

    # Delete the first two bogus rows
    #df.drop(df.index[0],inplace=True)
    
    # Find start and end indices of rays (a little obnoxious -- rays are stacked and marked with a "99999" in the TG column to separate)
    i9_inds = df[df['tg']==99999].index
    k1 = i9_inds[:-1] + 1
    k2 = i9_inds[1:]

    # Load damping file:
    damp = pd.read_csv(os.path.join(directory,damping_name),header=None,delim_whitespace=True)

    di9_inds = damp[damp[0]==99999].index

    df['power'] = damp.shift(i9_inds[0] + 1)
    df.scaled = False # Have we scaled the power yet, or is it just the damping vector?


    # Tidy up all the rays into a dictionary of dataframes:
    ray_dict = {}
    for x in xrange(len(k1)):
        #print int(df.iloc[k1[x]+1].LAT)
        ray_dict[(df.iloc[k1[x]].lat)] = df.iloc[k1[x]:k2[x]]
        ray_dict[(df.iloc[k1[x]].lat)].launch_lat  = df.iloc[k1[x]].lat
        ray_dict[(df.iloc[k1[x]].lat)].frequency  = frequency
        
    #print ray_dict.keys()
    return ray_dict
    #return df, k1, k2





if __name__ =='__main__':
    #df, k1, k2 = load_rayfile(directory='/shared/users/asousa/WIPP/WIPPy/python/',filename='newray200.dat',damping_name='d200.dat')
    df, k1, k2 = load_rayfile(directory='/shared/users/asousa/WIPP/WIPPy/python/',frequency = 200)
    
    print df