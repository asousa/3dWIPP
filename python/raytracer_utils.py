import numpy as np
import pandas as pd
import os
# from spacepy import coordinates as coord
# from spacepy.time import Ticktock



def read_rayfiles(directory, freq, latmin, latmax, lonmin, lonmax):
    '''read rayfiles from a directory '''

    out = []
    for root, dirs, files in os.walk(os.path.join(directory, "f_%d"%freq)):
        for f in files:
            if f.startswith('ray') and f.endswith('.ray'):
                row = (f.split(".")[0]).split("_")
                cfreq = float(row[1])
                clat  = float(row[2])
                clon  = float(row[3])
                
                if ( (cfreq == freq) and (clat >= latmin) and (clat <= latmax) and
                (clon >= lonmin) and (clon <= lonmax) ):
                    tmp = read_rayfile(os.path.join(root, f))

                    # check damping:
                    dpath = os.path.join(root, 'damp_%d_%d_%d.ray'%(cfreq, clat, clon))

                    dtmp = read_damp(dpath)

                    for ind in range(len(tmp)):
                        tmp[ind]['damping'] = dtmp[ind]['damping']
                    out.extend(tmp)
                    # out.extend(read_rayfile(os.path.join(root,f)))

    return out


def read_damp(dampfile):
    x = pd.read_csv(dampfile, delim_whitespace=True, header=None)

    raynums = np.unique(x[0])
    numrays = len(raynums)

    out = []

    for ii in range(numrays):
        tmp = x[x[0]==raynums[ii]]
        tmp.reset_index(drop=True)
        data = dict()
        data['time'] = tmp[1]
        data['damping'] = tmp[2]

        out.append(data)

    return out



def read_rayfile(rayfile):
    ''' Load output from Forest's raytracer'''
    x = pd.read_csv(rayfile,delim_whitespace=True, header=None)
    # % fields: raynum, stopcond, time, pos(3), vprel(3), vgrel(3), n(3),
    # % B0(3), Nspec, qs(Nspec), ms(Nspec), Ns(Nspec), nus(Nspec)

    raynums = np.unique(x[0])
    numrays = len(raynums)

    out = []
    for ii in range(numrays):
        tmp = x[x[0]==raynums[ii]]
        tmp.reset_index(drop=True)
        #data = pd.DataFrame(columns=['raynum','time','pos','vprel','vgrel','n','B0','qs','ms','Ns','nus'])
        data = dict()
        data['time'] = tmp[2]
        data['pos'] = tmp.loc[:,3:5]
        data['pos'].columns=['x','y','z']
        data['vprel'] = tmp.loc[:,6:8]
        data['vprel'].columns=['x','y','z']
        data['vgrel'] = tmp.loc[:,9:11]
        data['vgrel'].columns=['x','y','z']
        data['n'] = tmp.loc[:,12:14]
        data['n'].columns=['x','y','z']
        data['B0'] = tmp.loc[:,15:17]
        data['B0'].columns=['x','y','z']
        
        # #data = pd.DataFrame(columns=['raynum','time','pos','vprel','vgrel','n','B0','qs','ms','Ns','nus'])
        # data = dict()
        # data['time'] = tmp[2]

        # data['pos']   = coord.Coords(zip(tmp.loc[:,3], tmp.loc[:,4],  tmp.loc[:,5]), 'SM','car')

        # data['pos']   = data['pos'].
        # data['vprel'] = tmp.loc[:,6:8]
        # data['vprel'].columns=['x','y','z']
        # data['vgrel'] = tmp.loc[:,9:11]
        # data['vgrel'].columns=['x','y','z']
        # data['n'] = tmp.loc[:,12:14]
        # data['n'].columns=['x','y','z']
        # data['B0'] = tmp.loc[:,15:17]
        # data['B0'].columns=['x','y','z']
        
        



        

        data['w'] = tmp.iloc[0,18]       # Frequency
        data['Nspec'] = tmp.iloc[0,19]   # Number of species
        data['stopcond'] = tmp.iloc[0,1] # Stop condition
        data['qs'] =  tmp.loc[:,20 + 0*data['Nspec']:20 + 1*data['Nspec'] - 1]
        data['ms'] =  tmp.loc[:,20 + 1*data['Nspec']:20 + 2*data['Nspec'] - 1]
        data['Ns'] =  tmp.loc[:,20 + 2*data['Nspec']:20 + 3*data['Nspec'] - 1]
        data['nus']=  tmp.loc[:,20 + 3*data['Nspec']:20 + 4*data['Nspec'] - 1]

        # read damping data too, if we have it:
        if (tmp.shape[1] > 20+4*data['Nspec']):
            data['damping'] = tmp.loc[:,20+4*data['Nspec']]

            
        out.append(data)
    return out

def readdump(filename):
    ''' Reads the dump files from Forests raytracer
        % x - vector of x coordinates in meters
        % y - vector of y coordinates in meters
        % z - vector of z coordinates in meters
        % 
        % qs - array(ispecies, ix, iy, iz) of charges for each species
        % Ns - array(ispecies, ix, iy, iz) of number densities in m^-3 
        % Ms - array(ispecies, ix, iy, iz) of masses in kg
        % nus - array(ispecies, ix, iy, iz) of collision frequencies in s^-1
        % B0 - array(icomponent, ix,iy, iz) containing the background magnetic
    '''
    out = dict()
    
    f = open(filename,'r')
    line = f.readline().split()
    nspec = int(line[0])
    nx = int(line[1])
    ny = int(line[2])
    nz = int(line[3])
    
    line = f.readline().split()
    minx = float(line[0])
    maxx = float(line[1])
    miny = float(line[2])
    maxy = float(line[3])
    minz = float(line[4])
    maxz = float(line[5])
    
#     print maxz
    
    dat = [float(x) for x in f.readlines()]
#     print np.shape(dat)
#     print dat[0:10]
    dat = np.reshape(dat,[nspec*4 + 3, nx, ny, nz],order='f')
#     print np.shape(dat)
    
    out['qs'] = dat[0*nspec:(1*nspec),:,:,:]
    out['Ns'] = dat[1*nspec:(2*nspec),:,:,:]
    out['Ms'] = dat[2*nspec:(3*nspec),:,:,:]
    out['nus']= dat[3*nspec:(4*nspec),:,:,:]
    out['B0'] = dat[4*nspec:-1,:,:,:]
    
    out['x'] = np.linspace(minx, maxx, nx)
    out['y'] = np.linspace(miny, maxy, ny)
    out['z'] = np.linspace(minz, maxz, nz)
    
    return out