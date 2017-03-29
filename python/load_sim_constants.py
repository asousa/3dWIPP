# Load simulation constants from "consts.h" into Python.
# use:
# sc = load_sim_constants(path/to/const.h)


import numpy as np

# class consts_struct(object):
#     # pass
#     def __init__(self):
#         self.params = ['Q_EL','M_EL','E_EL','EPS0','C','Z0','R_E',
#                   'H_MAGNETO', 'H_IONO', 'P_A','P_B','H_E',
#                   'E_MIN', 'E_MAX', 'NUM_E','EALimS','EALimN','NUM_EA',
#                   'FREQ_STEP_SIZE','LAT_STEP_SIZE','LON_STEP_SIZE',
#                   'MAX_GROUND_DISTANCE','TIME_MAX','NUM_TIMES',
#                   'DAMPING_THRESH']

#         for p in self.params:
#             exec 'self.%s = None'%p


def load_sim_constants(directory):
    params = ['Q_EL','M_EL','E_EL','EPS0','C','Z0','R_E',
              'H_MAGNETO', 'H_IONO', 'P_A','P_B','H_E',
              'E_MIN', 'E_MAX', 'NUM_E','EALimS','EALimN','NUM_EA',
              'FREQ_STEP_SIZE','LAT_STEP_SIZE','LON_STEP_SIZE',
              'MAX_GROUND_DISTANCE','TIME_MAX','NUM_TIMES']


    sc = dict()
    with open(directory,'r+') as file:
        for line in iter(file):
            tmp = line.split()
            if len(tmp) >= 2:
                if tmp[0] == '#define':
                    if tmp[1] in params:
                        # print tmp
                        try:
                            sc[tmp[1]] = float(tmp[2])
                            # tmpstr = 'sc.%s = %s'%(tmp[1], tmp[2])
                            # exec tmpstr
                            # print "succeeded: ", tmpstr
                        except:
                            print "failed:", tmp

    # # Do the terms we know it won't pick up:

        
        sc['MU0'] = np.pi*4e-7

        sc['TIME_STEP'] = (1.0*((1.0*sc['TIME_MAX'])/sc['NUM_TIMES']))
        sc['E_EXP_BOT'] = np.log10(sc['E_MIN'])
        sc['E_EXP_TOP'] = np.log10(sc['E_MAX'])
        sc['DE_EXP'] = ( (sc['E_EXP_TOP'] - sc['E_EXP_BOT'])/sc['NUM_E'])

        # Generate energy and velocity arrays
        sc['E_tot_arr'] = pow(10,sc['E_EXP_BOT'] + sc['DE_EXP']*np.arange(0,sc['NUM_E']))
        sc['v_tot_arr'] = sc['C']*np.sqrt(1 - pow(sc['E_EL']/(sc['E_EL'] + sc['E_tot_arr']),2))
    
        # Cast some of the parameters to ints
        sc['NUM_TIMES'] = int(sc['NUM_TIMES'])
        sc['NUM_E'] = int(sc['NUM_E'])


    return sc
