import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt




def return_nonideal_fields(f0):
    v3x = f0['v3x'][0,:,:]
    v3y = f0['v3y'][0,:,:]
    v3z = f0['v3z'][0,:,:]

    bx = f0['bx'][0,:,:]
    by = f0['by'][0,:,:]
    bz = f0['bz'][0,:,:]
    ex = f0['ex'][0,:,:]#/bnorm                                                                                                                                                          
    ey = f0['ey'][0,:,:]#/bnorm                                                                                                                                                          
    ez = f0['ez'][0,:,:]#/bnorm                                                                                                                                                         
    e_vecfield = np.dstack((ex,ey,ez))
    b_vecfield = np.dstack((bx,by,bz))
    v_vecfield =np.dstack((v3x,v3y,v3z))
    vcrossb = np.cross(v_vecfield,b_vecfield)
    vcrossb_zcomp = vcrossb[:,:,2]
    vcrossb_ycomp = vcrossb[:,:,1]
    vcrossb_xcomp = vcrossb[:,:,0]
    nonidealz = ez + vcrossb_zcomp
    nonidealy = ey + vcrossb_ycomp
    nonidealx = ex + vcrossb_xcomp

    return nonidealx, nonidealy, nonidealz


def return_ideal_fields(f0):
    v3x = f0['v3x'][0,:,:]
    v3y = f0['v3y'][0,:,:]
    v3z = f0['v3z'][0,:,:]

    bx = f0['bx'][0,:,:]
    by = f0['by'][0,:,:]
    bz = f0['bz'][0,:,:]
    
    b_vecfield = np.dstack((bx,by,bz))
    v_vecfield =np.dstack((v3x,v3y,v3z))
    vcrossb = np.cross(v_vecfield,b_vecfield)
    ideal_zcomp = -vcrossb[:,:,2]
    ideal_ycomp = -vcrossb[:,:,1]
    ideal_xcomp = -vcrossb[:,:,0]
    
    return ideal_xcomp, ideal_ycomp, ideal_zcomp

