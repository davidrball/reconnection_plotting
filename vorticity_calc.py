import numpy as np



def return_vorticity(vx, vy): #return the z component of vorticity
    #vx and vy in code units, so x is across, y is up

    dvydx = np.roll(vy,1,axis=1)-np.roll(vy,-1,axis=1)
    dvxdy = np.roll(vx,1,axis=0)-np.roll(vx,-1,axis=0)
    
    vortz = dvydx - dvxdy
    return vortz
    #return dvxdy
    #return dvydx
