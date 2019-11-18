import matplotlib
matplotlib.use('Agg')

import numpy as np
import h5py
from temp_calc import temp_ions
import matplotlib.pyplot as plt

def beta_calc(T, bx, by, bz, densi,beta0, T0,ppc0):
    ni0 = ppc0/2.
    #return array of betas in same shape as T/b arrays
    #need array of ion temps, and b field strengths

    #beta init is just initialized beta
    #T init is just initial T

    #first find the upstream bfield values

    B0 = by[0][0] #just take grid cell far from center
    bx /= B0
    by /= B0
    bz /= B0

    bx_B02 = bx**2 + by**2 + bz**2
    
    beta = beta0*(1/bx_B02)*(T/T0)*(densi/ni0)
    
    return beta
    

    #
def beta_slice(file_string,beta0,T0,ppc0,xscan):
    myf = h5py.File(file_string,'r')
    bx = myf['bx'][0,:,:]
    xhlf = np.shape(bx)[1]/2
    new_xscan = xscan
    xscan=xhlf
    xlow = int(xhlf-xhlf)
    xup = int(xhlf+xhlf)


    bx = bx[:,xlow:xup]
    by = myf['by'][0,:,xlow:xup]
    bz = myf['bz'][0,:,xlow:xup]

    keni = myf['keni'][0,:,xlow:xup]
    densi = myf['densi'][0,:,xlow:xup]
    Ti = 2*temp_ions(myf,xhlf,xscan)
    #just return a slice
    beta = beta_calc(Ti, bx, by, bz, densi, beta0, T0, ppc0)[:,xhlf-new_xscan:xhlf+new_xscan]
    myf.close()
    #return np.mean(beta,axis=1)
    return beta
#test_str =  "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/flds.tot.020"
'''
test_str = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.7/output/flds.tot.020"
beta0 = 0.0033
T0 = .0005
ppc0 = 4

mybeta = beta_slice(test_str,beta0,T0,ppc0)

meanbeta = np.mean(mybeta,axis=1)
#meanTi = np.mean(Ti,axis=1)
plt.plot(np.log10(meanbeta))
#plt.ylim(-1,3)
plt.savefig('beta_slice_test.png')
'''
