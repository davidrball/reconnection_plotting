# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 11:14:44 2017

@author: davidrball
"""

# -*- coding: utf-8 -*-a
"""
Created on Thu Jul 28 11:21:48 2016

@author: davidrball
"""
import numpy as np
import h5py 
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 10})

def get_val_filename(input_string, fopen):
    str_list = list(fopen.keys())
    my_string = input_string
    if my_string in str_list:
        
        my_list = list(fopen.get(my_string))
        return my_list
    else:
        print("Value not stored, maybe you had a typo")
        print('your options are : ',list(fopen.keys()))
        
                
                
   
def temp_ions(fopen, xmid, xscan):
    #input params
    ppc0 = 4        
    keni = np.array(get_val_filename('keni',fopen))[0,:,:]

    #should make sure this xscan is appropriate across all timesteps
    xhlf = keni.shape[1]/2
    xhlf = xmid #put xhlf in by hand
    xmin = int(xhlf - xscan)
    xmax = int(xhlf + xscan)
    keni = keni[:,xmin:xmax]
    pdens_array = np.array(get_val_filename('pdens',fopen))[0,:,xmin:xmax]
    dens_array = np.array(get_val_filename('dens', fopen))[0,:,xmin:xmax]    
    bdens_array = np.array(get_val_filename('bdens',fopen))[0,:,xmin:xmax]
    
    
    #need to compute the bulk velocity and its associated gamma
    v3x = np.array(get_val_filename('v3xi',fopen))[0,:,xmin:xmax]
    v3y = np.array(get_val_filename('v3yi',fopen))[0,:,xmin:xmax]
    v3z = np.array(get_val_filename('v3zi', fopen))[0,:,xmin:xmax]
    #gamma array given by:
    vtot2 = v3x**2 + v3y**2 + v3z**2    
    gam = np.sqrt(1/(1-vtot2)) #lorentz factor associated with bulk velocity
    #print(np.max(gam))
    #should we identify recon region first?  Or just do it everywhere?
    #let's identify recon region
  
    
    Ti_array = np.zeros(np.shape(dens_array))
    for i in range(np.shape(dens_array)[0]):
        for j in range(np.shape(dens_array)[1]):
            ad_index = 4./3
                #do whatever you want here
            varx = dens_array[i][j] / (gam[i][j] * .5 * ppc0)
            varz = (keni[i][j]-gam[i][j]) * dens_array[i][j] / (.5 * ppc0)
                #pretty much just copying Lorenzo's script here
            maxabs = 1
            while maxabs >= 1e-2:
                vary = varz/(ad_index*gam[i][j]**2 - (ad_index-1))
                delgam1 = vary/varx * (ad_index-1) #dimensionless temp
                zeta1 = delgam1/(.24 + delgam1)
                ad_index1 = 1./3 * (5 - 1.121937*zeta1 + .18203*zeta1**2 - .96583*zeta1**3 + 2.32513*zeta1**4 - 2.39332*zeta1**5 + 1.07136*zeta1**6)
                    #computing the adiabatic index from temp
                maxabs = (ad_index-ad_index1) / ad_index
                    # this is the fractional change in the ad index
                ad_index = ad_index1
                    #print(delgam1, ad_index)
            Ti_array[i][j] = delgam1
    return Ti_array

