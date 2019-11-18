import numpy as np

import math                  
import numpy as np 
import matplotlib.pyplot as plt 
import h5py
from scipy import signal
from scipy.special import kn

def get_val_filename(input_string, fopen):
    str_list = list(fopen.keys())
    my_string = input_string
    if my_string in str_list:
        
        my_list = np.array(fopen.get(my_string))
        return my_list
    else:
        print("Value not stored, maybe you had a typo")


def marg_x(my_array):
    #my array is 2D, we want to marginalize along the "x" direction in the code

    xlen = np.shape(my_array)[1]
    ylen = np.shape(my_array)[0]
    oneD_dens = np.zeros(ylen)
    for i in range(xlen):
        oneD_dens += my_array[:,i]/xlen
    return oneD_dens


def decompose_spectra(myspece, myrspece, mygam):
    mylen = np.size(mygam)

    returnrspece = []
    #print(np.shape(myspece), np.shape(myrspece), np.shape(mygam))                           
    for i in range(mylen):

        if myspece[i] > 1.2*myrspece[i]:
            pass
        elif myspece[i] != 0 and mygam[i]*myspece[i] > 10:

        #elif myspece[i] > np.max(myspece)/100:                                              
            print(myspece[i], np.max(myspece))
            return np.array(mygam[i+4:]), np.array(myrspece[i+4:])
            break


def vecpot(myf,xwide0,istep):
    by = get_val_filename('by',myf)[0,:,:]
    bx = get_val_filename('bz',myf)[0,:,:]

    istep = 12

    dx = istep
    dy = istep

    #xmid = by.shape[1]/2
    #choose xmind based on by profile

    bysum = np.sum(by,axis=0)
    xmid = np.argmin(np.abs(bysum))
    xmin = int(xmid - xwide0)
    xmax = int(xmid + xwide0)
    
    bx = bx[:,xmin:xmax]
    by = by[:,xmin:xmax]
    myshape = bx.shape

    #create 2nd array to stack on top of eachother
    bxtot = np.zeros([myshape[0]*2,myshape[1]])
    bytot = np.zeros([myshape[0]*2,myshape[1]])
    for i in range(myshape[0]):
        bxtot[i,:] += bx[i,:]
        bxtot[i+myshape[0],:] += bx[i,:]
        bytot[i,:] += by[i,:]
        bytot[i+myshape[0],:] += by[i,:]
   
    sz = bxtot.shape
    mx00 = sz[0]
    my00 = sz[1]
    print('mx00, my00,' , mx00, my00)
    A0xtr = np.fft.fft2(bxtot)
    A0ytr = np.fft.fft2(bytot)

    Wxtr = np.zeros(A0xtr.shape)
    Wytr = np.zeros(A0ytr.shape)

    vec1 = np.zeros(A0xtr.shape)
    vec2 = np.zeros(A0ytr.shape)
    for i in range(0, int(mx00)):
  
        vec1[i,:] += (np.cos(2*np.pi * i / mx00) - 1.)/dx**2
    for j in range(0,my00):
        vec2[:,j] += (np.cos(2*np.pi*j/my00)-1.)/dy**2

    Wxtr = A0xtr / (2.*(vec1+vec2))
    Wytr = A0ytr / (2.*(vec1+vec2))

    Wxtr[0,0] = 0.
    Wytr[0,0] = 0.

    Wx = np.fft.ifft2(Wxtr)
    Wy = np.fft.ifft2(Wytr)

    mx00 = int(mx00/2)

    Wy1up = np.roll(Wy[0:mx00,:], -1, 1) #shifts 1 to the left
    Wy1dn = np.roll(Wy[0:mx00,:], 1, 1) #shifts 1 to the right
    Wx1up = np.roll(Wx[0:mx00,:], -1, 0)
    Wx1dn = np.roll(Wx[0:mx00,:], 1, 0)

    A1z = -(.5*(Wy1up - Wy1dn)/dx - .5*(Wx1up-Wx1dn)/dy)
    A1z[0,:] = 0
    A1z[:,0] = 0
    A1z[mx00-1,:] = 0
    A1z[:,my00-1]=0

    realvecpot = np.real(A1z)
    return realvecpot

def prtl_spec(gam_array, bin_num):
    gamma_array = gam_array - 1


    gamlist = list(gamma_array)

    gammin = np.min(gamma_array)
    gammax = np.max(gamma_array)
    gammin_log = np.log10(gammin)
    gammax_log = np.log10(gammax)

    hist_bins = np.logspace(gammin_log, gammax_log, bin_num)
    hist_array = np.zeros(bin_num)

    for i in range(bin_num-1):
        hist_low = hist_bins[i]
        hist_up = hist_bins[i+1]
        for loggam in gamlist:
            if hist_low < loggam < hist_up:
                hist_array[i] += 1
    return hist_bins, hist_array

def prtl_spec_fixedgam(gam_array, gammin, gammax, bin_num):
    gamma_array = gam_array - 1


    gamlist = list(gamma_array)

    #gammin = np.min(gamma_array)
    #gammax = np.max(gamma_array)
    gammin_log = np.log10(gammin)
    gammax_log = np.log10(gammax)

    hist_bins = np.logspace(gammin_log, gammax_log, bin_num)
    hist_array = np.zeros(bin_num)

    for i in range(bin_num-1):
        hist_low = hist_bins[i]
        hist_up = hist_bins[i+1]
        for loggam in gamlist:
            if hist_low < loggam < hist_up:
                hist_array[i] += 1
    return hist_bins, hist_array

def return_powerlaw(norm,p,x_array):
    y_array = norm * x_array**(-p+1)
    #print(x_array)
    #print(y_array)
    #myax.plot(x_array,y_array, color = col,linestyle="--", linewidth=2.0, label="Power law\n p="+str(p))
    return x_array,y_array

def vecpot2(flds,istep,comp):
   

    bx_arr = flds["bx"]
    by_arr = flds["by"]
   
 
    bx = bx_arr[0]
    by = by_arr[0]
   
 
    #Indexing convention opposite from IDL.
    mx = by.shape[1]
    my = by.shape[0]
 
    #print('mx, my',mx,my)
    xwide0 = 200000
    ybstep = 0
 
    #istep = 12.0 
    #comp = 3.0
 
    dx = 1.0*istep
    dy = 1.0*istep
 
    #Choose xmid, define center of box based off by.
    varmid = np.sum(by,axis=0)
    xmid = np.where(abs(varmid) == min(abs(varmid)))
    xmid = xmid[0][0]
 
    # Select region (in x) for computation.
    xwide = math.floor(xwide0*comp/istep) 
    xmin = max((xmid-xwide),2)
    xmax = min((xmid+xwide),mx-4)
 
    dlft = xmid-xmin
    drgt = xmax-xmid
    dchoice = min(dlft,drgt)
    xmin = xmid-dchoice
    xmax = xmid+dchoice
 
    #put xmin in by hand
    xmin = 5

    #print('xmin, xmax',xmin,xmax)
    #print('xlen : ',xmax-xmin)
    # Copy bx and by, removing z-dimension and creating periodic boundaries for
    # the Fourier transform in order to calculate potential.
    #print('ybstep : ',ybstep)
    #print('bx shape', bx.shape)
    varx = bx[ybstep:my-1-ybstep,xmin:xmax]
    #print('varx shape',varx.shape)
    A0x = np.concatenate((varx,np.fliplr(varx)),axis=1) 
 
    vary = by[ybstep:my-1-ybstep,xmin:xmax]
    A0y = np.concatenate((vary,np.fliplr(vary)),axis=1) 
 
    sz = A0x.shape
    mx00 = sz[1] 
    my00 = sz[0] 
    #print('mx00, my00', mx00, my00)
    #print('sz : ',sz)

    A0xtr = np.fft.fft2(A0x) #Fourier Transform
    A0ytr = np.fft.fft2(A0y) 
 
    Wxtr = A0xtr*0.0
    Wytr = A0ytr*0.0
 
    vec1 = A0xtr*0.0
    for i in range(0,mx00): 
        vec1[:,i] = (math.cos(float(2)*math.pi*i/(mx00*1.0))-1)/(dx**2)
 
    vec2 = A0ytr*0.0
    for j in range(0,my00):
        vec2[j,:] = (math.cos(float(2)*math.pi*j/(my00*1.0))-1)/(dy**2)
 
    Wxtr = np.divide(A0xtr,(2.0*(vec1+vec2))) 
    Wytr = np.divide(A0ytr,(2.0*(vec1+vec2)))
 
    Wxtr[0,0] = 0.0
    Wytr[0,0] = 0.0
 
    Wx = np.fft.ifft2(Wxtr) #Inverse Fourier Transform
    Wy = np.fft.ifft2(Wytr)
 
    mx00 = mx00/2
 
    A1z = np.empty((my00,mx00),dtype=float)
    Wx = np.array(Wx[:,0:mx00],dtype=float)
    Wy = np.array(Wy[:,0:mx00],dtype=float)
 
    Wy1up = np.roll(Wy,-1,axis=1) 
    Wy1dn = np.roll(Wy,1,axis=1)
    Wx1up = np.roll(Wx,-1,axis=0)
    Wx1dn = np.roll(Wx,1,axis=0)
    A1z = -(.5*(Wy1up-Wy1dn)/dx - .5*(Wx1up-Wx1dn)/dy)
 
    #Construct vector potential array
    A1z[:,0] = 0.0
    A1z[0,:] = 0.0
    A1z[:,mx00-1] = 0.0
    A1z[my00-1,:] = 0.0
    A1zraw = A1z
    A1z = np.empty((my,A1zraw[0].shape[0]),dtype=float)
    A1z[ybstep:my-ybstep-1,:] = A1zraw
    
    return np.real(A1z),xmin,xmax

def return_maxwell(theta, coeff,gamma_array):
    #coeff = 2.2e4                                                              
    gamma_array += 1
    beta_array = np.sqrt(1-1/gamma_array**2)
    dist_array = coeff*(gamma_array**2 * beta_array / (theta*kn(2,1/theta)))*np.exp(-gamma_array/theta)

    dist_array = coeff*(gamma_array**2 * np.float128(np.float128(beta_array))/np.float128((theta*kn(2,1/theta))))*np.exp(-gamma_array/theta)

    return dist_array


