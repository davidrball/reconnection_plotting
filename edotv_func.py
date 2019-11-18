import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
import numpy as np
from ideal_fields import return_ideal_fields

def return_edotv(myf,xscan):
    idealx, idealy, idealz = return_ideal_fields(myf)
    dens = myf['dens'][0,:,:]
    myshape = np.shape(dens)
    xlen = myshape[1]
    xhlf = xlen/2
    xlow = int(xhlf-xscan)
    xup = int(xhlf+xscan)
    dens = dens[:,xlow:xup]
    densi = myf['densi'][0,:,xlow:xup]
    dense = dens - densi 

    idealx = idealx[:,xlow:xup]
    idealy = idealy[:,xlow:xup]
    idealz = idealz[:,xlow:xup]
    ex = myf['ex'][0,:,xlow:xup]
    ey = myf['ey'][0,:,xlow:xup]
    ez = myf['ez'][0,:,xlow:xup]
    bx = myf['bx'][0,:,xlow:xup]
    by = myf['by'][0,:,xlow:xup]
    bz = myf['bz'][0,:,xlow:xup]
    u = myf['v3x'][0,:,xlow:xup]
    v = myf['v3y'][0,:,xlow:xup]
    w = myf['v3z'][0,:,xlow:xup]
    ui = myf['v3xi'][0,:,xlow:xup]
    vi = myf['v3yi'][0,:,xlow:xup]
    wi = myf['v3zi'][0,:,xlow:xup]
    mi=1836*densi
    me=1*dense
#here me is representing the mass density
    ue = ((mi+me)*u-mi*ui)/me
    ve = ((mi+me)*v-mi*vi)/me
    we = ((mi+me)*w-mi*wi)/me
    edotve = -(ex*ue + ey*ve + ez*we)
    bmag = np.sqrt(bx**2 + by**2 + bz**2)
    unitbx = bx/bmag
    unitby = by/bmag
    unitbz = bz/bmag
    emag = np.sqrt(ex**2 + ey**2 + ez**2)
    epar = ex*unitbx + ey*unitby + ez*unitbz
    eparx = epar*unitbx
    epary = epar*unitby
    eparz = epar*unitbz
    epardotve = -(eparx*ue + epary*ve + eparz*we)
    eidealdotve = -(idealx*ue + idealy*ve + idealz*we)
    eperpnonidealdotve = (edotve - epardotve) - eidealdotve

    return edotve, epardotve, eidealdotve


def return_epardotv_comps(myf,xscan):

    dens = myf['dens'][0,:,:]
    myshape = np.shape(dens)
    xlen = myshape[1]
    xhlf = xlen/2
    xlow = int(xhlf-xscan)
    xup = int(xhlf+xscan)
    dens = dens[:,xlow:xup]
    densi = myf['densi'][0,:,xlow:xup]
    dense = dens - densi 
    ex = myf['ex'][0,:,xlow:xup]
    ey = myf['ey'][0,:,xlow:xup]
    ez = myf['ez'][0,:,xlow:xup]
    bx = myf['bx'][0,:,xlow:xup]
    by = myf['by'][0,:,xlow:xup]
    bz = myf['bz'][0,:,xlow:xup]
    u = myf['v3x'][0,:,xlow:xup]
    v = myf['v3y'][0,:,xlow:xup]
    w = myf['v3z'][0,:,xlow:xup]
    ui = myf['v3xi'][0,:,xlow:xup]
    vi = myf['v3yi'][0,:,xlow:xup]
    wi = myf['v3zi'][0,:,xlow:xup]
    mi=1836*densi
    me=1*dense
#here me is representing the mass density                                                   
    ue = ((mi+me)*u-mi*ui)/me
    ve = ((mi+me)*v-mi*vi)/me
    we = ((mi+me)*w-mi*wi)/me

    bmag = np.sqrt(bx**2 + by**2 + bz**2)
    unitbx = bx/bmag
    unitby = by/bmag
    unitbz = bz/bmag
    emag = np.sqrt(ex**2 + ey**2 + ez**2)
    epar = ex*unitbx + ey*unitby + ez*unitbz
    eparx = epar*unitbx
    epary = epar*unitby
    eparz = epar*unitbz
    
    eparx *= -ue
    epary *= -ve
    eparz *= -we
    return eparx, epary, eparz
