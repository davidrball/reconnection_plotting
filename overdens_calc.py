import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
from tristan_funcs import recon_region
testfld = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.7/output/flds.tot.040"


#we want to calculate sheet thicknesses based on bdens

def return_avg_dens(file_string):
    myf = h5py.File(file_string)
    dens = np.array(myf['dens'][0,:,:])
    bdens = np.array(myf['bdens'][0,:,:])
    pdens = np.array(myf['pdens'][0,:,:])
    cury = recon_region(dens,pdens,bdens)
    dens -= bdens
    myshape = np.shape(cury)
    xlen = myshape[0]
    ylen = myshape[1]
    testarr = np.zeros((xlen,ylen))
    dens_list = []
    for i in range(xlen):
        for j in range(ylen):
            if cury[i][j]==1:      
               dens_list.append(dens[i][j])
    meandens = np.mean(np.array(dens_list))
    stddens = np.std(np.array(dens_list))
    #plt.imshow(testarr, origin='lower')
    #plt.savefig('testimg_dens.png')
    print('mean density in recon region : ',meandens)
    return meandens, stddens


                
#return_avg_dens(testfld)
