import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt

myt = 242
myt_str = '%03d'%myt
myt*=500*.15/100 #in units of 100 inverse omega p
myf = "../../tristan_acc-mec_Ez/bguide.3_wpar_zcomp/output/flds.tot."+myt_str
myf = h5py.File(myf,'r')


print('plot time:',myt)
dens = myf['dens'][0,:,:]

plt.imshow(dens,origin='lower')
plt.savefig('testingtime.png')
