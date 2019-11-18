import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})
plt.set_cmap("viridis")

time_string = "035"

bg0 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide0/output/flds.tot."+time_string
time_string = "035"

myf = h5py.File(bg0, 'r')





bx = myf['bx'][0,:,:]
by = myf['by'][0,:,:]



shape = np.shape(bx)
xsize = shape[1]
ysize = shape[0]

print(xsize)
xhlf = xsize / 2
xscan=100
xlow = int(xhlf-xscan)
xup = int(xhlf+xscan)
print(xlow, xup)
bx=bx[:,xlow:xup]
by=by[:,xlow:xup]

shape = np.shape(bx)
xsize = shape[0]
ysize = shape[1]



xarr = np.arange(0,xsize)
yarr = np.arange(0,ysize)

print(xsize, ysize)
plt.streamplot(yarr, xarr, by, bx)
plt.savefig('testing_streamplot.png')
