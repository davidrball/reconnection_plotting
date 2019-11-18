import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 8})
plt.set_cmap("viridis")

from timespec_func import plot_timespec


#define the simulation matrix, right now it's 3x3
#first index picks out guide field strength, 2nd is temp




#guide field 0, delgam=.00005, .0005, .005
#myname = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_run2/output/flds.tot."

myname = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/flds.tot."

fig, ax0  = plt.subplots(1,1)
t = 49
tstr = "%03d" % t

myname += tstr

print(myname)

xscan = 150
istep = 12
c_omp = 3


yext = xscan*istep/c_omp / 1000
           
mydens = h5py.File(myname,'r')['dens'][0,:,:]

bx = h5py.File(myname,'r')['bx'][0,:,:]
by = h5py.File(myname,'r')['by'][0,:,:]
myshape = np.shape(mydens)
xlen = myshape[1]
ylen = myshape[0]
yhlf = xlen/2.
ylow = int(yhlf - xscan)
yup = int(yhlf + xscan)

mydens = mydens[:,ylow:yup]


bx = bx[:,ylow:yup]
by = by[:,ylow:yup]

xext = xlen*istep/c_omp/1000
        #title = str(betaeff_list[myax_ind])
myim = ax0.imshow(mydens/4., vmin=0,vmax=5,extent = [-yext, yext, -xext/2, xext/2.],origin='lower')


myshape = np.shape(bx)
ylen = myshape[0]
xlen = myshape[1]


#xarr = np.linspace(0,xlen*100,xlen)
#yarr = np.linspace(0,ylen*100,ylen)
yarr = np.linspace(-xext/2,xext/2,ylen)
xarr = np.linspace(-yext,yext,xlen)

ax0.streamplot(xarr, yarr,bx,by,color="White",linewidth=0.5,arrowsize=.5, density=[1,2])
#ax0.quiver(xarr, yarr, by, bx)

#ax0.set_adjustable('box-forced')

#cbar_ax = fig.add_axes([.905,.145,.015,.7])
#cb= fig.colorbar(myim,cax=cbar_ax,orientation="vertical")
#cbar_ax.set_ylabel("Density ($N_{ppc}$ / $N_{ppc0}$)",fontsize=10)#, labelpad=-10)

#cbytick_obj = plt.getp(cb.ax.axes,'xticklabels')
#plt.setp(cbytick_obj,fontsize='x-small')

#cbytick_obj = plt.getp(cb.ax.axes,'xticklabels')
#plt.setp(cbytick_obj,fontsize='x-small')


xstr = "$x \; (1000 \; c/\omega_{p})$"
ystr = "$y \; (1000 \; c/\omega_{p})$"
ax0.set_ylabel(ystr)                                        
ax0.set_xlabel(xstr)
plt.savefig('testing_blines.png',dpi=300,bbox_inches='tight')
