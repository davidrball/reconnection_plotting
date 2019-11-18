import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.patches as patches


plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 12})
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
           
mydens = np.rot90(h5py.File(myname,'r')['dens'][0,:,:])

bx = -np.rot90(h5py.File(myname,'r')['bx'][0,:,:])
by = np.rot90(h5py.File(myname,'r')['by'][0,:,:])

myshape = np.shape(mydens)
xlen = myshape[1]
ylen = myshape[0]
yhlf = ylen/2.
ylow = int(yhlf - xscan)
yup = int(yhlf + xscan)

mydens = mydens[ylow:yup,:]


bx = bx[ylow:yup,:]
by = by[ylow:yup,:]

xext = xlen*istep/c_omp/1000
        #title = str(betaeff_list[myax_ind])
myim = ax0.imshow(mydens/4., vmin=0,vmax=5,extent = [-xext/2, xext/2, -yext, yext],origin='lower')


myshape = np.shape(bx)
ylen = myshape[0]
xlen = myshape[1]


#xarr = np.linspace(0,xlen*100,xlen)
#yarr = np.linspace(0,ylen*100,ylen)
xarr = np.linspace(-xext/2,xext/2,xlen)
yarr = np.linspace(-yext,yext,ylen)

ax0.streamplot(xarr, yarr,by,bx,color="White",linewidth=0.5,arrowsize=.5, density=[2,1])
#ax0.quiver(xarr, yarr, by, bx)

rectangle1 = patches.Rectangle((-.7,-.07),.1,.15,color="Red",fill=False,zorder=2)

rectangle2 = patches.Rectangle((.1,-.07),.15,.1,color="Cyan",fill=False,zorder=2)
ax0.add_patch(rectangle1)
ax0.add_patch(rectangle2)

ax0.scatter(.18,-.01,marker='x',color="Red",s=10)
ax0.scatter(-.66, 0,marker='x',color="Red",s=10)
#ax0.set_adjustable('box-forced')

#plt.colorbar()
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
plt.savefig('blines_withx.pdf',dpi=300,bbox_inches='tight')
plt.savefig('blines_withx.png',dpi=300,bbox_inches='tight')
