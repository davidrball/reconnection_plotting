import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
from matplotlib.colors import PowerNorm
from tristan_funcs import get_val_filename

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

def cut_field(arr, xscan):
    xhlf = np.shape(arr)[1]/2
    xlow = int(xhlf - xscan)
    xup = int(xhlf + xscan)
    return np.rot90(arr[:,xlow:xup])

def identify_central_xpoints(ez_over_vabxy):
    #first, identify middle of simulation
    shape = np.shape(ez_over_vabxy)
    xhlf = int(shape[0]/2)
    strip = ez_over_vabxy[xhlf,:]
    return strip

sigma=0.3
va = sigma/np.sqrt(1+sigma)



t = 50
tstr = str(t)
fld_base = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3_untriggered_stripe60/output/flds.tot.0"

if len(tstr)==1:
    fld_base += "0" + tstr
elif len(tstr)==2:
    fld_base += tstr


myf = h5py.File(fld_base,'r')

dens = get_val_filename('dens',myf)[0,:,:]
ez = get_val_filename('ez',myf)[0,:,:]

bx = get_val_filename('bx',myf)[0,:,:]
by = get_val_filename('by',myf)[0,:,:]
bz = get_val_filename('bz',myf)[0,:,:]

ey = get_val_filename('ey',myf)[0,:,:]
ex = get_val_filename('ex',myf)[0,:,:]


xscan = 150

dens = cut_field(dens, xscan)
ez = cut_field(ez, xscan)
ey = cut_field(ey, xscan)
ex = cut_field(ex, xscan)
bx = cut_field(bx, xscan)
by = cut_field(by, xscan)
bz = cut_field(bz, xscan)

ez_over_vabxy = ez /(va*np.sqrt(bx**2 + by**2))

mystrip = identify_central_xpoints(ez_over_vabxy)




#print(mystrip)

#plt.plot(mystrip)
#plt.ylim(-4,4)


istep = 12
c_omp = 3

yext = (np.shape(dens)[1]*istep / c_omp)/1000 /2
xext = (np.shape(dens)[0]*istep / c_omp)/1000 / 2



xarr = np.linspace(-yext,yext, np.size(mystrip))
#print(xarr)

fig, (ax1, ax2, ax3) = plt.subplots(3,1,sharex=True)

plt.set_cmap('viridis')
im1 = ax1.imshow(dens/4,origin='lower',vmin=0,vmax=5,extent=[-yext,yext,-xext,xext])
plt.set_cmap('RdBu')
im2 = ax2.imshow(ez_over_vabxy,origin='lower',vmin=-1,vmax=1,extent=[-yext,yext,-xext,xext])
#ax1.set_adjustable('box-forced')
#ax2.set_adjustable('box-forced')
'''
divider = make_axes_locatable(ax1)
cax = divider.append_axes('top', size='5%', pad=0.1)
fig.colorbar(im1, cax=cax, orientation = 'horizontal')

divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.1)
fig.colorbar(im2, cax=cax, orientation = 'vertical')
'''
cbar_ax = fig.add_axes([.905,.655,.015,.23])
cb= fig.colorbar(im1,cax=cbar_ax,orientation="vertical")              
#cbar_ax.set_xlabel("Density (particles per cell / $N_{ppc}$)",fontsize=14, labelpad=-4)

cbytick_obj = plt.getp(cb.ax.axes,'xticklabels')
plt.setp(cbytick_obj,fontsize='x-small')

cbar_ax = fig.add_axes([.905,.38,.015,.23])
cb= fig.colorbar(im2,cax=cbar_ax,orientation="vertical",ticks=[-1,-.5,0,.5, 1])
#cbar_ax.set_xlabel("Density (particles per cell / $N_{ppc}$)",fontsize=14, labelpad=-4)    

cbytick_obj = plt.getp(cb.ax.axes,'xticklabels')
plt.setp(cbytick_obj,fontsize='x-small')



#divider = make_axes_locatable(ax3)
#cax = divider.append_axes('right', size='5%', pad=0.1)


ax3.plot(xarr, mystrip)
ax3.set_ylim(-4,4)
ax3.set_ylabel('$E_{z}/v_{A}B_{xy}$')


ax2.set_ylabel('$E_{z}/v_{A}B_{xy}$')

ax1.set_ylabel('Density')



ax3.set_xlim(-yext,yext)


#ax3.set_aspect('auto')

ax3.set_xlabel('x ($1000\;c/\omega_{p}$)',size=14)


#plt.colorbar(im1, cax=ax1)
#plt.colorbar(im2, cax=ax2)


plt.savefig('untriggered_xpoints_stripe60.png',dpi=300,bbox_inches='tight')
plt.close() 

#myshape = np.shape(dens)

#for i in range(myshape[0]):
#    for j in range(myshape[1]):
#        tmpez = ez[i][j]
#        tmpbx = 
