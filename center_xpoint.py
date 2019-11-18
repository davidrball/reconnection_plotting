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

def cut_field(arr, xscan, yscan):
    xhlf = np.shape(arr)[1]/2
    yhlf = np.shape(arr)[0]/2                                                  
    xlow = int(xhlf - xscan)
    xup = int(xhlf + xscan)
    ylow = int(yhlf - yscan)
    yup = int(yhlf + yscan)
    return np.rot90(arr[ylow:yup,xlow:xup])





t = 30
tstr = str(t)
fld_base = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide2/output/flds.tot.0"

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



#we want to make a small square around the center of the box

xhlf = np.shape(dens)[1]/2
yhlf = np.shape(dens)[0]/2


#dimensions of square
xscan = 100
yscan = 100

xlow = int(xhlf - xscan)
xup = int(xhlf + xscan)
ylow = int(yhlf - yscan)
yup = int(yhlf + yscan)


dens = cut_field(dens, xscan,yscan)
ez = cut_field(ez, xscan,yscan)
ex = cut_field(ex, xscan, yscan)
ey = cut_field(ey,xscan,yscan)

bx = cut_field(bx, xscan,yscan)
by = cut_field(by, xscan,yscan)
bz = cut_field(bz, xscan,yscan)

EdotB = np.abs(ex*bx + ey*by + ez*bz)
ezbz = np.abs(ez*bz)

ez_bxy = np.abs(ez) / np.sqrt(bx**2 + by**2)


sig = .3
va = np.sqrt(sig / (sig+1))
ez_bxy /= va


print(np.shape(dens))
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,sharex='col',sharey='row')
im1 = ax1.imshow(dens/4.,origin='lower',vmin=0,vmax=5)
im2 = ax2.imshow(ez_bxy, origin='lower',cmap="Blues",vmin=0,vmax=2)
im3 = ax3.imshow(EdotB, origin='lower',cmap="Blues",vmin=0,vmax=.5)
im4 = ax4.imshow(bz, origin='lower',vmin=0,vmax=20)
print('maximum ez/bxy: ',np.max(ez_bxy))
plt.colorbar(im1, ax=ax1)
plt.colorbar(im2, ax=ax2)
plt.colorbar(im3, ax=ax3)
plt.colorbar(im4,ax=ax4)



ax1.set_title('Density')
ax2.set_title('$|E_{z}|/\\beta_{A}B_{xy}$')
ax3.set_title('$|E \cdot B|$')
ax4.set_title('$B_{z}$')

ax1.set_adjustable('box-forced')
ax2.set_adjustable('box-forced')
ax3.set_adjustable('box-forced')
ax4.set_adjustable('box-forced')



plt.savefig('bguide2_center.png',bbox_inches='tight', dpi=300)
plt.close()
myf.close()



