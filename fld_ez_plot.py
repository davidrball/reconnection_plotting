import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
from matplotlib.colors import PowerNorm

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

plt.set_cmap("viridis")


def get_val_filename(input_string,fopen):
    str_list = list(fopen.keys())
    my_string = input_string
    if my_string in str_list:
        my_list = np.array((fopen.get(my_string))[0])
        return my_list
    else:
        print("argument not found, your options are:")
        print(list(fopen.keys()))
        return 0

istep = 12 #downsampling                                                            
c_omp = 3 #cells per skin depth                                                     
t = 30
t_str = str(t)

#fld_base0 = "../../tristan_acc-mec_Ez/8k_run2/output/flds.tot.0"
#f0 = h5py.File(fld_base0 + t_str,'r')

#fld_base = "../../tristan_acc-mec_Ez/8k_untriggered_bguide.3/output/flds.tot.050"

#fld_base = "../../tristan-mp_reconnection/guide_fields/sig.3/mi1/lowbeta/bguide.3/output/flds.tot.014"

fld_base = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/flds.tot.046"

f0 = h5py.File(fld_base,'r')

dens = get_val_filename('dens',f0)
ez = get_val_filename('ez',f0)
bx = get_val_filename('bx',f0)
by = get_val_filename('by',f0)

sig = 0.3
va = np.sqrt(sig/ (1+sig))

ez_bxy = ez / (va*(bx**2 + by**2))


xscan = 50

x0len = np.shape(dens)[1]
x0hlf = x0len/2
x0low = int(x0hlf - xscan)
x0up = int(x0hlf+xscan)

dens = np.rot90(dens[:,x0low:x0up])
ez_bxy = np.rot90(ez_bxy[:,x0low:x0up])

xlen = 2*xscan
ylen = np.shape(dens)[1]

xext = xlen * istep / c_omp
yext = ylen * istep / c_omp

fig, (ax0, ax1) = plt.subplots(2,1,sharex=True)

im0 = ax0.imshow(dens/4., vmin=0, vmax=5, extent = [-yext/2, yext/2, -xext/2, xext/2],origin='lower')

im1 = ax1.imshow(ez_bxy, vmin=-2, vmax=2,extent = [-yext/2, yext/2, -xext/2, xext/2],origin='lower',cmap = 'RdBu')

cbar_ax = fig.add_axes([.125, .93, .77, .03])
cb = fig.colorbar(im0, cax=cbar_ax,orientation="horizontal")
cbar_ax.set_xlabel("Density (particles per cell / $N_{ppc}$",fontsize=14,labelpad=-45)

cbar_ax2 = fig.add_axes([.125, .4, .77, .015])
cb2 = fig.colorbar(im1, cax=cbar_ax2,orientation="horizontal")
cbar_ax2.set_xlabel("$E_{z}/v_{A}B_{xy}$",fontsize=12,labelpad=-36)


fig.savefig('guidefield_merge_ez.png',dpi=300,bbox_inches='tight')
