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
c_omp = 5 #cells per skin depth

sigma=.3
if sigma<1:
    sigma_str = str(sigma)[1:]
else:
    sigma_str = str(sigma)

delgam=.00016
delgam_str = str(delgam)


#flds_base0 = "../../tristan_acc-mec_Ez/testing_edotv/bguide.3_triggered_alwaystrack/output/flds.tot."

#flds_base0 = "../../tristan_testprt/tspec_as_rspectest/bguide.2_triggered/output/flds.tot."

#flds_base0 = "../../tristan_testprt/tristan_nopar/guo_test/output/flds.tot."
flds_base0 = "pleiades_out/interval1000_stride100_prtlrun/flds.tot."



tstart = 100
tfinal = 101

tlist = [2, 3, 4]

tlist = [85,90,95]

for t in range(1):
    t0_str = "%03d" % tlist[0]
    t1_str = "%03d" % tlist[1]
    t2_str = "%03d" % tlist[2]


    f0 = h5py.File(flds_base0+t0_str,'r')
    f1 = h5py.File(flds_base0+t1_str,'r')
    f2 = h5py.File(flds_base0+t2_str,'r')
    #f3 = h5py.File(fld_base3+t_str,'r')


    dens0 = get_val_filename('dens',f0)
    dens1 = get_val_filename('dens',f1)
    dens2 = get_val_filename('dens',f2)
    
    fig, (ax0, ax1, ax2) = plt.subplots(1,3,sharey=True)
    myvmax = 6

    im0 = ax0.imshow((dens0/4.), vmin=1, vmax=myvmax,origin='lower')
    im1 = ax1.imshow((dens1/4.),vmin=1,vmax=myvmax, origin='lower')
    im2 = ax2.imshow((dens2/4.),vmin=1,vmax=myvmax,origin='lower')
    #im3 = ax3.imshow((dens3/4.),vmin=0,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2],origin='lower')
#linear scale

    fig.savefig('testing_out.png',dpi=300,bbox_inches='tight')

    f0.close()
    f1.close()
    f2.close()
    #f3.close()
    plt.close()


