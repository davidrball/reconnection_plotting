import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
import numpy as np
from ideal_fields import return_ideal_fields
from edotv_func import return_edotv


plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 10})

#inputs
 
#mybase = "/home/u21/davidrball/david/guo_comparison/harris/bguide.1/ezoverbxy_run/output/flds.tot."
mybase = "/home/u21/davidrball/david/tristan_acc-mec_Ez/bguide.3_wpar_zcomp/gamprt_thresh50/output/flds.tot."


t0=2
tf=65

edotv_list = []
edotv_par_list = []
edotv_ideal_list = []

t_list = []

edotv_running_tot = 0
edotv_par_running_tot = 0
edotv_ideal_running_tot = 0


for t in range(t0, tf):
    print(t)
    tstr = '%03d'%t
    myf = h5py.File(mybase+tstr,'r')

    fig, (ax0, ax1, ax2, ax3) = plt.subplots(1,4)
    xscan = 150
    dens = myf['dens'][0,:,:]
    xlen = np.shape(dens)[1]
    xhlf = xlen/2
    xlow=xhlf-xscan
    xup=xhlf+xscan
    dens = dens[:,xlow:xup]
    edotve, epardotve, eidealdotve = return_edotv(myf,xscan)
   
    #vmax = np.max(np.abs(edotve))
    #vmin = -vmax

    vmax=.1
    vmin=-.1

    edotvtot = np.sum(edotve)
    epartot = np.sum(epardotve)
    eidealtot = np.sum(eidealdotve)
    
    edotv_running_tot += edotvtot
    edotv_par_running_tot += epartot
    edotv_ideal_running_tot += eidealtot

    edotv_list.append(edotv_running_tot)
    edotv_par_list.append(edotv_par_running_tot)
    edotv_ideal_list.append(edotv_ideal_running_tot)
    t_list.append(t)



    im0 = ax0.imshow(edotve, origin='lower',cmap="RdBu",vmin=vmin,vmax=vmax)
    im1 = ax1.imshow(epardotve,origin='lower',cmap='RdBu',vmin=vmin,vmax=vmax)
    im2 = ax2.imshow(eidealdotve,origin='lower',cmap='RdBu',vmin=vmin,vmax=vmax)
#im3 = ax3.imshow(eperpnonidealdotve,origin='lower',cmap='RdBu',vmin=vmin,vmax=vmax)
    im3 = ax3.imshow(dens,origin='lower',cmap='viridis',vmin=0,vmax=50)


    ax0.set_title('$q\\vec{E} \cdot \\vec{v_{e}}$')
    ax1.set_title('$q\\vec{E_{||}} \cdot \\vec{v_{e}}$')
    ax2.set_title('$q\\vec{E}_{ideal} \cdot \\vec{v_{e}}$')
#ax3.set_title('$q\\vec{E}_{\\bot,nonideal} \cdot \\vec{v_{e}}$')
    ax3.set_title('Density')
    ax1.set_yticks([])
    ax2.set_yticks([])
    ax3.set_yticks([])
#cbar_ax = fig.add_axes([.49,.19,.02,.65])
    cbar_ax2 = fig.add_axes([.98,.19,.02,.6])
#cb = fig.colorbar(im1,cax=cbar_ax)
    cb2 = fig.colorbar(im0,cax=cbar_ax2)
    fig.tight_layout()
    
    plt.savefig('edotv_comp_plots/bguide.3/'+tstr+'.png',dpi=300,bbox_inches='tight')
    plt.close()
    myf.close()

plt.plot(t_list, edotv_list,label='$W_{tot}$')
plt.plot(t_list, edotv_par_list,label='$W_{||}$')
plt.plot(t_list, edotv_ideal_list,label='$W_{ideal}$')
plt.legend(loc='upper left',frameon=False)
plt.ylabel('W')
plt.xlabel('Time')
plt.savefig('edotv_comp_plots/bguide.3/t_vs_W.png')
