import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
import numpy as np
from ideal_fields import return_ideal_fields
from edotv_func import return_edotv, return_epardotv_comps


plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 10})

#inputs
 
#mybase = "/home/u21/davidrball/david/guo_comparison/harris/bguide.1/ezoverbxy_run/output/flds.tot."
mybase = "/home/u21/davidrball/david/tristan_acc-mec_Ez/bguide.3_wpar_zcomp/gamprt_thresh50/output/flds.tot."


t0=2
tf=50

edotv_list = []
edotv_par_list = []
edotv_ideal_list = []

epardotvx_list = []
epardotvy_list = []
epardotvz_list = []


t_list = []

edotv_running_tot = 0
edotv_par_running_tot = 0
edotv_ideal_running_tot = 0

epardotvx_running_tot = 0
epardotvy_running_tot = 0
epardotvz_running_tot = 0


for t in range(t0, tf):
    print(t)
    tstr = '%03d'%t
    myf = h5py.File(mybase+tstr,'r')

    fig, (ax0, ax1, ax2, ax3) = plt.subplots(1,4)
    xscan = 50
    dens = myf['dens'][0,:,:]
    xlen = np.shape(dens)[1]
    xhlf = xlen/2
    xlow=xhlf-xscan
    xup=xhlf+xscan
    dens = dens[:,xlow:xup]
    edotve, epardotve, eidealdotve = return_edotv(myf,xscan)
   
    epardotvex, epardotvey, epardotvez = return_epardotv_comps(myf, xscan)

    vmax = np.max([np.max(np.abs(epardotvex)), np.max(np.abs(epardotvey)), np.max(np.abs(epardotvez))])
    vmin = -vmax

    edotvtot = np.sum(edotve)
    epartot = np.sum(epardotve)
    eidealtot = np.sum(eidealdotve)
    
    epardotvex_tot = np.sum(epardotvex)
    epardotvey_tot = np.sum(epardotvey)
    epardotvez_tot = np.sum(epardotvez)

    epardotvx_running_tot += epardotvex_tot
    epardotvy_running_tot += epardotvey_tot
    epardotvz_running_tot += epardotvez_tot

    epardotvx_list.append(epardotvx_running_tot)
    epardotvy_list.append(epardotvy_running_tot)
    epardotvz_list.append(epardotvz_running_tot)


    edotv_running_tot += edotvtot
    edotv_par_running_tot += epartot
    edotv_ideal_running_tot += eidealtot

    edotv_list.append(edotv_running_tot)
    edotv_par_list.append(edotv_par_running_tot)
    edotv_ideal_list.append(edotv_ideal_running_tot)
    t_list.append(t)



    im0 = ax0.imshow(epardotve, origin='lower',cmap="RdBu",vmin=vmin,vmax=vmax)
    im1 = ax1.imshow(epardotvex,origin='lower',cmap='RdBu',vmin=vmin,vmax=vmax)
    im2 = ax2.imshow(epardotvey,origin='lower',cmap='RdBu',vmin=vmin,vmax=vmax)
#im3 = ax3.imshow(eperpnonidealdotve,origin='lower',cmap='RdBu',vmin=vmin,vmax=vmax)
    im3 = ax3.imshow(epardotvez,origin='lower',cmap='RdBu',vmin=vmin,vmax=vmax)


    ax0.set_title('$q\\vec{E}_{||} \cdot \\vec{v_{e}}$')
    ax1.set_title('$qE_{||,x}v_{x}$')
    ax2.set_title('$qE_{||,y}v_{y}$')
#ax3.set_title('$q\\vec{E}_{\\bot,nonideal} \cdot \\vec{v_{e}}$')
    ax3.set_title('$qE_{||,z}v_{z}$')
    ax1.set_yticks([])
    ax2.set_yticks([])
    ax3.set_yticks([])
#cbar_ax = fig.add_axes([.49,.19,.02,.65])
    cbar_ax2 = fig.add_axes([.98,.19,.02,.6])
#cb = fig.colorbar(im1,cax=cbar_ax)
    cb2 = fig.colorbar(im0,cax=cbar_ax2)
    fig.tight_layout()
    
    plt.savefig('edotv_comp_plots/bguide.3/epar_comps/'+tstr+'.png',dpi=300,bbox_inches='tight')
    plt.close()
    myf.close()

plt.plot(t_list, epardotvx_list,label='$W_{||,x}$')
plt.plot(t_list, epardotvy_list,label='$W_{||,y}$')
plt.plot(t_list, epardotvz_list,label='$W_{||,z}$')
plt.legend(loc='upper left',frameon=False)
plt.ylabel('W')
plt.xlabel('Time')
plt.savefig('edotv_comp_plots/bguide.3/t_vs_Wpar_comps.png')
