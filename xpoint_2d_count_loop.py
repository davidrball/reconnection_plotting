import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches                               
from tristan_funcs import get_val_filename
from tristan_funcs_vecpot import vecpot2
from edge_detector import convolve
from localmin_func import localmin,localmax,localmin_2d, localmin_2d_bdens
import h5py

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

sigma=0.3
va = sigma/np.sqrt(1+sigma)



#fld_base = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam05/bguide.7/output/flds.tot."

fld_base = "../../tristan_acc-mec_Ez/8k_untriggered_bguide.3/output/flds.tot.0"
#fld_base = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/flds.tot.0"

#fld_base = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.7/output/flds.tot."


#if len(tstr)==1:
#    fld_base += "00" + tstr
#elif len(tstr)==2:
#    fld_base += "0"+tstr

t0 = 50
tf = 51

for t in range(t0,tf):
    t_string = str(t)
    myfld = fld_base + t_string

    myfld = h5py.File(myfld,'r')
    dens = get_val_filename('dens',myfld)[0,:,:]
    bdens = myfld['bdens'][0,:,:]

    vecpot = vecpot2(myfld,12,3)
#fig, (ax1, ax2, ax3) = plt.subplots(3,1,sharex=True)
    fig, (ax0) = plt.subplots(1,1,sharex=True)
#ax1.imshow(np.rot90(dens)/5,origin='lower',cmap='viridis',vmin=0,vmax=5)
#ax2.imshow(np.rot90(vecpot),origin='lower',cmap='jet')
#plt.savefig('vecpot_test.png')
#identify midplane 
    ax0.imshow(np.rot90(dens)/4.,origin='lower',vmin=0,vmax=5)
#ax1.imshow(np.rot90(vecpot),origin='lower')
#plt.colorbar(im1)
#plt.colorbar(im0)
#plt.colorbar(orientation='horizontal',label='$A_{z}$')
#plt.savefig('vecpot_fld.png',dpi=300,bbox_inches='tight')
    localmin_x,localmin_y = localmin_2d_bdens(vecpot, 3,bdens,1)

    print('num x points : ', len(localmin_x))
    print(localmin_x)
    print(localmin_y)
    xlen = np.shape(dens)[1]
    #print(xlen)
    for i in range(len(localmin_x)):
        circle = patches.Circle((localmin_y[i],xlen-localmin_x[i]+1),  color="Red",radius=.5)
        ax0.add_patch(circle)
    plt.savefig('plots/2d_xpoint/bguide_untriggered/'+t_string+'.png',dpi=300,bbox_inches='tight')
    myfld.close()
