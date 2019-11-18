import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from tristan_funcs import get_val_filename
from tristan_funcs_vecpot import vecpot2
from edge_detector import convolve
from localmin_func import localmin,localmax
import h5py

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

sigma=0.3
va = sigma/np.sqrt(1+sigma)



t = 25
tstr = str(t)

#fld_base = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/flds.tot.040"
fld_base = "../../tristan_acc-mec_Ez/8k_untriggered_bguide.3/output/flds.tot.050"
#fld_base = "../../tristan_acc-mec_Ez/8k_bguide0_untriggered_stride1_thresh2/output/flds.tot.024"

myfld = h5py.File(fld_base,'r')
dens = get_val_filename('dens',myfld)[0,:,:]
vecpot = vecpot2(myfld,12,3)
#fig, (ax1, ax2, ax3) = plt.subplots(3,1,sharex=True)
#fig, (ax0) = plt.subplots(1,1,sharex=True)
#ax1.imshow(np.rot90(dens)/5,origin='lower',cmap='viridis',vmin=0,vmax=5)
#ax2.imshow(np.rot90(vecpot),origin='lower',cmap='jet')
#plt.savefig('vecpot_test.png')
#identify midplane 
#ax0.imshow(np.rot90(dens)/4.,origin='lower',vmin=0,vmax=5)
#ax0.imshow(np.rot90(dens/4.),origin='lower',vmin=0,vmax=5)
xmid = np.shape(vecpot)[1]/2
vecpot_x_slice = vecpot[:,xmid]

xlen = np.shape(vecpot)[0]

my_localmin_list = localmin(vecpot_x_slice,15,15) #just the 
num_xpoints = len(my_localmin_list)
print('num of xpoints in main current layer' + str(num_xpoints))
localmin_list = []
for i in range(num_xpoints):
    yval = my_localmin_list[i]
    vecpot_yslice = vecpot[yval,:]
    vert_xpoints = localmin(vecpot_yslice,10,10)
    for j in range(len(vert_xpoints)):
        myx = vert_xpoints[j]
        myy = np.max(vecpot_yslice)
        plt.scatter(myx,myy,color="Red")
    

    #print(min(vecpot_yslice),max(vecpot_yslice))
    plt.plot(vecpot_yslice)
#plt.ylim(6000,8000)
plt.xlabel('$y$')
plt.ylabel('$A_{z}$')
plt.savefig('testing_vecpot_yslice_untriggered.png',dpi=300,bbox_inches='tight')

