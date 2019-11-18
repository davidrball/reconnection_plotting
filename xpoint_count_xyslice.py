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
fig, (ax0) = plt.subplots(1,1,sharex=True)
#ax1.imshow(np.rot90(dens)/5,origin='lower',cmap='viridis',vmin=0,vmax=5)
#ax2.imshow(np.rot90(vecpot),origin='lower',cmap='jet')
#plt.savefig('vecpot_test.png')
#identify midplane 
ax0.imshow(np.rot90(dens)/4.,origin='lower',vmin=0,vmax=5)
#ax0.imshow(dens/4.,origin='lower',vmin=0,vmax=5)
xmid = np.shape(vecpot)[1]/2
vecpot_x_slice = vecpot[:,xmid]

xlen = np.shape(vecpot)[0]
ylen = np.shape(vecpot)[1]

localmin_y_list = []
localmin_x_list = []
'''
for i in range(15,xlen-15):
    #take slices along y
    vecpot_y_slice = vecpot[i,:]
    tmplist = localmin(vecpot_y_slice,10,10)
    #print(len(tmplist))
    for yind in tmplist:
        localmin_y_list.append(ylen-yind)
        localmin_x_list.append(i)
'''


print('num of xpoints in vertical current layers : ', len(localmin_y_list))
my_localmin_list = localmin(vecpot_x_slice,15,15) #just the 
num_xpoints = len(my_localmin_list)
print('num of xpoints in main current layer' + str(num_xpoints))
for i in range(num_xpoints):
    circle = patches.Circle((my_localmin_list[i],xmid),fill=True,color="Red",radius=1)
    #circle = patches.Circle((xmid,my_localmin_list[i]),fill=True,color="Red",radius=1)
    ax0.add_patch(circle)
    
    #if you want to only search along y at x locations of xpoints in midplane
    
    yval = my_localmin_list[i]
    vecpot_yslice = vecpot[yval,:]
    vert_xpoints = localmin(vecpot_yslice,10,10)
    for j in range(len(vert_xpoints)):
        circle = patches.Circle((my_localmin_list[i],vert_xpoints[j]),fill=True,color="White",radius=1)
        ax0.add_patch(circle)
'''    
for i in range(1,xlen-2):
    vecpot_yslice = vecpot[i,:]
    vert_xpoints = localmin(vecpot_yslice,10,10)
    for j in range(len(vert_xpoints)):                                                
        circle = patches.Circle((i,vert_xpoints[j]),fill=True,color="White",radius=1) 
        ax0.add_patch(circle)  
'''

#print(localmin_x_list)
#print(localmin_y_list)
#for i in range(len(localmin_y_list)):
#    circle = patches.Circle((localmin_x_list[i],localmin_y_list[i]),fill=True,color="White",radius=.25)
    #circle = patches.Circle((localmin_y_list[i],localmin_x_list[i]),fill=True\,color="White",radius=.25)  
    #ax0.add_patch(circle)
plt.savefig('xpoint_count_untriggered.png',dpi=300,bbox_inches='tight')
