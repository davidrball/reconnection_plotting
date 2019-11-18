import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches                               
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tristan_funcs import get_val_filename
from tristan_funcs_vecpot import vecpot2
from edge_detector import convolve
from localmin_func import localmin,localmax,localmin_2d
import h5py

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

sigma=0.3
va = sigma/np.sqrt(1+sigma)



t = 45
tstr = str(t)
#fld_base = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam05/bguide.7/output/flds.tot."

#fld_base = "../../tristan_acc-mec_Ez/8k_untriggered_bguide.3/output/flds.tot.048"
fld_base = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/flds.tot.049"


#fld_base = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.7/output/flds.tot."


#if len(tstr)==1:
#    fld_base += "00" + tstr
#elif len(tstr)==2:
#    fld_base += "0"+tstr

myfld = h5py.File(fld_base,'r')
dens = get_val_filename('dens',myfld)[0,:,:]



vecpot = vecpot2(myfld,12,3)
vecpot_cut = np.rot90(vecpot)[150:-150,150:-150]
fig, axarr = plt.subplots(2,2)
#fig, (ax0) = plt.subplots(1,1,sharex=True)
#ax1.imshow(np.rot90(dens)/5,origin='lower',cmap='viridis',vmin=0,vmax=5)
#ax2.imshow(np.rot90(vecpot),origin='lower',cmap='jet')
#plt.savefig('vecpot_test.png')
#identify midplane 
ax0 = axarr[0,0]
ax1 = axarr[0,1]
ax2 = axarr[1,0]
ax3 = axarr[1,1]


im1 = ax0.imshow(np.rot90(dens)[150:-150,150:-150],origin='lower',vmin=0,vmax=20)

im0 = ax1.imshow(vecpot_cut,origin='lower',vmin=7000,vmax=8750)
ax1.set_title('$A_{z}$')
ax0.set_title("Density")

midx = np.shape(vecpot_cut)[0]/2
midy = np.shape(vecpot_cut)[1]/2 -166 

print('midx : ',midx)
print('midy : ',midy)

vecpot_yslice = vecpot_cut[midx,:] 
vecpot_xslice = vecpot_cut[:,midy]

print('shape vecpot slices : ')
print(np.shape(vecpot_xslice))
print(np.shape(vecpot_yslice))

slicesizex = np.size(vecpot_xslice)
slicesizey = np.size(vecpot_yslice)
my_xarr = np.arange(0,slicesizex)
my_yarr = np.arange(0,slicesizey)

xval_arr = midy*np.ones(slicesizex)
vertarr = np.linspace(6000,12000,slicesizex)

xpoint_val1 = 115
xpoint_val2 = 147
xpoint_val3 = 132


xval_arr1 = xpoint_val1*np.ones(slicesizey)
xval_arr2 = xpoint_val2*np.ones(slicesizey)
vertarr2 = np.linspace(7000,10000,slicesizey)

xval_arr3 = xpoint_val3*np.ones(slicesizey)

ax2.plot(xval_arr1, vertarr2,color="Red",linestyle='dashed')
ax2.plot(xval_arr2,vertarr2,color="Red",linestyle='dashed')
ax2.plot(xval_arr3, vertarr2,color="C1",linestyle='dashed')

ax3.plot(xval_arr, vertarr,color="Black",linestyle='dashed')
ax2.plot(my_xarr,vecpot_xslice)
ax3.plot(my_yarr,vecpot_yslice)
ax2.set_xlabel('y')
ax3.set_xlabel('x')
ax2.set_ylabel('$A_{z}$')
#ax3.set_ylim(5000,12000)
ax2.set_ylim(7000,9000)
ax3.set_yticks([])
#ax2.set_xlim(0,300)
#ax3.set_xlim(0,300)
divider1 = make_axes_locatable(ax1)
divider0 = make_axes_locatable(ax0)
#cax0 = divider0.append_axes("right",size="5%",pad=.05)
#cax1 = divider1.append_axes("right",size="5%",pad=.05)
#plt.colorbar(im1,cax=cax1)
#plt.colorbar(im0,cax=cax0)
#plt.colorbar(orientation='horizontal',label='$A_{z}$')
circle1 = patches.Circle((midy,xpoint_val1),color="Red",radius=2,fill=False)
circle2 = patches.Circle((midy,xpoint_val2),color="Red",radius=2,fill=False)
circle3 = patches.Circle((midy,xpoint_val3),color="C1",radius=2,fill=False)
ax0.add_patch(circle1)
ax0.add_patch(circle2)
ax0.add_patch(circle3)
plt.savefig('vecpot_fld_triggered.png',dpi=300,bbox_inches='tight')

'''
localmin_x,localmin_y = localmin_2d(vecpot, 5)

print('num x points : ', len(localmin_x))
print(localmin_x)
print(localmin_y)

for i in range(len(localmin_x)):
    circle = patches.Circle((localmin_y[i],localmin_x[i]), fill=True, color="Red",radius=1)
    ax0.add_patch(circle)
plt.savefig('2d_xpoint.png',dpi=300,bbox_inches='tight')
'''
