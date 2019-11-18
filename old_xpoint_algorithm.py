import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
                               
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



t = 47
tstr = str(t)
fld_base = "../../tristan_acc-mec_Ez/8k_untriggered_bguide.3/output/flds.tot."

#fld_base = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam05/bguide.7/output/flds.tot."

#fld_base = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam00005/output/flds.tot."

#fld_base = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.7/output/flds.tot."


if len(tstr)==1:
    fld_base += "00" + tstr
elif len(tstr)==2:
    fld_base += "0"+tstr

myfld = h5py.File(fld_base,'r')
dens = get_val_filename('dens',myfld)[0,:,:]
vecpot = vecpot2(myfld,12,3)


myshape = np.shape(dens)

convfac = (12/3.)/1000.
xext = myshape[1]*12 / 3. / 1000.
yext = myshape[0]*12 / 3. / 1000.


#fig, (ax1, ax2, ax3) = plt.subplots(3,1,sharex=True)
fig, (ax0) = plt.subplots(1,1)
ax0.imshow(dens,origin='lower',vmin=0,vmax=20,extent=[-xext/2,xext/2,-yext/2,yext/2])
#ax1.imshow(np.rot90(dens)/5,origin='lower',cmap='viridis',vmin=0,vmax=5)
#ax2.imshow(np.rot90(vecpot),origin='lower',cmap='jet')
#plt.savefig('vecpot_test.png')
#identify midplane 
xmid = np.shape(vecpot)[1]/2
vecpot_slice = vecpot[:,xmid]
print(np.shape(vecpot_slice))

#print(vecpot_slice[0],vecpot_slice[-1])



xarr = np.linspace(0,np.size(vecpot_slice),np.size(vecpot_slice))
#diff = (vecpot_roll_plus - vecpot_roll_minus)/2.


my_localmin_list = localmin(vecpot_slice,15,100)
my_localmax_list = localmax(vecpot_slice,15,100)
num_xpoints = len(my_localmin_list)
num_plasmoids = len(my_localmax_list)
print('num of xpoints' + str(num_xpoints))
print('num of plasmoids' + str(num_plasmoids))


ax0.set_xlabel('$x \; (1000 \; c/\omega_{p})$')
ax0.set_ylabel('$y \; (1000 \; c/\omega_{p})$')
for i in range(num_xpoints):
    myx = my_localmin_list[i]
    #find closest val in xarr
    #myarg = np.argmin(np.abs(myx - xarr))
    #myy = vecpot_slice[myarg]
    ax0.scatter(convfac*(xmid+2)-xext/2.,convfac*(myx)-yext/2.,color="Red",marker='x')
    #print(i, myx, myy)
'''
for i in range(num_plasmoids):
    myx = my_localmax_list[i]
    #find closest val in xarr
    myarg = np.argmin(np.abs(myx - xarr))
    myy = vecpot_slice[myarg]
    ax1.scatter(myx,myy,color="Blue")
    #print(i, myx, myy)
'''
plt.savefig('old_xpoint_algorithm_untrig_bguide.png',bbox_inches='tight')

