import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
from matplotlib.colors import PowerNorm
from vorticity_calc import return_vorticity


plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})


myfld = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide0_pressuretens/output/flds.tot.045"

f = h5py.File(myfld)
#print(f.keys())
tmpxx = np.rot90(f['tmpxx'][0,:,:])
tmpyy = np.rot90(f['tmpyy'][0,:,:])
tmpzz = np.rot90(f['tmpzz'][0,:,:])
#tmpxy = np.rot90(f['tmpxy'][0,:,:])
#tmpxz = np.rot90(f['tmpxz'][0,:,:])
#tmpyz = np.rot90(f['tmpyz'][0,:,:])



mymin = min(np.min(tmpxx),np.min(tmpyy),np.min(tmpzz),np.min(tmpxy),np.min(tmpyz),np.min(tmpxz))
mymax = max(np.max(tmpxx),np.max(tmpyy),np.max(tmpzz),np.max(tmpxy),np.max(tmpyz),np.max(tmpxz))

fig, (ax0, ax1, ax2, ax3, ax4, ax5) = plt.subplots(6,1,sharex=True)

ax0.imshow(tmpxx,origin='lower',vmin=mymin,vmax=mymax)
ax1.imshow(tmpyy,origin='lower',vmin=mymin,vmax=mymax)
ax2.imshow(tmpzz,origin='lower',vmin=mymin,vmax=mymax)
ax3.imshow(tmpxy,origin='lower',vmin=mymin,vmax=mymax)
ax4.imshow(tmpxz,origin='lower',vmin=mymin,vmax=mymax)
ax5.imshow(tmpyz,origin='lower',vmin=mymin,vmax=mymax)
plt.savefig('testing_tmp.png',dpi=300,bbox_inches='tight')






f.close()
