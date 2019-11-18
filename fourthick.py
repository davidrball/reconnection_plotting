import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
from matplotlib.colors import PowerNorm
from thickness_calc import return_thickness_list
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

'''
fld_base0 = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/flds.tot.0"
fld_base1 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/flds.tot.0"
fld_base2 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/flds.tot.0"
fld_base3 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.7/output/flds.tot.0"
'''

#highbeta
fld_base0 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide0/output/flds.tot.0"
fld_base1 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.1/output/flds.tot.0"
fld_base2 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.3/output/flds.tot.0"
fld_base3 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.7/output/flds.tot.0"


tstart = 15
tfinal = 45
tspan = tfinal - tstart
f0thicklist = []
f1thicklist = []
f2thicklist = []
f3thicklist = []
for t in range(tstart,tfinal):
    t_str = str(t)
    

    f0 = fld_base0+t_str
    f1 = fld_base1+t_str
    f2 = fld_base2+t_str
    f3 = fld_base3+t_str

    f0thicklistsub, f0mean, f0std = return_thickness_list(f0)
    f1thicklistsub, f1mean, f1std = return_thickness_list(f1)
    f2thicklistsub, f2mean, f2std = return_thickness_list(f2)
    f3thicklistsub, f3mean, f3std = return_thickness_list(f3)

    for thick in f0thicklistsub:
        f0thicklist.append(thick)
    for thick in f1thicklistsub:
        f1thicklist.append(thick)
    for thick in f2thicklistsub:
        f2thicklist.append(thick)
    for thick in f3thicklistsub:
        f3thicklist.append(thick)

    
minbin = min(min(f0thicklist), min(f1thicklist),min(f2thicklist),min(f3thicklist))    

maxbin = max(max(f0thicklist), max(f1thicklist),max(f2thicklist),max(f3thicklist))


binnum = 15
mybins = np.linspace(minbin,maxbin,binnum)
f0hist = np.zeros(binnum)
f1hist = np.zeros(binnum)
f2hist = np.zeros(binnum)
f3hist = np.zeros(binnum)
for i in range(binnum-1):
    low = mybins[i]
    up = mybins[i+1]
    for width in f0thicklist:
        if low < width < up:
            f0hist[i]+=1
    for width in f1thicklist:
        if low < width < up:
            f1hist[i]+=1
    for width in f2thicklist:
        if low < width < up:
            f2hist[i]+=1
    for width in f3thicklist:
        if low < width < up:
            f3hist[i]+=1
plt.plot(mybins,f0hist/tspan,label="$B_{g}=0$")
plt.plot(mybins,f1hist/tspan,label="$B_{g}=0.1$")
plt.plot(mybins,f2hist/tspan,label="$B_{g}=0.3$")
plt.plot(mybins,f3hist/tspan,label="$B_{g}=0.7$")
plt.legend(frameon=False,loc='upper right')
plt.xlabel('Sheet Thickness')
plt.ylabel('Time integrated counts')
plt.savefig('highbeta_thickhist.png',dpi=300,bbox_inches='tight')


