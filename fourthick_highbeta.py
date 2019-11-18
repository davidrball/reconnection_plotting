import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
from matplotlib.colors import PowerNorm
from thickness_calc import return_thickness_slice
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
fld_base4 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide1/output/flds.tot.0"

tstart = 25
tfinal = 45
tspan = tfinal - tstart
f0thicklist = []
f1thicklist = []
f2thicklist = []
f3thicklist = []
f4thicklist = []
for t in range(tstart,tfinal):
    t_str = str(t)
    

    f0 = fld_base0+t_str
    f1 = fld_base1+t_str
    f2 = fld_base2+t_str
    f3 = fld_base3+t_str
    f4 = fld_base4+t_str
    f0thicklistsub = return_thickness_slice(f0)
    f1thicklistsub = return_thickness_slice(f1)
    f2thicklistsub = return_thickness_slice(f2)
    f3thicklistsub = return_thickness_slice(f3)
    f4thicklistsub = return_thickness_slice(f4)
    for thick in f0thicklistsub:
        f0thicklist.append(thick)
    for thick in f1thicklistsub:
        f1thicklist.append(thick)
    for thick in f2thicklistsub:
        f2thicklist.append(thick)
    for thick in f3thicklistsub:
        f3thicklist.append(thick)
    for thick in f4thicklistsub:
        f4thicklist.append(thick)
    
minbin = min(min(f0thicklist), min(f1thicklist),min(f2thicklist),min(f3thicklist),min(f4thicklist))    

maxbin = max(max(f0thicklist), max(f1thicklist),max(f2thicklist),max(f3thicklist), max(f4thicklist))


binnum = 12
mybins = np.linspace(minbin,maxbin,binnum)
f0hist = np.zeros(binnum)
f1hist = np.zeros(binnum)
f2hist = np.zeros(binnum)
f3hist = np.zeros(binnum)
f4hist = np.zeros(binnum)
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
    for width in f4thicklist:
        if low < width < up:
            f4hist[i]+=1
plt.plot(mybins,f0hist/tspan,label="$B_{g}=0$")
plt.plot(mybins,f1hist/tspan,label="$B_{g}=0.1$")
plt.plot(mybins,f2hist/tspan,label="$B_{g}=0.3$")
plt.plot(mybins,f3hist/tspan,label="$B_{g}=0.7$")
plt.plot(mybins,f4hist/tspan,label="$B_{g}=1$")
plt.legend(frameon=False,loc='upper right')
plt.xlabel('Sheet Thickness')
plt.ylabel('Time integrated counts')
plt.savefig('highbeta_thickhist_slice.png',dpi=300,bbox_inches='tight')


