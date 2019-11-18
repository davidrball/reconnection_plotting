import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
from tristan_funcs import recon_region

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})


myfld = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.7_pressuretens/output/flds.tot.020"

f = h5py.File(myfld)
#print(f.keys())

#pretty sure these are pressures, not temperatures
pxx = np.rot90(f['tmpxx'][0,:,:])

ylen = pxx.shape[1]
print(ylen)
ylow = int(ylen/6)
yup = int(5*ylow)
pxx = pxx[:,ylow:yup]
pyy = np.rot90(f['tmpyy'][0,ylow:yup,:])
pzz = np.rot90(f['tmpzz'][0,ylow:yup,:])

#number of particles divided out in pxx / pyy
#so I think they're actually temperatures

bx = np.rot90(f['bx'][0,ylow:yup,:])
by = np.rot90(f['by'][0,ylow:yup,:])
bz = np.rot90(f['bz'][0,ylow:yup,:])

dens = np.rot90(f['dens'][0,ylow:yup,:])
bdens = np.rot90(f['bdens'][0,ylow:yup,:])
pdens = np.rot90(f['pdens'][0,ylow:yup,:])



cury = recon_region(dens,pdens,bdens)

btot = np.sqrt(bx**2 + by**2 + bz**2)
bx /= btot
by /= btot
bz /= btot


#ptot = np.sqrt(pxx**2 + pyy**2 + pzz**2)
#ppar = np.sqrt(bx**2 * pxx**2 + by**2 * pyy **2 + bz**2 * pzz**2)
#pperp = np.sqrt(ptot**2 - ppar**2)

ptot = pxx + pyy + pzz
ppar = np.abs(bx*pxx) + np.abs(by*pyy) + np.abs(bz*pzz)
pperp = (ptot - ppar)
#pperp = np.abs((1-bx)*pxx) + np.abs((1-by)*pyy) + np.abs((1-bz)*pzz)
ptot *= dens
ppar *= dens
pperp *= dens


trat = pperp / ppar
betapar = 8*np.pi*ppar/btot**2

myshape = np.shape(trat)

trat_list = []
betapar_list = []
for i in range(myshape[0]):
    for j in range(myshape[1]):
        if cury[i][j] == 1:
            trat_list.append(np.log10(trat[i][j]))
            betapar_list.append(np.log10(betapar[i][j]))

#trat = np.log10(trat.flatten())
#betapar = np.log10(betapar.flatten())
xarr = np.logspace(0,3,100)
yarr = 1 - 2/xarr
yarr2 = 1/xarr + 1


H, xedges, yedges = np.histogram2d(betapar_list, trat_list, bins=40)

plt.imshow(H, origin='low',norm=LogNorm(),extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]],aspect=4)
plt.plot(np.log10(xarr), np.log10(yarr), linestyle='dashed',color="Black")
#plt.plot(np.log10(xarr),np.log10(yarr2),linestyle='dashed',color="Black")
plt.ylabel('$\log(T_{\perp}/T_{||})$')
plt.xlabel('$\log(\\beta_{||})$')
plt.ylim(-.75,.5)
plt.xlim(-3,4)
plt.colorbar()
plt.savefig('bg.7_brazil.png',dpi=300,bbox_inches='tight')


'''
print(np.min(trat),np.max(trat))
plt.imshow(np.log10(betapar),origin='lower')
#plt.colorbar(label='$\log(T_{\perp}/T_{||})$')
plt.colorbar(label='$\\beta_{||}$')
plt.savefig('bg.7_betapar.png',dpi=300,bbox_inches='tight')
'''

f.close()
