import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
from matplotlib.colors import PowerNorm

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

plt.set_cmap("viridis")
from beta_calc import beta_slice

def get_val_filename(input_string,fopen):
    str_list = list(fopen.keys())
    my_string = input_string
    if my_string in str_list:   
        my_list = np.array((fopen.get(my_string))[0])
        return my_list
    else:
        print("argument not found, your options are:")
        print(list(fopen.keys()))
        return 0
        


istep = 12 #downsampling
c_omp = 3 #cells per skin depth

sigma=.3
if sigma<1:
    sigma_str = str(sigma)[1:]
else:
    sigma_str = str(sigma)

'''
fld_base0 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide0/output/flds.tot.0"
fld_base1 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.1/output/flds.tot.0"
fld_base2 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.3/output/flds.tot.0"
fld_base3 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.7/output/flds.tot.0"
'''

'''
fld_base0 = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/flds.tot.0"
fld_base1 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/flds.tot.0"
fld_base2 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/flds.tot.0"
fld_base3 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.7/output/flds.tot.0"
'''

#sigma.1 results

fld_base0 = "../../tristan-mp_reconnection/guide_fields/pleiades/sig.1/delgam0000165/bguide0/output/flds.tot.0"                                               
fld_base1 = "../../tristan-mp_reconnection/guide_fields/pleiades/sig.1/delgam0000165/bguide.1/output/flds.tot.0"
fld_base2 = "../../tristan-mp_reconnection/guide_fields/pleiades/sig.1/delgam0000165/bguide.3/output/flds.tot.0"
fld_base3 = "../../tristan-mp_reconnection/guide_fields/pleiades/sig.1/delgam0000165/bguide.7/output/flds.tot.0"



tstart = 50
tfinal = 51

mybeta = .0003
mytheta = 0.000165
ppc = 4
xslice = 100

for t in range(tstart,tfinal):
    t_str = str(t)


    f0 = h5py.File(fld_base0+t_str,'r')
    f1 = h5py.File(fld_base1+t_str,'r')
    f2 = h5py.File(fld_base2+t_str,'r')
    f3 = h5py.File(fld_base3+t_str,'r')


    by_ref = f0['by'][0][0][0] #just grab b field strength from corner of domain, y component is reconnecting
    bz_ref0 = f0['bz'][0][0][0]
    bz_ref1 = f1['bz'][0][0][0]
    bz_ref2 = f2['bz'][0][0][0]
    bz_ref3 = f3['bz'][0][0][0]
    dens_ref = 4. #ppc

    #print('max v3y : ',np.max(f0['v3y']))

    mi = 1836.
    me = 1.
    vy0 = (mi+me)*f0['v3y'][0,:,:] - mi*f0['v3yi'][0,:,:]
    vy1 = (mi+me)*f1['v3y'][0,:,:] - mi*f1['v3yi'][0,:,:]
    vy2 = (mi+me)*f2['v3y'][0,:,:] - mi*f2['v3yi'][0,:,:]
    vy3 = (mi+me)*f3['v3y'][0,:,:] - mi*f3['v3yi'][0,:,:]

    vy0 = np.abs(np.roll(vy0,2,axis=1)-np.roll(vy0,-1,axis=1))
    vy1 = np.abs(np.roll(vy1,2,axis=1)-np.roll(vy1,-1,axis=1))
    vy2 = np.abs(np.roll(vy2,2,axis=1)-np.roll(vy2,-1,axis=1))
    vy3 = np.abs(np.roll(vy3,2,axis=1)-np.roll(vy3,-1,axis=1))

    #vy0 = f0['v3y'][0,:,:]
    #vy1 = f1['v3y'][0,:,:]
    #vy2 = f2['v3y'][0,:,:]
    #vy3 = f3['v3y'][0,:,:]

    by0 = f0['by'][0,:,:]
    bx0 = f0['bx'][0,:,:]
    bz0 = f0['bz'][0,:,:]
    by1 = f1['by'][0,:,:]
    bx1 = f1['bx'][0,:,:]
    bz1 = f1['bz'][0,:,:]
    by2 = f2['by'][0,:,:]
    bx2 = f2['bx'][0,:,:]
    bz2 = f2['bz'][0,:,:]
    by3 = f3['by'][0,:,:]
    bx3 = f3['bx'][0,:,:]
    bz3 = f3['bz'][0,:,:]
    dens0 = f0['dens'][0,:,:]
    dens1 = f1['dens'][0,:,:]
    dens2 = f2['dens'][0,:,:]
    dens3 = f3['dens'][0,:,:]
    sig_init = 0.1
    guide0 = 0
    guide1 = .1
    guide2 = .3
    guide3 = .7
    sigtot0 = sig_init*(1+guide0**2)
    sigtot1 = sig_init*(1+guide1**2)
    sigtot2 = sig_init*(1+guide2**2)
    sigtot3 = sig_init*(1+guide3**2)

    sigfinal0 = sigtot0 * ((bx0**2 + by0**2 +bz0**2)/(by_ref**2 + bz_ref0**2))/(dens0/dens_ref)
    sigfinal1 = sigtot1 * ((bx1**2 + by1**2 +bz1**2)/(by_ref**2 + bz_ref1**2))/(dens1/dens_ref)
    sigfinal2 = sigtot2 * ((bx2**2 + by2**2 +bz2**2)/(by_ref**2 + bz_ref2**2))/(dens2/dens_ref)
    sigfinal3 = sigtot3 * ((bx3**2 + by3**2 +bz3**2)/(by_ref**2 + bz_ref3**2))/(dens3/dens_ref)
    #so we have our sigma values, the alfven speed is just sqrt(sigma / (1+sigma))
    #calculating these:
    #including factor of 2 for stability criterion
    va0 = 2*np.rot90(vy0/np.sqrt(sigfinal0 / (1+sigfinal0)))
    va1 = 2*np.rot90(vy1/np.sqrt(sigfinal1 / (1+sigfinal1)))
    va2 = 2*np.rot90(vy2/np.sqrt(sigfinal2 / (1+sigfinal2)))
    va3 = 2*np.rot90(vy3/np.sqrt(sigfinal3 / (1+sigfinal3)))


    #divide through by sigtot if you want to express sigma as fraction of upstream

    xscan = 100
    xhlf = np.shape(va0)[0]/2
    xup = int(xhlf + xscan)
    xlow = int(xhlf-xscan)
    va0 = va0[xlow:xup,:]
    va1 = va1[xlow:xup,:]
    va2 = va2[xlow:xup,:]
    va3 = va3[xlow:xup,:]

    #xlen = 2*xscan
    xlen = np.shape(va0)[0]
    ylen = np.shape(va0)[1]


    xext = xlen*istep / c_omp
    yext = ylen*istep / c_omp

    xext /= 1000
    yext /= 1000
    
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(4,1,sharex=True)
    #myvmax = 8
    myvmin = 0
    myvmax = 10
    im0 = ax0.imshow(va0,vmin=myvmin,vmax=myvmax, extent = [-yext/2,yext/2,-xext/2,xext/2],origin='lower')
    im1 = ax1.imshow(va1,vmin=myvmin,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2], origin='lower')
    im2 = ax2.imshow(va2,vmin=myvmin,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2],origin='lower')
    im3 = ax3.imshow(va3,vmin=myvmin,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2],origin='lower')
#linear scale

    ax0.set_adjustable('box-forced')
    ax1.set_adjustable('box-forced')
    ax2.set_adjustable('box-forced')
    ax3.set_adjustable('box-forced')
   

    ax3.set_xlabel('x ($1000\;c/\omega_{p}$)',size=14)
    #ax1.set_ylabel('y ($1000\;c/\omega_{p}$)',size=14)
    #ax2.set_ylabel('y ($1000\;c/\omega_{p}$)',size=14)
    #ax3.set_ylabel('y ($1000\;c/\omega_{p}$)',size=14)

#plt.colorbar()
#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([.81, .12, .05, .75])

    cbar_ax = fig.add_axes([.125,.93,.77,.03])
#cbar_ax = fig.add_axes([.28,.94,.46,.02])
#fig.colorbar(im3, orientation="horizontal",cax=cbar_ax,label="Density (particles per cell)")
    cb= fig.colorbar(im0,cax=cbar_ax,orientation="horizontal")#,label="Density (particles per cell)")
    cbar_ax.set_xlabel("$|\Delta v_{xe}|/v_{A}/2$",fontsize=14, labelpad=-45)

    cbytick_obj = plt.getp(cb.ax.axes,'xticklabels')
    plt.setp(cbytick_obj,fontsize='x-small')
#cb.ax.set_xticklabels([])
#cb.ax.tick_params(axis='x',direction='out',labeltop=True, top=True)

#cb.ax.set_xticklabels([0,4,8,12,16,20])
    #fig.savefig('../../tristan_plots/guide_flds/sig.3_delgam.0005/untriggered/'+t_str+'.png',dpi=300,bbox_inches='tight')
    fig.savefig('electron_KH_test_t50.png',dpi=300,bbox_inches='tight')

    plt.close()


