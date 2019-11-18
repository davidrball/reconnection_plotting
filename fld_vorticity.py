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

plt.set_cmap("RdBu")
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

    mi = 1836.
    me = 1.
    
    vy0 = np.rot90(f0['v3y'][0,:,:])
    vy1 = np.rot90(f1['v3y'][0,:,:])
    vy2 = np.rot90(f2['v3y'][0,:,:])
    vy3 = np.rot90(f3['v3y'][0,:,:])

    vx0 = np.rot90(f0['v3x'][0,:,:])
    vx1 = np.rot90(f1['v3x'][0,:,:])
    vx2 = np.rot90(f2['v3x'][0,:,:])
    vx3 = np.rot90(f3['v3x'][0,:,:])

    
    xscan = 100
    xhlf = np.shape(vy0)[0]/2
    xup = int(xhlf + xscan)
    xlow = int(xhlf-xscan)
    #print(xlow,xup)
    #print(np.shape(vx0))
    
    vy0 = vy0[xlow:xup,:]
    vy1 = vy1[xlow:xup,:]
    vy2 = vy2[xlow:xup,:]
    vy3 = vy3[xlow:xup,:]

    vx0 = vx0[xlow:xup,:]
    vx1 = vx1[xlow:xup,:]
    vx2 = vx2[xlow:xup,:]
    vx3 = vx3[xlow:xup,:]


    vort0 = return_vorticity(vx0,vy0)
    vort1 = return_vorticity(vx1,vy1)
    vort2 = return_vorticity(vx2,vy2)
    vort3 = return_vorticity(vx3,vy3)
    print(np.min(vort0),np.max(vort0))
    #print(np.shape(vx0))


    #xlen = 2*xscan
    xlen = np.shape(vx0)[0]
    ylen = np.shape(vx0)[1]


    xext = xlen*istep / c_omp
    yext = ylen*istep / c_omp

    xext /= 1000
    yext /= 1000
    
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(4,1,sharex=True)
    #myvmax = 8
    myvmin = -.1
    myvmax = .1
    im0 = ax0.imshow(vort0,vmin=myvmin,vmax=myvmax, extent = [-yext/2,yext/2,-xext/2,xext/2],origin='lower')
    im1 = ax1.imshow(vort1,vmin=myvmin,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2], origin='lower')
    im2 = ax2.imshow(vort2,vmin=myvmin,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2],origin='lower')
    im3 = ax3.imshow(vort3,vmin=myvmin,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2],origin='lower')
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
    #cbar_ax.set_xlabel("$|\Delta v_{xe}|/v_{A}/2$",fontsize=14, labelpad=-45)

    cbytick_obj = plt.getp(cb.ax.axes,'xticklabels')
    plt.setp(cbytick_obj,fontsize='x-small')
#cb.ax.set_xticklabels([])
#cb.ax.tick_params(axis='x',direction='out',labeltop=True, top=True)

#cb.ax.set_xticklabels([0,4,8,12,16,20])
    #fig.savefig('../../tristan_plots/guide_flds/sig.3_delgam.0005/untriggered/'+t_str+'.png',dpi=300,bbox_inches='tight')
    fig.savefig('vort_test.png',dpi=300,bbox_inches='tight')

    plt.close()


