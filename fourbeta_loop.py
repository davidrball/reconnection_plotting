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


fld_base0 = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/flds.tot.0"
fld_base1 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/flds.tot.0"
fld_base2 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/flds.tot.0"
fld_base3 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.7/output/flds.tot.0"



tstart = 30
tfinal = 31

mybeta = .003
mytheta = 0.0005
ppc = 4
xslice = 100

for t in range(tstart,tfinal):
    t_str = str(t)


    #f0 = h5py.File(fld_base0+t_str,'r')
    #f1 = h5py.File(fld_base1+t_str,'r')
    #f2 = h5py.File(fld_base2+t_str,'r')
    #f3 = h5py.File(fld_base3+t_str,'r')
    btot0 = beta_slice(fld_base0 + t_str,mybeta,mytheta,ppc,xslice)
    btot1= beta_slice(fld_base1 + t_str,mybeta,mytheta,ppc,xslice)
    btot2 = beta_slice(fld_base2 + t_str,mybeta,mytheta,ppc,xslice)
    btot3 = beta_slice(fld_base3 + t_str,mybeta,mytheta,ppc,xslice)




    btot0 = np.log10(np.rot90(btot0))
    btot1 = np.log10(np.rot90(btot1))
    btot2 = np.log10(np.rot90(btot2))
    btot3 = np.log10(np.rot90(btot3))
    xscan = 100
    xlen = 2*xscan
    ylen = np.shape(btot1)[1]


    xext = xlen*istep / c_omp
    yext = ylen*istep / c_omp

    xext /= 1000
    yext /= 1000
    
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(4,1,sharex=True)
    #myvmax = 8
    myvmin = -3
    myvmax = 2
    im0 = ax0.imshow(btot0,vmin=myvmin,vmax=myvmax, extent = [-yext/2,yext/2,-xext/2,xext/2],origin='lower')
    im1 = ax1.imshow(btot1,vmin=myvmin,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2], origin='lower')
    im2 = ax2.imshow(btot2,vmin=myvmin,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2],origin='lower')
    im3 = ax3.imshow(btot3,vmin=myvmin,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2],origin='lower')
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
    cbar_ax.set_xlabel("$\\beta$",fontsize=14, labelpad=-45)

    cbytick_obj = plt.getp(cb.ax.axes,'xticklabels')
    plt.setp(cbytick_obj,fontsize='x-small')
#cb.ax.set_xticklabels([])
#cb.ax.tick_params(axis='x',direction='out',labeltop=True, top=True)

#cb.ax.set_xticklabels([0,4,8,12,16,20])
    #fig.savefig('../../tristan_plots/guide_flds/sig.3_delgam.0005/untriggered/'+t_str+'.png',dpi=300,bbox_inches='tight')
    fig.savefig('betafld_lowbeta.png',dpi=300,bbox_inches='tight')
    '''
    plt.figure()
    ax = plt.gca()
    im = ax.imshow(np.rot90(dens),origin='lower',vmin=0,vmax=20, extent = [0,yext,0,xext])
    plt.annotate('Time ($L/V_{a}$)='+my_time_str ,(5, 20),color="white")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right",size="5%",pad=0.05)
    plt.colorbar(im,cax=cax,label="Density (particles per cell)")
    
    ax.set_xlabel('$c/\omega_{p}$',size=16)
    ax.set_ylabel('$c/\omega_{p}$',size=16)
    #plt.savefig('fld_plots/triggered_example/'+t_str+'.png',dpi=300,bbox_inches='tight')
    plt.savefig('fld_plots/sig.3/delgam00015/'+str(t)+'.png',dpi=300,bbox_inches='tight')
    plt.close()
    f.close()
    '''
    #f0.close()
    #f1.close()
    #f2.close()
    #f3.close()
    plt.close()


