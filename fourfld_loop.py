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

delgam=.00016
delgam_str = str(delgam)

fld_base0 = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/flds.tot.0"

fld_base1 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/flds.tot.0"

fld_base2 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/flds.tot.0"

fld_base3 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.7/output/flds.tot.0"

tstart = 30
tfinal = 31

for t in range(tstart,tfinal):
    t_str = str(t)


    f0 = h5py.File(fld_base0+t_str,'r')
    f1 = h5py.File(fld_base1+t_str,'r')
    f2 = h5py.File(fld_base2+t_str,'r')
    f3 = h5py.File(fld_base3+t_str,'r')

    bx0 = f0['bx'][0,:,:]
    by0 = f0['by'][0,:,:]
    bz0 = f0['bz'][0,:,:]
    btot0 = bx0**2 + by0**2 + bz0**2
    btot0 = bz0**2
    bx1 = f1['bx'][0,:,:]
    by1 = f1['by'][0,:,:]
    bz1 = f1['bz'][0,:,:]
    btot1 = bx1**2 + by1**2 + bz1**2
    btot1 = bz1**2
    bx2 = f2['bx'][0,:,:]
    by2 = f2['by'][0,:,:]
    bz2 = f2['bz'][0,:,:]
    btot2 = bx2**2 + by2**2 + bz2**2
    btot2 = bz2**2
    bx3 = f3['bx'][0,:,:]
    by3 = f3['by'][0,:,:]
    bz3 = f3['bz'][0,:,:]
    btot3 = bx3**2 + by3**2 + bz3**2
    btot3 = bz3**2
    xscan = 100

    x0len = np.shape(btot0)[1]
    x0hlf = x0len/2
    x0low = int(x0hlf - xscan)
    x0up = int(x0hlf+xscan)

    x1len = np.shape(btot1)[1]
    x1hlf = x1len/2
    x1low = int(x1hlf - xscan)
    x1up = int(x1hlf+xscan)

    x2len = np.shape(btot2)[1]
    x2hlf = x2len/2
    x2low = int(x2hlf - xscan)
    x2up = int(x2hlf+xscan)

    x3len = np.shape(btot3)[1]
    x3hlf = x3len/2
    x3low = int(x3hlf - xscan)
    x3up = int(x3hlf+xscan)




    btot0 = np.log10(np.rot90(btot0[:,x0low:x0up]))
    btot1 = np.log10(np.rot90(btot1[:,x1low:x1up]))
    btot2 = np.log10(np.rot90(btot2[:,x2low:x2up]))
    btot3 = np.log10(np.rot90(btot3[:,x3low:x3up]))

    xlen = 2*xscan
    ylen = np.shape(btot1)[1]


    xext = xlen*istep / c_omp
    yext = ylen*istep / c_omp

    xext /= 1000
    yext /= 1000
    '''
    L = ylen*istep
    c = .45
    sigma=.3
    interval = 2000
    va = np.sqrt(sigma / (1+sigma))
    alf_cross = L / (c*va*interval)
    my_time1 = (t1-bias) / alf_cross
    my_time_str1 = str(my_time1)[:4]
    my_time2 = (t2-bias) / alf_cross
    my_time_str2 = str(my_time2)[:4]
    my_time3 = (t3-bias) / alf_cross
    my_time_str3 = str(my_time3)[:4]
    '''
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(4,1,sharex=True)
    #myvmax = 8
    myvmin = -3
    myvmax = 2
    im0 = ax0.imshow(btot0,vmin=myvmin,vmax=myvmax, extent = [-yext/2,yext/2,-xext/2,xext/2],origin='lower')
    im1 = ax1.imshow(btot1,vmin=myvmin,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2], origin='lower')
    im2 = ax2.imshow(btot2,vmin=myvmin,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2],origin='lower')
    im3 = ax3.imshow(btot3,vmin=myvmin,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2],origin='lower')
#linear scale

    '''
    #nonlinear color scale
    im1 = ax1.imshow((dens1/4.),norm=PowerNorm(gamma=1,vmin=0,vmax=8),extent=[-yext/2,yext/2,-xext/2,xext/2])
    im2 = ax2.imshow((dens2/4.),norm=PowerNorm(gamma=1,vmin=0,vmax=8),extent=[-yext/2,yext/2,-xext/2,xext/2])
    im3 = ax3.imshow((dens3/4.),norm=PowerNorm(gamma=1,vmin=0,vmax=8),extent=[-yext/2,yext/2,-xext/2,xext/2])
    '''


#print(np.min(dens1))
#print(np.min(dens2))
#print(np.min(dens3))

#im1 = ax1.imshow(dens1/4.,norm=LogNorm,extent=[-yext/2,yext/2,-xext/2,xext/2])                     
#im2 = ax2.imshow(dens2/4.,norm=LogNorm,extent=[-yext/2,yext/2,-xext/2,xext/2])                     
#im3 = ax3.imshow(dens3/4.,norm=LogNorm,extent=[-yext/2,yext/2,-xext/2,xext/2])

#ax1.set_adjustable('box-forced')
#ax2.set_adjustable('box-forced')
#ax3.set_adjustable('box-forced')        

#ax1.set_yticks([-1,-.5,0,.5,1])
    ax0.set_adjustable('box-forced')
    ax1.set_adjustable('box-forced')
    ax2.set_adjustable('box-forced')
    ax3.set_adjustable('box-forced')
    #ax1.annotate("$t_{A}="+my_time_str1+"$",(-.3,.4),color="White",size=14, fontweight='bold')
    #ax2.annotate("$t_{A}="+my_time_str2+"$",(-.3,.4),color="White",size=14, fontweight='bold')
    #ax3.annotate("$t_{A}="+my_time_str3+"$",(-.3,.4),color="White",size=14, fontweight='bold')
    
    #ax0.annotate("$B_{g}/B_{0}=0.3$",(-.4,.25),color="White",size=14,fontweight='bold')
    #ax1.annotate("$B_{g}/B_{0}=0.3$",(-.4,.25),color="White",size=14,fontweight='bold')    
    #ax2.annotate("$B_{g}/B_{0}=0.7$" ,(-.4,.25),color="White",size=14,fontweight='bold')
    #ax3.annotate("$B_{g}/B_{0}=1$",(-.4,.25),color="White",size=14,fontweight='bold')
#rect = patches.Rectangle((-2.7,.35),.2,.2,facecolor="White")
#rect2 = patches.Rectangle((-2.7,.35),.2,.2,facecolor="White")
#rect3 = patches.Rectangle((-2.7,.35),.2,.2,facecolor="White")


#ax1.add_patch(rect)
#ax2.add_patch(rect2)
#ax3.add_patch(rect3)
#ax1.annotate("a)",(-2.7,.4),color="Black",size=14)
#ax2.annotate("b)",(-2.7,.4),color="Black",size=14)
#ax3.annotate("c)",(-2.7,.4),color="Black",size=14)



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
    cbar_ax.set_xlabel("magnetic energy density",fontsize=14, labelpad=-45)

    cbytick_obj = plt.getp(cb.ax.axes,'xticklabels')
    plt.setp(cbytick_obj,fontsize='x-small')
#cb.ax.set_xticklabels([])
#cb.ax.tick_params(axis='x',direction='out',labeltop=True, top=True)

#cb.ax.set_xticklabels([0,4,8,12,16,20])
    #fig.savefig('../../tristan_plots/guide_flds/sig.3_delgam.0005/untriggered/'+t_str+'.png',dpi=300,bbox_inches='tight')
    fig.savefig('density_lowbeta.png',dpi=300,bbox_inches='tight')
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
    f0.close()
    f1.close()
    f2.close()
    f3.close()
    plt.close()


