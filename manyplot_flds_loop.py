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

fld_base0 = "../../tristan-mp_reconnection/16k_triggered_finals/sig.1/delgam00015/output/flds.tot.060"

fld_base1 = "../../tristan-mp_reconnection/guide_fields/sig.1/simulation_matrix/delgam000165/bguide.1/output/flds.tot.060"

fld_base2 = "../../tristan-mp_reconnection/guide_fields/sig.1/simulation_matrix/delgam000165/bguide.3/output/flds.tot.0"


'''
#for untriggered sig.3 guide fields
fld_base0 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide0_untriggered/output/flds.tot.00"                               
fld_base1 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3_untriggered/output/flds.tot.00"                                             
fld_base2 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide1_untriggered/output/flds.tot.00" 
fld_base3 = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide2_untriggered/output/flds.tot.00"
'''



#for high beta
#fld_base0 = "../../tristan-mp_reconnection/16k_triggered_finals/sig1/delgam.2/output/flds.tot.0"
fld_base0 = "../../tristan-mp_reconnection/guide_fields/sig1/simulation_matrix/delgam0165/bguide.1/output/flds.tot.0"

fld_base1 = "../../tristan-mp_reconnection/guide_fields/sig1/simulation_matrix/delgam0165/bguide.3/output/flds.tot.0"     
fld_base2 = "../../tristan-mp_reconnection/guide_fields/sig1/simulation_matrix/delgam0165/bguide.7/output/flds.tot.0"     


#for sigma=0.3 highbeta
fld_base0 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.3/output/flds.tot.0"
fld_base1 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide1/output/flds.tot.0"
fld_base2 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide3/output/flds.tot.0"



fld_base0 = "../../tristan_acc-mec_Ez/8k_bguide0_triggered_stride1_thresh2/output/flds.tot."
fld_base1 = "../../tristan_acc-mec_Ez/8k_bguide.1_allprts/output/flds.tot."
fld_base2 = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/flds.tot."



 


mybase = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/"

t0 = 4
t1 = 8
t2 = 50

interval = 2000

t0str = "%03d" % t0
t1str = "%03d" % t1
t2str = "%03d" % t2

fld_base0 +=t0str
fld_base1 += t1str
fld_base2 += t2str

fld_base0 = mybase + "bguide.3_untriggered_stripe20/output/flds.tot."+t0str
fld_base1 = mybase + "bguide.3_untriggered/output/flds.tot."+t1str
fld_base2 = mybase + "bguide.3_untriggered_stripe60/output/flds.tot."+t2str

#fld_base0 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/boxsize_test/triggered/4k/output/flds.tot.030"
#fld_base1 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/flds.tot.060"
#fld_base2 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/flds.tot.070"


tstart = 17
tfinal = 18

for t in range(tstart,tfinal):
    t_str = "%03d"%t


    #f0 = h5py.File(fld_base0+t_str,'r')
    #f1 = h5py.File(fld_base1+t_str,'r')
    #f2 = h5py.File(fld_base2+t_str,'r')
    #f3 = h5py.File(fld_base3+t_str,'r')

    #for separate times
    f0 = h5py.File(fld_base0,'r')
    f1 = h5py.File(fld_base1,'r')
    f2 = h5py.File(fld_base2,'r')


    dens0 = get_val_filename('dens',f0)
    dens1 = get_val_filename('dens',f1)
    dens2 = get_val_filename('dens',f2)
    #dens3 = get_val_filename('dens',f3)


    xscan = 75
    
    print('shape is : ',np.shape(dens0))
    x0len = np.shape(dens0)[1]
    x0hlf = x0len/2
    x0low = int(x0hlf - xscan)
    x0up = int(x0hlf+xscan)
    print(x0len)
    x1len = np.shape(dens1)[1]
    x1hlf = x1len/2
    x1low = int(x1hlf - xscan)
    x1up = int(x1hlf+xscan)
    print(x1len)
    x2len = np.shape(dens2)[1]
    x2hlf = x2len/2
    x2low = int(x2hlf - xscan)
    x2up = int(x2hlf+xscan)
    print(x2len)
    #x3len = np.shape(dens3)[1]
    #x3hlf = x3len/2
    #x3low = int(x3hlf - xscan)
    #x3up = int(x3hlf+xscan)




    dens0 = np.rot90(dens0[:,x0low:x0up])
    dens1 = np.rot90(dens1[:,x1low:x1up])
    dens2 = np.rot90(dens2[:,x2low:x2up])
    #dens3 = np.rot90(dens3[:,x3low:x3up])

    xlen = 2*xscan
    ylen = np.shape(dens1)[1]


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
    fig, (ax0, ax1, ax2) = plt.subplots(3,1)#,sharex=True)
    myvmax = 5

    im0 = ax0.imshow((dens0/4.), vmin=0, vmax=myvmax, extent = [-yext/2,yext/2,-xext/2,xext/2],origin='lower')
    im1 = ax1.imshow((dens1/4.),vmin=0,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2], origin='lower')
    im2 = ax2.imshow((dens2/4.),vmin=0,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2],origin='lower')
    #im3 = ax3.imshow((dens3/4.),vmin=0,vmax=myvmax,extent=[-yext/2,yext/2,-xext/2,xext/2],origin='lower')
#linear scale


    #ax1.annotate("$t_{A}="+my_time_str1+"$",(-.3,.4),color="White",size=14, fontweight='bold')
    #ax2.annotate("$t_{A}="+my_time_str2+"$",(-.3,.4),color="White",size=14, fontweight='bold')
    #ax3.annotate("$t_{A}="+my_time_str3+"$",(-.3,.4),color="White",size=14, fontweight='bold')
    
    ax0.annotate("$B_{g}/B_{0}=0$",(-.2,.21),color="White",size=14,fontweight='bold')
    ax1.annotate("$B_{g}/B_{0}=0.1$",(-.2,.21),color="White",size=14,fontweight='bold')    
    ax2.annotate("$B_{g}/B_{0}=0.3$" ,(-.2,.21),color="White",size=14,fontweight='bold')
    #ax3.annotate("$B_{g}/B_{0}=1$",(-.4,.25),color="White",size=14,fontweight='bold')


    #ax0.annotate("$\Delta=6.66 \; c/\omega_{p}$", (-.5,.3),color="White",size=14,fontweight='bold')
    #ax1.annotate("$\Delta=13.33 \; c/\omega_{p}$", (-.5, .3), color="White",size=14,fontweight='bold')
    #ax2.annotate("$\Delta=20.0 \; c/\omega_{p}$", (-.5, .3), color="White",size=14, fontweight='bold')

    ax0.set_xticks([])
    ax1.set_xticks([])

    ax2.set_xlabel('x ($1000\;c/\omega_{p}$)',size=14)
    ax0.set_ylabel('y ($1000\;c/\omega_{p}$)',size=14)
    ax1.set_ylabel('y ($1000\;c/\omega_{p}$)',size=14)
    ax2.set_ylabel('y ($1000\;c/\omega_{p}$)',size=14)
    #ax3.set_ylabel('y ($1000\;c/\omega_{p}$)',size=14)

#plt.colorbar()
#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([.81, .12, .05, .75])

    cbar_ax = fig.add_axes([.125,.93,.77,.03])
#cbar_ax = fig.add_axes([.28,.94,.46,.02])
#fig.colorbar(im3, orientation="horizontal",cax=cbar_ax,label="Density (particles per cell)")
    cb= fig.colorbar(im1,cax=cbar_ax,orientation="horizontal")#,label="Density (particles per cell)")
    cbar_ax.set_xlabel("Density ($N_{ppc}$ / $N_{ppc0}$)",fontsize=14, labelpad=-45)

    cbytick_obj = plt.getp(cb.ax.axes,'xticklabels')
    plt.setp(cbytick_obj,fontsize='x-small')
#cb.ax.set_xticklabels([])
#cb.ax.tick_params(axis='x',direction='out',labeltop=True, top=True)

#cb.ax.set_xticklabels([0,4,8,12,16,20])
    #fig.savefig('../../tristan_plots/guide_flds/sig.3_delgam.0005/with_.1/'+t_str+'.png',dpi=300,bbox_inches='tight')

    #fig.savefig('diffguide_flds.pdf',dpi=300,bbox_inches='tight')
    #fig.savefig('diffguide_flds.png',dpi=300,bbox_inches='tight')
    fig.savefig('testing_thicktimes.png')

    f0.close()
    f1.close()
    f2.close()
    #f3.close()
    plt.close()


