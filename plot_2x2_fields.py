import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})
plt.set_cmap("viridis")

from timespec_func import plot_timespec


#define the simulation matrix, right now it's 3x3
#first index picks out guide field strength, 2nd is temp




#guide field 0, delgam=.00005, .0005, .005
#matrix00 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_run2/output/flds.tot."
matrix00 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/bguide.1_triggered_stride1_gamthresh50/output/flds.tot."
matrix01 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_bguide.3/output/flds.tot."


#matrix10 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_untriggered/output/flds.tot."
matrix10 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/bguide.1_untriggered_stride1_gamthresh50/output/flds.tot."


#matrix11 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_untriggered_bguide.3/output/flds.tot."
matrix11 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/bguide.3_untriggered_stride1_gamthresh50/output/flds.tot."


#matrix_namelist = [[matrix00, matrix01], [matrix10,matrix11]]

fig, ((ax0, ax1), (ax2, ax3))  = plt.subplots(2,2,sharex=True,sharey=True,figsize=(14,5))
plt.subplots_adjust(hspace=.01,wspace=.01)
axlist = [ax0, ax1, ax2, ax3]
mylen = len(axlist)
t = 12
tstr = '%03d' % t

matrix00 += tstr
matrix01 += tstr
matrix10 += tstr
matrix11 += tstr

#matrix_namelist = [[matrix00, matrix10], [matrix01,matrix11]]

#reordering to have high guide field at top
matrix_namelist = [[matrix01,matrix11],[matrix00,matrix10]]

xscan = 100
istep = 12
c_omp = 3


yext = xscan*istep/c_omp / 1000

for i in range(len(matrix_namelist)):
    for j in range(len(matrix_namelist[0])): #assuming that the matrix is square 
        myname = matrix_namelist[i][j]
        if myname == "pass":
            pass
        else:
            myax_ind = 2*i + j
            myax = axlist[myax_ind]
            mydens = np.rot90(h5py.File(myname,'r')['dens'][0,:,:])
            myshape = np.shape(mydens)
            xlen = myshape[1]
            ylen = myshape[0]
            yhlf = ylen/2.
            ylow = int(yhlf - xscan)
            yup = int(yhlf + xscan)
            mydens = mydens[ylow:yup,:]

            xext = xlen*istep/c_omp/1000
        #title = str(betaeff_list[myax_ind])
            myim = myax.imshow(mydens/4., origin='lower', vmin=0,vmax=5,extent = [-xext/2, xext/2, -yext, yext])
            myax.set_adjustable('box-forced')
        #myax.set_xlim(myxlim_low,myxlim_up)
            #myax.minorticks_off()
        #myax.set_yticks([])
        #myax.set_xticks([])
        #myax.set_xlim(1e0,1e4)
        #myax.set_xticks([1e0,1e1,1e2,1e3])
        #myax.set_yticks([1e2,1e4,1e6])
        #myfile.close()
cbar_ax = fig.add_axes([.125,.93,.775,.03])
cb= fig.colorbar(myim,cax=cbar_ax,orientation="horizontal")
cbar_ax.set_xlabel("Density ($N_{ppc}$ / $N_{ppc0}$)",fontsize=14,labelpad=-45)

cbytick_obj = plt.getp(cb.ax.axes,'xticklabels')
plt.setp(cbytick_obj,fontsize='x-small')

#cbar_ax = fig.add_axes([.905,.38,.015,.23])
#cb= fig.colorbar(im2,cax=cbar_ax,orientation="vertical",ticks=[-1,-.5,0,.5, 1])
#cbar_ax.set_xlabel("Density (particles per cell / $N_{ppc}$)",fontsize=14, labelpad=-4)

cbytick_obj = plt.getp(cb.ax.axes,'xticklabels')
plt.setp(cbytick_obj,fontsize='x-small')


#ax0.set_ylabel('$B_{g}/B_{0}=0$')
#ax3.set_ylabel('$B_{g}/B_{0}=0.1$')
#ax6.set_ylabel('$B_{g}/B_{0}=0.3$')

#ax6.set_xlabel('$\\theta=0.00005$')
#ax7.set_xlabel('$\\theta=0.0005$')
#ax8.set_xlabel('$\\theta=0.005$')
xstr = "$x \; (1000 \; c/\omega_{p})$"
ystr = "$y \; (1000 \; c/\omega_{p})$"
ax0.set_ylabel(ystr)                                        
ax2.set_ylabel(ystr)
ax2.set_xlabel(xstr)
ax3.set_xlabel(xstr)                                      

annx = -.4
anny=.3
ax0.annotate('$B_{g}/B_{0}=0.3$ Triggered',(annx,anny),color='White')
ax1.annotate('$B_{g}/B_{0}=0.3$ Untriggered',(annx,anny),color='White')
ax2.annotate('$B_{g}/B_{0}=0.1$ Triggered',(annx,anny),color='White')
ax3.annotate('$B_{g}/B_{0}=0.1$ Untriggered',(annx,anny),color='White')

#ax0.set_title('Triggered $B_{g}=0$')
#ax1.set_title('Triggered $B_{g}=0.3B_{0}$')
#ax2.set_title('Untriggered $B_{g}=0$')
#ax3.set_title('Untriggered $B_{g}=0.3B_{0}$')

plt.savefig('2x2_flds_reord.png',dpi=300,bbox_inches='tight')
plt.savefig('2x2_flds_reord.pdf',dpi=300,bbox_inches='tight')
















