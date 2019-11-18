import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 8})
plt.set_cmap("viridis")

from timespec_func import plot_timespec


#define the simulation matrix, right now it's 3x3
#first index picks out guide field strength, 2nd is temp




#guide field 0, delgam=.00005, .0005, .005
'''
matrix00 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_bguide0_triggered_stride1_thresh2/output/spect."
matrix01 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/spect."


matrix10 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_bguide0_untriggered_stride1_thresh2/output/spect."

matrix11 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_untriggered_bguide.3/output/spect."
'''

#accmec runs
matrix00 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/bguide.1_triggered_stride1_gamthresh50/output/spect."

matrix01 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/spect."
                   
matrix10 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/bguide.1_untriggered_stride1_gamthresh50/output/spect."
                   
matrix11 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/bguide.3_untriggered_stride1_gamthresh50/output/spect."




matrix_namelist = [[matrix00, matrix10], [matrix01,matrix11]]

fig, ((ax0, ax1), (ax2, ax3))  = plt.subplots(2,2,sharex=True,sharey=True)

axlist = [ax0, ax1, ax2, ax3]
mylen = len(axlist)
t0 = 5
tf = 65


p=-2.75
N1=6.5e10
N2=9e10
N3=10e10
N4=9.5e10
xarr = np.logspace(1,4,50)
yarr1 = N1*xarr**(p+1)
yarr2 = N2*xarr**(p+1)
yarr3 = N3*xarr**(p+1)
yarr4 = N4*xarr**(p+1)


yarr_list = [yarr1, yarr2, yarr3, yarr4]


for i in range(len(matrix_namelist)):
    for j in range(len(matrix_namelist[0])): #assuming that the matrix is square 
        myname = matrix_namelist[i][j]
        if myname == "pass":
            pass
        else:
            myax_ind = 2*i + j
            myax = axlist[myax_ind]
        #title = str(betaeff_list[myax_ind])
            if i==0 and j==0:
                plot_timespec_withlabel(myax,myname,t0,tf,5,2000,3.)
            else:
                plot_timespec(myax, myname, t0, tf, 5)
            myax.plot(xarr,yarr_list[myax_ind],color="Black",linestyle='--')
            myax.set_xscale('log')
            myax.set_yscale('log')
            myax.set_ylim(1e4,1e9)
            myax.set_xlim(1e1,1e4)
        #myax.set_xlim(myxlim_low,myxlim_up)
            myax.minorticks_off()
        #myax.set_yticks([])
        #myax.set_xticks([])
        #myax.set_xlim(1e0,1e4)
        #myax.set_xticks([1e0,1e1,1e2,1e3])
        #myax.set_yticks([1e2,1e4,1e6])
        #myfile.close()

#ax0.set_ylabel('$B_{g}/B_{0}=0$')
#ax3.set_ylabel('$B_{g}/B_{0}=0.1$')
#ax6.set_ylabel('$B_{g}/B_{0}=0.3$')

#ax6.set_xlabel('$\\theta=0.00005$')
#ax7.set_xlabel('$\\theta=0.0005$')
#ax8.set_xlabel('$\\theta=0.005$')
xstr = "$\gamma-1$"
ystr = "$(\gamma-1)dN/d\gamma$"
ax0.set_ylabel(ystr)                                        
ax2.set_ylabel(ystr)
ax2.set_xlabel(xstr)
ax3.set_xlabel(xstr)                                      
annx = 1e3
anny = 1e8

ax0.annotate('Triggered \n $Bg=0.1B_{0}$',(annx,anny),size=10)
ax1.annotate('Untriggered \n $Bg=0.1B_{0}$',(annx,anny),size=10)
ax2.annotate('Triggered \n $Bg=0.3B_{0}$',(annx,anny),size=10)
ax3.annotate('Untriggered \n $Bg=0.3B_{0}$',(annx,anny),size=10)

plt.savefig('2x2_timespec.png',dpi=300,bbox_inches='tight')
plt.savefig('2x2_timespec.pdf',dpi=300,bbox_inches='tight')
















