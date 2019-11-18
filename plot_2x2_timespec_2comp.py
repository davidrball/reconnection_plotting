import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy import special
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 12})
plt.set_cmap("viridis")

from timespec_func import plot_timespec, plot_timespec_label


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




#matrix_namelist = [[matrix00, matrix10], [matrix01,matrix11]]
matrix_namelist = [[matrix01,matrix11],[matrix00,matrix10]]


fig, ((ax0, ax1), (ax2, ax3))  = plt.subplots(2,2,sharex=True,sharey=True)

axlist = [ax0, ax1, ax2, ax3]
mylen = len(axlist)
t0 = 15
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


yarr_list = [yarr3, yarr4, yarr1, yarr2]

def MJ(N,theta,gamma):
    beta = np.sqrt(1-1/gamma**2)
    bess = special.kn(2,1/theta)
    tmp = gamma**2 * beta * np.exp(-gamma/theta) / (theta * bess)
    return N*tmp

def xpoint(N,cutoff,gamma,falloff):
    return (N*gamma**-1)*np.exp(-(gamma/cutoff)**falloff)

Nxpoint_untrig = 8e5
Nxpoint_trig = Nxpoint_untrig /3.

mycut = 450

xpoint_spec_untrig = xpoint(Nxpoint_untrig, mycut, xarr,1)
xpoint_spec_trig = xpoint(Nxpoint_trig,mycut,xarr,1)

NMJ_untrig = 1.7e7
theta_untrig=17
MJ_untrig_spec = MJ(NMJ_untrig,theta_untrig,xarr)

NMJ_trig = 1.7e7
theta_trig=14
MJ_trig_spec = MJ(NMJ_trig, theta_trig,xarr)

for i in range(len(matrix_namelist)):
    for j in range(len(matrix_namelist[0])): #assuming that the matrix is square 
        myname = matrix_namelist[i][j]
        if myname == "pass":
            pass
        else:
            myax_ind = 2*i + j
            myax = axlist[myax_ind]
        #title = str(betaeff_list[myax_ind])

            if i==1 and j==1:
                plot_timespec_label(myax,myname,t0,tf,5,2000,3.)
                ax3.legend(loc='lower left',frameon=False,prop={'size':6},ncol=2)
            else:
                plot_timespec(myax, myname, t0, tf, 5)
            myax.plot(xarr,yarr_list[myax_ind],color="Black",linestyle='--')
            myax.set_xscale('log')
            myax.set_yscale('log')
            myax.set_ylim(1e2,2e8)
            myax.set_xlim(1e0,1e4)
            #myax.set_xticks([1,1e1,1e2,1e3])
        #myax.set_xlim(myxlim_low,myxlim_up)
            myax.minorticks_off()

            if myax_ind ==0:
                myax.plot(xarr, xarr*xpoint_spec_trig,color="Green",linestyle='dashed',linewidth=.5,label='$dN/d\gamma = N\gamma^{-1}e^{-\gamma/\gamma_{c}}$')
                myax.plot(xarr, xarr*MJ_trig_spec,color="Orange",linestyle='dashed',linewidth=.5,label='$dN/d\gamma=$ MJ($\\theta=17$)')
                myax.legend(loc='lower left', frameon=False,prop={'size':8})
            elif myax_ind ==1:
                myax.plot(xarr,xarr*xpoint_spec_untrig,color="Green",linestyle='dashed',linewidth=.5,label='$dN/d\gamma = 3N\gamma^{-1}e^{-\gamma/\gamma_{c}}$')
                myax.plot(xarr,xarr*MJ_untrig_spec,color="Orange",linestyle='dashed',linewidth=.5,label='$dN/d\gamma=$ MJ($\\theta=14$)')
                myax.legend(loc='lower left',frameon=False,prop={'size':8})
        #myax.set_yticks([])
        #myax.set_xticks([])
        #myax.set_xlim(1e0,1e4)
        myax.set_xticks([1e0,1e1,1e2,1e3])
        #myax.set_yticks([1e2,1e4,1e6])
        #myfile.close()



#ax0.set_ylabel('$B_{g}/B_{0}=0$')
#ax3.set_ylabel('$B_{g}/B_{0}=0.1$')
#ax6.set_ylabel('$B_{g}/B_{0}=0.3$')

#ax6.set_xlabel('$\\theta=0.00005$')
#ax7.set_xlabel('$\\theta=0.0005$')
#ax8.set_xlabel('$\\theta=0.005$')

ax3.annotate('Time ($\omega_{p}^{-1}$)',(2.5,.8e5),size=10)


plt.subplots_adjust(wspace=.1,hspace=.1)
xstr = "$\gamma-1$"
ystr = "$(\gamma-1)dN/d\gamma$"
ax0.set_ylabel(ystr)                                        
ax2.set_ylabel(ystr)
ax2.set_xlabel(xstr)
ax3.set_xlabel(xstr)                                      
annx = 6e2
anny = 1e7

ax0.annotate('Triggered \n $B_{g}=0.3B_{0}$',(annx,anny),size=10)
ax1.annotate('Untriggered \n $B_{g}=0.3B_{0}$',(annx,anny),size=10)
ax2.annotate('Triggered \n $B_{g}=0.1B_{0}$',(annx,anny),size=10)
ax3.annotate('Untriggered \n $B_{g}=0.1B_{0}$',(annx,anny),size=10)

plt.savefig('2x2_timespec_2comp_reord.png',dpi=300,bbox_inches='tight')
plt.savefig('2x2_timespec_2comp_reord.pdf',dpi=300,bbox_inches='tight')
















