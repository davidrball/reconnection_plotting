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
matrix00 = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam00005/output/spect."
matrix01 = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/spect."
matrix02 =  "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam005/output/spect."


matrix10 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam00005/bguide.1/output/spect."
matrix11 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/spect."
matrix12 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.1/output/spect."

matrix20 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam00005/bguide.3/output/spect."
matrix21 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/spect."
matrix22 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.3/output/spect."


matrix_namelist = [[matrix00, matrix01, matrix02], [matrix10,matrix11,matrix12], [matrix20,matrix21,matrix22]]


fig, ((ax0, ax1, ax2), (ax3, ax4, ax5), (ax6, ax7, ax8)) = plt.subplots(3,3,sharex=True,sharey=True)

axlist = [ax0, ax1, ax2, ax3, ax4, ax5, ax6 ,ax7, ax8]
betaeff_list = [.0003, .003, .03, .01, .013, .043, .09, .093, .12]

mylen = len(axlist)
t0 = 5
tf = 50

for i in range(len(matrix_namelist)):
    for j in range(len(matrix_namelist[0])): #assuming that the matrix is square 
        myname = matrix_namelist[i][j]
        
        myax_ind = 3*i + j
        myax = axlist[myax_ind]
        #title = str(betaeff_list[myax_ind])
        
        plot_timespec(myax, myname, t0, tf, 5)
        myxlim_low = 10**j * .1
        myxlim_up = 10**j * 1e4
        #print(myxlim_low, myxlim_up)
        #myax.set_title(title,size=10)
        myax.set_xscale('log')
        myax.set_yscale('log')
        myax.set_ylim(1e2,1e9)
        myax.set_xlim(myxlim_low,myxlim_up)
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
ax3.set_ylabel(ystr)                                      
ax6.set_ylabel(ystr)                                      
ax6.set_xlabel(xstr)                                      
ax7.set_xlabel(xstr)                                       
ax8.set_xlabel(xstr)  


plt.savefig('testing_matrix_plot_timespec.png',dpi=300,bbox_inches='tight')

















