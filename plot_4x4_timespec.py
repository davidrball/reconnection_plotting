import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 10})
plt.set_cmap("viridis")

from tristan_funcs import get_val_filename
from timespec_func import plot_timespec


#define the simulation matrix, right now it's 3x3
#first index picks out guide field strength, 2nd is temp




#guide field 0, delgam=.00005, .0005, .005
matrix00 = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam00005/output/spect."
matrix01 = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/spect."
matrix02 =  "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam005/output/spect."
matrix03 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide0/output/spect."


matrix10 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam00005/bguide.1/output/spect."
matrix11 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/spect."
matrix12 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.1/output/spect."
matrix13  = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.1/output/spect." 

#matrix13 = "pass"
#assign unfinished sims as "pass"


matrix20 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam00005/bguide.3/output/spect."
matrix21 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/spect."
matrix22 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.3/output/spect."
matrix23 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.3/output/spect." 

#matrix23 = "pass"

matrix30 = "pass"
matrix31 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.7/output/spect." 
matrix32 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.7/output/spect."

matrix33 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.7/output/spect."



matrix_namelist = [[matrix00, matrix01, matrix02,matrix03], [matrix10,matrix11,matrix12,matrix13], [matrix20,matrix21,matrix22,matrix23],[matrix30,matrix31,matrix32,matrix33]]



fig, ((ax0, ax1, ax2, ax3), (ax4, ax5, ax6, ax7), (ax8, ax9, ax10, ax11), (ax12, ax13, ax14, ax15)) = plt.subplots(4,4,sharex=True,sharey=True)

axlist = [ax0, ax1, ax2, ax3, ax4, ax5, ax6 ,ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15]
betaeff_list = [.0003, .003, .03, .3, .01, .013, .043, .34,.09, .093, .12, .42, .49, .493, .52, .82]

mylen = len(axlist)

t0=10
tf = 45

xarr = np.logspace(1.5,4,10)
p_list = [-1.9,-1.8,-3.5,-6.2]
N_list = [1.5e11,1e11,3e15,8e25]
yarr_list = []

for i in range(len(p_list)):
    myp = p_list[i]
    myN = N_list[i]
    yarr = myN*xarr**myp
    yarr_list.append(yarr)


for i in range(len(matrix_namelist)):
    for j in range(len(matrix_namelist[0])): #assuming that the matrix is square 
        if matrix_namelist[i][j] == "pass":
            pass
        else:
            myname = matrix_namelist[i][j]
            myax_ind = 4*i + j
            myax = axlist[myax_ind]

            plot_timespec(myax, myname, t0, tf, 5)
            #myxlim_low = 10**j * .1
            #myxlim_up = 10**j * 1e4
            #print(myxlim_low, myxlim_up)                                                 
            #myax.set_title(title,size=10)                                                
            myax.set_xscale('log')
            myax.set_yscale('log')
            myax.set_ylim(1e2,1e9)
            myax.set_xlim(1e1,1e4)
            #myax.set_xlim(myxlim_low,myxlim_up)
            myax.minorticks_off()
           
    
            myax.plot(xarr,yarr_list[j],linestyle='--',color="Black")
          
xstr = "$\gamma-1$"
ystr = "$(\gamma-1)dN/d\gamma$"
ax0.set_ylabel(ystr)
ax4.set_ylabel(ystr)
ax8.set_ylabel(ystr)
ax12.set_ylabel(ystr)
ax12.set_xlabel(xstr)
ax13.set_xlabel(xstr)
ax14.set_xlabel(xstr)
ax15.set_xlabel(xstr)
plt.savefig('testing_4x4_timespec.png',dpi=300,bbox_inches='tight')

















