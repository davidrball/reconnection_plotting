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



#define the simulation matrix, right now it's 3x3
#first index picks out guide field strength, 2nd is temp


time_string = "045"

#guide field 0, delgam=.00005, .0005, .005
matrix00 = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam00005/output/spect."+time_string
matrix01 = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/spect."+time_string
matrix02 =  "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam005/output/spect."+time_string
matrix03 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide0/output/spect." + time_string


matrix10 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam00005/bguide.1/output/spect."+time_string
matrix11 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/spect."+time_string
matrix12 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.1/output/spect."+time_string
matrix13  = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.1/output/spect." + time_string

#matrix13 = "pass"
#assign unfinished sims as "pass"


matrix20 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam00005/bguide.3/output/spect."+time_string
matrix21 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/spect."+time_string
matrix22 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.3/output/spect."+time_string
matrix23 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.3/output/spect." + time_string

#matrix23 = "pass"

matrix30 = "pass"

matrix31 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.7/output/spect." + time_string
matrix32 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.7/output/spect."+time_string
matrix33 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.7/output/spect."+time_string

#matrix33 = "pass"

matrix_namelist = [[matrix00, matrix01, matrix02,matrix03], [matrix10,matrix11,matrix12,matrix13], [matrix20,matrix21,matrix22,matrix23],[matrix30,matrix31,matrix32,matrix33]]


fig, ((ax0, ax1, ax2, ax3), (ax4, ax5, ax6, ax7), (ax8, ax9, ax10, ax11), (ax12, ax13, ax14, ax15)) = plt.subplots(4,4,sharex=True,sharey=True)

axlist = [ax0, ax1, ax2, ax3, ax4, ax5, ax6 ,ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15]
betaeff_list = [.0003, .003, .03, .3, .01, .013, .043, .34,.09, .093, .12, .42, .49, .493, .52, .82]

mylen = len(axlist)


for i in range(len(matrix_namelist)):
    for j in range(len(matrix_namelist[0])): #assuming that the matrix is square 
        if matrix_namelist[i][j] == "pass":
            pass
        else:
             myfile = h5py.File(matrix_namelist[i][j],'r')
             gam = np.array(get_val_filename('gamma', myfile))
             myspec = np.array(get_val_filename('rspece',myfile))[:,0,:]
             summed_spec = np.zeros(np.size(gam))
             for k in range(np.size(summed_spec)):
                 for l in range(np.shape(myspec)[1]):
                     summed_spec[k] += myspec[k][l]

             myax_ind = 4*i + j
             myax = axlist[myax_ind]
             title = str(betaeff_list[myax_ind])
             myax.plot(gam,gam*summed_spec)
             myax.set_adjustable('box-forced')
             #myax.set_title(title,size=10)
             myax.set_xscale('log')
             myax.set_yscale('log')
             myax.set_ylim(1e2,1e7)
             myax.minorticks_off()
        #myax.set_yticks([])                                                                
        #myax.set_xticks([])                                                                
             myax.set_xlim(1e0,1e4)
             myax.set_xticks([1e0,1e1,1e2,1e3])
             myax.set_yticks([1e2,1e4,1e6])
             myfile.close()


ax0.set_ylabel('$B_{g}/B_{0}=0$')
ax4.set_ylabel('$B_{g}/B_{0}=0.1$')
ax8.set_ylabel('$B_{g}/B_{0}=0.3$')
ax12.set_ylabel('$B_{g}/B_{0}=0.7$')
ax12.set_xlabel('$\\theta=0.00005$')
ax13.set_xlabel('$\\theta=0.0005$')
ax14.set_xlabel('$\\theta=0.005$')
ax15.set_xlabel('$\\theta=0.05$')
plt.savefig('testing_4x4_spec.png',dpi=300,bbox_inches='tight')

















