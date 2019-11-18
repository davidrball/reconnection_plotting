import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})
plt.set_cmap("viridis")

from tristan_funcs import get_val_filename



#define the simulation matrix, right now it's 3x3
#first index picks out guide field strength, 2nd is temp


time_string = "025"

#guide field 0, delgam=.00005, .0005, .005
matrix00 = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam00005/output/flds.tot."+time_string
matrix01 = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/flds.tot."+time_string
matrix02 =  "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam005/output/flds.tot."+time_string


matrix10 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam00005/bguide.1/output/flds.tot."+time_string
matrix11 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/flds.tot."+time_string
matrix12 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.1/output/flds.tot."+time_string

matrix20 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam00005/bguide.3/output/flds.tot."+time_string
matrix21 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/flds.tot."+time_string
matrix22 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.3/output/flds.tot."+time_string


matrix_namelist = [[matrix00, matrix01, matrix02], [matrix10,matrix11,matrix12], [matrix20,matrix21,matrix22]]


fig, ((ax0, ax1, ax2), (ax3, ax4, ax5), (ax6, ax7, ax8)) = plt.subplots(3,3,sharex=True,sharey=True)

axlist = [ax0, ax1, ax2, ax3, ax4, ax5, ax6 ,ax7, ax8]
betaeff_list = [.0003, .003, .03, .01, .013, .043, .09, .093, .12]

mylen = len(axlist)


for i in range(len(matrix_namelist)):
    for j in range(len(matrix_namelist[0])): #assuming that the matrix is square 
        myfile = h5py.File(matrix_namelist[i][j],'r')
        mydens = np.rot90(get_val_filename('dens',myfile)[0])

        myax_ind = 3*i + j
        myax = axlist[myax_ind]
        title = str(betaeff_list[myax_ind])
        myax.imshow(mydens/4,vmin=0,vmax=5,origin='lower')
        myax.set_adjustable('box-forced')
        myax.set_title(title)
        myax.set_yticks([])
        myax.set_xticks([])
        myax.minorticks_off()
        myfile.close()

ax0.set_ylabel('$B_{g}/B_{0}=0$')
ax3.set_ylabel('$B_{g}/B_{0}=0.1$')
ax6.set_ylabel('$B_{g}/B_{0}=0.3$')

ax6.set_xlabel('$\\theta=0.00005$')
ax7.set_xlabel('$\\theta=0.0005$')
ax8.set_xlabel('$\\theta=0.005$')
plt.savefig('testing_matrix_plot.png',dpi=300,bbox_inches='tight')

















