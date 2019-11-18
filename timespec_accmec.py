import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 8})
plt.set_cmap("viridis")
plt.rcParams['figure.figsize'] = 2.4, 8



from timespec_func import plot_timespec

plt.rcParams.update({'font.size':6})
#define the simulation matrix, right now it's 3x3
#first index picks out guide field strength, 2nd is temp




#guide field 0, delgam=.00005, .0005, .005
matrix00 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_run2/output/spect."
matrix01 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_bguide.3/output/spect."


matrix10 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_untriggered/output/spect."

matrix11 = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_untriggered_bguide.3/output/spect."

matrix_namelist = [matrix00, matrix10, matrix01,matrix11]

fig, (ax0, ax1, ax2, ax3)  = plt.subplots(4,1,sharex=True)
#fig = plt.figure()#figsize=(1.6,5.33))
#ax0 = fig.add_subplot(4,1,1,adjustable='box',aspect=.4)
#ax1 = fig.add_subplot(4,1,2,adjustable='box',aspect=.4)
#ax2 = fig.add_subplot(4,1,3,adjustable='box',aspect=.4)
#ax3 = fig.add_subplot(4,1,4,adjustable='box',aspect=.4)


#ax0 = fig.add_subplot(4,1,1,aspect=.4)
#ax1 = fig.add_subplot(4,1,2,aspect=.4)
#ax2 = fig.add_subplot(4,1,3,aspect=.4)
#ax3 = fig.add_subplot(4,1,4,aspect=.4)#


axlist = [ax0, ax1, ax2, ax3]
mylen = len(axlist)
t0 = 3
tf = 35


p=-2.5
N=2e10
xarr = np.logspace(1,4,50)
yarr = N*xarr**(p+1)



for i in range(len(matrix_namelist)):
    myname = matrix_namelist[i] 
    myax = axlist[i]
        #title = str(betaeff_list[myax_ind])
            
    plot_timespec(myax, myname, t0, tf, 4)
    myax.plot(xarr,yarr,color="Black",linestyle='--')
    myax.set_xscale('log')
    myax.set_yscale('log')
    myax.set_ylim(1e2,1e9)
    myax.set_xlim(1e0,1e4)
    myax.set_aspect(.4, adjustable='box-forced')
        #myax.set_xlim(myxlim_low,myxlim_up)
    myax.minorticks_off()
  
    myax.set_yticks([1e2,1e4,1e6,1e8])
    if i!=3:
        myax.set_xticks([])
        #myax.set_yticks([])
        #myax.set_xticks([])
        #myax.set_xlim(1e0,1e4)
        #myax.set_xticks([1e0,1e1,1e2,1e3])
        #myax.set_yticks([1e2,1e4,1e6])
        #myfile.close()

ax1.set_ylabel('$(\gamma-1) dN/d\gamma$',size=10)
ax1.yaxis.set_label_coords(-.2,-.15)
#ax3.set_ylabel('$B_{g}/B_{0}=0.1$')
#ax6.set_ylabel('$B_{g}/B_{0}=0.3$')

#ax6.set_xlabel('$\\theta=0.00005$')
#ax7.set_xlabel('$\\theta=0.0005$')
#ax8.set_xlabel('$\\theta=0.005$')
xstr = "$\gamma-1$"
ystr = "$(\gamma-1)dN/d\gamma$"
#ax0.set_ylabel(ystr)                                        
#ax2.set_ylabel(ystr)
#ax2.set_xlabel(xstr)
ax3.set_xlabel(xstr,size=10)                                      
#ax0.set_title('Triggered Bguide=0')
#ax1.set_title('Triggered Bguide=0.3')
#ax2.set_title('Untriggered Bguide=0')
#ax3.set_title('Untriggered Bguide=0.3')

plt.savefig('4x1_timespec.png',dpi=300, bbox_inches='tight')
#plt.savefig('testing_4x1.png')
















