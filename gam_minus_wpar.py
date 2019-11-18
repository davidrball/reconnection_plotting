
import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

#inputs
#prtl_base = "../../tristan_acc-mec_Ez/testing_edotv/untriggered/output/prtl.tot."
prtl_base = "../../tristan_acc-mec_Ez/testing_edotv/bguide.1/output/prtl.tot."
savename = "bguide.1_triggered"
tfinal = 65
gamthresh=50 #don't include cold thermals


tfinal_str = "%03d" % tfinal
namefinal = prtl_base + tfinal_str
myf_final = h5py.File(namefinal,'r')
gammae_final = np.array(myf_final['gammae'])

initsize = np.size(gammae_final)


select = gammae_final > gamthresh
myind = np.nonzero(select)
gammae_final = gammae_final[myind]
wpar_final = np.abs(np.array(myf_final['edotve']))[myind]*2/.45
#not many negativs, just take abs to make us able to plot in log, take this away as needed
finalsize = np.size(gammae_final)

print('init size : ', initsize)
print('final size : ', finalsize)

mingam = np.min(gammae_final)
maxgam = np.max(gammae_final)


gam_min_wpar = gammae_final - wpar_final

#plt.scatter(gammae_final,gam_min_wpar)

#plt.xlim(200,1e4)
#plt.ylim(200,1e4)
#plt.xscale('log')
#plt.yscale('log')
#plt.savefig('gam_min_wpar_test.png')




plt.hist2d(np.log10(gammae_final),np.log10(gam_min_wpar),bins=40,norm=LogNorm())
#plt.hist2d(wpar_frac, np.log10(gammae_final),bins=40,norm=LogNorm()),range=[[0,1],[np.log10(mingam),np.log10(maxgam)]])
plt.xlabel('$\gamma$')
plt.ylabel('$\gamma - W_{||}$')
#plt.xlim(0,1)
plt.colorbar(label="$N_{e}$")
plt.savefig('testing_gam_min.png',dpi=300,bbox_inches='tight')
plt.close()

