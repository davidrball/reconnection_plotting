import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tristan_funcs import get_val_filename, list_types
import h5py
import numpy as np
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 15})


#myprtl = "../../tristan_acc-mec_Ez/8k_untriggered_bguide.3/output/prtl.tot.035"
#myprtl = "../../tristan_acc-mec_Ez/test2/output/prtl.tot.004"
#myfld = "../../tristan_acc-mec/untrig_test/output/flds.tot.004"
#myprtl = "../../tristan_acc-mec/untrig_test/output/prtl.tot.003"
#myprtl = "../../tristan_acc-mec_Ez/8k_bguide.3_allprts/output/prtl.tot.001"

myprtl = "../../tristan_acc-mec_Ez/testing_edotv/bguide.1_untriggered/output/prtl.tot.065"

myprtl = h5py.File(myprtl,'r')

gammae = np.array(get_val_filename('gammae',myprtl))

tcse = np.array(get_val_filename('tcse',myprtl))
ycse = np.array(get_val_filename('ycse',myprtl))
ez_bxye = np.array(get_val_filename('ezbxye',myprtl))


sigma= 0.3
va = np.sqrt(sigma / (1+sigma))

#count number of prtls w/ tcse = 1
mycount =0
#for i in range(np.size(gammae)):
#    if tcse[i] == 1:
#        mycount +=1
#print('number of tcse=1 : ', mycount)

tmax = np.max(tcse)
print('tmax : ' , tmax)
logmax = np.log10(tmax)

prtnum = np.size(tcse)

#print('prtnum : ', prtnum)
istep = 12.
c_omp = 3.
my0 = 8160
#ycse *= istep / c_omp
col_array = tcse


xext = my0 / c_omp

#plotting ez_bxy vs y
sc = plt.scatter(ycse, ez_bxye/va,c=col_array,cmap='magma', s = .2) 
plt.ylim(-5,5)
plt.xlim(0, xext)
plt.colorbar(sc,label='$t_{cs}/t_{max}$')
plt.ylabel('$|E_{z}|/v_{A}B_{xy}$')
plt.xlabel('$y_{cs} \; (c/\omega_{p})$')
plt.savefig('8k_bguide.1_untriggered_ycs_ez_allprts.png',dpi=300,bbox_inches='tight')
plt.close()

sc = plt.scatter(ez_bxye/va, gammae, c=col_array,cmap='magma', s=.2)
plt.colorbar(sc,label='$t_{cs}/t_{max}$')
plt.xlim(-5,5)
plt.yscale('log')
plt.xlabel('$|E_{z}|/v_{A}B_{xy}$')
plt.ylabel('$\gamma$')
plt.savefig('8k_bguide.1_untriggered_ez_gam_allprts.png',bbox_inches='tight',dpi=300)
plt.close()
 


ho_arr = np.linspace(0,xext,10)
vert_arr = np.ones(10)*1836*sigma / 2.
sc = plt.scatter(ycse, gammae, c=col_array, cmap='magma',s=.2)     
plt.xlabel('$y_{cs} \; (c/\omega_{p})$')
plt.xlim(0,xext)
plt.plot(ho_arr, vert_arr,color="Black",linestyle='--',label='$\sigma_{e}/2.$')
plt.legend(frameon=False)
plt.colorbar(sc,label='$t_{cs}/t_{max}$')
#plt.title('Triggered Bguide=0.3')
plt.ylabel('$\gamma$')
plt.yscale('log')
plt.savefig('8k_bguide.1_untriggered_ycs_gam_allprts.png',bbox_inches='tight',dpi=300)     
plt.close()
myprtl.close()

