
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
prtl_base = "../../tristan_acc-mec_Ez/testing_edotv/output/prtl.tot."
savename = "bguide.3_triggered"


interval=2000
tscan = 3 #go tscan output timesteps past when the particle is accelerated at tcs to assess the wpar work
#just go grab last timestep for identifying tcs's

tfinal = 65
tfinal_str = "%03d" % tfinal
namefinal = prtl_base + tfinal_str


myf_final = h5py.File(namefinal,'r')

gammae_final = myf_final['gammae']
myind_final = myf_final['inde']
myproc_final = myf_final['proce']
tcse_final = myf_final['tcse']


file_dict = {} 
for t in range(1,tfinal+1):
    tstr = "%03d" % t
    myfile = h5py.File(prtl_base + tstr,'r')
    file_dict[tstr]=myfile


#so now we have all the files open and can call h5py files from this dictionary via file_dict[tstr]

numprts = np.size(gammae_final)

wpar_list = []
finalgam_list = []
midgam_list = []
print('iterating through : ',numprts)

nonzero_tcs_prtls = np.sum(np.array(tcse_final)>0)
print('nonzero tcs prtls : ',nonzero_tcs_prtls)
count = 0
for i in range(numprts):
    #iterate through final particles
    if tcse_final[i] != 0: #only grab particles with nonzero tcs
        myt = tcse_final[i]//interval + tscan
        if myt > tfinal:#probably shouldn't be too many of these, but worth checking
            myt = tfinal
        myt_str = "%03d" % myt
        myf = file_dict[myt_str]
        myproc = myproc_final[i]
        myind = myind_final[i]

        #now need to find particle with proper index and proc
        tmpproc = np.array(myf['proce'])
        tmpind = np.array(myf['inde'])
        wpar = np.array(myf['edotve'])*2 / .45
        tmpgam = myf['gammae']
        numprt_tmp = np.size(tmpproc)
        for j in range(numprt_tmp):
            if tmpproc[j]==myproc and tmpind[j]==myind:
                #then we've found the correct particle at the time we want and we just need to grab its wpar
                wpar_list.append(wpar[j])
                midgam_list.append(tmpgam[j])
                finalgam_list.append(gammae_final[i])
                count+=1
                break
        print(count/nonzero_tcs_prtls, ' of the way through')

print('done!')
wpar_frac = np.array(wpar_list)/np.array(midgam_list)

mingam = min(finalgam_list)
maxgam = max(finalgam_list)

xarr = np.linspace(0,4000)

plt.hist2d(wpar_list,midgam_list,bins=40,norm=LogNorm())                     
plt.xlabel('$W_{||}$')                                                         
plt.ylabel('$\gamma_{final}$')                                                 
plt.plot(xarr,xarr,linestyle='dashed',color="Orange")                          
plt.colorbar(label="$N_{e}$")                                                  
plt.savefig(savename+'_wpar_midgam.png',dpi=300,bbox_inches='tight')  
plt.close()   


'''
plt.hist2d(wpar_list,finalgam_list,bins=40,norm=LogNorm())
plt.xlabel('$W_{||}$')
plt.ylabel('$\gamma_{final}$')
plt.plot(xarr,xarr,linestyle='dashed',color="Orange")
plt.colorbar(label="$N_{e}$")
plt.savefig(savename+'_wpar_gam_hist_tscan3.png',dpi=300,bbox_inches='tight')
plt.close()


plt.hist2d(wpar_frac, finalgam_list,bins=40,norm=LogNorm(),range=[[-2,2],[mingam,maxgam]])
plt.xlabel('$W_{||}/\gamma_{tcs}$')
plt.ylabel('$\gamma_{final}$')
plt.xlim(-2,2)
plt.colorbar(label="$N_{e}$")
plt.savefig(savename+'_wparfrac_tscan3.png',dpi=300,bbox_inches='tight')
plt.close()
'''

