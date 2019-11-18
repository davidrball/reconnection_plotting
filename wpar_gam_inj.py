import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})


savename = "bguide.3_triggered_injhist_offset0_highgamthresh"
prtl_base = "/home/u21/davidrball/david/tristan_acc-mec_Ez/bguide.3_wpar_zcomp/gamprt_thresh50/output/prtl.tot."
tfinal = 20 #start at early times for testing
interval = 2000.
tcs_offset = 0 #how many output timesteps to go after injection timestep to calculate the quantities

gamthresh=500 #start really high here for testing purposes
final_str = "%03d" % tfinal
namefinal = prtl_base + final_str
myf_final = h5py.File(namefinal,'r')
gammae_final = np.array(myf_final['gammae'])
initsize = np.size(gammae_final)
select = gammae_final > gamthresh
myind = np.nonzero(select)
gammae_final = gammae_final[myind]
wpar_final = np.abs(np.array(myf_final['edotve']))[myind]*2/.45
tcs_final = np.round(np.array(myf_final['tcse'])[myind]/interval) #integer division by interval
ind_final = np.array(myf_final['inde'])[myind]
proc_final = np.array(myf_final['proce'])[myind]


finalsize = np.size(gammae_final)
print('init size : ', initsize)
print('final size : ', finalsize)
gam_list = []
wpar_list = []
ind_list = []
proc_list = []
#let's see if there's enough RAM to fit all of this...  I think there should be? huge master lists containing all the relevant particle properties

for t in range(1,tfinal+1+tcs_offset):
    myt = "%03d" % t
    print(myt)
    myf = prtl_base+ myt
    myf = h5py.File(myf,'r')
    tmp_gamma = np.array(myf['gammae'])

    
    #we should be able to trim these ones too since we are only considering particles once they've been injected
    gamthresh2 = 500
    select2 =  tmp_gamma > gamthresh2
    myind2 = np.nonzero(select2)
    tmp_gamma = tmp_gamma[myind2]
    gam_list.append(tmp_gamma)
    wpar_list.append(np.array(myf['edotve'])[myind2]*2/.45)
    ind_list.append(np.array(myf['inde'])[myind2])
    proc_list.append(np.array(myf['proce'])[myind2])

print('done saving info')
#ok now our master lists contain all the info we need, first index is the output timestep, second index will grab particle

#now iterate through final particles, find their injection time and grab relevant quantities

gam_plot_list = []
wpar_plot_list = []

print('iterating through particles now')

for i in range(finalsize):
    print(float(i)/finalsize)
    #i corresponds to index of initial particles
    myproc = proc_final[i]
    myind = ind_final[i]
    tcs =int(tcs_final[i]) #tcs should be an integer corresponding to the output time closest to injection
    

    #tcs is related to the index we call in the the lists by just subtracting 1
    tcs = tcs-1 + tcs_offset
    

    #print('proc, ind', myproc, myind)

    #now we go grab the relevant quantities at the right time
    print('tcs : ',tcs)
    my_ind_list = ind_list[tcs]
    my_proc_list = proc_list[tcs]
    my_gam_list = gam_list[tcs]
    my_wpar_list = wpar_list[tcs]

    #now we need to iterate through these lists and find index that matches myproc, myind
    tmpsize = len(my_ind_list)
    #print('testing sizes, should all be the same')
    #print(tmpsize, np.size(gam_list[tcs]), np.size(wpar_list[tcs]))
    
    for j in range(tmpsize):
        if my_ind_list[j]==myind and my_proc_list[j]==myproc:
            print(my_proc_list[j], my_ind_list[j])
            print('index ' , j)

            mygam = my_gam_list[j]
            mywpar = my_wpar_list[j]
           
            gam_plot_list.append(mygam)
            wpar_plot_list.append(mywpar)
            break

plt.scatter(wpar_plot_list, gam_plot_list)
plt.xscale('log')
plt.yscale('log')

plt.savefig('testing_inj.png')

            
                    
    
