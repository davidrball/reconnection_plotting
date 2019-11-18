# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 11:50:18 2017

@author: dball
"""

# -*- coding: utf-8 -*-a
"""
Created on Thu Jul 28 11:21:48 2016

@author: davidrball
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})



                
                
#t = 49
#for sigma=.3, one alfven crossing time corresponds to 24 output timesteps, let's go a bit under this to be conservative
t_final=200

t_string_final = '%03d' % t_final
c_omp = 3 
istep = 12



prtl_base = "../../tristan_acc-mec_Ez/bguide.3_wpar_zcomp/output/prtl.tot."



min_part_num = 1
max_part_num = 2
for part_num in range(min_part_num,max_part_num):
    #fileName = bguide1_base + t_string
    prtl_filename = prtl_base + t_string_final
    print(prtl_filename)
    f_prtl = h5py.File(prtl_filename, "r")

    #identify highest energy particle at end of run, save its proc and ind
    gammae = np.array(f_prtl['gammae'])
    inde = np.array(f_prtl['inde'])
    proce = np.array(f_prtl['proce'])
    max_E_ind = np.argmax(gammae)

    if part_num > 1:
        for i in range(1,part_num):
        
            gammae[max_E_ind] = 0
            max_E_ind = np.argmax(gammae)

    max_E_ind = np.argmax(gammae)
#if you want to track particles lower down the energy list
#gammae[max_E_ind] = 0
    print('highest gamma we have:')
    print(gammae[max_E_ind])

#max_E_ind = np.argmax(gammae) # 2nd most energetic particle
    max_E_proc_val = proce[max_E_ind]
    max_E_ind_val = inde[max_E_ind]



    print("index , proc")
    print(max_E_ind_val, max_E_proc_val)

    gamma_list = []
    t_list = []
    wpar_list = []
    for t in range(1,t_final,1): #in some cases the high E particle is injected later in time, so it'
        t_list.append(t)    
        t_string = '%03d' % t
        prtl_filename = prtl_base + t_string
        f_prtl = h5py.File(prtl_filename, "r")
        inde = np.array(f_prtl['inde'])
        proce = np.array(f_prtl['proce'])
        tmp_proc = np.abs(proce - max_E_proc_val)
        tmp = np.abs(inde - max_E_ind_val)
    
    #beginning of new stuff to handle processor index as well
    
        for i in range(np.size(tmp)):
        #print(i)
            if tmp[i] == 0 and tmp_proc[i]==0:
                
                tmp_ind = i
                print(inde[tmp_ind]) # should match above printed index
        
                gamma = np.array(f_prtl['gammae'])[tmp_ind]
                gamma_list.append(gamma-1)
                
                proc = np.array(f_prtl['proce'])[tmp_ind]
                
            
                tmpwpar = np.array(f_prtl['edotve'])[tmp_ind]*2/.45
                wpar_list.append(tmpwpar)
                break
            
            elif i==np.size(tmp)-1:
                tmp_ind = i
                print("particle not injected yet")
                
                gamma_list.append(0)
                wpar_list.append(0)
    tarr = np.array(t_list)*6/10. #in units of 100 w_p ^-1            
    plt.plot(tarr, gamma_list,color="C2",label='$\gamma$')
    plt.plot(tarr,wpar_list,color="C5",label="$W_{||,z}$",linestyle='--')
    plt.yscale('log')
    plt.ylim(1e-1,1e4)
    plt.xlabel('Time ($100 \; \omega_{p}^{-1}$)')
    plt.ylabel('$\gamma$')
    plt.legend(loc="upper left",frameon=False)
    plt.savefig('prtls_zcomp/E'+str(part_num)+'_wparz_gam.pdf',dpi=300,bbox_inches='tight')
                
