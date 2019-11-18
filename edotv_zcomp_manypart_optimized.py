
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
t_final=260

t_string_final = '%03d' % t_final
c_omp = 3 
istep = 12

tskip = 20

#prtl_base = "../../tristan_acc-mec_Ez/bguide.3_wpar_zcomp/output/prtl.tot."
prtl_base = "../../tristan_acc-mec_Ez/testing_edotv/bguide.3_triggered_alwaystrack_hightimeres/output/prtl.tot."


col_list = []


#prt_list = [3,57,50010]

prt_list = [3,57, 99000]

#use 1 and 15 for PCS and MCS prtls

numprts = len(prt_list)

ind_list = []
proc_list = []
gam0_list = []

for i in range(numprts):
    col_string = "C"+str(i)
    col_list.append(col_string)



print('len col list : ', col_list)

mycount = 0
t_list = []
gamma_list = [] #will have sublists corresponding to each particle
wpar_list = []
for part_num in prt_list:
    
    #fileName = bguide1_base + t_string
    prtl_filename = prtl_base + t_string_final
    #print(prtl_filename)
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
    print('gammas')
    print(gammae[max_E_ind])

#max_E_ind = np.argmax(gammae) # 2nd most energetic particle
    max_E_proc_val = proce[max_E_ind]
    max_E_ind_val = inde[max_E_ind]
    ind_list.append(max_E_ind_val)
    proc_list.append(max_E_proc_val)

    gam0_list.append([])
    gamma_list.append([])
    wpar_list.append([])
    

for t in range(1,t_final,tskip): #in some cases the high E particle is injected later in time, so it'
    print(t)
    t_list.append(t)    
    t_string = '%03d' % t
    prtl_filename = prtl_base + t_string
    f_prtl = h5py.File(prtl_filename, "r")
    inde = np.array(f_prtl['inde'])
    proce = np.array(f_prtl['proce'])
    
    

        
    for i in range(np.size(inde)):
        #print(i)
        tmpcount = 0
        for j in range(numprts): 
            if inde[i] == ind_list[j] and proce[i]==proc_list[j]:
                
                tmp_ind = i
                #print(inde[tmp_ind]) # should match above printed index
                
                gamma = np.array(f_prtl['gammae'])[tmp_ind]
                tmpwpar = np.array(f_prtl['edotve'])[tmp_ind]*2/.45

                if len(gam0_list[j])==0:
                    gam0_list[j].append(gamma)
                    
                
                gamma_list[j].append(gamma)
                wpar_list[j].append(tmpwpar)
                tmpcount += 1 #prtl found
                if tmpcount == numprts:
                    break #found all the particles at this timestep, now stop

            
    #if prtl isnt injected, need to add a 0 to its list
    for j in range(numprts):
        if len(gamma_list[j]) == len(t_list)-1:
            gamma_list[j].append(0)
            wpar_list[j].append(0)
    f_prtl.close()

    #tarr = np.array(t_list)*6/10. #in units of 100 w_p ^-1            
tarr = np.array(t_list)*(500*.15)/100.



print('gam0 list')
print(gam0_list)
for j in range(numprts):
    
    my_gam0 = gam0_list[j]
    gam0_arr = np.ones(np.size(tarr))*my_gam0
    plt.plot(tarr, np.array(gamma_list[j])-gam0_arr,color=col_list[j])
    plt.plot(tarr, wpar_list[j],color=col_list[j],linestyle='--')



mysig = .3
sige = mysig*1836
sige /= 2.

yarr = sige*np.ones(np.size(tarr))



plt.plot(0,0,color="Black",label='$\Delta \gamma$')
plt.plot(0,0,color="Black",label='$W_{||,z}$',linestyle='--')
plt.plot(tarr,yarr,color="Magenta",label='$\sigma_{e}/2$')




plt.yscale('log')
plt.xlim(0,tarr[-1])
plt.ylim(1e0,1e4)
plt.xlabel('Time ($100 \; \omega_{p}^{-1}$)')
plt.ylabel('$\Delta \gamma$')
plt.legend(loc="upper left",frameon=False,prop={'size':12})

plt.savefig('prtltrack_test_opt_morelow.png',dpi=300,bbox_inches='tight')
plt.savefig('prtltrack_test_opt_morelow.pdf',dpi=300,bbox_inches='tight')
