
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
t_final=130

t_string_final = '%03d' % t_final
c_omp = 3 
istep = 12

tskip = 1

#prtl_base = "../../tristan_acc-mec_Ez/bguide.3_wpar_zcomp/output/prtl.tot."
prtl_base = "../../tristan_acc-mec_Ez/testing_edotv/bguide.3_triggered_alwaystrack_hightimeres/output/prtl.tot."


prt_step = 200
min_part_num = 1
max_part_num = 1000


col_list = []
count = 0

prt_list = []

for part_num in range(min_part_num, max_part_num+1, prt_step):
    prt_list.append(part_num)

#putting in our prt_list by hand
#prt_list = [3,5,7,9,11,100,1000]


prt_list = [1,10,20,100,1000,10000]
#use 1 and 15 for PCS and MCS prtls


#number 5 is a good example of merger prtl

#prt_list = [1,100,10001]


for part_num in prt_list:
    col_string = "C"+str(count)
    col_list.append(col_string)
    count += 1


print('len col list : ', col_list)

mycount = 0
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
    print('highest gamma we have:')
    print(gammae[max_E_ind])

#max_E_ind = np.argmax(gammae) # 2nd most energetic particle
    max_E_proc_val = proce[max_E_ind]
    max_E_ind_val = inde[max_E_ind]



    #print("index , proc")
    #print(max_E_ind_val, max_E_proc_val)

    gamma_list = []
    t_list = []
    wpar_list = []
    gamma0=0
    for t in range(1,t_final,tskip): #in some cases the high E particle is injected later in time, so it'
        print(t)
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
                #print(inde[tmp_ind]) # should match above printed index
                
                gamma = np.array(f_prtl['gammae'])[tmp_ind]
                if gamma0==0:
                    gamma0=gamma
                gamma_list.append(gamma-gamma0)
                print('gamma0 is : ',gamma0)
                proc = np.array(f_prtl['proce'])[tmp_ind]
                
            
                tmpwpar = np.array(f_prtl['edotve'])[tmp_ind]*2/.45
                wpar_list.append(tmpwpar)
                break
            
            elif i==np.size(tmp)-1:
                tmp_ind = i
                print("particle not injected yet")
                
                gamma_list.append(0)
                wpar_list.append(0)

    #tarr = np.array(t_list)*6/10. #in units of 100 w_p ^-1            
    tarr = np.array(t_list)*(500*.15)/100.
    print('mycount : ', mycount)
    #if mycount == 0: #only put label on first prtls so we don't get lots of labels
    #    plt.plot(tarr, gamma_list,color=col_list[mycount],label='$\gamma$')
    #    plt.plot(tarr,wpar_list,color=col_list[mycount],label="$W_{||,z}$",linestyle='--')

    plt.plot(tarr, gamma_list,color=col_list[mycount])
    plt.plot(tarr,wpar_list,color=col_list[mycount],linestyle='--')
    mycount += 1

mysig = .3
sige = mysig*1836
sige /= 2.

yarr = sige*np.ones(np.size(tarr))



plt.plot(0,0,color="Black",label='$\Delta \gamma$')
plt.plot(0,0,color="Black",label='$W_{||,z}$',linestyle='--')
plt.plot(tarr,yarr,color="Magenta",label='$\sigma_{e}/2$')




plt.yscale('log')
plt.xlim(0,tarr[-1])
plt.ylim(1e-1,1e4)
plt.xlabel('Time ($100 \; \omega_{p}^{-1}$)')
plt.ylabel('$\Delta \gamma$')
plt.legend(loc="upper left",frameon=False,prop={'size':12})
#plt.savefig('manypart_wparz'+str(prt_step)+'.png',dpi=300,bbox_inches='tight')
#plt.savefig('manypart_wparz'+str(prt_step)+'.pdf',dpi=300,bbox_inches='tight')

#plt.savefig('manypart_wparz_handselect3.png',dpi=300,bbox_inches='tight')         
#plt.savefig('manypart_wparz_handselect3.pdf',dpi=300,bbox_inches='tight') 
#plt.savefig('manypart_wparz_deltagamma_log.png',dpi=300,bbox_inches='tight')
#plt.savefig('manypart_wparz_deltagamma_log.pdf',dpi=300,bbox_inches='tight')
#plt.savefig('testing_moreprts.png',dpi=300,bbox_inches='tight')
plt.savefig('prtltrack_test.png',dpi=300,bbox_inches='tight')
plt.savefig('prtltrack_test.pdf',dpi=300,bbox_inches='tight')
