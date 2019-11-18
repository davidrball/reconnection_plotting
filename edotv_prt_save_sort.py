# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 11:50:18 2017

@author: dball
"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 10})

#inputs

#inputs                                                                                     
t_final=605
mysig = .3
prtl_base = "../../tristan_acc-mec_Ez/untriggered_bguide_hightimeres/output/prtl.tot."     
#prtl_base = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/prtl.to\t."                                                                                         
#prtl_base = "../../tristan_acc-mec_Ez/bguide.3_triggered_hightimeres_nothresh/prtl.tot."
min_part_num = 1000
#max_part_num = 1000 #currently using nummax
partskip = 100
interval = 200

#derived quantities                                                                         
t_final_string = '%03d' % t_final
final_file_string = prtl_base + t_final_string
f_prtl_num = h5py.File(final_file_string,'r')

gammae_final = np.array(f_prtl_num['gammae'])
print('sorting final arrays')
sort_arr = gammae_final.argsort()
gammae_sorted = gammae_final[sort_arr][::-1]
#print(gammae_sorted[0:20])                                                                 
inde_sorted = np.array(f_prtl_num['inde'])[sort_arr][::-1]
proce_sorted = np.array(f_prtl_num['proce'])[sort_arr][::-1]
tcse_sorted = np.array(f_prtl_num['tcse'])[sort_arr][::-1]//interval
print('done sorting')

ycs_nonzero = f_prtl_num['ycse']
print('total particles : ', np.size(ycs_nonzero))
prtl_num = np.count_nonzero(ycs_nonzero)
print('nonzero ycs part num : ',prtl_num)
f_prtl_num.close()

my_nummax = prtl_num/4.

logmin = np.log10(min_part_num)
#print('logmin : ', logmin)                                                                 
logmax = np.log10(my_nummax)
logsteps = 100
partnum_list = []
#for part_num in range(min_part_num,max_part_num,partskip):                                 
#for part_num in np.logspace(logmin,logmax,num=logsteps):
for part_num in range(min_part_num,my_nummax,partskip):
    my_part_num = int(part_num)
    while my_part_num in partnum_list:
        print('logstep too small, adding 1 to partnum')
        my_part_num += 1
    partnum_list.append(my_part_num)


    prtl_filename = prtl_base + t_final_string
    print(prtl_filename)
    f_prtl = h5py.File(prtl_filename, "r")

    finalgam = gammae_sorted[my_part_num]
    my_inde = inde_sorted[my_part_num]
    my_proce = proce_sorted[my_part_num]
    my_tcs = tcse_sorted[my_part_num]
    print('tcs : ', my_tcs)
    print('gamma : ',finalgam)



    target = open('edotv_prtls/bguide.3_untriggered_hightimeres_latetime/withtcs/'+str(part_num) + '.txt','w')
    target.write('proc, ind, tcs,t , u , v , w , ex , ey , ez , bx , by , b\
z , gamma')
    target.write('\n')

    print("index , proc")
    #print(max_E_ind_val, max_E_proc_val)
    target.write(str(my_inde))
    target.write('\n')
    target.write(str(my_proce))
    target.write('\n')
    target.write(str(my_tcs))
    target.write('\n')
    gamma_list = []
    t_list = []
    vx_list = []
    vy_list = []
    vz_list = []
    ex_list = []
    ey_list = []
    ez_list = []
    bx_list = []
    by_list = []
    bz_list = []
    f_prtl.close()
    for t in range(1,t_final):
        t_list.append(t)    
        t_string = '%03d' % t
        prtl_filename = prtl_base + t_string
        f_prtl = h5py.File(prtl_filename, "r")
        inde = f_prtl['inde']
        proce = f_prtl['proce']
        tmp_proc = np.abs(proce - my_proce)
        tmp = np.abs(inde - my_inde)
    
        for i in range(np.size(tmp)):
        #print(i)
            if tmp[i] == 0 and tmp_proc[i]==0:
                
                tmp_ind = i
                print(inde[tmp_ind]) # should match above printed index
        
                gamma = f_prtl['gammae'][tmp_ind]
                gamma_list.append(gamma)
                
                vx = f_prtl["ue"][tmp_ind]
                vy = f_prtl["ve"][tmp_ind]
                vz = f_prtl["we"][tmp_ind]
                vx_list.append(vx)
                vy_list.append(vy)
                vz_list.append(vz)
                bx = f_prtl['bxe'][tmp_ind]
                by = f_prtl['bye'][tmp_ind]
                bz = f_prtl['bze'][tmp_ind]
                bx_list.append(bx)
                by_list.append(by)
                bz_list.append(bz)
                ex = f_prtl['exe'][tmp_ind]
                ey = f_prtl['eye'][tmp_ind]
                ez = f_prtl['eze'][tmp_ind]
                ex_list.append(ex)
                ey_list.append(ey)
                ez_list.append(ez)
                


            #took out exception handling for particles not being saved at the timestep, may need to put this back in
    #target.write(str(ind) + ' , ' + str(proc))
    #target.write('\n')
    target.write(str(t_list))
    target.write('\n')
    target.write(str(vx_list))
    target.write('\n')
    target.write(str(vy_list))
    target.write('\n')
    target.write(str(vz_list))
    target.write('\n')
    target.write(str(ex_list))
    target.write('\n')
    target.write(str(ey_list))
    target.write('\n')
    target.write(str(ez_list))
    target.write('\n')
    target.write(str(bx_list))
    target.write('\n')
    target.write(str(by_list))
    target.write('\n')
    target.write(str(bz_list))
    target.write('\n')
    target.write(str(gamma_list))
    target.close()
