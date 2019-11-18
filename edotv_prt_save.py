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
t_final=605
mysig = .3
prtl_base = "../../tristan_acc-mec_Ez/untriggered_bguide_hightimeres/output/prtl.tot."
#prtl_base = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/prtl.tot."
min_part_num = 14
max_part_num = 20000

partspace = 100

#derived quantities
alfven_v = np.sqrt(mysig / (1+mysig))
t_final_string = '%03d' % t_final

for part_num in range(min_part_num,max_part_num,partspace):
    prtl_filename = prtl_base + t_final_string
    print(prtl_filename)
    f_prtl = h5py.File(prtl_filename, "r")

    #identify highest energy particle at end of run, save its proc and ind
    gammae = np.array(f_prtl['gammae'])
    inde = f_prtl['inde']
    proce = f_prtl['proce']
    max_E_ind = np.argmax(gammae)

    if part_num > 1:
        for i in range(1,part_num):
        
            gammae[max_E_ind] = 0 
            max_E_ind = np.argmax(gammae)

    max_E_ind = np.argmax(gammae)
#if you want to track particles lower down the energy list
#gammae[max_E_ind] = 0
    print(max_E_ind)

#max_E_ind = np.argmax(gammae) # 2nd most energetic particle
    max_E_proc_val = proce[max_E_ind]
    max_E_ind_val = inde[max_E_ind]

    target = open('edotv_prtls/bguide.3_untriggered_hightimeres_latetime/'+str(part_num) + '.txt','w')
    target.write('proc, ind, t , u , v , w , ex , ey , ez , bx , by , b\
z , gamma')
    target.write('\n')

    print("index , proc")
    print(max_E_ind_val, max_E_proc_val)
    target.write(str(max_E_ind_val))
    target.write('\n')
    target.write(str(max_E_proc_val))
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
        tmp_proc = np.abs(proce - max_E_proc_val)
        tmp = np.abs(inde - max_E_ind_val)
    
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
