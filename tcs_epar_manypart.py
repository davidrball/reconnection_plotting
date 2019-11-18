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
from scipy.integrate import trapz
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 10})

#inputs
t_final=605
mysig = .3
prtl_base = "../../tristan_acc-mec_Ez/untriggered_bguide_hightimeres/output/prtl.tot."
#prtl_base = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/prtl.tot."
min_part_num = 1
max_part_num = 3
partskip = 1
interval = 200
tscan_low = 4 #output timesteps
tscan_up = 15


#derived quantities
alfven_v = np.sqrt(mysig / (1+mysig))
t_final_string = '%03d' % t_final

finalgam_list = []
epar_frac_list = []



f_prtl_num = h5py.File(prtl_base + t_final_string,'r')
prtl_num = np.size(f_prtl_num['gammae'])
print('tot prtl num : ',prtl_num)

#we can be smarter, first calculate the number of particles with nonzero ycs
ycs_nonzero = f_prtl_num['ycse']
prtl_num = np.count_nonzero(ycs_nonzero)
print('nonzero ycs part num : ',prtl_num)


f_prtl_num.close()
#gammae_final = np.array(f_prtl_final['gammae'])
#gammae = np.array(f_prtl['gammae'])
#inde = f_prtl_final['inde']
#proce = f_prtl_final['proce']
#tcse = np.array(f_prtl_final['tcse'])//interval
#max_E_ind = np.argmax(gammae_final)


#f_prtl_final.close()


my_nummax = prtl_num

logmin = np.log10(min_part_num)
#print('logmin : ', logmin)
logmax = np.log10(my_nummax-10)
logsteps = 20


#ok let's add sortin functionality, just presort the gamma list here, then indices will be in order of lorentz factor and we don't have to go through and zero out shit for every particle


partnum_list = []
#for part_num in range(min_part_num,max_part_num,partskip):
for part_num in np.logspace(logmin,logmax,num=logsteps):
    my_part_num = int(part_num)
    print(partnum_list)
    while my_part_num in partnum_list: #something here isn't working
        print('logstep too small, adding 1 to partnum : ')
        my_part_num +=1
    partnum_list.append(my_part_num)

    prtl_filename = prtl_base + t_final_string
    print('part num: ',my_part_num)
    f_prtl = h5py.File(prtl_filename, "r")

    #identify highest energy particle at end of run, save its proc and ind
    gammae = np.array(f_prtl['gammae'])
    inde = f_prtl['inde']
    proce = f_prtl['proce']
    tcse = np.array(f_prtl['tcse'])//interval 
    max_E_ind = np.argmax(gammae)


    #this step seems to take a long time, maybe try to figure out a better way to do it
    #maybe we should just presort the array
    if my_part_num > 1:
        for i in range(1,my_part_num):
        
            gammae[max_E_ind] = 0 #zero out all gamma entries above part_num's
            max_E_ind = np.argmax(gammae)
    print('finished zeroing process')
    max_E_ind = np.argmax(gammae)
    finalgam = np.max(gammae)
    finalgam_list.append(finalgam)
#if you want to track particles lower down the energy list
#gammae[max_E_ind] = 0
    #print(max_E_ind)

#max_E_ind = np.argmax(gammae) # 2nd most energetic particle
    max_E_proc_val = proce[max_E_ind]
    max_E_ind_val = inde[max_E_ind]


    tcs = tcse[max_E_ind]
    print('tcs : ', tcs)
    tlow = int(tcs - tscan_low)
    tup = int(tcs + tscan_up)



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
    #f_prtl.close()


    #for an individual particle, need to select time window around tcs
    for t in range(tlow,tup):
        t_list.append(t)
        print(t)
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
                #print(inde[tmp_ind]) # should match above printed index
        
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
                f_prtl.close()
                break
    

    gamma_arr =np.array(gamma_list)
    vx_arr = np.array(vx_list)
    vy_arr = np.array(vy_list)
    vz_arr = np.array(vz_list)
    bx_arr = np.array(bx_list)
    by_arr = np.array(by_list)
    bz_arr = np.array(bz_list)
    ex_arr = np.array(ex_list)
    ey_arr = np.array(ey_list)
    ez_arr = np.array(ez_list)

    vx_arr /= gamma_arr
    vy_arr /= gamma_arr
    vz_arr /= gamma_arr


    edotv_x = -vx_arr*ex_arr
    edotv_y = -vy_arr*ey_arr
    edotv_z = -vz_arr*ez_arr

    edotv_tot = edotv_x + edotv_y + edotv_z


    bmag = np.sqrt(bx_arr**2 + by_arr**2 + bz_arr**2)
    unitbx = bx_arr / bmag
    unitby = by_arr / bmag
    unitbz = bz_arr / bmag

    upar = vx_arr*unitbx + vy_arr*unitby + vz_arr*unitbz
    epar = ex_arr*unitbx + ey_arr*unitby + ez_arr*unitbz
    epar_upar = -upar*epar
    eperp_uperp = edotv_tot - epar_upar

    edotv_par_cumulative = []
    edotv_perp_cumulative = []
    edotv_tot_cumulative = []

    for i in range(np.size(eperp_uperp)):
        tmp_par = trapz(epar_upar[0:i])
        tmp_perp = trapz(eperp_uperp[0:i])
        tmp_tot = trapz(edotv_tot[0:i])
        edotv_par_cumulative.append(tmp_par)
        edotv_tot_cumulative.append(tmp_tot)
        edotv_perp_cumulative.append(tmp_perp)

    
    epar_frac = edotv_par_cumulative[-1]/edotv_tot_cumulative[-1]
    print('fraction of total work done by parallel E field : ',  epar_frac)

    epar_frac_list.append(epar_frac)

plt.scatter(epar_frac_list,finalgam_list,color="Blue")
plt.xlabel('$W_{||}/W_{tot}$')
plt.ylabel('$\gamma_{final}$')
plt.xlim(0,1)
plt.yscale('log')
plt.savefig('wpar_gam_scattertest.png',dpi=300,bbox_inches='tight')


    
    #ax.plot(t_list,gamma_list-gamma_list[0],color="Black",label='$\gamma$')
    
