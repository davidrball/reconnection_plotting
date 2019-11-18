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
max_part_num = 10000
part_jump = 1000
interval = 200
tscan_low = 4 #output timesteps
tscan_up = 15


#derived quantities
alfven_v = np.sqrt(mysig / (1+mysig))
t_final_string = '%03d' % t_final

for part_num in range(min_part_num,max_part_num,part_jump):
    prtl_filename = prtl_base + t_final_string
    print(prtl_filename)
    f_prtl = h5py.File(prtl_filename, "r")

    #identify highest energy particle at end of run, save its proc and ind
    gammae = np.array(f_prtl['gammae'])
    inde = f_prtl['inde']
    proce = f_prtl['proce']
    tcse = np.array(f_prtl['tcse'])//interval 
    max_E_ind = np.argmax(gammae)

    if part_num > 1:
        for i in range(1,part_num):
        
            gammae[max_E_ind] = 0 
            max_E_ind = np.argmax(gammae)

    max_E_ind = np.argmax(gammae)
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
    f_prtl.close()


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
                

                break
    gamma_arr = np.array(gamma_list)
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
        #tmp_par = trapz(epar_upar[0:i])
        #tmp_perp = trapz(eperp_uperp[0:i])
        #tmp_tot = trapz(edotv_tot[0:i])
        tmp_par = np.sum(epar_upar[0:i])
        tmp_perp = np.sum(eperp_uperp[0:i])
        tmp_tot = np.sum(edotv_tot[0:i])

        edotv_par_cumulative.append(tmp_par)
        edotv_tot_cumulative.append(tmp_tot)
        edotv_perp_cumulative.append(tmp_perp)

    mynorm = (gamma_list[-1]-gamma_list[0]) / edotv_tot_cumulative[-1]
    print('mynorm : ',mynorm)

    fig, ax = plt.subplots()

    epar_frac = edotv_par_cumulative[-1]/edotv_tot_cumulative[-1]
    print('fraction of total work done by parallel E field : ',  epar_frac)

    
    #ax.plot(t_list,gamma_list-gamma_list[0],color="Black",label='$\gamma$')
    
    #plotting cumulative 
    ax.plot(t_list,mynorm*np.array(edotv_par_cumulative),color="Red",label="$\int v_{||}E_{||}dt$")
    ax.plot(t_list,mynorm*np.array(edotv_perp_cumulative),color="Blue",label="$\int v_{\\bot}E_{\\bot}dt$")

    ax.plot(t_list,mynorm*(np.array(edotv_par_cumulative)+np.array(edotv_perp_cumulative)),color="Black",linestyle='--')
    ax.set_ylabel('$\Delta \gamma$')
    ax.set_xlabel('Time (output timesteps)')
    ax.legend(loc='upper left')

    #ax.set_yscale('log')
    #ax2 = ax.twinx()
    #ax2.plot(t_list,epar_upar,color="Blue")
    #ax2.plot(t_list,eperp_uperp,color="Red")
    #ax2.plot(t_list,edotv_tot,color="Black")
    plt.savefig('tcs_zoom/untriggered_bguide/'+str(part_num)+'.png',dpi=300,bbox_inches='tight')
