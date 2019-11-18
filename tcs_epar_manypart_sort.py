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
#prtl_base = "../../tristan_acc-mec_Ez/untriggered_bguide_hightimeres/output/prtl.tot."
#prtl_base = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/prtl.tot."
prtl_base = "../../tristan_acc-mec_Ez/bguide.3_triggered_hightimeres_nothresh/prtl.tot."

min_part_num = 1
#max_part_num = 20000
partskip = 50
interval = 200
tscan_low = 4 #output timesteps
tscan_up = 15

#derived quantities
alfven_v = np.sqrt(mysig / (1+mysig))
t_final_string = '%03d' % t_final
final_file_string = prtl_base + t_final_string
finalgam_list = []
epar_frac_list = []
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
for part_num in range(min_part_num,max_part_num,partskip):
for part_num in np.logspace(logmin,logmax,num=logsteps):
    my_part_num = int(part_num)
    print(partnum_list)
    while my_part_num in partnum_list: #something here isn't working
        print('logstep too small, adding 1 to partnum : ')
        my_part_num +=1
    partnum_list.append(my_part_num)
    prtl_filename = prtl_base + t_final_string
    print('part num: ',my_part_num)
    
    #because we've already sorted, we can now just read off the index and processor

    finalgam = gammae_sorted[my_part_num]
    my_inde = inde_sorted[my_part_num]
    my_proce = proce_sorted[my_part_num]
    my_tcs = tcse_sorted[my_part_num]
    print('tcs : ', my_tcs)
    print('gamma : ',finalgam)
        
    

    tlow = int(my_tcs - tscan_low)
    tup = int(my_tcs + tscan_up)

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
    #for an individual particle, need to select time window around tcs
    for t in range(tlow,tup):
        t_list.append(t)
        print(t)
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
    finalgam_list.append(finalgam)
plt.scatter(epar_frac_list,finalgam_list,color="Blue",s=.5)
plt.xlabel('$W_{||}/W_{tot}$')
plt.ylabel('$\gamma_{final}$')
plt.xlim(0,1)
plt.yscale('log')
plt.savefig('wpar_gam_sort_100point_test_logspace_triggered.png',dpi=300,bbox_inches='tight')


    
    #ax.plot(t_list,gamma_list-gamma_list[0],color="Black",label='$\gamma$')
    
