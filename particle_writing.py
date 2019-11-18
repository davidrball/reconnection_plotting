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
plt.rcParams.update({'font.size': 10})
def get_val(input_string, file):
    str_list = list(file.keys())
    my_string = input_string
    if my_string in str_list:   
        my_list = np.array(file.get(my_string))
        return my_list
    else:
        print("argument not found, your options are:")
        print(list(file.keys()))
        return 0
             
def get_val_fld(input_string, file):
    str_list = list(file.keys())
    my_string = input_string
    if my_string in str_list:   
        #my_list = np.array(list(file.get(my_string))[0])
        my_list = np.array(file.get(my_string))[0]
        return my_list
    else:
        print("argument not found, your options are:")
        print(list(file.keys()))
        return 0
                
                
def plot_val(input_string, file):
    my_array = get_val(input_string, file)
    if type(my_array) == np.ndarray:
        #plt.imshow(np.log(my_array[:,:]), origin='lower')
        plt.figure(figsize=(12,8))
        #plt.imshow(np.log10((np.rot90(my_array))), origin='lower' ,vmin=0,vmax=2)
        plt.imshow((np.rot90(my_array))[:,:], origin='lower')
        plt.colorbar()
        print(my_array.shape)
    else:
        pass
                    
def plotx_val(input_string, file,my_scan):
    my_array = get_val(input_string,file)
    if type(my_array) == np.ndarray:
        shape = my_array.shape
        xlen = shape[1]
        ylen = shape[0]
                            
        xhlf = int(xlen/2.0)
        scan = my_scan
        xup = xhlf+scan
        xlow = xhlf-scan
                            
         #plt.imshow(np.log(my_array[:,:]), origin='lower')
        #plt.figure(figsize=(6.5*2,4*2))
        #plt.imshow(np.log10((np.rot90(np.abs(my_array[:,xlow:xup])))), origin='lower' ''',vmin=0,vmax=70''', vmin=0, vmax = 2)
       
        plt.imshow((np.rot90(my_array[:,xlow:xup])), origin='lower' ,vmin=0,vmax=20)        
        #plt.colorbar(label="Log(Dens)")
        #plt.xscale("log")
        #print(my_array.shape)
    else:
        pass

#t = 49
#for sigma=.3, one alfven crossing time corresponds to 24 output timesteps, let's go a bit under this to be conservative
t_final=25 #corresponds to .8 of an alfven crossing time

t_string = str(t_final)
c_omp = 3 
istep = 12



min_part_num = 4
max_part_num = 1000
step = 10
for part_num in range(min_part_num,max_part_num, step):
    #opening the file we are writing to named after particle energy number
    target = open('particle_txt/sig.3_delgam0005/bguide.1_earlytime/'+str(part_num)+'.txt','w')
    target.write('proc, ind, t , x , y , u , v , w , ex , ey , ez , bx , by , bz , gamma')
    target.write('\n')
    prtl_base = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/prtl.tot."

    if len(t_string) == 2:
        prtl_base += "0"
    elif len(t_string)==1:
        prtl_base += "00"
    
    prtl_filename = prtl_base + t_string
    print(prtl_filename)
    f_prtl = h5py.File(prtl_filename, "r")

    #identify highest energy particle at end of run, save its proc and ind
    gammae = get_val("gammae",f_prtl)
    inde = get_val("inde", f_prtl)
    proce = get_val("proce",f_prtl)
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
    

    print("index , proc")
    print(max_E_ind_val, max_E_proc_val)
    gamma_list = []
    x_list = []
    y_list = []
    u_list = []
    v_list = []
    w_list = []
    t_list = []
    ex_list = []
    ey_list = []
    ez_list = []
    bx_list = []
    by_list = []
    bz_list = []

    for t in range(1,t_final): #in some cases the high E particle is injected later in time, so it'
        t_list.append(t)    
        t_string = str(t)
        prtl_base ="../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/prtl.tot."


        if len(t_string) == 2:
            prtl_base += "0"
        elif len(t_string)==1:
            prtl_base += "00"
    
        prtl_filename = prtl_base + t_string
        f_prtl = h5py.File(prtl_filename, "r")
        inde = get_val("inde",f_prtl)
        proce = get_val("proce", f_prtl)
        tmp_proc = np.abs(proce - max_E_proc_val)
        tmp = np.abs(inde - max_E_ind_val)
    
    #beginning of new stuff to handle processor index as well
    
        for i in range(np.size(tmp)):
        #print(i)
            if tmp[i] == 0 and tmp_proc[i]==0:
                
                tmp_ind = i
                print(inde[tmp_ind]) # should match above printed index
        
                

                gamma = get_val("gammae",f_prtl)[tmp_ind]
                gamma_list.append(gamma)
                x_pos = get_val("xe",f_prtl)[tmp_ind]
                x_list.append(x_pos)
                y_pos = get_val("ye",f_prtl)[tmp_ind]
                y_list.append(y_pos)
                proc = get_val("proce",f_prtl)[tmp_ind]
                ind = get_val("inde",f_prtl)[tmp_ind]
                tmpu = get_val("ue",f_prtl)[tmp_ind]
                tmpv = get_val("ve",f_prtl)[tmp_ind]
                tmpw = get_val("we", f_prtl)[tmp_ind]
                u_list.append(np.abs(tmpu))
                v_list.append(np.abs(tmpv))
                w_list.append(np.abs(tmpw))
                bx = get_val('bxe',f_prtl)[tmp_ind]
                by = get_val('bye',f_prtl)[tmp_ind]
                bz = get_val('bze',f_prtl)[tmp_ind]
                ex = get_val('exe',f_prtl)[tmp_ind]
                ey = get_val('eye',f_prtl)[tmp_ind]
                ez = get_val('eze',f_prtl)[tmp_ind]
                bx_list.append(bx)
                by_list.append(by)
                bz_list.append(bz)
                ex_list.append(ex)
                ey_list.append(ey)
                ez_list.append(ez)
                f_prtl.close()
                break
            
            elif i==np.size(tmp)-1:
                tmp_ind = i
                print("particle not injected yet")
    
    target.write(str(ind) + ' , ' + str(proc))
    target.write('\n')
    target.write(str(t_list))
    target.write('\n')
    target.write(str(x_list))
    target.write('\n')
    target.write(str(y_list))
    target.write('\n')
    target.write(str(u_list))
    target.write('\n')
    target.write(str(v_list))
    target.write('\n')
    target.write(str(w_list))
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
