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
                
                
#t = 49
#for sigma=.3, one alfven crossing time corresponds to 24 output timesteps, let's go a bit under this to be conservative
t_final=605 #corresponds to .8 of an alfven crossing time

mysig = .3
alfven_v = np.sqrt(mysig / (1+mysig))



t_string_final = '%03d' % t_final
c_omp = 3 
istep = 12
#fld_base = "../../tristan_acc-mec_Ez/triggered_bguide_hightimeres_begintrack/output/flds.tot."
#prtl_base ="../../tristan_acc-mec_Ez/triggered_bguide_hightimeres_begintrack/output/prtl.tot."

fld_base = "../../tristan_acc-mec_Ez/untriggered_bguide_hightimeres/output/flds.tot."
prtl_base = "../../tristan_acc-mec_Ez/untriggered_bguide_hightimeres/output/prtl.tot."



min_part_num = 100
max_part_num = 110
for part_num in range(min_part_num,max_part_num):
    

    
    #fileName = bguide1_base + t_string
    fld_filename = fld_base + t_string_final
    prtl_filename = prtl_base + t_string_final

    print(fld_filename)
    print(prtl_filename)
    f_fld = h5py.File(fld_filename,  "r")
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

    plt.set_cmap("viridis")
    gamma_list = []
    t_list = []
    u_list = []
    v_list = []
    w_list = []
    #t_max = 350
    
    EoverB_list = []
    for t in range(20,610,10): #in some cases the high E particle is injected later in time, so it'
        t_list.append(t)    
        t_string = '%03d' % t
        

        fld_filename = fld_base + t_string
        prtl_filename = prtl_base + t_string

        f_fld = h5py.File(fld_filename,  "r")
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
                gamma_list.append(gamma-1)
                x_pos = get_val("xe",f_prtl)[tmp_ind]
                y_pos = get_val("ye",f_prtl)[tmp_ind]
                proc = get_val("proce",f_prtl)[tmp_ind]
            #tmpu = get_val("ue",f_prtl)[tmp_ind]
            #tmpv = get_val("ve",f_prtl)[tmp_ind]
            #tmpw = get_val("we", f_prtl)[tmp_ind]
            #u_list.append(np.abs(tmpu))
            #v_list.append(np.abs(tmpv))
            #w_list.append(np.abs(tmpw))
            
                bx = get_val('bxe',f_prtl)[tmp_ind]
                by = get_val('bye',f_prtl)[tmp_ind]

                ez = get_val('eze',f_prtl)[tmp_ind]
            
                EoverB_list.append(ez/(np.sqrt(bx**2+by**2)*alfven_v))

            #print("proc")
            #print(proc)
            #print("x, y")
            #print(x_pos, y_pos) #I think in units of computational cells?
                
            
                intx = int(x_pos)
                inty = int(y_pos)
            
                intx /= istep
                inty /= istep        
            
            #print("on our grid:")
            #print(intx,inty)
            #fld information
                dens = get_val_fld("dens", f_fld)
            #plt.imshow(dens,origin='lower')
            #fig,ax1 = plt.subplots(1)
                ax1 = plt.subplot2grid((1,2),(0,0),rowspan=1,colspan=1)
                ax2 = plt.subplot2grid((1,2),(0,1),colspan=1)
                xhlf = np.shape(dens)[1]/2        
                xscan = 100
                xlow = int(xhlf - xscan)
                xup = int(xhlf + xscan)
                xext = istep*2*xscan / c_omp # in electron skin depths
                yext = istep * np.shape(dens)[0] /c_omp
                ax1.imshow(dens[:,xlow:xup],origin='lower',vmax=20, extent=[0,xext,0,yext])
                ax1.set_xlabel('x $(c/\omega_{p})$',size=14)
                ax1.set_ylabel('y $(c/\omega_{p})$',size=14)
                
                    
                ax2.plot(t_list,gamma_list,color="Red",label='$\gamma$')
                ax2.scatter(t,gamma-1,color="Red")
                ax2.set_xlim(0,t_final)
                ax2.set_ylim(1e-1,1e5)
                ax2.set_yscale('log')
                interval=1000
                
                ax2.set_xlabel('Time $(30\omega_{p}^{-1})$')
                ax2.set_ylabel('$\gamma-1$')
                ax2.legend(loc='upper right',frameon=False)

                ax3 = ax2.twinx()
            #ax3.plot(t_list,u_list,label="u="+str(tmpu)[0:4],color="Blue",linewidth=.5)
            #ax3.plot(t_list,v_list,label="v="+str(tmpv)[0:4],color="Green",linewidth=.5)
            #ax3.plot(t_list,w_list,label="w="+str(tmpw)[0:4],color="Cyan",linewidth=.5)
                ax3.plot(t_list, EoverB_list, color = "Blue",linewidth=.5,label='$E_{z}/v_{A}B_{xy}$')
                ax3.legend(loc='lower right',frameon=False)
                ax3.set_ylim(-10,10)
                ax3.set_ylabel('$E_{z}/(v_{A} B_{xy})$')
                ax2.set_xlim(0,t_final)
            #ax3.set_yscale('log')
            #ax3.set_xlim(1,t_max)
            #ax3.set_ylim(1,1e4)
            #plt.legend(loc='upper left',prop={'size':6})
            #have to reindex the position of particle since we've shifted coordinates by narrowing window
                intx -= xlow
            
            #now need to convert position of particle to skin depths:
                intx *= istep/c_omp            
                inty *= istep/c_omp
            #circle = patches.Circle((intx,inty),radius=gamma/50 + 5,fill=True,color='red') #for having radius grow
                circle = patches.Circle((intx,inty),radius=20,fill=False,color='red')
                plt.tight_layout()
                ax1.add_patch(circle)
            #plt.savefig("/home/dball/tristan_out/real_runs/movie_plots/particle_tracking/untriggered/betap1/E"+str(int(part_num)) + "_particle/"+t_string+".png",bbox_inches='tight',dpi=400)
                tstr_save = str(t)
                plt.savefig('hightimeres_prts/untriggered_bguide/E'+ str(part_num)+'_' +tstr_save+'.png',bbox_inches='tight',dpi=300)
                plt.close()
                f_fld.close()
                f_prtl.close()
                break
            
            elif i==np.size(tmp)-1:
                tmp_ind = i
                print("particle not injected yet")
                gamma = get_val("gammae",f_prtl)[tmp_ind]
                gamma_list.append(gamma)
                x_pos = get_val("xe",f_prtl)[tmp_ind]
                y_pos = get_val("ye",f_prtl)[tmp_ind]
                proc = get_val("proce",f_prtl)[tmp_ind]
            #tmpu = get_val("ue",f_prtl)[tmp_ind]                                                                                                                                        
            #tmpv = get_val("ve",f_prtl)[tmp_ind]                                                                                                                                        
            #tmpw = get_val("we", f_prtl)[tmp_ind]                                                                                                                                       
            #u_list.append(np.abs(tmpu))                                                                                                                                                 
            #v_list.append(np.abs(tmpv))                                                                                                                                                 
            #w_list.append(np.abs(tmpw))                                                                                                                                                 

                bx = get_val('bxe',f_prtl)[tmp_ind]
                by = get_val('bye',f_prtl)[tmp_ind]

                ez = get_val('eze',f_prtl)[tmp_ind]

                EoverB_list.append(ez/(alfven_v*np.sqrt(bx**2+by**2)))
