import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tristan_funcs import get_val_filename
import h5py
import numpy as np
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 15})

myfld_base = "../../tristan_acc-mec_Ez/8k_untriggered/output/flds.tot."

prtl_final = "../../tristan_acc-mec_Ez/8k_untriggered/output/prtl.tot.035"
prtl_base = "../../tristan_acc-mec_Ez/8k_untriggered/output/prtl.tot."

f_prtl_final = h5py.File(prtl_final, 'r')



tcse = np.array(get_val_filename('tcse',f_prtl_final))
ycse = np.array(get_val_filename('ycse',f_prtl_final))
gammae = np.array(get_val_filename('gammae',f_prtl_final))
ez_bxy = np.array(get_val_filename('ezbxye',f_prtl_final))

prtnum = np.size(tcse)

gam_max_final = np.max(gammae)

f_prtl_final.close()
final_len = np.size(tcse)

sig=0.3
va = np.sqrt(sig/(sig+1))
istep = 12
c_omp = 3
interval = 2000.

t_start = 35
t_stop = 36




fld_scan = 120

yext = fld_scan * istep / (1000*c_omp)
tot_gam_list = []
tot_ycs_list = []
tot_tcs_list = []
for t in range(t_start,t_stop):
    t_string = str(t)
    if len(t_string) == 2:
        fld_base =myfld_base+ "0"
      #  myprtl = prtl_base + "0"
    elif len(t_string)==1:  
        fld_base =myfld_base+ "00"
     #   myprtl = prtl_base + "00"
    fld_base += t_string
    #myprtl += t_string

    f_fld = h5py.File(fld_base,'r')
    #f_prtl = h5py.File(myprtl,'r')
    dens = np.rot90(get_val_filename('dens',f_fld)[0,:,:])
    ymid = np.shape(dens)[0]/2
    ylow = int(ymid - fld_scan)
    yup = int(ymid + fld_scan)
    dens = dens[ylow:yup,:]
    #dens = np.fliplr(dens)
    newmid = np.shape(dens)[0]/2
    xext = np.shape(dens)[1]*istep / (1000*c_omp)
    #fig = plt.figure()
    fig, (ax,ax2) = plt.subplots(2,1,sharex=True)
    #ax = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
    #ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3,sharex=ax)

    im = ax.imshow(dens/4.,origin='lower',cmap='viridis',vmin=1,vmax=6, extent= [-xext/2, xext/2, -yext, yext])
    ax2.set_xlabel('$x \; (1000 \; c/\omega_{p})$')
    ax.set_ylabel('$y \; (1000 \; c/\omega_{p})$')
    ax2.set_ylabel('$\gamma$')
    ax2.set_xlim(-xext/2, xext/2)
    #fig.subplots_adjust(hspace=0)
    cbar_ax = fig.add_axes([.125, .935, .77, .03])
    cb = fig.colorbar(im,cax=cbar_ax, orientation="horizontal")
    cbar_ax.set_xlabel("Density (particles per cell / $N_{ppc})$",fontsize=14,labelpad=-45)
    #fld info plotted
    
    #now loop through particles, for any with tcs/interval between t-1 and t, add plot them over the fld quantities

    f_fld.close()
    count = 0
    gam_list = []
    ycs_list = []
    tcs_list = []
    ez_bxy_list = []
    x_list = []
    tlow = t-.5
    tup = t+.5
    print(tlow, tup)
    for i in range(prtnum):
        mycs = ycse[i] #in units of cells
        mytcs = tcse[i]
        #mygam = gammae[i]
        if mycs != 0 and mytcs != 1:
            #in actual computational steps
            mytcs /= interval
            
            #print(tlow, mytcs, tup)
            if tlow < mytcs <= tup:
                #myx = xe[i]/istep

                count +=1
                #mycs /= istep
                mycs /= (c_omp * 1000)
                mycs -= xext/2
                ycs_list.append(mycs)
                tcs_list.append(mytcs/t_stop)
                #print(mycs)
                #x_list.append(myx-ylow)              
                my_ezbxy = np.abs(ez_bxy[i])/va
                #print(my_ezbxy)
                mygam = gammae[i]
                loggam = np.log10(mygam)
                gam_list.append(mygam)
                ez_bxy_list.append(my_ezbxy)
                #ax.scatter(mycs, -.25,s=((mygam/200.)**2.5), c=my_ezbxy,vmin=0,vmax =4 , cmap='YlOrRd')
    #loggam_arr = np.log10(np.array(gam_list))
    size_arr = (np.array(gam_list)/200)**2.5


    im2 = ax.scatter(ycs_list, -.25*np.ones(np.size(ycs_list)), s=size_arr, c=ez_bxy_list,vmin=0,vmax=4,cmap='YlOrRd')

    tot_gam_list.append(gam_list)
    tot_ycs_list.append(ycs_list)
    tot_tcs_list.append(tcs_list)
    
    for index in range(len(tot_gam_list)):
        sub_gam_list = tot_gam_list[index]
        sub_ycs_list = tot_ycs_list[index]
        sub_tcs_list = tot_tcs_list[index]
        ax2.scatter(sub_ycs_list,sub_gam_list, c=sub_tcs_list,cmap='plasma',marker='o',s=1,vmin=0,vmax=1)
    ax2.set_yscale('log')
    ax2.set_ylim(5e1,5e3)
    cbar_ax2 = fig.add_axes([.92, .53, .02, .35])
    cb2 = fig.colorbar(im2,cax=cbar_ax2, orientation="vertical")
    cbar_ax2.set_ylabel("$E_{z}/B_{xy}$",fontsize=14)
    cbar_ax2.yaxis.set_label_position('right')
    cb2.set_ticks([0,1,2,3,4])
    #plt.scatter(ycs_list, x_list, s=1, c=loggam_arr,vmin=0,vmax=np.log10(np.max(gam_final)),cmap='viridis')
    
    if count == 0:
        pass
    else:
        print('particles plotted : ',count)
        print('max gam : ', max(gam_list))
        print('max of all particles : ', gam_max_final)
 
    
    #plt.savefig('fld_acc_plots/untriggered/'+t_string+'.png',dpi=300,bbox_inches='tight')
    plt.savefig('testing_fld_prt_acc.png',dpi=300,bbox_inches='tight')
    plt.close()
