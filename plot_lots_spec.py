import matplotlib
matplotlib.use('Agg') #for elgato
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from scipy.special import kn 
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 18})

from tristan_funcs import get_val_filename, decompose_spectra, return_rspec, return_spece_withcut

prebase = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/"

#reference case
spec0 = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/spect."

#triggered guide fields

spec0 = "../../tristan_acc-mec_Ez/8k_bguide0_triggered_stride1_thresh2/output/spect."
spec1 = "../../tristan_acc-mec_Ez/8k_bguide.1_allprts/output/spect."
spec2 = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/spect."
spec_base_list = [spec0,spec1,spec2]
#spec_base_list = [spec0, spec1, spec2, spec3]
#L_list = [16000,16000,16000,16000]
L_list = [8000,8000,8000]
bguide_list = [0,.1,.3]#,.7,1]
color_list = ["C1", "C2", "C3"]#, "C4"]

listlen = len(spec_base_list)

t_start = 40
t_end = 41

for t in range(t_start,t_end):
    mybase_list = []
    for specbase in spec_base_list:
        mybase_list.append(specbase)
    
    t_str = str(t)

    if len(t_str) == 1: 
        for i in range(len(mybase_list)):
            mybase_list[i] += "00"                                                         
    if len(t_str) == 2:
        for i in range(len(mybase_list)):
            mybase_list[i] += "0"                                                           
    for i in range(len(mybase_list)):
        mybase_list[i] += t_str

    file_list = []
    for mybase in mybase_list:
        print(mybase)
        file_list.append(h5py.File(mybase,'r'))
    
    #if you want to just input your own files
    #myf1 = prebase + "/tristan-mp_reconnection/16k_triggered_finals/sig1/delgam.2/output/spect.045"
    #myf2 = prebase + "/tristan-mp_reconnection/outflow/sig1/delgam.2/output/spect.090"
    #file_list = [h5py.File(myf1,"r"),h5py.File(myf2,"r")]

    

    mygam_list = []
    myrspece_list = []

    gam, ion, lec = return_rspec(file_list)
    totgam, totion, totlec = return_spece_withcut(file_list,L_list)
     
    #print(len(totgam),len(totion),len(totlec))
    #print(len(gam), len(lec), len(ion))


    
    for i in range(len(lec)):
        mygam, myrspece = decompose_spectra(totlec[i],lec[i],gam[i])
        mygam_list.append(mygam)
        myrspece_list.append(myrspece)
    
    '''
    for i in range(len(lec)):
        plt.plot(totgam[i],totgam[i]*totlec[i], color=color_list[i], linewidth=0.5)
        plt.plot(gam[i],gam[i]*lec[i],color=color_list[i],linewidth=2,label='$B_{g}/B_{0}=$'+str(bguide_list[i]))
    '''
    for i in range(len(lec)):
        my_totgam = totgam[i]
        my_rspecgam = mygam_list[i]
        my_totlec = totlec[i]
        my_rspec = myrspece_list[i]
        plt.plot(my_totgam, my_totgam*my_totlec, color=color_list[i],linewidth=.5)
        plt.plot(my_rspecgam, my_rspecgam*my_rspec,color=color_list[i],linewidth=2,label='$B_{g}/B_{0}=$'+str(bguide_list[i]))

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\gamma-1$')
    plt.ylabel('$(\gamma-1)dN/d\gamma$')
    plt.ylim(1e2,1e9)
    plt.xlim(1e0,1e4)
    #plt.title('$\sigma=1 \; \; \\beta=.15$')
    #plt.title('$\sigma=0.3 \; \; \\beta=0.006$')
    #plt.legend(loc='upper left', prop={'size':8})
    plt.legend(frameon=False,prop={'size':14})
    
    plt.savefig('bg_spec_compare.pdf',dpi=300,bbox_inches='tight')
    #plt.savefig('plots/spectra_movies/sig.3_delgam0005_triggered/'+t_str+'.png',dpi=300,bbox_inches='tight')
    plt.close()
    for file in file_list:
        file.close()
