import matplotlib
matplotlib.use('Agg') #for elgato
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from scipy.special import kn 
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 10})

from tristan_funcs import get_val_filename, decompose_spectra, return_rspec, return_spece_withcut

prebase = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/"

#reference case
spec0 = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/spect."

#triggered guide fields
spec1 = prebase+"bguide.3/output/spect."
spec2 = prebase + "bguide.7/output/spect."
spec3 = prebase + "bguide1/output/spect."

spec_base_list = [spec0, spec1, spec2, spec3]
L_list = [16000,16000,16000,16000]




sigma = 0.3
bguide_list = [0,.3,.7,1]
color_list = ["C1", "C2", "C3", "C4"]

bguide_array = np.array(bguide_list)
sigma_g_array = (bguide_array**2)*sigma

time_multiplier = np.sqrt((sigma+1+sigma_g_array)/(sigma+1))




listlen = len(spec_base_list)

t_start = 5
t_end = 70

for t in range(t_start,t_end):
    base_time = t
    base_time_string = str(t)
    mybase_list = []
    for specbase in spec_base_list:
        mybase_list.append(specbase)
    for i in range(len(mybase_list)):
        my_time_multiplier = time_multiplier[i]
        mytime = int(t*my_time_multiplier)
        t_str = str(mytime)

        if len(t_str) == 1: 
            mybase_list[i] += "00"                                                         
        if len(t_str) == 2:
            mybase_list[i] += "0"                                                          
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


    '''
    for i in range(len(lec)):
        mygam, myrspece = decompose_spectra(totlec[i],lec[i],gam[i])
        mygam_list.append(mygam)
        myrspece_list.append(myrspece)
    '''
    for i in range(len(lec)):
        plt.plot(totgam[i],totgam[i]*totlec[i], color=color_list[i], linewidth=0.5)
        plt.plot(gam[i],gam[i]*lec[i],color=color_list[i],linewidth=2,label='$B_{g}/B_{0}=$'+str(bguide_list[i]))


    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\gamma-1$')
    plt.ylabel('$(\gamma-1)dN/d\gamma$')
    plt.ylim(1e2,1e9)
    plt.xlim(1e0,1e4)
    plt.legend(frameon=False)
    
    plt.savefig('plots/spectra_movies/sig.3_delgam0005_triggered/time_adjusted/'+base_time_string+'.png',dpi=300,bbox_inches='tight')
    plt.close()
    for file in file_list:
        file.close()
