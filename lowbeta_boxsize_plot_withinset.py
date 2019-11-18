from tristan_funcs import decompose_spectra
from tristan_funcs import return_spece_withcut
from sklearn.linear_model import LinearRegression
import matplotlib
matplotlib.use('Agg') #for elgato
import numpy as np
import h5py 
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 18})

#annotate_str = "Triggered $B_{g}=0.3B_{0}$"
#savename = "triggered_bguide_3_boxsize_lesspoints"
 
annotate_str = "Triggered $B_{g}=0.1B_{0}$"
savename = "triggered_bguide_1_boxsize"
   
def plot_in_stripe_lecs(f, col):
    sheet_lecs = np.array(f["speceb"])
  
    
    sheet_lecs = np.rot90(sheet_lecs)
    sheet_spect_array = np.zeros(sheet_lecs.shape[1])
    for i in range(sheet_lecs.shape[0]):
        sheet_spect_array += sheet_lecs[i]      
    
    
    lec_array = np.array(f["rspece"])
    lec_array = lec_array[:,0,:]
    spect_array = np.zeros(lec_array.shape[0])
    for i in range(lec_array.shape[1]):
        spect_array += lec_array[:,i]
    #print(spect_array.shape)
    plotarray = np.array(f["gamma"])
    #plt.plot(plotarray, spect_array*plotarray**2, color = "black") #for comparing what happens when you subtract out ptls initialized in sheet
    spect_array = spect_array - sheet_spect_array#subtracting out particles initialized in current sheet
    #only really makes a difference at early times when high gamma is dominated by particles initialized in sheet
    plt.plot(plotarray,spect_array*plotarray, color = col)
    plt.yscale("log")
   #plt.ylim(1,1e8)
    plt.xscale("log")
   # plt.xlim(1,1e2)
    plt.xlabel('$\gamma-1$')
    plt.ylabel('($\gamma-1$) $N_{\gamma}$')

def plot_in_stripe_ions(f,col):
   
  
    ion_array = np.array(f["rspecp"])
    ion_array = ion_array[:,0,:]
    spect_array = np.zeros(ion_array.shape[0])
    for i in range(ion_array.shape[1]):
        spect_array += ion_array[:,i]
    #print(spect_array.shape)
   
    plotarray = np.array(f["gamma"])
    #plt.plot(plotarray,spect_array*plotarray**2, color = "black")
    #spect_array = spect_array - sheet_spect_array#subtracting out particles initialized in current sheet
    plt.plot(plotarray,spect_array*plotarray, color = col)
    plt.yscale("log")
    #plt.ylim(1,1e8)
    plt.xscale("log")
    #plt.xlim(1,1e2)
    plt.xlabel('$\gamma-1$')
    plt.ylabel('($\gamma-1$) d$N$/d$\gamma$')

def plot_in_stripe_lecs_filename(col, fopen1, my_label):
 
    lec_array = np.array(fopen1["rspece"])
    lec_array = lec_array[:,0,:]
    spect_array = np.zeros(lec_array.shape[0])
    for i in range(lec_array.shape[1]):
        spect_array += lec_array[:,i]
    #print(spect_array.shape)
    plotarray = np.array(fopen1["gamma",fopen1])
   
    plt.plot(plotarray,spect_array*plotarray, color = col, label=my_label)
    plt.yscale("log")
    #plt.ylim(1,1e8)
    plt.xscale("log")
    plt.xlim(5,10000)
    plt.xlabel('$\gamma-1$')
    plt.ylabel('($\gamma-1$) d$N$d$\gamma$')

def plot_in_stripe_ions_filename(col,fopen1, my_label):
    ion_array = np.array(fopen1["rspecp"])
    ion_array = ion_array[:,0,:]
    spect_array = np.zeros(ion_array.shape[0])
    for i in range(ion_array.shape[1]):
        spect_array += ion_array[:,i]
    #print(spect_array.shape)
   
    plotarray = np.array(fopen1["gamma"])
    #plt.plot(plotarray,spect_array*plotarray**2, color = "black")
    #spect_array = spect_array - sheet_spect_array#subtracting out particles initialized in current sheet
    plt.plot(plotarray,spect_array*plotarray, color = col, label = my_label)
    plt.yscale("log")
    #plt.ylim(1,2e10)
    plt.xscale("log")
    #plt.xlim(0,.01)
    plt.xscale("log")
    plt.xlabel('$\gamma-1$')
    plt.ylabel('$\left(\gamma-1\\right)$ d$N$/d$\gamma$')

#plot_in_stripe_ions_filename("red", f, "normal ppc" )
#plot_in_stripe_lecs_filename("blue",f2,"high ppc")
    
def return_rspece(file_list):
    out_lec_list = []    
    out_ion_list = []
    out_gam_list = []
    for file in file_list:
        dgam = file['dgam']
        rspece = np.array(file["rspece"])
        rspecp = np.array(file["rspecp"])
        gamarray = np.array(file["gamma"])
        rspece = rspece[:,0,:]
        rspecp = rspecp[:,0,:]
        lec_plot_array = np.zeros(rspece.shape[0])
        ion_plot_array = np.zeros(rspecp.shape[0])
        
        for i in range(rspecp.shape[1]):
            ion_plot_array += rspecp[:,i]
        for i in range(rspece.shape[1]):
            lec_plot_array += rspece[:,i]
        
        out_gam_list.append(gamarray)
        out_lec_list.append(lec_plot_array/dgam)
        out_ion_list.append(ion_plot_array/dgam)
    return out_gam_list, out_ion_list, out_lec_list


def return_spece(file_list):
    out_lec_list = []    
    out_ion_list = []
    out_gam_list = []
    for file in file_list:
        dgam = file["dgam"]
        rspece = np.array(file["spece"])
        rspecp = np.array(file["specp"])
        
        sheetlecs = np.array(file['speceb'])
        sheetions = np.array(file['specpb'])
        
        rspece -= sheetlecs
        rspecp -= sheetions
        gamarray = np.array(file["gamma"])
        #print(rspece.shape)
        #print(sheetlecs.shape)
        #rspece = rspece[:,0,:]
        #rspecp = rspecp[:,0,:]
        lec_plot_array = np.zeros(rspece.shape[0])
        ion_plot_array = np.zeros(rspecp.shape[0])
        

        xlen = rspecp.shape[1]
        #print(xlen, xlow, xup)
        for i in range(xlen):
            ion_plot_array += rspecp[:,i]
        for i in range(xlen):
            lec_plot_array += rspece[:,i]
        
        out_gam_list.append(gamarray)
        out_lec_list.append(lec_plot_array/dgam)
        out_ion_list.append(ion_plot_array/dgam)
    return out_gam_list, out_ion_list, out_lec_list
    


#for untriggered stuff 
'''
mybase = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/"

k4_ext = mybase + "bguide.3_4k_untriggered/output/spect.020"
k8_ext = mybase + "bguide.3_8k_untriggered/output/spect.040"
k16_ext = mybase +"bguide.3_16k_longtime/output/spect.072"

k4_file = h5py.File(k4_ext,'r')
k8_file = h5py.File(k8_ext,'r')
k16_file = h5py.File(k16_ext,'r')

file_list = [k4_file, k8_file, k16_file]
'''
L_list = [4080., 8160., 16320.]


#all for the triggered case

'''
mybase = "/home/u21/davidrball/david/tristan_acc-mec_Ez/boxsize_test/triggered/"
k2_ext = "2k/output/spect.015"
k4_ext = "4k/output/spect.030"
k6_ext = "6k/output/spect.045"

#8k and 16 ones are different
k8_file = "/home/u21/davidrball/david/tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/spect.060"
k16_file = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/spect.070"


k2_file = h5py.File(mybase + k2_ext,'r')
k4_file = h5py.File(mybase + k4_ext,'r')
k6_file = h5py.File(mybase + k6_ext,'r')
k8_file = h5py.File(k8_file,'r')
k16_file = h5py.File(k16_file,'r')
'''
#for lots of files
#file_list = [k2_file, k4_file, k6_file, k8_file,k16_file]
#L_list = [2040., 4080., 6120., 8160., 16320.]
#file_list = [k4_file, k8_file, k16_file]
L_list = [4080., 8160., 16320.]
#end of stuff for triggered case


#for bguide .1 untriggered
k4_file = h5py.File("/home/u21/davidrball/david/tristan_acc-mec_Ez/boxsize_test/triggered/bguide.1/untriggered_4k/output/spect.032",'r')
k8_file = h5py.File("/home/u21/davidrball/david/tristan_acc-mec_Ez/bguide.1_untriggered/output/spect.055",'r')
k16_file = h5py.File("/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1_untriggered/output/spect.075",'r')
file_list = [k4_file, k8_file,k16_file]
#L_list = [4080./3., 8160./3.,16320./3.]

L_label_list = ['1,360', '2,720','5,440']


#for bguide.1 triggered
k16_file = h5py.File("/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/spect.070",'r')
k8_file = h5py.File("/home/u21/davidrball/david/tristan_acc-mec_Ez/bguide.1_triggered_stride1_gamthresh50/output/spect.060",'r')
k4_file = h5py.File("/home/u21/davidrball/david/tristan_acc-mec_Ez/boxsize_test/triggered/bguide.1/4k/output/spect.030",'r')

file_list = [k4_file, k8_file, k16_file]
L_list = [4080., 8160., 16320.]


norm_list = []
ref_norm = L_list[-1]
power=2
for L in L_list:
    myfrac = ref_norm/L
    mynorm = myfrac**power
    norm_list.append(mynorm)


gam, ion, lec = return_rspece(file_list)
#totgam, totion, totlec = return_spece(file_list)
totgam, totion, totlec = return_spece_withcut(file_list, L_list)
mygam_list = []
myrspece_list = []

for i in range(len(lec)):
    print(i)
    mygam, myrspece = decompose_spectra(totlec[i],lec[i],gam[i])
    mygam_list.append(mygam)
    myrspece_list.append(myrspece)

#because we don't evolve to infintie time, largest simulation is compared at an "earlier" time in terms of alfven crossing, but its shape has stopped changing.  Just give it an extra boost by hand so thermal peaks line up
#hand_norm_list = [.8,1,1.2]
hand_norm_list=[1,1,1.5]
xarr = np.linspace(1e2,2e3)


'''
#keep these values for each set of boxsize sims
p_list = [-1., -1.6, -2.2, -2.5, -3.]#p list for triggered bg=0.3
N_arr = [5e9,1.5e11,4e12,2e13,2e14]#N_list for triggered bg=0.3
'''
#also for triggered bguide=0.3
#p_list = [-1.6, -2.5, -3.]
#N_arr =[1.5e11, 4e13, 2e14]

#lists for untriggered bguide.1
p_list = [-1.8, -1.9, -2.0]
N_arr = [2e11, 5e11, 1.5e12]

gamlow = 200 #set lower limit of power law
gamhigh = 800 #be conservative so we don't hit exp cutoff for small boxes

p_list = []
N_list = [1e10,1e10,1e10]


print('looping through ' + str(len(totgam)) +' sims')

for i in range(len(totgam)):
    gamlow_arg = np.argmin(np.abs(totgam[i] - gamlow))
    gamup_arg = np.argmin(np.abs(totgam[i]-gamhigh))


    X = np.log10(totgam[i])[gamlow_arg:gamup_arg].reshape(-1,1)
    Y = np.log10(totlec[i])[gamlow_arg:gamup_arg]
    
    #print(X)

    #print(X.shape)
    #print(Y.shape)
    
    #X = np.reshape(X,1,-1)
    #doing the fitting
    model = LinearRegression().fit(X,Y)
    
    print('testing model score : ', model.score(X,Y))
    print('fit values are : ', model.coef_)
    p_list.append(-model.coef_)
    #N_list.append(model.coef_[0])



#ok this is getting silly let's just automate the fits


for i in range(len(totgam)):
    print(i)
    mycol = "C"+str(i)
    plt.plot(totgam[i],hand_norm_list[i]*norm_list[i]*totgam[i]*totlec[i],label=L_label_list[i],color=mycol,linewidth=1)
    plt.plot(mygam_list[i], hand_norm_list[i]*norm_list[i]*mygam_list[i]*myrspece_list[i],color=mycol,linewidth=3)


    #for plotting lines over the top to do rough fit
    #my_yarr = N_arr[i]*xarr**p_list[i]
    #plt.plot(xarr, my_yarr, color = mycol, linestyle='dashed')
    #print(np.min(my_yarr),np.max(my_yarr))
plt.annotate(annotate_str,(13,4e8))

plt.yscale('log')
plt.xscale('log')
plt.ylim(1e2,1e9)
plt.xlim(1e1,1e4)

plt.legend(loc='upper right',prop={'size':16},frameon=False, title='$L_{x}$ ($c/\omega_{p}$)')
#plt.xlabel('$(m_{s}/m_{e})(\gamma-1$)')
plt.tick_params(axis='y',which='minor')
plt.xlabel('$\gamma-1$')
plt.ylabel('$(\gamma-1)dN/d\gamma$')
for file in file_list:
    file.close()


#inset plot
L_list_thousand = [1.36, 2.72, 5.44]

a = plt.axes([.25,.25,.3,.3])
plt.rcParams.update({'font.size': 10})
plt.scatter(L_list_thousand, p_list)

plt.xscale('log')
plt.xlabel('$L_{x}\; (1000 \; c/\omega_{p})$')
plt.ylabel('Power-law index')
plt.yticks([2,2.5,3,3.5,4])
####plt.xlim(1,10)
a.set_xticks([1,10])

plt.minorticks_off()



plt.savefig(savename+'.png',dpi=300, bbox_inches='tight')
plt.savefig(savename+'.pdf',dpi=300, bbox_inches='tight')
plt.close()
