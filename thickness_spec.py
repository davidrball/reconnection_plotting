from tristan_funcs import decompose_spectra
from tristan_funcs import return_spece_withcut
import matplotlib
matplotlib.use('Agg') #for elgato
import numpy as np
import h5py 
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 18})

annotate_str = "Untriggered $B_{g}=0.3B_{0}$"
savename = "untriggered_bguide.3_thickness_withtrig"
    
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


mybase = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/"
thick20_ext = mybase + "bguide.3_untriggered_stripe20/output/spect.060"
thick40_ext = mybase + "bguide.3_untriggered/output/spect.060"
thick60_ext = mybase + "bguide.3_untriggered_stripe60/output/spect.097"

triggered_ext = mybase + "bguide.3/output/spect.053"


thick20_file = h5py.File(thick20_ext,'r')
thick40_file = h5py.File(thick40_ext,'r')
thick60_file = h5py.File(thick60_ext,'r')

triggered_file = h5py.File(triggered_ext,'r')

file_list = [thick20_file, thick40_file, thick60_file, triggered_file]

thick_list = [20, 40, 60] #in units of cells I think, so divide by c_omp to go to skin depths
thick_list_comp = []
div3 = lambda a : str(a/3)[:4]
for i in range(len(thick_list)):
    thick_list_comp.append(div3(thick_list[i]))

thick_list_comp.append('Triggered')


L_list = [16000, 16000, 16000, 16000]
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

xarr = np.linspace(1e2,2e3)


'''
#keep these values for each set of boxsize sims
p_list = [-1., -1.6, -2.2, -2.5, -3.]#p list for triggered bg=0.3
N_arr = [5e9,1.5e11,4e12,2e13,2e14]#N_list for triggered bg=0.3
'''

p_list = [-1.6, -2.5, -3.]
N_arr =[1.5e11, 2e13, 2e14]

#lists for untriggered
#p_list = [-1.95, -2., -2.]
#N_arr = [1.5e12, 4e11, 2e10]


for i in range(len(totgam)):
    print(i)
    mycol = "C"+str(i)
    mylin = 'solid'
    if i==3:
        mycol="Black"
        #mylin = 'dashed'
    plt.plot(totgam[i],totgam[i]*totlec[i],label=thick_list_comp[i],color=mycol,linewidth=1,linestyle=mylin)
    plt.plot(mygam_list[i],mygam_list[i]*myrspece_list[i],color=mycol,linewidth=3,linestyle=mylin)


    #for plotting lines over the top to do rough fit
    #my_yarr = N_arr[i]*xarr**p_list[i]
    #plt.plot(xarr, my_yarr, color = mycol, linestyle='dashed')
    #print(np.min(my_yarr),np.max(my_yarr))
plt.annotate(annotate_str,(7,4e8))

plt.yscale('log')
plt.xscale('log')
plt.ylim(1e2,1e9)
plt.xlim(5,1e4)

plt.legend(loc='upper right',prop={'size':14},frameon=False, title='$\Delta$ ($c/\omega_{p}$)')
#plt.xlabel('$(m_{s}/m_{e})(\gamma-1$)')
plt.tick_params(axis='y',which='minor')
plt.xlabel('$\gamma-1$')
plt.ylabel('$(\gamma-1)dN/d\gamma$')
for file in file_list:
    file.close()


#inset plot
true_p_list = [] #convert negative p with extra factor of gamma to actual distribution function p

for p in p_list:
    true_p_list.append(-p + 1)
 
'''
a = plt.axes([.25,.25,.3,.3])
plt.rcParams.update({'font.size': 10})
plt.scatter(L_list_thousand, true_p_list)
plt.xlabel('$L_{x}\; (1000 \; c/\omega_{p})$')
plt.ylabel('Power-law index')
plt.yticks([2,2.5,3,3.5,4])

plt.xscale('log')
plt.xticks([4,10])
'''
plt.savefig(savename+'.png',dpi=300, bbox_inches='tight')
plt.savefig(savename+'.pdf',dpi=300, bbox_inches='tight')
plt.close()
