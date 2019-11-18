import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from scipy.special import kn 

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 16})


def get_val_filename(input_string, fopen):
    str_list = list(fopen.keys())
    my_string = input_string
    if my_string in str_list:
        
        my_list = fopen.get(my_string)
        return my_list
    else:
        print("Value not stored, maybe you had a typo")



def plot_powerlaw(norm,p,xlow,xup):
    x_array = np.linspace(xlow,xup,100)
    y_array = norm * x_array**(-p+1)
    #print(x_array)
    #print(y_array)
    plt.plot(x_array,y_array, color = "Blue",linestyle="--", linewidth=2.0, label="Power law\n p="+str(p))



def decompose_spectra(myspece, myrspece, mygam):
    mylen = np.size(mygam)

    returnrspece = []
    #print(np.shape(myspece), np.shape(myrspece), np.shape(mygam))              
    for i in range(mylen):

        if myspece[i] > 2*myrspece[i]:
            pass
        elif myspece[i] != 0 and mygam[i]*myspece[i] > 1000:
        
        #elif myspece[i] > np.max(myspece)/100:
            #print(myspece[i], np.max(myspece))                         
            return np.array(mygam[i:]), np.array(myrspece[i:])
            break






def plot_in_stripe_ions_filename(col, fopen1):
   
  
    ion_array = np.array(get_val_filename("rspecp", fopen1))
    ion_array = ion_array[:,0,:]
    spect_array = np.zeros(ion_array.shape[0])
    tot_ion_array = np.array(get_val_filename('specp',fopen1))
    sheet_ion_array = np.array(get_val_filename('specpb',fopen1))
    dgam = np.array(get_val_filename('dgam',fopen1))
    tot_ion_array -= sheet_ion_array

    xscan=40
    xhlf = np.shape(tot_ion_array)[1]/2
    xlow = int(xhlf-xscan)
    xup = int(xhlf+xscan)
    spect_ion_array = np.zeros(tot_ion_array.shape[0])
    for i in range(xlow,xup):
        spect_ion_array+=tot_ion_array[:,i]


    for i in range(ion_array.shape[1]):
        spect_array += ion_array[:,i]
    #print(spect_array.shape)
    
    

    plotarray = np.array(get_val_filename("gamma",fopen1))
    
    mygam, my_rspece = decompose_spectra(spect_ion_array, spect_array, plotarray)
    #plt.plot(plotarray,spect_array*plotarray**2, color = "black")
    #spect_array = spect_array - sheet_spect_array#subtracting out particles initialized in current sheet
    plt.plot(1836*mygam,mygam*my_rspece/dgam, color = col, linewidth=1.5)
    plt.plot(1836*plotarray,plotarray*spect_ion_array/dgam, color=col,linewidth=0.5) 
    plt.yscale("log")
    #plt.ylim(1,1e8)
    plt.xscale("log")
def plot_in_stripe_lecs_filename(col, fopen1, my_label):
    tot_lec_array = np.array(get_val_filename("spece",fopen1))
    sheet_lec_array = np.array(get_val_filename("speceb",fopen1))
    tot_lec_array = tot_lec_array - sheet_lec_array
    print(np.shape(tot_lec_array))
    tot_lec_array_spect = np.zeros(tot_lec_array.shape[0])
    lec_array = np.array(get_val_filename("rspece",fopen1))
    dgam = get_val_filename("dgam",fopen1)
    lec_array = lec_array[:,0,:]
    spect_array = np.zeros(lec_array.shape[0])
    print(tot_lec_array.shape, lec_array.shape)
    for i in range(lec_array.shape[1]):
        spect_array += lec_array[:,i]
    


    xscan = 40
    xhlf = tot_lec_array.shape[1]/2
    xup=int(xhlf+xscan)
    xlow=int(xhlf-xscan)
    for i in range(xlow,xup):
        tot_lec_array_spect += tot_lec_array[:,i]

    




    #print(spect_array.shape)
    plotarray = np.array(get_val_filename("gamma",fopen1))
    print(type(tot_lec_array_spect))
    print(type(spect_array))
    print(type(plotarray))
    
    mygam, my_rspece = decompose_spectra(tot_lec_array_spect, spect_array, plotarray)


    if my_label != "0":
        plt.plot(mygam, mygam*my_rspece/dgam, color = col, label=my_label,linewidth=1.5)
        plt.plot(plotarray, plotarray*tot_lec_array_spect/dgam, color=col, linewidth=.5)
    
    if my_label == "0":
        plt.plot(mygam,mygam*my_rspece/dgam, color = col,linewidth=1.5)
        plt.plot(plotarray, plotarray*tot_lec_array_spect/dgam, color=col, linewidth=.5)

    plt.yscale("log")
    #plt.ylim(1,1e8)
    plt.xscale("log")
    
    plt.xlabel('$(m_{s}/m_{e})(\gamma-1)$',size=18)
    plt.ylabel('($\gamma-1$) d$N$/d$\gamma$',size=18)

t_start = 10
t_end  = 35
step =5


for t in range(t_start,t_end+1,step):
    print(t)
    
    #myf = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/spect."
    myf = "../../tristan_acc-mec_Ez/8k_bguide.3/output/spect."
    #myf = "../../tristan-mp_reconnection/guide_fields/sig1/delgam.2/bguide2/output/spect."
    t_string = str(t)
    if t<10:
        myf+="00"
    elif t>=10 and t<100:
        myf+="0"
    elif t>=100:
        pass
    myf += t_string
    print(myf)
    f = h5py.File(myf,'r')
        
    #timestrings
    #t is in output intervals, converting this to real time:
    #.15*interval = plasma times (\omega_{p}^{-1}
    #L=16000.
    L=16000.
    c=.45
    interval=2000.
    sigma=0.3
    va=np.sqrt(sigma/(sigma+1))

    L_over_Va = ((L/c)/interval)/va

    current_crossing_time = t/L_over_Va

    time_string = str(current_crossing_time)[:4]
    
    #if int(t_string[0])%2 == 0:
    #    time_string = "0"
    #for testing other colors:
    tcol = t/t_end
    if tcol < 1:
        plot_in_stripe_lecs_filename((1-tcol,1-tcol,tcol),f,time_string) #I like these colors
    elif tcol ==1:
        plot_in_stripe_lecs_filename((1,0,0),f,time_string)
    if t_end - step<t < t_end+step:
        plot_in_stripe_ions_filename("Cyan",f)
        print('last timestep triggered')
xlow = 1e0
xup = 1e4
yup = 1e9
plt.xlim(xlow,xup)
plt.ylim(1e2,1e9)

xarr = np.linspace(1e2,1e3, 100)
p=-2.7
C=1e11
yarr = C*xarr**(p+1)
plt.plot(xarr, yarr, color="Black",linestyle='--')
#plt.annotate('$p='+str(p)+'$',(.3e5,.2e8))

#plot_maxwell(7, .7e7, "Blue")
#plt.annotate('$t / t_{A}$',(xup*.27,yup*.52))
#plt.annotate('$t / t_{A}$',(3.3e2, 1.e7),size=18)
#plt.annotate('$\sigma=1$',(1e3,2e8),size=18)
#plt.annotate('$\\beta=0.16$',(1e3,.6e8),size=18)
#plt.annotate('$B_{g}/B_{0}=0$',(1e3,1e7))
#plt.annotate('Electrons',(1.5e3, .4e9),size=18)
#plt.title('$B_{g}/B_{0}=0.1$')
#plt.legend(loc=9,bbox_to_anchor=(1.15,.95))
#plt.rcParams.update({'font.size': 12})

theta = .00005
mime = 1836
thetae = theta*mime
#plot_maxwell_lecs(thetae,6e8,"Orange")
#plot_maxwell_ions_nonrel(theta,5e10,"Black")
#plt.legend(frameon=False, loc='lower left', prop={'size':14})
plt.savefig('sig.3_bguide.3_triggered.png',dpi=300,bbox_inches='tight')
#plt.savefig('sig1_delgam.2_timespec.png', dpi=300, bbox_inches='tight')
#plt.savefig('sig1_delgam.2_timespec.pdf',bbox_inches='tight',dpi=300)
#plt.savefig('/home/dball/beta006_larger_colorspec.png', bbox_inches='tight')
#plt.savefig("/home/dball/tristan_out/real_runs/movie_plots/outflow_test/spectra/outfow_spectra.png",bbox_inches='tight')
#plt.savefig('testing_colorspec_xslice.png')
plt.close()
f.close()
