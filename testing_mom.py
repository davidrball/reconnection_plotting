import matplotlib
matplotlib.use('Agg') #for elgato
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.special import kn
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 10})

from tristan_funcs import get_val_filename

spec_base = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.3/output/momentum.0"
#spec_base = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide1/output/momentum.0"
#spec_base = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam02/output/momentum.0"

t0 = 5
tf = 30

px_std_list = []
py_std_list = []
pz_std_list = []
t_list = []


for t in range(t0, tf):
    t_str = str(t)
    if len(t_str)==1:
        spec0 = spec_base + "0"+t_str
    elif len(t_str)==2:
        spec0 = spec_base + t_str

    t_list.append(t)
    #spec0 = spec_base + str(t)
    myf = h5py.File(spec0,'r')
    px = np.array(get_val_filename('pxelogsp',myf))#[:,0:10]
    py = get_val_filename('pyelogsp',myf)#[:,0:10]
    pz = get_val_filename('pzelogsp',myf)#[:,0:10]

    pxb = np.array(get_val_filename('pxeblogsp',myf))
    print(np.shape(pxb))
    px -= pxb

    pxbin = np.array(get_val_filename('pxbin',myf))
    pybin = np.array(get_val_filename('pybin',myf))
    pzbin = np.array(get_val_filename('pzbin',myf))
    #testbin = get_val_filename('nonsense',myf)
    #print(pxbin)
    px_spec = np.zeros(np.size(pxbin))
    py_spec = np.zeros(np.size(pybin))
    pz_spec = np.zeros(np.size(pzbin))

    xsize = np.shape(px)[1]
    #print(xsize)
    xhlf = xsize / 2
    xindex = int(xhlf + 10)
    
    for i in range(np.size(px_spec)):
        px_spec[i] += np.mean(px[i,xindex-1:xindex+1])
        py_spec[i] += np.mean(py[i,xindex-1:xindex+1])
        pz_spec[i] += np.mean(pz[i,xindex-1:xindex+1])
        #print(px_spec[i])
    
    '''
    print(np.shape(px))
    for i in range(np.size(px_spec)):
        px_spec += (px[:,-10])
        py_spec += (py[:,-10])
        pz_spec += (pz[:,-10])
        #print(px_spec[i])
    '''
    print(t)
    #print(np.shape(px[0,-10:]))
    #print(np.shape(px))
    #print(np.shape(pxbin))

    #print(np.shape(px_spec))

    #calculating the dispersion
    #this is px*abs(px)*dN/dpx
    px_sum = np.sum(px_spec*np.abs(pxbin)*pxbin)
    py_sum = np.sum(py_spec*np.abs(pybin)*pybin)
    pz_sum = np.sum(pz_spec*np.abs(pzbin)*pzbin) 


    #print(pxbin)
    px_spec_norm = px_spec*np.abs(pxbin) 
    py_spec_norm = py_spec*np.abs(pybin)
    pz_spec_norm = pz_spec*np.abs(pzbin)


    px_spec_norm_sum = np.sum(px_spec_norm)
    py_spec_norm_sum = np.sum(py_spec_norm)
    pz_spec_norm_sum = np.sum(pz_spec_norm)
    #so avg is px_sum / px_spec_norm_sum
    xmean = px_sum / px_spec_norm_sum
    ymean = py_sum / py_spec_norm_sum
    zmean = pz_sum / pz_spec_norm_sum

    xdiff = (pxbin - xmean)**2
    #print(xdiff)
    ydiff = (pybin - ymean)**2
    zdiff = (pzbin - zmean)**2

    xdiff = np.abs(pxbin-xmean)**2

    xdiff *= px_spec_norm
    ydiff *= py_spec_norm
    zdiff *= pz_spec_norm
    #print(px_spec_norm)
    #print(np.sum(xdiff))

    #px_std = np.sqrt(np.sum(xdiff)/px_spec_norm_sum)
    px_std = np.sqrt(np.sum(xdiff)/px_spec_norm_sum)
    py_std = np.sqrt(np.sum(ydiff)/py_spec_norm_sum)
    pz_std = np.sqrt(np.sum(zdiff)/pz_spec_norm_sum)
    #print('mean x : ', xmean)  
    #print('px std : ',px_std)
    #print('py std : ', py_std)
    #print('pz std : ',pz_std)
    #print(pxbin)
    px_std_list.append(px_std)
    py_std_list.append(py_std)
    pz_std_list.append(pz_std)
    myf.close()
'''
print('mean x : ', xmean)
#binlist.pop(argmax)
plt.plot(pxbin,np.abs(pxbin)*px_spec,color="Blue")
ho_arr = np.ones(10)*xmean
vert_arr = np.linspace(0,np.max(px_spec)/100000,10)

ho_arr_low = np.ones(10)*(xmean-px_std)
ho_arr_up = np.ones(10)*(xmean+px_std)

plt.plot(ho_arr, vert_arr,color="Black",linewidth=3.)
plt.plot(ho_arr_low, vert_arr, color="Black",linestyle='--')
plt.plot(ho_arr_up, vert_arr, color="Black",linestyle='--')


#plt.plot(pybin, py_spec*abs(pybin),color="Green")
ho_arr_y = np.ones(10)*ymean
ho_arr_y_low = np.ones(10)*(ymean-py_std)
ho_arr_y_up = np.ones(10)*(ymean+py_std)

#plt.plot(ho_arr_y, vert_arr,color="Black",linewidth=3.)
#plt.plot(ho_arr_y_low, vert_arr, color="Black",linestyle='--')
#plt.plot(ho_arr_y_up, vert_arr, color="Black",linestyle='--')



#plt.xscale('log')
plt.xlim(-1000,1000)
plt.yscale('log')
plt.xlabel('$p_{x}$')
plt.ylabel('dN/dpx')

plt.savefig('px_spec_allx_t0.png')
'''


plt.plot(t_list, px_std_list,color="Blue",label="$p_{x}$")
plt.plot(t_list, py_std_list,color="Green",label="$p_{y}$")
plt.plot(t_list, pz_std_list,color="Orange",label="$p_{z}$")

rstep_arr = np.arange(6,max(t_list),6)
arrsize = np.size(rstep_arr)
rstep_up = np.ones(arrsize)*np.max(px_std_list)

plt.scatter(rstep_arr, rstep_up,color="Black",label="right clean")

#plt.yscale('log')
plt.legend()
plt.ylabel('Dispersion')
plt.xlabel('Time (output intervals)')
plt.savefig('delgam005_bguide.3_t_vs_disp.png',bbox_inches='tight',dpi=300)
plt.close()  


'''
    plt.plot(pxbin, px_spec, color="Blue",label="$p_{x}$")
    plt.plot(pybin, py_spec, color="Green",label="$p_{y}$")
    plt.plot(pzbin, pz_spec, color="Red", label="$p_{z}$")

    plt.xlabel('$p_{i}$')
    plt.ylabel('$dN/dp$')
    plt.xscale('log')
    plt.yscale('log')
    #plt.ylim(1e3,1e6)
    plt.legend()
    plt.xlim(1e-2,1e2)
    plt.savefig('momspec/'+str(t)+'.png',bbox_inches='tight')
    myf.close()
    plt.close()
 '''

#test = get_val_filename('nonsense',myf)

