import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
                               
from tristan_funcs import get_val_filename
from tristan_funcs_vecpot import vecpot2
from edge_detector import convolve
from localmin_func import localmin, localmax
import h5py

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

sigma=.3
va = np.sqrt(sigma)/np.sqrt(1+sigma)
interval = 2000
istep = 12

scan = 10

t0 = 5
tf = 45


matrixbase = "00"

fld_base = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.7/output/flds.tot."

#fld_base = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam00005/output/flds.tot."

#fld_base = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig1/delgam.2/output/flds.tot."

#fld_base = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig1/delgam.2/bguide.3/output/flds.tot."

xpoint_list = []
plasmoid_list = []
for t in range(t0, tf):
    tstr = str(t)



    if len(tstr)==1:
        myfld = fld_base + "00" + tstr
    elif len(tstr)==2:
        myfld = fld_base + "0"+tstr

    

    myfld = h5py.File(myfld,'r')
    dens = get_val_filename('dens',myfld)[0,:,:]
    #bdens = get_val_filename('bdens',myfld)[0,:,:]
    vecpot = vecpot2(myfld,12,3)
    fig, (ax0, ax1) = plt.subplots(2,1,sharex=True)
 
    ylen = np.shape(vecpot)[0]
    #offset /= istep
    ylenhlf = ylen/2 

    time_offset = 4
    time_delay = 3
    myedge = max(int(ylenhlf  - (va * .45 * interval *(t-time_offset) / istep/time_delay)),scan+1)
    print('myedge : ', myedge)
    xmid = np.shape(vecpot)[1]/2

    #bdens_slice = bdens[:,xmid]
    

    vecpot_slice = vecpot[:,xmid]
    xarr = np.linspace(0,np.size(vecpot_slice),np.size(vecpot_slice))
    my_localmin_list = localmin(vecpot_slice,scan,myedge)
    my_localmax_list = localmax(vecpot_slice,scan,myedge)
    num_xpoints = len(my_localmin_list)
    num_plasmoids = len(my_localmax_list)
    xpoint_list.append(num_xpoints)
    plasmoid_list.append(num_plasmoids)
    print('num of xpoints based on localmin \n ' + str(num_xpoints))
    print('num of plasmoids \n ' + str(num_plasmoids))
    vert_array = np.linspace(0,np.max(vecpot_slice),20)
    ho_array1 = np.ones(20)*myedge
    ho_array2 = np.ones(20)*(ylen - myedge)
    ax0.imshow(np.rot90(dens),origin='lower',cmap='viridis',vmin=0,vmax=20)
    for i in range(num_xpoints):
        myx = my_localmin_list[i]
        myarg = np.argmin(np.abs(myx - xarr))
        myy = vecpot_slice[myarg]
        ax1.scatter(myx,myy,color="Red")
    for i in range(num_plasmoids):
        myx = my_localmax_list[i]
        myarg = np.argmin(np.abs(myx - xarr))
        myy = vecpot_slice[myarg]
        ax1.scatter(myx,myy,color="Blue")
        
        
    ax1.set_ylabel('$A_{z}$')
    ax1.plot(xarr,vecpot_slice)
    ax1.plot(ho_array1, vert_array,color="Black",linestyle='--')
    ax1.plot(ho_array2, vert_array,color="Black",linestyle='--')
    ax1.set_ylim(0,np.median(vecpot_slice)*1.5)
    plt.savefig('xpoint_plots/sig.3_delgam.005_bguide.7/'+tstr+".png",bbox_inches='tight')
    plt.close()
print('mean xpoints:')
print(np.mean(np.array(xpoint_list)[10:45]+np.array(plasmoid_list)[10:45]))
