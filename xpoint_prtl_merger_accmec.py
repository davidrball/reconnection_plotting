import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

from tristan_funcs import get_val_filename, return_spec
from tristan_funcs_vecpot import vecpot2
from localmin_func import localmin, localmax, localmin_wrap
import h5py

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

sigma=0.3
va = np.sqrt(sigma)/np.sqrt(1+sigma)

scan = 10
interval = 2000.
istep = 12
c_omp = 3


tscan= .5

xpoint_scan = tscan * interval * .45 / istep
xpoint_scan = 25

myprtl = "../../tristan_acc-mec_Ez/highbeta/triggered_guide/output/prtl.tot.065"
#myprtl = "../../tristan_acc-mec_Ez/8k_bguide.1_allprts/output/prtl.tot.045"
myprtl = h5py.File(myprtl,'r')


tcse = np.array(get_val_filename('tcse',myprtl))/interval
ycse = np.array(get_val_filename('ycse',myprtl))/istep
gammae = np.array(get_val_filename('gammae',myprtl))
ez_bxy = -np.array(get_val_filename('ezbxye',myprtl))/va


prtnum = np.size(tcse)



t0 = 5
tf = 50

xpoint_prtl_ycs = []
xpoint_prtl_gam = []

else_prtl_ycs = []
else_prtl_gam = []

merger_prtl_ycs = []
merger_prtl_gam = []

for t in range(t0,tf):
    myfld = "../../tristan_acc-mec_Ez/highbeta/triggered_guide/output/flds.tot."
    t_str = str(t)
    print(t_str)
    if len(t_str)==1:
        fld_base = myfld + "00"
    elif len(t_str)==2:
        fld_base = myfld+"0"
    fld_base += t_str
    myfld = h5py.File(fld_base,'r')
    dens = np.rot90(get_val_filename('dens',myfld)[0,:,:])

    vecpot = vecpot2(myfld,istep,c_omp)
    #print(np.shape(vecpot))
    xmid = np.shape(vecpot)[1]/2
    vecpot_slice = vecpot[:,xmid]
   
    lowcount = 0
    for index in range(np.size(vecpot_slice)):
        if vecpot_slice[index]==0:
            pass
        if vecpot_slice[index] != 0:
            lowcount = index
            break
    vecpot_slice = vecpot_slice[lowcount:]
    #print('trimmed ' + str(lowcount) + ' zeros from bottom')
    #print('vecpot slice first 10 , ', vecpot_slice[0:10])
    highcount = 0
    for index in range(np.size(vecpot_slice)):
        modified_index = np.size(vecpot_slice)-index-1
        if vecpot_slice[modified_index]==0:
            pass
        if vecpot_slice[modified_index]!=0:
            highcount = modified_index
            break
    #print('trimmed ' + str(index) +' zeros from top')
    #print('highcount index : ', highcount)
    vecpot_slice = vecpot_slice[:modified_index]

    #print('vecpot slice last 10 , ', vecpot_slice[-10:])
    #print(np.shape(vecpot_slice))
    ylen = np.shape(vecpot)[0]
    ylenhlf = ylen/2
    time_offset = 4
    time_delay = 3
    myedge = max(int(ylenhlf - (va*.45*interval*(t-time_offset)/istep / time_delay)),scan+1)
    #myedge = 20
    #need to write a localmin func that wraps around and takes care of boundaries, didn't need this for triggered simulation, need it for untriggered
    #xpoint_list = localmin(vecpot_slice,scan, myedge)
    xpoint_list = localmin_wrap(vecpot_slice,scan)
    num_xpoints = len(xpoint_list)
    print('num xpoints ',num_xpoints)

    tlow = t-tscan
    tup = t+tscan
    
    time_ycs_xpoint = []
    time_ycs_else = []
    time_gam_xpoint = []
    time_gam_else = []
    time_ycs_merger = []
    time_gam_merger = []
    for i in range(prtnum):
        myt = tcse[i]
        if tlow < myt < tup and gammae[i]>500:#particle accelerated near this fld timestep
            #now go through xpoints and correlate particle ycs to xpoint positio     
            myycs = ycse[i]
            for j in range(num_xpoints):
                xpoint_loc = xpoint_list[j]
                xpoint_low = xpoint_loc - xpoint_scan
                xpoint_up = xpoint_loc + xpoint_scan
                #print(xpoint_low, xpoint_up)

                #let's demand that the particles associated with the xpoint are causally connected
                #c is .45 cells / timestep
                #delta x / delta t < .45 (as long as delta x and t are in the right units, cells and timesteps)

                deltax = np.abs(xpoint_loc - myycs) # in units of downsampled cells
                deltax *= istep #in units of cells
                deltat = np.abs(myt - t)*interval #in units of computational steps

                deltax_deltat = deltax / deltat

                myezbxy = ez_bxy[i]
                if xpoint_low < myycs < xpoint_up:
                #if deltax_deltat < .45: #causally connected particles
                    #print(i)
                    if myezbxy > 0:
                        xpoint_prtl_ycs.append(myycs)
                        xpoint_prtl_gam.append(gammae[i])
                        time_ycs_xpoint.append(myycs)
                        time_gam_xpoint.append(gammae[i])
                    elif myezbxy < 0:
                        merger_prtl_ycs.append(myycs)
                        merger_prtl_gam.append(gammae[i])
                        time_ycs_merger.append(myycs)
                        time_gam_merger.append(gammae[i])
                if xpoint_loc < xpoint_scan: #handling boundaries
                    overflow = np.abs(xpoint_loc-xpoint_scan)
                    upperlim = ylen - overflow
                    if myycs > upperlim:
                        if myezbxy > 0:
                            xpoint_prtl_ycs.append(myycs)
                            xpoint_prtl_gam.append(gammae[i])
                            time_ycs_xpoint.append(myycs)
                            time_gam_xpoint.append(gammae[i])
                        elif myezbxy < 0:
                            merger_prtl_ycs.append(myycs)
                            merger_prtl_gam.append(gammae[i])
                            time_ycs_merger.append(myycs)
                            time_gam_merger.append(gammae[i])
                if xpoint_loc+xpoint_scan > ylen:
                    overflow = np.abs(xpoint_loc + xpoint_scan - ylen)
                    lowerlim = overflow
                    if myycs < lowerlim:
                        if myezbxy > 0:
                            xpoint_prtl_ycs.append(myycs)
                            xpoint_prtl_gam.append(gammae[i])
                            time_ycs_xpoint.append(myycs)
                            time_gam_xpoint.append(gammae[i])
                        elif myezbxy < 0:
                            merger_prtl_ycs.append(myycs)
                            merger_prtl_gam.append(gammae[i])
                            time_ycs_merger.append(myycs)
                            time_gam_merger.append(gammae[i])
            if myycs in xpoint_prtl_ycs or myycs in merger_prtl_ycs:
                pass
            else:
                #print(i)
                if myezbxy > 0:
                    else_prtl_ycs.append(myycs)
                    else_prtl_gam.append(gammae[i])
                    time_ycs_else.append(myycs)
                    time_gam_else.append(gammae[i])
                elif myezbxy < 0:
                    merger_prtl_ycs.append(myycs)
                    merger_prtl_gam.append(gammae[i])
                    time_ycs_merger.append(myycs)
                    time_gam_merger.append(gammae[i])
    #if we don't want to plot each timestep
    
    print('num xpoint prtls : ' ,len(xpoint_prtl_ycs))
    print('num else prtls : ', len(else_prtl_ycs))
    print('num merger prtls : ', len(merger_prtl_ycs))
    vert_array = np.logspace(0,3,10)#only works if looking at a single timestep
    fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)

    for xpoint_loc in xpoint_list:
        ho_array = np.ones(10)*xpoint_loc
        ax2.plot(ho_array,vert_array,color="Black",linestyle='--')

        #for plotting all prtls over time
    #ax2.scatter(else_prtl_ycs, else_prtl_gam,color="Red")
    #ax2.scatter(xpoint_prtl_ycs, xpoint_prtl_gam,color="Blue")
    
        #only plot prtls injected in this timestep
    ax2.scatter(time_ycs_else, time_gam_else,color="Red",s=10)
    ax2.scatter(time_ycs_xpoint,time_gam_xpoint,color="Blue",s=10)
    ax2.scatter(time_ycs_merger, time_gam_merger, color="Orange",s=10)
    ax2.set_xlim(0,np.size(vecpot_slice))
    ax2.set_ylabel('$\gamma$')
    ax2.set_yscale('log')
    ax2.set_ylim(100,4e3)
    ax1.imshow(dens/4,vmin=0,vmax=5,origin='lower')
#print(np.shape(dens))
    filepath = 'xpoint_prtl/triggered_bguide.1/'
    plt.savefig(filepath+t_str+ '.png')
    plt.close()
    
    #gives locations of xpoints in downsampled cells
    #now need to identify particles with tcs closest to this timestep, and correlate them with xpoints
    myfld.close()
myprtl.close()
gamlow = min(min(xpoint_prtl_gam),min(else_prtl_gam))
gamup = max(max(xpoint_prtl_gam),max(else_prtl_gam))
bin_num = 50

histbins, xpoint_hist = return_spec(np.array(xpoint_prtl_gam),gamlow,gamup,bin_num)
histbins2, else_hist = return_spec(np.array(else_prtl_gam),gamlow,gamup,bin_num)
histbins3, merger_hist = return_spec(np.array(merger_prtl_gam),gamlow,gamup,bin_num)

plt.plot(histbins, histbins*xpoint_hist,color="Blue",label="Xpoint")
plt.plot(histbins3, histbins3*merger_hist,color="Orange",label="Merger") 
plt.plot(histbins2, histbins2*else_hist,color="Red",label="Other")
plt.xlim(500,1e4)
plt.legend(frameon=False)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$\gamma$')
plt.ylabel('$\gamma dN/d\gamma$')
spect_filepath = 'bguide.1_triggered_'
plt.savefig(spect_filepath+'spect.png',dpi=300,bbox_inches='tight')
