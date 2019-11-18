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

scan = 6
interval = 2000.
istep = 12
c_omp = 3


tscan= .5

xpoint_scan = tscan * interval * .45 / istep
xpoint_scan = 20

myprtl = "../../tristan_acc-mec_Ez/8k_bguide0_triggered_stride1_thresh2/output/prtl.tot.045"
#myprtl = "../../tristan_acc-mec_Ez/8k_bguide.1_allprts/output/prtl.tot.045"
myprtl = h5py.File(myprtl,'r')


tcse = np.array(get_val_filename('tcse',myprtl))
ycse = np.array(get_val_filename('ycse',myprtl))/c_omp / 1000
gammae = np.array(get_val_filename('gammae',myprtl))
ez_bxy = -np.array(get_val_filename('ezbxye',myprtl))/va


prtnum = np.size(tcse)



t0 = 5
tf = 25
tstep=5
xpoint_prtl_ycs = []
xpoint_prtl_gam = []

else_prtl_ycs = []
else_prtl_gam = []

merger_prtl_ycs = []
merger_prtl_gam = []

for t in range(t0,tf,tstep):
    myfld = "../../tristan_acc-mec_Ez/8k_bguide0_triggered_stride1_thresh2/output/flds.tot."
    t_str = str(t)
    print(t_str)
    if len(t_str)==1:
        fld_base = myfld + "00"
    elif len(t_str)==2:
        fld_base = myfld+"0"
    fld_base += t_str
    myfld = h5py.File(fld_base,'r')

    
    dens = np.rot90(get_val_filename('dens',myfld)[0,:,:])
    dens_shape = np.shape(dens)

    

    xlen = dens_shape[1]
    xext = xlen * istep / c_omp /1000

    print('xext : ' ,xext)
    yext = dens_shape[0] * istep / c_omp / 1000

    vecpot = vecpot2(myfld,istep,c_omp)
    #print(np.shape(vecpot))
    xmid = np.shape(vecpot)[1]/2
    vecpot_slice = vecpot[:,xmid]
    
    bdens_slice = myfld['bdens'][0,:,dens_shape[0]/2]
    

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
    xpoint_list_ind = localmin_wrap(vecpot_slice,scan)

    #remove any xpoints where bdens is high
    print('original xpoint list length : ', len(xpoint_list_ind))
    #print('mean bdens : ', np.mean(bdens_slice))

    xpoint_dellist = []
    #print('max xpoint and bdens len')
    #print(max(xpoint_list_ind), np.size(bdens_slice))
    for index in range(len(xpoint_list_ind)):
        xpoint_loc = xpoint_list_ind[index]
        if bdens_slice[xpoint_loc] > 1:
            #print('bdens criterion triggered')
            #print(xpoint_loc)
            xpoint_dellist.append(index)
    xpoint_list_ind = np.array(xpoint_list_ind)
    xpoint_list_ind = np.delete(xpoint_list_ind, xpoint_dellist)
    xpoint_list = xpoint_list_ind * istep / c_omp / 1000 - xext/2.
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
    dellist = []
    prtnum = np.size(tcse)
    for i in range(prtnum):
        myt = tcse[i]
        
        if tlow < myt/interval < tup and gammae[i] > 400:#particle accelerated near this fld timestep
            #now go through xpoints and correlate particle ycs to xpoint positio           
            dellist.append(i)
            myycs = ycse[i] - xext/2.
            myezbxy = ez_bxy[i]
            for j in range(num_xpoints):
                xpoint_loc = xpoint_list[j]#in units of 1000 skin depths
                deltat = np.abs(myt - tlow*interval) #in units of computational steps
            
                #assume xpoint doesn't move
                deltax = xpoint_loc - myycs
   
                #putting this back into computational units
                deltax *= 1000 *c_omp #back in terms of cells

                
                deltax_deltat = np.abs(deltax / deltat)
                #print(deltax_deltat)
                alfven_speed = va*.45 
                if deltax_deltat <= alfven_speed:
                    #print(deltax, deltat, deltax_deltat)
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



                #print(ylen_compunits)
                '''
                if xpoint_loc < -xext/2. + xpoint_scan_causal:
                    #print(xpoint_loc)
                    overflow = xext/2. + xpoint_loc-xpoint_scan_causal
                    
                    upperlim = xext/2. + overflow
                    #print('upper lim : ',  upperlim)
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
                 
                
                elif xpoint_loc+xpoint_scan_causal > xext/2.:
                    overflow = np.abs(xpoint_loc + xpoint_scan_causal - xext/2.) #how much you've overflows by

                    #print('overflow : ', overflow)
                    lowerlim = -xext/2. + overflow
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
                
            
                '''
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
    
    #delete the entries that we've already added to lists
    gammae = np.delete(gammae,dellist)
    ycse = np.delete(ycse, dellist)
    tcse = np.delete(tcse, dellist)
    ez_bxy = np.delete(ez_bxy, dellist)
    #if we don't want to plot each timestep
    
    print('num xpoint prtls : ' ,len(xpoint_prtl_ycs))
    print('num else prtls : ', len(else_prtl_ycs))
    print('num merger prtls : ', len(merger_prtl_ycs))
    vert_array = np.logspace(3.8,4,10)#only works if looking at a single timestep
    fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)
    
    #commenting out x-point plotting
    for xpoint_loc in xpoint_list:
        ho_array = np.ones(10)*xpoint_loc
        ax2.plot(ho_array,vert_array,color="Cyan",linestyle='--')

        #for plotting all prtls over time
    #ax2.scatter(else_prtl_ycs, else_prtl_gam,color="Red")
    #ax2.scatter(xpoint_prtl_ycs, xpoint_prtl_gam,color="Blue")
    
        #only plot prtls injected in this timestep
    ax2.scatter(np.array(time_ycs_else), time_gam_else,color="Red",s=1)
    ax2.scatter(np.array(time_ycs_xpoint),time_gam_xpoint,color="Blue",s=1)
    ax2.scatter(np.array(time_ycs_merger), time_gam_merger, color="Orange",s=1)
    ax2.set_xlim(-xext/2.,xext/2.)
    ax2.set_ylabel('$\gamma$')
    ax2.set_yscale('log')
    ax2.set_ylim(400,1e4)
    ax1.imshow(dens/4,vmin=0,vmax=5,extent=[-xext/2, xext/2, -yext/2, yext/2],origin='lower')
    ax1.set_ylabel('$y \; (1000 \; c/\omega_{p})$')
    ax2.set_xlabel('$x \; (1000 \; c/\omega_{p})$')
#print(np.shape(dens))
    filepath = 'paper_plots/triggered_bguide0_snapshot_alfven'
    plt.savefig(filepath+t_str+ '.png',bbox_inches='tight',dpi=300)
    plt.close()
    
    #gives locations of xpoints in downsampled cells
    #now need to identify particles with tcs closest to this timestep, and correlate them with xpoints
    myfld.close()
myprtl.close()
