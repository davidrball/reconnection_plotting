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
stride=1

t0 = 1
tf = 45

tscan= .5

xpoint_scan = tscan * interval * .45 / istep
xpoint_scan = 25

#myprtl = "../../tristan_acc-mec_Ez/8k_bguide.3/output/prtl.tot.035"

#mybase = "../../tristan_acc-mec_Ez/8k_untriggered_bguide.3/output/"
mybase = "../../tristan_acc-mec_Ez/8k_bguide0_triggered_stride1_thresh2/output/"
tf_str = "%03d" % tf

myprtl = mybase + "prtl.tot." + tf_str

myspec = mybase + "spect." + tf_str

#myprtl = "../../tristan_acc-mec_Ez/8k_bguide.1_allprts/output/prtl.tot.045"
myprtl = h5py.File(myprtl,'r')
myspec = h5py.File(myspec,'r')

tcse = np.array(get_val_filename('tcse',myprtl))
ycse = np.array(get_val_filename('ycse',myprtl))/c_omp / 1000
gammae = np.array(get_val_filename('gammae',myprtl))-1
ez_bxy = -np.array(get_val_filename('ezbxye',myprtl))/va

#print('max gammae from prtls : ',np.max(gammae))
#print('min gammae from prtls : ',np.min(gammae))
prtnum = np.size(tcse)



xpoint_prtl_ycs = []
xpoint_prtl_gam = []

else_prtl_ycs = []
else_prtl_gam = []

merger_prtl_ycs = []
merger_prtl_gam = []

for t in range(t0,tf):
    #myfld = "../../tristan_acc-mec_Ez/8k_bguide.3/output/flds.tot."
    myfld = mybase + "flds.tot."
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
    #vdens = myfld['bdens'][0,:,:]
    xlen = dens_shape[1]
    xext = xlen * istep / c_omp /1000

    yext = dens_shape[0] * istep / c_omp / 1000

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
    xpoint_list = np.array(localmin_wrap(vecpot_slice,scan))*istep / c_omp / 1000 - xext/2.
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
        
        if tlow < myt/interval < tup:# and gammae[i] > 200:#particle accelerated near this fld timestep
            #now go through xpoints and correlate particle ycs to xpoint positio           
            dellist.append(i)
            myycs = ycse[i] - xext/2. 
            for j in range(num_xpoints):
                xpoint_loc = xpoint_list[j]#in units of 1000 skin depths
                
                #for a fixed scan difference
                #xpoint_low = xpoint_loc - xpoint_scan
                #xpoint_up = xpoint_loc + xpoint_scan
                
                #deltat = np.abs(myt - t*interval) #in units of computational steps

                deltax = np.abs(xpoint_loc - myycs)
                deltat = tscan*interval
                #the maximum distance a particle could have traveled in this time in terms of computational units is:
                speed_adjust = 2. #if we want the connection to x-point to be slower than speed of light
                xpoint_scan_causal = .45*deltat/speed_adjust #in computational units
                xpoint_scan_causal /= (c_omp * 1000) #into 1000 c_omp
                #print(xpoint_scan_causal)
                xpoint_low = xpoint_loc - xpoint_scan_causal
                xpoint_up = xpoint_loc + xpoint_scan_causal
                #print(xpoint_low,xpoint_up)
                myezbxy = ez_bxy[i]
                deltax_deltat = deltax / deltat
                #if xpoint_low < myycs < xpoint_up:
                #print(deltax_deltat)
                #if deltax_deltat < .45: #causally connected particles
                if deltax < xpoint_scan_causal:
                
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
                if xpoint_loc < xpoint_scan_causal: #handling boundaries
                    overflow = np.abs(xpoint_loc-xpoint_scan_causal)
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
                if xpoint_loc+xpoint_scan_causal > ylen:
                    overflow = np.abs(xpoint_loc + xpoint_scan_causal - ylen)
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
    
    #delete the entries that we've already added to lists
    gammae = np.delete(gammae,dellist)
    ycse = np.delete(ycse, dellist)
    tcse = np.delete(tcse, dellist)
    ez_bxy = np.delete(ez_bxy, dellist)
    #if we don't want to plot each timestep
    
    print('num xpoint prtls : ' ,len(xpoint_prtl_ycs))
    print('num else prtls : ', len(else_prtl_ycs))
    print('num merger prtls : ', len(merger_prtl_ycs))
    
    myfld.close()
myprtl.close()

#calculating spectra
#gamlow = min(min(xpoint_prtl_gam),min(else_prtl_gam))
#gamup = max(max(xpoint_prtl_gam),max(else_prtl_gam))
#bin_num = 50


#grabbing gams from spectra
specgam = np.array(myspec['gamma'])
dgam = np.array(myspec['dgam'])[0]
#print('dgam : ',10**dgam)

bin_num=200
gamlow = np.min(specgam)
gamup = np.max(specgam)
rspec = np.array(myspec['spece'])#[:,:]
#print('shape of spece : ',np.shape(rspec))
speceb = np.array(myspec['speceb'])#[:,:]

rspec -= speceb
#print('size of rspec : ',np.shape(rspec))
rspec_sum = np.zeros(np.shape(rspec)[0])
#print(np.size(rspec_sum))

for i in range(np.size(rspec_sum)):
    rspec_sum[i] += np.sum(rspec[i,:])
    



#print('min from spec gam : ', np.min(specgam))
corr=1.
plt.plot(specgam, specgam*rspec_sum/corr,color="Black",label='Total Spectra')


histbins, xpoint_hist = return_spec(np.array(xpoint_prtl_gam),gamlow,gamup,bin_num)
histbins2, else_hist = return_spec(np.array(else_prtl_gam),gamlow,gamup,bin_num)
histbins3, merger_hist = return_spec(np.array(merger_prtl_gam),gamlow,gamup,bin_num)


#plt.plot(histbins, histbins*xpoint_hist,color="Blue",label="Xpoint")
#plt.plot(histbins3, histbins3*merger_hist,color="Orange",label="Merger") 
#plt.plot(histbins2, histbins2*else_hist,color="Red",label="Other")
totprtspec = xpoint_hist + merger_hist + else_hist

plt.plot(histbins, stride*xpoint_hist,color="Blue",label="Xpoint")          
plt.plot(histbins, stride*merger_hist,color="Orange",label="Merger")      
plt.plot(histbins, stride*else_hist,color="Red",label="Other")             
#plt.plot(histbins, stride*totprtspec,color="Black",linestyle='--')

#print('specgam minmax : ', np.min(specgam),np.max(specgam))
#print('histbin minmax : ', np.min(histbins), np.max(histbins))
#print(histbins)
#print(specgam)

sigma_e= 1836*.3
plt.xlim(sigma_e/10.,1e4)
plt.legend(frameon=False)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$\gamma-1$')
plt.ylabel('$(\gamma-1) dN/d\gamma$')

spect_filepath = 'triggered_bguide0_stride1_'
plt.savefig(spect_filepath+'spect.png',dpi=300,bbox_inches='tight')
myspec.close()
