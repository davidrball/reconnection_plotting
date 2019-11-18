import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

from tristan_funcs import get_val_filename, return_spec
from tristan_funcs_vecpot import vecpot2
#from localmin_func import localmin, localmax, localmin_wrap
from new_xpoint_func import return_xpoints
import h5py

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

t0 = 48
tf = 49

tf_prtl=65
tf_prtl = '%03d'%tf_prtl

sigma=0.3
va = np.sqrt(sigma)/np.sqrt(1+sigma)

interval = 2000.
istep = 12
c_omp = 3


tscan= .5

#for now i think this makes the most sense to do with our xcs ycs runs since we are identifying x-points in 2d
mybase = "../../accmec_withy/bguide.3_triggered/output/"

myprtl = mybase + "prtl.tot."+tf_prtl
myprtl = h5py.File(myprtl,'r')

fld_base = mybase + "flds.tot."

tcse = np.array(myprtl['tcse'])
ycse = np.array(myprtl['ycse'])/c_omp / 1000
xcse = np.array(myprtl['xcse'])/c_omp / 1000 #let's think about smartest units to use here... I don't actually see a reason to convert them in this way
gammae = np.array(myprtl['gammae'])
ez_bxy = -np.array(myprtl['ezbxye'])/va

prtnum = np.size(tcse)



xpoint_prtl_ycs = []
xpoint_prtl_gam = []

else_prtl_ycs = []
else_prtl_gam = []

merger_prtl_ycs = []
merger_prtl_gam = []

#parameters for xpoint detection
edge=4
smoothlen=1
myroll=2
bdens_tolerence=4
zero_tolerence=2

grid_to_plot_coord = istep/c_omp / 1000. #converting output field coords to plotting units


for t in range(t0,tf):
    #myfld = "../../tristan_acc-mec_Ez/8k_untriggered_bguide.3_stride1/output/flds.tot."
    t_str = '%03d' % t
    myfld = fld_base + t_str
    print(t_str)
    myfld = h5py.File(myfld,'r')
    vecpot = vecpot2(myfld,istep,c_omp)
    #print(np.shape(vecpot))
    xpoint_loc, leftcut, rightcut, bottomcut, topcut = return_xpoints(myfld, edge,smoothlen,myroll,bdens_tolerence,zero_tolerence)

    dens = myfld['dens'][0,bottomcut:topcut,leftcut:rightcut]
    #dens = np.rot90(get_val_filename('dens',myfld)[0,:,:])                                      
    dens_shape = np.shape(dens)



    xlen = dens_shape[1]
    xext = xlen * grid_to_plot_coord
    yext = dens_shape[0] * grid_to_plot_coord
    vecpot = vecpot2(myfld,istep,c_omp)



    print('xext : ', xext)
    print('yext : ', yext)

    #print(xpoint_loc)
    #print(leftcut, rightcut, bottomcut, topcut)
    plt.imshow(np.rot90(dens), origin='lower',vmin=0,vmax=20,extent=[-yext/2,yext/2,-xext/2,xext/2])
    
    for sublist in xpoint_loc:
        x = sublist[0]*grid_to_plot_coord
        y = sublist[1]*grid_to_plot_coord
        plt.scatter(x-yext/2,y-xext/2,color="red",s=10,marker='x')
    plt.savefig('testing_newfunc.png',bbox_inches='tight',origin='lower')


    
    
    num_xpoints = len(xpoint_loc)
   
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
    '''
    '''
    vert_array = np.logspace(0,3,10)#only works if looking at a single timestep
    fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)


    #commenting out x-point plotting
    for xpoint_loc in xpoint_list:
        ho_array = np.ones(10)*xpoint_loc
        ax2.plot(ho_array,vert_array,color="Black",linestyle='--')

        #for plotting all prtls over time
    #ax2.scatter(else_prtl_ycs, else_prtl_gam,color="Red")
    #ax2.scatter(xpoint_prtl_ycs, xpoint_prtl_gam,color="Blue")
    
        #only plot prtls injected in this timestep
    ax2.scatter(np.array(time_ycs_else), time_gam_else,color="Red",s=1)
    ax2.scatter(np.array(time_ycs_xpoint),time_gam_xpoint,color="Blue",s=1)
    ax2.scatter(np.array(time_ycs_merger)-xext/2., time_gam_merger, color="Orange",s=1)
    ax2.set_xlim(-xext/2.,xext/2.)
    ax2.set_ylabel('$\gamma$')
    ax2.set_yscale('log')
    ax2.set_ylim(400,1e4)
    ax1.imshow(dens/4,vmin=0,vmax=5,extent=[-xext/2, xext/2, -yext/2, yext/2],origin='lower')
    ax1.set_ylabel('$y \; (1000 \; c/\omega_{p})$')
    ax2.set_xlabel('$x \; (1000 \; c/\omega_{p})$')
#print(np.shape(dens))
    filepath = 'xpoint_prtl/stride1/triggered_bguide.3/'
    plt.savefig(filepath+t_str+ '.png',bbox_inches='tight',dpi=300)
    plt.close()
    
    '''
    #gives locations of xpoints in downsampled cells
    #now need to identify particles with tcs closest to this timestep, and correlate them with xpoints
    myfld.close()
myprtl.close()

'''
'''
#calculating spectra
gamlow = min(min(xpoint_prtl_gam),min(else_prtl_gam))
gamup = max(max(xpoint_prtl_gam),max(else_prtl_gam))
bin_num = 50

histbins, xpoint_hist = return_spec(np.array(xpoint_prtl_gam),gamlow,gamup,bin_num)
histbins2, else_hist = return_spec(np.array(else_prtl_gam),gamlow,gamup,bin_num)
histbins3, merger_hist = return_spec(np.array(merger_prtl_gam),gamlow,gamup,bin_num)

plt.plot(histbins, histbins*xpoint_hist,color="Blue",label="Xpoint")
plt.plot(histbins3, histbins3*merger_hist,color="Orange",label="Merger") 
plt.plot(histbins2, histbins2*else_hist,color="Red",label="Other")
plt.xlim(200,1e4)
plt.legend(frameon=False)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$\gamma$')
plt.ylabel('$\gamma dN/d\gamma$')
spect_filepath = 'untriggered_bguide.3_stride1_'
plt.savefig(spect_filepath+'spect.png',dpi=300,bbox_inches='tight')
'''
