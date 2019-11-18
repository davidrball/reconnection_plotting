import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import h5py
from matplotlib.colors import LogNorm
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 8})
plt.set_cmap("viridis")

#file name
mybase = "/home/u21/davidrball/david/tristan_acc-mec_Ez/bguide.1_triggered_stride1_gamthresh50/output/"

prtl_base = mybase + "prtl.tot."
spec_base = mybase + "spect."
fld_base = mybase + "flds.tot."


tstart = 6
finalt = 58
finalt_str = "%03d" % finalt
final_prtl = h5py.File(prtl_base+finalt_str,'r')

final_ycs = np.array(final_prtl['ycse'])
final_gam = np.array(final_prtl['gammae'])
final_tcs = np.array(final_prtl['tcse'])
final_wpar = np.array(final_prtl['edotve'])*2/.45


#plotting params
istep = 12
c_omp = 3
interval = 2000
conv = istep / c_omp / 1000

#just grabbing grid dimensions to rescale final_ycs
testdens = h5py.File(fld_base + finalt_str,'r')['dens'][0,:,:]
testshape = testdens.shape
ylen = testshape[0] # zero because it's not rotated
yext = ylen*conv/2

final_ycs /= (c_omp*1000)
final_ycs -= yext

#plotting stuff for 2d hists
mingam = 50
maxgam = np.max(final_gam)
mingamlog = np.log10(mingam)
maxgamlog = np.log10(maxgam)
minwpar = .1
maxwpar = np.max(final_wpar)
minwparlog = np.log10(minwpar)
maxwparlog = np.log10(maxwpar)



binnum = 40

wparbins = np.linspace(minwparlog, maxwparlog,binnum)
gambins = np.linspace(mingamlog, maxgamlog, binnum)

xarr = np.linspace(-1,4,binnum)
print('about to start loop')

for t in range(tstart,finalt):
    tstr = "%03d" % t
    print(tstr)
    myfld = h5py.File(fld_base + tstr,'r')
    myspec = h5py.File(spec_base + tstr,'r')
    myprtl = h5py.File(prtl_base + tstr,'r')

    xscan = 50
    dens = np.rot90(myfld['dens'][0,:,:])/4.
    
    shape = dens.shape
    xlen = shape[0]
    xhlf = xlen/2
    xlow = int(xhlf - xscan)
    xup = int(xhlf + xscan)
    dens = dens[xlow:xup,:]


    xext = xscan*conv
    #defining our axes
    '''
    ax0 = plt.subplot2grid((4,4),(0,0),colspan=4)
    ax1 = plt.subplot2grid((4,4),(1,0),colspan=4)
    ax2 = plt.subplot2grid((4,4),(2,0),colspan=2,rowspan=2)
    ax3 = plt.subplot2grid((4,4),(2,2),colspan=2,rowspan=2)
    '''
    #trying something out for spacing

    ax0 = plt.subplot2grid((8,8),(0,0),colspan=8,rowspan=2)                          
    ax1 = plt.subplot2grid((8,8),(2,0),colspan=8,rowspan=2)                          
    ax2 = plt.subplot2grid((8,8),(4,0),colspan=3,rowspan=4)                          
    ax3 = plt.subplot2grid((8,8),(4,5),colspan=3,rowspan=4)     


    #fld plotting
    ax0.imshow(dens,origin='lower',vmin=1,vmax=5,extent=[-yext,yext,-xext,xext])
    ax1.set_xticks([])

    #prt acc plotting
    myt = interval*t #in same units now as tcs
    tscan = 500
    bool_arr1 = final_tcs < myt+tscan
    bool_arr2 = final_tcs > myt-tscan
    bool_arr = bool_arr1 * bool_arr2
    
    my_ycs = final_ycs*bool_arr
    my_gam = final_gam*bool_arr

    ax1.scatter(my_ycs, my_gam,s=.1)
    ax1.set_yscale('log')
    ax1.set_xlim(-yext,yext)
    ax1.set_ylim(1e2,1e4)
    ax0.set_ylabel('$y$ ($1000 \; c/\omega_{p}$)')
    ax1.set_ylabel('$\gamma_{\\rm{f}}$')
    ax0.set_xlabel('$x$ ($1000 \; c/\omega_{p}$)')
    ax0.xaxis.tick_top()
    ax0.xaxis.set_label_position('top')


    #particle stuff for histogram
    gammae = np.array(myprtl['gammae'])
    wpar = np.array(myprtl['edotve'])*2/.45 
    logwpar = np.log10(wpar)
    loggam = np.log10(gammae)
    ax2.hist2d(logwpar, loggam, bins = (wparbins, gambins),norm=LogNorm())
    ax2.plot(xarr,xarr,linestyle='--',color="orange")
    ax2.set_xlabel('$\log(W_{||,z})$')
    ax2.set_ylabel('$\log(\gamma)$')

    #spectra
    rspece = np.array(myspec['rspece'])
    mygam = np.array(myspec['gamma'])
    rspece = rspece[:,0,:]
    lec = np.zeros(rspece.shape[0])
    for i in range(rspece.shape[1]):
        lec += rspece[:,i]
    ax3.plot(mygam, mygam*lec,color="C1")
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlim(10,1e4)
    ax3.set_ylim(1e2,1e7)
    ax3.set_xlabel('$\gamma-1$')
    ax3.set_ylabel('$(\gamma-1)dN/d\gamma$')
        


    #plt.tight_layout()


    plt.savefig('accmec_movie/bguide.1_triggered/'+tstr[1:]+'.png',dpi=300,bbox_inches='tight')


    myfld.close()
    myspec.close()
    myprtl.close()
