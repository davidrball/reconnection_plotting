
import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

#inputs
#prtl_base = "../../tristan_acc-mec_Ez/testing_edotv/untriggered/output/prtl.tot."
#prtl_base = "../../tristan_acc-mec_Ez/bguide.1_wpar_zcomp/output/prtl.tot."
#prtl_base = "../../tristan_acc-mec_Ez/bguide.3_wpar_zcomp/output/prtl.tot."
#prtl_base = "../../tristan_acc-mec_Ez/bguide.3_untriggered_stride1_gamthresh50/output/prtl.tot."

#our new run with the right output params
#prtl_base = "../../tristan_acc-mec_Ez/bguide.3_wpar_zcomp/gamprt_thresh50/output/prtl.tot."

prtl_base = "../../tristan_acc-mec_Ez/highbeta_example/output/prtl.tot."

#prtl_base = "pleiades_out/prtl.tot."

savename = "bguide3_triggered_final"
tfinal = 5
gamthresh=7000 #don't include cold thermals
interval = 13000 #for labeling with time
time = .15*interval*tfinal
timestr = str(time)[:-2]
labelstr = '$t='+timestr+' \; \omega_{p}^{-1}$'


tfinal_str = "%03d" % tfinal
namefinal = prtl_base + tfinal_str
myf_final = h5py.File(namefinal,'r')
gammae_final = np.array(myf_final['gammae'])

initsize = np.size(gammae_final)


select = gammae_final > gamthresh
myind = np.nonzero(select)
gammae_final = gammae_final[myind]
wpar_final = np.abs(np.array(myf_final['edotve']))[myind]*2/.45
#not many negativs, just take abs to make us able to plot in log, take this away as needed
finalsize = np.size(gammae_final)

print('init size : ', initsize)
print('final size : ', finalsize)

mingam = gamthresh
maxgam = np.max(gammae_final)
mingamlog = np.log10(mingam)
maxgamlog = np.log10(maxgam)


minwpar = .1#np.min(wpar_final)
maxwpar = np.max(wpar_final)

minwparlog = np.log10(minwpar)
maxwparlog = np.log10(maxwpar)

#print(minwpar,maxwpar)
wpar_frac = wpar_final / gammae_final

xarr = np.linspace(0,np.max(gammae_final),50)
#plt.hist2d(wpar_final,gammae_final,bins=40,norm=LogNorm())

#plt.hist2d(np.log10(wpar_final),np.log10(gammae_final),bins=40,norm=LogNorm())
binnum=40
wparbins = np.linspace(minwpar,maxwpar,binnum)
gambins = np.linspace(mingam,maxgam,binnum)
#wparbins = np.logspace(minwparlog, maxwparlog,binnum)
#gambins = np.logspace(mingamlog,maxgamlog,binnum)
wparbins = np.linspace(minwparlog, maxwparlog,binnum)
gambins = np.linspace(np.log10(mingam),np.log10(maxgam),binnum)



print(wparbins)
print(gambins)



myhist, xedges, yedges, image = plt.hist2d(np.log10(np.abs(wpar_final)), np.log10(gammae_final),bins=(wparbins,gambins),norm=LogNorm(vmin=100,vmax=1e5))
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim(1e2,2e3)
#print('type of myhist : ', type(myhist))
print('shape of myhist : ',np.shape(myhist))
#print('gambins : ', gambins)
histshape = np.shape(myhist) 
mean_list = []
#ok time to calculate the median values
for i in range(binnum-2):
    gamnum = myhist[i,:]
    gam_times_gamnum = gamnum*10**gambins[:-1]
    gamnum = np.nan_to_num(gamnum)
    gam_times_gamnum = np.nan_to_num(gam_times_gamnum)
    meangam = np.sum(gam_times_gamnum)/np.sum(gamnum)
    

    #mymean = np.mean(gamnum)
    mean_list.append(meangam)

ax=plt.gca()
print('mean list : ',mean_list)
#for i in range(len(mean_list)):
    #patch = patches.Ellipse((wparbins[i]+.05,np.log10(mean_list[i])),.1,.075,facecolor='red')
    #ax.add_patch(patch)


sig_e = 1836*.3/2.
sig_e = np.log10(sig_e)
sig_e_arr = sig_e*np.ones(np.size(xarr))
ax.plot(xarr, sig_e_arr, color="magenta",label="$\Delta\gamma=\sigma_{e}/2$")
ax.plot(wparbins[1:-1]-.05,np.log10(mean_list),color="cyan",label='$\\bar{\Delta \gamma}$')    

#ax = plt.gca()
#patch = patches.Circle((1,2),1,facecolor="Red")
#ax.add_patch(patch)


plt.xlabel('$\log(W_{|| \; , \; z})$')
plt.ylabel('$\log(\Delta \gamma)$')
#plt.xlabel('$W_{||,z}$')
#plt.ylabel('$\gamma$')
plt.plot(xarr,xarr,linestyle='dashed',color="Orange",label='$\Delta \gamma=W_{||,z}$')
plt.colorbar(label="$N_{e}$")
plt.legend(frameon=False,loc="upper left")
#plt.savefig('testing_logscale_hist.png')
plt.ylim(1.75,4)
plt.annotate(labelstr, (2.,3.85))
plt.xlim(0,3.5)
#plt.savefig(savename+'2dhist_withmean'+str(tfinal)+'.pdf',dpi=300,bbox_inches='tight')
plt.savefig('highbeta_2dhist.png',dpi=300,bbox_inches='tight')
plt.close()


#logwpar = np.log10(wpar_frac)


'''
plt.hist2d(wpar_frac, gammae_final,bins=40,norm=LogNorm(),range=[[0,1],[mingam,maxgam]])
#plt.hist2d(wpar_frac, np.log10(gammae_final),bins=40,norm=LogNorm()),range=[[0,1],[np.log10(mingam),np.log10(maxgam)]])
plt.xlabel('$W_{||}/\gamma_{tcs}$')
plt.ylabel('$\gamma_{final}$')
plt.xlim(0,1)
plt.colorbar(label="$N_{e}$")
plt.savefig(savename+'_wparfrac_gam_tothist.png',dpi=300,bbox_inches='tight')
plt.close()
'''
