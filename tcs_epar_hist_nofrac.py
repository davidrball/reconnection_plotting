import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz
from matplotlib.colors import LogNorm
from tristan_funcs import load_particle_edotv_withtcs
import os



plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

#inputs

#mypath = "edotv_prtls/bguide.3_untriggered_hightimeres_latetime/tcs_cut/"
mypath = "edotv_prtls/bguide.3_triggered_hightimeres_nothresh/tcs_cut/"



#open all .txt files in directory pointed to by filepath
filelist = os.listdir(mypath)
#print(filelist)


#for testing purposes
#filelist = ["edotv_prtls/bguide.3_untriggered_hightimeres_latetime/1.txt"]

wpar_list = []
finalgam_list = []

for filename in filelist:
    fullname = mypath + filename
    #testing
    #print(fullname)
    
    filesize = os.stat(fullname).st_size
    if filesize == 0:
        pass
    else:
    

        d = load_particle_edotv_withtcs(fullname)
        tcs = d['tcs']
    
        tlow=0 #change this if you want to narrow window that calculation is done in
        tup = -1

        #finalgam = d['gamma'][-1]
    #indexing will only work as long as we save prt info from the first timestep on, which we do, but just something to be careful about not to change
        t_arr = np.array(d['time'])[tlow:tup]
        u_arr = np.array(d['u'])[tlow:tup]
        v_arr = np.array(d['v'])[tlow:tup]
        w_arr = np.array(d['w'])[tlow:tup]
        bx_arr = np.array(d['bx'])[tlow:tup]
        by_arr = np.array(d['by'])[tlow:tup]
        bz_arr = np.array(d['bz'])[tlow:tup]
        ex_arr = -np.array(d['ex'])[tlow:tup]
        ey_arr = -np.array(d['ey'])[tlow:tup]
        ez_arr = -np.array(d['ez'])[tlow:tup]
        gamma_arr = np.array(d['gamma'])[tlow:tup]
        finalgam = d['finalgam']
        #print('tcs in output units is : ',tcs)


        bmag = np.sqrt(bx_arr**2 + by_arr**2 + bz_arr**2)

        unitbx = bx_arr / bmag
        unitby = by_arr / bmag
        unitbz = bz_arr / bmag
        
    #unit vector in b direction is (unitbx, unitby, unitbz)                                

    #dividing momenta by gamma to turn them into velocities
        u_arr /= gamma_arr
        v_arr /= gamma_arr
        w_arr /= gamma_arr


        edotv_tot = u_arr*ex_arr + v_arr*ey_arr + w_arr*ez_arr
    
    

        upar = u_arr*unitbx + v_arr*unitby + w_arr*unitbz                                   
        epar = ex_arr*unitbx + ey_arr*unitby + ez_arr*unitbz                                   
        epar_upar = upar*epar                                                                 
        eperp_uperp = edotv_tot - epar_upar
        edotv_par_cumulative = []                                              
        edotv_perp_cumulative = []                                             
        edotv_tot_cumulative = []                                              
        for i in range(np.size(eperp_uperp)):                                  
            tmp_par = trapz(epar_upar[0:i])                                    
            tmp_perp = trapz(eperp_uperp[0:i])                                 
            tmp_tot = trapz(edotv_tot[0:i])                                    
            edotv_par_cumulative.append(tmp_par)                               
            edotv_tot_cumulative.append(tmp_tot)                               
            edotv_perp_cumulative.append(tmp_perp)

    
    #testing
    #plt.plot(t_arr, edotv_par_cumulative,color="Blue")
    #plt.plot(t_arr, edotv_perp_cumulative,color="Red")
    #plt.plot(t_arr, edotv_tot_cumulative,color="Black")
    #plt.savefig('testing_reading.png')
    
    
        wpar_frac = 200*edotv_par_cumulative[-1]*2/.45
        #print('wpar_frac : ', wpar_frac)
        #print('finalgam : ', finalgam)
        wpar_list.append(wpar_frac)
        finalgam_list.append(finalgam)



maxgam = max(finalgam_list)
mingam = min(finalgam_list)

maxwpar = max(wpar_list)
minwpar = min(wpar_list)

binnum = 40
#wpar_bins = np.linspace(-2,1,binnum)
wpar_bins = np.linspace(minwpar, maxwpar,40)
gambins = np.linspace(min(finalgam_list),max(finalgam_list)+1,binnum)


 
#plt.scatter(wpar_list,finalgam_list)
#plt.xlim(0,1.1)
#plt.savefig('testing_scatter.png')

#problem is that there are some artificially low values of wpar (down to -60 ish), so when we cut up the bins this skews the entire distribution




#not putting bins in by hand
#H,xedges,yedges=np.histogram2d(wpar_list,finalgam_list,bins=40)
  
#H,xedges,yedges = np.histogram2d(wpar_list,finalgam_list,bins=[wpar_bins,gambins])
#H,xedges,yedges = np.histogram2d(wpar_list,finalgam_list,bins=40,range=[[0,1],[mingam,maxgam]])  
#plt.yscale('log')

#print(H)
#extent_list = [int(wpar_bins[0]),int(wpar_bins[-1]),int(gambins[0]),int(gambins[-1])]
#extent_list = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
#print('bins : ', extent_list)

plt.hist2d(wpar_list,finalgam_list,bins=40,range=[[0,maxwpar],[mingam,maxgam]],norm=LogNorm(),cmin=1,cmax=2001,alpha=1)

plt.xlim(0,maxwpar)
plt.ylim(0,maxgam)


xarr = np.linspace(0,maxgam,40)
plt.plot(xarr,xarr,linestyle='dashed',color="Black")

#log scale on gamma
#plt.hist2d(wpar_list, np.log10(finalgam_list),bins=40,range=[[0,1],[np.log10(mingam),np.log10(maxgam)]],norm=LogNorm())

plt.colorbar(label='$N_{e}$')
#plt.imshow(H,origin='lower',extent=extent_list,aspect="auto")#,extent=extent_list,aspect="auto")
#plt.yscale('log')

plt.xlabel('$W_{||}$')
plt.ylabel('$\gamma_{final}$')
#plt.ylabel('$\log{\gamma_{final}}$')

plt.savefig('wpar_nofrac.png')

plt.close()
