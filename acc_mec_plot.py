import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tristan_funcs import get_val_filename, list_types
import h5py
import numpy as np
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 15})


#myprtl = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/prtl.tot.036"

myprtl = "../../accmec_withy/bguide.3_triggered_highthresh/output/prtl.tot.045"

myprtl = h5py.File(myprtl,'r')

xcse = np.array(myprtl['xcse'])
#print('min max xcse')
#print(np.min(np.nonzero(xcse)),np.max(xcse))
#print(np.max(xcse))

gammae = np.array(get_val_filename('gammae',myprtl))

tcse = np.array(get_val_filename('tcse',myprtl))
ycse = np.array(get_val_filename('ycse',myprtl))
ez_bxye = np.array(get_val_filename('ezbxye',myprtl))
xe = np.array(myprtl['xe'])
print('minmax of xe : ', np.min(xe), np.max(xe))

c = .45
rstep_first = 12000
rstep_interval = 2000
rstep_jump = 2000
u_sh = .1
btheta = 0.3
sigma=0.3
tcse_jump = (tcse - rstep_first)//rstep_interval+1 #how many jumps there are at this tcs
tcse_jump[tcse_jump<0] = 0
#print(tcse)
#print(tcse_jump)

xcse_adjusted = xcse - c*tcse


corspeed = 1-u_sh*np.sqrt(sigma/(1+sigma*(1+btheta**2)))
print('corrected speed : ', corspeed)
xcse_jump = xcse_adjusted + c*corspeed*rstep_jump*tcse_jump

#xcse_adjusted = xcse - .45*tcse


#print(np.max(ycse), np.min(ycse))
#print(np.max(tcse), np.min(tcse))


sigma= 0.3
va = np.sqrt(sigma / (1+sigma))

 
tmax = np.max(tcse)

logmax = np.log10(tmax)

prtnum = np.size(tcse)

#print('prtnum : ', prtnum)
istep = 12.
c_omp = 3.
my0 = 8160
#ycse *= istep / c_omp



mygam_list = []
mytcs_list = []
myycs_list = []
myxcs_list = []

myezbxy_list = []
prt_count = 0
for i in range(prtnum):
    mygam = gammae[i]
    myycs = ycse[i]/c_omp 
    mytcs = tcse[i]
    myxcs = xcse_jump[i]/c_omp# - 1500/c_omp
    #mycol = np.log10(myt)/logmax
    #mycol = (0, mycol, mycol)
    if myycs != 0:# and mygam > 350: #i think when the boundaries move they inject particles with tcs=1, making noisy plots
        #mycol = (myt-tmin)/tmax                                                   
        #mycol = (0, mycol, mycol) 
        #mycol = np.log10(myt)/logmax
        #mycol = (0, mycol, mycol)
        #plt.scatter(myycs, mygam, color=mycol)
        mygam_list.append(gammae[i])
        mytcs_list.append(mytcs)
        myycs_list.append(myycs)
        myxcs_list.append(myxcs)
        prt_count += 1
        myezbxy_list.append(ez_bxye[i]/va)
#print('max ycs' , max(myycs_list))
tcs_array = np.array(mytcs_list)
tcs_min = np.min(tcs_array)
tcs_max = np.max(tcs_array)

myezbxy_array = np.array(myezbxy_list)

#col_array = (tcs_array - tcs_min)/tcs_max
col_array = tcs_array / tcs_max
print('prt count : ',prt_count) 

xext = my0 / c_omp

#plotting ez_bxy vs y
sc = plt.scatter(myycs_list, -myezbxy_array,c=col_array,cmap='magma',s=.1) 
plt.ylim(-5,5)
plt.xlim(0, xext)
plt.clim(-5,5)
plt.colorbar(sc,label='$t_{cs}/t_{max}$')
plt.ylabel('$E_{z}/v_{A}B_{xy}$')
plt.xlabel('$y_{cs} \; (c/\omega_{p})$')
plt.savefig('ycs_ez_bguide.3_triggered_accmec.png',dpi=300,bbox_inches='tight')
plt.close()
'''
sc = plt.scatter(-myezbxy_array, mygam_list, c=col_array,cmap='magma',s=.1)
plt.colorbar(sc,label='$t_{cs}/t_{max}$')
plt.xlim(-5,5)
plt.yscale('log')
plt.xlabel('$E_{z}/v_{A}B_{xy}$')
plt.ylabel('$\gamma$')
plt.savefig('ez_gam_8k_untriggered_bguide0.png',bbox_inches='tight',dpi=300)
plt.close()
'''


ho_arr = np.linspace(0,xext,10)
vert_arr = np.ones(10)*1836*sigma / 2.
sc = plt.scatter(myycs_list, mygam_list, c=col_array, cmap='magma',s=.1)     
plt.xlabel('$y_{cs} \; (c/\omega_{p})$')
plt.xlim(0,xext)
#plt.plot(ho_arr, vert_arr,color="Black",linestyle='--',label='$\sigma_{e}/2.$')
plt.legend(frameon=False)
plt.colorbar(sc,label='$t_{cs}/t_{max}$')
#plt.title('Triggered Bguide=0.3')
plt.ylabel('$\gamma$')
plt.yscale('log')
plt.savefig('ycs_gam_bguide.3_triggered_highthresh_accmec.png',bbox_inches='tight',dpi=300)     
plt.close()



sc = plt.scatter(myxcs_list, mygam_list, c=col_array, cmap='magma',s=.1)
plt.xlabel('$x_{cs} \; (c/\omega_{p})$')
#plt.xlim(0,3000/c_omp)
#plt.plot(ho_arr, vert_arr,color="Black",linestyle='--',label='$\sigma_{e}/2.$')                                                                              
plt.legend(frameon=False)
plt.colorbar(sc,label='$t_{cs}/t_{max}$')
#plt.title('Triggered Bguide=0.3')                                              
plt.ylabel('$\gamma$')
plt.yscale('log')
plt.savefig('xcs_gam_bguide.3_triggered_highthresh_accmec.png',bbox_inches='tight',dpi=300)
plt.close()

mygam_col_array = np.log10(np.array(mygam_list))
ezbxy_array = np.array(myezbxy_list)
#ezbxy_array[ezbxy_array>1]=1
#ezbxy_array[ezbxy_array<-1]=-1
sc = plt.scatter(myycs_list, myxcs_list, c=col_array, cmap='magma',s=.1)
plt.xlabel('$x_{cs} \; (c/\omega_{p})$')
#plt.xlim(0,3000/c_omp)
#plt.plot(ho_arr, vert_arr,color="Black",linestyle='--',label='$\sigma_{e}/2.$')                                                                             
plt.xlim(0,xext)
#plt.ylim(-1500,1500)
plt.legend(frameon=False)
plt.colorbar(sc,label='$\log{\gamma}$')
#plt.title('Triggered Bguide=0.3')                                                         
plt.ylabel('$y_{cs} \; (c/\omega_{p})$')
#plt.yscale('log')
plt.savefig('xcs_ycs_bguide_triggered_highthresh_accmec.png',bbox_inches='tight',dpi=300)
plt.close()



myprtl.close()


