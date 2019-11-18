import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tristan_funcs import get_val_filename, list_types
import h5py
import numpy as np
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 15})

myprtl_base = "../../tristan_acc-mec/test2/output/prtl.tot."
myfld_base = "../../tristan_acc-mec/test2/output/flds.tot."

sig=0.3
va = np.sqrt(sig/(sig+1))
istep = 12
c_omp = 3
interval = 1000


t0 = 1
tf = 4

tf_str = str(tf)

if len(tf_str) == 1:
    myprtl = myprtl_base + "00" + tf_str
elif len(tf_str)==2:
    myprtl = myprtl_base + "0" + tf_str


myprtl = h5py.File(myprtl,'r')

#loop for fld arrays

ez_bxy_list = []

plt.set_cmap('viridis')
#colors for ez_bxy plots
tsteps = tf - t0

colors = plt.cm.viridis(np.linspace(0,1,tf*interval+1))
mycount = 0

fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)

lap0 = t0*interval
for t in range(t0,tf+1):
     lap = t*interval
     #lap_index = lap-lap0
     #print(lap_index)
     print(lap)
     t_str = str(t)
     print(t_str)
     if len(t_str) == 1:
         myfld = myfld_base + "00" + t_str
     elif len(t_str)==2:
         myfld = myfld_base + "0" + t_str
     myfld = h5py.File(myfld,'r')
     ez = get_val_filename('ez',myfld)[0,:,:]
     bx = get_val_filename('bx',myfld)[0,:,:]
     by = get_val_filename('by',myfld)[0,:,:]
     #in the first loop, identify size of array in y-direction
     if t==t0:
         mysizey = np.shape(ez)[0]
         yarr = np.linspace(0,mysizey*istep/c_omp,mysizey)
     mysizex = np.shape(ez)[1]
     xhlf = int(mysizex/2)
     
     ez = ez[:,xhlf]
     bx = bx[:,xhlf]
     by = by[:,xhlf]
     ez_over_vabxy = ez / (va*np.sqrt(bx**2 + by**2))
     ez_bxy_list.append(ez_over_vabxy)

     
     ax2.plot(yarr,ez_over_vabxy, color=colors[lap])
     #mycount += 1

#plt.savefig('test_ezlist.png')

gammae = np.array(get_val_filename('gammae',myprtl))
tcse = np.array(get_val_filename('tcse',myprtl))
ycse = np.array(get_val_filename('ycse',myprtl))


tmax = np.max(tcse)

prtnum = np.size(tcse)
istep = 12.
c_omp = 3.
my0 = 1000
#ycse *= istep / c_omp



mygam_list = []
mytcs_list = []
myycs_list = []


for i in range(prtnum):
    #mygam = gammae[i]
    myycs = ycse[i]/c_omp 
    #myt = tcse[i]
    #mycol = np.log10(myt)/logmax
    #mycol = (0, mycol, mycol)
    if myycs != 0:
        #mycol = (myt-tmin)/tmax                                                   
        #mycol = (0, mycol, mycol) 
        #mycol = np.log10(myt)/logmax
        #mycol = (0, mycol, mycol)
        #plt.scatter(myycs, mygam, color=mycol)
        mygam_list.append(gammae[i])
        mytcs_list.append(tcse[i])
        myycs_list.append(myycs)

#print('max ycs' , max(myycs_list))
tcs_array = np.array(mytcs_list)
tcs_min = np.min(tcs_array)
tcs_max = np.max(tcs_array)

#col_array = (tcs_array - tcs_min)/tcs_max
col_array = tcs_array /tcs_max

sc = ax1.scatter(myycs_list, mygam_list,c=col_array,cmap='viridis',vmin=0,vmax=1)
#plt.colorbar(sc,label='$t_{cs}/t_{max}$')
ax1.set_yscale('log')
ax2.set_xlabel('$y_{cs} \; (c/\omega_{pe})$')
ax1.set_ylabel('$\gamma$')
ax1.set_ylim(1e1,1e3)
ax1.set_xlim(0,my0/c_omp)
ax2.set_ylabel('$E_{z}/v_{A}B_{xy}$')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([.85,.15,.05,.7])
fig.colorbar(sc, cax=cbar_ax, label='$t_{cs}/t_{max}$')

plt.savefig('acc_mec_withfld.png', bbox_inches= 'tight', dpi=300)


 
myprtl.close()
myfld.close()
