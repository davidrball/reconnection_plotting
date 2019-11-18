
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
import matplotlib.patches as patches

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})



mybase = "../../tristan_acc-mec_Ez/testing_edotv/bguide.3_triggered_alwaystrack/output/"


tselect = 43
gamthresh = 2000

interval = 2000
tscan = 5
istep = 12
c_omp = 3

t0 = 1
tf = 40



tstr = '%03d' % tselect
#print(tstr)
myprt = mybase + "prtl.tot." + tstr
#identify particles we want to look at
myprt = h5py.File(myprt, 'r')


gammae = np.array(myprt['gammae'])
ind  = np.array(myprt['inde'])
proc = np.array(myprt['proce'])
tcs = np.array(myprt['tcse'])
ezbxy = np.array(myprt['ezbxye'])



#go through tcs and pick out particles accelerated near the time we want
prtnum = np.size(gammae)
bool_array = gammae > gamthresh 


myind = np.nonzero(bool_array)

select_ind = ind[myind]
select_proc = proc[myind]
select_len = np.size(select_ind)
select_tcs = tcs[myind]
select_ezbxy = ezbxy[myind]

maxtcs = np.max(select_tcs)
mintcs = np.min(select_tcs)
print('min tcs : ',mintcs)
#print(select_proc)
print('selecting ',np.size(select_ind),' particles')

myprt.close()


t_list = []
gam_list = []
wpar_list = []


for i in range(np.size(select_ind)):
    gam_list.append([])
    wpar_list.append([])
for t in range(t0, tf):
    t_list.append(t)
    tstr = '%03d' % t
    print(tstr)
    myprt = mybase + "prtl.tot." + tstr
    myprt = h5py.File(myprt,'r')
    ind = np.array(myprt['inde'])
    proc = np.array(myprt['proce'])
    gammae = np.array(myprt['gammae'])-1
    wpar = np.array(myprt['edotve'])*2/.45
    mylen = np.size(ind)

    time_wpar_list = []
    time_gam_list = []
    #fig,(ax1,ax2) = plt.subplots(1,2)
    
    for i in range(select_len):
        #print(myind, myproc)
        myind = select_ind[i]
        myproc = select_proc[i]
        mytcs = select_tcs[i]
        myezbxy = select_ezbxy[i]
        for j in range(mylen):
            if ind[j]==myind and proc[j]==myproc:
                #not sure if order is preserved for tracking indidual particles through time here but let's try
                tcol = (mytcs - mintcs)/(maxtcs - mintcs)
                #print(tcol)
                mygam = gammae[j]
                mywpar = wpar[j]
                if myezbxy > 0: 
                    mymark = "s"
                    scattercol = "Cyan"
                elif myezbxy < 0:
                    mymark = "o"
                    scattercol = "Orange"

                if mywpar > 0:#let's see what happens here...
                    gam_list[i].append(mygam)
                    wpar_list[i].append(mywpar)
                    time_wpar_list.append(mywpar)
                    time_gam_list.append(mygam)


                plt.plot(wpar_list[i],gam_list[i],color=(tcol,0,1-tcol),linewidth=0.7,zorder=1)
                plt.scatter(mywpar, mygam,color=scattercol,s=4,zorder=2,marker=mymark)
                break
    #plt.scatter(time_wpar_list,time_gam_list,color="Orange",s=3,zorder=2)
    myprt.close()
    xarr = np.linspace(0,5000)
    plt.plot(xarr, xarr,linestyle='dashed',color="Black")
    plt.xlabel('$W_{||}$')
    plt.ylabel('$\gamma-1$')
    plt.xlim(1e-3,5000)
    plt.ylim(1e-3,5000)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('edotv_gam_timetrack/triggered_bguide.3/alwaystrack_abs/'+tstr+'.png',dpi=300,bbox_inches='tight')
    plt.close()   

    

