
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
import matplotlib.patches as patches
from tristan_funcs import return_spec


plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 10})





#mybase = "../../tristan_acc-mec_Ez/8k_untriggered_bguide.3/output/"
mybase = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/"
#mybase = "../../tristan_acc-mec_Ez/8k_bguide.3/output/"


tselect = 48
interval = 2000
tscan = 2
istep = 12
c_omp = 3

t0 = 40
tf = 52



tstr = '%03d' % tselect
#print(tstr)

myfld = mybase + "flds.tot." + tstr
myprt = mybase + "prtl.tot." + tstr

#define time interval we look in
tlow = (tselect - tscan)*interval
tup = (tselect + tscan)*interval

#identify particles we want to look at
myprt = h5py.File(myprt, 'r')
myfld = h5py.File(myfld,'r')
#print(dict(myprt).keys())

tcs = np.array(myprt['tcse'])
ind = np.array(myprt['inde'])
ycs = np.array(myprt['ycse'])
proc = np.array(myprt['proce'])
ez_bxy = np.array(myprt['ezbxye'])

x = np.array(myprt['xe'])
y = np.array(myprt['ye'])

#go through tcs and pick out particles accelerated near the time we want
prtnum = np.size(tcs)
bool_array = tlow < tcs 
bool_array2 = tcs < tup
bool_array3 = ez_bxy > 0  #>0 for merge, < 0 for xpoint


#lazily getting max / min yvals
ymax = np.max(y)
#if we only want particles with ycs near central
#yhlf = ymax / 2.
#put yhlf where you want point we scan away from
yfrac_hlf = 1.8/3.
yhlf = ymax * yfrac_hlf


yscan_frac = 1/15.
yscan = ymax * yscan_frac
ylow = yhlf - yscan
yup = yhlf + yscan

bool_array4 = ycs < yup
bool_array5 = ycs > ylow

#can also put constraint on ycs



#bool_array_comb = bool_array*bool_array2*bool_array3
#including in criteria on y position
bool_array_comb = bool_array*bool_array2*bool_array3#*bool_array4*bool_array5

myind = np.nonzero(bool_array_comb)

select_ind = ind[myind]
select_proc = proc[myind]
select_len = np.size(select_ind)


#print(select_proc)
print('selecting ',np.size(select_ind),' particles')

myprt.close()
myfld.close()
'''
x = x[myind]
y = y[myind]
mylen = np.size(x)
#testing the overplotting to make sure our criteria are working
dens = myfld['dens'][0,:,:]
fig, ax = plt.subplots()
ax.imshow(dens/4.,origin='lower',vmin=0,vmax=4)
print(mylen)

for i in range(mylen):
    myx = x[i] / istep #in output units, let's just do that for now
    myy = y[i] / istep
    circle = patches.Circle((myx,myy),radius=.2,fill=True,color='red',alpha=.5)
    ax.add_patch(circle)

plt.savefig('testing_select.png',dpi=300,bbox_inches='tight')
'''
#now we loop through time steps, identify particles based off their ID, and plot just these particles we want
t_list = []
gam_list = []
EdotB_list = []

avg_gam_list = []
avg_EdotB_list = []

for i in range(np.size(select_ind)):
    gam_list.append([])
    EdotB_list.append([])

for t in range(t0, tf):
    t_list.append(t)
    tstr = '%03d' % t
    print(tstr)
    myfld = mybase + "flds.tot." + tstr
    myprt = mybase + "prtl.tot." + tstr
    myfld = h5py.File(myfld,'r')
    myprt = h5py.File(myprt,'r')
    dens = myfld['dens'][0,:,:]
    ind = np.array(myprt['inde'])
    proc = np.array(myprt['proce'])
    gammae = np.array(myprt['gammae'])
    mylen = np.size(ind)
    x = myprt['xe']
    y = myprt['ye']

    exe = np.array(myprt['exe'])
    eye = np.array(myprt['eye'])
    eze = np.array(myprt['eze'])
    

    bxe = np.array(myprt['bxe'])
    bye = np.array(myprt['bye'])
    bze = np.array(myprt['bze'])

    EdotB = exe*bxe + eye*bye + eze*bze

    #fig,(ax1,ax2) = plt.subplots(1,2)
    ax1 = plt.subplot2grid((1,2),(0,0),rowspan=1,colspan=1)
    ax2 = plt.subplot2grid((1,2),(0,1),colspan=1)

    xscan = 150
    xlen = np.shape(dens)[1]
    xhlf = xlen/2
    xup = int(xhlf + xscan)
    xlow = int(xhlf - xscan)
    dens = dens[:,xlow:xup]
    xlow_compunits = xlow *istep / c_omp / 1000

    xext = istep*xscan / c_omp / 1000
    yext = istep * np.shape(dens)[0]/c_omp / 2. / 1000


    ax1.imshow(dens/4.,vmin=0,vmax=5,origin='lower', extent = [-xext, xext, -yext, yext])

  
    ytot = 2*yext
    rect_loc = yfrac_hlf * ytot - yext - yscan_frac*ytot

    rect_h = yscan_frac * ytot * 2

  
    #rectangle = patches.Rectangle((-xext, rect_loc),2*xext, rect_h,color="White",fill=False)
    #ax1.add_patch(rectangle)
    #show rectangle around region we're selecting
  

    ax1.set_xlabel('$x \; (1000 \; c/\omega_{p})$')
    ax1.set_ylabel('$y \; (1000 \; c/\omega_{p})$')
    ax3 = ax2.twiny()
    #print('iterating through ', mylen)
    #print('shape of ind : ',np.shape(ind)) 
    #indin = np.in1d(ind, select_ind) #returns boolean in shape time_ind 
    #procin = np.in1d(proc, select_proc)

    #final_select = indin*procin #this line doesn't seem to be doing what I want it, seeing duplicate indices which should be zeroed out by the processor number..

    #print(select_ind)
    #print(select_proc)

    #print(ind[final_select])
    #print(proc[final_select])
    #print('picked out ',np.sum(final_select), ' particles')
    #now find indices
    #final_select = np.nonzero(final_select)
    #print(np.size(final_select))
    
    #need to first loop through selected indices, and then search for indeces in the larger array where ind[i] and proc[i] match the pre-specified ones.
    avg_gam = 0
    avg_EdotB = 0
    spectra_list = []
    edotb_spectra_list = []
    for i in range(select_len):
        #print(myind, myproc)
        myind = select_ind[i]
        myproc = select_proc[i]
        for j in range(mylen):
            if ind[j]==myind and proc[j]==myproc:
                #not sure if order is preserved for tracking indidual particles through time here but let's try
                mygam = gammae[j]
                avg_gam += mygam
                myEdotB = EdotB[j]
                avg_EdotB += myEdotB
                gam_list[i].append(mygam)
                EdotB_list[i].append(myEdotB)
                spectra_list.append(mygam)
                edotb_spectra_list.append(myEdotB)
                #ax2.plot(t_list, gam_list[i], color="Red",linewidth=0.7)
                #ax2.scatter(t,mygam,color="Red",s=1)
                #ax3.plot(t_list, np.abs(EdotB_list[i]),color="Blue",linewidth=0.3)

                #print(myind, myproc)
        #now loop through prtls, should we zero out all particle w/o either correct proc or ind numbers first from array to save time?  Let's see how long this takes
                #mycount +=1
        #if ind[i] in select_ind and proc[i] in select_proc:
            #print(ind[i],proc[i])
                #myx = x[j] / istep #in output units, let's just do that for now    
                #myy = y[j] / istep   
                
                #let's put myx and myy in the units we want to present

                #commenting out overplotting particles
                '''
                myx = x[j] / istep - xlow
                myx *= istep / c_omp / 1000 
                myx -= xext
                
                myy = y[j] / c_omp / 1000 - yext
                
                #print(myx, myy)
                circle = patches.Circle((myx,myy),radius=.002,fill=True,color='red')

                #myx = x[j]/c_omp/1000 #into skin depths
                #myy = x[j]/c_omp/1000

                #print(myx, myy)
                #circle = patches.Circle((myx-xlow_compunits-xext,myy-yext),radius=.0002,fill=True,color='red') 
                ax1.add_patch(circle)
                '''
                break
    
    #make spectra out of spectra list
    if len(spectra_list) > 0:
        bin_num = 35
        #gamlow = min(spectra_list)
        #gamup = max(spectra_list)
        gamlow = 1
        gamup = 5e3
        histbins, merger_hist = return_spec(np.array(spectra_list),gamlow,gamup,bin_num)

    
    ax3.hist(edotb_spectra_list,bins='auto',histtype='step',label='$E \cdot B$')
    #ax2.hist(spectra_list,color="C2",histtype='step',label="Particle Spectra", bins = "auto")
    ax2.plot(histbins, merger_hist,color="C1",label="Particle Spectrum")


    #avg_gam /= select_len
    #avg_EdotB /= select_len
    #avg_gam_list.append(avg_gam)
    #avg_EdotB_list.append(avg_EdotB)
    #ax2.plot(t_list, avg_gam_list,color="Cyan",label='$\\bar{\gamma}$')
    #ax3.plot(t_list, np.abs(avg_EdotB_list), color="Blue",label='$|\overline{E \cdot B}|$')
    ax2.legend(loc = 'upper left',frameon=False)
    ax3.legend(loc= 'upper right',frameon=False)
    #print('plotted ' ,mycount,' particles')
    #ax2.set_xlabel('Time')
    #ax2.set_xlim(t0,tf)
    ax2.set_xlabel('$\gamma$')
    ax3.set_xlabel('$E \cdot B$')
    #ax3.set_ylabel('$|E \cdot B|$')
    ax3.set_xlim(-1,1)
    #ax2.plot(histbins, merger_hist,color="C1",label="Particle Spectra")

    #ax2.set_ylabel('$\gamma$')
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_xlim(1,1e4)
    ax2.set_ylim(1,1e2)
    #ax3.set_ylim(0,.25)
    plt.savefig('select_plots_spec/triggered_bguide/'+tstr+'.png',dpi=300,bbox_inches='tight')
    plt.close()
    myfld.close()
    myprt.close()
    

