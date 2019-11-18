import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

from tristan_funcs import get_val_filename, return_spec
from tristan_funcs_vecpot import vecpot2
from localmin_func import localmin_wrap
#from localmin_func import localmin, localmax, localmin_wrap                                
from new_xpoint_func import return_xpoints
import h5py

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

plt.rcParams['figure.figsize'] = 10,4


#bguide triggered
mybase = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/"
#mybase = "../../tristan_acc-mec_Ez/bguide.3_untriggered_stride1_gamthresh50/output/"

#mybase = "../../tristan_acc-mec_Ez/bguide.1_untriggered_stride1_gamthresh50/output/"

#mybase = "/home/u21/davidrball/david/tristan_acc-mec_Ez/triggered_bguide_hightimeres_begintrack/output/"


spec_binary=0 #whether you want to plot / save spectra
second_xpoint_col="red" #different color for secondary x-point finding just to diagnose problems

fldbase = mybase +"flds.tot."
prtlbase = mybase + "prtl.tot."

#savename = "xpoint_plots/hightimeres/"
#savename = "final_snapshot_bguide1_untriggered"
savename = "final_plots/bguide1_untriggered_snapshot"


#don't need this if we're not doing the spectra saving
spect_filepath = "bguide.3_triggered_hightimeres"
filepath = "classification_txt/bguide.3_triggered_hightimeres.txt"  

t0=9
tf = 65
#t_list = range(t0,tf)
t_list = [9,33,46]
#t_list = [7,20,30]
#sim params
interval=2000
istep=12
c_omp=3
va=.48
c=.45

xscan = (interval/2)*va*c / (c_omp*1000)
xscan = .072 #alfven connectivity in plot units (1000 c/omp) for case of interval=2000




#xscan=.32 #c connectivity in case of high time res

grid_to_plot_coord = istep/c_omp/1000.


#xpoint finding params
edge=4
smoothlen=0
myroll=1
bdens_tolerence=4
zero_tolerence=3

#first we loop through all fld files, make lists of x-points
time_xpoint_list = []
xpoint_spec_list = []
merger_spec_list = []
else_spec_list = []


tfinal_str = "%03d" % tf
prtl_final = prtlbase + tfinal_str
myf_prtl = h5py.File(prtl_final,'r')

gamma = np.array(myf_prtl['gammae'])

ycs = np.array(myf_prtl['ycse'])/(1000*c_omp) #3 cells per c_omp



tcs = np.array(myf_prtl['tcse'])
#convert tcs down to output times
tcs = np.round(tcs / float(interval),decimals=0)
ezbxy = np.array(myf_prtl['ezbxye'])
prtnump = np.size(tcs)
print('total particles : ',prtnump)

scan=5

myf_prtl.close()
#t_list = [t_list[0]]


fig,((ax0,ax1,ax2),(ax3,ax4,ax5)) = plt.subplots(2,3,sharex=True)

#fig, ax = plt.subplots(2,3)
#ax0, ax1, ax2, ax3, ax4, ax5=ax.flatten()
axlist = [[ax0, ax1],[ax2,ax3],[ax4,ax5]]
fldaxlist = [ax0,ax1,ax2]
scataxlist = [ax3,ax4,ax5]


mycount = 0
for t in t_list:
    print(t)
    tstr = "%03d" % t
    myfld = fldbase + tstr
    myfld = h5py.File(myfld,'r')
    #max_list = axlist[mycount]
    #fldax = max_list[0]
    #scatax = max_list[1]
    fldax = fldaxlist[mycount]
    scatax = scataxlist[mycount]


    xpoint_loc, leftcut, rightcut, bottomcut, topcut = return_xpoints(myfld,edge,smoothlen,myroll,bdens_tolerence,zero_tolerence)
    
    #using a secondary method to ensure we don't miss xpoints
    vecpot,junk1,junk2 = vecpot2(myfld,istep,c_omp)
    xmid = np.shape(vecpot)[1]/2
    vecpot_slice = vecpot[:,xmid]
    bdens = myfld['bdens'][0,:,:]
    xmid = np.shape(bdens)[1]/2
    bdens_slice = bdens[:,xmid]
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
    vecpot_slice = vecpot_slice[:modified_index]
    xpoint_list_ind = localmin_wrap(vecpot_slice,scan)
 



    xpoint_dellist = []
    for index in range(len(xpoint_list_ind)):
        xpoint_loc_tmp = xpoint_list_ind[index]
        if bdens_slice[xpoint_loc_tmp] > bdens_tolerence:
            #print('bdens criterion triggered')                                             
            #print(xpoint_loc)                                                              
            xpoint_dellist.append(index)


        #just for plotting, doesn't make an actual difference for spectra
  
    xpoint_list_ind = np.array(xpoint_list_ind)
    xpoint_list_ind = np.delete(xpoint_list_ind, xpoint_dellist)

    print('number of secondary x-points : ', np.size(xpoint_list_ind))
    
    #plotting for sanity check
    #fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)
    plt.subplots_adjust(hspace=0,wspace=0)
    dens = myfld['dens'][0,bottomcut:topcut,leftcut:rightcut]
    dens_shape = np.shape(dens)
    xlen = dens_shape[1]
    xext = xlen * grid_to_plot_coord
    yext = dens_shape[0] * grid_to_plot_coord
    fldax.imshow(np.rot90(dens), origin='lower',vmin=0,vmax=20,extent=[-yext/2,yext/2,-xext/2,xext/2])
    
    xpoint_list2 = xpoint_list_ind * istep / c_omp / 1000 - yext/2.

    #just for plotting local min xpoints which won't account for inclined current layer
    '''
    if mycount ==1:
        for j in range(len(xpoint_list2)):
            myx = xpoint_list2[j]
            if myx < 0:
                myy=.02
        #if myx < -.2: 
        #    myy=-.015
            elif myx > 0:
                myy=-.01
            elif myx == 0:
                myy=0        
            fldax.scatter(myx, myy,color=second_xpoint_col,marker='x',s=10)
    
     else:
        fldax.scatter(xpoint_list2-.01, np.zeros(np.size(xpoint_list2)),color=second_xpoint_col,marker='x',s=10)
    #print('xext : ',xext, 'yext : ',yext)
    '''
    fldax.scatter(xpoint_list2-.01, np.zeros(np.size(xpoint_list2)),color=second_xpoint_col,marker='x',s=10)  


    
    print('NUMBER OF 2D XPOINTS FOUND : ', np.shape(xpoint_loc))
    for sublist in xpoint_loc:

        x = sublist[0]*grid_to_plot_coord-yext/2
        y = sublist[1]*grid_to_plot_coord-xext/2
        print(x, y)
        #if x < 1:#just for plotting temporarily
        #    ax1.scatter(x,y,color="red",s=10,marker='x')
        #ax1.scatter(x,-y,color="white",s=10,marker='x')
    xpoint_xvals = []
    #for sublist in xpoint_loc:
    #    xpoint_xvals.append(sublist[0]*grid_to_plot_coord-yext/2)

    
    for xpoint2 in xpoint_list2:
        xpoint_xvals.append(xpoint2)
    
    print('xpoint vals')
    print(xpoint_xvals)
    myind = tcs==t #boolean array of prtls that were accelerated near this timestep                                                                          
    myind = np.nonzero(myind)#indices of prtls accelerated                     
    #slicing off the only particles we care about at this timestep             
    tmp_ycs = ycs[myind]-yext/2
    tmp_gamma = gamma[myind]
    tmp_ezbxy = ezbxy[myind]

    #unfortunately need to do loop to associate prtls w/ xpoints etc.
    prtnum = np.size(tmp_ycs)
    col_list = []
    ycs_xpoint = []
    ycs_merger = []
    ycs_else = []
    gam_xpoint = []
    gam_merger = []
    gam_else = []

    #to handle boundaries let's use a little hack
    #loop through xpoints, if xpoint - xscan < -yext/2, just add an xpoint at yext + overflow where overflow is how much past the boundary it went
    addlist = []
    for xpoint_xval in xpoint_xvals:
        if xpoint_xval -xscan < -yext/2:
            overflow = abs(-yext/2 - (xpoint_xval-xscan))
            addlist.append(yext/2 + overflow)

    for tmpval in addlist:
        xpoint_xvals.append(tmpval)

    for i in range(prtnum):
        my_ycs = tmp_ycs[i]     
        counted = 0
        for xpoint_xval in xpoint_xvals:
            if np.abs(xpoint_xval - my_ycs) < xscan:
                if tmp_ezbxy[i]<0:
                    ycs_xpoint.append(my_ycs)
                    gam_xpoint.append(tmp_gamma[i])
                    counted =1
                elif tmp_ezbxy[i]>0:
                    ycs_merger.append(my_ycs)
                    gam_merger.append(tmp_gamma[i])
                    counted = 1
                break
            '''
            if xpoint_xval - xscan < -yext/2: #handling boundaries
                #find how far over the boundary xscan goes:
                overflow = yext/2 - np.abs(xpoint_xval-xscan)
                upperlim = yext/2 - overflow
                if my_ycs > upperlim:
                    if tmp_ezbxy[i]<0:
                        ycs_xpoint.append(my_ycs)
                        gam_xpoint.append(tmp_gamma[i])
                        counted =1
                        break
                    elif tmp_ezbxy[i]>0:
                        ycs_merger.append(my_ycs)
                        gam_merger.append(tmp_gamma[i])
                        counted = 1
                        break
            '''
        if counted ==0:
            ycs_else.append(my_ycs)
            gam_else.append(tmp_gamma[i])
            
        
    #scatax.scatter(ycs_else,gam_else,s=1,c='Red',label="Other")
    mysize = .25
    scatax.scatter(ycs_xpoint,gam_xpoint,s=mysize,c='Blue',label="X-point (PCS)") 
    scatax.scatter(ycs_merger,gam_merger,s=mysize,c='Orange',label="X-point (MCS)") 
    scatax.scatter(ycs_else,gam_else,s=mysize,c='Red',label="Other")

    ax3.legend(frameon=False,prop={'size':8},loc="upper left")
    scatax.set_yscale('log')
    scatax.set_ylim(50,1e4)
    scatax.set_xlim(-yext/2,yext/2)


    tann_str = str(round(t*.15*2000))
    annstr = "t="+tann_str+' $\omega_{p}^{-1}$'
    fldax.annotate(annstr,(-.4,.5),color="White")
    fldax.set_xticks([])
    scatax.set_xlabel('$x \; (1000 \; c/\omega_{p})$')
    #scatax.set_ylabel('$\gamma_{\\rm{f}}$')
    scatax.set_xticks([-1,-.5,0,.5,1])
    if mycount==0:
        fldax.set_ylabel('$y \; (1000 \; c/\omega_{p})$')
        fldax.set_yticks([-.5, 0, .5])
        scatax.set_ylabel('$\gamma_{\\rm{f}}$')
        scatax.set_yticks([1e2,1e3,1e4])
    if mycount != 0:
        fldax.set_yticks([])
        scatax.set_yticks([])

    #plt.savefig(savename+tstr+'.png',bbox_inches='tight',dpi=300)
    #plt.savefig(savename+tstr+'.pdf',bbox_inches='tight',dpi=300)
    #plt.close()
    
    #add to total spectra lists
    
    for gam in gam_xpoint:
        xpoint_spec_list.append(gam)
    for gam in gam_merger:
        merger_spec_list.append(gam)
    for gam in gam_else:
        else_spec_list.append(gam)
    myfld.close()
    mycount += 1


plt.savefig('bguide3_accmec_singleplot.png',dpi=300,bbox_inches='tight')
plt.savefig('bguide3_accmec_singleplot.pdf',dpi=300,bbox_inches='tight')

