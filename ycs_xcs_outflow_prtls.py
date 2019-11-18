import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tristan_funcs import get_val_filename
import h5py
import numpy as np
import matplotlib.patches as patches

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 15})

myfld_base = "../../accmec_withy/bguide.3_triggered/output/flds.tot."

prtl_final = "../../accmec_withy/bguide.3_triggered/output/prtl.tot.020"
prtl_base = "../../accmec_withy/bguide.3_triggered/output/prtl.tot."

f_prtl_final = h5py.File(prtl_final, 'r')



tcse = np.array(get_val_filename('tcse',f_prtl_final))
ycse = np.array(get_val_filename('ycse',f_prtl_final))
xcse = np.array(get_val_filename('xcse',f_prtl_final))
gammae = np.array(get_val_filename('gammae',f_prtl_final))
ez_bxy = -np.array(get_val_filename('ezbxye',f_prtl_final))

c=.45
rstep_first = 12000
rstep_jump = 2000
rstep_interval = 2000
u_sh = .1
btheta=0.3
sigma=0.3

tcse_jump = (tcse-rstep_first)//rstep_interval+1
tcse_jump[tcse-rstep_first==0]=1
tcse_jump[tcse_jump<0]=0

#xcse_adjusted = xcse - c*tcse
corspeed = 1-u_sh*np.sqrt(sigma/(1+sigma*(1+btheta**2)))
print('corspeed : ',corspeed)
#xcse_jump = xcse_adjusted + c*corspeed*rstep_jump*tcse_jump
prtnum = np.size(tcse)

gam_max_final = np.max(gammae)

f_prtl_final.close()
final_len = np.size(tcse)

sig=0.3
va = np.sqrt(sig/(sig+1))
istep = 12
c_omp = 3
interval = 2000.

convfac = istep/c_omp / 1000.

t_start = 15
t_stop = 16

fld_scan = 100

yext = fld_scan * istep / (1000*c_omp)

tscan = .5

c_connect = va*tscan*interval * c #this will be in computational cells

c_connect *= convfac #converting computational units to plotting units
print('c connect : ',c_connect)

for t in range(t_start,t_stop):
    t_string = str(t)
    if len(t_string) == 2:
        fld_base =myfld_base+ "0"
      #  myprtl = prtl_base + "0"
    elif len(t_string)==1:  
        fld_base =myfld_base+ "00"
     #   myprtl = prtl_base + "00"
    fld_base += t_string
    #myprtl += t_string

    f_fld = h5py.File(fld_base,'r')
    #f_prtl = h5py.File(myprtl,'r')
    dens = np.rot90(get_val_filename('dens',f_fld)[0,:,:])
    ymid = np.shape(dens)[0]/2
    ylow = int(ymid - fld_scan)
    yup = int(ymid + fld_scan)
    dens = dens[ylow:yup,:]
    #dens = np.fliplr(dens)
    newmid = np.shape(dens)[0]/2
    xext = np.shape(dens)[1]*istep / (1000*c_omp)
    
    fig, ax = plt.subplots(1)


    im = ax.imshow(dens/4.,origin='lower',cmap='viridis',vmin=1,vmax=4, extent= [-xext/2, xext/2, -yext, yext])
    ax.set_xlabel('$x \; (1000 \; c/\omega_{p})$')
    ax.set_ylabel('$y \; (1000 \; c/\omega_{p})$')


    #[.125,.7,.77,.03]
    cbar_ax = fig.add_axes([.125, .7, .77, .03])
    cb = fig.colorbar(im,cax=cbar_ax, orientation="horizontal")
    cbar_ax.set_xlabel("Density (particles per cell / $N_{ppc})$",fontsize=14,labelpad=-45)
    #fld info plotted
    
    #now loop through particles, for any with tcs/interval between t-1 and t, add plot them over the fld quantities

    f_fld.close()
    count = 0
    gam_list = []
    ycs_list = []
    xcs_list = []
    ez_bxy_list = []
    x_list = []
    tlow = t-tscan
    tup = t+tscan
    print(tlow, tup)

    
    #expansion = interval*t*c

    for i in range(prtnum):
        mycs = ycse[i] #in units of cells
        mytcs = tcse[i]
        #myxcs = xcse_jump[i]
        mygam = gammae[i]
        myxcs = xcse[i]
        #lap = mytcs
        lap = mytcs
        if mycs != 0 and mytcs != 1 and 90<mygam<110:
            #in actual computational steps
            mytcs /= interval
            
            #print(tlow, mytcs, tup)
            if tlow < mytcs <= tup:
                #myx = xe[i]/istep
                #print('tcse : ',tcse[i])
                #print('tcse jump : ', tcse_jump[i])
                
                count +=1
                #mycs /= istep
                mycs /= (c_omp * 1000)
                mycs -= xext/2
                
                #account for moving injector
                #print(lap)
                if lap < rstep_first:
                    jumpsteps = 0
                    #print('jump0')
                elif lap == rstep_first:
                    jumpsteps=1
                    #print('jump1')
                elif lap > rstep_first:
                    jumpsteps = lap // rstep_jump 
                

                jumpdistance = jumpsteps*rstep_interval
                #myxcs -= c*mytcs*interval
                #myxcs -= c*corspeed*jumpdistance
                expansion = c*corspeed*mytcs*interval / (c_omp*1000)
                #print('expansion : ',expansion)
                contraction = c*corspeed*jumpdistance / (c_omp*1000)
                #print('contraction : ', contraction)
                #print('expansion - contraction : ', expansion - contraction)
                #xoffset = expansion-contraction
                #xoffset = 2.46 - (t-49)*(1-corspeed)*2.46
                xoffset=2.455
                #print('working offset : ', xoffset)
                myxcs /= (c_omp*1000)
                myxcs -= yext
                myxcs -= xoffset
                if mycs < -xext/2.:
                    print('triggered minus')
                    mycs = -xext / 2+.01
                if mycs > xext/2.:
                    print('triggered plus')
                    mycs = xext/2.-.01
                ycs_list.append(mycs)
                xcs_list.append(myxcs)
                #print(mycs)
                #x_list.append(myx-ylow)              
                my_ezbxy = (ez_bxy[i])/va
                #print(my_ezbxy)
                mygam = gammae[i]
                loggam = np.log10(mygam)
                gam_list.append(mygam)
                ez_bxy_list.append(my_ezbxy)
                #ax.scatter(mycs, -.25,s=((mygam/200.)**2.5), c=my_ezbxy,vmin=0,vmax =4 , cmap='YlOrRd')
    #loggam_arr = np.log10(np.array(gam_list))
    size_arr = (np.array(gam_list)/350)**3

    #disable scattering
    
    
    im2 = ax.scatter(ycs_list, xcs_list, s=.15, c=ez_bxy_list,vmin=-1,vmax=1,cmap='coolwarm')
    
    #im3 = ax.scatter(0,0,s=100,edgecolors='white',facecolors='none',linestyle='--')

    circle = patches.Circle((-.66,.017),fill=False,color="White",radius=c_connect)
    #ax.add_patch(circle)
    cbar_ax2 = fig.add_axes([.92, .34, .02, .31])
    cb2 = fig.colorbar(im2,cax=cbar_ax2, orientation="vertical")
    cbar_ax2.set_ylabel("$E_{z}/B_{xy}$",fontsize=14)
    cbar_ax2.yaxis.set_label_position('right')
    cb2.set_ticks([-1,-.5,0,.5,1])
    
    #plt.scatter(ycs_list, x_list, s=1, c=loggam_arr,vmin=0,vmax=np.log10(np.max(gam_final)),cmap='viridis')
    
    if count == 0:
        pass
    else:
        print('particles plotted : ',count)
        print('max gam : ', max(gam_list))
        print('max of all particles : ', gam_max_final)
 
    
    #plt.savefig('fld_acc_plots/'+t_string,dpi=300,bbox_inches='tight')
    #plt.savefig('prtl_acc_xy_noscatter.png',dpi=300,bbox_inches='tight')
    plt.savefig('prtl_outflow_acc.png',dpi=300,bbox_inches='tight')
    plt.close()
