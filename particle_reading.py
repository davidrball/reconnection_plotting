import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from tristan_funcs import get_val_filename, load_particle
import h5py
import matplotlib.patches as patches
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 10})



partnumstart = 6000
partnumstop = 6002
mypath = 'particle_txt/sig.3_delgam0005/bguide.3_earlytime/'
c_omp = 3
istep = 12

sigma=.3
va = np.sqrt(sigma/(1+sigma))

t_start = 1
t_stop = 25
fld_base = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/flds.tot."

for partnum in range(partnumstart, partnumstop):
    d = load_particle(mypath, partnum)
    t_list = d['time']
    #print(t_list)
    x_list = np.array(d['x'])
    y_list = np.array(d['y'])
    u_list = np.array(d['u'])
    v_list = np.array(d['v'])
    w_list = np.array(d['w'])
    bx_list = np.array(d['bx'])
    by_list = np.array(d['by'])
    bz_list = np.array(d['bz'])
    ex_list = np.array(d['ex'])
    ey_list = np.array(d['ey'])
    ez_list = np.array(d['ez'])
    gamma_list = np.array(d['gamma'])
    #print(gamma_list)
    maxgam = max(gamma_list)
    mingam = min(gamma_list)


    #print(t_list)

    for t in range(t_start, t_stop):
        tstr = str(t)
        print(tstr)
        if len(tstr) == 2:
            myfld = fld_base + "0"
        elif len(tstr) == 1:
            myfld = fld_base + "00"
        myfld += tstr
        print(myfld)
        f_fld = h5py.File(myfld)
        
        #pick out prtl values at the correct time:
        for index in range(len(t_list)):
            if t_list[index] == t:
                i = index
                break
        #to extract particle values at this time, use [index]
        x = d['x'][i]
        y = d['y'][i]
        u = d['u'][i]
        v = d['v'][i]
        w = d['w'][i]
        bx = d['bx'][i]
        by = d['by'][i]
        bz = d['bz'][i]
        ex = d['ex'][i]
        ey = d['ey'][i]
        ez = d['ez'][i]
        gamma = d['gamma'][i]

        gamcol = np.log10(gamma/maxgam) / np.log10(mingam/maxgam)
        print(gamcol)

        intx = int(x)
        inty = int(y)
        intx /= istep
        inty /= istep

        dens = np.array(get_val_filename("dens",f_fld))[0]
        
        
        fig = plt.figure()
        ax1 = plt.subplot2grid((1,2),(0,0),rowspan=1,colspan=1)
        ax2 = plt.subplot2grid((1,2),(0,1),colspan=1)
        xhlf = np.shape(dens)[1]/2
        xscan = 100
        xlow = int(xhlf - xscan)
        xup = int(xhlf + xscan)
        xext = istep*2*xscan / c_omp # in electron skin depths                      
        yext = istep * np.shape(dens)[0] /c_omp
        ax1.imshow(dens[:,xlow:xup],origin='lower',vmax=20, extent=[0,xext,0,yext])
        ax1.set_xlabel('$c/\omega_{p}$',size=14)
        ax1.set_ylabel('$c/\omega_{p}$',size=14)


        ax2.plot(t_list[:i],gamma_list[:i],color="Red",label="$\gamma$")
        #ax2.scatter(t,gamma,color="Red")
        plt.legend(loc='lower left',prop={'size':12},frameon=False)

        ax2.set_ylim(1e0,1e4)
        ax2.set_yscale('log')
        ax2.set_xlim(0,70)
        ax2.set_xlabel('Time $(300\omega_{p}^{-1})$')
        ax2.set_ylabel('$\gamma$')
        ax3 = ax2.twinx()
        ax3.plot(t_list[:i], ez_list[:i]/(va*np.sqrt(bx_list[:i]**2+by_list[:i]**2)), color = "Blue",label="$E_{z}/B_{xy}$")
        ax3.set_ylim(-3,3)
        #ax3.set_ylim(0,4)
        ax3.set_ylabel('$E_{z}/(V_{a} B_{xy})$')
        ax3.set_xlim(t_start-1, t_stop+1)
        intx -= xlow
 
        intx *= istep/c_omp
        inty *= istep/c_omp
          
        circle = patches.Circle((intx,inty),radius=20,fill=False,color=(gamcol, 1-gamcol,1))
        plt.tight_layout()
        ax1.add_patch(circle)
        plt.legend(loc='lower right',prop={'size':12},frameon=False)
        plt.savefig('plots/particle_plots/sig.3_delgam0005/bguide.3_earlytime/'+'E'+str(partnum)+'_'+str(t)+'.png',bbox_inches='tight',dpi=300)
        f_fld.close()
        plt.close()

