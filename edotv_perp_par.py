import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from tristan_funcs import get_val_filename, load_particle_edotv
from scipy.integrate import trapz
import h5py
import matplotlib.patches as patches
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

c=.45
interval = 2000

partnumstart = 1
partnumstop = 2
mypath = 'edotv_prtls/bguide.3_untriggered_hightimeres_latetime/'
#mypath = 'edotv_prtls/bguide.3_triggered_hightimeres_nothresh/'

for partnum in range(partnumstart, partnumstop):
    d = load_particle_edotv(mypath, partnum)
    t_list = d['time']
    #print(t_list)
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


    bmag = np.sqrt(bx_list**2 + by_list**2 + bz_list**2)
    
    unitbx = bx_list / bmag
    unitby = by_list / bmag
    unitbz = bz_list / bmag
                
    #unit vector in b direction is (unitbx, unitby, unitbz)


    u_list /= gamma_list
    v_list /= gamma_list
    w_list /= gamma_list

    upar_list = u_list*unitbx + v_list*unitby + w_list*unitbz #u dot b
    epar_list = ex_list*unitbx + ey_list*unitby + ez_list*unitbz

    epar_upar = -upar_list*epar_list


    len_t = np.size(t_list)
    len_vlist = np.size(u_list)
    diff = len_t - len_vlist
    newt_list = t_list[diff:]
    edotv_list = []
    edotv_list_x = []
    edotv_list_y = []
    edotv_list_z = []
    edotv_x = -u_list*ex_list
    edotv_y = -v_list*ey_list
    edotv_z = -w_list*ez_list

    edotv_tot = edotv_x + edotv_y + edotv_z

    eperp_uperp = edotv_tot - epar_upar

    deltat = 1000 #inverse plasma frequencies, but what unit of time is appropriate?

    edotv_cumulative = []
    for i in range(np.size(edotv_x)):
        tmp = np.sum(edotv_tot[0:i])*deltat
    #tmp = trapz(edotv_tot[0:i])
        edotv_cumulative.append(tmp)

#mynorm = gamma_list[-1]/edotv_cumulative[-1]
#print('mynorm : ',mynorm)
    '''
    fig,ax1 = plt.subplots()
    ax2=ax1.twinx()
    mynorm = 1
    ax1.plot(newt_list,mynorm*np.array(edotv_cumulative),color="Red")
    ax2.plot(newt_list,gamma_list-1,color="Black")
    ax2.set_yscale('log')
    plt.savefig('testing_cumulative.png',dpi=300,bbox_inches='tight')
    '''
    firsttime = 1
    lasttime = 300


    #mynorm = gamma_list[-1] / (np.sum(eperp_uperp) + np.sum(epar_upar))
    mynorm = 50

    fig, ax1 = plt.subplots()
    #ax1.plot(newt_list[:lasttime],edotv_x[:lasttime],color="Blue", label='$q_{e}v_{x}E_{x}$')
    #ax1.plot(newt_list[:lasttime],edotv_y[:lasttime],color="Red", label='$q_{e}v_{y}E_{y}$')
    #ax1.plot(newt_list[:lasttime],edotv_z[:lasttime],color="Green", label='$q_{e}v_{z}E_{z}$')
    ax1.plot(newt_list[firsttime:lasttime],mynorm*eperp_uperp[firsttime:lasttime],color="Blue",label="$E_{\\bot}v_{\\bot}$")
    ax1.plot(newt_list[firsttime:lasttime],mynorm*epar_upar[firsttime:lasttime],color="Red",label="$E_{||}v_{||}$")
    ax1.plot(newt_list[firsttime:lasttime],mynorm*edotv_tot[firsttime:lasttime],color="Black", label = '$q_{e}\\vec{v} \\cdot \\vec{E}$',linewidth=.7,linestyle='dotted')
    #ax1.set_ylim(0,.2)

    ax2 = ax1.twinx()
    ax2.plot(newt_list[firsttime:lasttime],gamma_list[firsttime:lasttime],color="Magenta",label='$\gamma$')
    ax2.set_yscale('log')
    #ax2.set_ylim(1e0,1e4)
    ax1.legend(loc='upper left')
    ax2.legend(loc='lower right')
    ax1.set_ylabel('$\\vec{v}\\cdot\\vec{E}$')
    ax2.set_ylabel('$\gamma$')
    ax1.set_xlabel('Time (output timesteps)')
    plt.savefig('vdote_plots/untriggered_latetime/perp_par/'+str(partnum)+'.png',bbox_inches='tight',dpi=300)
    #plt.savefig('testing_eparvpar2.png',dpi=300,bbox_inches='tight')

'''
fig, ax1 = plt.subplots()                                                                       
ax1.plot(newt_list,u_list,color="Blue", label='$q_{e}v_{x}E_{x}$')                             
ax1.plot(newt_list,v_list,color="Red", label='$q_{e}v_{y}E_{y}$')                              
ax1.plot(newt_list,w_list,color="Green", label='$q_{e}v_{z}E_{z}$') 
plt.savefig('testing_velocities.png')
'''


'''
totmom = np.sqrt(u_list**2 + v_list**2 + w_list**2)
plt.plot(t_list[diff:],u_list*ex_list,color="Red")
plt.plot(t_list[diff:],v_list*ey_list,color="Blue")
plt.plot(t_list[diff:],w_list*ez_list,color="Green")
#plt.plot(t_list[diff:],gamma_list,color="Black")
#plt.plot(t_list[diff:],totmom,color="Magenta")
plt.savefig('test_read_work.png')
plt.close()

'''
