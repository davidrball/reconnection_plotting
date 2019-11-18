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

partnumstart = 2
partnumstop = 15
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

    edotv_par_cumulative = []
    edotv_perp_cumulative = []
    edotv_tot_cumulative = []
    for i in range(np.size(eperp_uperp)):
        tmp_par = trapz(epar_upar[0:i])
    #tmp = trapz(edotv_tot[0:i])
        tmp_perp = trapz(eperp_uperp[0:i])
        edotv_par_cumulative.append(tmp_par)
        edotv_tot_cumulative_tmp = trapz(edotv_tot[0:i])
        edotv_perp_cumulative.append(tmp_perp)
        edotv_tot_cumulative.append(edotv_tot_cumulative_tmp)

    mynorm = gamma_list[-1]/edotv_tot_cumulative[-1]
    print('mynorm : ',mynorm)
    

    '''
    #for two separate axes
    fig,ax1 = plt.subplots()
    ax2=ax1.twinx()
    #mynorm = 1
    ax1.plot(newt_list,mynorm*np.array(edotv_par_cumulative),color="Red",label="$\int v_{||}E_{||}dt$")
    ax1.plot(newt_list,mynorm*np.array(edotv_perp_cumulative),color="Blue",label="$\int v_{\\bot}E_{\\bot}dt$")
    ax2.plot(newt_list,gamma_list-1,color="Black",label='$\gamma-1$')
    ax2.set_yscale('log')
    ax1.legend(loc='upper left')
    ax2.legend(loc='lower right')
    ax1.set_xlabel('Time (output timesteps)')
    ax1.set_ylabel('$\\vec{E} \\cdot \\vec{v}$')
    '''
    fig, ax1 = plt.subplots()
    ax1.plot(newt_list,mynorm*np.array(edotv_par_cumulative),color="Red",label="$\int v_{||}E_{||}dt$")                                                       
    ax1.plot(newt_list,mynorm*np.array(edotv_perp_cumulative),color="Blue",label="$\int v_{\\bot}E_{\\bot}dt$")                                               
    ax1.plot(newt_list,gamma_list,color="Black",label='$\gamma-1$')         
    #ax1.set_yscale('log')                                                     
    ax1.legend(loc='upper left')                                              
    #ax2.legend(loc='lower right')                                             
    ax1.set_xlabel('Time (output timesteps)')                                 
    #ax1.set_ylabel('$\\vec{E} \\cdot \\vec{v}$')     
    ax1.set_ylabel('$\Delta \gamma$')
    plt.savefig('vdote_plots/untriggered_latetime/cumulative/'+str(partnum)+'.png',dpi=300,bbox_inches='tight')
    
