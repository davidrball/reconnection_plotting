import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from tristan_funcs import get_val_filename, load_particle
import h5py
import matplotlib.patches as patches
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 18})



partnumstart = 1
partnumstop = 641
stride = 10


mypath = 'particle_txt/sig.3_delgam0005/bguide.1_earlytime/'
c_omp = 3
istep = 12

sigma=.3
va = np.sqrt(sigma/(1+sigma))

t_start = 1
t_stop = 25
fld_base = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/flds.tot."

y_acc_list = []
gam_final_list = []


for partnum in range(partnumstart, partnumstop,stride):
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
    ez_over_vabxy = ez_list / (va*np.sqrt(by_list**2 + bx_list**2))

#plt.plot(t_list, x_list,label='x') #need to reindex x values
#plt.plot(t_list, y_list, label='y')
#plt.plot(t_list, gamma_list,label='$\gamma$')

    gam0 = gamma_list[0]
    thresh = 10
    gamthresh = thresh*gam0
    gam_final = gamma_list[-1]
    for i in range(np.size(gamma_list)):
        mygam = gamma_list[i]
        if mygam > gamthresh:
            mymax = i
            break
    y_acc_list.append(y_list[mymax])
    gam_final_list.append(gam_final)
    #mymax given by the first time particle's lorentz factor exceeds 10x its initial values
    #at this point we want to record its y location along the sheet

plt.scatter(y_acc_list, gam_final_list)

plt.title('$\sigma=0.3$ $\\beta=0.003$ $B_{g}/B_{0}=0.1$')
plt.xlabel('$Y_{0}$')
plt.ylabel('$\gamma_{final}$')
plt.yscale('log')
plt.ylim(1e2,1e4)
plt.savefig('xpoint_acc_test_bguide.1.png',bbox_inches='tight')
