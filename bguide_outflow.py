import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
from matplotlib.colors import PowerNorm
from tristan_funcs import outflow

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})


outflow_list = [[],[],[],[],[],[]]
t_list = []

tstart = 1
tfinal = 24

#for high beta
fldbase = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide"
b0 = "0/output/flds.tot.0"
b1 = ".3/output/flds.tot.0"
b2 = ".5/output/flds.tot.0"
b3 = ".7/output/flds.tot.0"
bone = "1/output/flds.tot.0"
bg3 = "3/output/flds.tot.0"
'''
#for low beta
fldbase = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/"

fldbase0 = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/flds.tot.0"

b1 = "bguide.1/output/flds.tot.0"
b2 = "bguide.3/output/flds.tot.0"
b3 = "bguide.7/output/flds.tot.0"
'''
for t in range(tstart,tfinal):
    print(t)
    fld_base0 = fldbase + b0            
    #fld_base0 = fldbase0
    fld_base1 = fldbase + b1
    fld_base2 = fldbase + b2
    fld_base3 = fldbase + b3
    fld_base4 = fldbase + bone
    fld_base5 = fldbase + bg3

    fld_base_list = [fld_base0, fld_base1, fld_base2, fld_base3,fld_base4, fld_base5]
    t_list.append(t)


    t_str = str(t)
    file_list = []


    if len(t_str) == 1:
        for index in range(len(fld_base_list)):
            fld_base_list[index] += "0" + t_str
    elif len(t_str) == 2:
        for index in range(len(fld_base_list)):
            fld_base_list[index] +=  t_str

    print(fld_base_list)
    for fld_base in fld_base_list:
        file_list.append(h5py.File(fld_base,'r'))

                
    for index in range(len(fld_base_list)):
        myoutflow = outflow(file_list[index])
        outflow_list[index].append(myoutflow)
#print(outflow_list)

bg_b0_point3 = .3
bg_b0_2 = 2
bg_b0_1 = 1


sig = 0.3
va_all = np.sqrt(sig/(1+sig))

Va_point3_bg0 = np.sqrt(sig/(1+sig))
Va_point3_bgpoint3 = np.sqrt(sig / (1+sig+sig*bg_b0_point3**2))
Va_point3_bg2 = np.sqrt(sig / (1+sig+sig*bg_b0_2**2))
Va_point3_bg1 = np.sqrt(sig / (1+sig+sig*bg_b0_1**2))

Va_point3_bg0 = np.sqrt(sig/(1+sig))
Va_point3_bgpoint3 = np.sqrt(sig / (1+sig))
Va_point3_bg2 = np.sqrt(sig / (1+sig))
Va_point3_bg1 = np.sqrt(sig / (1+sig))





plt.title('$\sigma=0.3 \; \\beta=0.3$')
plt.plot(t_list, np.array(outflow_list[0])/va_all,label='$B_{g}/B_{0}=0$')
plt.plot(t_list, np.array(outflow_list[1])/va_all,label='$B_{g}/B_{0}=0.3$')
plt.plot(t_list,np.array(outflow_list[2])/va_all ,label='$B_{g}/B_{0}=0.5$')
plt.plot(t_list,np.array(outflow_list[3])/va_all,label='$B_{g}/B_{0}=0.7$')
plt.plot(t_list, np.array(outflow_list[4])/va_all,label='$B_{g}/B_{0}=1$')
plt.plot(t_list, np.array(outflow_list[5])/va_all,label='$B_{g}/B_{0}=3$')
plt.legend()
plt.xlim(1,45)
#plt.ylabel('Outflow Rate /$\sqrt{\sigma / (1+\sigma + \sigma_{g})}$')
plt.ylabel('Outflow rate ($v_{A}$)')
plt.xlabel('Output Timesteps')
plt.savefig('outflow_highbeta_6flds.png')


