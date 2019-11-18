
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib.colors import LogNorm

from matplotlib.ticker import MultipleLocator

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

#sigpoint3_fp = open('sigmapoint3_inflow_rates.txt','r')
sigpoint3_fp = open('sigma.3_guide_inflow_highbeta_bg1.txt','r')
sigpoint3_inflow_list = []
for line in sigpoint3_fp.readlines():
    mylist = line.split(' , ')
    floatlist = []
    for i in range(len(mylist)):
        if i==0:
            floatlist.append(mylist[i])
        else:
            
            floatlist.append(float(mylist[i]))
    #filename = mylist[0]
    #inflow = mylist[1:]
    
    sigpoint3_inflow_list.append(floatlist)
#print(inflow_list[0][1:])
mylen = len(sigpoint3_inflow_list[0][1:])
print('mylen = :'+str(mylen))

sigmapoint3_bguide0_inflow = np.array(sigpoint3_inflow_list[0][1:])
sigmapoint3_bguidepoint3_inflow = np.array(sigpoint3_inflow_list[1][1:])
sigmapoint3_bguidepoint7_inflow = np.array(sigpoint3_inflow_list[2][1:])
sigmapoint3_bguide1_inflow = np.array(sigpoint3_inflow_list[3][1:])
sigmapoint3_bguide2_inflow = np.array(sigpoint3_inflow_list[4][1:])
Va_point3 = np.sqrt(.3/1.3)


bg_b0_point3 = .3
bg_b0_point7 = .7
bg_b0_1 = 1


sig = 0.3
va_all = np.sqrt(sig/(1+sig))
Va_point3_bgpoint3 = np.sqrt(sig / (1+sig+sig*bg_b0_point3**2))
Va_point3_bgpoint7 = np.sqrt(sig / (1+sig+sig*bg_b0_point7**2))
Va_point3_bg1 = np.sqrt(sig / (1+sig+sig*bg_b0_1**2))

print(Va_point3)
print(Va_point3_bgpoint3)
print(Va_point3_bgpoint7)
print(Va_point3_bg1)

#plt.plot(sigma3_delgam05_inflow/Va_3, label='$\\beta=0.025$')                     
#plt.plot(sigma3_delgam005_inflow/Va_3, label='$\\beta=0.003$')                    
#plt.plot(sigma3_delgam0005_inflow/Va_3, label='$\\beta=0.0003$')
size1 = np.size(sigmapoint3_bguide0_inflow)
size2 = np.size(sigmapoint3_bguidepoint3_inflow)
size3 = np.size(sigmapoint3_bguidepoint7_inflow)
size4 = np.size(sigmapoint3_bguide1_inflow)
size5 = np.size(sigmapoint3_bguide2_inflow)
xarr1 = (np.linspace(0,size1,size1)*300)/1000
xarr2 = (np.linspace(0,size2,size1)*300)/1000
xarr3 = (np.linspace(0,size3,size1)*300)/1000
xarr4 = (np.linspace(0,size4,size1)*300)/1000
xarr5 = (np.linspace(0,size5,size1)*300)/1000
print(size1, np.size(xarr1))

 
plt.plot(xarr1,sigmapoint3_bguide0_inflow/va_all, label='$B_{g}/B_{0}=0$')
plt.plot(xarr2,sigmapoint3_bguidepoint3_inflow/va_all, label='$B_{g}/B_{0}=0.1$')
plt.plot(xarr3,sigmapoint3_bguidepoint7_inflow/va_all, label='$B_{g}/B_{0}=0.3$')
plt.plot(xarr4, sigmapoint3_bguide1_inflow/va_all,label='$B_{g}/B_{0}=0.7$')
plt.plot(xarr5, sigmapoint3_bguide2_inflow/va_all,label='$B){g}/B_{0}=1$')
plt.annotate('$\sigma=0.3$',(1,.13),size=18)
plt.annotate('$\\beta=0.3$',(1,.12),size=18)
plt.xlabel('Time ($1000 \; \omega_{p}^{-1}$)',fontsize=16)
plt.ylabel('Inflow Rate ($v_{A}$)',fontsize=16)
#plt.xlim(0,50)
#plt.rcParams.update({'font.size': 10})
plt.legend(prop={'size':14},loc="upper left",frameon=False)
plt.ylim(0,.08)
plt.xlim(.5,13)
plt.savefig('plots/guide_inflows_highbeta.png', dpi=300,bbox_inches='tight')
plt.close()
