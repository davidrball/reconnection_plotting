import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})


t_offset = 5
offset_str = str(t_offset)
#filename = "bguide.3_triggered_injhist_offset"+offset_str + "_gamthresh300_tfinal60.txt"
filename = "bguide.1_triggered_injhist_offset"+offset_str  +"_gamthresh300_tfinal60.txt"

path = "wpar_gam_txt/"
filename = path+filename

myfile = open(filename,'r')
outlist = []
linesum = 0
for line in myfile.readlines():
    
    linesum+=1
    mylist = line.split(' , ')
    floatlist = []
    print('len list : ',len(mylist))

    for i in range(len(mylist)-1):
        floatlist.append(float(mylist[i]))

    outlist.append(floatlist)


wpar_list = outlist[0]
gam_list = outlist[1]

wparlog = np.log10(np.abs(wpar_list))
gamlog = np.log10(gam_list)

wparmin = np.min(wparlog)
wparmin=2.3

wparmax = np.max(wparlog)
wparmax = 3.6

gammin = np.min(gamlog)
gammax = np.max(gamlog)
gammin=2.3
gammax=3.6


wpar_bins = np.linspace(wparmin,wparmax,40)
gambins=np.linspace(gammin,gammax,40)

plt.hist2d(wparlog, gamlog, bins=(wpar_bins,gambins),norm=LogNorm())
plt.xlabel('$\log(W_{|| \; , \; z})$')
plt.ylabel('$\log(\gamma)$')
xarr= np.linspace(-3,4,10)
plt.plot(xarr,xarr,linestyle='dashed',color="Orange")
plt.colorbar(label="$N_{e}$")


offset_str_comp = str(t_offset*300 + 150)

plt.annotate('$t=t_{inj}+'+offset_str_comp+' \omega_{p}^{-1}$',(wparmin+.1,gammax-.1),color="White")



savename = "bguide.1_triggered_wpar_inj_offset"+offset_str+"_log"
plt.savefig(savename + '.png',dpi=300,bbox_inches='tight')
plt.savefig(savename + '.pdf',dpi=300,bbox_inches='tight')
plt.close()
