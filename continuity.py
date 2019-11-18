import matplotlib
matplotlib.use('Agg')                                                         
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.special import kn
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

from inflow_test import mean_inflow
from tristan_funcs import outflow
from thickness_calc import return_thickness_slice
from overdens_calc import return_avg_dens


bg0 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide0/output/flds.tot.040"

bg3 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.3/output/flds.tot.040"
bg7 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.7/output/flds.tot.040"
bg1 = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide1/output/flds.tot.040"


def cont_check(file_string,ppc0):
    myf = h5py.File(file_string,'r')
    dens = myf['dens'][0,:,:]
    mylen = np.shape(dens)[0]/2
    #print(mylen)
    thicklist = return_thickness_slice(file_string)
    meandens, stddens = return_avg_dens(file_string)
    myoutflow = outflow(myf)
    myinflow = mean_inflow(myf)
    myf.close
    outlist = [np.mean(thicklist),myoutflow,meandens]
    inlist = [mylen, myinflow, ppc0]
    cont_out = myoutflow*meandens*np.mean(thicklist)
    cont_in = myinflow*ppc0*mylen
    print('inflow: ',cont_in)
    print('outflow: ',cont_out)
    return cont_out, cont_in, outlist, inlist

bg0_out, bg0_in, bg0_outlist, bg0_inlist = cont_check(bg0,4)
bg3_out, bg3_in, bg3_outlist, bg3_inlist = cont_check(bg3,4)
bg7_out, bg7_in, bg7_outlist, bg7_inlist = cont_check(bg7,4)
bg1_out, bg1_in, bg1_outlist, bg1_inlist = cont_check(bg1,4)

bglist = [0, .3, .7, 1]

norm0_0 = (bg0_outlist[0]/bg0_inlist[0])
norm0_1 = bg0_outlist[1]/bg0_inlist[1]
norm0_2 = (bg0_outlist[2]/bg0_inlist[2])
norm0_0 = 1
norm0_1 = 1
norm0_2 = 1

plt.scatter(bglist[0],bg0_outlist[0]/bg0_inlist[0]/norm0_0,color="Blue",label='$w/L$')
plt.scatter(bglist[0],bg0_outlist[1]/bg0_inlist[1]/norm0_1,color="Red",label='$v_{out}/v_{in}$')
plt.scatter(bglist[0],bg0_outlist[2]/bg0_inlist[2]/norm0_2,color="Orange",label="$n_{sheet}/n_{0}$")

plt.scatter(bglist[1],bg3_outlist[0]/bg3_inlist[0]/norm0_0,color="Blue")
plt.scatter(bglist[1],bg3_outlist[1]/bg3_inlist[1]/norm0_1,color="Red")
plt.scatter(bglist[1],bg3_outlist[2]/bg3_inlist[2]/norm0_2,color="Orange")

plt.scatter(bglist[2],bg7_outlist[0]/bg7_inlist[0]/norm0_0,color="Blue")
plt.scatter(bglist[2],bg7_outlist[1]/bg7_inlist[1]/norm0_1,color="Red")
plt.scatter(bglist[2],bg7_outlist[2]/bg7_inlist[2]/norm0_2,color="Orange")

plt.scatter(bglist[3],bg1_outlist[0]/bg1_inlist[0]/norm0_0,color="Blue")
plt.scatter(bglist[3],bg1_outlist[1]/bg1_inlist[1]/norm0_1,color="Red")
plt.scatter(bglist[3],bg1_outlist[2]/bg1_inlist[2]/norm0_2,color="Orange")
#plt.yscale('log')
plt.ylim(0,2)
plt.legend(frameon=False,loc='upper left')
plt.xlabel('$B_{g}/B_{0}$')
plt.ylabel('Outflow Inflow Ratio')
plt.savefig('cont_ratios.png',dpi=300,bbox_inches='tight')



'''
plt.scatter(bglist[0],bg0_out,marker='o',color='C1',label='$n_{sheet}v_{out}w$')
plt.scatter(bglist[0],bg0_in,color='C1',marker='x',label='$n_{0}v_{in}L$')
plt.scatter(bglist[1],bg3_out,color='C1',marker='o')
plt.scatter(bglist[1],bg3_in,color='C1',marker='x')
plt.scatter(bglist[2],bg7_out,color='C1',marker='o')
plt.scatter(bglist[2],bg7_in,color='C1',marker='x')
plt.scatter(bglist[3],bg1_out,color='C1',marker='o')
plt.scatter(bglist[3],bg1_in,color='C1',marker='x')
plt.legend(frameon=False)
plt.ylabel('Mass flux')
plt.xlabel('$B_{g}/B_{0}$')
plt.savefig('cont_test.png',dpi=300,bbox_inches='tight')
plt.close()
'''





