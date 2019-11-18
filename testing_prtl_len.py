import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tristan_funcs import get_val_filename, list_types
import h5py
import numpy as np
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 10})

#myprtl = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/deglam005/output/prtl.tot.030"

myprtl = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam005/output/prtl.tot.030"

#myprtl = "../../tristan_acc-mec/correct_gam/output/prtl.tot.001"
#myprtl = "../../tristan_acc-mec/comm_test/output/prtl.tot.001"
myprtl = h5py.File(myprtl,'r')

#myfld = "../../tristan_acc-mec/test/output/flds.tot.009"
#myfld = h5py.File(myfld,'r')

#list_types(myprtl)

inde = np.array(get_val_filename('inde',myprtl))
indi = np.array(get_val_filename('indi',myprtl))
gammae = np.array(get_val_filename('gammae',myprtl))
gammai = np.array(get_val_filename('gammai',myprtl))

xe = np.array(get_val_filename('xe',myprtl))
ue = np.array(get_val_filename('ue',myprtl))

exe = np.array(get_val_filename('exe',myprtl))
we = np.array(get_val_filename('we',myprtl))

#tcsi = np.array(get_val_filename('tcsi',myprtl))
#ycsi = np.array(get_val_filename('ycsi',myprtl))

bsynce = np.array(get_val_filename('bsynce',myprtl))
bsynci = np.array(get_val_filename('bsynci',myprtl))


print(np.size(inde), np.size(gammae), np.size(xe), np.size(ue), np.size(exe),np.size(we),np.size(bsynce))

print(np.size(bsynci))
