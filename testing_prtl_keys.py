from tristan_funcs import get_val_filename, list_types
import h5py
import numpy as np
myprtl = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/deglam005/output/prtl.tot.030"
#myprtl = "../../tristan_acc-mec/test/output/prtl.tot.003"
myprtl = h5py.File(myprtl)

#list_types(myprtl)

tcs = np.array(get_val_filename('tcs',myprtl))
ycs = np.array(get_val_filename('ycs',myprtl))

boolarray = ycs == 0

print(tcs,ycs)
print(False in boolarray)
myprtl.close()

