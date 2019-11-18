import matplotlib
matplotlib.use('Agg') #for elgato
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.special import kn
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 10})

from tristan_funcs import get_val_filename

#spec_base = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide2/output/momentum.0"
#spec_base = "../../tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide1/output/momentum.0"
spec_base = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam005/output/momentum.0"

t0 = 40
tf = 41

px_std_list = []
py_std_list = []
pz_std_list = []
t_list = []


for t in range(t0, tf):
    t_str = str(t)
    if len(t_str)==1:
        spec0 = spec_base + "0"+t_str
    elif len(t_str)==2:
        spec0 = spec_base + t_str

    t_list.append(t)
    #spec0 = spec_base + str(t)
    myf = h5py.File(spec0,'r')
    px = np.array(get_val_filename('pxelogsp',myf))#[:,0:10]
    py = get_val_filename('pyelogsp',myf)#[:,0:10]
    pz = get_val_filename('pzelogsp',myf)#[:,0:10]

    pxb = np.array(get_val_filename('pxeblogsp',myf))
    print(np.shape(pxb))
    px -= pxb

    pxbin = np.array(get_val_filename('pxbin',myf))
    pybin = np.array(get_val_filename('pybin',myf))
    pzbin = np.array(get_val_filename('pzbin',myf))
    #testbin = get_val_filename('nonsense',myf)
    #print(pxbin)
    px_spec = np.zeros(np.size(pxbin))
    py_spec = np.zeros(np.size(pybin))
    pz_spec = np.zeros(np.size(pzbin))

    xsize = np.shape(px)[1]

    ysize = np.shape(px)[0]
    for i in range(xsize):
        plt.plot(pxbin, np.abs(pxbin)*px[:,i])
plt.yscale('log')
plt.ylim(1e2,1e6)
plt.xlim(-10,10)
plt.savefig('mom_all_slice.png')
