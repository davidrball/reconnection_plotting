import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt

from tristan_funcs import get_val_filename, decompose_spectra, return_spece_withcut, return_rspec


#plt.rcParams['mathtext.fontset'] = 'stix'
#plt.rcParams['font.family'] = 'STIXGeneral'
#plt.rcParams.update({'font.size': 16})


#give it the ax you want to plot, filebase, init step, final step, stride
def plot_timespec(ax, filebase, t0, tf, stride):
    for t in range(t0, tf+1, stride):
        myf = filebase
        t_string = str(t)
        if t<10:
            myf+="00"
        elif t>=10 and t<100:
            myf+="0"
        elif t>=100:
            pass
        myf += t_string
        print(myf)
        f = h5py.File(myf,'r')        
        gam_spece, ion_spece, lec_spece = return_spece_withcut([f], [16000])
        gam_spece = gam_spece[0]
        ion_spece = ion_spece[0]
        lec_spece = lec_spece[0]
        gam_rspec, ion_rspec, lec_rspec = return_rspec([f])
        gam_rspec = gam_rspec[0]
        ion_rspec = ion_rspec[0]
        lec_rspec = lec_rspec[0]
        f.close()
        #gam_array = np.array(get_val_filename("gamma",f))
        gam_array = gam_spece
        decomp_gam, decomp_rspec = decompose_spectra(lec_spece,lec_rspec,gam_array)
        tcol = t/tf
        if tcol < 1:
            tcol_rgb = (1-tcol,1-tcol,tcol)
            ax.plot(decomp_gam, decomp_gam*decomp_rspec, color = tcol_rgb, linewidth = 1.1)
            ax.plot(gam_array,gam_array*lec_spece, color=tcol_rgb, linewidth=.5)
        elif tcol==1:
            ax.plot(decomp_gam, decomp_gam*decomp_rspec, color = "Red", linewidth = 1.5)
            ax.plot(gam_array,gam_array*lec_spece, color="Red", linewidth=.5)

def plot_timespec_label(ax, filebase, t0, tf,stride,interval,comp):
    for t in range(t0, tf+1, stride):
        myf = filebase
        t_string = str(t)
        if t<10:
            myf+="00"
        elif t>=10 and t<100:
            myf+="0"
        elif t>=100:
            pass
        myf += t_string
        print(myf)

        t_plasmafreq = t*interval*.45/comp
        tlabel_str = str(t_plasmafreq)[:-2]
        f = h5py.File(myf,'r')
        gam_spece, ion_spece, lec_spece = return_spece_withcut([f], [16000])
        gam_spece = gam_spece[0]
        ion_spece = ion_spece[0]
        lec_spece = lec_spece[0]
        gam_rspec, ion_rspec, lec_rspec = return_rspec([f])
        gam_rspec = gam_rspec[0]
        ion_rspec = ion_rspec[0]
        lec_rspec = lec_rspec[0]
        f.close()
        #gam_array = np.array(get_val_filename("gamma",f))                                                                                                                               
        gam_array = gam_spece
        decomp_gam, decomp_rspec = decompose_spectra(lec_spece,lec_rspec,gam_array)
        tcol = t/tf
        if tcol < 1:
            tcol_rgb = (1-tcol,1-tcol,tcol)
            ax.plot(decomp_gam, decomp_gam*decomp_rspec, color = tcol_rgb, linewidth = 1.1,label=tlabel_str)
            ax.plot(gam_array,gam_array*lec_spece, color=tcol_rgb, linewidth=.5)
        elif tcol==1:
            ax.plot(decomp_gam, decomp_gam*decomp_rspec, color = "Red", linewidth = 1.5,label=tlabel_str)
            ax.plot(gam_array,gam_array*lec_spece, color="Red", linewidth=.5)
        plt.legend(loc='lower left', frameon=False)


#testing
'''
mybase = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/spect."

tstart = 15
t_end = 50

fig, ax1 = plt.subplots(1,1)

plot_timespec(ax1, mybase,tstart, t_end,5)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(1e-1,1e4)
ax1.set_ylim(1e4,1e9)
plt.savefig('testing_timespec.png')
'''
