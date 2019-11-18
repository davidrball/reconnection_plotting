import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tristan_funcs import get_val_filename, return_rspec, return_spec
import h5py
import numpy as np

plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams.update({'font.size':15})

prtl_final = "../../tristan_acc-mec_Ez/8k_run2/output/prtl.tot.035"
f_fld = "../../tristan_acc-mec_Ez/8k_run2/output/flds.tot.035"
spec_final = "../../tristan_acc-mec_Ez/8k_run2/output/spect.020"

f_prtl = h5py.File(prtl_final,'r')
f_fld = h5py.File(f_fld,'r')
f_spect = h5py.File(spec_final,'r')

specgam, specion, speclec = return_rspec([f_spect])
dens = np.rot90(get_val_filename('dens',f_fld)[0,:,:])
c_omp = 3
ycs_conv = c_omp * 1000 #just to convert ycse from cells to the units we show in the fld plots for easier comparison

istep = 12
xext = np.shape(dens)[1]*istep / ycs_conv



ycse = np.array(get_val_filename('ycse',f_prtl))
gammae = np.array(get_val_filename('gammae',f_prtl))
tcse = np.array(get_val_filename('tcse',f_prtl))
prtnum = np.size(gammae)

#make a list of gammas from the central region and from elsewhere

xpoint_gams = []
else_gams = []

xlow = -.5
xup = .5
 
for i in range(prtnum):
    if tcse[i] != 0:
        myycse = ycse[i] #in units of cells
        myycse /= ycs_conv
        myycse -= xext/2.
        #now myycse is in the same units as we show in the other plots w/ flds
        if xlow < myycse < xup:
            xpoint_gams.append(gammae[i])
        else:
            else_gams.append(gammae[i])


xpoint_gams = np.array(xpoint_gams)
else_gams = np.array(else_gams)
gamup = max(np.max(xpoint_gams),np.max(else_gams))
gamlow = min(np.min(xpoint_gams),np.min(else_gams))
bin_num = 50

histbins, xpoint_hist = return_spec(xpoint_gams,gamlow,gamup,bin_num)
histbins2, else_hist = return_spec(else_gams,gamlow,gamup,bin_num)

stride = 2
plt.plot(histbins, histbins*xpoint_hist,color="Blue",label='$|y_{cs}|<0.5$')
plt.plot(histbins2, histbins2*else_hist,color="Red",label='$|y_{cs}|>0.5$')
plt.plot(specgam[0],specgam[0]*speclec[0]/stride,color="Black",linestyle='dashed',label='total spectrum')
plt.legend(frameon=False)
plt.xscale('log')
plt.yscale('log')
plt.ylim(100,2e7)
plt.xlim(gamlow,2*gamup)
plt.xlabel('$\gamma$')
plt.ylabel('$\gamma dN/d\gamma$')
plt.title('Triggered Bguide=0')
#plt.ylim(10,3e3)
plt.savefig('triggered_bguide0_prt_spec.png')


f_prtl.close()
f_fld.close()
f_spect.close()
