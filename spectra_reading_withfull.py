import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 16})

mybase = "../../tristan_acc-mec_Ez/bguide.1_untriggered_stride1_gamthresh50/output/spect.065"
#mybase = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/spect.060"
#mybase = "../../tristan_acc-mec_Ez/bguide.3_untriggered_stride1_gamthresh50/output/spect.060"
#mybase = "../../tristan_acc-mec_Ez/bguide.1_triggered_stride1_gamthresh50/output/spect.065"

myspec = h5py.File(mybase,'r')

rspec = np.array(myspec['rspece'])
rspecshape = np.shape(rspec)

mygam = np.array(myspec['gamma'])
dgam = np.array(myspec['dgam'])
print(dgam)
gambins = np.size(mygam)
myrspec = np.zeros(gambins)

for j in range(gambins):
    myrspec[j]+= np.sum(rspec[j][:])
myspec.close()
filepath = "classification_txt/bguide.1_untriggered_spec.txt"

target = open(filepath,'r')
target.readline() #reading histbins
histbins = target.readline().split(' , ')[:-2]
target.readline() #reading xpoint spec
xpoint_spec = target.readline().split(' , ')[:-2]
target.readline() #reading merger spec
merger_spec = target.readline().split(' , ')[:-2]
target.readline() #reading other spec
other_spec = target.readline().split(' , ')[:-2]

#convert lists of strings to lists of floats
histbins_f = []
xpoint_spec_f = []
merger_spec_f = []
other_spec_f = []
mylen = len(histbins)

for i in range(mylen):
    histbins_f.append(float(histbins[i]))
    xpoint_spec_f.append(float(xpoint_spec[i]))
    merger_spec_f.append(float(merger_spec[i]))
    other_spec_f.append(float(other_spec[i]))

xpoint_spec_f = np.array(xpoint_spec_f)
merger_spec_f = np.array(merger_spec_f)
other_spec_f = np.array(other_spec_f)
histbins = np.array(histbins_f)

xarr = 275*np.ones(10)
yarr = np.linspace(1e1,1e7,10)

plt.plot(histbins,xpoint_spec_f,color="Blue",label='X-point (PCS)')
plt.plot(histbins,merger_spec_f,color="Orange",label='X-point (MCS)')
plt.plot(histbins,other_spec_f,color="Red",label='Other')
plt.plot(xarr, yarr, color="Cyan",linestyle="dashed",label='$\sigma_{e}/2$')

plt.plot(mygam, mygam*myrspec,color="Black",linestyle='--',label='Total Spectrum')



plt.legend(frameon=False,prop={'size':12})
plt.xscale('log')
plt.yscale('log')
plt.ylim(10,1e7)
plt.xlim(10,1e4)
plt.xlabel('$\gamma$')
plt.ylabel('$\gamma dN/d\gamma$')
plt.savefig('bguide.1_untriggered_withspec.pdf',dpi=300,bbox_inches='tight')
plt.savefig('bguide.1_untriggered_withspec.png',dpi=300,bbox_inches='tight')
