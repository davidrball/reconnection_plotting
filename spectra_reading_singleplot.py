import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 12})

bg1_ut = "../../tristan_acc-mec_Ez/bguide.1_untriggered_stride1_gamthresh50/output/spect.065"
bg3_t = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/spect.060"
bg3_ut = "../../tristan_acc-mec_Ez/bguide.3_untriggered_stride1_gamthresh50/output/spect.060"
bg1_t = "../../tristan_acc-mec_Ez/bguide.1_triggered_stride1_gamthresh50/output/spect.065"

baselist = [bg3_t, bg3_ut, bg1_t,bg1_ut]
filepathbase = "classification_txt/"
bg1_ut_filepath = filepathbase + "bguide.1_untriggered_spec.txt"
bg1_t_filepath = filepathbase + "bguide.1_triggered_spec.txt"
bg3_ut_filepath = filepathbase + "bguide.3_untriggered_spec.txt"
bg3_t_filepath = filepathbase + "bguide.3_triggered_spec.txt"

bg1_str = "$B_{g}=0.1B_{0}$"
bg3_str = "$B_{g}=0.3B_{0}$"


filepathlist = [bg3_t_filepath, bg3_ut_filepath, bg1_t_filepath, bg1_ut_filepath]

fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, sharex=True,sharey=True)
axlist = [ax0,ax1,ax2,ax3]
for k in range(len(baselist)):
    myax = axlist[k]
    mybase = baselist[k]
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
    filepath = filepathlist[k]

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

    myax.plot(histbins,xpoint_spec_f/dgam,color="Blue",label='X-point (PCS)')
    myax.plot(histbins,merger_spec_f/dgam,color="Orange",label='X-point (MCS)')
    myax.plot(histbins,other_spec_f/dgam,color="Red",label='Other')
    myax.plot(xarr, yarr, color="Cyan",linestyle="dashed",label='$\sigma_{e}/2$')
    
    myax.plot(mygam, mygam*myrspec/dgam,color="Black",linestyle='--',label='Total Spectrum')


    if k==0:
        myax.legend(frameon=False,prop={'size':6})

    if k==0 or k==1:
        myax.annotate(bg3_str,(15,200))
    elif k==2 or k==3:
        myax.annotate(bg1_str, (15,200))

    if k==0 or k==2:
        myax.annotate('Triggered',(15,750))
    elif k==1 or k==3:
        myax.annotate('Untriggered',(15,750))

    myax.set_xscale('log')
    myax.set_yscale('log')
    myax.set_ylim(100,1e8)
    myax.set_xlim(10,1e4)
    if k==2 or k==3:
        myax.set_xlabel('$\gamma-1$')
    if k==0 or k==2:
        myax.set_ylabel('$(\gamma-1) dN/d\gamma$')
    #myax.set_yticks([1e1,1e3,1e5])
plt.subplots_adjust(wspace=.1, hspace=.1)
plt.savefig('spectra_fourplot.pdf',dpi=300,bbox_inches='tight')
plt.savefig('spectra_fourplot.png',dpi=300,bbox_inches='tight')
