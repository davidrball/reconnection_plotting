
import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

#inputs

#prtl_base = "../../tristan_acc-mec_Ez/testing_edotv/bguide.1_untriggered/output/prtl.tot."
prtl_base = "../../tristan_acc-mec_Ez/testing_edotv/untriggered/output/prtl.tot."

#20 to 40 seem messed up, not sure why...
tlow=  1
tup = 50
tstep = 1




mymax = 0

for t in range(tlow,tup,tstep):
    fig, ax1 = plt.subplots(1)
    tstr = '%03d'%t
    filename = prtl_base + tstr
    print(filename)
    f_prtl = h5py.File(filename,'r')
    gammae = np.array(f_prtl['gammae'])
    edotv = np.array(f_prtl['edotve'])
    inde = np.array(f_prtl['inde'])
    proce = np.array(f_prtl['proce'])

    numprtls = np.size(gammae)


    print('maxgam : ', np.max(gammae))

    tmpmax = np.max(gammae)
    mymax = max([mymax,tmpmax])

    edotv *= 2 #because the ex0 saved is actually ex/2 because mover does 2 episodes of E field acceleration

    edotv /= .45 #need factor of c?  not sure about this one

    tcol = (t-tlow)/(tup-tlow)
    tcol_rgb = (tcol/4,tcol/2,tcol)

    #print('min max edotve')
    #print(np.min(edotve))
    #print(np.max(edotve))
    ax1.scatter(edotv, gammae,s=1)#,label="$t="+tstr[1:]+"$")
    #ax1.scatter(inde,gammae,s=.5,label=tstr)
    f_prtl.close()

    xarr = np.linspace(1e-1,1e4,50)
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.plot(xarr, xarr, linestyle='dashed',color="Red",label="$\gamma=W_{||}$")
    ax1.legend(loc="upper left",frameon=False,ncol=1,prop={'size':12})
    
#ax1.set_ylim(0,1e3)
    
#ax1.set_xlim(0,1e6)
#ax1.set_ylim(0,mymax)
    
    ax1.set_xlim(1e-1,1e4)
    ax1.set_ylim(20,1e4)
    
    ax1.set_xlabel('$W_{||}$')
    ax1.set_ylabel('$\gamma$')
    plt.savefig('wpar_vid/bguide.3_untriggered/'+tstr[1:]+'.png',dpi=300,bbox_inches='tight')
    plt.close()
