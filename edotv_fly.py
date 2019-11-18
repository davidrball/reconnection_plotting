
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

prtl_base = "../../tristan_acc-mec_Ez/testing_edotv/output/prtl.tot."
#prtl_base = "../../tristan_acc-mec_Ez/testing_edotv/recon_rgn/output/prtl.tot."

#20 to 40 seem messed up, not sure why...
tlow=  10
tup = 31
tstep = 20
maxind = 1.2e9 #taken empirically from first timestep, where particles that are initialized in sheet are actually saved in the output, 1 per proc (should be ~340 of them), let's just delete them out here
minind = 1e5


fig, ax1 = plt.subplots(1)

mymax = 0

for t in range(tlow,tup,tstep):
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

    #coldprtls = gammae < 20.
    #gammae*coldprtls #only look at particles accelerated out of thermal pool


    #goodprtls = inde<maxind

    #print(badprtls)
    #badprtlsmin = inde < minind
    #numbad = np.sum(badprtls) + np.sum(badprtlsmin)
    #print('num bad : ',numbad)
    #gammae = gammae*goodprtls
    #gammae = gammae*badprtlsmin
    tmpmax = np.max(gammae)
    mymax = max([mymax,tmpmax])

    edotv *= 2 #because the ex0 saved is actually ex/2 because mover does 2 episodes of E field acceleration

    edotv /= .45 #need factor of c?  not sure about this one

    tcol = (t-tlow)/(tup-tlow)
    tcol_rgb = (tcol/4,tcol/2,tcol)

    #print('min max edotve')
    #print(np.min(edotve))
    #print(np.max(edotve))
    ax1.scatter(edotv, gammae,s=1,label="$t="+tstr[1:]+"$")
    #ax1.scatter(inde,gammae,s=.5,label=tstr)
    f_prtl.close()

xarr = np.linspace(0,mymax,50)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.plot(xarr, xarr, linestyle='dashed',color="Red",label="$\gamma=W_{||}$")
ax1.legend(loc="upper left",frameon=False,ncol=1,prop={'size':12})

#ax1.set_ylim(0,1e3)

#ax1.set_xlim(0,1e6)
#ax1.set_ylim(0,mymax)

ax1.set_xlim(20,1e4)
ax1.set_ylim(200,1e4)

ax1.set_xlabel('$W_{||}$')
ax1.set_ylabel('$\gamma$')
plt.savefig('edotv_fly_recon_newtest.png')

