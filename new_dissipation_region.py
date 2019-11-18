import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
from matplotlib.colors import PowerNorm

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

plt.set_cmap("viridis")



mybase = "../../tristan_acc-mec_Ez/testing_edotv/bguide.3_triggered_alwaystrack/output/"

flds_base0 = mybase + "flds.tot."
param = mybase + "param."


tstart = 10
tfinal = 11





for t in range(tstart,tfinal):
    t0_str = "%03d" % t
    f0 = h5py.File(flds_base0+t0_str,'r')
    fparam = h5py.File(param + t0_str,'r')
    
    bdens = f0['bdens'][0,:,:]
    bdensi = f0['bdensi'][0,:,:]
    dens = f0['dens'][0,:,:]#-bdens
    densi = f0['densi'][0,:,:]#-bdensi
    dense = dens-densi
    qi = np.array(fparam['qi'])[0]
    mi = np.array(fparam['mi'])[0]
    me = np.array(fparam['me'])[0]
    mass_ratio = mi/me
    #print('qi, mi, me')
    #print(qi, mi, me)

    keni = f0['keni'][0,:,:]
    ken = f0['ken'][0,:,:]
    kene = mass_ratio*(ken - keni) #think this should be avg kinetic energy
    #add one to make it actual lorentz factor instead of gamma -1
    kene +=1 
    
    print('max kene : ', np.max(kene))
    print('min kene : ', np.min(kene))
    #ok so I think what's actually being saved is q/m, so qe 

    

    rhoc = qi*(densi - dense)

    #print(qi)
    
    v3x = f0['v3x'][0,:,:]
    v3xi = f0['v3xi'][0,:,:]
    v3xe = (v3x * (me*dense + mi*densi) - v3xi*mi*densi)/(me*dense)

    v3y = f0['v3y'][0,:,:]
    v3yi = f0['v3yi'][0,:,:]
    v3ye = (v3y * (me*dense + mi*densi) - v3yi*mi*densi)/(me*dense)


    v3z = f0['v3z'][0,:,:]
    v3zi = f0['v3zi'][0,:,:]
    v3ze = (v3z * (me*dense + mi*densi) - v3zi*mi*densi)/(me*dense)
    '''
    #for plotting components of electron velocity
    fig, (ax0, ax1, ax2) = plt.subplots(1,3,sharey=True)
    im0 = ax0.imshow(v3xe,origin='lower',vmin=-1,vmax=1)
    im1 = ax1.imshow(v3ye,origin='lower',vmin=-1,vmax=1)
    im2 = ax2.imshow(v3ze,origin='lower',vmin=-1,vmax=1)

    cbar_ax2 = fig.add_axes([.92, .34, .02, .31])
    cb2 = fig.colorbar(im2,cax=cbar_ax2, orientation="vertical")
    ax0.set_title('v3xe')
    ax1.set_title('v3ye')
    ax2.set_title('v3ze')

    plt.savefig('v3xe_test.png',bbox_inches='tight')
    '''
    bx = f0['bx'][0,:,:]
    by = f0['by'][0,:,:]
    bz = f0['bz'][0,:,:]

    ex = f0['ex'][0,:,:]
    ey = f0['ey'][0,:,:]
    ez = f0['ez'][0,:,:]

    jx = -f0['jx'][0,:,:]
    jy = -f0['jy'][0,:,:]
    jz = -f0['jz'][0,:,:]

    

    #print(list(f0.keys()))
    #print('param keys')
    #print(list(fparam.keys()))
    #qi = np.array(fparam['qi'])[0]

     

    #construct arrays of vectors
    j_vecfield = np.dstack((jx,jy,jz))
    e_vecfield = np.dstack((ex,ey,ez))
    b_vecfield = np.dstack((bx,by,bz))
    ve_vecfield =np.dstack((v3xe,v3ye,v3ze))
    
    vcrossb = np.cross(ve_vecfield,b_vecfield)
    
    xlen,ylen,veclen=np.shape(j_vecfield)
    D=np.zeros((xlen,ylen))
    
    for i in range(xlen):
        for j in range(ylen):
            jvec = j_vecfield[i,j]
            vevec = ve_vecfield[i,j]
            bvec = b_vecfield[i,j]
            evec = e_vecfield[i,j]
            myrhoc = rhoc[i,j]
            tmp1 = evec + np.cross(vevec,bvec)
            tmp2 = np.dot(tmp1,jvec)
            tmp3 = myrhoc * np.dot(evec,vevec)
            tmp4 = tmp2 - tmp3

            D[i,j]=tmp4
    D *= kene
    #D = np.dot(j_vecfield,e_vecfield)

    #testing vector field shapes
    #print(np.shape(b_vecfield))
    #print(np.shape(bx))
    #print(np.shape(vcrossb))

    #plt.imshow(D,origin='lower',vmin=-.005,vmax=.005)
    xhlf = ylen / 2
    xscan = 50
    xlow = int(xhlf-xscan)
    xup = int(xhlf+xscan)
    print(xlen, ylen)
    print(xlow, xup)
    plt.imshow(D[200:-200,xlow:xup],origin='lower',vmin=-.1,vmax=.1)
    plt.colorbar(label='$D_{e}$')
    plt.savefig('diss_test.png')
    
    f0.close()
    fparam.close()
