import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from tristan_funcs_vecpot import vecpot2
import h5py
from scipy.ndimage.filters import gaussian_filter
plt.set_cmap('RdBu')



plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 15})

t=47
t_str = '%03d' % t
#fldbase = "../../tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam00005/output/flds.tot."

#fldbase = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/flds.tot."

fldbase = "../../tristan_acc-mec_Ez/8k_untriggered_bguide.3/output/flds.tot."

#fldbase = "../../tristan_acc-mec_Ez/8k_bguide0_untriggered_stride1_thresh2/output/flds.tot."

fldbase += t_str
myfld = h5py.File(fldbase,'r')
dens = myfld['dens'][0,:,:]
bdens = myfld['bdens'][0,:,:]
edge = 5
smoothlen = 1

vecpot = vecpot2(myfld,12,3)[edge:-1*edge,:] #clip off top and bottom to remove numerical issue at boundary with vector potential calc
vecpot_smooth = gaussian_filter(vecpot,smoothlen)




#just getting bdens to be the same shape as vecpot
mx = bdens.shape[1]
bdens = bdens[:,2:mx-5]
bdens = bdens[edge:-1*edge,:]



print(np.shape(vecpot),np.shape(bdens))

xhlf = np.shape(vecpot)[1]/2

#print(xhlf)

#vecpot_slice = vecpot[:,xhlf]
#vecpot_smooth_slice = vecpot_smooth[:,xhlf]

#plt.plot(vecpot_slice,color="Red")
#plt.plot(vecpot_smooth_slice,color="Blue")
#plt.savefig('vecpot_slices_smooth.png')
#looks decent now let's do the differencing

#0 axis is up along y
#1 axis is in x-direction

#basic central difference:
myroll = 2




dady = (np.roll(vecpot_smooth,myroll,axis=0) -np.roll(vecpot_smooth,-1*myroll,axis=0))/myroll
dadx = (np.roll(vecpot_smooth,myroll,axis=1) - np.roll(vecpot_smooth,-1*myroll,axis=1))/myroll
d2ad2x = (np.roll(vecpot_smooth,myroll,axis=1)-2*vecpot_smooth+np.roll(vecpot_smooth,-1*myroll,axis=1))/myroll**2
d2ad2y = (np.roll(vecpot_smooth,myroll,axis=0)-2*vecpot_smooth+np.roll(vecpot_smooth,-1*myroll,axis=0))/myroll**2
d2adydx = (np.roll(dady,myroll,axis=1) - np.roll(dady,-1*myroll,axis=1))/myroll
d2adxdy = (np.roll(dadx,myroll,axis=0) - np.roll(dadx,-1*myroll,axis=0))/myroll





#first_deriv_tot = dady+dadx


#first_deriv_tot = first_deriv_tot[:,smoothlen+1:-smoothlen-2] #cut off stuff from rolling at end to avoid analyzing nonsense
dady = dady[:,smoothlen+1:-smoothlen-2]
dadx = dadx[:,smoothlen+1:-smoothlen-2]
d2ad2x = d2ad2x[:,smoothlen+1:-smoothlen-2]
d2ad2y = d2ad2y[:,smoothlen+1:-smoothlen-2]
d2adydx = d2adydx[:,smoothlen+1:-smoothlen-2]
d2adxdy = d2adxdy[:,smoothlen+1:-smoothlen-2]
#so first we're gonna have to loop through the first derivative and identify zeros.  

bdens = bdens[:,smoothlen+1:-smoothlen-2]


convfac = 12/3./1000
xext = np.shape(dens)[1]*convfac
yext = np.shape(dens)[0]*convfac


bdens_tolerence = 4
zero_tolerence = 2


bdens_loc =  bdens<bdens_tolerence
dadx_zero_loc = np.abs(dadx)<zero_tolerence
dady_zero_loc = np.abs(dady)<zero_tolerence

first_deriv_zero_loc = dadx_zero_loc*dady_zero_loc*bdens_loc


#first_deriv_zero_loc = np.abs(first_deriv_tot)<zero_tolerence #locations where the first derivative is sufficiently close to zero (sufficiently being controlled by parameter "zero tolerence")

num_deriv_zero = np.sum(first_deriv_zero_loc)
print('number of zero derivative points : ',num_deriv_zero)

#now we loop through null points and ask the questions about 2nd derivatives
ylen,xlen = np.shape(dady)
#print(xlen,ylen)

#testarr = np.zeros((ylen,xlen))
#print(np.shape(testarr))

saddlepoint_arr = np.zeros((ylen,xlen))
saddle_count = 0
xpoint_loc_list = []

for i in range(ylen):
    for j in range(xlen):
        

        #only need to do operations if its a null point
        if first_deriv_zero_loc[i][j]==1:
            #testarr[i][j] += 1 #indexing is correct
            
            H00 = d2ad2x[i][j] #values of the hessian matrix
            H11 = d2ad2y[i][j]
            H01 = d2adydx[i][j]

            #coeffs of quadratic equations solving for eigenvals
            a=1
            b=H00+H11
            c=H00*H11 - H01**2

            #eigenvals
            lambdaplus = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
            lambdaminus = (-b - np.sqrt(b**2 -4*a*c))/(2*a)

            #print('eigenvals : ',lambdaplus,lambdaminus)
            #print('eigenvec : ',-b,a-lambdaplus)


            if lambdaplus*lambdaminus < 0: #if they have opposite signs
            #if True:
                xpoint_loc_list.append([i,j])
                saddle_count += 1
                saddlepoint_arr[i][j]+=1
print('saddle point count : ',saddle_count)


plt.set_cmap('viridis')
#testing various quantities along the way
fig, ax0 = plt.subplots(1,1)#,sharey=True)
ax0.imshow(dens/4.,origin='lower',vmin=0,vmax=5,extent=[-xext/2, xext/2, -yext/2,yext/2])
for i in range(saddle_count):
    ax0.scatter((xpoint_loc_list[i][1]+4)*convfac - xext/2,(xpoint_loc_list[i][0]+1)*convfac - yext/2,marker='x',s=10,color="Red")
ax0.set_xlabel('$x \; (1000 c/\omega_{p})$')
ax0.set_ylabel('$y \; (1000 c/\omega_{p})$')

#ax0.set_title('$\\frac{\partial A_{z}}{\partial x} + \\frac{\partial A_{z}}{\partial y}$')
#ax1.imshow(saddlepoint_arr,origin='lower')
#ax1.set_title('Saddle Points')

#plt.imshow(first_deriv_tot,origin='lower',vmin=-150,vmax=150)#,vmin=-5,vmax=5)
#plt.colorbar()
plt.savefig('new_xpoint_algorithm_untrig_bguide.png',dpi=300,bbox_inches='tight')

#plt.imshow(vecpot_smooth,origin='lower')
#plt.savefig('vecpot_smoothed.png',dpi=300,bbox_inches='tight')
#plt.close()
#plt.imshow(vecpot,origin='lower')
#plt.savefig('vecpot_orig.png',dpi=300,bbox_inches='tight')
#plt.close()


myfld.close()
