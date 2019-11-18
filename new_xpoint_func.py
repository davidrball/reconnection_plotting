import numpy as np
from tristan_funcs_vecpot import vecpot2
import h5py
from scipy.ndimage.filters import gaussian_filter

def return_xpoints(myfld,edge,smoothlen,myroll,bdens_tolerence,zero_tolerence): #value of 5 is fine for edge, 1 or 2 for myroll (distance over which derivatives are calculated

    #bdens_tolerence makes it so we don't care about stuff that happens where the particles are all just particles initialized in the current sheet, 4 works fine
    
    #zero_tolerence is how close the derivatives have to be to 0 in order to count, 2 is fine for this

    #one thing we could implement is a higher order construction of the derivative, might be worth looking into

    leftcut = 0 #adding up all the cuts to the arrays we make so we can return them at the end so we can make sure we line up the locations properly
    rightcut = 0

    bdens = myfld['bdens'][0,:,:]

    vecpot,xmin,xmax = vecpot2(myfld,12,3) #clip off top and bottom to remove numerical issue at boundary with vector potential calc
    vecpot = vecpot[edge:-1*edge,:]
    vecpot_smooth = gaussian_filter(vecpot,smoothlen)
#just getting bdens to be the same shape as vecpot
    #print('bdens shape : ', bdens.shape)
    #print('vecpot shape : ', vecpot.shape)
    mx = bdens.shape[1]
    mxvec = vecpot_smooth.shape[1]
    mxdiff = np.abs(mx-mxvec)
    #print('mx, mxvec',mx,mxvec)
    #print('mxdiff',mxdiff)
    #print('init shapes')
    #print(vecpot_smooth.shape,bdens.shape)
    

    myx = np.shape(bdens)[1]
    bdens = bdens[:,xmin:xmax]
    bdens = bdens[edge:-1*edge,:]

    #print(np.shape(vecpot_smooth),np.shape(bdens))

    xhlf = np.shape(vecpot)[1]/2
    #leftcut += edge
    #rightcut += edge
    topcut = -edge
    bottomcut = edge
    leftcut += xmin
    rightcut = (xmax-myx)
    #0 axis is up along y
    #1 axis is in x-direction

    #basic central difference:
    dady = (np.roll(vecpot_smooth,myroll,axis=0) -np.roll(vecpot_smooth,-1*myroll,axis=0))/myroll
    dadx = (np.roll(vecpot_smooth,myroll,axis=1) - np.roll(vecpot_smooth,-1*myroll,axis=1))/myroll
    d2ad2x = (np.roll(vecpot_smooth,myroll,axis=1)-2*vecpot_smooth+np.roll(vecpot_smooth,-1*myroll,axis=1))/myroll**2
    d2ad2y = (np.roll(vecpot_smooth,myroll,axis=0)-2*vecpot_smooth+np.roll(vecpot_smooth,-1*myroll,axis=0))/myroll**2
    d2adydx = (np.roll(dady,myroll,axis=1) - np.roll(dady,-1*myroll,axis=1))/myroll
    d2adxdy = (np.roll(dadx,myroll,axis=0) - np.roll(dadx,-1*myroll,axis=0))/myroll


    #print('shapes before smoothlen cut', dady.shape, bdens.shape)

    #cutting off stuff that rolled over the boundary and results in nonsense
    dady = dady[:,smoothlen+1:-smoothlen-2]
    dadx = dadx[:,smoothlen+1:-smoothlen-2]
    d2ad2x = d2ad2x[:,smoothlen+1:-smoothlen-2]
    d2ad2y = d2ad2y[:,smoothlen+1:-smoothlen-2]
    d2adydx = d2adydx[:,smoothlen+1:-smoothlen-2]
    d2adxdy = d2adxdy[:,smoothlen+1:-smoothlen-2]
#so first we're gonna have to loop through the first derivative and identify zeros.  
    bdens = bdens[:,smoothlen+1:-smoothlen-2]

    leftcut += smoothlen+1
    rightcut += -(smoothlen+2)

    bdens_loc =  bdens<bdens_tolerence
    dadx_zero_loc = np.abs(dadx)<zero_tolerence
    dady_zero_loc = np.abs(dady)<zero_tolerence

    #print(np.shape(bdens_loc),np.shape(dadx),np.shape(dady))
    first_deriv_zero_loc = dadx_zero_loc*dady_zero_loc*bdens_loc


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
    #remember that the returned locations are with respect to our new cut array
    return xpoint_loc_list, leftcut, rightcut, bottomcut, topcut
