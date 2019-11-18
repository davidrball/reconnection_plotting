import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
from tristan_funcs import recon_region, recon_region_slice
#testfld = "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.7/output/flds.tot.040"


#we want to calculate sheet thicknesses based on bdens

#calculates thickness given a single slice of the mixing parameter across the sheet (just at single location)
def slice_thickness(cury_slice):
    mylen = np.size(cury_slice)
    lowind = 0
    highind = 0
    #loop from the beginning
    for i in range(mylen):
        if cury_slice[i] == 1:
            lowind = i
            break

    #loop from end
    for i in range(mylen-1,0,-1):
        if cury_slice[i] == 1:
            highind = i
            break
    thickness = highind-lowind
    return thickness

def return_thickness_list(file_string):
    myf = h5py.File(file_string)


    dens = np.array(myf['dens'][0,:,:])
    print(np.shape(dens))
    #let's say we only want the central third of the box, doing this operation in code units (so y is along the sheet)
    ylen = np.shape(dens)[0]
    ylow = int(ylen/4.)
    yup = 3*ylow
    #this gives central half the box (from 1/4 to 3/4)
    #ylow = 0
    #yup = ylen
    dens = dens[ylow:yup]
    bdens = np.array(myf['bdens'][0,ylow:yup,:])
    pdens = np.array(myf['pdens'][0,ylow:yup,:])
    cury = recon_region_slice(dens,pdens,bdens)

    myshape = np.shape(cury)
    xlen = myshape[0]
    ylen = myshape[1]

    print(np.shape(cury))
    #cury[Y_val,:] gives the profile sliced along x, need to just identify lowest value w/ a 1 in it, highest value w/ a 1 in it
    thickness_list = []
    #testslice = cury[0,:]
    print(ylen)
    for i in range(xlen):
        myslice = cury[i,:]
        thickness = slice_thickness(myslice)
        thickness_list.append(thickness)
    thickarr = np.array(thickness_list)
    meanthick = np.mean(thickarr)
    stdthick = np.std(thickarr)
    #print(meanthick, stdthick)
    return thickness_list, meanthick, stdthick
#thickness_list = return_thickness_list(testfld)
#plt.hist(thickness_list,bins=40)
#plt.savefig('thick_hist.png')
#plt.close()

def return_thickness_slice(file_string):
    myf = h5py.File(file_string)
    dens = np.array(myf['dens'][0,:,:])
    #print(np.shape(dens))
    #let's say we only want a slice at +- 1/3, doing this operation in code units (so y is along the sheet)                                                                 
    ylen = np.shape(dens)[0]
    ylow = int(ylen/3.)
    yup = 2*ylow
    #this gives central half the box (from 1/4 to 3/4)                                      
    #ylow = 0                                                                               
    #yup = ylen                                                                             
    densup = dens[yup]
    bdensup = np.array(myf['bdens'][0,yup,:])
    pdensup = np.array(myf['pdens'][0,yup,:])
    #print(np.shape(densup))
    #print(np.shape(bdensup))

    curyup = recon_region_slice(densup,pdensup,bdensup)
    
    #print('got to here')
    denslow = dens[ylow]
    bdenslow = np.array(myf['bdens'][0,ylow,:])
    pdenslow = np.array(myf['pdens'][0,ylow,:])
    curylow = recon_region_slice(denslow,pdenslow,bdenslow)



    mysize = np.size(curylow)
 
    thickness_list = []
    #testslice = cury[0,:]                                                                  
    #print(ylen)
    mysliceup = curyup
    thicknessup = slice_thickness(mysliceup)
    thickness_list.append(thicknessup)
    myslicelow = curylow
    thicknesslow = slice_thickness(myslicelow)
    thickness_list.append(thicknesslow)

    thickarr = np.array(thickness_list)
    print('thicknesses : ', thickness_list)
    return thickness_list

