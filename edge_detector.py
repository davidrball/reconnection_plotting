import numpy as np


#define the convolution kernel
#mykernel = np.array([-1, -1, -1, 8, -1, -1, -1])




#apply the convolution
#just 1d for our purposes
def convolve(image, kernel,edge_scan):
    imagelen = np.size(image)
    new_image = np.zeros(imagelen)


    kernellen = np.size(kernel)
    kernelmid = kernellen/2
    #use edge scan to define upper and lower bounds of convolution
    #edge scan must be greater than kernellen/2
    for i in range(edge_scan, imagelen-edge_scan):
        #i is the point we're applying the convolution to
        imageslice = image[i-kernelmid:i+kernelmid]
        kernel_times_image = imageslice * kernel
        mysum = np.sum(kernel_times_image)
        new_image[i] = mysum
    return new_image

#testing
#testimage = np.ones(100)
#test_convolve = convolve(testimage, mykernel, 10)
#print(test_convolve)
    
