import numpy as np

#brute force method of identifying local mins

#testarr = np.array([1,2,3,4,5,4,3,2,1,2,3,4,5,6,7,6,5,6,7,8,9,10])


def localmin(array, stencil_half_size,edge): #only for 1D arrays, edge must be > stencil_half_size
    
    local_min_list = []
    for i in range(edge,np.size(array)-edge):
        myval = array[i]
        arrayslice = array[i-stencil_half_size:i+stencil_half_size+1]
        #print(myval, arrayslice)
        bool_array = myval <= arrayslice
        if False in bool_array:
            pass
        else:
            #print(i)
            local_min_list.append(i)
    return local_min_list

#localmin_list = localmin(testarr,2, 3)

#print(localmin_list)


def localmax(array, stencil_half_size,edge): 
    local_min_list = []
    for i in range(edge,np.size(array)-edge):
        myval = array[i]
        arrayslice = array[i-stencil_half_size:i+stencil_half_size+1]
        #print(myval, arrayslice)
        bool_array = myval >= arrayslice
        if False in bool_array:
            pass
        else:
            #print(i)
            local_min_list.append(i)
    return local_min_list
        
def localmin_wrap(array, stencil_half_size): #only for 1D arrays, edge mus\t be > stencil_half_size
    local_min_list = []


    for i in range(0,np.size(array)):
        if i - stencil_half_size < 0:
            #print('lower boundary triggered')
            newlist = []
            underflow = i - stencil_half_size
            #print('underflow ' , underflow)
            underflow_array = array[np.size(array)+underflow:]
            #print('size of underflow array : ', np.size(underflow_array))            
            for j in range(np.size(underflow_array)):
                #if underflow_array[j] > 0 :
                newlist.append(underflow_array[j])
            reg_array = array[:i+stencil_half_size+1]
            for j in range(np.size(reg_array)):
                #if reg_array[j] > 0:
                newlist.append(reg_array[j])
            arrayslice = np.array(newlist)
            #print('underflow arrayslice ' , np.size(arrayslice))
        elif i + stencil_half_size > np.size(array):
            #print('upper boundary triggered')
            newlist = []
            overflow = i+ stencil_half_size - np.size(array)
            overflow_array = array[:overflow]
            reg_array = array[i-stencil_half_size:]
            for j in range(np.size(reg_array)):
                if reg_array[j]>0:
                    newlist.append(reg_array[j])
            for j in range(np.size(overflow_array)):
                if overflow_array[j] > 0:

                    newlist.append(overflow_array[j])
            arrayslice = np.array(newlist)
        elif i-stencil_half_size > 0 and i+stencil_half_size < np.size(array):
            arrayslice = array[i-stencil_half_size:i+stencil_half_size+1]
        #print('size of arrayslice ' , np.size(arrayslice))
        #print(myval, arrayslice)
        myval = array[i]
        bool_array = myval <= arrayslice
        if False in bool_array:
            pass
        else:
            #print(i)
            local_min_list.append(i)
    return local_min_list        


def localmin_2d(array, scan):
    localmin_x_list = []
    localmin_y_list = []

    myshape = np.shape(array)
    
    for i in range(scan,myshape[0]-scan):
        for j in range(scan,myshape[1]-scan):
            
            myvecpot = array[i][j]
            #define subarray
            subarray = array[i-scan:i+scan,j-scan:j+scan]
            #print(np.shape(subarray))
            #now test if myvecpot is the minima of this array
            #print(myvecpot, np.min(subarray))
            #if myvecpot == np.amin(subarray):
            #    localmin_x_list.append(j)
            #    localmin_y_list.append(i)
            bool_array = myvecpot <= subarray
            if False in bool_array:
                pass
            else:
            #print(i)                                                          
                localmin_x_list.append(j)
                localmin_y_list.append(i)
    return localmin_x_list, localmin_y_list

def localmin_2d_bdens(array, scan,bdens,thresh):
    localmin_x_list = []
    localmin_y_list = []

    myshape = np.shape(array)

    for i in range(scan,myshape[0]-scan):
        for j in range(scan,myshape[1]-scan):

            myvecpot = array[i][j]
            #define subarray                                                                                                                                                             
            subarray = array[i-scan:i+scan+1,j-scan:j+scan+1]
            #print(np.shape(subarray))                                                                                                                                                   
            #now test if myvecpot is the minima of this array                                                                                                                            
            #print(myvecpot, np.min(subarray))                                                                                                                                           
            #if myvecpot == np.amin(subarray):                                                                                                                                           
            #    localmin_x_list.append(j)                                                                                                                                               
            #    localmin_y_list.append(i)                                                                                                                                               
            bool_array = myvecpot <= subarray
            if False in bool_array:
                pass
            elif bdens[i][j] < thresh:
            #print(i)                                                                                                                                                                    
                localmin_x_list.append(j)
                localmin_y_list.append(i)
    return localmin_x_list, localmin_y_list
