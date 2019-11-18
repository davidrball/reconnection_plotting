import matplotlib
#matplotlib.use('Agg')
import numpy as np
import h5py 
import matplotlib.pyplot as plt
from scipy.special import kn 
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'


def get_val(input_string):
    str_list = list(f.keys())
    my_string = input_string
    if my_string in str_list:   
        my_list = np.array(f.get(my_string))[0]
        return my_list
    else:
        print("argument not found, your options are:")
        print(list(f.keys()))
        return 0
        
def get_val_filename(input_string, fopen):
    str_list = list(fopen.keys())
    my_string = input_string
    if my_string in str_list:
        
        my_list = np.array(fopen.get(my_string))
        return my_list
    else:
        print("Value not stored, maybe you had a typo")
        
    
    
def mean_inflow(file):
    x_vel1 = np.array(get_val_filename("v3x",file))
    #print('shape of xvel1 before anything : ', x_vel1.shape)
    #print('shape of x velocity array : ', np.shape(x_vel1))
    x_vel1 = (x_vel1[0,:,:])
    #print('shape of xvel' , np.shape(x_vel1))
    #print('shape of xvel1 : ', x_vel1.shape)
    #plt.imshow(x_vel1, vmin = -.1, vmax = .1, origin = 'lower') 
    
    ytot = np.shape(x_vel1)[0]    
    yhlf = ytot/2
    yscan = ytot/4
    
    ymin = int(yhlf-yscan)
    ymax = int(yhlf+yscan)
    
    xtot = np.shape(x_vel1)[1]
    xhlf = xtot/2
    xscan = 100 #changing this to be a flat number since xtot will change throught he sim
    xmin = int(xhlf - xscan)
    xmax = int(xhlf + xscan)
    smallarray = x_vel1[ymin:ymax,xmin:xmax]
    #print("small array shape: ", smallarray.shape)
    #average along y direction first
    int_smallarray = np.zeros(smallarray.shape[1])

    for i in range(ymax-ymin):
        int_smallarray += smallarray[i,:]
    int_smallarray /= float(ymax-ymin) #average velocity of slice in the y direction as a function of x
    
    absarray = np.abs(int_smallarray)
    mean_inflow = np.mean(absarray)
    print('mean inflow : ',mean_inflow)
    return(mean_inflow)
    #plt.imshow(smallarray, origin = 'lower')
    #xmin_array1 = x_vel1[ymin:ymax,0,xmin:xmax]
    #mean_inflow1 = np.mean(xmin_array1)
    #mean_inflow_list1.append(mean_inflow1)

    '''
    x_vel2 = np.rot90(get_val_filename("v3x",fld_lapfile))

 
    xmin_array2 = x_vel2[ymin:ymax,0,xmin:xmax]
    mean_inflow2 = np.mean(xmin_array2)
    mean_inflow_list2.append(mean_inflow2)
    '''

'''
t0 = 1
tf = 45

base2 = "/output/flds.tot."
fldbase = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide"
b0 = "0"+base2
b1 = ".1"+base2
b2 = ".3"+base2
b3 = ".7"+base2
b4 = "1"+base2
f0 = fldbase + b0
f1 = fldbase + b1
f2 = fldbase + b2
f3 = fldbase + b3
f4 = fldbase + b4
file_list = [f0, f1, f2, f3,f4]


#file_list = ["../../tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/flds.tot.","../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/flds.tot.", "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/flds.tot.", "../../tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.7/output/flds.tot."]




inflow_list = []
file_list_len = len(file_list)
for i in range(file_list_len):
    inflow_list.append([])
target=open('sigma.3_guide_inflow_highbeta_bg1.txt','w')


for t in range(t0,tf):
    h5_list = []
    tmp_file_list = []
    #tmp_file_list=file_list #reset the tmp_file_list so we can stuff to it
    for filebase in file_list:
        tmp_file_list.append(filebase)
    t_string = str(t)
    print('t string length' , len(t_string))
    if len(t_string) == 2:
         for i in range(file_list_len):
             tmp_file_list[i] += "0"
    elif len(t_string)==1:
         for i in range(file_list_len):
             tmp_file_list[i] += "00"
    
    for filestr in tmp_file_list:
        filestr += t_string
        print(filestr)
        h5_list.append(h5py.File(filestr))


    for i in range(file_list_len):
        inflow_list[i].append(mean_inflow(h5_list[i]))
for i in range(file_list_len):
    target.write(file_list[i])
    for inflow_rate in inflow_list[i]:
        target.write(' , ' + str(inflow_rate))
    target.write('\n')
    
target.close()
'''
