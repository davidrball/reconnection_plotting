import numpy as np




def outflow(filename):
    y_vel1 = np.array(get_val_filename('v3y',filename))[0,:,:]
    ymax = y_vel1.shape[0]
    xmax = y_vel1.shape[1]
    xhlf = int(xmax/2)
    yhlf=int(ymax/2)
    outflow_profile1 = y_vel1[:,xhlf]
    bdens = np.array(get_val_filename('bdens',filename))
    bdens = bdens[0,:,:]
    intbdens = np.zeros(ymax)
    for i in range(xmax):
        intbdens += bdens[:,i]
    xhlf = int(xmax/2)
    #very spiky outflow, let's just take the third highest val or something                 
    outflow_x = range(ymax)
    bdens_up = 0
    bdens_low = 0
    for i in range(yhlf,ymax):
        if intbdens[i] > 20:
            bdens_up = i
            break
    for i in range(yhlf,0,-1):
        if intbdens[i] > 20:
            bdens_low = i
            break
    #print('bdens low, bdens up')
    #print(bdens_low, bdens_up)

    for i in range(ymax):
        if i < bdens_low or i > bdens_up:
            outflow_profile1[i] = 0
    abs_outflow = np.abs(outflow_profile1)
    highest_num = 4
    myavg = 0
    for i in range(highest_num):
        tmp = np.argmax(abs_outflow)
        myavg += np.max(abs_outflow)
        abs_outflow[tmp] = 0
        #zero out the 'highest num' most extreme                                            
    myavg /= highest_num
    outflow_rate = myavg
    print ('outflow rate : ', outflow_rate)
    return outflow_rate



def marg_x(my_array):
    #my array is 2D, we want to marginalize along the "x" direction in the code

    xlen = np.shape(my_array)[1]
    ylen = np.shape(my_array)[0]
    oneD_dens = np.zeros(ylen)
    for i in range(xlen):
        oneD_dens += my_array[:,i]/xlen
    return oneD_dens

def get_val_filename(input_string, fopen):
    str_list = list(fopen.keys())
    my_string = input_string
    if my_string in str_list:

        my_list = fopen.get(my_string)
        return my_list
    else:
        print("Value not stored, maybe you had a typo")
        print(str_list)
def list_types(fopen):
    str_list = list(fopen.keys())
    for mykey in str_list:
        print(type(np.array(fopen.get(mykey))[0]))


def decompose_spectra(myspece, myrspece, mygam):
    mylen = np.size(mygam)
    returnrspece = []
    print(np.shape(myspece), np.shape(myrspece), np.shape(mygam))
    #let's try making it so it can only occur after the peak
    argmax = np.argmax(myspece)
    #print('done')
    for i in range(argmax,mylen):
        if myspece[i] > 1.7*myrspece[i]:
            #print('done')
            pass
        elif myspece[i] != 0:# and mygam[i]*myspece[i] > 10:
        #elif myspece[i] > np.max(myspece)/100:
            #print('TRIGGERED')
            print(myspece[i], np.max(myspece))
            return np.array(mygam[i+4:]), np.array(myrspece[i+4:])
            break
def return_spece_withcut(file_list, my0_list):
    out_lec_list = []
    out_ion_list = []
    out_gam_list = []
    for file_index in range(len(file_list)):
        file = file_list[file_index]
        dgam = get_val_filename("dgam",file)
        spece = np.array(get_val_filename("spece",file))
        specp = np.array(get_val_filename("specp",file)) 
        sheetlecs = np.array(get_val_filename('speceb',file))
        sheetions = np.array(get_val_filename('specpb',file))
        spece -= sheetlecs
        specp -= sheetions
        gamarray = np.array(get_val_filename('gamma',file))
        lec_plot_array=np.zeros(spece.shape[0])
        ion_plot_array=np.zeros(specp.shape[0])

        #spectra sliced by 100s of cells
        #my0_list should be in thousands (4000, 8000 etc.)
        xsize = spece.shape[1]
        xmid = int(xsize/2.)


        my0 = my0_list[file_index]
        my0_downsampled = my0 / 100. #now in units that spectra is sliced in
        F = .2 #to either side of the sheet
        scan = int(F*my0_downsampled)
        
        xlow = int(xmid-scan)
        xup = int(xmid+scan)
        
        for i in range(xlow,xup):
            ion_plot_array += specp[:,i]
            lec_plot_array += spece[:,i]
        out_gam_list.append(gamarray)
        out_lec_list.append(lec_plot_array/dgam)
        out_ion_list.append(ion_plot_array/dgam)
    return out_gam_list, out_ion_list, out_lec_list    

        

def recon_region(dens, pdens, bdens):
    myshape = np.shape(dens)
    frac = pdens / (dens - bdens)
    cury = np.zeros(myshape)
    for i in range(myshape[0]):
        for j in range(myshape[1]):
            if frac[i][j] > .01 and frac[i][j]<.99:
                cury[i][j] +=1
    return cury

def recon_region_slice(dens,pdens,bdens):
    mysize = np.size(dens)
    frac = pdens / (dens - bdens)
    cury = np.zeros(mysize)
    for i in range(mysize):
        if frac[i]>0.1 and frac[i]<.99:
            cury[i] += 1
    return cury


def load_particle(path, partnum):
    target = open(path+str(partnum)+ '.txt','r')
    target.readline()
    target.readline()
    dict_list = []
    d = {'time':0, 'x':0, 'y':0, 'u':0, 'v':0, 'w':0, 'ex':0, 'ey':0, 'ez':0, 'bx':0, 'by':0, 'bz':0,'gamma':0}
    for line in target.readlines():
        line=line[1:-2]
        #print(line)                                                                        
        mylist = line.split(', ')
        tmplist = []
        for item in mylist:
            tmplist.append(float(item))
            #print(mylist)                                                                  
        dict_list.append(tmplist)
    d["time"] = dict_list[0]
    d["x"] = dict_list[1]
    d["y"] = dict_list[2]
    d["u"] = dict_list[3]
    d["v"] = dict_list[4]
    d["w"] = dict_list[5]
    d["ex"] = dict_list[6]
    d["ey"] = dict_list[7]
    d["ez"] = dict_list[8]
    d["bx"] = dict_list[9]
    d["by"] = dict_list[10]
    d["bz"] = dict_list[11]
    d["gamma"] = dict_list[12]
    #print(d["gamma"])
    target.close()
    return d

def load_particle_edotv(path, partnum):
    target = open(path+str(partnum)+ '.txt','r')
    target.readline()#get past header
    target.readline()#index
    target.readline()#proc
    #target.readline()
    dict_list = []
    d = {'time':0, 'u':0, 'v':0, 'w':0, 'ex':0, 'ey':0, 'ez':0, 'bx':0, 'by':0, 'bz':0,'gamma':0}
    for line in target.readlines():
        #print(line)
        line=line[1:-2]
        mylist = line.split(', ')
        tmplist = []
        for item in mylist:
            tmplist.append(float(item))
            #print(mylist)                                                                \
                                                                                            
        dict_list.append(tmplist)
    d["time"] = dict_list[0]
    d["u"] = dict_list[1]
    d["v"] = dict_list[2]
    d["w"] = dict_list[3]
    d["ex"] = dict_list[4]
    d["ey"] = dict_list[5]
    d["ez"] = dict_list[6]
    d["bx"] = dict_list[7]
    d["by"] = dict_list[8]
    d["bz"] = dict_list[9]
    d["gamma"] = dict_list[10]
    #print(d["gamma"])                                                                      
    target.close()
    return d

def load_particle_edotv_withtcs(fullname):
    target=open(fullname,'r')
    target.readline()#get past header                                                       
    target.readline()#index                                                                 
    target.readline()#proc                                                                  
    mytcs = float(target.readline()) #tcs
    #target.readline()                                                            
    finalgam = float(target.readline())

    dict_list = []
    d = {'time':0, 'u':0, 'v':0, 'w':0, 'ex':0, 'ey':0, 'ez':0, 'bx':0, 'by':0, 'bz':0,'gamma':0,'tcs':0,'finalgam':0}
    for line in target.readlines():
        #print(line)                                                                        
        line=line[1:-2]
        mylist = line.split(', ')
        tmplist = []
        for item in mylist:
            tmplist.append(float(item))
            #print(mylist)                                                                \
                                                                     

        dict_list.append(tmplist)
    d["tcs"] = mytcs
    d["finalgam"] = finalgam
    d["time"] = dict_list[0]
    d["u"] = dict_list[1]
    d["v"] = dict_list[2]
    d["w"] = dict_list[3]
    d["ex"] = dict_list[4]
    d["ey"] = dict_list[5]
    d["ez"] = dict_list[6]
    d["bx"] = dict_list[7]
    d["by"] = dict_list[8]
    d["bz"] = dict_list[9]
    d["gamma"] = dict_list[10]                                                                                            
    target.close()
    return d



def return_rspec(file_list):
    out_lec_list = []
    out_ion_list = []
    out_gam_list = []
    for file in file_list:
        dgam = get_val_filename('dgam',file)
        rspece = np.array(get_val_filename("rspece",file))
        rspecp = np.array(get_val_filename("rspecp",file))
        gamarray = np.array(get_val_filename("gamma",file))
        #print(np.shape(rspece),np.shape(rspecp),np.shape(gamarray))                        
        rspece = rspece[:,0,:]
        rspecp = rspecp[:,0,:]
        lec_plot_array = np.zeros(rspece.shape[0])
        ion_plot_array = np.zeros(rspecp.shape[0])

        for i in range(rspecp.shape[1]):
            ion_plot_array += rspecp[:,i]
        for i in range(rspece.shape[1]):
            lec_plot_array += rspece[:,i]

        out_gam_list.append(gamarray)
        out_lec_list.append(lec_plot_array/dgam)
        out_ion_list.append(ion_plot_array/dgam)
    return out_gam_list, out_ion_list, out_lec_list



def return_spec(gam_array, gammin, gammax, bin_num):
    gamma_array = gam_array - 1


    gamlist = list(gamma_array)
    tmpgamlist = list(gamma_array)
    #gammin = np.min(gamma_array)
    #gammax = np.max(gamma_array)
    gammin_log = np.log10(gammin)
    gammax_log = np.log10(gammax)

    hist_bins = np.logspace(gammin_log, gammax_log, bin_num)
    hist_array = np.zeros(bin_num)

    for i in range(bin_num-1):
        hist_low = hist_bins[i]
        hist_up = hist_bins[i+1]
        #print('bin ' +str(i))
        #print('len gamlist: '+str(len(gamlist)))
        #del_list = []
        #for j in range(len(gamlist)):
            #print(j)
            #loggam = gamlist[j]
        for loggam in gamlist:
            if hist_low < loggam < hist_up:
                hist_array[i] += 1
                #del_list.append(j)
        #for index in del_list:
            #del gamlist[index]
    return hist_bins, hist_array
