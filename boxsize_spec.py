import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
import numpy as np
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 15})

#inputs
t0=45
tf=46

mybase = "../../tristan_acc-mec_Ez/boxsize_test/triggered/"
endbase = "/output/spect."

twok = mybase + "2k" + endbase
fourk = mybase + "4k" + endbase
sixk = mybase + "6k" + endbase
eightk = "../../tristan_acc-mec_Ez/8k_bguide.3_triggered_stride1_thresh2/output/spect."


size_list = [2., 4., 6., 8.]
filestr_list = [twok,fourk,sixk,eightk]
reference_size = 8.
offset = 5.
col_list = ["C1","C2","C3","C4"]
label_list = []
for i in range(len(size_list)):
    label_list.append("Box Length : "+str(int(size_list[i])) + "k cells")


norm_arr = np.array(size_list)/reference_size

listlen = len(filestr_list)

for myt in range(t0, tf):
    print(myt)
    t_list = []
    for i in range(listlen):
        mynorm = norm_arr[i]
        adjusted_t = int(myt * mynorm + offset)
        adjusted_t_str = "%03d"%adjusted_t
        t_list.append(adjusted_t_str)
    #print(t_list)
    for i in range(listlen):
        myfile_str = filestr_list[i] + t_list[i]
        print(myfile_str)
        f_spec = h5py.File(myfile_str)
        mycol = col_list[i]
        mylabel = label_list[i]
    

        mygam = np.array(f_spec['gamma'])
        rspec=np.array(f_spec['rspece'])[:,0,:]
        tspec=np.array(f_spec['tspece'])[:,0,:]
        gambins,spacebins = np.shape(tspec)
        mytspec = np.zeros(gambins)
        myrspec = np.zeros(gambins)
        for j in range(gambins):
            myrspec[j]+=np.sum(rspec[j][:])
            mytspec[j]+=np.sum(tspec[j][:])
        plt.plot(mygam,mygam*myrspec,color=mycol,label=mylabel)
        #plt.plot(mygam,mygam*mytspec,color=mycol,label=mylabel)
        f_spec.close()
plt.xscale('log')
plt.yscale('log')
    
plt.xlabel('$\gamma-1$')
plt.ylabel('$(\gamma-1) dN/d\gamma$')
    
plt.xlim(1e0,1e4)
plt.ylim(1e1,1e8)
    
plt.legend(frameon=False,prop={'size':12})
plt.savefig('testing_boxsize_spec.png',bbox_inches='tight',dpi=300)
