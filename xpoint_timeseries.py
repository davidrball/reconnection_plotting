import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

from tristan_funcs import get_val_filename
from tristan_funcs_vecpot import vecpot2
from localmin_func import localmin
import h5py

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})

sigma=0.3
va = np.sqrt(sigma)/np.sqrt(1+sigma)

t_start = 5
t_final = 45

scan = 10
interval = 2000
istep = 12
matrix00 = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam00005/output/flds.tot."
matrix01 = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/flds.tot."
matrix02 =  "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam005/output/flds.tot."
matrix03 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide0/output/flds.tot." 


matrix10 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam00005/bguide.1/output/flds.tot."
matrix11 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/flds.tot."
matrix12 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.1/output/flds.tot."
matrix13  = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.1/output/flds.tot."

#assign unfinished sims as "pass"                                                           


matrix20 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam00005/bguide.3/output/flds.tot."
matrix21 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/flds.tot."
matrix22 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.3/output/flds.tot."
matrix23 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.3/output/flds.tot."

#matrix23 = "pass"                                                                          

matrix30 = "pass"
matrix31 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.7/output/flds.tot."
matrix32 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.7/output/flds.tot."
matrix33 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.7/output/flds.tot."



#matrix_namelist = [[matrix00, matrix01, matrix02], [matrix10,matrix11,matrix12], [matrix20,matrix21,matrix22]]

matrix_namelist = [[matrix00, matrix01, matrix02,matrix03], [matrix10,matrix11,matrix12,matrix13], [matrix20,matrix21,matrix22,matrix23],[matrix30,matrix31,matrix32,matrix33]]




xpoint_countlist = [[[],[],[],[]], [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]

for i in range(len(matrix_namelist)):
    for j in range(len(matrix_namelist[0])): #assuming that the matrix is square
        if matrix_namelist[i][j] == "pass":
            pass
        else:
            print(i,j)
            myname = matrix_namelist[i][j]
            for t in range(t_start, t_final):
                t_str = str(t)
                print(t_str)
                if len(t_str) == 1:
                    tmpname = myname + "00" + t_str
                elif len(t_str)==2:
                    tmpname = myname + "0" + t_str
                myf = h5py.File(tmpname,'r')
                vecpot = vecpot2(myf,12,3)
                xmid = np.shape(vecpot)[1]/2
                vecpot_slice = vecpot[:,xmid]

                ylen = np.shape(vecpot)[0]
                ylenhlf = ylen/2
                time_offset = 4
                time_delay = 3
                myedge = max(int(ylenhlf - (va * .45 * interval * (t-time_offset)/istep / time_delay)), scan+1)

                my_localmin_list = localmin(vecpot_slice,scan,myedge)
                num_xpoints=len(my_localmin_list)
                xpoint_countlist[i][j].append(num_xpoints)
                myf.close()
for i in range(len(matrix_namelist)):
    for j in range(len(matrix_namelist[0])):
        if matrix_namelist[i][j] == "pass":
            pass
        else:
            my_timeseries = xpoint_countlist[i][j]
            name = 'xpoint_txt/edge/matrix' + str(i) + str(j)
            target = open(name,'w')
            target.write(str(my_timeseries))
            target.close()

        '''
        if i==0:
            col = "Blue"
        elif i==1:
            col = "Red"
        elif i==2:
            col = "Orange"
        if j ==0:
            mylin = "solid"
        elif j==1:
            mylin = "dashed"
        elif j==2:
            mylin = "dotted"
        
        plt.plot(xpoint_countlist[i][j],color=col,linestyle=mylin)
        '''
#plt.savefig('xpoint_time_count.png',bbox_inches='tight')


