import matplotlib
matplotlib.use('Agg')
import numpy as np
import h5py
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})
plt.set_cmap("viridis")

from tristan_funcs import get_val_filename
from overdens_calc import return_avg_dens
bguide_list = [0,0.003,0.006,0.009,.1,.103,.106,.109,.3,.303,.306,.309,.7,.703,.706,.709]

#define the simulation matrix, right now it's 3x3
#first index picks out guide field strength, 2nd is temp
t0 = 20
tf = 22

dens_superlist = [[[],[],[],[]],[[],[],[],[]],[[],[],[],[]],[[],[],[],[]]]
std_superlist = dens_superlist
for t in range(t0,tf):
    #guide field 0, delgam=.00005, .0005, .005
    time_string = '%03d' % t
    print(time_string)
    matrix00 = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam00005/output/flds.tot."+time_string
    matrix01 = "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/flds.tot."+time_string
    matrix02 =  "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam005/output/flds.tot."+time_string
    matrix03 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide0/output/flds.tot." + time_string


    matrix10 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam00005/bguide.1/output/flds.tot."+time_string
    matrix11 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/flds.tot."+time_string
    matrix12 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.1/output/flds.tot."+time_string
    matrix13  = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.1/output/flds.tot." + time_string

#assign unfinished sims as "pass"


    matrix20 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam00005/bguide.3/output/flds.tot."+time_string
    matrix21 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/flds.tot."+time_string
    matrix22 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.3/output/flds.tot."+time_string
    matrix23 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.3/output/flds.tot." + time_string

#matrix23 = "pass"

    matrix30 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam00005/bguide.7/output/flds.tot."+time_string
    matrix31 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.7/output/flds.tot." + time_string
    matrix32 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam005/bguide.7/output/flds.tot."+time_string
    matrix33 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.7/output/flds.tot." + time_string

    matrix_namelist = [[matrix00, matrix01, matrix02,matrix03], [matrix10,matrix11,matrix12,matrix13], [matrix20,matrix21,matrix22,matrix23],[matrix30,matrix31,matrix32,matrix33]]

    

    for i in range(len(matrix_namelist)):
        for j in range(len(matrix_namelist[0])): #assuming that the matrix is squar                
            print(i,j)
            myind = 4*i + j
            myname = matrix_namelist[i][j]
            mydens,stddens = return_avg_dens(myname)
            if j==0:
                col = "Blue"
            elif j==1:
                col = "Green"                                                 
            elif j==2:
                col = "Red"
            elif j==3:
                col = "Orange"
            if i ==0:                                                    
                mymark = 'o'
            elif i==1:
                mymark = 'x'
            elif i==2:
                mymark = 'D'
            elif i==3 :
                mymark = "*"

            #plt.errorbar(bguide_list[myind], mydens,yerr=stddens,marker=mymark, color=col)
            #plt.scatter(bguide_list[myind],mydens,marker=mymark,color=col)
            dens_superlist[i][j].append(mydens)
for i in range(4):
    for j in range(4):
        myind = 4*i + j
        if j==0:
            col = "Blue"
        elif j==1:
            col = "Green"
        elif j==2:
            col = "Red"
        elif j==3:
            col = "Orange"
        if i ==0:
            mymark = 'o'
        elif i==1:
            mymark = 'x'
        elif i==2:
            mymark = 'D'
        elif i==3 :
            mymark = "*"
        plt.scatter(bguide_list[myind],np.mean(np.array(dens_superlist[i][j])),marker=mymark,color=col)
#plt.ylabel('Mean Thickness (downsampled cells)')
plt.ylabel('Average Density (particles per cell)')
plt.xlabel('$B_{g}/B_{0}$')
plt.xlim(-.2,1)
plt.savefig('sig.3_time_avg_dens.png',bbox_inches='tight',dpi=300)
plt.close()









