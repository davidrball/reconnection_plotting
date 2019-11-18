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
from thickness_calc import return_thickness_list

bguide_offset_list = [0,0.003,0.006,0.009,.1,.103,.106,.109,.3,.303,.306,.309,.7,.703,.706,.709]
bguide_list = [0,0,0,0,.1,.1,.1,.1,.3,.3,.3,.3,.7,.7,.7,.7]
#define the simulation matrix, right now it's 3x3
#first index picks out guide field strength, 2nd is temp
t0 = 10
tf = 45

thick_superlist = [[[],[],[],[],[]],[[],[],[],[]],[[],[],[],[]],[[],[],[],[]]]

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
            tmp_thick_list, meanthick, stdthick = return_thickness_list(myname)
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

            #plt.errorbar(bguide_list[myind], meanthick ,yerr=stdthick,marker=mymark, color=col)
            #plt.scatter(bguide_list[myind],stdthick,marker=mymark,color=col)
            thick_superlist[i][j].append(meanthick)
for i in range(len(matrix_namelist)):
    for j in range(len(matrix_namelist[0])):
         tmp_thick_list, meanthick, stdthick = return_thickness_list(myname)
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
         myind = 4*i + j
         plt.scatter(bguide_list[myind], np.mean(np.array(thick_superlist[i][j])), marker=mymark, color=col)
#plt.ylabel('Mean Thickness (downsampled cells)')
plt.ylabel('Mean thickness')
plt.xlabel('$B_{g}/B_{0}$')
plt.xlim(-.2,1)
plt.savefig('sig.3_thickness_timeavg.png',bbox_inches='tight',dpi=300)
plt.close()









