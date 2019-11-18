import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 16})

#betaeff_list = [.0003,.003,.03,.01,.013,.043,.09,.093,.12]
#betaeff_list = [.0003, .003, .03, .25, .01, .013, .043, .34,.09, .093, .12, .42, .49, .493, .52, .82]

#betaeff_list = [.0003, .003, .03, .25, .01, .01, .02, .1, .09, .09, .0965,.155, .49, .49, .494, .52 ]

betaeff_list = [.0003, .003, .03, .25, .01, .013, .043, .33, .083, .085, .1, .38, .33, .33, .35, .55]

bguide_list = [0,0,0,0,.1,.1,.1,.1,.3,.3,.3,.3,.7,.7,.7,.7,1,1,1,1]

actual_betalist = [.0003, .003, .03, .3, .0003, .003, .03, .3, .0003, .003, .03, .3, .0003, .003, .03, .3, .0003, .003, .03, .3, .0003, .003, .03, .3]

avg_list = []


for i in range(5):
    for j in range(5):
        print(i, j)
        if i==4 and j!=3 or j==4:
            pass
        else:
            myname = "plasmoid_xpoint_txt/"+str(i)+str(j)
            target = open(myname,'r')
            mylist = target.readline()[1:-1].split(', ')
            for k in range(len(mylist)):
                mylist[k] = int(mylist[k])
            myarr = np.array(mylist[5:45])
            mymean = np.mean(myarr)
            avg_list.append(mymean)
            myind = 4*i + j
            print(mymean)
            #print(mylist)
        
            if j==0:                                                           
                col = "Blue"                                                  
            elif j==1:                                                         
                col = "Green"                                                  
            elif j==2:                                                         
                col = "Red"
            elif j==3:
                col = "Orange"
            if i ==0:                               
            #mylin = "solid"
                mymark = 'o'
            elif i==1:                                                     
            #mylin = "dashed"
                mymark = 'x'
            elif i==2:                                                      
            #mylin = "dotted"
                mymark = 'D'
            elif i==3 :
                mymark = "*"
            elif i==4:
                mymark = "s"
        #plt.scatter(actual_betalist[myind],mymean,marker=mymark,color=col)  
        #plt.scatter(betaeff_list[myind],mymean,marker=mymark,color=col)
            plt.scatter(bguide_list[myind],mymean,marker=mymark,color=col)
#plt.xlim(10,50)
#plt.savefig('xpoint_time_count.png',bbox_inches='tight')

#plt.xscale('log')
plt.ylabel('$\overline{N}_{xpoint+plasmoid}$')
#plt.xlabel('$\\beta_{eff}$')
#plt.xlim(1e-4,1)
plt.xlabel('$B_{g}/B_{0}$')
plt.xlim(-.2,1)
plt.ylim(0,20)
plt.savefig('sig.3_bguide_vs_xpoint.png',bbox_inches='tight',dpi=300)
plt.close()

'''
for i in range(len(betaeff_list)):
    if i <= 3 or i%4 ==0:
        plt.scatter(betaeff_list[i],avg_list[i], color=actual_betalist[i])


plt.xscale('log')
plt.ylabel('$\overline{N}_{xpoint+plasmoid}$')
plt.xlabel('$\\beta_{eff}$')
plt.xlim(1e-4,1)
plt.ylim(0,20)
plt.savefig('0ind_betaeff_test.png',bbox_inches='tight',dpi=300)
'''
