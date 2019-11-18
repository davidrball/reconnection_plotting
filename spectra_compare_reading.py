import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 16})

filepath = "classification_txt/bguide.3_triggered_spec.txt"
filepath2 = "classification_txt/bguide.3_untriggered_spec.txt"


target = open(filepath,'r')
target.readline() #reading histbins
histbins = target.readline().split(' , ')[:-2]
target.readline() #reading xpoint spec
xpoint_spec = target.readline().split(' , ')[:-2]
target.readline() #reading merger spec
merger_spec = target.readline().split(' , ')[:-2]
target.readline() #reading other spec
other_spec = target.readline().split(' , ')[:-2]

#convert lists of strings to lists of floats
histbins_f = []
xpoint_spec_f = []
merger_spec_f = []
other_spec_f = []
mylen = len(histbins)

for i in range(mylen):
    histbins_f.append(float(histbins[i]))
    xpoint_spec_f.append(float(xpoint_spec[i]))
    merger_spec_f.append(float(merger_spec[i]))
    other_spec_f.append(float(other_spec[i]))

xpoint_spec_f = np.array(xpoint_spec_f)
merger_spec_f = np.array(merger_spec_f)
other_spec_f = np.array(other_spec_f)
histbins = np.array(histbins_f)
target.close()

target2 = open(filepath2,'r')
target2.readline() #reading histbins                                                   
histbins2 = target2.readline().split(' , ')[:-2]
target2.readline() #reading xpoint spec
xpoint_spec2 = target2.readline().split(' , ')[:-2]
target2.readline() #reading merger spec                                                
merger_spec2 = target2.readline().split(' , ')[:-2]
target2.readline() #reading other spec                                 
other_spec2 = target2.readline().split(' , ')[:-2]

#convert lists of strings to lists of floats                                                                                                                                             
histbins_f2 = []
xpoint_spec_f2 = []
merger_spec_f2 = []
other_spec_f2 = []
mylen2 = len(histbins2)



for i in range(mylen2):
    histbins_f2.append(float(histbins2[i]))
    xpoint_spec_f2.append(float(xpoint_spec2[i]))
    merger_spec_f2.append(float(merger_spec2[i]))
    other_spec_f2.append(float(other_spec2[i]))

xpoint_spec_f2 = np.array(xpoint_spec_f2)
merger_spec_f2 = np.array(merger_spec_f2)
other_spec_f2 = np.array(other_spec_f2)
histbins2 = np.array(histbins_f2)

gamarg = np.argmin(np.abs(histbins-1000))
gamarg2 = np.argmin(np.abs(histbins2-1000))
print('gam args : ', gamarg, gamarg2)


xpoint_gamarg = xpoint_spec_f[gamarg]
merger_gamarg = merger_spec_f[gamarg]
ratio = xpoint_gamarg / merger_gamarg


xpoint_gamarg2 = xpoint_spec_f2[gamarg2]
merger_gamarg2 = merger_spec_f2[gamarg2]
ratio2 = xpoint_gamarg2 / merger_gamarg2
print('triggered xpoint ratio : ',ratio)
print('untriggered xpoint ratio : ', ratio2)


trig_untrig_ratio = ratio / ratio2
print('triggered vs untriggered ratio : ', trig_untrig_ratio)





plt.plot(histbins,histbins*xpoint_spec_f,color="Blue",label='X-point')
plt.plot(histbins,histbins*merger_spec_f,color="Orange",label='Merger')
plt.plot(histbins,histbins*other_spec_f,color="Red",label='Other')

plt.plot(histbins2,histbins2*xpoint_spec_f2,color="Blue",label='Untrig X-point',linestyle='--')
plt.plot(histbins2,histbins2*merger_spec_f2,color="Orange",label='Untrig Merger',linestyle='--')
plt.plot(histbins2,histbins2*other_spec_f2,color="Red",label='Other',linestyle='--')




plt.legend(frameon=False,prop={'size':12})
plt.xscale('log')
plt.yscale('log')
#plt.ylim(1e5,1e8)
#plt.xlim(100,1e4)
plt.xlabel('$\gamma$')
plt.ylabel('$\gamma dN/d\gamma$')
plt.savefig('spectra_compare.png',dpi=300,bbox_inches='tight')
