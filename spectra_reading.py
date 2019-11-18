import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 16})

filepath = "classification_txt/bguide.3_triggered_hightimeres.txt"

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

plt.plot(histbins,xpoint_spec_f,color="Blue",label='X-point')
plt.plot(histbins,merger_spec_f,color="Orange",label='Merger')
plt.plot(histbins,other_spec_f,color="Red",label='Other')
plt.legend(frameon=False,prop={'size':12})
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e0,1e4)
plt.xlim(100,1e4)
plt.xlabel('$\gamma$')
plt.ylabel('$\gamma dN/d\gamma$')
plt.savefig('bguide.3_triggered_hightimeres_test.png',dpi=300,bbox_inches='tight')
