import matplotlib 
matplotlib.use('Agg')

import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})


#beta_eff_list = [.0003, .003, .03, .01, .013, .043, .09, .093, .12]

beta_eff_list =[.0001, .001, .01, .01, .011, .02, .09, .0905, .096] 

plasmoid_list = [8,12,8,7,9,3,1,1,8]

plt.scatter(beta_eff_list,plasmoid_list)
plt.xlabel('$\\beta_{eff}$')
plt.ylabel('$N_{xpoint}$')
plt.xscale('log')
plt.xlim(1e-4,.2)
plt.savefig('beta_eff_xpoints.png',dpi=200,bbox_inches='tight')
