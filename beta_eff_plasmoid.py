import matplotlib 
matplotlib.use('Agg')

import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 14})


beta_eff_list = [.0003, .003, .03, .01, .013, .043, .09, .093, .12]
plasmoid_list = [7,6,5,5,8,2,0,0,6]

plt.scatter(beta_eff_list,plasmoid_list)
plt.xlabel('$\\beta_{eff}$')
plt.ylabel('$N_{plasmoid}$')
plt.xscale('log')
plt.xlim(1e-4,.2)
plt.savefig('beta_eff_plasmoids.png',dpi=200,bbox_inches='tight')
