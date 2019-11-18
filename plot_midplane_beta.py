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
from temp_calc import temp_ions
from beta_calc import beta_slice

#define the simulation matrix, right now it's 3x3
#first index picks out guide field strength, 2nd is temp


time_string = "035"
'''
#low beta case
bg0 =  "/home/u21/davidrball/david/tristan-mp_reconnection/16k_triggered_finals/sig.3/delgam0005/output/flds.tot."+time_string
bg1 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.1/output/flds.tot."+time_string
bg3 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.3/output/flds.tot."+time_string
bg7 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/delgam0005/bguide.7/output/flds.tot."+time_string
'''

bg0 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide0/output/flds.tot."+time_string
bg1 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.1/output/flds.tot."+time_string
bg3 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.3/output/flds.tot."+time_string
bg7 = "/home/u21/davidrball/david/tristan-mp_reconnection/guide_fields/sig.3/simulation_matrix/delgam05/bguide.7/output/flds.tot."+time_string



delgam = .05
beta = .25

beta_slice0 = beta_slice(bg0, beta, delgam, 4, 10)
beta_slice1 = beta_slice(bg1, beta, delgam, 4, 10)
beta_slice3 = beta_slice(bg3, beta, delgam, 4, 10)
beta_slice7 = beta_slice(bg7, beta, delgam, 4, 10)



plt.plot(np.log10(beta_slice0),color="Blue",label="$B_{g}=0$")
plt.plot(np.log10(beta_slice1), color="Green",label="$B_{g}=0.1$")
plt.plot(np.log10(beta_slice3), color="Red",label="$B_{g}=0.3$")
plt.plot(np.log10(beta_slice7), color="Orange",label="$B_{g}=0.7$")
plt.ylim(-1,3)
plt.legend(frameon=False)
plt.xlabel('x')
plt.ylabel('$log\\beta$')
plt.savefig('beta_slices_beta.3.png')
















