from photochem import EvoAtmosphere
import numpy as np
from matplotlib import pyplot as plt

pc = EvoAtmosphere(
    "zahnle_earth.yaml", 
    "input/settings_100pc.yaml", 
    "input/Sun_0.0Ga.txt", 
    "100pc_PT_profile/100%PAL_Sun_0.0Ga.txt" 
)
pc.var.verbose = 0 
pc.var.atol = 1e-23 


pl = pc.production_and_loss('O2', pc.wrk.usol)
sol = pc.mole_fraction_dict()
pressure_hpa = sol['pressure'] / 1e2   #Pa to hPa

# indiex for the O2 photolysis 

idx_o_o = -1
idx_o_o1d = -1

for i, rx in enumerate(pl.loss_rx):
    if rx == 'O2 + hv => O + O':
        idx_o_o = i
    elif rx == 'O2 + hv => O + O1D':
        idx_o_o1d = i
    
# every O2 loss gives 2 odd ox 
ox_prod_cm3 = 2 * pl.loss[:, idx_o_o] + 2 * pl.loss[:, idx_o_o1d]

# to molecules/m3/s
ox_prod_m3 = ox_prod_cm3 * 1e6


plt.rcParams.update({'font.size': 13})
fig, ax = plt.subplots(1, 1, figsize=[8, 6])

ax.plot(ox_prod_m3, pressure_hpa)

ax.set_xscale('log')
ax.set_yscale('log')
ax.invert_yaxis()
ax.grid(alpha=0.4)

ax.set_xlim(1e3, 1e14)
ax.set_ylim(10e3,10e-5)

ax.set_ylabel('Pressure (hPa)')
ax.set_xlabel('Production (molecules/m$^3$/s)')


plt.tight_layout()
plt.show()
