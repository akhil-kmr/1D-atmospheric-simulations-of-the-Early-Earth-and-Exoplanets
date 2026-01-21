from photochem import EvoAtmosphere
import numpy as np
from matplotlib import pyplot as plt

plt.rcParams.update({'font.size': 15})

pc = EvoAtmosphere(
    "zahnle_earth.yaml", 
    "input/settings_100pc.yaml", 
    "input/Sun_0.0Ga.txt", 
    "100pc_PT_profile/100%PAL_Sun_0.0Ga.txt" 
)
pc.var.verbose = 0 
pc.var.atol = 1e-23 

pl = pc.production_and_loss('O3', pc.wrk.usol)

# top X reactions and their integrated loss rates
loss = pl.loss_rx[:15]
loss_rates = pl.integrated_loss[:15]
production = pl.production_rx[:15]
production_rates = pl.integrated_production[:15]

print("---Production ---")
for rx, rate in zip(production, production_rates):
    print(f"{rx} - {rate:.2e} molecules/cm^2")

print("\n--- Loss ---")
for rx, rate in zip(loss, loss_rates):
    print(f"{rx} - {rate:.2e} molecules/cm^2")

#my method

# calculating o3 column
ozone_index = pc.dat.species_names.index('O3') 
#  no. density (molecules/cm^3)
ozone_density = pc.wrk.usol[ozone_index, :] 
altitudes = pc.var.z
# density over alt to get column (molecules/cm^2)
column_density = np.trapz(ozone_density, altitudes)   
dobson_units = column_density / 2.6867e16 # Convert to DU
print(f"Ozone column: {dobson_units:.2f} Dobson Units")

#wogans method
ozone_index = pc.dat.species_names.index('O3') 
ozone_density = pc.wrk.usol[ozone_index, :] 
dz = pc.var.z[1] - pc.var.z[0]
column_density = np.sum(ozone_density * dz)
dobson_units = column_density / 2.6867e16 # Convert to DU
print(f"Ozone column Wogan: {dobson_units:.2f} Dobson Units")

#plotting
sol = pc.mole_fraction_dict()

fig1, ax1 = plt.subplots(1, 1, figsize=[6, 5]) # Production
fig2, ax2 = plt.subplots(1, 1, figsize=[6, 5]) # Loss
fig3, ax3 = plt.subplots(1, 1, figsize=[6, 5]) # Ozone Profile

# production plot
for i in range(1): 
    ax1.plot(pl.production[:, i], sol['pressure']/1e6, label=pl.production_rx[i])

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.invert_yaxis()
ax1.grid(alpha=0.4)
ax1.set_ylabel('Pressure (bar)')
ax1.set_xlabel('NO Production')
ax1.legend(ncol=1, bbox_to_anchor=(0.5, 1.05), loc='lower center')

# loss plot
for i in range(1): 
    ax2.plot(pl.loss[:, i], sol['pressure']/1e6, '--', label=pl.loss_rx[i])

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.invert_yaxis()
ax2.grid(alpha=0.4)
ax2.set_ylabel('Pressure (bar)')
ax2.set_xlabel('NO Loss')
ax2.legend(ncol=1, bbox_to_anchor=(0.5, 1.05), loc='lower center')

# o3 profile

ax3.plot(sol['O3'], sol['pressure']/1e6, color='green', label='O3')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.invert_yaxis()
ax3.grid(alpha=0.4)
ax3.set_ylabel('Pressure (bar)')
ax3.set_xlabel('Ozone Mole Fraction')
ax3.set_title(f'Ozone Profile DU)')
ax3.legend(loc='best')


plt.show()