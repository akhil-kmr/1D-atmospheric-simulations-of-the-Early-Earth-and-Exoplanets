#run a photochemical simulation of the earth 
import numpy as np
from matplotlib import pyplot as plt
from IPython.display import clear_output
plt.rcParams.update({'font.size': 15})
from photochem import EvoAtmosphere

pc = EvoAtmosphere(
    "zahnle_earth.yaml", # Reactions & thermodynamics (not changed ever)
    "input/settings.yaml", # Settings
    "input/Sun_0.0Ga.txt", # The flux
    "input/WACCM6_profile_100pc.txt" # PT profile
)
pc.var.verbose = 0 # Turn off printing.
pc.var.atol = 1e-23 # integration absolute tolerance

#this is if you want live plotting
'''
from utils import plot_atmosphere
def find_steady_state(pc, plot=True):
    pc.initialize_robust_stepper(pc.wrk.usol) 
    while True:
        if plot:
            clear_output(wait=True)
            fig, ax = plot_atmosphere(pc, species = ['H2O','N2','O2','CO2','O3','NO','CH4','H2Oaer'])
            plt.show()
        for i in range(10):
            give_up, converged = pc.robust_step()
            if give_up or converged:
                break
        if give_up or converged:
            break
    return converged

find_steady_state(pc)
'''

converged = pc.find_steady_state()
print(converged)

surf, top = pc.gas_fluxes() # all in (molecules/cm2/s)
print(surf['N2'])
print(surf['O2'])
print(surf['CO2'])
print(surf['CH4'])
print(surf['CO'])
print(surf['N2O'])
print(surf['OCS'])
print(surf['H2'])

ozone_index = pc.dat.species_names.index('O3') #find index of O3 in the species list
ozone_density = pc.wrk.usol[ozone_index, :] #retrieving number density (molecules/cm3)
altitudes = pc.var.z
column_density = np.trapezoid(ozone_density, altitudes) #trapezoidal rule approximation to integrate the number density
dobson_units = column_density / 2.6867e16  #loschmidt constant

print(f"Ozone column: {dobson_units:.2f} Dobson Units")

#save to wherever you want
#pc.out2atmosphere_txt('100pc_PT_profile/zenithangle/100%PAL_Sun_0.0_zenithangle_58.txt',overwrite=True)