"""
Early Earth 1D model comparisons — main script.

Loads model output, runs diagnostics, and generates figures.
All functions and shared constants: early_earth_lib.py (imported below).

File layout
-----------
1. early_earth_lib.py  — imports, paths, O3 budget, plotting, I/O, column integrals
2. Early_Earth.py      — data loading, then analysis cells (#%%) in execution order

Run this file cell-by-cell in Spyder, or execute the whole script.
"""

from early_earth_lib import thick_axes, Read_O3_Run, V_O3_col_z_trapz
from early_earth_lib import _o3_loss_fraction_figure
from early_earth_lib import *
import matplotlib.pyplot as plt
from early_earth_lib import Photochem_O3_col, plot_o3_nox_hox_fraction_vulcan_vs_photo
from early_earth_lib import waccm_o3_catalytic_budget, waccm_has_o3_loss_diagnostics
import pandas as pd
import xarray as xr
import numpy as np

#%% import file for pressure and gaussian weights
File =  "/Users/gregcooke/H_escape/b.e21.BWma1850.f19_g17.PC_b.SSPO.016.cam.h0.0320-0320.nc" #file name

Gw_file = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

gw = Gw_file.gw.values
Pressure = Gw_file.lev.values
#%%
def LWAV(DS, time = False):
    try:
        DS = DS.mean(dim='lon')
    except:
        DS = DS
    if (time == False):
        try:
            DS = DS.mean(dim='time')
        except:
            DS = DS
    DS = gw*DS #multiply by gaussian weights
    try:
        DS = np.sum(DS,axis=1)/gw.sum() #weighted mean
    except:
        DS = np.sum(DS,axis=0)/gw.sum() #weighted mean 
    return DS #return a numpy array that has been modified above

#%% User notes
"""
Instructions for user to go here
"""

figure = plt.figure(figsize = (8,5))
lw = 2

plt.title('Ozone profiles', fontsize = 15, weight = 'bold')

#%% Read in VULCAN files

#%% VULCAN PCb runs
PCb_V, PCb_V_spec = Read_O3_Run(file_path='PCb_1e12s_482SZA_WBC_WPT_1rtol.vul')
PCb_50pc_V, PCb_50pc_V_spec = Read_O3_Run(file_path='PCb_50pc_o2_1e12s_482SZA_WBC_WPT_1rtol.vul')
PCb_10pc_V, PCb_10pc_V_spec = Read_O3_Run(file_path='PCb_10pc_o2_1e12s_482SZA_WBC_WPT_1rtol.vul')
PCb_5pc_V, PCb_5pc_V_spec = Read_O3_Run(file_path='PCb_5pc_o2_1e12s_482SZA_WBC_WPT_1rtol.vul')
PCb_1pc_V, PCb_1pc_V_spec = Read_O3_Run(file_path='PCb_1pc_o2_1e12s_482SZA_WBC_WPT_1rtol.vul')
PCb_05pc_V, PCb_1pc_V_spec = Read_O3_Run(file_path='PCb_0.5pc_o2_1e12s_482SZA_WBC_WPT_1rtol.vul')

#%% VULCAN regular oxygen cases

One50_pc_V_60SZA, One50_pc_spec_60SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_60SZA_WPT_1rtol.vul')
One50_pc_V_482SZA, One50_pc_spec_482SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_48.2SZA_WPT_1rtol.vul')
One50_pc_V_45SZA, One50_pc_spec_45SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_45SZA_WPT_1rtol.vul')

PI_V_60SZA, PI_spec_60SZA = Read_O3_Run(file_path='Earth_1e12s_60SZA_WPT_1rtol.vul')
PI_V_482SZA, PI_spec_482SZA = Read_O3_Run(file_path='Earth_1e12s_48.2SZA_WPT_1rtol.vul')
PI_V_45SZA, PI_spec_45SZA = Read_O3_Run(file_path='Earth_1e12s_45SZA_WPT_1rtol.vul')

Fifty_pc_V_60SZA, Fifty_pc_spec_60SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_60SZA_WPT_1rtol.vul')
Fifty_pc_V_482SZA, Fifty_pc_spec_482SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_48.2SZA_WPT_1rtol.vul')
Fifty_pc_V_45SZA, Fifty_pc_spec_45SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Ten_pc_V_60SZA, Ten_pc_spec_60SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_60SZA_WPT_1rtol.vul')
Ten_pc_V_482SZA, Ten_pc_spec_482SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_48.2SZA_WPT_1rtol.vul')
Ten_pc_V_482SZA_T, Ten_pc_spec_482SZA_T = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_48.2SZA_WPT_Temp_1rtol.vul')
Ten_pc_V_45SZA, Ten_pc_spec_45SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Five_pc_V_60SZA, Five_pc_spec_60SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_60SZA_WPT_1rtol.vul')
Five_pc_V_482SZA, Five_pc_spec_482SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_48.2SZA_WPT_1rtol.vul')
Five_pc_V_45SZA, Five_pc_spec_45SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_45SZA_WPT_1rtol.vul')

One_pc_V_60SZA, One_pc_spec_60SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_60SZA_WPT_1rtol.vul')
One_pc_V_482SZA, One_pc_spec_482SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_48.2SZA_WPT_1rtol.vul')
One_pc_V_482SZA_T, One_pc_spec_482SZA_T = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_48.2SZA_WPT_Temp_1rtol.vul')
One_pc_V_45SZA, One_pc_spec_45SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Z5_pc_V_60SZA, Z5_pc_spec_60SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_60SZA_WPT_1rtol.vul')
Z5_pc_V_482SZA, Z5_pc_spec_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_WPT_1rtol.vul')
Z5_pc_V_45SZA, Z5_pc_spec_45SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Z1_pc_V_60SZA, Z1_pc_spec_60SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_60SZA_WPT_1rtol.vul')
Z1_pc_V_482SZA, Z1_pc_spec_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_WPT_1rtol.vul')
Z1_pc_V_482SZA_T, Z1_pc_spec_482SZA_T = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_WPT_Temp_1rtol.vul')
Z1_pc_V_45SZA, Z1_pc_spec_45SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_45SZA_WPT_1rtol.vul')

#%% VULCAN cases with different methane
One50_pc_V_10xCH4_482SZA, One50_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_48.2SZA_10xCH4_WPT_1rtol.vul')
One50_pc_V_5xCH4_482SZA, One50_pc_spec_5xCH4_482SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_48.2SZA_5xCH4_WPT_1rtol.vul')
One50_pc_V_05xCH4_482SZA, One50_pc_spec_05xCH4_482SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_48.2SZA_0.5xCH4_WPT_1rtol.vul')
One50_pc_V_01xCH4_482SZA, One50_pc_spec_01xCH4_482SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_48.2SZA_0.1xCH4_WPT_1rtol.vul')

One50_pc_V_10xCH4f_482SZA, One50_pc_spec_10xCH4f_482SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_48.2SZA_10xCH4flux_WPT_1rtol.vul')
One50_pc_V_5xCH4f_482SZA, One50_pc_spec_5xCH4f_482SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_48.2SZA_5xCH4flux_WPT_1rtol.vul')
One50_pc_V_1xCH4f_482SZA, One50_pc_spec_1xCH4f_482SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_48.2SZA_1xCH4flux_WPT_1rtol.vul')
One50_pc_V_05xCH4f_482SZA, One50_pc_spec_05xCH4f_482SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_48.2SZA_0.5xCH4flux_WPT_1rtol.vul')
One50_pc_V_01xCH4f_482SZA, One50_pc_spec_01xCH4f_482SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_48.2SZA_0.1xCH4flux_WPT_1rtol.vul')

PI_V_10xCH4f_482SZA, PI_pc_spec_10xCH4f_482SZA = Read_O3_Run(file_path='Earth_1e12s_48.2SZA_10xCH4flux_WPT_1rtol.vul')
PI_V_5xCH4f_482SZA, PI_pc_spec_5xCH4f_482SZA = Read_O3_Run(file_path='Earth_1e12s_48.2SZA_5xCH4flux_WPT_1rtol.vul')
PI_V_1xCH4f_482SZA, PI_pc_spec_1xCH4f_482SZA = Read_O3_Run(file_path='Earth_1e12s_48.2SZA_1xCH4flux_WPT_1rtol.vul')
PI_V_05xCH4f_482SZA, PI_pc_spec_05xCH4f_482SZA = Read_O3_Run(file_path='Earth_1e12s_48.2SZA_0.5xCH4flux_WPT_1rtol.vul')
PI_V_01xCH4f_482SZA, PI_pc_spec_01xCH4f_482SZA = Read_O3_Run(file_path='Earth_1e12s_48.2SZA_0.1xCH4flux_WPT_1rtol.vul')

PI_V_10xCH4_482SZA, PI_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_1e12s_48.2SZA_10xCH4_WPT_1rtol.vul')
PI_V_5xCH4_482SZA, PI_pc_spec_5xCH4_482SZA = Read_O3_Run(file_path='Earth_1e12s_48.2SZA_5xCH4_WPT_1rtol.vul')
PI_V_05xCH4_482SZA, PI_pc_spec_05xCH4_482SZA = Read_O3_Run(file_path='Earth_1e12s_48.2SZA_0.5xCH4_WPT_1rtol.vul')
PI_V_01xCH4_482SZA, PI_pc_spec_01xCH4_482SZA = Read_O3_Run(file_path='Earth_1e12s_48.2SZA_0.1xCH4_WPT_1rtol.vul')

Fifty_pc_V_10xCH4_482SZA, Fifty_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_48.2SZA_10xCH4_WPT_1rtol.vul')
Fifty_pc_V_5xCH4_482SZA, Fifty_pc_spec_5xCH4_482SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_48.2SZA_5xCH4_WPT_1rtol.vul')
Fifty_pc_V_05xCH4_482SZA, Fifty_pc_spec_05xCH4_482SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_48.2SZA_0.5xCH4_WPT_1rtol.vul')
Fifty_pc_V_01xCH4_482SZA, Fifty_pc_spec_01xCH4_482SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_48.2SZA_0.1xCH4_WPT_1rtol.vul')

Fifty_pc_V_10xCH4f_482SZA, Fifty_pc_spec_10xCH4f_482SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_48.2SZA_10xCH4flux_WPT_1rtol.vul')
Fifty_pc_V_5xCH4f_482SZA, Fifty_pc_spec_5xCH4f_482SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_48.2SZA_5xCH4flux_WPT_1rtol.vul')
Fifty_pc_V_1xCH4f_482SZA, Fifty_pc_spec_1xCH4f_482SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_48.2SZA_1xCH4flux_WPT_1rtol.vul')
Fifty_pc_V_05xCH4f_482SZA, Fifty_pc_spec_05xCH4f_482SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_48.2SZA_0.5xCH4flux_WPT_1rtol.vul')
Fifty_pc_V_01xCH4f_482SZA, Fifty_pc_spec_01xCH4f_482SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_48.2SZA_0.1xCH4flux_WPT_1rtol.vul')

Ten_pc_V_10xCH4_482SZA, Ten_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_48.2SZA_10xCH4_WPT_1rtol.vul')
Ten_pc_V_5xCH4_482SZA, Ten_pc_spec_5xCH4_482SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_48.2SZA_5xCH4_WPT_1rtol.vul')
Ten_pc_V_05xCH4_482SZA, Ten_pc_spec_05xCH4_482SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_48.2SZA_0.5xCH4_WPT_1rtol.vul')
Ten_pc_V_01xCH4_482SZA, Ten_pc_spec_01xCH4_482SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_48.2SZA_0.1xCH4_WPT_1rtol.vul')

Ten_pc_V_10xCH4f_482SZA, Ten_pc_spec_10xCH4f_482SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_48.2SZA_10xCH4flux_WPT_1rtol.vul')
Ten_pc_V_5xCH4f_482SZA, Ten_pc_spec_5xCH4f_482SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_48.2SZA_5xCH4flux_WPT_1rtol.vul')
Ten_pc_V_1xCH4f_482SZA, Ten_pc_spec_1xCH4f_482SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_48.2SZA_1xCH4flux_WPT_1rtol.vul')
Ten_pc_V_05xCH4f_482SZA, Ten_pc_spec_05xCH4f_482SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_48.2SZA_0.5xCH4flux_WPT_1rtol.vul')
Ten_pc_V_01xCH4f_482SZA, Ten_pc_spec_01xCH4f_482SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_48.2SZA_0.1xCH4flux_WPT_1rtol.vul')

Five_pc_V_10xCH4_482SZA, Five_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_48.2SZA_10xCH4_WPT_1rtol.vul')
Five_pc_V_5xCH4_482SZA, Five_pc_spec_5xCH4_482SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_48.2SZA_5xCH4_WPT_1rtol.vul')
Five_pc_V_05xCH4_482SZA, Five_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_48.2SZA_0.5xCH4_WPT_1rtol.vul')
Five_pc_V_01xCH4_482SZA, Five_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_48.2SZA_0.1xCH4_WPT_1rtol.vul')

Five_pc_V_10xCH4f_482SZA, Five_pc_spec_10xCH4f_482SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_48.2SZA_10xCH4flux_WPT_1rtol.vul')
Five_pc_V_5xCH4f_482SZA, Five_pc_spec_5xCH4f_482SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_48.2SZA_5xCH4flux_WPT_1rtol.vul')
Five_pc_V_1xCH4f_482SZA, Five_pc_spec_1xCH4f_482SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_48.2SZA_1xCH4flux_WPT_1rtol.vul')
Five_pc_V_05xCH4f_482SZA, Five_pc_spec_05xCH4f_482SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_48.2SZA_0.5xCH4flux_WPT_1rtol.vul')
Five_pc_V_01xCH4f_482SZA, Five_pc_spec_01xCH4f_482SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_48.2SZA_0.1xCH4flux_WPT_1rtol.vul')

One_pc_V_10xCH4_482SZA, One_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_48.2SZA_10xCH4_WPT_1rtol.vul')
One_pc_V_5xCH4_482SZA, One_pc_spec_5xCH4_482SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_48.2SZA_5xCH4_WPT_1rtol.vul')
One_pc_V_05xCH4_482SZA, One_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_48.2SZA_0.5xCH4_WPT_1rtol.vul')
One_pc_V_01xCH4_482SZA, One_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_48.2SZA_0.1xCH4_WPT_1rtol.vul')

One_pc_V_10xCH4f_482SZA, One_pc_spec_10xCH4f_482SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_48.2SZA_10xCH4flux_WPT_1rtol.vul')
One_pc_V_5xCH4f_482SZA, One_pc_spec_5xCH4f_482SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_48.2SZA_5xCH4flux_WPT_1rtol.vul')
One_pc_V_1xCH4f_482SZA, One_pc_spec_1xCH4f_482SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_48.2SZA_1xCH4flux_WPT_1rtol.vul')
One_pc_V_05xCH4f_482SZA, One_pc_spec_05xCH4f_482SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_48.2SZA_0.5xCH4flux_WPT_1rtol.vul')
One_pc_V_01xCH4f_482SZA, One_pc_spec_01xCH4f_482SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_48.2SZA_0.1xCH4flux_WPT_1rtol.vul')

Z5_pc_V_10xCH4_482SZA, Z5_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_10xCH4_WPT_1rtol.vul')
Z5_pc_V_5xCH4_482SZA, Z5_pc_spec_5xCH4_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_5xCH4_WPT_1rtol.vul')
Z5_pc_V_05xCH4_482SZA, Z5_pc_spec_05xCH4_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_0.5xCH4_WPT_1rtol.vul')
Z5_pc_V_01xCH4_482SZA, Z5_pc_spec_01xCH4_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_0.1xCH4_WPT_1rtol.vul')

Z5_pc_V_10xCH4f_482SZA, Z5_pc_spec_10xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_10xCH4flux_WPT_1rtol.vul')
Z5_pc_V_5xCH4f_482SZA, Z5_pc_spec_5xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_5xCH4flux_WPT_1rtol.vul')
Z5_pc_V_1xCH4f_482SZA, Z5_pc_spec_1xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_1xCH4flux_WPT_1rtol.vul')
Z5_pc_V_05xCH4f_482SZA, Z5_pc_spec_05xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_0.5xCH4flux_WPT_1rtol.vul')
Z5_pc_V_01xCH4f_482SZA, Z5_pc_spec_01xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_0.1xCH4flux_WPT_1rtol.vul')

Z1_pc_V_10xCH4_482SZA, Z1_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_10xCH4_WPT_1rtol.vul')
Z1_pc_V_5xCH4_482SZA, Z1_pc_spec_5xCH4_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_5xCH4_WPT_1rtol.vul')
Z1_pc_V_05xCH4_482SZA, Z1_pc_spec_05xCH4_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_0.5xCH4_WPT_1rtol.vul')
Z1_pc_V_01xCH4_482SZA, Z1_pc_spec_01xCH4_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_0.1xCH4_WPT_1rtol.vul')

Z1_pc_V_10xCH4f_482SZA, Z1_pc_spec_10xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_10xCH4flux_WPT_1rtol.vul')
Z1_pc_V_5xCH4f_482SZA, Z1_pc_spec_5xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_5xCH4flux_WPT_1rtol.vul')
Z1_pc_V_1xCH4f_482SZA, Z1_pc_spec_1xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_1xCH4flux_WPT_1rtol.vul')
Z1_pc_V_05xCH4f_482SZA, Z1_pc_spec_05xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_0.5xCH4flux_WPT_1rtol.vul')
Z1_pc_V_01xCH4f_482SZA, Z1_pc_spec_01xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_0.1xCH4flux_WPT_1rtol.vul')

#%% Cross sections and checks
# Cross sections and the wavelength grid are in photos.pdat (read by READPHOTO),
# not in OUTPUT_PLOT.dat (that file is only Z, P, photolysis rates J, etc.).

KASTING_PHOTOS_PDAT = (
    '/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Kasting_1D_model/'
    'Old sims/Kasting model/DATA/photos.pdat'
)

data_blocks = read_kasting_file(KASTING_PHOTOS_PDAT)
df_flux = data_blocks['flux_block']
df_cross = data_blocks['cross_sections']

df_flux[['Wave_Min', 'Wave_Max']] = (
    df_flux['Range_A'].str.split('-', expand=True).astype(float)
)
# Wavelength at bin centre [nm] (photos.pdat uses Angstrom)
wav_nm_flux = (df_flux['Wave_Min'] + df_flux['Wave_Max']) / 20.0
n_cross = len(df_cross)
wav_nm_cross = wav_nm_flux.iloc[:n_cross].to_numpy()

print("--- Wavelength & flux (first 10 bins) ---")
print(df_flux[['Int', 'Wave_Min', 'Wave_Max', 'Flux']].head(10))

print("\n--- Cross sections (first 10 bins) ---")
print(df_cross[['Int', 'O2', 'H2O', 'CO2', 'N2O']].head(10))

# O2 in the SR bands: use Ozone1 + Ozone2 from the flux table (not the O2 column,
# which is zero there because O2 uses correlated-k, not a single sigma per bin).
sigma_o2_sr = df_flux['Ozone1'].iloc[:n_cross] + df_flux['Ozone2'].iloc[:n_cross]

species_plots = [
    ('H2O', df_cross['H2O'], None),
    ('N2O', df_cross['N2O'], None),
    ('CO2', df_cross['CO2'], None),
    ('HO2', df_cross['HO2'], None),
    ('O2 (Ozone1+Ozone2, SR grid)', sigma_o2, None),
]

for title, sigma, _ in species_plots:
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(wav_nm_cross, sigma, color='teal', lw=2, label='Kasting photos.pdat')
    # Optional overlay with VULCAN (requires PI_V_482SZA loaded above)
    if 'PI_V_482SZA' in globals():
        vul_bins = PI_V_482SZA['variable']['bins']
        cj = PI_V_482SZA['variable']['cross_J']
        if title.startswith('H2O'):
            ax.plot(vul_bins, cj['H2O', 1], color='m', alpha=0.7, label='VULCAN branch 1')
            ax.plot(vul_bins, cj['H2O', 2], color='m', ls='--', alpha=0.7, label='VULCAN branch 2')
        elif title.startswith('N2O'):
            ax.plot(vul_bins, cj['N2O', 1], color='m', alpha=0.7, label='VULCAN')
        elif title.startswith('O2'):
            ax.plot(
                vul_bins, cj['O2', 1] + cj['O2', 2],
                color='black', alpha=0.7, label='VULCAN O2 (branches 1+2)',
            )
    ax.set_yscale('log')
    ax.set_xlim(1, 300)
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('Cross section [cm²]')
    ax.set_title(title)
    ax.legend(loc='best')
    thick_axes()
    plt.tight_layout()
#%% Read in Kasting simulations

path = '/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Kasting_1D_model/'
K_150pc_J = pd.read_csv(
    path+"150pc/SZA_48.2/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_150pc = K_150pc_J['PRESS']/1000
K_150pc_J_45SZA = pd.read_csv(
    path+"150pc/SZA_45/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_150pc_45SZA = K_150pc_J_45SZA['PRESS']/1000
K_150pc_J_60SZA = pd.read_csv(
    path+"150pc/SZA_60/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_150pc_60SZA = K_150pc_J_60SZA['PRESS']/1000
z = K_150pc_J['Z']/1000



K_100pc_J_8G = pd.read_csv(
    path+"100pc/8point/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_100pc_8G = K_100pc_J_8G['PRESS']/1000

K_100pc_J = pd.read_csv(
    path+"100pc/SZA_48.2/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_100pc = K_100pc_J['PRESS']/1000
K_100pc_J_45SZA = pd.read_csv(
    path+"100pc/SZA_45/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_100pc_45SZA = K_100pc_J_45SZA['PRESS']/1000
K_100pc_J_60SZA = pd.read_csv(
    path+"100pc/SZA_60/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_100pc_60SZA = K_100pc_J_60SZA['PRESS']/1000
z = K_100pc_J['Z']/1000

K_10pc_J_8G = pd.read_csv(
    path+"10pc/8point/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_10pc_8G = K_10pc_J_8G['PRESS']/1000

K_10pc_J = pd.read_csv(
    path+"10pc/SZA_48.2/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_10pc = K_10pc_J['PRESS']/1000
K_10pc_J_45SZA = pd.read_csv(
    path+"10pc/SZA_45/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_10pc_45SZA = K_10pc_J_45SZA['PRESS']/1000
K_10pc_J_60SZA = pd.read_csv(
    path+"10pc/SZA_60/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_10pc_60SZA = K_10pc_J_60SZA['PRESS']/1000
z = K_10pc_J['Z']/1000


K_1pc_J_8G = pd.read_csv(
    path+"1pc/8point/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_1pc_8G = K_1pc_J_8G['PRESS']/1000

K_1pc_J = pd.read_csv(
    path+"1pc/SZA_48.2/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_1pc = K_1pc_J['PRESS']/1000
K_1pc_J_45SZA = pd.read_csv(
    path+"1pc/SZA_45/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_1pc_45SZA = K_1pc_J_45SZA['PRESS']/1000
K_1pc_J_60SZA = pd.read_csv(
    path+"1pc/SZA_60/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_1pc_60SZA = K_1pc_J_60SZA['PRESS']/1000
z = K_1pc_J['Z']/1000

K_01pc_J_8G = pd.read_csv(
    path+"0.1pc/8point/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_01pc_8G = K_01pc_J_8G['PRESS']/1000

K_01pc_J = pd.read_csv(
    path+"0.1pc/SZA_48.2/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_01pc = K_1pc_J['PRESS']/1000
K_01pc_J_45SZA = pd.read_csv(
    path+"0.1pc/SZA_45/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_01pc_45SZA = K_1pc_J_45SZA['PRESS']/1000
K_01pc_J_60SZA = pd.read_csv(
    path+"0.1pc/SZA_60/OUTPUT_PLOT.dat", delim_whitespace=True, engine="python")
p_01pc_60SZA = K_1pc_J_60SZA['PRESS']/1000
z = K_1pc_J['Z']/1000

#%% J rates plot












#%% Read in Atmos files


import pandas as pd


'''
Now read in horribly formatted photolysis data
'''


file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/100pc/"
Atmos_100pc = pd.read_csv(file_path+'SZA_48.2/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_100pc_01xCH4 = pd.read_csv(file_path+'../Methane_Perturbations/0.1x/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_100pc_05xCH4 = pd.read_csv(file_path+'../Methane_Perturbations/0.5x/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_100pc_1xCH4 = pd.read_csv(file_path+'../Methane_Perturbations/1x/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_100pc_2xCH4 = pd.read_csv(file_path+'../Methane_Perturbations/2x/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_100pc_3xCH4 = pd.read_csv(file_path+'../Methane_Perturbations/3x/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_100pc_4xCH4 = pd.read_csv(file_path+'../Methane_Perturbations/4x/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_100pc_5xCH4 = pd.read_csv(file_path+'../Methane_Perturbations/5x/PTZ_mixingratios_out.dist', delim_whitespace=True)

Atmos_100pc_45SZA = pd.read_csv(file_path+'SZA_45/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_100pc_60SZA = pd.read_csv(file_path+'SZA_60/PTZ_mixingratios_out.dist', delim_whitespace=True)
file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/10pc/"
Atmos_10pc = pd.read_csv(file_path+'SZA_48.2/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_10pc_45SZA = pd.read_csv(file_path+'SZA_45/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_10pc_60SZA = pd.read_csv(file_path+'SZA_60/PTZ_mixingratios_out.dist', delim_whitespace=True)
file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/1pc/"
Atmos_1pc = pd.read_csv(file_path+'SZA_48.2/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_1pc_45SZA = pd.read_csv(file_path+'SZA_45/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_1pc_60SZA = pd.read_csv(file_path+'SZA_60/PTZ_mixingratios_out.dist', delim_whitespace=True)
file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/0.1pc/"
Atmos_01pc = pd.read_csv(file_path+'SZA_48.2/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_01pc_45SZA = pd.read_csv(file_path+'SZA_45/PTZ_mixingratios_out.dist', delim_whitespace=True)
Atmos_01pc_60SZA = pd.read_csv(file_path+'SZA_60/PTZ_mixingratios_out.dist', delim_whitespace=True)

file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/100pc/SZA_48.2/out.O2prates"
Atmos_XS_100pc_Z, Atmos_XS_100pc_JO2 = Atmos_photo(file_path)
file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/10pc/SZA_48.2/out.O2prates"
Atmos_XS_10pc_Z, Atmos_XS_10pc_JO2 = Atmos_photo(file_path)
file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/1pc/SZA_48.2/out.O2prates"
Atmos_XS_1pc_Z, Atmos_XS_1pc_JO2= Atmos_photo(file_path)
file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/0.1pc/SZA_48.2/out.O2prates"
Atmos_XS_01pc_Z, Atmos_XS_01pc_JO2= Atmos_photo(file_path)

file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/100pc/SZA_48.2/out.O2prates"
Atmos_XS_100pc_Z, Atmos_XS_100pc_JO2 = Atmos_photo(file_path)
file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/10pc/SZA_48.2/out.O2prates"
Atmos_XS_10pc_Z, Atmos_XS_10pc_JO2 = Atmos_photo(file_path)
file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/1pc/SZA_48.2/out.O2prates"
Atmos_XS_1pc_Z, Atmos_XS_1pc_JO2= Atmos_photo(file_path)
file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/0.1pc/SZA_48.2/out.O2prates"
Atmos_XS_01pc_Z, Atmos_XS_01pc_JO2= Atmos_photo(file_path)

#%% Read in Phtoochem files

# Photochem Ozone column 
path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/150pc/"
Photo_150pc = pd.read_csv(path+"Earth_150pc_48.2.txt", delim_whitespace=True)
'''
Photo_150pc_01CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_150pc_48.2_0.1x_methane.txt", delim_whitespace=True)
Photo_150pc_05CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_0.5x_methane.txt", delim_whitespace=True)
Photo_150pc_1CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_1x_methane.txt", delim_whitespace=True)
Photo_150pc_2CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_2x_methane.txt", delim_whitespace=True)
Photo_150pc_3CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_3x_methane.txt", delim_whitespace=True)
Photo_150pc_4CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_4x_methane.txt", delim_whitespace=True)
Photo_150pc_5CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_5x_methane.txt", delim_whitespace=True)
'''
path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/100pc/"
Photo_PI = pd.read_csv(path+"Earth_100pc_48.2.txt", delim_whitespace=True)
Photo_PI_60SZA = pd.read_csv(path+"Earth_100pc_60.txt", delim_whitespace=True)
Photo_PI_45SZA = pd.read_csv(path+"Earth_100pc_45.txt", delim_whitespace=True)
Photo_PI_01CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_0.1x_methane.txt", delim_whitespace=True)
Photo_PI_05CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_0.5x_methane.txt", delim_whitespace=True)
Photo_PI_1CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_1x_methane.txt", delim_whitespace=True)
Photo_PI_2CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_2x_methane.txt", delim_whitespace=True)
Photo_PI_3CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_3x_methane.txt", delim_whitespace=True)
Photo_PI_4CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_4x_methane.txt", delim_whitespace=True)
Photo_PI_5CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_5x_methane.txt", delim_whitespace=True)


path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/10pc/"
Photo_10pc = pd.read_csv(path+'Earth_10pc_48.2.txt', delim_whitespace=True)
Photo_10pc_60SZA = pd.read_csv(path+'Earth_10pc_60.txt', delim_whitespace=True)
Photo_10pc_45SZA = pd.read_csv(path+'Earth_10pc_45.txt', delim_whitespace=True)
Photo_10pc_01CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_10pc_48.2_0.1x_methane.txt", delim_whitespace=True)
Photo_10pc_05CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_10pc_48.2_0.5x_methane.txt", delim_whitespace=True)
Photo_10pc_1CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_10pc_48.2_1x_methane.txt", delim_whitespace=True)
Photo_10pc_2CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_10pc_48.2_2x_methane.txt", delim_whitespace=True)
Photo_10pc_3CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_10pc_48.2_3x_methane.txt", delim_whitespace=True)
Photo_10pc_4CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_10pc_48.2_4x_methane.txt", delim_whitespace=True)
Photo_10pc_5CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_10pc_48.2_5x_methane.txt", delim_whitespace=True)

path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/1pc/"
Photo_1pc = pd.read_csv(path+'Earth_1pc_48.2.txt', delim_whitespace=True)
Photo_1pc_60SZA = pd.read_csv(path+'Earth_1pc_60.txt', delim_whitespace=True)
Photo_1pc_45SZA = pd.read_csv(path+'Earth_1pc_45.txt', delim_whitespace=True)
Photo_1pc_01CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_1pc_48.2_0.1x_methane.txt", delim_whitespace=True)
Photo_1pc_05CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_1pc_48.2_0.5x_methane.txt", delim_whitespace=True)
Photo_1pc_1CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_1pc_48.2_1x_methane.txt", delim_whitespace=True)
Photo_1pc_2CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_1pc_48.2_2x_methane.txt", delim_whitespace=True)
Photo_1pc_3CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_1pc_48.2_3x_methane.txt", delim_whitespace=True)
Photo_1pc_4CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_1pc_48.2_4x_methane.txt", delim_whitespace=True)
Photo_1pc_5CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_1pc_48.2_5x_methane.txt", delim_whitespace=True)

path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/0.1pc/"
Photo_01pc = pd.read_csv(path+'Earth_0.1pc_48.2.txt', delim_whitespace=True)
Photo_01pc_60SZA = pd.read_csv(path+'Earth_0.1pc_60.txt', delim_whitespace=True)
Photo_01pc_45SZA = pd.read_csv(path+'Earth_0.1pc_45.txt', delim_whitespace=True)
Photo_01pc_01CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_0.1pc_48.2_0.1x_methane.txt", delim_whitespace=True)
Photo_01pc_05CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_0.1pc_48.2_0.5x_methane.txt", delim_whitespace=True)
Photo_01pc_1CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_0.1pc_48.2_1x_methane.txt", delim_whitespace=True)
Photo_01pc_2CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_0.1pc_48.2_2x_methane.txt", delim_whitespace=True)
Photo_01pc_3CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_0.1pc_48.2_3x_methane.txt", delim_whitespace=True)
Photo_01pc_4CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_0.1pc_48.2_4x_methane.txt", delim_whitespace=True)
Photo_01pc_5CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_0.1pc_48.2_5x_methane.txt", delim_whitespace=True)

#%% color choices

color_Kasting = 'teal'
color_VULCAN = 'magenta'
color_Photochem = 'blue'
color_Atmos = 'darkorange'
color_WACCM = 'k'


blues = plt.colormaps.get_cmap('Blues')
darkblue = blues(0.8)
blue = blues(0.6)
lightblue = blues(0.4)

greens = plt.colormaps.get_cmap('Greens')
darkgreen = greens(0.8)
green = greens(0.6)
lightgreen = greens(0.4)

purples = plt.colormaps.get_cmap('Purples')
darkpurple = purples(1.0)
purple = purples(0.75)
lightpurple = purples(0.4)

PI_color = 'black'; 
One_50_color = 'grey'
Fifty_pc_color = darkgreen; Ten_pc_color = green; Five_pc_color = lightgreen
One_pc_color = darkblue; Zero5_pc_color = blue; Zero1_pc_color = lightblue; 
verylightblue = blues(0.2)
YS_color = 'm'; 
#%% import file for pressure and gaussian weights


File =  "/Users/gregcooke/H_escape/b.e21.BWma1850.f19_g17.PC_b.SSPO.016.cam.h0.0320-0320.nc" #file name

Gw_file = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

gw = Gw_file.gw.values
Pressure = Gw_file.lev.values




#%%


#%%

k_B = 1.380649e-23  # J/K

path = "/Users/gregcooke/H_escape/Earth_PI_SRB_scat.cam.h1.0001-01-01-01800.nc"
ds = xr.open_dataset(path, decode_times=False)
    
#%% SRB

SRB = xr.open_dataset('/Users/gregcooke/H_escape/Earth_PI_SRB_no_scat.cam.h1.0001-01-01-01800.nc',decode_times=False) #open the file and decode time as false
SRB_scat = xr.open_dataset('/Users/gregcooke/H_escape/Earth_PI_SRB_scat.cam.h1.0001-01-01-01800.nc',decode_times=False) #open the file and decode time as false
SRB_scat2 = xr.open_dataset('/Users/gregcooke/H_escape/Earth_PI_SRB_test_scat.cam.h1.0001-01-01-01800.nc',decode_times=False) #open the file and decode time as false

print(LWAV(SRB.jo2_a+SRB.jo2_b))
print(LWAV(SRB_scat.jo2_a+SRB_scat.jo2_b))


plt.figure()
plt.plot(LWAV(SRB.jo2_a+SRB.jo2_b), SRB.lev)
'''
plt.plot(LWAV(SRB.jo2_a), SRB.lev)
plt.plot(LWAV(SRB.jo2_b), SRB.lev)
'''
plt.plot(LWAV(SRB_scat.jo2_a+SRB_scat.jo2_b), SRB_scat.lev)
plt.plot(LWAV(SRB_scat2.jo2_a+SRB_scat2.jo2_b), SRB_scat2.lev)
'''
plt.plot(LWAV(SRB_scat.jo2_a), SRB_scat.lev)
plt.plot(LWAV(SRB_scat.jo2_b), SRB_scat.lev)
'''
plt.ylim(1e3, 1e-6); plt.yscale('log')
#plt.xlim(0,2)
plt.xscale('log')
plt.xlim(1e-13, 1e-7)

#%% define ozone column calculation function


#%% K stars
path = '/Users/gregcooke/H_escape/'

File = path+"b.e21.BWma1850.f19_g17.PI.no_obliq.001.cam.h0.0021.nc"  
Pre_0obq = xr.open_dataset(File,decode_times=False) #open the file and decode time as false  

File = path+"b.e21.BWma1850.f19_g17.10pc_o2.no_obliq.001.cam.h0.0044.nc"  
Ten_0obq = xr.open_dataset(File,decode_times=False) #open the file and decode time as false  

File = path+"b.e21.BWma1850.f19_g17.1pc_o2.no_obliq.001.cam.h0.0061.nc"  
One_0obq = xr.open_dataset(File,decode_times=False) #open the file and decode time as false  

File = path+"b.e21.BWma1850.f19_g17.0.1pc_o2.no_obliq.001.cam.h0.0066.nc"  
Zero1_0obq = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

#%% Read in low oxygen cases
path = '/Users/gregcooke/H_escape/'
File =  path+"b.e21.BWma1850.f19_g17.baseline.cam.h2.0014-0017.zm.nc" #file name
Pre_pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 150% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.150pc_o2.002.cam.h2.0034-0037.zm.nc" #file name
One50_pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 50% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.50pc_o2.002.cam.h2.0030-0033.zm.nc" #file name
Fifty_pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 10% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.10pc_o2.001.cam.h2.0034-0037.zm.nc" #file name
Ten_pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
#Ten_pc_h2 = Ten_pc_h2.mean(dim = 'time')

# 5% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.5pc_o2.002.cam.h2.0032-0035.zm.nc" #file name
Five_pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 1% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.1pc_o2.cam.h2.0028-0031.zm.nc" #file name
One_pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 0.5% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.0.5pc_o2.001.cam.h2.0046-0049.zm.nc" #file name
Zero5_pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 0.1% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.0.1pc_o2.001.cam.h2.0080-0083.zm.nc" #file name
Zero1_pc_h2 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
 

File = path+"b.e21.BWma1850.f19_g17.baseline.cam.h0.0014-0017.nc"
PI = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

File =  "/Users/gregcooke/CESM_data/Earth_baseline_PI.cam.h0.0001.nc" #file name
Pre_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# Pre-industrial baseliine run (PI; 100% PAL of O2)
File =  "/Users/gregcooke/H_escape/b.e21.BWma1850.f19_g17.baseline.cam.h0.0013-0017.nc" #file name
Pre_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 50% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.50pc_o2.002.cam.h0.0031-0034.nc" #file name
Fifty_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 10% PAL of O2
File = "/Users/gregcooke/H_escape/b.e21.BWma1850.f19_g17.10pc_o2.001.cam.h0.0034-0037.nc" #file name
Ten_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 5% PAL of O2
File =  path+"b.e21.BWma1850.f19_g17.5pc_o2.002.cam.h0.0032-0035.nc" #file name
Five_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 1% PAL of O2
File = "/Users/gregcooke/H_escape/b.e21.BWma1850.f19_g17.1pc_o2.cam.h0.0028-031.nc" #file name
One_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 0.5% PAL of O2
File = path+"b.e21.BWma1850.f19_g17.0.5pc_o2.001.cam.h0.0046-0049.nc" #file name
Zero5_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

# 0.1% PAL of O2
File =  "/Users/gregcooke/H_escape/b.e21.BWma1850.f19_g17.0.1pc_o2.001.cam.h0.0080-0083.nc" #file name
Zero1_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

File =  path+"b.e21.BWma1850.f19_g17.150pc_o2.002.cam.h0.0034-0037.nc" #file name
One50_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

#%% P-T profile plot
alpha = 0.2
plt.figure(figsize = (14,7))

plt.plot(LWAV(Pre_pc_h0.T), Pre_pc_h0.lev, lw = lw, ls = '-', color = 'black', label = '100% PAL')
T_min = np.empty(0); T_max = np.empty(0)
for i in range(len(LWAV(Pre_pc_h0.T))):
    T_min = np.append(T_min, Pre_pc_h0.T.isel(lev = i).min())
    T_max = np.append(T_max, Pre_pc_h0.T.isel(lev = i).max())
plt.fill_betweenx(Pre_pc_h0.lev, T_min, T_max, color = 'black', alpha = alpha)

plt.plot(LWAV(Ten_pc_h0.T), Pre_pc_h0.lev, lw = lw, ls = '-', color = Ten_pc_color, label = '10% PAL')
T_min = np.empty(0); T_max = np.empty(0)
for i in range(len(LWAV(Ten_pc_h0.T))):
    T_min = np.append(T_min, Ten_pc_h0.T.isel(lev = i).min())
    T_max = np.append(T_max, Ten_pc_h0.T.isel(lev = i).max())
plt.fill_betweenx(Ten_pc_h0.lev, T_min, T_max, color = Ten_pc_color, alpha =  alpha)

plt.plot(LWAV(One_pc_h0.T), Pre_pc_h0.lev, lw = lw, ls = '-', color = One_pc_color, label = '1% PAL')
T_min = np.empty(0); T_max = np.empty(0)
for i in range(len(LWAV(Ten_pc_h0.T))):
    T_min = np.append(T_min, One_pc_h0.T.isel(lev = i).min())
    T_max = np.append(T_max, One_pc_h0.T.isel(lev = i).max())
plt.fill_betweenx(One_pc_h0.lev, T_min, T_max, color = One_pc_color, alpha =  alpha)

plt.plot(LWAV(Zero1_pc_h0.T), Pre_pc_h0.lev, lw = lw, ls = '-', color = Zero1_pc_color, label = '0.1% PAL')
T_min = np.empty(0); T_max = np.empty(0)
for i in range(len(LWAV(Ten_pc_h0.T))):
    T_min = np.append(T_min, Zero1_pc_h0.T.isel(lev = i).min())
    T_max = np.append(T_max, Zero1_pc_h0.T.isel(lev = i).max())
plt.fill_betweenx(One_pc_h0.lev, T_min, T_max, color = Zero1_pc_color, alpha =  alpha)

plt.ylim(1e3, 1e-4); plt.yscale('log')
plt.xlabel('Temperature [K]', fontsize = 15, weight = 'bold')
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.xlim(145, 300)
thick_axes(top = True)
plt.legend(loc = (0.71, 0.06), fontsize = 15, 
           frameon = False, ncol = 1)
plt.savefig('/Users/gregcooke/python_output/P-T_profiles', dpi = 400, bbox_inches = 'tight')
#%% Now calculate ozone columns

One50_col = O3_col(One50_pc_h0.mean(dim='time'), time = False, O2_mr = 0.21*1.5, lon = True, g = 9.81)
Pre_col = O3_col(Pre_pc_h0.mean(dim='time'), time = False, O2_mr = 0.21, lon = True, g = 9.81)
Fifty_col = O3_col(Fifty_pc_h0.mean(dim='time'), time = False, O2_mr = 0.21/2, lon = True, g = 9.81)
Ten_col = O3_col(Ten_pc_h0.mean(dim='time'), time = False, O2_mr = 0.21/10, lon = True, g = 9.81)
Five_col = O3_col(Five_pc_h0.mean(dim='time'), time = False, O2_mr = 0.21/20, lon = True, g = 9.81)
One_col = O3_col(One_pc_h0.mean(dim='time'), time = False, O2_mr = 0.21/100, lon = True, g = 9.81)
Zero5_col = O3_col(Zero5_pc_h0.mean(dim='time'), time = False, O2_mr = 0.21/200, lon = True, g = 9.81)
Zero1_col = O3_col(Zero1_pc_h0.mean(dim='time'), time = False, O2_mr = 0.21/1000, lon = True, g = 9.81)
Pre_0obq_col = O3_col(Pre_0obq.mean(dim='time'), time = False, O2_mr = 0.21, lon = True, g = 9.81)
Ten_0obq_col = O3_col(Ten_0obq.mean(dim='time'), time = False, O2_mr = 0.21/10, lon = True, g = 9.81)
One_0obq_col = O3_col(One_0obq.mean(dim='time'), time = False, O2_mr = 0.21/100, lon = True, g = 9.81)
Zero1_0obq_col = O3_col(Zero1_0obq.mean(dim='time'), time = False, O2_mr = 0.21/1000, lon = True, g = 9.81)

#%% Winds - look at polar vortex if you can

plt.figure()
Pre_pc_h0.U.isel(lev = 44).plot(cmap = 'RdBu_r', vmax = 40, vmin = -40)

plt.figure()
Zero1_pc_h0.U.isel(lev = 44).plot(cmap = 'RdBu_r', vmax = 40, vmin = -40)

plt.figure()
Pre_pc_h0.V.isel(lev = 44).plot(cmap = 'RdBu_r', vmax = 5, vmin = -5)

plt.figure()
Zero1_pc_h0.V.isel(lev = 44).plot(cmap = 'RdBu_r', vmax = 5, vmin = -5)

#%% Tropospheric polar vortex

lev = 500
geopotential_height = Pre_pc_h0['Z3'].sel(lev=lev, method='nearest')
location = np.arange(4000, 6000, 50); contour = 20000

# Select a specific time step if you have a time dimension
geopotential_slice = geopotential_height.isel(time=0)

# Create a figure and axes with a polar projection
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())

# Set the extent of the plot (e.g., poleward of 30°N)
ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
ax.coastlines()
ax.gridlines()

# Plot geopotential height as filled contours
geopotential_slice.plot.contourf(
    ax=ax, transform=ccrs.PlateCarree(), levels=location,
    cmap='viridis', cbar_kwargs={'label': 'Geopotential Height (m)'}
)

# Plot the 50 hPa geopotential height contour that defines the edge of the vortex
# The exact value (e.g., 19000 m) is often chosen based on climatology or a specific definition.
# A common choice is to plot the 19,000 m contour for the Northern Hemisphere winter.
geopotential_slice.plot.contour(
    ax=ax, transform=ccrs.PlateCarree(), levels=[contour], colors='white', linewidths=2
)

plt.title('Polar Vortex at ' + str(lev) + ' hPa')
plt.show()

geopotential_height = Zero1_pc_h0['Z3'].sel(lev=lev, method='nearest')
location = np.arange(4000, 6000, 50); contour = 20000

# Select a specific time step if you have a time dimension
geopotential_slice = geopotential_height.isel(time=0)

# Create a figure and axes with a polar projection
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())

# Set the extent of the plot (e.g., poleward of 30°N)
ax.set_extent([-180, 180, 30, 90], crs=ccrs.PlateCarree())
ax.coastlines()
ax.gridlines()

# Plot geopotential height as filled contours
geopotential_slice.plot.contourf(
    ax=ax, transform=ccrs.PlateCarree(), levels=location,
    cmap='viridis', cbar_kwargs={'label': 'Geopotential Height (m)'}
)

# Plot the 50 hPa geopotential height contour that defines the edge of the vortex
# The exact value (e.g., 19000 m) is often chosen based on climatology or a specific definition.
# A common choice is to plot the 19,000 m contour for the Northern Hemisphere winter.
geopotential_slice.plot.contour(
    ax=ax, transform=ccrs.PlateCarree(), levels=[contour], colors='white', linewidths=2
)

plt.title('Polar Vortex at ' + str(lev) + ' hPa')
plt.show()

#%% Cross sections
file_path = "/Users/gregcooke/Aoshuang_2_CorrK4_clean/DATA/photos.pdat"

with open(file_path, "r") as f:
    lines = f.readlines()

# Print lines 50–100 (Python is 0-indexed, so line 50 = index 49)
for line in lines[0:110]:
    print(line.rstrip())
        
#%%



file_path = "/Users/gregcooke/Aoshuang_2_CorrK4_clean/DATA/photos.pdat"

int_list = []
flux_list = []
ozone1_list = []
ozone2_list = []
ch3cho_list = []
ch3coch3_list = []
range_min_list = []
range_max_list = []

skipped_lines = []

with open(file_path, "r") as f:
    lines = f.readlines()

# Skip the first two header lines
for i, line in enumerate(lines[2:], start=3):
    line = line.strip()
    if not line:
        continue
    
    parts = line.split()
    
    # Skip lines that don't have enough columns
    if len(parts) < 7:
        skipped_lines.append((i, line))
        continue
    
    try:
        int_val = int(parts[0])
        range_min, range_max = map(float, parts[1].split('-'))
        flux_val = float(parts[2])
        ozone1_val = float(parts[3])
        ozone2_val = float(parts[4])
        ch3cho_val = float(parts[5])
        ch3coch3_val = float(parts[6])
        
        int_list.append(int_val)
        range_min_list.append(range_min)
        range_max_list.append(range_max)
        flux_list.append(flux_val)
        ozone1_list.append(ozone1_val)
        ozone2_list.append(ozone2_val)
        ch3cho_list.append(ch3cho_val)
        ch3coch3_list.append(ch3coch3_val)
        
    except ValueError:
        skipped_lines.append((i, line))

# Convert to NumPy arrays
int_arr      = np.array(int_list)
flux_arr     = np.array(flux_list)
ozone1_arr   = np.array(ozone1_list)
ozone2_arr   = np.array(ozone2_list)
ch3cho_arr   = np.array(ch3cho_list)
ch3coch3_arr = np.array(ch3coch3_list)
range_min    = np.array(range_min_list)
range_max    = np.array(range_max_list)

print(f"Parsed {len(int_arr)} valid lines, skipped {len(skipped_lines)} lines.")
# Optional: print skipped lines
for ln_num, ln_text in skipped_lines[:5]:
    print(f"Skipped line {ln_num}: {ln_text}")

wav = (range_min+range_max)/2
plt.figure(figsize = (9,5))
plt.plot(wav/10, ozone1_arr)
plt.plot(wav/10, ozone2_arr)
thick_axes()
plt.yscale('log')
plt.ylabel('Cross section')
plt.xlabel('Wavelength [nm]')

#%% Kasting photo files


file_path = "/Users/gregcooke/Aoshuang_2_CorrK4_clean/DATA/photos.pdat"

int_list = []
o2_list = []
h2o_list = []
co2_list = []
ho2_list = []
n2o_list = []
hcl_list = []
ccl2f2_list = []
chclf2_list = []
ch3cl_list = []
ch3ccl3_list = []

parsing_table = False
last_int = 0

with open(file_path, "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        
        # Detect the start of the "Assorted cross sections" table
        if line.startswith("Int #") and "O2" in line:
            parsing_table = True
            last_int = 0
            continue
        if parsing_table and line.startswith("="):
            continue  # skip separator line
        
        if parsing_table:
            parts = line.split()
            if len(parts) < 11:
                continue  # skip malformed lines
            
            try:
                int_val = int(parts[0])
                # Stop if the Int # resets (new table starts)
                if int_val <= last_int:
                    break
                last_int = int_val
                
                int_list.append(int_val)
                o2_list.append(float(parts[1]))
                h2o_list.append(float(parts[2]))
                co2_list.append(float(parts[3]))
                ho2_list.append(float(parts[4]))
                n2o_list.append(float(parts[5]))
                hcl_list.append(float(parts[6]))
                ccl2f2_list.append(float(parts[7]))
                chclf2_list.append(float(parts[8]))
                ch3cl_list.append(float(parts[9]))
                ch3ccl3_list.append(float(parts[10]))
            except ValueError:
                continue

# Convert to NumPy arrays
int_arr     = np.array(int_list)
o2_arr      = np.array(o2_list)
h2o_arr     = np.array(h2o_list)
co2_arr     = np.array(co2_list)
ho2_arr     = np.array(ho2_list)
n2o_arr     = np.array(n2o_list)
hcl_arr     = np.array(hcl_list)
ccl2f2_arr  = np.array(ccl2f2_list)
chclf2_arr  = np.array(chclf2_list)
ch3cl_arr   = np.array(ch3cl_list)
ch3ccl3_arr = np.array(ch3ccl3_list)

print("CO2 array length:", len(co2_arr))
print("First 5 CO2 values:", co2_arr[:5])

#%% Solar file from different codes comparison 

file_path = "/Users/gregcooke/Aoshuang_2_CorrK4_clean/DATA/faruv_sun_3d.pdat"

L_list = []
SFX_list = []

with open(file_path, "r") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("=") or line.lower().startswith("note"):
            continue  # skip separators and note lines
        parts = line.split()
        if len(parts) < 2:
            continue  # skip lines without enough data
        try:
            L_list.append(float(parts[0]))
            SFX_list.append(float(parts[1]))
        except ValueError:
            continue  # skip lines that don't contain numbers

# Convert to NumPy arrays
FarUV_3D_Wav_arr = np.array(L_list)
FarUV_3D_SFX_arr = np.array(SFX_list)

file_path = "/Users/gregcooke/Aoshuang_2_CorrK4_clean/DATA/faruv_sun.pdat"

L_list = []
SFX_list = []

with open(file_path, "r") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("=") or line.lower().startswith("note"):
            continue  # skip separators and note lines
        parts = line.split()
        if len(parts) < 2:
            continue  # skip lines without enough data
        try:
            L_list.append(float(parts[0]))
            SFX_list.append(float(parts[1]))
        except ValueError:
            continue  # skip lines that don't contain numbers

# Convert to NumPy arrays
FarUV_Wav_arr = np.array(L_list)
FarUV_SFX_arr = np.array(SFX_list)



FarUV_SFX_3D_arr_corrected = convert_photons_to_energy(FarUV_3D_Wav_arr/10, FarUV_3D_SFX_arr)
FarUV_SFX_arr_corrected = convert_photons_to_energy(FarUV_Wav_arr/10, FarUV_SFX_arr)
Flux_arr_corrected = convert_photons_to_energy(wav/10, flux_arr)

VULCAN_Sun = pd.read_csv('/Users/gregcooke/VULCAN/atm/stellar_flux/VPL_solar.txt', 
                         skiprows = 1,  delim_whitespace = True,
                         names = ['Wav', 'Flux'])

VULCAN_Sun_EUV = pd.read_csv('/Users/gregcooke/VULCAN/atm/stellar_flux/Gueymard_solar.txt', 
                         skiprows = 1,  delim_whitespace = True,
                         names = ['Wav', 'Flux'])
VULCAN_Sun['Flux'] = VULCAN_Sun['Flux']*1e4*1e-7
VULCAN_Sun_EUV['Flux'] = VULCAN_Sun_EUV['Flux']*1e4*1e-7


solar_dir = '/Users/gregcooke/stellar_files/'
solar_file = 'SolarForcingCMIP6piControl_c160921.nc'#read in file
ds = xr.open_dataset(solar_dir+solar_file)#attach file to dataset
ssi = ds['ssi'].isel(time=0) #define dataset from file

WACCM_Flux = ssi.values
WACCM_wavelength = ssi.wavelength.values

#normalise all spectra

WACCM_TSI = np.trapz(WACCM_Flux[165:760], WACCM_wavelength[165:760])/1000
Kasting_TSI = np.trapz(Flux_arr_corrected, wav/10)/1000
Kasting_factor = Kasting_TSI / WACCM_TSI
VULCAN_TSI = np.trapz(VULCAN_Sun['Flux'][88:175], VULCAN_Sun['Wav'][88:175])/1000
VULCAN_factor = VULCAN_TSI / WACCM_TSI

plt.figure(figsize = (10,5))
plt.plot(WACCM_wavelength, WACCM_Flux, color = 'orange', lw = 2, label = 'WACCM6')
#plt.plot(VULCAN_Sun['Wav'][88:175], VULCAN_Sun['Flux'][88:175]/VULCAN_factor, color = 'b', lw = 2, label = 'VULCAN')
plt.plot(VULCAN_Sun_EUV['Wav'], VULCAN_Sun_EUV['Flux']/VULCAN_factor, color = 'b', lw = 2)
plt.plot(wav/10, Flux_arr_corrected/Kasting_factor, color = 'black', lw = 2, label = 'Kasting')
#plt.plot(wav/10, 3*Flux_arr_corrected/Kasting_factor, color = 'black', lw = 2, label = 'Kasting')
#plt.plot(FarUV_Wav_arr/10, FarUV_SFX_arr_corrected/Kasting_factor, color = 'black', lw = 2)#, label = 'Kasting 1D model')
#plt.plot(FarUV_3D_Wav_arr/10, FarUV_SFX_3D_arr_corrected/Kasting_factor, color = 'black', lw = 2)#, label = 'Kasting 1D model')

plt.yscale('log')
plt.xlim(100,870)
plt.ylim(1e-3, 1e4)
thick_axes()
plt.legend(loc = 0, fontsize = 15, frameon = False)
plt.ylabel('Top of atmosphere irradiance\nper unit wavelength [W m'+r'$^\mathbf{-2}$'+' nm'+r'$^\mathbf{-1}$'+']',fontsize=15,color='k', weight = 'bold')
plt.xlabel('Wavelength [nm]', fontsize = 15, weight = 'bold')

#%% Youn Sun

File = '/Users/gregcooke/stellar_files/Sun_0.0Ga.txt'
Sun_0_0 = pd.read_csv(File, delim_whitespace=True, 
                      skiprows=2, names=['Wav', 'Flux'])

File = '/Users/gregcooke/stellar_files/Sun_0.5Ga.txt'
Sun_0_5 = pd.read_csv(File, delim_whitespace=True, 
                      skiprows=2, names=['Wav', 'Flux'])


File = '/Users/gregcooke/stellar_files/Sun_1.0Ga.txt'
Sun_1_0 = pd.read_csv(File, delim_whitespace=True, 
                      skiprows=2, names=['Wav', 'Flux'])

File = '/Users/gregcooke/stellar_files/Sun_2.0Ga.txt'
Sun_2_0 = pd.read_csv(File, delim_whitespace=True, 
                      skiprows=2, names=['Wav', 'Flux'])
File = '/Users/gregcooke/stellar_files/Sun_2.4Ga.txt'
Sun_2_4 = pd.read_csv(File, delim_whitespace=True, 
                      skiprows=2, names=['Wav', 'Flux'])

plt.figure(figsize = (10,5))

plt.plot(WACCM_wavelength, WACCM_Flux, color = 'orange', lw = 2, label = 'WACCM6')
plt.plot(Sun_0_0['Wav'], Sun_0_0['Flux'], color = 'black', lw = 2, label = 'Sun 0.0 Ga')
plt.plot(Sun_0_5['Wav'], Sun_0_5['Flux'], color = 'grey', lw = 2, label = 'Sun 0.5 Ga')
plt.plot(Sun_1_0['Wav'], Sun_1_0['Flux'], color = 'b', lw = 2, label = 'Sun 1.0 Ga')
plt.plot(Sun_2_0['Wav'], Sun_2_0['Flux'], color = 'cyan', lw = 2, label = 'Sun 2.0 Ga')
plt.plot(Sun_2_4['Wav'], Sun_2_4['Flux'], color = 'g', lw = 2, label = 'Sun 2.4 Ga')

#plt.yscale('log')
plt.xlim(10,1000)
#plt.ylim(1e-3)
#%% Cross section comparison

def load_photochem_photodiss_xs(species):
    """Effective photodissociation cross section [cm^2] from photochem_clima_data HDF5."""
    import h5py
    from photochem_clima_data import DATA_DIR

    filename = f"{DATA_DIR}/xsections/{species}.h5"
    with h5py.File(filename, "r") as f:
        w_nm = f["wavelengths"][:].astype(float)
        sigma = f["photodissociation"][:].astype(float)
        if "photodissociation-qy" in f:
            qy_group = f["photodissociation-qy"]
            w_qy = qy_group["wavelengths"][:].astype(float)
            sigma_eff = np.zeros_like(sigma)
            for branch in qy_group.keys():
                if branch == "wavelengths":
                    continue
                qy = np.interp(w_nm, w_qy, qy_group[branch][:].astype(float))
                sigma_eff += sigma * qy
            return w_nm, sigma_eff
        return w_nm, sigma


atmos_path = '/Users/gregcooke/atmos/PHOTOCHEM/DATA/XSECTIONS/'
csv_a_O2 = pd.read_csv(atmos_path+'O2/Yoshino92.XS.dat', delim_whitespace=True, skiprows=4, names=['Wav', 'XS'])
csv_a_O2 = pd.read_csv(atmos_path+'O2/O2_alinc.dat', delim_whitespace=True, skiprows=7, names=['Wav', 'XS'])
csv_a_O3 = pd.read_csv(atmos_path+'O3/O3.XS.dat', delim_whitespace=True, skiprows=4, names=['Wav', 'XS'])

# path = '/Users/gregcooke/K2_18_XS/'
# v_path = '/Users/gregcooke/VULCAN/thermo/photo_cross/'
# csv_a_H2O = pd.read_csv(path+'H2O.XS.dat', delim_whitespace=True, skiprows=4, names=['Wav', 'XS'])
# csv_a_CO2 = pd.read_csv(path+'CO2.XS.dat', delim_whitespace=True, skiprows=4, names=['Wav', 'XS'])
# csv_v_H2O = pd.read_csv(v_path+'H2O/H2O_cross.csv', skiprows=1, names=['Wav', 'XS', 'XS1', 'XS2'])
# csv_v_CO2 = pd.read_csv(v_path+'CO2/CO2_cross.csv', skiprows=1, names=['Wav', 'XS', 'XS1', 'XS2'])

w_pc_O2, xs_pc_O2 = load_photochem_photodiss_xs('O2')
w_pc_O3, xs_pc_O3 = load_photochem_photodiss_xs('O3')

wav = (range_min+range_max)/2
plt.figure(figsize = (11,5))
plt.plot(wav/10, ozone1_arr, 'teal', label = 'O3 (Kasting)')
plt.plot(wav/10, ozone2_arr, 'teal')
O3 = PI_V_482SZA['variable']['cross_J']['O3',1]+PI_V_482SZA['variable']['cross_J']['O3',2]
plt.plot(PI_V_482SZA['variable']['bins'], O3, 'm', lw=2,ls = '--', label = 'O3 (VULCAN)')
O2 = PI_V_482SZA['variable']['cross_J']['O2',1]+PI_V_482SZA['variable']['cross_J']['O2',2]
plt.plot(PI_V_482SZA['variable']['bins'], O2, 'm', lw=2,ls = '--', label = 'O2 (VULCAN)')
plt.plot((wav/10)[:35], o2_arr, color = 'teal', lw=2,label = 'O2 (Kasting)')
plt.plot(csv_a_O2['Wav']/10, csv_a_O2['XS'], color = 'orange',lw=2, ls = ':', label = 'O2 (atmos)')
plt.plot(csv_a_O3['Wav'], csv_a_O3['XS'], color = 'orange', lw=2,ls = ':', label = 'O2 (atmos)')
plt.plot(w_pc_O2, xs_pc_O2, color='b', ls='-.', lw=2, label='O2 (Photochem)')
plt.plot(w_pc_O3, xs_pc_O3, color='b', ls='-.', lw=2, label='O3 (Photochem)')
# plt.plot((wav/10)[:35], n2o_arr, color = 'purple', label = 'N2O')
# plt.plot((wav/10)[:35], co2_arr, color = 'g', lw = 2, ls = '-', label = 'CO2')
# plt.plot(csv_a_CO2['Wav']/10, csv_a_CO2['XS'], color = 'g', lw = 2, ls = '--')
# plt.plot(csv_v_CO2['Wav'], csv_v_CO2['XS'], color = 'g', lw = 2, ls = ':')
# plt.plot((wav/10)[:35], h2o_arr, color = 'b', lw = 2, ls = '-', label = 'H2O')
# plt.plot(csv_a_H2O['Wav']/10, csv_a_H2O['XS'], color = 'b', lw = 2, ls = '--')
# plt.plot(csv_v_H2O['Wav'], csv_v_H2O['XS'], color = 'b', lw = 2, ls = ':')
thick_axes()
plt.ylim(1e-24,1e-17)
plt.legend(loc = 0, fontsize = 15, frameon = False)
plt.xlim(150,250)
plt.yscale('log')
plt.ylabel('Cross section', fontsize = 15, weight = 'bold')
plt.xlabel('Wavelength [nm]', fontsize = 15, weight = 'bold')

#%%
# Original (pressure-model) grid
z_old = LWAV(Pre_pc_h0.Z3).values / 1000.0   # km
o3_old = LWAV(Pre_pc_h0.O3).values           # mixing ratio

# New (photochem) grid
z_new = Photo_PI['alt'].values                   # km

# Interpolate mixing ratio onto new altitude grid
o3_interp = np.interp(
    z_new,
    np.flip(z_old),
    np.flip(o3_old),
    left=np.nan,
    right=np.nan
)

plt.figure()
plt.plot(LWAV(Pre_pc_h0.O3), LWAV(Pre_pc_h0.Z3)/1000)
plt.plot(o3_interp, Photo_PI['alt'].values, color = 'r')
plt.plot(Photo_PI['O3'], Photo_PI['alt'].values, color = 'grey')
plt.xscale('log'); plt.xlim(1e-9, 1e-4)
plt.ylim(0, 140)

P_new_WACCM = Photo_PI
P_new_WACCM['O3'] = o3_interp

print(o3_trapz(P_new_WACCM))
print(o3_mass_conserving(P_new_WACCM))
print(o3_mass_conserving_trapz(P_new_WACCM, g=9.81, MO2=0.21))

#%% Truncation error

# Example usage
errors = layerwise_truncation_error(P_new_WACCM)


plt.figure()
plt.plot(errors['pct_error'], errors['z_mid'])
plt.xlabel('Layerwise % error (height vs pressure)')
plt.ylabel('Altitude (km)')
plt.grid(True)
plt.show()


# Example usage
total_error = total_pct_error(errors)
print(f"Total weighted % error: {total_error:f}%")



# Example usage
cumulative_column_and_error(errors)


n_m3 = P_new_WACCM['O3'].to_numpy() * P_new_WACCM['den'].to_numpy() * 1e6
z_m = P_new_WACCM['alt'].to_numpy() * 1e3
d2n_dz2 = np.gradient(np.gradient(n_m3, z_m), z_m)

plt.figure()
plt.plot(d2n_dz2, z_m/1e3)
plt.xlabel(r'$d^2 n_{O_3}/dz^2$ [molecules/m^5]')
plt.ylabel('Altitude [km]')
plt.grid(True)
plt.show()


#%%



print(o3_trapz(Photo_PI))
print(o3_trapz(Photo_10pc))
print(o3_trapz(Photo_1pc))
print(o3_trapz(Photo_01pc))

Photo_PI_col_mid = Photochem_O3_col_mid(Photo_PI, g=9.81, MO2 = 0.21)
Photo_PI_col = Photochem_O3_col(Photo_PI, g=9.81, MO2 = 0.21)
print(Photo_PI_col)
Photo_10pc_col = Photochem_O3_col(Photo_10pc, g=9.81, MO2 = 0.21*0.1)
print(Photo_10pc_col)
Photo_1pc_col = Photochem_O3_col(Photo_1pc, g=9.81, MO2 = 0.21*0.01)
print(Photo_1pc_col)
Photo_01pc_col = Photochem_O3_col(Photo_01pc, g=9.81, MO2 = 0.21*0.001)
print(Photo_01pc_col)
#%%
# VULCAN Ozone column 
'''
def V_O3_col_P_sum(DS):
    species = DS['variable']['species']
    O3 = DS['variable']['ymix'][:,species.index('O3')]
    p = DS['atm']['pico']/10
    mu = Earth_V_1pc['atm']['mu'] * 1.66e-27
    g = 9.81
    
    dp = np.empty(0)
    for i in range(len(p)-1):
        dp_new = p[i+1] - p[i]
        dp = np.append(dp, dp_new)
        
    O3_col = (O3 * dp / (g * mu)).sum() / 2.6867e20
    
    return -O3_col

def V_O3_col_P_trapz(DS):
    species = DS['variable']['species']
    O3 = DS['variable']['ymix'][:,species.index('O3')]
    p = DS['atm']['pco']/10
    mu = Earth_V_1pc['atm']['mu'] * 1.66e-27
    g = 9.81
    O3_gas = O3  / mu
        
    O3_col = (np.trapz(O3_gas, p)/(g)) / 2.6867e20
    
    return -O3_col
'''


'''
E_WPT_V_60SZA, E_WPT_V_spec_60SZA = Read_O3_Run(file_path='Earth_1e12s_60SZA_WPT_1rtol.vul')
E_WPT_V_45SZA, E_WPT_V_spec_45SZA = Read_O3_Run(file_path='Earth_1e12s_45SZA_WPT_1rtol.vul')

E_WPT_V_70pc_58SZA, Earth_V_spec_70pc_58SZA = Read_O3_Run(file_path='Earth_70pc_o2_1e12s_58SZA_WPT_1rtol.vul')
E_WPT_V_80pc_58SZA, Earth_V_spec_80pc_58SZA = Read_O3_Run(file_path='Earth_80pc_o2_1e12s_58SZA_WPT_1rtol.vul')
E_WPT_V_90pc_58SZA, Earth_V_spec_90pc_58SZA = Read_O3_Run(file_path='Earth_90pc_o2_1e12s_58SZA_WPT_1rtol.vul')
E_WPT_V_110pc_60SZA, Earth_V_spec_90pc_58SZA = Read_O3_Run(file_path='Earth_110pc_o2_1e12s_58SZA_WPT_1rtol.vul')
E_WPT_V_110pc_noD_60SZA, Earth_V_spec_90pc_58SZA = Read_O3_Run(file_path='Earth_110pc_o2_1e12s_60SZA_WPT_1rtol_no_Di_factor.vul')


E_WPT_V_125pc_58SZA, Earth_V_spec_125pc_58SZA = Read_O3_Run(file_path='Earth_125pc_o2_1e12s_58SZA_WPT_1rtol.vul')
E_WPT_V_140pc_58SZA, Earth_V_spec_140pc_58SZA = Read_O3_Run(file_path='Earth_140pc_o2_1e12s_58SZA_WPT_1rtol.vul')

Earth_V_60SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_60SZA)
Earth_V_58SZA_70pc_O3_col = V_O3_col_z_trapz(E_WPT_V_70pc_58SZA)
Earth_V_58SZA_80pc_O3_col = V_O3_col_z_trapz(E_WPT_V_80pc_58SZA)
Earth_V_58SZA_90pc_O3_col = V_O3_col_z_trapz(E_WPT_V_90pc_58SZA)
Earth_V_58SZA_125pc_O3_col = V_O3_col_z_trapz(E_WPT_V_125pc_58SZA)
Earth_V_58SZA_140pc_O3_col = V_O3_col_z_trapz(E_WPT_V_140pc_58SZA)
Earth_V_60SZA_110pc_O3_col = V_O3_col_z_trapz(E_WPT_V_110pc_60SZA)
Earth_V_60SZA_110pc_O3_col_noD_60SZA = V_O3_col_z_trapz(E_WPT_V_110pc_noD_60SZA)

E_WPT_V_150pc_60SZA, Earth_V_spec_150pc_60SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_60SZA_WPT_1rtol.vul')
E_WPT_V_150pc_45SZA, Earth_V_spec_150pc_45SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Earth_V_150pc, Earth_V_spec_150pc = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_58SZA_1rtol.vul')
Earth_V_150pc_45SZA, Earth_V_spec_150pc_45SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_45SZA_1rtol.vul')


Earth_V, Earth_V_spec = Read_O3_Run(file_path='Earth_1e12s_low_rtol.vul')
Earth_V_45SZA, Earth_V_spec_45SZA = Read_O3_Run(file_path='Earth_1e12s_45SZA_1rtol.vul')

E_WPT_V_50pc_60SZA, E_WPT_V_50pc_spec_60SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_60SZA_WPT_1rtol.vul')
E_WPT_V_50pc_45SZA, E_WPT_V_50pc_spec_45SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Earth_V_50pc, Earth_V_spec_50pc = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_58SZA_1rtol.vul')
Earth_V_50pc_45SZA, Earth_V_spec_50pc_45SZA = Read_O3_Run(file_path='Earth_50pc_o2_1e12s_45SZA_1rtol.vul')

E_WPT_V_10pc_60SZA, E_WPT_V_10pc_spec_60SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_60SZA_WPT_1rtol.vul')
E_WPT_V_10pc_45SZA, E_WPT_V_10pc_spec_45SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Earth_V_10pc, Earth_V_spec_10pc = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_low_rtol.vul')
Earth_V_10pc_45SZA, Earth_V_spec_10pc_45SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_45SZA_1rtol.vul')

E_WPT_V_5pc_60SZA, E_WPT_V_5pc_spec_60SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_60SZA_WPT_1rtol.vul')
E_WPT_V_5pc_45SZA, E_WPT_V_5pc_spec_45SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Earth_V_5pc, Earth_V_spec_5pc = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_58SZA_1rtol.vul')
Earth_V_5pc_45SZA, Earth_V_spec_5pc_45SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_45SZA_1rtol.vul')

E_WPT_V_1pc_60SZA, E_WPT_V_1pc_spec_60SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_60SZA_WPT_1rtol.vul')
E_WPT_V_1pc_45SZA, E_WPT_V_1pc_spec_45SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Earth_V_1pc, Earth_V_spec_1pc = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_low_rtol.vul')
Earth_V_1pc_45SZA, Earth_V_spec_1pc_45SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_45SZA_1rtol.vul')

E_WPT_V_05pc_60SZA, E_WPT_V_05pc_spec_60SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_60SZA_WPT_1rtol.vul')
E_WPT_V_05pc_45SZA, E_WPT_V_05pc_spec_45SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Earth_V_05pc, Earth_V_spec_05pc = Read_O3_Run(file_path='Earth_05pc_o2_1e12s_58SZA_1rtol.vul')
Earth_V_05pc_45SZA, Earth_V_spec_05pc_45SZA = Read_O3_Run(file_path='Earth_05pc_o2_1e12s_45SZA_1rtol.vul')

E_WPT_V_01pc_60SZA, E_WPT_V_01pc_spec_60SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_60SZA_WPT_1rtol.vul')
E_WPT_V_01pc_45SZA, E_WPT_V_01pc_spec_45SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Earth_V_01pc, Earth_V_spec_01pc = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_low_rtol.vul')
Earth_V_01pc, Earth_V_spec_01pc = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_60SZA_WPT_1rtol.vul')
Earth_V_01pc_45SZA, Earth_V_spec_01pc_45SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_45SZA_1rtol.vul')

E_WPT_V_60_01xCH4, E_WPT_V_60_01xCH4_spec = Read_O3_Run(file_path='Earth_0.1xCH4_1e12s_60SZA_WPT_1rtol.vul')
E_WPT_V_60_05xCH4, E_WPT_V_60_05xCH4_spec = Read_O3_Run(file_path='Earth_0.5xCH4_1e12s_60SZA_WPT_1rtol.vul')
E_WPT_V_60_1xCH4, E_WPT_V_60_1xCH4_spec = Read_O3_Run(file_path='Earth_1xCH4_1e12s_60SZA_WPT_1rtol.vul')
E_WPT_V_60_5xCH4, E_WPT_V_60_01xCH4_spec = Read_O3_Run(file_path='Earth_5xCH4_1e12s_60SZA_WPT_1rtol.vul')
E_WPT_V_60_10xCH4, E_WPT_V_60_05xCH4_spec = Read_O3_Run(file_path='Earth_10xCH4_1e12s_60SZA_WPT_1rtol.vul')

E_WPT_V_45_01xCH4, E_WPT_V_45_01xCH4_spec = Read_O3_Run(file_path='Earth_0.1xCH4_1e12s_45SZA_WPT_1rtol.vul')
E_WPT_V_45_05xCH4, E_WPT_V_45_05xCH4_spec = Read_O3_Run(file_path='Earth_0.5xCH4_1e12s_45SZA_WPT_1rtol.vul')
E_WPT_V_45_1xCH4, E_WPT_V_45_1xCH4_spec = Read_O3_Run(file_path='Earth_1xCH4_1e12s_45SZA_WPT_1rtol.vul')
E_WPT_V_45_5xCH4, E_WPT_V_45_5xCH4_spec = Read_O3_Run(file_path='Earth_5xCH4_1e12s_45SZA_WPT_1rtol.vul')
E_WPT_V_45_10xCH4, E_WPT_V_45_10xCH4_spec = Read_O3_Run(file_path='Earth_10xCH4_1e12s_45SZA_WPT_1rtol.vul')
E_WPT_V_45_20xCH4, E_WPT_V_45_20xCH4_spec = Read_O3_Run(file_path='Earth_20xCH4_1e12s_45SZA_WPT_1rtol.vul')
E_WPT_V_45_50xCH4, E_WPT_V_45_50xCH4_spec = Read_O3_Run(file_path='Earth_50xCH4_1e12s_45SZA_WPT_1rtol.vul')


P_test = False
trapz = True
if (P_test == True):
    if (trapz == True):
        Earth_V_58SZA_150pc_O3_col = V_O3_col_P_trapz(Earth_V_150pc)
    else:
        Earth_V_58SZA_150pc_O3_col = V_O3_col_P_sum(Earth_V_150pc)
else:
    if (trapz == True):
        Earth_V_58SZA_150pc_O3_col = V_O3_col_z_trapz(Earth_V_150pc)
        Earth_V_45SZA_150pc_O3_col = V_O3_col_z_trapz(Earth_V_150pc_45SZA)
        Earth_V_58SZA_O3_col = V_O3_col_z_trapz(Earth_V)
        Earth_V_45SZA_O3_col = V_O3_col_z_trapz(Earth_V_45SZA)
        Earth_V_58SZA_50pc_O3_col = V_O3_col_z_trapz(Earth_V_50pc)
        Earth_V_45SZA_50pc_O3_col = V_O3_col_z_trapz(Earth_V_50pc_45SZA)
        Earth_V_58SZA_10pc_O3_col = V_O3_col_z_trapz(Earth_V_10pc)
        Earth_V_45SZA_10pc_O3_col = V_O3_col_z_trapz(Earth_V_10pc_45SZA)
        Earth_V_58SZA_5pc_O3_col = V_O3_col_z_trapz(Earth_V_5pc)
        Earth_V_45SZA_5pc_O3_col = V_O3_col_z_trapz(Earth_V_5pc_45SZA)
        Earth_V_58SZA_1pc_O3_col = V_O3_col_z_trapz(Earth_V_1pc)
        Earth_V_45SZA_1pc_O3_col = V_O3_col_z_trapz(Earth_V_1pc_45SZA)
        Earth_V_58SZA_05pc_O3_col = V_O3_col_z_trapz(Earth_V_05pc)
        Earth_V_45SZA_05pc_O3_col = V_O3_col_z_trapz(Earth_V_05pc_45SZA)
        Earth_V_58SZA_01pc_O3_col = V_O3_col_z_trapz(Earth_V_01pc)
        Earth_V_45SZA_01pc_O3_col = V_O3_col_z_trapz(Earth_V_01pc_45SZA)
        
        E_WPT_V_01pc_45SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_01pc_45SZA)
        E_WPT_V_05pc_45SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_05pc_45SZA)
        E_WPT_V_1pc_45SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_1pc_45SZA)
        E_WPT_V_5pc_45SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_5pc_45SZA)
        E_WPT_V_10pc_45SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_10pc_45SZA)
        E_WPT_V_50pc_45SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_50pc_45SZA)
        E_WPT_V_45SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_45SZA)
        E_WPT_V_150pc_45SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_150pc_45SZA)
        
        E_WPT_V_01pc_60SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_01pc_60SZA)
        E_WPT_V_05pc_60SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_05pc_60SZA)
        E_WPT_V_1pc_60SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_1pc_60SZA)
        E_WPT_V_5pc_60SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_5pc_60SZA)
        E_WPT_V_10pc_60SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_10pc_60SZA)
        E_WPT_V_50pc_60SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_50pc_60SZA)
        E_WPT_V_60SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_60SZA)
        E_WPT_V_150pc_60SZA_O3_col = V_O3_col_z_trapz(E_WPT_V_150pc_60SZA)
        
        E_WPT_V_45_01xCH4_O3_col = V_O3_col_z_trapz(E_WPT_V_45_01xCH4)
        E_WPT_V_45_05xCH4_O3_col = V_O3_col_z_trapz(E_WPT_V_45_05xCH4)
        E_WPT_V_45_1xCH4_O3_col = V_O3_col_z_trapz(E_WPT_V_45_1xCH4)
        E_WPT_V_45_5xCH4_O3_col = V_O3_col_z_trapz(E_WPT_V_45_5xCH4)
        E_WPT_V_45_10xCH4_O3_col = V_O3_col_z_trapz(E_WPT_V_45_10xCH4)
        E_WPT_V_45_20xCH4_O3_col = V_O3_col_z_trapz(E_WPT_V_45_20xCH4)
        E_WPT_V_45_50xCH4_O3_col = V_O3_col_z_trapz(E_WPT_V_45_50xCH4)
        
        E_WPT_V_60_01xCH4_O3_col = V_O3_col_z_trapz(E_WPT_V_60_01xCH4)
        E_WPT_V_60_05xCH4_O3_col = V_O3_col_z_trapz(E_WPT_V_60_05xCH4)
        E_WPT_V_60_1xCH4_O3_col = V_O3_col_z_trapz(E_WPT_V_60_1xCH4)
        E_WPT_V_60_5xCH4_O3_col = V_O3_col_z_trapz(E_WPT_V_60_5xCH4)
        E_WPT_V_60_10xCH4_O3_col = V_O3_col_z_trapz(E_WPT_V_60_10xCH4)
        
    else:
        Earth_V_58SZA_150pc_O3_col = V_O3_col_z_sum(Earth_V_150pc)
'''
#%% Ozone column variation with different methane

'''
VULCAN columns
'''

One50_col = V_O3_col_z_trapz(One50_pc_V_482SZA)
One50_10xCH4_col = V_O3_col_z_trapz(One50_pc_V_10xCH4_482SZA)
One50_5xCH4_col = V_O3_col_z_trapz(One50_pc_V_5xCH4_482SZA)
One50_05xCH4_col = V_O3_col_z_trapz(One50_pc_V_05xCH4_482SZA)
One50_01xCH4_col = V_O3_col_z_trapz(One50_pc_V_01xCH4_482SZA)

One50_10xCH4f_col = V_O3_col_z_trapz(One50_pc_V_10xCH4f_482SZA)
One50_5xCH4f_col = V_O3_col_z_trapz(One50_pc_V_5xCH4f_482SZA)
One50_1xCH4f_col = V_O3_col_z_trapz(One50_pc_V_1xCH4f_482SZA)
One50_05xCH4f_col = V_O3_col_z_trapz(One50_pc_V_05xCH4f_482SZA)
One50_01xCH4f_col = V_O3_col_z_trapz(One50_pc_V_01xCH4f_482SZA)

PI_01xCH4f_col = V_O3_col_z_trapz(PI_V_01xCH4f_482SZA)
PI_05xCH4f_col = V_O3_col_z_trapz(PI_V_05xCH4f_482SZA)
PI_1xCH4f_col = V_O3_col_z_trapz(PI_V_1xCH4f_482SZA)
PI_5xCH4f_col = V_O3_col_z_trapz(PI_V_5xCH4f_482SZA)
PI_10xCH4f_col = V_O3_col_z_trapz(PI_V_10xCH4f_482SZA)

PI_col = V_O3_col_z_trapz(PI_V_482SZA)
PI_10xCH4_col = V_O3_col_z_trapz(PI_V_10xCH4_482SZA)
PI_5xCH4_col = V_O3_col_z_trapz(PI_V_5xCH4_482SZA)
PI_05xCH4_col = V_O3_col_z_trapz(PI_V_05xCH4_482SZA)
PI_01xCH4_col = V_O3_col_z_trapz(PI_V_01xCH4_482SZA)

Fifty_col = V_O3_col_z_trapz(Fifty_pc_V_482SZA)
Fifty_10xCH4_col = V_O3_col_z_trapz(Fifty_pc_V_10xCH4_482SZA)
Fifty_5xCH4_col = V_O3_col_z_trapz(Fifty_pc_V_5xCH4_482SZA)
Fifty_05xCH4_col = V_O3_col_z_trapz(Fifty_pc_V_05xCH4_482SZA)
Fifty_01xCH4_col = V_O3_col_z_trapz(Fifty_pc_V_01xCH4_482SZA)

Fifty_1xCH4f_col = V_O3_col_z_trapz(Fifty_pc_V_1xCH4f_482SZA)
Fifty_10xCH4f_col = V_O3_col_z_trapz(Fifty_pc_V_10xCH4f_482SZA)
Fifty_5xCH4f_col = V_O3_col_z_trapz(Fifty_pc_V_5xCH4f_482SZA)
Fifty_05xCH4f_col = V_O3_col_z_trapz(Fifty_pc_V_05xCH4f_482SZA)
Fifty_01xCH4f_col = V_O3_col_z_trapz(Fifty_pc_V_01xCH4f_482SZA)

Ten_col = V_O3_col_z_trapz(Ten_pc_V_482SZA)
Ten_col_T = V_O3_col_z_trapz(Ten_pc_V_482SZA_T)
Ten_10xCH4_col = V_O3_col_z_trapz(Ten_pc_V_10xCH4_482SZA)
Ten_5xCH4_col = V_O3_col_z_trapz(Ten_pc_V_5xCH4_482SZA)
Ten_05xCH4_col = V_O3_col_z_trapz(Ten_pc_V_05xCH4_482SZA)
Ten_01xCH4_col = V_O3_col_z_trapz(Ten_pc_V_01xCH4_482SZA)

Ten_1xCH4f_col = V_O3_col_z_trapz(Ten_pc_V_1xCH4f_482SZA)
Ten_10xCH4f_col = V_O3_col_z_trapz(Ten_pc_V_10xCH4f_482SZA)
Ten_5xCH4f_col = V_O3_col_z_trapz(Ten_pc_V_5xCH4f_482SZA)
Ten_05xCH4f_col = V_O3_col_z_trapz(Ten_pc_V_05xCH4f_482SZA)
Ten_01xCH4f_col = V_O3_col_z_trapz(Ten_pc_V_01xCH4f_482SZA)

Five_col = V_O3_col_z_trapz(Five_pc_V_482SZA)
Five_10xCH4_col = V_O3_col_z_trapz(Five_pc_V_10xCH4_482SZA)
Five_5xCH4_col = V_O3_col_z_trapz(Five_pc_V_5xCH4_482SZA)
Five_05xCH4_col = V_O3_col_z_trapz(Five_pc_V_05xCH4_482SZA)
Five_01xCH4_col = V_O3_col_z_trapz(Five_pc_V_01xCH4_482SZA)

Five_10xCH4f_col = V_O3_col_z_trapz(Five_pc_V_10xCH4f_482SZA)
Five_5xCH4f_col = V_O3_col_z_trapz(Five_pc_V_5xCH4f_482SZA)
Five_1xCH4f_col = V_O3_col_z_trapz(Five_pc_V_1xCH4f_482SZA)
Five_05xCH4f_col = V_O3_col_z_trapz(Five_pc_V_05xCH4f_482SZA)
Five_01xCH4f_col = V_O3_col_z_trapz(Five_pc_V_01xCH4f_482SZA)

One_col = V_O3_col_z_trapz(One_pc_V_482SZA)
One_col_T = V_O3_col_z_trapz(One_pc_V_482SZA_T)
One_10xCH4_col = V_O3_col_z_trapz(One_pc_V_10xCH4_482SZA)
One_5xCH4_col = V_O3_col_z_trapz(One_pc_V_5xCH4_482SZA)
One_05xCH4_col = V_O3_col_z_trapz(One_pc_V_05xCH4_482SZA)
One_01xCH4_col = V_O3_col_z_trapz(One_pc_V_01xCH4_482SZA)

One_10xCH4f_col = V_O3_col_z_trapz(One_pc_V_10xCH4f_482SZA)
One_5xCH4f_col = V_O3_col_z_trapz(One_pc_V_5xCH4f_482SZA)
One_1xCH4f_col = V_O3_col_z_trapz(One_pc_V_1xCH4f_482SZA)
One_05xCH4f_col = V_O3_col_z_trapz(One_pc_V_05xCH4f_482SZA)
One_01xCH4f_col = V_O3_col_z_trapz(One_pc_V_01xCH4f_482SZA)

Z5_col = V_O3_col_z_trapz(Z5_pc_V_482SZA)
Z5_10xCH4_col = V_O3_col_z_trapz(Z5_pc_V_10xCH4_482SZA)
Z5_5xCH4_col = V_O3_col_z_trapz(Z5_pc_V_5xCH4_482SZA)
Z5_05xCH4_col = V_O3_col_z_trapz(Z5_pc_V_05xCH4_482SZA)
Z5_01xCH4_col = V_O3_col_z_trapz(Z5_pc_V_01xCH4_482SZA)

Z5_10xCH4f_col = V_O3_col_z_trapz(Z5_pc_V_10xCH4f_482SZA)
Z5_5xCH4f_col = V_O3_col_z_trapz(Z5_pc_V_5xCH4f_482SZA)
Z5_1xCH4f_col = V_O3_col_z_trapz(Z5_pc_V_1xCH4f_482SZA)
Z5_05xCH4f_col = V_O3_col_z_trapz(Z5_pc_V_05xCH4f_482SZA)
Z5_01xCH4f_col = V_O3_col_z_trapz(Z5_pc_V_01xCH4f_482SZA)

Z1_col = V_O3_col_z_trapz(Z1_pc_V_482SZA)
Z1_col_T = V_O3_col_z_trapz(Z1_pc_V_482SZA_T)
Z1_10xCH4_col = V_O3_col_z_trapz(Z1_pc_V_10xCH4_482SZA)
Z1_5xCH4_col = V_O3_col_z_trapz(Z1_pc_V_5xCH4_482SZA)
Z1_05xCH4_col = V_O3_col_z_trapz(Z1_pc_V_05xCH4_482SZA)
Z1_01xCH4_col = V_O3_col_z_trapz(Z1_pc_V_01xCH4_482SZA)

Z1_10xCH4f_col = V_O3_col_z_trapz(Z1_pc_V_10xCH4f_482SZA)
Z1_5xCH4f_col = V_O3_col_z_trapz(Z1_pc_V_5xCH4f_482SZA)
Z1_1xCH4f_col = V_O3_col_z_trapz(Z1_pc_V_1xCH4f_482SZA)
Z1_05xCH4f_col = V_O3_col_z_trapz(Z1_pc_V_05xCH4f_482SZA)
Z1_01xCH4f_col = V_O3_col_z_trapz(Z1_pc_V_01xCH4f_482SZA)

# Data mapping to make looping possible
# Format: (Title, Methane_Data_List, Column_Data_List, has_fixed_flux)
plot_configs = [
    ('150% PAL', [0.1, 0.5, 1, 5], [One50_01xCH4_col, One50_05xCH4_col, One50_col, One50_5xCH4_col], True),
    ('PI', [0.1, 0.5, 1, 5, 10], [PI_01xCH4_col, PI_05xCH4_col, PI_col, PI_5xCH4_col, PI_10xCH4_col], True),
    ('50% PAL', [0.1, 0.5, 1, 5, 10], [Fifty_01xCH4_col, Fifty_05xCH4_col, Fifty_col, Fifty_5xCH4_col, Fifty_10xCH4_col], True),
    ('10% PAL', [0.1, 0.5, 1, 5, 10], [Ten_01xCH4_col, Ten_05xCH4_col, Ten_col, Ten_5xCH4_col, Ten_10xCH4_col], True),
    ('5% PAL', [0.1, 0.5, 1, 5, 10], [Five_01xCH4_col, Five_05xCH4_col, Five_col, Five_5xCH4_col, Five_10xCH4_col], True),
    ('1% PAL', [0.1, 0.5, 1, 5, 10], [One_01xCH4_col, One_05xCH4_col, One_col, One_5xCH4_col, One_10xCH4_col], True),
    ('0.5% PAL', [0.1, 0.5, 1, 5, 10], [Z5_01xCH4_col, Z5_05xCH4_col, Z5_col, Z5_5xCH4_col, Z5_10xCH4_col], True),
    ('0.1% PAL', [0.1, 0.5, 1, 5, 10], [Z1_01xCH4_col, Z1_05xCH4_col, Z1_col, Z1_5xCH4_col, Z1_10xCH4_col], True)
]

'''
Photochem columns
'''

Photo_PI_col = Photochem_O3_col(Photo_PI, g=9.81, MO2 = 0.21)
Photo_PI_01CH4_col = Photochem_O3_col(Photo_PI_01CH4, g=9.81, MO2 = 0.21)
Photo_PI_05CH4_col = Photochem_O3_col(Photo_PI_05CH4, g=9.81, MO2 = 0.21)
Photo_PI_1CH4_col = Photochem_O3_col(Photo_PI_1CH4, g=9.81, MO2 = 0.21)
Photo_PI_2CH4_col = Photochem_O3_col(Photo_PI_2CH4, g=9.81, MO2 = 0.21)
Photo_PI_3CH4_col = Photochem_O3_col(Photo_PI_3CH4, g=9.81, MO2 = 0.21)
Photo_PI_4CH4_col = Photochem_O3_col(Photo_PI_4CH4, g=9.81, MO2 = 0.21)
Photo_PI_5CH4_col = Photochem_O3_col(Photo_PI_5CH4, g=9.81, MO2 = 0.21)

Photo_10pc_col = Photochem_O3_col(Photo_10pc, g=9.81, MO2 = 0.021)
Photo_10pc_01CH4_col = Photochem_O3_col(Photo_10pc_01CH4, g=9.81, MO2 = 0.021)
Photo_10pc_05CH4_col = Photochem_O3_col(Photo_10pc_05CH4, g=9.81, MO2 = 0.021)
Photo_10pc_1CH4_col = Photochem_O3_col(Photo_10pc_1CH4, g=9.81, MO2 = 0.021)
Photo_10pc_2CH4_col = Photochem_O3_col(Photo_10pc_2CH4, g=9.81, MO2 = 0.021)
Photo_10pc_3CH4_col = Photochem_O3_col(Photo_10pc_3CH4, g=9.81, MO2 = 0.021)
Photo_10pc_4CH4_col = Photochem_O3_col(Photo_10pc_4CH4, g=9.81, MO2 = 0.021)
Photo_10pc_5CH4_col = Photochem_O3_col(Photo_10pc_5CH4, g=9.81, MO2 = 0.021)

Photo_1pc_col = Photochem_O3_col(Photo_1pc, g=9.81, MO2 = 0.0021)
Photo_1pc_01CH4_col = Photochem_O3_col(Photo_1pc_01CH4, g=9.81, MO2 = 0.0021)
Photo_1pc_05CH4_col = Photochem_O3_col(Photo_1pc_05CH4, g=9.81, MO2 = 0.0021)
Photo_1pc_1CH4_col = Photochem_O3_col(Photo_1pc_1CH4, g=9.81, MO2 = 0.0021)
Photo_1pc_2CH4_col = Photochem_O3_col(Photo_1pc_2CH4, g=9.81, MO2 = 0.0021)
Photo_1pc_3CH4_col = Photochem_O3_col(Photo_1pc_3CH4, g=9.81, MO2 = 0.0021)
Photo_1pc_4CH4_col = Photochem_O3_col(Photo_1pc_4CH4, g=9.81, MO2 = 0.0021)
Photo_1pc_5CH4_col = Photochem_O3_col(Photo_1pc_5CH4, g=9.81, MO2 = 0.0021)

Photo_01pc_col = Photochem_O3_col(Photo_01pc, g=9.81, MO2 = 0.00021)
Photo_01pc_01CH4_col = Photochem_O3_col(Photo_01pc_01CH4, g=9.81, MO2 = 0.00021)
Photo_01pc_05CH4_col = Photochem_O3_col(Photo_01pc_05CH4, g=9.81, MO2 = 0.00021)
Photo_01pc_1CH4_col = Photochem_O3_col(Photo_01pc_1CH4, g=9.81, MO2 = 0.00021)
Photo_01pc_2CH4_col = Photochem_O3_col(Photo_01pc_2CH4, g=9.81, MO2 = 0.00021)
Photo_01pc_3CH4_col = Photochem_O3_col(Photo_01pc_3CH4, g=9.81, MO2 = 0.00021)
Photo_01pc_4CH4_col = Photochem_O3_col(Photo_01pc_4CH4, g=9.81, MO2 = 0.00021)
Photo_01pc_5CH4_col = Photochem_O3_col(Photo_01pc_5CH4, g=9.81, MO2 = 0.00021)

'''
Plot figure
'''

fig = plt.figure(figsize=(22, 10))
gs = gridspec.GridSpec(2, 4, hspace=0.15, wspace=0.25)

for i, (title, meth_vals, col_vals, has_flux) in enumerate(plot_configs):
    ax = fig.add_subplot(gs[i // 4, i % 4])
    
    # Plot standard Fixed MR
    methane_array = np.array(meth_vals) * 0.8e-6
    ax.plot(methane_array, col_vals, lw=lw,  marker='s', ls='', color='m', label='VULCAN')
    
    # Special case for PI Fixed Flux
    if has_flux:
        # Assuming these variables are defined in your global scope
        if (i == 0):
            flux_methane = [
                One50_pc_V_10xCH4f_482SZA['variable']['ymix'][:,One50_pc_spec_10xCH4f_482SZA.index('CH4')][0],
                One50_pc_V_5xCH4f_482SZA['variable']['ymix'][:,One50_pc_spec_5xCH4f_482SZA.index('CH4')][0],
                One50_pc_V_1xCH4f_482SZA['variable']['ymix'][:,One50_pc_spec_1xCH4f_482SZA.index('CH4')][0],
                One50_pc_V_05xCH4f_482SZA['variable']['ymix'][:,One50_pc_spec_05xCH4f_482SZA.index('CH4')][0],
                One50_pc_V_01xCH4f_482SZA['variable']['ymix'][:,One50_pc_spec_01xCH4f_482SZA.index('CH4')][0]
            ]
            flux_cols = [One50_10xCH4f_col, One50_5xCH4f_col, One50_1xCH4f_col, One50_05xCH4f_col, One50_01xCH4f_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='m', marker='s', ls='')#, label='VULCAN')

        if (i == 1):
            flux_methane = [
                PI_V_01xCH4f_482SZA['variable']['ymix'][:,PI_pc_spec_01xCH4f_482SZA.index('CH4')][0],
                PI_V_05xCH4f_482SZA['variable']['ymix'][:,PI_pc_spec_05xCH4f_482SZA.index('CH4')][0],
                PI_V_1xCH4f_482SZA['variable']['ymix'][:,PI_pc_spec_1xCH4f_482SZA.index('CH4')][0],
                PI_V_5xCH4f_482SZA['variable']['ymix'][:,PI_pc_spec_5xCH4f_482SZA.index('CH4')][0],
                PI_V_10xCH4f_482SZA['variable']['ymix'][:,PI_pc_spec_10xCH4f_482SZA.index('CH4')][0]
            ]
            flux_cols = [PI_01xCH4f_col, PI_05xCH4f_col, PI_1xCH4f_col, PI_5xCH4f_col, PI_10xCH4f_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='m', marker='s', ls='')#, label='Fixed Flux, VULCAN')
            
            flux_methane = [
                Photo_PI_01CH4['CH4'][0],Photo_PI_05CH4['CH4'][0],Photo_PI_1CH4['CH4'][0],
                Photo_PI_2CH4['CH4'][0], Photo_PI_3CH4['CH4'][0], Photo_PI_4CH4['CH4'][0],
                Photo_PI_5CH4['CH4'][0]
            ]
            flux_cols = [Photo_PI_01CH4_col, Photo_PI_05CH4_col, Photo_PI_1CH4_col, Photo_PI_2CH4_col,
                         Photo_PI_3CH4_col, Photo_PI_4CH4_col, Photo_PI_5CH4_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='b', marker='s', ls='', label='Photochem')
            
            flux_cols = [251.93, 256.11,259.7,264.68,268.17,270.79,272.83]
            ax.plot(0.808e-6*np.array([0.1,.5,1,2,3,4,5]), flux_cols, lw=lw, color=color_Atmos, marker='s', ls='', label='Atmos')
            
        if (i == 2):
            flux_methane = [
                Fifty_pc_V_01xCH4f_482SZA['variable']['ymix'][:,Fifty_pc_spec_01xCH4f_482SZA.index('CH4')][0],
                Fifty_pc_V_05xCH4f_482SZA['variable']['ymix'][:,Fifty_pc_spec_05xCH4f_482SZA.index('CH4')][0],
                Fifty_pc_V_1xCH4f_482SZA['variable']['ymix'][:,Fifty_pc_spec_1xCH4f_482SZA.index('CH4')][0],
                Fifty_pc_V_5xCH4f_482SZA['variable']['ymix'][:,Fifty_pc_spec_5xCH4f_482SZA.index('CH4')][0],
                Fifty_pc_V_10xCH4f_482SZA['variable']['ymix'][:,Fifty_pc_spec_10xCH4f_482SZA.index('CH4')][0]
                ]
        
            flux_cols = [Fifty_01xCH4f_col, Fifty_05xCH4f_col, Fifty_1xCH4f_col, Fifty_5xCH4f_col, Fifty_10xCH4f_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='m', marker='s', ls='', label='Fixed Flux, VULCAN')


            
        if (i == 3):
            flux_methane = [
                Ten_pc_V_01xCH4f_482SZA['variable']['ymix'][:,Ten_pc_spec_01xCH4f_482SZA.index('CH4')][0],
                Ten_pc_V_05xCH4f_482SZA['variable']['ymix'][:,Ten_pc_spec_05xCH4f_482SZA.index('CH4')][0],
                Ten_pc_V_1xCH4f_482SZA['variable']['ymix'][:,Ten_pc_spec_1xCH4f_482SZA.index('CH4')][0],
                Ten_pc_V_5xCH4f_482SZA['variable']['ymix'][:,Ten_pc_spec_5xCH4f_482SZA.index('CH4')][0],
                Ten_pc_V_10xCH4f_482SZA['variable']['ymix'][:,Ten_pc_spec_10xCH4f_482SZA.index('CH4')][0]
                ]
            flux_cols = [Ten_01xCH4f_col, Ten_05xCH4f_col, Ten_1xCH4f_col, Ten_5xCH4f_col, Ten_10xCH4f_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='m', marker='s', ls='', label='Fixed Flux, VULCAN')
        
            flux_methane = [
                Photo_10pc_01CH4['CH4'][0],Photo_10pc_05CH4['CH4'][0],Photo_10pc_1CH4['CH4'][0],
                Photo_10pc_2CH4['CH4'][0], Photo_10pc_3CH4['CH4'][0], Photo_10pc_4CH4['CH4'][0],
                Photo_10pc_5CH4['CH4'][0]
            ]
            flux_cols = [Photo_10pc_01CH4_col, Photo_10pc_05CH4_col, Photo_10pc_1CH4_col, Photo_10pc_2CH4_col,
                         Photo_10pc_3CH4_col, Photo_10pc_4CH4_col, Photo_10pc_5CH4_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='b', marker='s', ls='', label='Photochem')
            
        if (i == 4):
            flux_methane = [
                Five_pc_V_01xCH4f_482SZA['variable']['ymix'][:,Five_pc_spec_01xCH4f_482SZA.index('CH4')][0],
                Five_pc_V_05xCH4f_482SZA['variable']['ymix'][:,Five_pc_spec_05xCH4f_482SZA.index('CH4')][0],
                Five_pc_V_1xCH4f_482SZA['variable']['ymix'][:,Five_pc_spec_1xCH4f_482SZA.index('CH4')][0],
                Five_pc_V_5xCH4f_482SZA['variable']['ymix'][:,Five_pc_spec_5xCH4f_482SZA.index('CH4')][0],
                Five_pc_V_10xCH4f_482SZA['variable']['ymix'][:,Five_pc_spec_10xCH4f_482SZA.index('CH4')][0]
                ]
            flux_cols = [Five_01xCH4f_col, Five_05xCH4f_col, Five_1xCH4f_col, Five_5xCH4f_col, Five_10xCH4f_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='m', marker='s', ls='', label='Fixed Flux, VULCAN')
        
        if (i == 5):
            flux_methane = [
                One_pc_V_01xCH4f_482SZA['variable']['ymix'][:,One_pc_spec_01xCH4f_482SZA.index('CH4')][0],
                One_pc_V_05xCH4f_482SZA['variable']['ymix'][:,One_pc_spec_05xCH4f_482SZA.index('CH4')][0],
                One_pc_V_1xCH4f_482SZA['variable']['ymix'][:,One_pc_spec_1xCH4f_482SZA.index('CH4')][0],
                One_pc_V_5xCH4f_482SZA['variable']['ymix'][:,One_pc_spec_5xCH4f_482SZA.index('CH4')][0],
                One_pc_V_10xCH4f_482SZA['variable']['ymix'][:,One_pc_spec_10xCH4f_482SZA.index('CH4')][0]
                ]
            flux_cols = [One_01xCH4f_col, One_05xCH4f_col, One_1xCH4f_col, One_5xCH4f_col, One_10xCH4f_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='m', marker='s', ls='', label='Fixed Flux, VULCAN')
        
            flux_methane = [
                Photo_1pc_01CH4['CH4'][0],Photo_1pc_05CH4['CH4'][0],Photo_1pc_1CH4['CH4'][0],
                Photo_1pc_2CH4['CH4'][0], Photo_1pc_3CH4['CH4'][0], Photo_1pc_4CH4['CH4'][0],
                Photo_1pc_5CH4['CH4'][0]
            ]
            flux_cols = [Photo_1pc_01CH4_col, Photo_1pc_05CH4_col, Photo_1pc_1CH4_col, Photo_1pc_2CH4_col,
                         Photo_1pc_3CH4_col, Photo_1pc_4CH4_col, Photo_1pc_5CH4_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='b', marker='s', ls='', label='Photochem')
            
        
        if (i == 6):
            flux_methane = [
                Z5_pc_V_01xCH4f_482SZA['variable']['ymix'][:,Z5_pc_spec_01xCH4f_482SZA.index('CH4')][0],
                Z5_pc_V_05xCH4f_482SZA['variable']['ymix'][:,Z5_pc_spec_05xCH4f_482SZA.index('CH4')][0],
                Z5_pc_V_1xCH4f_482SZA['variable']['ymix'][:,Z5_pc_spec_1xCH4f_482SZA.index('CH4')][0],
                Z5_pc_V_5xCH4f_482SZA['variable']['ymix'][:,Z5_pc_spec_5xCH4f_482SZA.index('CH4')][0],
                Z5_pc_V_10xCH4f_482SZA['variable']['ymix'][:,Z5_pc_spec_10xCH4f_482SZA.index('CH4')][0]
                ]
            flux_cols = [Z5_01xCH4f_col, Z5_05xCH4f_col, Z5_1xCH4f_col, Z5_5xCH4f_col, Z5_10xCH4f_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='m', marker='s', ls='', label='Fixed Flux, VULCAN')
        
        if (i == 7):
            flux_methane = [
                Z1_pc_V_01xCH4f_482SZA['variable']['ymix'][:,Z1_pc_spec_01xCH4f_482SZA.index('CH4')][0],
                Z1_pc_V_05xCH4f_482SZA['variable']['ymix'][:,Z1_pc_spec_05xCH4f_482SZA.index('CH4')][0],
                Z1_pc_V_1xCH4f_482SZA['variable']['ymix'][:,Z1_pc_spec_1xCH4f_482SZA.index('CH4')][0],
                Z1_pc_V_5xCH4f_482SZA['variable']['ymix'][:,Z1_pc_spec_5xCH4f_482SZA.index('CH4')][0],
                Z1_pc_V_10xCH4f_482SZA['variable']['ymix'][:,Z1_pc_spec_10xCH4f_482SZA.index('CH4')][0]
                ]
            flux_cols = [Z1_01xCH4f_col, Z1_05xCH4f_col, Z1_1xCH4f_col, Z1_5xCH4f_col, Z1_10xCH4f_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='m', marker='s', ls='', label='Fixed Flux, VULCAN')
            
            flux_methane = [
                Photo_01pc_01CH4['CH4'][0],Photo_01pc_05CH4['CH4'][0],Photo_01pc_1CH4['CH4'][0],
                Photo_01pc_2CH4['CH4'][0], Photo_01pc_3CH4['CH4'][0], Photo_01pc_4CH4['CH4'][0],
                Photo_01pc_5CH4['CH4'][0]
            ]
            flux_cols = [Photo_01pc_01CH4_col, Photo_01pc_05CH4_col, Photo_01pc_1CH4_col, Photo_01pc_2CH4_col,
                         Photo_01pc_3CH4_col, Photo_01pc_4CH4_col, Photo_01pc_5CH4_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='b', marker='s', ls='', label='Photochem')
            

    # Formatting
    ax.set_xscale('log')
    ax.set_title(title, fontsize=15, weight='bold')
    thick_axes(ax) # Pass ax if your function supports it, else call as you did
    
    # Only label outer axes to keep it clean
    if i >= 4: # Bottom row
        ax.set_xlabel(f'Mixing ratio of CH$_{4}$ at surface', fontsize=12, weight='bold')
    if i % 4 == 0: # Leftmost column
        ax.set_ylabel('Ozone column [DU]', fontsize=12, weight='bold')
    
    # Place legend only on the first or second plot to avoid clutter
    if i == 1:
        ax.legend(loc=(0), fontsize=15, frameon = False)

plt.savefig('/Users/gregcooke/python_output/Varying_CH4_O3_col.png')


#%% TRAPPIST-1 e files
TP1e_01pc_MS_60SZA, TP1e_01pc_MS_60SZA_spec = Read_O3_Run(file_path='TP1e_MS_0.1pc_o2_1e12s_60SZA_WBC_WPT_1rtol.vul')
#TP1e_01pc_W21_60SZA, TP1e_01pc_W21_60SZA_spec = Read_O3_Run(file_path='TP1e_W21_01pc_o2_1e12s_60SZA_WBC_WPT_1rtol.vul')
#TP1e_01pc_P19_60SZA, TP1e_01pc_P19_60SZA_spec = Read_O3_Run(file_path='TP1e_P19_01pc_o2_1e12s_60SZA_WBC_WPT_1rtol.vul')
TP1e_01pc_MS_60SZA_col = V_O3_col_z_trapz(TP1e_01pc_MS_60SZA); #TP1e_01pc_W21_60SZA_col = V_O3_col_z_trapz(TP1e_01pc_W21_60SZA); #TP1e_01pc_P19_60SZA_col = V_O3_col_z_trapz(TP1e_01pc_P19_60SZA)

TP1e_1pc_MS_60SZA, TP1e_1pc_MS_60SZA_spec = Read_O3_Run(file_path='TP1e_MS_1pc_o2_1e12s_60SZA_WBC_WPT_1rtol.vul')
TP1e_1pc_W21_60SZA, TP1e_1pc_W21_60SZA_spec = Read_O3_Run(file_path='TP1e_W21_1pc_o2_1e12s_60SZA_WBC_WPT_1rtol.vul')
TP1e_1pc_P19_60SZA, TP1e_1pc_P19_60SZA_spec = Read_O3_Run(file_path='TP1e_P19_1pc_o2_1e12s_60SZA_WBC_WPT_1rtol.vul')
TP1e_1pc_MS_60SZA_col = V_O3_col_z_trapz(TP1e_1pc_MS_60SZA); TP1e_1pc_W21_60SZA_col = V_O3_col_z_trapz(TP1e_1pc_W21_60SZA); TP1e_1pc_P19_60SZA_col = V_O3_col_z_trapz(TP1e_1pc_P19_60SZA)

TP1e_10pc_MS_60SZA, TP1e_10pc_MS_60SZA_spec = Read_O3_Run(file_path='TP1e_MS_10pc_o2_1e12s_60SZA_WBC_WPT_1rtol.vul')
TP1e_10pc_W21_60SZA, TP1e_10pc_W21_60SZA_spec = Read_O3_Run(file_path='TP1e_W21_10pc_o2_1e12s_60SZA_WBC_WPT_1rtol.vul')
TP1e_10pc_P19_60SZA, TP1e_10pc_P19_60SZA_spec = Read_O3_Run(file_path='TP1e_P19_10pc_o2_1e12s_60SZA_WBC_WPT_1rtol.vul')
TP1e_10pc_MS_60SZA_col = V_O3_col_z_trapz(TP1e_10pc_MS_60SZA); TP1e_10pc_W21_60SZA_col = V_O3_col_z_trapz(TP1e_10pc_W21_60SZA); TP1e_10pc_P19_60SZA_col = V_O3_col_z_trapz(TP1e_10pc_P19_60SZA)

TP1e_MS_60SZA, TP1e_MS_60SZA_spec = Read_O3_Run(file_path='TP1e_MS_1e12s_60SZA_WBC_WPT_1rtol.vul')
TP1e_W21_60SZA, TP1e_W21_60SZA_spec = Read_O3_Run(file_path='TP1e_W21_1e12s_60SZA_WBC_WPT_1rtol.vul')
TP1e_P19_60SZA, TP1e_P19_60SZA_spec = Read_O3_Run(file_path='TP1e_P19_1e12s_60SZA_WBC_WPT_1rtol.vul')
TP1e_MS_60SZA_col = V_O3_col_z_trapz(TP1e_MS_60SZA); TP1e_W21_60SZA_col = V_O3_col_z_trapz(TP1e_W21_60SZA); TP1e_P19_60SZA_col = V_O3_col_z_trapz(TP1e_P19_60SZA)

plt.figure()
O2 = [1, 0.1, 0.01, 0.001]   
cols = TP1e_MS_60SZA_col, TP1e_10pc_MS_60SZA_col, TP1e_1pc_MS_60SZA_col, TP1e_01pc_MS_60SZA_col
plt.plot(O2, cols, marker = 's', color = 'black')
O2 = [1, 0.1, 0.01] #, 0.001] 
cols = TP1e_P19_60SZA_col, TP1e_10pc_P19_60SZA_col, TP1e_1pc_P19_60SZA_col
plt.plot(O2, cols, marker = 's', color = 'orange')
O2 = [1, 0.1, 0.01]#, 0.001]  
cols = TP1e_W21_60SZA_col, TP1e_10pc_W21_60SZA_col, TP1e_1pc_W21_60SZA_col
plt.plot(O2, cols, marker = 's', color = 'b')
plt.yscale('log')
plt.xscale('log')

plt.figure(figsize = (10,6))
gs = gridspec.GridSpec(2, 2)
plt.subplot(gs[0,0])
plt.plot(TP1e_MS_60SZA['variable']['ymix'][:,TP1e_MS_60SZA_spec.index('O3')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'black', ls = '--', label = 'O'+sub(3))
plt.plot(TP1e_P19_60SZA['variable']['ymix'][:,TP1e_MS_60SZA_spec.index('O3')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'black', ls = '--', label = 'O'+sub(3))
plt.plot(LWAV(Pre_pc_h0.O3), Pre_pc_h0.lev/1e3, color = 'black', lw = 2)
plt.xscale('log')
plt.xlim(1e-9, 1e-3)
plt.ylim(1e0, 1e-7)
plt.yscale('log')


#%% Ozone column plot

#Photo_150pc_60SZA_col = calculate_column_height_method(Photo_150pc_60SZA)
Photo_PI_60SZA_col = calculate_column_height_method(Photo_PI_60SZA)
Photo_10pc_60SZA_col = calculate_column_height_method(Photo_10pc_60SZA)
Photo_1pc_60SZA_col = calculate_column_height_method(Photo_1pc_60SZA)
Photo_01pc_60SZA_col = calculate_column_height_method(Photo_01pc_60SZA)

#Photo_150pc_45SZA_col = calculate_column_height_method(Photo_150pc_45SZA)
Photo_PI_45SZA_col = calculate_column_height_method(Photo_PI_45SZA)
Photo_10pc_45SZA_col = calculate_column_height_method(Photo_10pc_45SZA)
Photo_1pc_45SZA_col = calculate_column_height_method(Photo_1pc_45SZA)
Photo_01pc_45SZA_col = calculate_column_height_method(Photo_01pc_45SZA)

alpha = 0.25
photochem_60 = [Photo_01pc_60SZA_col, 106.851, Photo_1pc_60SZA_col, 191.3, Photo_10pc_60SZA_col, 265.73, Photo_PI_60SZA_col, 253.26]
photochem_45 = [Photo_01pc_45SZA_col, 125.3, Photo_1pc_45SZA_col, 230.58, Photo_10pc_45SZA_col, 346.74, Photo_PI_45SZA_col, 365.15]

o2_conc_less = [0.001, 0.01, 0.1, 1, 1.5]
o2_conc_photo = [0.001, 0.01, 0.1, 1]
kasting_45sza = [29.89765538, 106.2858206, 247.7112021, 353.0815035, 352.3557871]
kasting_60sza = [28.87458132, 93.61369557, 205.1246744, 264.4175661, 256.4346855]

#Atmos
A_SZA_60 = [80.38, 130.55, 148.79, 175.8, 191.42, 199.49, 204.83, 200.47]
A_SZA_45 = [92.6, 160, 178.54, 217.52, 239.41, 261.56, 272.85, 272.57]

'''
160 IS WRIBG FOR ATMOS 0.5 at 60 SZA need to change
'''

#VULCAN

One50_pc_V_45SZA_O3_col = V_O3_col_z_trapz(One50_pc_V_45SZA)
PI_V_45SZA_O3_col = V_O3_col_z_trapz(PI_V_45SZA)
Fifty_pc_V_45SZA_O3_col = V_O3_col_z_trapz(Fifty_pc_V_45SZA)
Ten_pc_V_45SZA_O3_col = V_O3_col_z_trapz(Ten_pc_V_45SZA)
One_pc_V_45SZA_O3_col = V_O3_col_z_trapz(One_pc_V_45SZA)
Five_pc_V_45SZA_O3_col = V_O3_col_z_trapz(Five_pc_V_45SZA)
Z5_pc_V_45SZA_O3_col = V_O3_col_z_trapz(Z5_pc_V_45SZA)
Z1_pc_V_45SZA_O3_col = V_O3_col_z_trapz(Z1_pc_V_45SZA)

One50_pc_V_60SZA_O3_col = V_O3_col_z_trapz(One50_pc_V_60SZA)
PI_V_60SZA_O3_col = V_O3_col_z_trapz(PI_V_60SZA)
Fifty_pc_V_60SZA_O3_col = V_O3_col_z_trapz(Fifty_pc_V_60SZA)
Ten_pc_V_60SZA_O3_col = V_O3_col_z_trapz(Ten_pc_V_60SZA)
One_pc_V_60SZA_O3_col = V_O3_col_z_trapz(One_pc_V_60SZA)
Five_pc_V_60SZA_O3_col = V_O3_col_z_trapz(Five_pc_V_60SZA)
Z5_pc_V_60SZA_O3_col = V_O3_col_z_trapz(Z5_pc_V_60SZA)
Z1_pc_V_60SZA_O3_col = V_O3_col_z_trapz(Z1_pc_V_60SZA)

V_SZA_60 = [Z1_pc_V_60SZA_O3_col, Z5_pc_V_60SZA_O3_col, One_pc_V_60SZA_O3_col, Five_pc_V_60SZA_O3_col,
            Ten_pc_V_60SZA_O3_col, Fifty_pc_V_60SZA_O3_col, PI_V_60SZA_O3_col, One50_pc_V_60SZA_O3_col]
V_SZA_45 = [Z1_pc_V_45SZA_O3_col, Z5_pc_V_45SZA_O3_col, One_pc_V_45SZA_O3_col, Five_pc_V_45SZA_O3_col,
            Ten_pc_V_45SZA_O3_col, Fifty_pc_V_45SZA_O3_col, PI_V_45SZA_O3_col, One50_pc_V_45SZA_O3_col]

o2_conc = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 1.5]
WACCM = [LWAV(Zero1_col), LWAV(Zero5_col), LWAV(One_col), LWAV(Five_col), LWAV(Ten_col), LWAV(Fifty_col), LWAV(Pre_col), LWAV(One50_col)]
WACCM_max = [Zero1_col.max(), Zero5_col.max(), One_col.max(), Five_col.max(), Ten_col.max(), Fifty_col.max(), Pre_col.max(), One50_col.max()]
WACCM_min = [Zero1_col.min(), Zero5_col.min(), One_col.min(), Five_col.min(), Ten_col.min(), Fifty_col.min(), Pre_col.min(), One50_col.min()]

plt.figure(figsize = (12,6))
plt.fill_between(o2_conc, photochem_60, photochem_45, alpha = alpha, color = 'b', label = 'Photochem, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, V_SZA_45, V_SZA_60, alpha = alpha, color = 'm', label = 'VULCAN, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, WACCM_min, WACCM_max, alpha = alpha, color = 'black', label = 'WACCM6; Cooke et al. (2022)', lw = 4)
plt.fill_between(o2_conc, A_SZA_60, A_SZA_45, alpha = alpha, color = 'darkorange', label = 'Atmos 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc_less, kasting_60sza, kasting_45sza, alpha = alpha, color = 'teal', label = 'Kasting 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.plot(o2_conc, WACCM, color = 'black')
plt.xscale('log')
thick_axes(top = True)
#plt.title('Simulations with Cooke et al. (2022) WACCM6 boundary conditions', fontsize = 15, weight = 'bold')
'''
if (P_test == True and trapz == True):
    plt.title('P method, trapz', fontsize = 15, weight = 'bold')
    save = '_P_trapz'
if (P_test == True and trapz == False):
    plt.title('P method, sum', fontsize = 15, weight = 'bold')
    save = '_P_sum'
if (P_test == False and trapz == True):
    plt.title('Z method, trapz', fontsize = 15, weight = 'bold')
    save = '_Z_trapz'
if (P_test == False and trapz == False):
    plt.title('Z method, sum', fontsize = 15, weight = 'bold')
    save = '_Z_sum'
'''
plt.title('')
plt.xlim(1e-3, 1.5)
markersize = 7
plt.ylim(0,380)
plt.plot([0.001, 1], [18, 330], marker = 's', ls = '', color = 'teal', markersize = markersize, label = 'Kasting 1D model; Ji et al. (2024)')
plt.plot([1],[292], color = 'black', marker = 'o', markersize = markersize)
plt.legend(loc = 0, fontsize = 15, frameon = False)
plt.xlabel('Oxygen concentration [PAL]', fontsize = 15, weight = 'bold')
plt.ylabel('O'+sub(3)+' column [DU]', fontsize = 15, weight = 'bold')
#plt.savefig('/Users/gregcooke/python_output/Ozone_column_vs_o2_curve'+save+'.png', dpi = 200, bbox_inches = 'tight')
plt.savefig('/Users/gregcooke/python_output/Ozone_column_vs_o2_curve.png', dpi = 200, bbox_inches = 'tight')
#%%

# 1. Create a common grid for the 1D models to find the min/max envelope
o2_common = [1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 0.5, 1, 1.5] # Grid from 1e-3 to ~1.5 PAL

# Interpolate all 1D models onto the common grid
photo60_interp = np.interp(o2_common, o2_conc, photochem_60)
photo45_interp = np.interp(o2_common, o2_conc, photochem_45)
vulcan45_interp = np.interp(o2_common, o2_conc, V_SZA_45)
vulcan60_interp = np.interp(o2_common, o2_conc, V_SZA_60)
atmos45_interp = np.interp(o2_common, o2_conc, A_SZA_45)
atmos60_interp = np.interp(o2_common, o2_conc, A_SZA_60)
kasting45_interp = np.interp(o2_common, o2_conc_less, kasting_45sza)
kasting60_interp = np.interp(o2_common, o2_conc_less, kasting_60sza)

# 2. Calculate the collective 1D envelope
all_1d_stack = np.vstack([photo45_interp, photo60_interp,  
                          vulcan45_interp, vulcan60_interp, 
                          atmos45_interp, atmos60_interp, 
                          kasting45_interp, kasting60_interp])

one_d_min = np.min(all_1d_stack, axis=0)
one_d_max = np.max(all_1d_stack, axis=0)

# --- Updated Plotting Function ---

# --- Execute Grid ---
plt.figure(figsize=(14, 18))
gs = gridspec.GridSpec(3, 4, hspace=0.3, wspace=0.3)
gs.update(hspace = 0.18)
gs.update(wspace = 0.18)

# Row 0: Unified Comparison
ax0 = plt.subplot(gs[0, :])
ax0.set_xlabel('Oxygen mixing ratio [PAL]', fontsize=15, weight='bold')
plot_panel(ax0, highlight='unified_1d', show_legend=True)
ax0.set_title("O"+sub(2)+"-O"+sub(3)+" curve between 1D models and WACCM6", fontsize=15, weight='bold')
# Row 1 & 2: Breakdown
ax1 = plt.subplot(gs[1, 0:2]); plot_panel(ax1, highlight='kasting', show_legend=True, show_ylabel=False)
ax2 = plt.subplot(gs[1, 2:4]); plot_panel(ax2, highlight='photochem', show_legend=True, show_ylabel=False)
thick_axes(top=True, labelleft = False)
ax3 = plt.subplot(gs[2, 0:2]); plot_panel(ax3, highlight='atmos', show_legend=True)
ax4 = plt.subplot(gs[2, 2:4]); plot_panel(ax4, highlight='vulcan', show_legend=True, show_ylabel=False)
ax3.set_xlabel('Oxygen mixing ratio [PAL]', fontsize=15, weight='bold')
ax4.set_xlabel('Oxygen mixing ratio [PAL]', fontsize=15, weight='bold')
thick_axes(top=True, labelleft = False)
plt.savefig('/Users/gregcooke/python_output/Ozone_column_vs_o2_curve.png', dpi = 200, bbox_inches = 'tight')


#%% Big panel money plot


# --- Figure Creation ---
plt.figure(figsize=(14, 18))
gs = gridspec.GridSpec(3, 4, hspace=0.3, wspace=0.3)

# 1. Top Panel: THE ALL-MODEL RANGE
ax0 = plt.subplot(gs[0, :])
plot_panel(ax0, highlight='all', show_legend=True)
ax0.set_title("Comparison: All 1D Models vs WACCM6", fontsize=16, weight='bold')

# 2. Second Row: Highlights
ax1 = plt.subplot(gs[1, 0:2]); plot_panel(ax1, highlight='kasting', show_legend=True)
ax1.set_title("WACCM6 vs Kasting 1D")

ax2 = plt.subplot(gs[1, 2:4]); plot_panel(ax2, highlight='photochem', show_legend=True, show_ylabel=False)
ax2.set_title("WACCM6 vs Photochem")

# 3. Third Row: Highlights
ax3 = plt.subplot(gs[2, 0:2]); plot_panel(ax3, highlight='atmos', show_legend=True)
ax3.set_title("WACCM6 vs Atmos")

ax4 = plt.subplot(gs[2, 2:4]); plot_panel(ax4, highlight='vulcan', show_legend=True, show_ylabel=False)
ax4.set_title("WACCM6 vs VULCAN")



#%%


plt.figure()
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('O3')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'black', ls = '-', label = 'O'+sub(3))
plt.plot(np.flip(LWAV(Pre_pc_h0.O3)), np.flip(Pre_pc_h0.lev/1e3), lw = lw, color = 'black', ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('H2O')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'blue', ls = '-', label = 'O'+sub(3))
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('CO2')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'green', ls = '-', label = 'O'+sub(3))
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('O2')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'grey', ls = '-', label = 'O'+sub(3))
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('N2O')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'purple', ls = '-', label = 'O'+sub(3))
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('CH4')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'orange', ls = '-', label = 'O'+sub(3))
plt.xscale('log'); plt.xlim(1e-9,1)
plt.text(1e-5, 1e-4, str(int(Earth_V_58SZA_O3_col)) + ' DU')
plt.yscale('log'); plt.ylim(1,1e-7)

plt.figure()
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('O3')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'black', ls = '-', label = 'O'+sub(3))
plt.plot(Earth_V_45SZA['variable']['ymix'][:,Earth_V_spec.index('O3')], Earth_V_45SZA['atm']['pco']/1e6, lw = lw, color = 'black', ls = ':', label = 'O'+sub(3))
plt.plot(np.flip(LWAV(Pre_pc_h0.O3)), np.flip(Pre_pc_h0.lev/1e3), lw = lw, color = 'black', ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V_10pc['variable']['ymix'][:,Earth_V_spec_10pc.index('O3')], Earth_V_10pc['atm']['pco']/1e6, lw = lw, color = Ten_pc_color, ls = '-', label = 'O'+sub(3))
plt.plot(np.flip(LWAV(Ten_pc_h0.O3)), np.flip(Ten_pc_h0.lev/1e3), lw = lw, color = Ten_pc_color, ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V_1pc['variable']['ymix'][:,Earth_V_spec_1pc.index('O3')], Earth_V_1pc['atm']['pco']/1e6, lw = lw, color = One_pc_color, ls = '-', label = 'O'+sub(3))
plt.plot(np.flip(LWAV(One_pc_h0.O3)), np.flip(One_pc_h0.lev/1e3), lw = lw, color = One_pc_color, ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V_01pc['variable']['ymix'][:,Earth_V_spec_01pc.index('O3')], Earth_V_01pc['atm']['pco']/1e6, lw = lw, color = Zero1_pc_color, ls = '-', label = 'O'+sub(3))
plt.plot(np.flip(LWAV(Zero1_pc_h0.O3)), np.flip(Zero1_pc_h0.lev/1e3), lw = lw, color = Zero1_pc_color, ls = '--', label = 'O'+sub(3))

plt.xscale('log'); plt.xlim(1e-9,1e-4)
plt.text(1e-5, 1e-4, str(int(Earth_V_58SZA_O3_col)) + ' DU')
plt.text(1e-5, 3e-4, str(int(LWAV(Pre_col))) + ' DU')
plt.text(1e-5, 1e-3, str(int(Earth_V_58SZA_10pc_O3_col)) + ' DU', color = Ten_pc_color)
plt.text(1e-5, 3e-3, str(int(LWAV(Ten_col))) + ' DU', color = Ten_pc_color)
plt.text(1e-5, 1e-2, str(int(Earth_V_58SZA_1pc_O3_col)) + ' DU', color = One_pc_color)
plt.text(1e-5, 3e-2, str(int(LWAV(One_col))) + ' DU', color = One_pc_color)
plt.text(1e-5, 1e-1, str(int(Earth_V_58SZA_01pc_O3_col)) + ' DU', color = Zero1_pc_color)
plt.text(1e-5, 3e-1, str(int(LWAV(Zero1_col))) + ' DU', color = Zero1_pc_color)
plt.yscale('log'); plt.ylim(1,1e-7)

#%%


One50_pc_V_482SZA_O3_col = V_O3_col_z_trapz(One50_pc_V_482SZA)
PI_V_482SZA_O3_col = V_O3_col_z_trapz(PI_V_482SZA)
Fifty_pc_V_482SZA_O3_col = V_O3_col_z_trapz(Fifty_pc_V_482SZA)
Ten_pc_V_482SZA_O3_col = V_O3_col_z_trapz(Ten_pc_V_482SZA)
Ten_pc_V_482SZA_O3_col_T = V_O3_col_z_trapz(Ten_pc_V_482SZA_T)
One_pc_V_482SZA_O3_col = V_O3_col_z_trapz(One_pc_V_482SZA)
One_pc_V_482SZA_O3_col_T = V_O3_col_z_trapz(One_pc_V_482SZA_T)
Five_pc_V_482SZA_O3_col = V_O3_col_z_trapz(Five_pc_V_482SZA)
Z5_pc_V_482SZA_O3_col = V_O3_col_z_trapz(Z5_pc_V_482SZA)
Z1_pc_V_482SZA_O3_col = V_O3_col_z_trapz(Z1_pc_V_482SZA)
Z1_pc_V_482SZA_O3_col_T = V_O3_col_z_trapz(Z1_pc_V_482SZA_T)

One50_pc_V_45SZA_O3_col = V_O3_col_z_trapz(One50_pc_V_45SZA)
Fifty_pc_V_45SZA_O3_col = V_O3_col_z_trapz(Fifty_pc_V_45SZA)
Ten_pc_V_45SZA_O3_col = V_O3_col_z_trapz(Ten_pc_V_45SZA)
One_pc_V_45SZA_O3_col = V_O3_col_z_trapz(One_pc_V_45SZA)
Five_pc_V_45SZA_O3_col = V_O3_col_z_trapz(Five_pc_V_45SZA)
Z5_pc_V_45SZA_O3_col = V_O3_col_z_trapz(Z5_pc_V_45SZA)
Z1_pc_V_45SZA_O3_col = V_O3_col_z_trapz(Z1_pc_V_45SZA)

One50_pc_V_60SZA_O3_col = V_O3_col_z_trapz(One50_pc_V_60SZA)
Fifty_pc_V_60SZA_O3_col = V_O3_col_z_trapz(Fifty_pc_V_60SZA)
Ten_pc_V_60SZA_O3_col = V_O3_col_z_trapz(Ten_pc_V_60SZA)
One_pc_V_60SZA_O3_col = V_O3_col_z_trapz(One_pc_V_60SZA)
Five_pc_V_60SZA_O3_col = V_O3_col_z_trapz(Five_pc_V_60SZA)
Z5_pc_V_60SZA_O3_col = V_O3_col_z_trapz(Z5_pc_V_60SZA)
Z1_pc_V_60SZA_O3_col = V_O3_col_z_trapz(Z1_pc_V_60SZA)

markersize = 5

figure = plt.figure(figsize = (8,5))
#WACCM results
plt.plot([1.5, 1.5, 1.5], [One50_col.min(), LWAV(One50_col), One50_col.max()], marker = 's', color = 'crimson', markersize = markersize, lw = 0, label = 'WACCM6')
plt.plot([1, 1, 1], [Pre_col.min(), LWAV(Pre_col), Pre_col.max()], marker = 's', color = 'crimson', markersize = markersize, lw = 0)
plt.plot([0.5, 0.5, 0.5], [Fifty_col.min(), LWAV(Fifty_col), Fifty_col.max()], marker = 's', color = 'crimson', markersize = markersize, lw = 0)
plt.plot([0.1, 0.1, 0.1], [Ten_col.min(), LWAV(Ten_col), Ten_col.max()], marker = 's', color = 'crimson', markersize = markersize, lw = 0)
plt.plot([0.05, 0.05, 0.05], [Five_col.min(), LWAV(Five_col), Five_col.max()], marker = 's', color = 'crimson', markersize = markersize, lw = 0)
plt.plot([0.01, 0.01, 0.01], [One_col.min(), LWAV(One_col), One_col.max()], marker = 's', color = 'crimson', markersize = markersize, lw = 0)
plt.plot([5e-3, 5e-3, 5e-3], [Zero5_col.min(), LWAV(Zero5_col), Zero5_col.max()], marker = 's', color = 'crimson', markersize = markersize, lw = 0)
plt.plot([1e-3, 1e-3, 1e-3], [Zero1_col.min(), LWAV(Zero1_col), Zero1_col.max()], marker = 's', color = 'crimson', markersize = markersize, lw = 0)

#VULCAN results 48 degree solar zenith angle
plt.plot([1.5], One50_pc_V_482SZA_O3_col, marker = 'o', color = 'black', markersize = markersize, lw = 0, label = 'VULCAN SZA = 48.2'+r'$^\circ$')
plt.plot([1], PI_V_482SZA_O3_col, marker = 's', color = 'black', markersize = markersize, lw = 2)
plt.plot([0.5], Fifty_pc_V_482SZA_O3_col, marker = 'o', color = 'black', markersize = markersize, lw = 2)
plt.plot([0.1], Ten_pc_V_482SZA_O3_col, marker = 'o', color = 'black', markersize = markersize, lw = 2)
plt.plot([0.05], Five_pc_V_482SZA_O3_col, marker = 'o', color = 'black', markersize = markersize, lw = 2)
plt.plot([0.01], One_pc_V_482SZA_O3_col, marker = 'o', color = 'black', markersize = markersize, lw = 2)
plt.plot([0.005], Z5_pc_V_482SZA_O3_col, marker = 'o', color = 'black', markersize = markersize, lw = 2)
plt.plot([0.001], Z1_pc_V_482SZA_O3_col, marker = 'o', color = 'black', markersize = markersize, lw = 2)

plt.legend(loc = 0, fontsize = 15, frameon = False)

plt.xscale('log')
plt.xlim(5e-4, 2)
plt.ylabel('Ozone column depth [DU]', fontsize = 15, weight = 'bold')
plt.xlabel('Oxygen mixing ratio [PAL]', fontsize = 15, weight = 'bold')
thick_axes()
plt.savefig('/Users/gregcooke/Downloads/VULCAN_O2_O3_test.png', dpi = 400)

#%% Proxima Centauri stuff WACCM
path = '/Users/gregcooke/H_escape/'

File = path+"b.e21.BWma1850.f19_g17.PC_b.SSPO.016.cam.h0.0320-0320.nc"
PCb = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.PC_b.SSPL.012.cam.h0.0143-0160.nc"
PCb_SPL = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.PC_b.10pc_o2.001.cam.h0.0328.nc"
PCb_10pc = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.PC_b.1pc_o2.002.cam.h0.0338.nc"
PCb_1pc = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
File = path+"b.e21.BWma1850.f19_g17.PC_b.0.1pc_o2.001.cam.h0.0320.nc"
PCb_01pc = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

PCb_W_O3_col = O3_col(PCb, time = False, O2_mr = 0.21, lon = True, g = 9.81)
PCb_W_10pc_O3_col = O3_col(PCb_10pc, time = False, O2_mr = 0.21, lon = True, g = 9.81)
PCb_W_1pc_O3_col = O3_col(PCb_1pc, time = False, O2_mr = 0.21, lon = True, g = 9.81)
PCb_W_01pc_O3_col = O3_col(PCb_01pc, time = False, O2_mr = 0.21, lon = True, g = 9.81)

PCB_W_cols = [LWAV(PCb_W_01pc_O3_col), LWAV(PCb_W_1pc_O3_col), 
              LWAV(PCb_W_10pc_O3_col), LWAV(PCb_W_O3_col)]

PCb_V_482SZA_O3_col = V_O3_col_z_trapz(PCb_V)
PCb_50pc_V_482SZA_O3_col = V_O3_col_z_trapz(PCb_50pc_V)
PCb_10pc_V_482SZA_O3_col = V_O3_col_z_trapz(PCb_10pc_V)
PCb_1pc_V_482SZA_O3_col = V_O3_col_z_trapz(PCb_1pc_V)
#PCb_01pc_V_482SZA_O3_col = V_O3_col_z_trapz(PCb_01pc_V)
PCB_V_cols = [PCb_1pc_V_482SZA_O3_col, PCb_10pc_V_482SZA_O3_col, PCb_V_482SZA_O3_col]


path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/Proxima_Centauri/"
PCb_photo_data = pd.read_csv(path+'Proxima_100pc_48.2.txt', delim_whitespace=True)
PCb_10pc_photo_data = pd.read_csv(path+'Proxima_10pc_48.2.txt', delim_whitespace=True)
PCb_1pc_photo_data = pd.read_csv(path+'Proxima_1pc_48.2.txt', delim_whitespace=True)
PCb_01pc_photo_data = pd.read_csv(path+'Proxima_0.1pc_48.2.txt', delim_whitespace=True)

Photo_PCb_col = Photochem_O3_col(PCb_photo_data, g=9.81, MO2 = 0.21)
Photo_PCb_10pc_col = Photochem_O3_col(PCb_10pc_photo_data, g=9.81, MO2 = 0.21)
Photo_PCb_1pc_col = Photochem_O3_col(PCb_1pc_photo_data, g=9.81, MO2 = 0.21)
Photo_PCb_01pc_col = Photochem_O3_col(PCb_01pc_photo_data, g=9.81, MO2 = 0.21)

PCB_P_cols = [Photo_PCb_01pc_col, Photo_PCb_1pc_col, Photo_PCb_10pc_col, Photo_PCb_col]

markersize = 9
plt.figure(figsize = (9,5))
plt.plot([0.001, 0.01, 0.1, 1], PCB_W_cols, color = 'black', marker = 's', markersize = markersize, ls = '', label = 'WACCM6')
plt.plot([0.01, 0.1, 1], PCB_V_cols, color = 'm', marker = 'o', markersize = markersize, ls = '', label = 'VULCAN')
plt.plot([0.001, 0.01, 0.1, 1], PCB_P_cols, color = 'b', marker = 'v', markersize = markersize, ls = '', label = 'Photochem')
plt.xscale('log')
plt.yscale('log')
thick_axes(top = True)
plt.legend(loc=(1.0,0.76), fontsize = 15,
           frameon = False, ncol = 1)
plt.title('Proxima Centauri b simulations', fontsize = 15, weight = 'bold')
plt.ylabel('O'+sub(3)+' column [DU]', fontsize = 15, weight = 'bold')
plt.xlabel('O'+sub(2)+' mixing ratio [PAL]', fontsize = 15, weight = 'bold')
plt.savefig('/Users/gregcooke/python_output/PCb_ozone_cols.png', bbox_inches = 'tight', dpi = 400)

#%% Create new VULCAN files

DS = Fifty_pc_h0
# --------------------------------------------------
# Input / output files
# --------------------------------------------------
infile = "/Users/gregcooke/VIH_cases/atm/atm_Earth_Jan_Kzz.txt"
outfile = "/Users/gregcooke/VIH_cases/atm/atm_Earth_PCb_0.1pc_o2_WACCM_Kzz.txt"
#outfile = "/Users/gregcooke/Downloads/atm_Earth_PCb_PI_WACCM_Kzz.txt"
outfile = "/Users/gregcooke/VIH_cases/atm/atm_Earth_50pc_o2_WACCM_Kzz_Z3.txt"

# --------------------------------------------------
# Load original Kzz profile
# --------------------------------------------------
data = np.loadtxt(
    infile,
    skiprows=3
)


P_old = data[:, 0]    # dyne/cm^2
Kzz_old = data[:, 2]  # cm^2/s

# --------------------------------------------------
# Load WACCM pressure & temperature
# --------------------------------------------------
P_new = np.flip(np.array(DS.lev/1e3) )    # dyne/cm^2
T_new = np.flip(np.array(LWAV(DS.T)) )     # K
Z3_new = np.flip(np.array(LWAV(DS.Z3)) )     # K

# --------------------------------------------------
# Ensure monotonic ordering for interpolation
# (np.interp requires increasing x)
# --------------------------------------------------
sort_idx = np.argsort(P_old)
P_old_sorted = P_old[sort_idx]
Kzz_old_sorted = Kzz_old[sort_idx]

# --------------------------------------------------
# Log–log interpolation of Kzz
# --------------------------------------------------
logP_old = np.log10(P_old_sorted)
logKzz_old = np.log10(Kzz_old_sorted)

logP_new = np.log10(P_new)

logKzz_new = np.interp(
    logP_new,
    logP_old,
    logKzz_old,
    left=logKzz_old[0],
    right=logKzz_old[-1]
)

Kzz_new = 10.0**logKzz_new

# --------------------------------------------------
# Write output file
# --------------------------------------------------
header = (
    "# (bar)           (K)      (cm2/s)     (m)\n"
    "# Pressure        Temp     Kzz     Alt\n"
)

with open(outfile, "w") as f:
    f.write(header)
    for p, t, k, z in zip(P_new, T_new, Kzz_new, Z3_new):
        f.write(f"{p:12.5E} {t:8.2f} {k:10.3E} {z:10.3E}\n")

print(f"Wrote {outfile}")


#%% Proxima Centauri stuff VULCAN

'''
Should do TRAPPIST-1 e as well?
'''


PCb_V_45SZA, PCb_spec_V_45SZA = Read_O3_Run(file_path='PCb_1e12s_45SZA_WPT_1rtol.vul')
PCb_V_60SZA, PCb_spec_V_60SZA = Read_O3_Run(file_path='PCb_1e12s_60SZA_WBC_WPT_1rtol.vul')
PCb_10pc_V_45SZA, PCb_10pc_spec_V_45SZA = Read_O3_Run(file_path='PCb_10pc_o2_1e12s_45SZA_WPT_1rtol.vul')
PCb_10pc_V_60SZA, PCb_10pc_spec_V_60SZA = Read_O3_Run(file_path='PCb_10pc_o2_1e12s_60SZA_WBC_WPT_1rtol.vul')
PCb_1pc_V_45SZA, PCb_1pc_spec_V_45SZA = Read_O3_Run(file_path='PCb_1pc_o2_1e12s_45SZA_WPT_1rtol.vul')
PCb_1pc_V_60SZA, PCb_1pc_spec_V_60SZA = Read_O3_Run(file_path='PCb_1pc_o2_1e12s_60SZA_WBC_WPT_1rtol.vul')
PCb_01pc_V_45SZA, PCb_01pc_spec_V_45SZA = Read_O3_Run(file_path='PCb_01pc_o2_1e12s_45SZA_WPT_1rtol.vul')
PCb_01pc_V_60SZA, PCb_01pc_spec_V_60SZA = Read_O3_Run(file_path='PCb_0.1pc_o2_1e12s_60SZA_WBC_WPT_1rtol.vul')


PCb_V_45SZA_O3_col = V_O3_col_z_trapz(PCb_V_45SZA)
PCb_V_60SZA_O3_col = V_O3_col_z_trapz(PCb_V_60SZA)
PCb_10pc_V_45SZA_O3_col = V_O3_col_z_trapz(PCb_10pc_V_45SZA)
PCb_10pc_V_60SZA_O3_col = V_O3_col_z_trapz(PCb_10pc_V_60SZA)
PCb_1pc_V_45SZA_O3_col = V_O3_col_z_trapz(PCb_1pc_V_45SZA)
PCb_1pc_V_60SZA_O3_col = V_O3_col_z_trapz(PCb_1pc_V_60SZA)
PCb_01pc_V_45SZA_O3_col = V_O3_col_z_trapz(PCb_01pc_V_45SZA)
PCb_01pc_V_60SZA_O3_col = V_O3_col_z_trapz(PCb_01pc_V_60SZA)

plt.plot()
ls = '-'
#species = data['variable']['species']
ls = ':'
plt.plot(PCb_01pc_V_45SZA['variable']['ymix'][:,PCb_01pc_spec_V_45SZA.index('O3')], PCb_01pc_V_45SZA['atm']['pco']/1e6, lw = lw, color = Zero1_pc_color, ls = ls, label = 'O'+sub(3))
plt.plot(PCb_1pc_V_45SZA['variable']['ymix'][:,PCb_1pc_spec_V_45SZA.index('O3')], PCb_1pc_V_45SZA['atm']['pco']/1e6, lw = lw, color = One_pc_color, ls = ls, label = 'O'+sub(3))
plt.plot(PCb_10pc_V_45SZA['variable']['ymix'][:,PCb_10pc_spec_V_45SZA.index('O3')], PCb_10pc_V_45SZA['atm']['pco']/1e6, lw = lw, color = Ten_pc_color, ls = ls, label = 'O'+sub(3))
plt.plot(PCb_V_45SZA['variable']['ymix'][:,PCb_spec_V_45SZA.index('O3')], PCb_1pc_V_45SZA['atm']['pco']/1e6, lw = lw, color = 'black', ls = ls, label = 'O'+sub(3))
plt.xlim(1e-9,1e-4)
plt.xscale('log')
plt.yscale('log')
plt.ylim(1,1e-6)

plt.figure(figsize = (9,5))
markersize = 10; capsize = 7; lw = 3
M_cols_TL = 387
plt.errorbar([1], M_cols_TL, yerr=[[M_cols_TL-262],[1466-M_cols_TL]], capsize = capsize,  lw = lw, marker = 's', markersize = markersize, color = 'g')
M_cols_SOR32 = 731
plt.errorbar([1], M_cols_SOR32, yerr=[[M_cols_SOR32-613],[878-M_cols_SOR32]], capsize = capsize,  lw = lw, marker = '*', markersize = markersize, color = 'b')
M3D_cols_TL = 20
plt.errorbar([0.01], M3D_cols_TL, yerr=[[M3D_cols_TL-0],[63-M3D_cols_TL]], capsize = capsize,  lw = lw, marker = 's', markersize = markersize, color = 'g')
W_cols = np.array([LWAV(PCb_W_01pc_O3_col).values, LWAV(PCb_W_1pc_O3_col).values, LWAV(PCb_W_10pc_O3_col).values, LWAV(PCb_W_O3_col).values])
W_min_cols = np.array([PCb_W_01pc_O3_col.min().values, PCb_W_1pc_O3_col.min().values, PCb_W_10pc_O3_col.min().values, PCb_W_O3_col.min().values])
W_max_cols = np.array([PCb_W_01pc_O3_col.max().values, PCb_W_1pc_O3_col.max().values, PCb_W_10pc_O3_col.max().values, PCb_W_O3_col.max().values])
plt.errorbar([0.001, 0.01, 0.1, 1], W_cols, yerr=[W_cols-W_min_cols,W_max_cols-W_cols], capsize = capsize,  lw = lw, marker = 's', markersize = markersize, color = 'black')
V_cols = [PCb_01pc_V_45SZA_O3_col, PCb_1pc_V_45SZA_O3_col, PCb_10pc_V_45SZA_O3_col, PCb_V_45SZA_O3_col]
plt.plot([0.001, 0.01, 0.1, 1], V_cols, color = 'magenta', lw = lw, marker = 'o', markersize = markersize)
V_cols = [PCb_01pc_V_60SZA_O3_col, PCb_1pc_V_60SZA_O3_col, PCb_10pc_V_60SZA_O3_col, PCb_V_60SZA_O3_col]
plt.plot([0.001, 0.01, 0.1, 1], V_cols, color = 'pink', lw = lw, marker = 'o', markersize = markersize)
plt.ylim(5,2000)
plt.yscale('log')
plt.xscale('log')
thick_axes(labelleft=True, top = True, right = True)

#%% 
GJ551 = pd.read_csv('/Users/gregcooke/MUSCLES/sflux-GJ551_0.0_albedo.txt',
                       delim_whitespace = True, skiprows = 1, names = ['Wav','Flux'])


hc = 1.98644582E-9 # Planck constant times the light speed (erg nm)
au = 1.49597871E13  # Astronomical Unit (cm)
r_sun = 6.957E10 # solar radius (cm)
r_jup = 7.1492E9 # Jupiter equatorial radius (cm)
ag0 = 0 # the asymmetry factor in RT (0 for isotropic scattering)
r_star = 0.1542 # stellar radius in solar radius
Rp = 6.681E8 # Planetary radius (cm) (for computing gravity)
orbit_radius = 0.04856  # planet-star distance in A.U.

plt.figure()
flux = GJ551['Flux']*((r_star*r_sun)/(au*orbit_radius))**2
plt.plot(GJ551['Wav'], flux, color = 'black')

x = np.trapz(flux, GJ551['Wav'])/1000

solar_dir = '/Users/gregcooke/stellar_files/'
solar_file = 'SolarForcingCMIP6piControl_c160921.nc'#read in file
ds = xr.open_dataset(solar_dir+solar_file)#attach file to dataset
ssi = ds['ssi'].isel(time=0) #define dataset from file

WACCM_Flux = ssi.values
WACCM_wavelength = ssi.wavelength.values

plt.plot(WACCM_wavelength, WACCM_Flux)

y = np.trapz(WACCM_Flux, WACCM_wavelength)/1000

plt.yscale('log')
plt.xscale('log')
plt.xlim(1,10000)

print(x)
print(y)

GJ551 = pd.read_csv('/Users/gregcooke/MUSCLES/sflux-GJ551_0.0_albedo.txt',
                       delim_whitespace = True, skiprows = 1, names = ['Wav','Flux'])
GJ551_W = pd.read_csv('/Users/gregcooke/MUSCLES/PCb_WACCM_0.0_albedo.txt',
                       delim_whitespace = True, skiprows = 1, names = ['Wav','Flux'])


hc = 1.98644582E-9 # Planck constant times the light speed (erg nm)
au = 1.49597871E13  # Astronomical Unit (cm)
r_sun = 6.957E10 # solar radius (cm)
r_jup = 7.1492E9 # Jupiter equatorial radius (cm)
ag0 = 0 # the asymmetry factor in RT (0 for isotropic scattering)
r_star = 0.1542 # stellar radius in solar radius
Rp = 6.681E8 # Planetary radius (cm) (for computing gravity)
orbit_radius = 	0.04848  # planet-star distance in A.U.

plt.figure()

flux = GJ551['Flux']*((r_star*r_sun)/(au*orbit_radius))**2
plt.plot(GJ551['Wav'], flux, color = 'black')
flux = GJ551_W['Flux']*((r_star*r_sun)/(au*orbit_radius))**2
plt.plot(GJ551_W['Wav'], flux, color = 'm')

#x = np.trapz(flux, GJ551['Wav'])/1000

solar_dir = '/Users/gregcooke/stellar_files/'
solar_file = 'SolarForcingCMIP6piControl_c160921.nc'#read in file
ds = xr.open_dataset(solar_dir+solar_file)#attach file to dataset
ssi = ds['ssi'].isel(time=0) #define dataset from file

WACCM_Flux = ssi.values
WACCM_wavelength = ssi.wavelength.values

plt.plot(WACCM_wavelength, WACCM_Flux)

y = np.trapz(WACCM_Flux, WACCM_wavelength)/1000

plt.yscale('log')
plt.xscale('log')
plt.xlim(1,10000)

print(x)
print(y)

#%% Convert .nc TRAPPIST-1 e files to txt files

plt.figure()
solar_dir = '/Users/gregcooke/stellar_files/'
solar_file = 'TRAPPIST1_flux_at_e_P19.nc'#read in file
ds = xr.open_dataset(solar_dir+solar_file)#attach file to dataset
ssi = ds['ssi'].isel(time=0) #define dataset from file
P19_Flux = ssi.values
WACCM_wavelength = ssi.wavelength.values
plt.plot(WACCM_wavelength, WACCM_Flux, color = 'black')
plt.plot(WACCM_wavelength, P19_Flux, color = 'darkorange')
solar_file = 'TRAPPIST1_flux_at_e_W21.nc'#read in file
ds = xr.open_dataset(solar_dir+solar_file)#attach file to dataset
ssi = ds['ssi'].isel(time=0) #define dataset from file
W21_Flux = ssi.values
WACCM_wavelength = ssi.wavelength.values
plt.plot(WACCM_wavelength, W21_Flux, color = 'navy')

flux = GJ551['Flux']*((r_star*r_sun)/(au*orbit_radius))**2
plt.plot(GJ551['Wav'], flux, color = 'black')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1,10000)

#%% Write to files
# WL(nm)         Flux(ergs/cm**2/s/nm)
hc = 1.98644582E-9 # Planck constant times the light speed (erg nm)
au = 1.49597871E13  # Astronomical Unit (cm)
r_sun = 6.957E10 # solar radius (cm)
r_jup = 7.1492E9 # Jupiter equatorial radius (cm)
ag0 = 0 # the asymmetry factor in RT (0 for isotropic scattering)
r_star =	 0.1192 # stellar radius in solar radius
Rp = 6.681E8 # Planetary radius (cm) (for computing gravity)
orbit_radius = 	0.02925

W21_new_flux = (W21_Flux / ((((r_star*r_sun)) / (au*orbit_radius))**2))
new_str = '# WL(nm)\t Flux(ergs/cm**2/s/nm)\n'

for i in range(len(W21_new_flux)):
    new_str += '{:<12}'.format(WACCM_wavelength[i]) + "{:>12.2E}".format(float(W21_new_flux[i])) + '\n'
albedo = 0.0
path = '/Users/gregcooke/MUSCLES/'
out_name = 'TP-1e_W21'
with open(path+out_name+'_'+str(albedo)+'_albedo.txt', 'w+') as f: f.write(new_str)

P19_new_flux = (P19_Flux / ((((r_star*r_sun)) / (au*orbit_radius))**2))
new_str = '# WL(nm)\t Flux(ergs/cm**2/s/nm)\n'

for i in range(len(P19_new_flux)):
    new_str += '{:<12}'.format(WACCM_wavelength[i]) + "{:>12.2E}".format(float(P19_new_flux[i])) + '\n'
albedo = 0.0
path = '/Users/gregcooke/MUSCLES/'
out_name = 'TP-1e_P19'
with open(path+out_name+'_'+str(albedo)+'_albedo.txt', 'w+') as f: f.write(new_str)

solar_dir = '/Users/gregcooke/stellar_files/'
solar_file = 'Proxima_Centauri_at_b.nc'#read in file
ds = xr.open_dataset(solar_dir+solar_file)#attahc file to dataset
ssi = ds['ssi'].isel(time=0) #define dataset from file

'''
PCb file
'''
hc = 1.98644582E-9 # Planck constant times the light speed (erg nm)
au = 1.49597871E13  # Astronomical Unit (cm)
r_sun = 6.957E10 # solar radius (cm)
r_jup = 7.1492E9 # Jupiter equatorial radius (cm)
ag0 = 0 # the asymmetry factor in RT (0 for isotropic scattering)
r_star =	 0.1542 # stellar radius in solar radius
Rp = 6.681E8 # Planetary radius (cm) (for computing gravity)
orbit_radius = 	0.04848

PC_b_flux = ssi.values

df = ssi.to_dataframe(name='Solar flux (mW/m^2/nm)').reset_index()

# 2. Select only the columns you want
# Check if your coordinate is named 'wavelength' or 'wvl' and adjust if necessary
df = df[['wavelength', 'Solar flux (mW/m^2/nm)']]

# 3. Rename columns for the header
df.columns = ['Wavelength (nm)', 'Solar flux (mW/m^2/nm)']

# 4. Save to text file with specific formatting
output_path = '/Users/gregcooke/python_output/Proxima_Centauri_at_b.txt'

# We use sep='\t' or spaces and float_format for that clean look
df.to_csv(output_path, 
          sep='\t', 
          index=False, 
          float_format='%.6e', # Ensures scientific notation like your example
          quoting=None)

print(f"File saved successfully to {output_path}")

PC_b_wavelength = ssi.wavelength.values

PCb_new_flux = (PC_b_flux / ((((r_star*r_sun)) / (au*orbit_radius))**2))
new_str = '# WL(nm)\t Flux(ergs/cm**2/s/nm)\n'

for i in range(len(PCb_new_flux)):
    new_str += '{:<12}'.format(WACCM_wavelength[i]) + "{:>12.2E}".format(float(PCb_new_flux[i])) + '\n'
albedo = 0.0
path = '/Users/gregcooke/MUSCLES/'
out_name = 'PCb_WACCM'
with open(path+out_name+'_'+str(albedo)+'_albedo.txt', 'w+') as f: f.write(new_str)

#%% Is there a photon flux factor that is missing in the conversion? Check your script and Shamis script
GJ551 = pd.read_csv('/Users/gregcooke/MUSCLES/sflux-TP-1e_0.0_albedo.txt',
                       delim_whitespace = True, skiprows = 1, names = ['Wav','Flux'])
flux = GJ551['Flux']
plt.plot(GJ551['Wav'], flux, color = 'black')
GJ551 = pd.read_csv('/Users/gregcooke/MUSCLES/TP-1e_P19_0.0_albedo.txt',
                       delim_whitespace = True, skiprows = 1, names = ['Wav','Flux'])
flux = GJ551['Flux']
plt.plot(GJ551['Wav'], flux, color = 'orange')
GJ551 = pd.read_csv('/Users/gregcooke/MUSCLES/TP-1e_W21_0.0_albedo.txt',
                       delim_whitespace = True, skiprows = 1, names = ['Wav','Flux'])
flux = GJ551['Flux']
plt.plot(GJ551['Wav'], flux, color = 'blue')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1,10000)


#%% Checks for Earth, TRAPPIST-1 e, and TOI-1468 c

GJ551 = pd.read_csv('/Users/gregcooke/VIH_cases/atm/stellar_flux/Gueymard_solar.txt',
                       delim_whitespace = True, skiprows = 1, names = ['Wav','Flux'])
hc = 1.98644582E-9 # Planck constant times the light speed (erg nm)
au = 1.49597871E13  # Astronomical Unit (cm)
r_sun = 6.957E10 # solar radius (cm)
r_jup = 7.1492E9 # Jupiter equatorial radius (cm)
ag0 = 0 # the asymmetry factor in RT (0 for isotropic scattering)
r_star = 0.1542 # stellar radius in solar radius
Rp = 6.681E8 # Planetary radius (cm) (for computing gravity)
orbit_radius = 0.04856  # planet-star distance in A.U.
r_star = 1
orbit_radius = 1
flux = GJ551['Flux']*((r_star*r_sun)/(au*orbit_radius))**2
wav = GJ551['Wav']

print(np.trapz(flux, wav)/1e3)

GJ551 = pd.read_csv('/Users/gregcooke/VIH_cases/atm/stellar_flux/TP-1e_W21_0.0_albedo.txt',
                       delim_whitespace = True, skiprows = 1, names = ['Wav','Flux'])

hc = 1.98644582E-9 # Planck constant times the light speed (erg nm)
au = 1.49597871E13  # Astronomical Unit (cm)
r_sun = 6.957E10 # solar radius (cm)
r_jup = 7.1492E9 # Jupiter equatorial radius (cm)
ag0 = 0 # the asymmetry factor in RT (0 for isotropic scattering)
r_star =	 0.1192 # stellar radius in solar radius
Rp = 6.681E8 # Planetary radius (cm) (for computing gravity)
orbit_radius = 	0.02925
flux = GJ551['Flux']*((r_star*r_sun)/(au*orbit_radius))**2
wav = GJ551['Wav']

print(np.trapz(flux, wav)/1e3)

GJ551 = pd.read_csv('/Users/gregcooke/MUSCLES/sflux-TOI-1468_0.0_albedo.txt',
                       delim_whitespace = True, skiprows = 1, names = ['Wav','Flux'])
'''
GJ551['Flux'] = GJ551['Flux']/ 1.3190196360217292

GJ551.to_csv('/Users/gregcooke/MUSCLES/sflux-TOI-1468_0.0_albedo.txt',
               index = False,        sep = '\t')
'''
hc = 1.98644582E-9 # Planck constant times the light speed (erg nm)
au = 1.49597871E13  # Astronomical Unit (cm)
r_sun = 6.957E10 # solar radius (cm)
r_jup = 7.1492E9 # Jupiter equatorial radius (cm)
ag0 = 0 # the asymmetry factor in RT (0 for isotropic scattering)
r_star =	 0.3714 # stellar radius in solar radius
Rp = 13.54E8 # Planetary radius (cm) (for computing gravity)
orbit_radius = 0.08608
flux = GJ551['Flux']*((r_star*r_sun)/(au*orbit_radius))**2
wav = GJ551['Wav']

print(np.trapz(flux, wav)/1e3)

GJ551['Flux'] = flux

GJ551.to_csv('/Users/gregcooke/photochem/examples/TOI1468c/sflux-TOI-1468_0.0_albedo.txt',
               index = False,        sep = '\t')


#%% O2 photolysis Earth

File = "/Users/gregcooke/H_escape/b.e21.BWma1850.f19_g17.baseline.H2O_SRB_updated_xs.001.cam.h1.0001-01.nc" #file name
base = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
jo2_base= base.jo2_a + base.jo2_b
ox_prod_base = (base.jo2_a + base.jo2_b)*n_dens(base, 'O2')
File = "/Users/gregcooke/H_escape/b.e21.BWma1850.f19_g17.10pc_o2.LBC_fixed_fluxes_SRB_update_XS.001.cam.h1.0057-0060.nc" #file name
Ten_pc_new = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
jo2_ten = Ten_pc_new.jo2_a + Ten_pc_new.jo2_b
ox_prod_ten = (Ten_pc_new.jo2_a + Ten_pc_new.jo2_b)*n_dens(Ten_pc_new, 'O2')
File = "/Users/gregcooke/H_escape/b.e21.BWma1850.f19_g17.1pc_o2.LBC_fixed_fluxes.006.cam.h1.0053-0056.nc" #file name
One_pc_new = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
jo2_one = One_pc_new.jo2_a + One_pc_new.jo2_b
ox_prod_one = (One_pc_new.jo2_a + One_pc_new.jo2_b)*n_dens(One_pc_new, 'O2')
File = "/Users/gregcooke/H_escape/b.e21.BWma1850.f19_g17.0.1pc_o2.LBC_fixed_fluxes_SRB_update_XS.004.cam.h1.0080-0083.nc" #file name
Zero1_pc_new = xr.open_dataset(File,decode_times=False) #open the file and decode time as false
jo2_z1 = Zero1_pc_new.jo2_a + Zero1_pc_new.jo2_b
ox_prod_z1 = (Zero1_pc_new.jo2_a + Zero1_pc_new.jo2_b)*n_dens(Zero1_pc_new, 'O2')


int_jo2_base = int_O2_photo(base, O2 = 0.21, g = 9.81)
int_jo2_ten = int_O2_photo(Ten_pc_new, O2 = 0.021, g = 9.81)
int_jo2_one = int_O2_photo(One_pc_new, O2 = 0.021, g = 9.81)
int_jo2_z1 = int_O2_photo(Zero1_pc_new, O2 = 0.021, g = 9.81)
plt.figure()
int_jo2_base.plot()
int_jo2_ten.plot()
int_jo2_one.plot()
int_jo2_z1.plot()
#plt.yscale('log')
plt.ylim(6e15, 0.35e17)

plt.figure(figsize = (8,6))
ax1 = plt.gca()
idx = LWAV(ox_prod_ten).argmax()
ox_prod_ten.isel(lev = idx).mean(dim = 'lon').plot(color=Ten_pc_color, lw = lw)
plt.axhline(LWAV(ox_prod_ten).max(), color=Ten_pc_color, ls = ':')
idx = LWAV(ox_prod_one).argmax()
ox_prod_one.isel(lev = idx).mean(dim = 'lon').plot(color=One_pc_color, lw = lw)
plt.axhline(LWAV(ox_prod_one).max(), color=One_pc_color, ls = ':')
idx = LWAV(ox_prod_z1).argmax()
ox_prod_z1.isel(lev = idx).mean(dim = 'lon').plot(color=Zero1_pc_color, lw = lw)
plt.axhline(LWAV(ox_prod_z1).max(), color=Zero1_pc_color, ls = ':')
plt.xlim(-90,90)
plt.title('Ox production versus ozone column with latitude', fontsize = 15, weight = 'bold')
plt.xlabel('Latitude ['+sup('\circ')+']', fontsize = 15, weight = 'bold')
plt.ylabel('Max Ox production [molecules m-3 s-1]', fontsize = 15, weight = 'bold')
thick_axes(top = True, labelleft = True)
plt.axvline(39);plt.axvline(-40)

ax2 = ax1.twinx()
Ten_col.mean(dim = 'lon').plot(color=Ten_pc_color, lw = lw, ls = '--')
One_col.mean(dim = 'lon').plot(color=One_pc_color, lw = lw, ls = '--')
Zero1_col.mean(dim = 'lon').plot(color=Zero1_pc_color, lw = lw, ls = '--')
thick_axes(top = True, labelleft = False)

plt.figure(figsize = (18,12))
gs = gridspec.GridSpec(2,2)
plt.subplot(gs[0,1])
ax1 = plt.gca()
idx = LWAV(ox_prod_ten).argmax()
ox_prod_ten.isel(lev = idx).mean(dim = 'lon').plot(color=Ten_pc_color, lw = lw)
plt.axhline(LWAV(ox_prod_ten).max(), color=Ten_pc_color, ls = ':')
plt.xlim(-90,90)
plt.suptitle('Ox production versus ozone column with latitude', fontsize = 15, weight = 'bold', y = 0.95)
plt.xlabel('Latitude ['+sup('\circ')+']', fontsize = 15, weight = 'bold')
plt.ylabel('Max Ox production [molecules m-3 s-1]', fontsize = 15, weight = 'bold')
thick_axes(top = True, labelleft = True)
plt.axvline(39);plt.axvline(-40)

ax2 = ax1.twinx()
Ten_col.mean(dim = 'lon').plot(color=Ten_pc_color, lw = lw, ls = '--')
plt.axhline(LWAV(Ten_col), color=Ten_pc_color, ls = 'dashdot')
thick_axes(top = True, labelleft = False)

plt.subplot(gs[1,0])
ax1 = plt.gca()
idx = LWAV(ox_prod_one).argmax()
ox_prod_one.isel(lev = idx).mean(dim = 'lon').plot(color=One_pc_color, lw = lw)
plt.axhline(LWAV(ox_prod_one).max(), color=One_pc_color, ls = ':')
plt.xlim(-90,90)
plt.xlabel('Latitude ['+sup('\circ')+']', fontsize = 15, weight = 'bold')
plt.ylabel('Max Ox production [molecules m-3 s-1]', fontsize = 15, weight = 'bold')
thick_axes(top = True, labelleft = True)
plt.axvline(39);plt.axvline(-40)

ax2 = ax1.twinx()
One_col.mean(dim = 'lon').plot(color=One_pc_color, lw = lw, ls = '--')
plt.axhline(LWAV(One_col), color=One_pc_color, ls = 'dashdot')
thick_axes(top = True, labelleft = False)

plt.subplot(gs[1,1])
ax1 = plt.gca()
idx = LWAV(ox_prod_z1).argmax()
ox_prod_z1.isel(lev = idx).mean(dim = 'lon').plot(color=Zero1_pc_color, lw = lw)
plt.axhline(LWAV(ox_prod_z1).max(), color=Zero1_pc_color, ls = ':')
plt.xlim(-90,90)
plt.xlabel('Latitude ['+sup('\circ')+']', fontsize = 15, weight = 'bold')
plt.ylabel('Max Ox production [molecules m-3 s-1]', fontsize = 15, weight = 'bold')
thick_axes(top = True, labelleft = True)
plt.axvline(39);plt.axvline(-40)

ax2 = ax1.twinx()
Zero1_col.mean(dim = 'lon').plot(color=Zero1_pc_color, lw = lw, ls = '--')
plt.axhline(LWAV(Zero1_col), color=Zero1_pc_color, ls = 'dashdot')
thick_axes(top = True, labelleft = False)

'''

Now plot the integrated jo2 at different latitudes. I think you have a 
function for it in a different file
in
'''

#%% O2 photolysis Earth

path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/Proxima_Centauri/"
Photo_PCb_100pc_txt = pd.read_csv(path+"Proxima_100pc_48.2.txt", delim_whitespace=True)
Photo_PCb_10pc_txt = pd.read_csv(path+"Proxima_10pc_48.2.txt", delim_whitespace=True)
Photo_PCb_1pc_txt = pd.read_csv(path+"Proxima_1pc_48.2.txt", delim_whitespace=True)
Photo_PCb_01pc_txt = pd.read_csv(path+"Proxima_0.1pc_48.2.txt", delim_whitespace=True)

plt.figure(figsize = (14,10))
plt.plot(LWAV(PCb.O3), PCb.lev, color = 'black')
plt.plot(Photo_PCb_100pc_txt['O3'], Photo_PCb_100pc_txt['press']*1e3, ls = ':', color = 'black', lw = lw)
O3 = PCb_V['variable']['ymix'][:,PCb_V_spec.index('O3')]
plt.plot(O3, PCb_V['atm']['pco']/1e3, ls = '--', color = Ten_pc_color, lw = lw)

plt.plot(LWAV(PCb_10pc.O3), PCb.lev, color = Ten_pc_color, lw = lw)
O3 = PCb_10pc_V['variable']['ymix'][:,PCb_V_spec.index('O3')]
plt.plot(O3, PCb_10pc_V['atm']['pco']/1e3, ls = '--', color = Ten_pc_color, lw = lw)
plt.plot(Photo_PCb_10pc_txt['O3'], Photo_PCb_10pc_txt['press']*1e3, ls = ':', color = Ten_pc_color, lw = lw)

plt.plot(LWAV(PCb_1pc.O3), PCb.lev, color = One_pc_color, lw = lw)
O3 = PCb_1pc_V['variable']['ymix'][:,PCb_V_spec.index('O3')]
plt.plot(O3, PCb_1pc_V['atm']['pco']/1e3, ls = '--', color = Ten_pc_color, lw = lw)
plt.plot(Photo_PCb_1pc_txt['O3'], Photo_PCb_1pc_txt['press']*1e3, ls = ':', color = One_pc_color, lw = lw)

plt.plot(LWAV(PCb_01pc.O3), PCb.lev, color = Zero1_pc_color, lw = lw)
plt.plot(Photo_PCb_01pc_txt['O3'], Photo_PCb_01pc_txt['press']*1e3, ls = ':', color = Zero1_pc_color, lw = lw)
plt.yscale('log')
plt.xscale('log');  thick_axes(top = True);
plt.xlim(1e-10, 1e-4)
plt.ylim(1e3, 1e-5)

ylim = (1e3, 1e-3)
xlim = (1e11, 6e12)

PCb_photo, p_PCb_photo = compute_ox_production(
    mech_file=P_path+"Old sims/zahnle_earth.yaml",
    settings_file=P_path+"Proxima_Centauri/input/settings_100pc.yaml",
    flux_file=P_path+"Proxima_Centauri/input/Proxima_Centauri.txt",
    pt_file=P_path+"Proxima_Centauri/Proxima_100pc_48.2.txt",
    atol=1e-23,
    verbose=0
)

PCb_photo_10pc, p_PCb_photo_10pc = compute_ox_production(
    mech_file=P_path+"Old sims/zahnle_earth.yaml",
    settings_file=P_path+"Proxima_Centauri/input/settings_10pc.yaml",
    flux_file=P_path+"Proxima_Centauri/input/Proxima_Centauri.txt",
    pt_file=P_path+"Proxima_Centauri/Proxima_10pc_48.2.txt",
    atol=1e-23,
    verbose=0
)

PCb_photo_1pc, p_PCb_photo_1pc = compute_ox_production(
    mech_file=P_path+"Old sims/zahnle_earth.yaml",
    settings_file=P_path+"Proxima_Centauri/input/settings_1pc.yaml",
    flux_file=P_path+"Proxima_Centauri/input/Proxima_Centauri.txt",
    pt_file=P_path+"Proxima_Centauri/Proxima_1pc_48.2.txt",
    atol=1e-23,
    verbose=0
)

PCb_photo_01pc, p_PCb_photo_01pc = compute_ox_production(
    mech_file=P_path+"Old sims/zahnle_earth.yaml",
    settings_file=P_path+"Proxima_Centauri/input/settings_0.1pc.yaml",
    flux_file=P_path+"Proxima_Centauri/input/Proxima_Centauri.txt",
    pt_file=P_path+"Proxima_Centauri/Proxima_0.1pc_48.2.txt",
    atol=1e-23,
    verbose=0
)

plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)
gs.update(wspace = 0.2, hspace = 0.1)
plt.subplot(gs[0,0])
plt.title('100% PAL  ', fontsize = 15, weight = 'bold', y = 0.9, loc = 'right')
plt.plot(LWAV(prox_ox_W(PCb)), base.lev, color = 'black', label = 'WACCM6', lw = 2)
plt.plot(PCb_photo, p_PCb_photo/10, color = 'b', label = 'Photochem', lw = 2)
plt.plot(prox_ox_V(PCb_V, PCb_V_spec)*3/8, PCb_V['atm']['pco']/1e3, color = 'm', lw = 2, label = 'VULCAN',)
plt.yscale('log')
plt.xscale('log');  thick_axes(top = True);
plt.xlim(1e10, 5e12)
plt.ylim(1e3, 1e-5)
plt.legend(loc = (0), handlelength = 0.5, ncol=1, 
           fontsize = 15, frameon = False, columnspacing = 0.5)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[0,1])
plt.title('10% PAL  ', fontsize = 15, weight = 'bold', y = 0.9, loc = 'right')
plt.plot(LWAV(prox_ox_W(PCb_10pc)), base.lev, color = 'black', label = 'WACCM6', lw = 2)
plt.plot(PCb_photo_10pc, p_PCb_photo_10pc/10, color = 'b', label = 'Photochem', lw = 2)
plt.plot(prox_ox_V(PCb_10pc_V, PCb_V_spec)*3/8, PCb_10pc_V['atm']['pco']/1e3, color = 'm', lw = 2, label = 'VULCAN',)
plt.yscale('log')
plt.xscale('log');  thick_axes(top = True);
plt.xlim(1e10, 5e12)
plt.ylim(1e3, 1e-5)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')


plt.subplot(gs[1,0])
plt.title('1% PAL  ', fontsize = 15, weight = 'bold', y = 0.9, loc = 'right')
plt.plot(LWAV(prox_ox_W(PCb_1pc)), base.lev, color = 'black', label = 'WACCM6', lw = 2)
plt.plot(PCb_photo_1pc, p_PCb_photo_1pc/10, color = 'b', label = 'Photochem', lw = 2)
plt.plot(prox_ox_V(PCb_1pc_V, PCb_V_spec)*3/8, PCb_1pc_V['atm']['pco']/1e3, color = 'm', lw = 2, label = 'VULCAN',)
plt.yscale('log')
plt.xscale('log');  thick_axes(top = True);
plt.xlim(1e10, 5e12)
plt.ylim(1e3, 1e-5)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[1,1])
plt.title('0.1% PAL  ', fontsize = 15, weight = 'bold', y = 0.9, loc = 'right')
plt.plot(LWAV(prox_ox_W(PCb_01pc)), base.lev, color = 'black', label = 'WACCM6', lw = 2)
plt.plot(PCb_photo_01pc, p_PCb_photo_01pc/10, color = 'b', label = 'Photochem', lw = 2)
#plt.plot(prox_ox_V(PCb_01pc_V, PCb_V_spec)*3/8, PCb_1pc_V['atm']['pco']/1e3, color = 'm', lw = 2, label = 'VULCAN',)
plt.yscale('log')
plt.xscale('log');  thick_axes(top = True);
plt.xlim(1e10, 5e12)
plt.ylim(1e3, 1e-5)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

ox_prod_100pc_photo, p_100pc_photo = compute_ox_production(
    mech_file=P_path+"Old sims/zahnle_earth.yaml",
    settings_file=P_path+"Old sims/input/settings_100pc.yaml",
    flux_file=P_path+"Old sims/input/Sun_0.0Ga.txt",
    pt_file=P_path+"100pc/Earth_100pc_48.2.txt",
    atol=1e-23,
    verbose=0
)

ox_prod_10pc_photo, p_10pc_photo = compute_ox_production(
    mech_file=P_path+"Old sims/zahnle_earth.yaml",
    settings_file=P_path+"Old sims/input/settings_100pc.yaml",
    flux_file=P_path+"Old sims/input/Sun_0.0Ga.txt",
    pt_file=P_path+"10pc/Earth_10pc_48.2.txt",
    atol=1e-23,
    verbose=0
)

ox_prod_1pc_photo, p_1pc_photo = compute_ox_production(
    mech_file=P_path+"Old sims/zahnle_earth.yaml",
    settings_file=P_path+"Old sims/input/settings_100pc.yaml",
    flux_file=P_path+"Old sims/input/Sun_0.0Ga.txt",
    pt_file=P_path+"1pc/Earth_1pc_48.2.txt",
    atol=1e-23,
    verbose=0
)

ox_prod_01pc_photo, p_01pc_photo = compute_ox_production(
    mech_file=P_path+"Old sims/zahnle_earth.yaml",
    settings_file=P_path+"Old sims/input/settings_100pc.yaml",
    flux_file=P_path+"Old sims/input/Sun_0.0Ga.txt",
    pt_file=P_path+"0.1pc/Earth_0.1pc_48.2.txt",
    atol=1e-23,
    verbose=0
)

Gauss = True

plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)
gs.update(wspace = 0.1, hspace = 0.1)
plt.subplot(gs[0,0])
plt.title('100% PAL  ', fontsize = 15, weight = 'bold', y = 0.9, loc = 'right')
plt.plot(LWAV(prox_ox_W(base)), base.lev, color = 'black', label = 'WACCM6', lw = 2)
if (Gauss == True):
    plt.plot(prox_ox_K(K_100pc_J_8G), p_100pc_8G, color = 'teal', lw = 2, ls = '--')
plt.plot(prox_ox_K(K_100pc_J), p_100pc, color = 'teal', label = 'Kasting', lw = 2)
plt.plot(prox_ox_V(PI_V_482SZA, PI_spec_482SZA)*3/8, PI_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN', lw = 2)
plt.plot(ox_prod_100pc_photo, p_100pc_photo/10, color = 'b', label = 'Photochem', lw = 2)
plt.plot(Atmos_XS_100pc_JO2*Atmos_dens(Atmos_100pc)*Atmos_100pc["O2"], Atmos_100pc["PRESS"]*1e3, color = 'darkorange', lw = 2, label = 'Atmos')
plt.yscale('log')
plt.xscale('log');  thick_axes(top = True);
plt.xlim(xlim)
plt.ylim(ylim)
plt.legend(loc = (0.01, 0.3), handlelength = 0.5, ncol=1, 
           fontsize = 15, frameon = False, columnspacing = 0.5)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[0,1])
plt.title('10% PAL  ', fontsize = 15, weight = 'bold', y = 0.9, loc = 'right')
plt.plot(LWAV(prox_ox_W(Ten_pc_new)), base.lev, color = 'black', lw = 2)
if (Gauss == True):
    plt.plot(prox_ox_K(K_10pc_J_8G), p_10pc_8G, color = 'teal', lw = 2, ls = '--')
plt.plot(prox_ox_K(K_10pc_J), p_10pc, color = 'teal', lw = 2)
plt.plot(prox_ox_V(Ten_pc_V_482SZA, Ten_pc_spec_482SZA)*3/8, Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = 2)
plt.plot(ox_prod_10pc_photo, p_10pc_photo/10, color = 'b', label = 'Photochem', lw = 2)
plt.plot(Atmos_XS_10pc_JO2*Atmos_dens(Atmos_10pc)*Atmos_10pc["O2"], Atmos_10pc["PRESS"]*1e3, color = 'darkorange', lw = 2, label = 'Atmos')
plt.yscale('log')
plt.xscale('log');  thick_axes(top = True, labelleft = False)
plt.xlim(xlim)
plt.ylim(ylim)

plt.subplot(gs[1,0])
plt.title('1% PAL  ', fontsize = 15, weight = 'bold', y = 0.9, loc = 'right')
plt.plot(LWAV( prox_ox_W(One_pc_new)), base.lev, color = 'black', lw = 2)
if (Gauss == True):
    plt.plot(prox_ox_K(K_1pc_J_8G), p_1pc_8G, color = 'teal', lw = 2, ls = '--')
plt.plot(prox_ox_K(K_1pc_J), p_1pc, color = 'teal', lw = 2)
plt.plot(prox_ox_V(One_pc_V_482SZA, One_pc_spec_482SZA)*3/8, One_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = 2)
plt.plot(ox_prod_1pc_photo, p_1pc_photo/10, color = 'b', label = 'Photochem', lw = 2)
plt.plot(Atmos_XS_1pc_JO2*Atmos_dens(Atmos_1pc)*Atmos_1pc["O2"], Atmos_1pc["PRESS"]*1e3, color = 'darkorange', lw = 2, label = 'Atmos')
plt.yscale('log')
plt.xscale('log');  thick_axes(top = True)
plt.xlim(xlim)
plt.ylim(ylim)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.xlabel('O'+sub('x')+' production [molecules m'+sup(-3)+' s'+sup(-1)+']', fontsize = 15, weight = 'bold')

plt.subplot(gs[1,1])
plt.title('0.1% PAL  ', fontsize = 15, weight = 'bold', y = 0.9, loc = 'right')
plt.plot(LWAV( prox_ox_W(Zero1_pc_new)), base.lev, color = 'black', lw = 2)
plt.plot(prox_ox_K(K_01pc_J), p_01pc, color = 'teal', lw = 2)
if (Gauss == True):
    plt.plot(prox_ox_K(K_01pc_J_8G), p_01pc_8G, color = 'teal', lw = 2, ls = '--')
plt.plot(prox_ox_V(Z1_pc_V_482SZA, Z1_pc_spec_482SZA)*3/8, Z1_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = 2)
plt.plot(ox_prod_01pc_photo, p_01pc_photo/10, color = 'b', label = 'Photochem', lw = 2)
plt.plot(Atmos_XS_01pc_JO2*Atmos_dens(Atmos_01pc)*Atmos_01pc["O2"], Atmos_01pc["PRESS"]*1e3, color = 'darkorange', lw = 2, label = 'Atmos')
plt.yscale('log')
plt.xscale('log');  thick_axes(top = True, labelleft = False)
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel('O'+sub('x')+' production [molecules m'+sup(-3)+' s'+sup(-1)+']', fontsize = 15, weight = 'bold')

plt.savefig('/Users/gregcooke/python_output/Ox_production_all_models.png')

#%% Ozone density all models

xlim = (5e16, 6e18)
ylim = (1e3, 1e-1)

plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)
gs.update(wspace = 0.1, hspace = 0.1)
plt.subplot(gs[0,0])
plt.title('100% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(n_dens(Pre_pc_h0, "O3")), Pre_pc_h0.lev, color = 'black', label = 'WACCM6')
plt.plot(PI_V_482SZA['variable']['y'][:,PI_spec_482SZA.index('O3')]*1e6, PI_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN')
plt.plot(Photo_PI["O3"]*P_dens(Photo_PI), Photo_PI["press"]*1e3, color = 'b', label = 'Photochem', lw = 2)
plt.plot(Atmos_100pc["O3"]*Atmos_dens(Atmos_100pc), Atmos_100pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = 2)
plt.plot(K_100pc_J["FO3"]*K_100pc_J["DEN"]*1e6, K_100pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting', lw = 2)
plt.yscale('log')
plt.xlim(xlim)
plt.ylim(ylim)
plt.xscale('log');  thick_axes(top = True)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[0,1])
plt.title('10% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(n_dens(Ten_pc_h0, "O3")), Pre_pc_h0.lev, color = 'black', lw = lw)
plt.plot(Ten_pc_V_482SZA['variable']['y'][:,Ten_pc_spec_482SZA.index('O3')]*1e6, Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(Photo_10pc["O3"]*P_dens(Photo_10pc), Photo_10pc["press"]*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_10pc["O3"]*Atmos_dens(Atmos_10pc), Atmos_10pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_10pc_J["FO3"]*K_10pc_J["DEN"]*1e6, K_10pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting', lw = 2)
plt.yscale('log')
plt.xlim(xlim)
plt.ylim(ylim)
plt.xscale('log');  thick_axes(top = True, labelleft = False)

plt.subplot(gs[1,0])
plt.title('1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(n_dens(One_pc_h0, "O3")), Pre_pc_h0.lev, color = 'black', lw = lw, label = 'WACCM6')
plt.plot(One_pc_V_482SZA['variable']['y'][:,One_pc_spec_482SZA.index('O3')]*1e6, Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN', lw = lw)
plt.plot(Photo_1pc["O3"]*P_dens(Photo_1pc), Photo_1pc["press"]*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_1pc["O3"]*Atmos_dens(Atmos_1pc), Atmos_1pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_1pc_J["FO3"]*K_1pc_J["DEN"]*1e6, K_1pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting', lw = 2)
plt.yscale('log')
plt.xlim(xlim)
plt.ylim(ylim)
plt.legend(loc = (0.17, 0.6), handlelength = 1,ncol = 2, 
           fontsize = 15, frameon = False, columnspacing = 0.5)
plt.xscale('log');  thick_axes(top = True)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.xlabel('O'+sub(3)+' number density [molecules m'+sup(-3)+']', fontsize = 15, weight = 'bold')

plt.subplot(gs[1,1])
plt.title('0.1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(n_dens(Zero1_pc_h0, "O3")), Pre_pc_h0.lev, color = 'black', lw = lw)
plt.plot(Z1_pc_V_482SZA['variable']['y'][:,Z1_pc_spec_482SZA.index('O3')]*1e6, Z1_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(Photo_01pc["O3"]*P_dens(Photo_01pc), Photo_01pc["press"]*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_01pc["O3"]*Atmos_dens(Atmos_01pc), Atmos_01pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_01pc_J["FO3"]*K_01pc_J["DEN"]*1e6, K_01pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting', lw = 2)
plt.yscale('log')
plt.xscale('log');  thick_axes(top = True, labelleft = False)
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel('O'+sub(3)+' number density [molecules m'+sup(-3)+']', fontsize = 15, weight = 'bold')

plt.savefig('/Users/gregcooke/python_output/O3_density_all_models.png')

#%% NOx comparison

plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)
gs.update(wspace = 0.1, hspace = 0.1)
plt.subplot(gs[0,0])
plt.title('100% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Pre_pc_h2.NOX), Pre_pc_h0.lev, color = 'black', label = 'WACCM6', lw = lw)
N = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('N')]
NO = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('NO')]
NO2 = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('NO2')]
NOX = N+NO+NO2
plt.plot(NOX, PI_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN', lw = lw)
plt.plot(P_NOX(Photo_PI), Photo_PI['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_NOX(Atmos_100pc), Atmos_100pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(Kasting_NOX(K_100pc_J), K_100pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-5)
plt.xscale('log');  thick_axes(top = True)
plt.xlim(1e-12, 1e-6)
plt.legend(loc = (0.03, 0.4), handlelength = 1,
           fontsize = 15, frameon = False)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[0,1])
plt.title('10% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Ten_pc_h2.NOX), Pre_pc_h0.lev, color = 'black', lw = lw)
N = Ten_pc_V_482SZA['variable']['ymix'][:,Ten_pc_spec_482SZA.index('N')]
NO = Ten_pc_V_482SZA['variable']['ymix'][:,Ten_pc_spec_482SZA.index('NO')]
NO2 = Ten_pc_V_482SZA['variable']['ymix'][:,Ten_pc_spec_482SZA.index('NO2')]
NOX = N+NO+NO2
plt.plot(Atmos_NOX(Atmos_10pc), Atmos_10pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(NOX, Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(P_NOX(Photo_10pc), Photo_10pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Kasting_NOX(K_10pc_J), K_10pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-5)
plt.xscale('log');  thick_axes(top = True, labelleft = False)
plt.xlim(1e-12, 1e-6)

plt.subplot(gs[1,0])
plt.title('1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(One_pc_h2.NOX), Pre_pc_h0.lev, color = 'black', lw = lw)
N = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('N')]
NO = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('NO')]
NO2 = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('NO2')]
NOX = N+NO+NO2
plt.plot(NOX, One_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(Atmos_NOX(Atmos_1pc), Atmos_1pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(P_NOX(Photo_1pc), Photo_1pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Kasting_NOX(K_1pc_J), K_1pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-5)
plt.xscale('log');  thick_axes(top = True)
plt.xlim(1e-12, 1e-6)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.xlabel('NO'+sub('\\rm x')+' mixing ratio', fontsize = 15, weight = 'bold')

plt.subplot(gs[1,1])
plt.title('0.1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Zero1_pc_h2.NOX), Pre_pc_h0.lev, color = 'black', lw = lw)
N = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('N')]
NO = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('NO')]
NO2 = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('NO2')]
NOX = N+NO+NO2
plt.plot(NOX, Z1_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(P_NOX(Photo_01pc), Photo_01pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_NOX(Atmos_01pc), Atmos_01pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(Kasting_NOX(K_01pc_J), K_01pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-5)
plt.xscale('log');  thick_axes(top = True, labelleft = False)
plt.xlim(1e-12, 1e-6)
plt.xlabel('NO'+sub('\\rm x')+' mixing ratio', fontsize = 15, weight = 'bold')

plt.savefig('/Users/gregcooke/python_output/NOX_all_models.png')

#%% HOX comparison
xlim = (5e-13, 1e-6)
plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)
gs.update(wspace = 0.1, hspace = 0.1)
plt.subplot(gs[0,0])
plt.title('100% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Pre_pc_h2.HOX), Pre_pc_h0.lev, color = 'black', label = 'WACCM6', lw = lw)
H = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('H')]
OH = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('OH')]
HO2 = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('HO2')]
H2O2 = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('H2O2')]
HOX = H+OH+HO2+H2O2
plt.plot(HOX, PI_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN', lw = lw)
plt.plot(P_HOX(Photo_PI), Photo_PI['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_HOX(Atmos_100pc), Atmos_100pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(Kasting_HOX(K_100pc_J), K_100pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-3)
plt.xscale('log');  thick_axes(top = True)
plt.xlim(xlim)
plt.legend(loc = (0.03, 0.4), handlelength = 1,
           fontsize = 15, frameon = False)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[0,1])
plt.title('10% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Ten_pc_h2.HOX), Pre_pc_h0.lev, color = 'black', lw = lw)
H = Ten_pc_V_482SZA['variable']['ymix'][:,Ten_pc_spec_482SZA.index('H')]
OH = Ten_pc_V_482SZA['variable']['ymix'][:,Ten_pc_spec_482SZA.index('OH')]
HO2 = Ten_pc_V_482SZA['variable']['ymix'][:,Ten_pc_spec_482SZA.index('HO2')]
H2O2 = Ten_pc_V_482SZA['variable']['ymix'][:,Ten_pc_spec_482SZA.index('H2O2')]
HOX = H+OH+HO2+H2O2
plt.plot(HOX, Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(P_HOX(Photo_10pc), Photo_10pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_HOX(Atmos_10pc), Atmos_10pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(Kasting_HOX(K_10pc_J), K_10pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Atmos', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-3)
plt.xscale('log');  thick_axes(top = True, labelleft = False)
plt.xlim(xlim)

plt.subplot(gs[1,0])
plt.title('1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(One_pc_h2.HOX), Pre_pc_h0.lev, color = 'black', lw = lw)
H = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('H')]
OH = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('OH')]
HO2 = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('HO2')]
H2O2 = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('H2O2')]
HOX = H+OH+HO2+H2O2
plt.plot(HOX, One_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(P_HOX(Photo_1pc), Photo_1pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_HOX(Atmos_1pc), Atmos_1pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(Kasting_HOX(K_1pc_J), K_1pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-3)
plt.xscale('log');  thick_axes(top = True)
plt.xlim(xlim)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.xlabel('HO'+sub('\\rm x')+' mixing ratio', fontsize = 15, weight = 'bold')

plt.subplot(gs[1,1])
plt.title('0.1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Zero1_pc_h2.HOX), Pre_pc_h0.lev, color = 'black', lw = lw)
H = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('H')]
OH = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('OH')]
HO2 = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('HO2')]
H2O2 = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('H2O2')]
HOX = H+OH+HO2+H2O2
plt.plot(HOX, Z1_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(P_HOX(Photo_01pc), Photo_01pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_HOX(Atmos_01pc), Atmos_01pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(Kasting_HOX(K_01pc_J), K_01pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-3)
plt.xscale('log');  thick_axes(top = True, labelleft = False)
plt.xlim(xlim)
plt.xlabel('HO'+sub('\\rm x')+' mixing ratio', fontsize = 15, weight = 'bold')

plt.savefig('/Users/gregcooke/python_output/HOX_all_models.png')

#%% Water vapour profile comparison
ylim = 1e-4 
plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)
gs.update(wspace = 0.1, hspace = 0.1)
plt.subplot(gs[0,0])
plt.title('100% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Pre_pc_h0.H2O), Pre_pc_h0.lev, color = 'black', label = 'WACCM6', lw = lw)
plt.plot(PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('H2O')], PI_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN')
plt.plot(Photo_PI["H2O"], Photo_PI['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_100pc["H2O"], Atmos_100pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_100pc_J["FH2O"], K_100pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting', lw = lw)
plt.yscale('log')
plt.ylim(1e3, ylim)
plt.xscale('log');  thick_axes(top = True)
plt.xlim(1e-7, 1e-4)
plt.legend(loc = (0.01, 0.01), handlelength = 1,
           fontsize = 15, frameon = False)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[0,1])
plt.title('10% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Ten_pc_h0.H2O), Ten_pc_h0.lev, color = 'black', lw = lw)
plt.plot(Ten_pc_V_482SZA['variable']['ymix'][:,Ten_pc_spec_482SZA.index('H2O')], Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN', lw = lw)
plt.plot(Photo_10pc["H2O"], Photo_10pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_10pc["H2O"], Atmos_10pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_10pc_J["FH2O"], K_10pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting', lw = lw)
plt.yscale('log')
plt.ylim(1e3, ylim)
plt.xscale('log');  thick_axes(top = True, labelleft = False)
plt.xlim(1e-7, 1e-4)

plt.subplot(gs[1,0])
plt.title('1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(One_pc_h0.H2O), One_pc_h0.lev, color = 'black', lw = lw)
plt.plot(One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('H2O')], One_pc_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN', lw = lw)
plt.plot(Photo_1pc["H2O"], Photo_1pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_1pc["H2O"], Atmos_1pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_1pc_J["FH2O"], K_1pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-3)
plt.xscale('log');  thick_axes(top = True)
plt.xlim(1e-7, ylim)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.xlabel('H'+sub(2)+'O mixing ratio', fontsize = 15, weight = 'bold')

plt.subplot(gs[1,1])
plt.title('0.1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Zero1_pc_h0.H2O), Zero1_pc_h0.lev, color = 'black', lw = lw)
plt.plot(Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('H2O')], Z1_pc_V_60SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN', lw = lw)
plt.plot(Photo_01pc["H2O"], Photo_01pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_01pc["H2O"], Atmos_01pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_01pc_J["FH2O"], K_01pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting', lw = lw)
plt.yscale('log')
plt.ylim(1e3, ylim)
plt.xscale('log');  thick_axes(top = True, labelleft = False)
plt.xlim(1e-7, 1e-4)
plt.xlabel('H'+sub(2)+'O mixing ratio', fontsize = 15, weight = 'bold')

plt.savefig('/Users/gregcooke/python_output/H2O_all_models.png')
#%% Integrated Ox production calculation


'''
Need to figure out the factors here. 
Are the photolysis rates outputted depending on the diurnal averaging factor or not?
'''

Vulcan_factor = 2

# --- Compute VULCAN integrated JO2 ---
int_JO2_PI = vulcan_integrated_JO2(PI_V_482SZA)
int_JO2_Ten = vulcan_integrated_JO2(Ten_pc_V_482SZA)
int_JO2_One = vulcan_integrated_JO2(One_pc_V_482SZA)
int_JO2_Zero1 = vulcan_integrated_JO2(Z1_pc_V_482SZA)
# --- Existing cumulative calculations ---
cum_int_base = cum_int_O2_photo(base, O2=0.21, g=9.81)
cum_int_ten = cum_int_O2_photo(Ten_pc_new, O2=0.021, g=9.81)
cum_int_one = cum_int_O2_photo(One_pc_new, O2=0.0021, g=9.81)
cum_int_zero1 = cum_int_O2_photo(Zero1_pc_new, O2=0.00021, g=9.81)
# Photochem
int_ox_prod_100pc_photo = photochem_integrated_JO2(ox_prod_100pc_photo)
int_ox_prod_10pc_photo = photochem_integrated_JO2(ox_prod_10pc_photo)
int_ox_prod_1pc_photo = photochem_integrated_JO2(ox_prod_1pc_photo)
int_ox_prod_01pc_photo = photochem_integrated_JO2(ox_prod_01pc_photo)

# Atmos and Kasting: cumulative column integral from TOA downward.
# profile: molecules m^-3 s^-1, z_cm: altitude grid in cm.
def cumulative_from_toa(profile, z_cm):
    arr = np.asarray(profile, dtype=float)
    z_m = np.asarray(z_cm, dtype=float) / 100.0
    dz_m = np.abs(np.gradient(z_m))
    layer_prod = arr * dz_m
    return layer_prod[::-1].cumsum()[::-1]

atmos_ox_prod_100pc = Atmos_XS_100pc_JO2 * Atmos_dens(Atmos_100pc) * Atmos_100pc["O2"]
atmos_ox_prod_10pc = Atmos_XS_10pc_JO2 * Atmos_dens(Atmos_10pc) * Atmos_10pc["O2"]
atmos_ox_prod_1pc = Atmos_XS_1pc_JO2 * Atmos_dens(Atmos_1pc) * Atmos_1pc["O2"]
atmos_ox_prod_01pc = Atmos_XS_01pc_JO2 * Atmos_dens(Atmos_01pc) * Atmos_01pc["O2"]

int_ox_prod_100pc_atmos = cumulative_from_toa(atmos_ox_prod_100pc, Atmos_XS_100pc_Z)
int_ox_prod_10pc_atmos = cumulative_from_toa(atmos_ox_prod_10pc, Atmos_XS_10pc_Z)
int_ox_prod_1pc_atmos = cumulative_from_toa(atmos_ox_prod_1pc, Atmos_XS_1pc_Z)
int_ox_prod_01pc_atmos = cumulative_from_toa(atmos_ox_prod_01pc, Atmos_XS_01pc_Z)

kasting_ox_prod_100pc = prox_ox_K(K_100pc_J)
kasting_ox_prod_10pc = prox_ox_K(K_10pc_J)
kasting_ox_prod_1pc = prox_ox_K(K_1pc_J)
kasting_ox_prod_01pc = prox_ox_K(K_01pc_J)

int_ox_prod_100pc_kasting = cumulative_from_toa(kasting_ox_prod_100pc, K_100pc_J["Z"])
int_ox_prod_10pc_kasting = cumulative_from_toa(kasting_ox_prod_10pc, K_10pc_J["Z"])
int_ox_prod_1pc_kasting = cumulative_from_toa(kasting_ox_prod_1pc, K_1pc_J["Z"])
int_ox_prod_01pc_kasting = cumulative_from_toa(kasting_ox_prod_01pc, K_01pc_J["Z"])

# --- Plot ---
plt.figure(figsize=(14,10))
gs = gridspec.GridSpec(2,2)

xlim = (5e15, 2e17)

plt.subplot(gs[0,0])
# WACCM6
plt.plot(cum_int_base, base.lev, color='k', label='WACCM6', lw = lw)
# VULCAN case
plt.plot(int_JO2_PI, PI_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw, label='VULCAN')
plt.plot(int_ox_prod_100pc_photo, Photo_PI['press']*1e3, color = 'b', lw = lw, label='Photochem')
plt.plot(int_ox_prod_100pc_atmos, Atmos_100pc["PRESS"]*1e3, color='darkorange', lw=lw, label='Atmos')
plt.plot(int_ox_prod_100pc_kasting, K_100pc_J["PRESS"]/1e3, color='teal', lw=lw, label='Kasting')
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.legend(loc = 0, fontsize = 15, frameon = False)
thick_axes(top = True)
plt.yscale('log'); plt.ylim(1e3, 1e-5)
plt.xscale('log'); plt.xlim(xlim)

plt.subplot(gs[0,1])
plt.plot(int_JO2_Ten, Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(int_ox_prod_10pc_photo, Photo_10pc['press']*1e3, color = 'b', lw = lw, label='Photochem')
plt.plot(cum_int_ten, base.lev, color = 'black', lw = lw)
plt.plot(int_ox_prod_10pc_atmos, Atmos_10pc["PRESS"]*1e3, color='darkorange', lw=lw)
plt.plot(int_ox_prod_10pc_kasting, K_10pc_J["PRESS"]/1e3, color='teal', lw=lw)
thick_axes(top = True)
plt.yscale('log'); plt.ylim(1e3, 1e-5)
plt.xscale('log');  plt.xlim(xlim)

plt.subplot(gs[1,0])
plt.plot(int_JO2_One, One_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(cum_int_one, base.lev, color = 'black', lw = lw)
plt.plot(int_ox_prod_1pc_photo, Photo_1pc['press']*1e3, color = 'b', lw = lw, label='Photochem')
plt.plot(int_ox_prod_1pc_atmos, Atmos_1pc["PRESS"]*1e3, color='darkorange', lw=lw)
plt.plot(int_ox_prod_1pc_kasting, K_1pc_J["PRESS"]/1e3, color='teal', lw=lw)
thick_axes(top = True)
plt.yscale('log'); plt.ylim(1e3, 1e-5)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.xscale('log');  plt.xlim(xlim)
plt.xlabel('Integrated O$_2$ photolysis', fontsize = 15, weight = 'bold')

plt.subplot(gs[1,1])
plt.plot(int_JO2_Zero1, Z1_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(int_ox_prod_01pc_photo, Photo_01pc['press']*1e3, color = 'b', lw = lw, label='Photochem')
plt.plot(cum_int_zero1, base.lev, color = 'black', lw = lw)
plt.plot(int_ox_prod_01pc_atmos, Atmos_01pc["PRESS"]*1e3, color='darkorange', lw=lw)
plt.plot(int_ox_prod_01pc_kasting, K_01pc_J["PRESS"]/1e3, color='teal', lw=lw)

thick_axes(top = True)
plt.yscale('log'); plt.ylim(1e3, 1e-5)
plt.xscale('log'); plt.xlim(xlim)
plt.xlabel('Integrated O$_2$ photolysis', fontsize = 15, weight = 'bold')

#thick_axes(top = True)
plt.tight_layout()
plt.savefig('/Users/gregcooke/python_output/Integrated_Ox_Prod.png', dpi = 400, bbox_inches = 'tight')
#%% Other photolysis pots

'''
There is a factor of 10 from somewhere in VULCAN?
VUCLAN multipied by 2x2 because 2 O atoms produed and
photolysis rates are divided by 2 due to diurnal factor
But shouldn't that be taken into account on the nightside?'
'''


base_o2_dens = n_dens(base, var = 'O2')
jo2 = base.jo2_a + base.jo2_b
PI_photo = 2*LWAV(base_o2_dens  * jo2)

Ten_pc_o2_dens = n_dens(Ten_pc_new, var = 'O2')
jo2 = Ten_pc_new.jo2_a + Ten_pc_new.jo2_b
Ten_photo = LWAV(2* Ten_pc_o2_dens  * jo2)

One_pc_o2_dens = n_dens(One_pc_new, var = 'O2')
jo2 = One_pc_new.jo2_a + One_pc_new.jo2_b
One_photo = LWAV(2* One_pc_o2_dens  * jo2)

Zero1_pc_o2_dens = n_dens(Zero1_pc_new, var = 'O2')
jo2 = Zero1_pc_new.jo2_a + Zero1_pc_new.jo2_b
Zero1_photo = LWAV(2* Zero1_pc_o2_dens  * jo2)

'''
Begin figure comparing photolysis rates
'''

ylim = 1e-5
plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)
plt.subplot(gs[0, 0])

plt.plot(PI_photo, base.lev, lw = lw, color = 'black')

lw = 2
V_PI_photo = (
    PI_pc_V_60SZA['variable']['J_sp']['O2',0] +
    PI_pc_V_60SZA['variable']['J_sp']['O2',1] +
    PI_pc_V_60SZA['variable']['J_sp']['O2',2]
)
species = PI_pc_V_60SZA['variable']['species']
O2_dens_PI = PI_pc_V_60SZA['variable']['y'][:,species.index('O2')] * 1e6
V_PI_prod_ox = Vulcan_factor*V_PI_photo*O2_dens_PI
plt.plot(V_PI_prod_ox, PI_pc_V_60SZA['atm']['pco']/1e3, lw = lw, ls = ':', color = 'black')
thick_axes()
plt.ylim(1e3,ylim)
plt.xlim(1e10,1e13)
plt.yscale('log')
plt.xscale('log')

plt.subplot(gs[0, 1])
plt.plot(Ten_photo, base.lev, lw = lw, color = Ten_pc_color)
V_10pc_photo = Ten_pc_V_60SZA['variable']['J_sp']['O2', 0] + Ten_pc_V_60SZA['variable']['J_sp']['O2', 1] + Ten_pc_V_60SZA['variable']['J_sp']['O2', 2]
species = Ten_pc_V_60SZA['variable']['species']
O2_dens_10pc = Ten_pc_V_60SZA['variable']['y'][:,species.index('O2')] * 1e6
plt.plot(Vulcan_factor*V_10pc_photo*O2_dens_10pc, Ten_pc_V_60SZA['atm']['pco']/1e3, lw = lw, ls = ':', color = Ten_pc_color)
thick_axes()
plt.ylim(1e3,ylim)
plt.xlim(1e10,1e13)
plt.yscale('log')
plt.xscale('log')

plt.subplot(gs[1, 0])
plt.plot(One_photo, base.lev, lw = lw, color = One_pc_color)
V_1pc_photo = One_pc_V_60SZA['variable']['J_sp']['O2', 0] + One_pc_V_60SZA['variable']['J_sp']['O2', 1] + One_pc_V_60SZA['variable']['J_sp']['O2', 2]
species = One_pc_V_60SZA['variable']['species']
O2_dens_1pc = One_pc_V_60SZA['variable']['y'][:,species.index('O2')] * 1e6
plt.plot(Vulcan_factor*V_1pc_photo*O2_dens_1pc, One_pc_V_60SZA['atm']['pco']/1e3, lw = lw, ls = ':', color = One_pc_color)
thick_axes()
plt.ylim(1e3,ylim)
plt.xlim(1e10,1e13)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Production of O'+sub('\\rm x')+' [molecules m'+sup(-3)+' s'+sup(-1)+']', fontsize = 15, weight = 'bold')

plt.subplot(gs[1, 1])
plt.plot(Zero1_photo, base.lev, lw = lw, color = Zero1_pc_color)
V_01pc_photo = Z1_pc_V_60SZA['variable']['J_sp']['O2', 0] + Z1_pc_V_60SZA['variable']['J_sp']['O2', 1] + Z1_pc_V_60SZA['variable']['J_sp']['O2', 2]
species = Z1_pc_V_60SZA['variable']['species']
O2_dens_01pc = Z1_pc_V_60SZA['variable']['y'][:,species.index('O2')] * 1e6
V_z1_prod_ox = Vulcan_factor*V_01pc_photo*O2_dens_01pc
plt.plot(V_z1_prod_ox, Z1_pc_V_60SZA['atm']['pco']/1e3, lw = lw, ls = ':', color = Zero1_pc_color)

plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.xlabel('Production of O'+sub('\\rm x')+' [molecules m'+sup(-3)+' s'+sup(-1)+']', fontsize = 15, weight = 'bold')

thick_axes()
plt.ylim(1e3,ylim)
plt.xlim(1e10,1e13)
plt.yscale('log')
plt.xscale('log')
plt.savefig('/Users/gregcooke/python_output/VULCAN_Ox_Prod.png', dpi = 400, bbox_inches = 'tight')


'''
Interolate photolysis rates onto the same pressure grid so I can
divide one by the other
'''

'''
Interpolate photolysis rates onto the same pressure grid so I can
divide one by the other
'''


# Target pressure grid (CAM)
p_cam = base.lev.values

# --- PI case ---
p_v = PI_pc_V_60SZA['atm']['pco']/1e3
vulcan_PI_interp = interp1d(p_v, Vulcan_factor*V_PI_photo*O2_dens_PI,
                            bounds_error=False, fill_value=np.nan)

V_PI_on_cam = vulcan_PI_interp(p_cam)
ratio_PI = PI_photo / V_PI_on_cam


# --- 10% case ---
p_v = Ten_pc_V_60SZA['atm']['pco']/1e3
species = Ten_pc_V_60SZA['variable']['species']
O2_dens_10pc = Ten_pc_V_60SZA['variable']['y'][:,species.index('O2')] * 1e5
V_10_prod = Vulcan_factor * V_10pc_photo * O2_dens_10pc

vulcan_10_interp = interp1d(p_v, V_10_prod,
                           bounds_error=False, fill_value=np.nan)

V_10_on_cam = vulcan_10_interp(p_cam)
ratio_10 = Ten_photo / V_10_on_cam


# --- 1% case ---
p_v = One_pc_V_60SZA['atm']['pco']/1e3
species = One_pc_V_60SZA['variable']['species']
O2_dens = One_pc_V_60SZA['variable']['y'][:,species.index('O2')] * 1e5
V_1_prod = Vulcan_factor * V_1pc_photo * O2_dens

vulcan_1_interp = interp1d(p_v, V_1_prod,
                          bounds_error=False, fill_value=np.nan)

V_1_on_cam = vulcan_1_interp(p_cam)
ratio_1 = One_photo / V_1_on_cam


# --- 0.1% case ---
p_v = Z1_pc_V_60SZA['atm']['pco']/1e3
species = Z1_pc_V_60SZA['variable']['species']
O2_dens = Z1_pc_V_60SZA['variable']['y'][:,species.index('O2')] * 1e5
V_01_prod = Vulcan_factor * V_01pc_photo * O2_dens

vulcan_01_interp = interp1d(p_v, V_01_prod,
                           bounds_error=False, fill_value=np.nan)

V_01_on_cam = vulcan_01_interp(p_cam)
ratio_01 = Zero1_photo / V_01_on_cam

'''
Cumulative photolysis integrated

'''


plt.figure(figsize=(6,5))
plt.plot(ratio_PI, p_cam, color='k')
plt.plot(ratio_10, p_cam, color=Ten_pc_color)
plt.plot(ratio_1, p_cam, color=One_pc_color)
plt.plot(ratio_01, p_cam, color=Zero1_pc_color)

plt.yscale('log')
plt.xscale('log')
plt.ylim(1e3,1e-2)
plt.xlim(0.1, 50)

plt.axvline(1, ls = '--', color = 'm')

plt.xlabel('CAM / VULCAN Ox production', fontsize = 15, weight = 'bold')
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
thick_axes()

#%% O2 photolysis and temperature

plt.figure(figsize = (14,5))
gs = gridspec.GridSpec(1, 3)
plt.subplot(gs[0,0])
prox_ox_W(base).mean(dim='lon').plot(levels = 11)
plt.ylim(100, 0.1)
plt.yscale('log')

plt.subplot(gs[0,1])
n_dens(base, "O").mean(dim='lon').plot(norm = colors.LogNorm(vmax = 5e15, vmin = 1e14))
plt.ylim(100, 0.1)
plt.yscale('log')

plt.subplot(gs[0,2])
base.T.mean(dim='lon').plot(vmin = 200, vmax = 270, cmap = 'coolwarm', levels = 11)
plt.ylim(100, 0.1)
plt.yscale('log')

plt.figure(figsize = (14,7))
gs = gridspec.GridSpec(1, 2)
plt.subplot(gs[0,0])
prox_ox_W(Ten_pc_new).mean(dim='lon').plot(levels = 11)
plt.ylim(100, 0.1)
plt.yscale('log')

plt.subplot(gs[0,1])
Ten_pc_new.T.mean(dim='lon').plot(vmin = 170, vmax = 230, cmap = 'coolwarm', levels = 11)
plt.ylim(100, 0.1)
plt.yscale('log')

plt.figure(figsize = (14,7))
gs = gridspec.GridSpec(1, 2)
plt.subplot(gs[0,0])
prox_ox_W(One_pc_new).mean(dim='lon').plot(levels = 11)
plt.ylim(100, 0.1)
plt.yscale('log')

plt.subplot(gs[0,1])
One_pc_new.T.mean(dim='lon').plot(vmin = 170, vmax = 230, cmap = 'coolwarm', levels = 11)
plt.ylim(100, 0.1)
plt.yscale('log')
#%% O2 photolysis PCb
PCb_o2_dens = n_dens(PCb, var = 'O2')
jo2 = PCb.jo2_a + PCb.jo2_b
PI_PCb_photo = LWAV(2* PCb_o2_dens  * jo2)

PCb_01pc_o2_dens = n_dens(PCb_01pc, var = 'O2')
jo2 = PCb.jo2_a + PCb.jo2_b
PCb_01pc_photo = LWAV(2* PCb_01pc_o2_dens  * jo2)

plt.figure(figsize = (7,5))
plt.plot(LWAV(jo2), base.lev, lw = lw, color = 'black')
thick_axes()
plt.ylim(1e3,1e-5)
plt.xscale('log')
plt.yscale('log')
#plt.xlim(1e10,1e13)

plt.figure(figsize = (7,5))

plt.plot(PI_PCb_photo , base.lev, lw = lw, color = 'black')
#plt.plot(Ten_photo, base.lev, lw = lw, color = Ten_pc_color)
#plt.plot(One_photo, base.lev, lw = lw, color = One_pc_color)
plt.plot(Zero1_photo, base.lev, lw = lw, color = Zero1_pc_color)

plt.title('PCb simulations', fontsize = 15, weight = 'bold')
lw = 2
V_PCb_PI_photo = PCb_V_60SZA['variable']['J_sp']['O2', 0] + PCb_V_60SZA['variable']['J_sp']['O2', 1] + PCb_V_60SZA['variable']['J_sp']['O2', 2]
species = PCb_V_60SZA['variable']['species']
O2_dens = PCb_V_60SZA['variable']['y'][:,species.index('O2')] * 1e5
plt.plot(Vulcan_factor*V_PCb_PI_photo*O2_dens, PCb_V_60SZA['atm']['pco']/1e3, lw = lw, ls = ':', color = 'black')

V_PCb_01pc_photo = PCb_01pc_V_60SZA['variable']['J_sp']['O2', 0] + PCb_01pc_V_60SZA['variable']['J_sp']['O2', 1] + PCb_01pc_V_60SZA['variable']['J_sp']['O2', 2]
species = PCb_01pc_V_60SZA['variable']['species']
O2_dens = PCb_01pc_V_60SZA['variable']['y'][:,species.index('O2')] * 1e5
plt.plot(Vulcan_factor*V_PCb_01pc_photo*O2_dens, PCb_01pc_V_60SZA['atm']['pco']/1e3, lw = lw, ls = ':', color = Zero1_pc_color)

'''
V_10pc_photo = Ten_pc_V_60SZA['variable']['J_sp']['O2', 0] + Ten_pc_V_60SZA['variable']['J_sp']['O2', 1] + Ten_pc_V_60SZA['variable']['J_sp']['O2', 2]
species = Ten_pc_V_60SZA['variable']['species']
O2_dens = Ten_pc_V_60SZA['variable']['y'][:,species.index('O2')] * 1e5
plt.plot(Vulcan_factor*V_10pc_photo*O2_dens, Ten_pc_V_60SZA['atm']['pco']/1e3, lw = lw, ls = ':', color = Ten_pc_color)

V_1pc_photo = One_pc_V_60SZA['variable']['J_sp']['O2', 0] + One_pc_V_60SZA['variable']['J_sp']['O2', 1] + One_pc_V_60SZA['variable']['J_sp']['O2', 2]
species = One_pc_V_60SZA['variable']['species']
O2_dens = One_pc_V_60SZA['variable']['y'][:,species.index('O2')] * 1e5
plt.plot(Vulcan_factor*V_1pc_photo*O2_dens, One_pc_V_60SZA['atm']['pco']/1e3, lw = lw, ls = ':', color = One_pc_color)

V_01pc_photo = Z1_pc_V_60SZA['variable']['J_sp']['O2', 0] + Z1_pc_V_60SZA['variable']['J_sp']['O2', 1] + Z1_pc_V_60SZA['variable']['J_sp']['O2', 2]
species = Z1_pc_V_60SZA['variable']['species']
O2_dens = Z1_pc_V_60SZA['variable']['y'][:,species.index('O2')] * 1e5
plt.plot(Vulcan_factor*V_01pc_photo*O2_dens, Z1_pc_V_60SZA['atm']['pco']/1e3, lw = lw, ls = ':', color = Zero1_pc_color)
'''
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.xlabel('Production of O'+sub('\\rm x')+' [molecules m'+sup(-3)+' s'+sup(-1)+']', fontsize = 15, weight = 'bold')

thick_axes()
plt.ylim(1e3,1e-5)
plt.xlim(1e10,1e13)
plt.yscale('log')
plt.xscale('log')
plt.savefig('/Users/gregcooke/python_output/VULCAN_Ox_Prod_PCb.png', dpi = 400, bbox_inches = 'tight')

#%% Hycean

plt.figure(figsize = (8,5))
TOI, TOI_spec = Read_O3_Run(file_path='TOI_1468c_M50_10bar_1e12s.vul')

H2 = TOI['variable']['ymix'][:,TOI_spec.index('H2')]
CH4 = TOI['variable']['ymix'][:,TOI_spec.index('CH4')]
NH3 = TOI['variable']['ymix'][:,TOI_spec.index('NH3')]
CO = TOI['variable']['ymix'][:,TOI_spec.index('CO')]
CO2 = TOI['variable']['ymix'][:,TOI_spec.index('CO2')]

plt.plot(H2, TOI['atm']['pco']/1e6, lw = lw)
plt.plot(CH4, TOI['atm']['pco']/1e6, lw = lw)
plt.plot(NH3, TOI['atm']['pco']/1e6, lw = lw)
plt.plot(CO, TOI['atm']['pco']/1e6, lw = lw)
plt.plot(CO2, TOI['atm']['pco']/1e6, lw = lw)
plt.ylim(10,1e-8)
plt.yscale('log')
plt.xscale('log')
plt.xlim(1e-10,1)
thick_axes()








#%% PCb JO2 and WACCM6 JO2

'''
Should include integrated photolysis rate too!
'''

xlim=(1e9, 1e13)
ylim=(1e3, 1e-5)
'''
plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)

plt.subplot(gs[0, 0])
PCb_W_JO2 = prox_ox_W(PCb)
PCb_V_JO2 = prox_ox_V(PCb_V, PCb_spec_V)
plt.plot(LWAV(PCb_W_JO2), PCb.lev, lw = lw, color = 'black', label = 'WACCM6')
plt.plot(PCb_V_JO2, PCb_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)


plt.subplot(gs[0, 1])
PCb_W_JO2 = prox_ox_W(PCb_10pc)
#PCb_V_JO2 = prox_ox_V(PCb_10pc_V, PCb_10pc_spec_V)
plt.plot(LWAV(PCb_W_JO2), PCb.lev, lw = lw, color = 'black', label = 'WACCM6')
plt.plot(PCb_V_JO2, PCb_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)


plt.subplot(gs[1, 0])
PCb_W_JO2 = prox_ox_W(PCb_1pc)
#PCb_V_JO2 = prox_ox_V(PCb_1pc_V, PCb_1pc_spec_V)
plt.plot(LWAV(PCb_W_JO2), PCb.lev, lw = lw, color = 'black', label = 'WACCM6')
plt.plot(PCb_V_JO2, PCb_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)

plt.subplot(gs[1, 1])
PCb_W_JO2 = prox_ox_W(PCb_01pc)
#PCb_V_JO2 = prox_ox_V(PCb_01pc_V, PCb_01pc_spec_V)
plt.plot(LWAV(PCb_W_JO2), PCb.lev, lw = lw, color = 'black', label = 'WACCM6')
plt.plot(PCb_V_JO2, PCb_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)
'''

'''
Now plot integrated photolysis rate
'''

int_JO2_PI = vulcan_integrated_JO2(PCb_V)
int_JO2_Ten = vulcan_integrated_JO2(PCb_10pc_V)
int_JO2_One = vulcan_integrated_JO2(PCb_1pc_V)
#int_JO2_Zero1 = vulcan_integrated_JO2(PCb_01pc_V)
# --- Existing cumulative calculations ---
cum_int_base = cum_int_O2_photo(PCb, O2=0.21, g=9.81)
cum_int_ten = cum_int_O2_photo(PCb_10pc, O2=0.021, g=9.81)
cum_int_one = cum_int_O2_photo(PCb_1pc, O2=0.0021, g=9.81)
#cum_int_zero1 = cum_int_O2_photo(PCb_01pc, O2=0.00021, g=9.81)

xlim = 1e14, 1e17
plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)

plt.subplot(gs[0, 0])
plt.plot(int_JO2_PI, PCb_V['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.plot(cum_int_base, PCb.lev, color = 'black', lw =lw, label= 'WACCM6')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')


plt.subplot(gs[0, 1])
plt.plot(int_JO2_Ten, PCb_10pc_V['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.plot(cum_int_ten, PCb.lev, color = 'black', lw =lw, label= 'WACCM6')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)


plt.subplot(gs[1, 0])
plt.plot(int_JO2_One, PCb_1pc_V['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.plot(cum_int_one, PCb.lev, color = 'black', lw =lw, label= 'WACCM6')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)
plt.xlabel('Integrated O2 photolysis', fontsize = 15, weight = 'bold')
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[1, 1])
#plt.plot(int_JO2_Zero1, PCb_01pc_V['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.plot(cum_int_zero1, PCb.lev, color = 'black', lw =lw, label= 'WACCM6')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)
plt.xlabel('Integrated O2 photolysis', fontsize = 15, weight = 'bold')

#%% J rates calculation check with aflux



J_rates = calculate_J_rates(PI_V_60SZA)

plt.figure()
plt.plot(J_rates, PI_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.plot(PI_V_60SZA['variable']['J_sp']['O2',1], PI_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
xlim = (1e-12,1e-3)
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)


#%% UV spectra

solar_dir = '/Users/gregcooke/stellar_files/'
solar_file = 'Proxima_Centauri_at_b.nc'#read in file
ds = xr.open_dataset(solar_dir+solar_file)#attahc file to dataset
ssi = ds['ssi'].isel(time=0) #define dataset from file

PC_b_flux = ssi.values
PC_b_wavelength = ssi.wavelength.values

solar_dir = '/Users/gregcooke/stellar_files/'
solar_file = 'SolarForcingCMIP6piControl_c160921.nc'#read in file
ds = xr.open_dataset(solar_dir+solar_file)#attahc file to dataset

ssi = ds['ssi'].isel(time=0) #define dataset from file

WACCM_Flux = ssi.values
WACCM_wavelength = ssi.wavelength.values


plt.figure(figsize = (9*1.5,5*1.2))
markercolor = 'whitesmoke'
#alpha = 0.75
alpha = 0.9
plt.plot(WACCM_wavelength, WACCM_Flux/1000, lw = 2, color = 'black', label = 'Sun at Earth', zorder = 3, alpha = alpha)
plt.plot(PC_b_wavelength, PC_b_flux/1000, lw = 2, color = '#32d1a7', label = 'Proxima Centauri at b', zorder = 0, alpha = alpha)
plt.ylabel('Top of atmosphere irradiance\nper unit wavelength [W m'+r'$^\mathbf{-2}$'+' nm'+r'$^\mathbf{-1}$'+']',fontsize=15,color='k', weight = 'bold')
plt.xscale('log')
plt.yscale('log')
plt.tick_params(labelsize = 15, length = 4, width = 2, right = True, colors = 'k')#, labelbottom = False) #change the size of the ticks and labels  
plt.tick_params(which = 'minor', labelsize = 15, length = 2, width = 1, right = True, colors = 'k')#, labelbottom = False) #change the size of the ticks and labels  
plt.xlim(100,5e3)#WACCM_wavelength.max())

plt.legend(frameon = False,  labelcolor='linecolor', loc = (4), fontsize = 15)

plt.axvspan(100, 400, alpha = 0.25, facecolor = 'slategrey', edgecolor = 'slategrey', lw = 3)
plt.ylim(1e-6,10)
plt.yticks([1e-6, 1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1])
plt.text(200,0.7e-5,'UV',color='#434343',fontsize=20, horizontalalignment = 'center', weight = 'bold')
plt.text(121,0.2,'L-'+r'$\alpha$',color='#434343',fontsize=15, horizontalalignment = 'center', weight = 'bold')
plt.xticks([100, 200, 400, 600, 1000, 2000, 4000], [100, 200, 400, 600, 1000, 2000, 4000])

plt.xlabel('Wavelength [nm]', fontsize = 15,color='k', weight = 'bold')
thick_axes(top = True)

plt.savefig('/Users/gregcooke/python_output/PC_spectra.png', dpi = dpi, bbox_inches = 'tight')#, transparent = True)



#%% Spectra compared to cross sections
XSfactor = 1e5
plt.figure()
plt.plot(PI_pc_V_482SZA['variable']['bins'], XSfactor * (PI_pc_V_482SZA['variable']['cross_J']['O2',1]+PI_pc_V_482SZA['variable']['cross_J']['O2',2]), color = 'black')
plt.plot(PI_pc_V_482SZA['variable']['bins'], PI_pc_V_482SZA['variable']['cross_J']['O3',1]+PI_pc_V_482SZA['variable']['cross_J']['O3',2], color = 'm')
plt.ylim(1e-21, 1e-16)
#plt.plot((df_flux['Wave_Min']+df_flux['Wave_Max'])[:35]/20, df_cross['O2'], color = 'teal')
plt.xlim(150, 250)
plt.yscale('log')


plt.figure()
plt.plot(PCb_V_60SZA['variable']['y'][:,PCb_spec_V_60SZA.index('O3')] * 1e6, PCb_V_60SZA['atm']['pco']/1e3, color = 'black')
plt.plot(LWAV(n_dens(PCb, "O3")), PCb.lev, color = 'black', ls = '--')
plt.plot(PCb_10pc_V_60SZA['variable']['y'][:,PCb_10pc_spec_V_60SZA.index('O3')] * 1e6, PCb_10pc_V_60SZA['atm']['pco']/1e3, color = Ten_pc_color)
plt.plot(LWAV(n_dens(PCb_10pc, "O3")), PCb.lev, color = Ten_pc_color, ls = '--')
plt.plot(PCb_1pc_V_60SZA['variable']['y'][:,PCb_1pc_spec_V_60SZA.index('O3')] * 1e6, PCb_1pc_V_60SZA['atm']['pco']/1e3, color = One_pc_color)
plt.plot(LWAV(n_dens(PCb_1pc, "O3")), PCb.lev, color = One_pc_color, ls = '--')
plt.plot(PCb_01pc_V_60SZA['variable']['y'][:,PCb_01pc_spec_V_60SZA.index('O3')] * 1e6, PCb_01pc_V_60SZA['atm']['pco']/1e3, color = Zero1_pc_color)
plt.plot(LWAV(n_dens(PCb_01pc, "O3")), PCb.lev, color = Zero1_pc_color, ls = '--')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e3, 1e-5)
plt.xlim(1e12, 1e19)


plt.figure()
plt.plot(PCb_V_60SZA['variable']['y'][:,PCb_spec_V_60SZA.index('O3')] * 1e6, PCb_V_60SZA['atm']['pco']/1e3, color = 'black')
O3 = LWAV_spec(n_dens(PCb, "O3"), time = False, latitude = np.arange(47-5, 47+5))
plt.plot(O3, PCb.lev, color = 'black', ls = '--')
plt.plot(PCb_10pc_V_60SZA['variable']['y'][:,PCb_10pc_spec_V_60SZA.index('O3')] * 1e6, PCb_10pc_V_60SZA['atm']['pco']/1e3, color = Ten_pc_color)
O3 = LWAV_spec(n_dens(PCb_10pc, "O3"), time = False, latitude = np.arange(47-5, 47+5))
plt.plot(O3, PCb.lev, color = Ten_pc_color, ls = '--')
plt.plot(PCb_1pc_V_60SZA['variable']['y'][:,PCb_1pc_spec_V_60SZA.index('O3')] * 1e6, PCb_1pc_V_60SZA['atm']['pco']/1e3, color = One_pc_color)
O3 = LWAV_spec(n_dens(PCb_1pc, "O3"), time = False, latitude = np.arange(47-5, 47+5))
plt.plot(O3, PCb.lev, color = One_pc_color, ls = '--')
plt.plot(PCb_01pc_V_60SZA['variable']['y'][:,PCb_01pc_spec_V_60SZA.index('O3')] * 1e6, PCb_01pc_V_60SZA['atm']['pco']/1e3, color = Zero1_pc_color)
O3 = LWAV_spec(n_dens(PCb_01pc, "O3"), time = False, latitude = np.arange(47-5, 47+5))
plt.plot(O3, PCb.lev, color = Zero1_pc_color, ls = '--')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e3, 1e-5)
plt.xlim(1e12, 1e19)



'''


Is it something to do with diffuse radiation?
Plot the cumulated integrated ozone column compared to the
cumulated integrated O2 column
I think ozone is self shielding in the Herzberg continuum because
of the low flux
This then affects intgrated O2 photolysis unlike the early Earth cases.
Self shielding doesn't take place in 3D mdoels because
ozone accumulates on the nightside
And yet, 3D models end up having overall more ozone? Are you sure about this?
'''

plt.figure()
plt.plot(XSfactor*PCb_V_60SZA['variable']['y'][:,PCb_spec_V_60SZA.index('O3')] * 1e6, PCb_V_60SZA['atm']['pco']/1e3, color = 'black')
plt.plot(PCb_V_60SZA['variable']['y'][:,PCb_spec_V_60SZA.index('O2')] * 1e6, PCb_V_60SZA['atm']['pco']/1e3, color = 'black', ls = '--')
plt.plot(XSfactor*PCb_10pc_V_60SZA['variable']['y'][:,PCb_10pc_spec_V_60SZA.index('O3')] * 1e6, PCb_10pc_V_60SZA['atm']['pco']/1e3, color = Ten_pc_color)
plt.plot(PCb_10pc_V_60SZA['variable']['y'][:,PCb_10pc_spec_V_60SZA.index('O2')] * 1e6, PCb_10pc_V_60SZA['atm']['pco']/1e3, color = Ten_pc_color, ls = '--')
plt.plot(XSfactor*PCb_1pc_V_60SZA['variable']['y'][:,PCb_1pc_spec_V_60SZA.index('O3')] * 1e6, PCb_1pc_V_60SZA['atm']['pco']/1e3, color = One_pc_color)
plt.plot(PCb_1pc_V_60SZA['variable']['y'][:,PCb_1pc_spec_V_60SZA.index('O2')] * 1e6, PCb_1pc_V_60SZA['atm']['pco']/1e3, color = One_pc_color, ls = '--')
plt.plot(XSfactor*PCb_01pc_V_60SZA['variable']['y'][:,PCb_01pc_spec_V_60SZA.index('O3')] * 1e6, PCb_01pc_V_60SZA['atm']['pco']/1e3, color = Zero1_pc_color)
plt.plot(PCb_01pc_V_60SZA['variable']['y'][:,PCb_01pc_spec_V_60SZA.index('O2')] * 1e6, PCb_01pc_V_60SZA['atm']['pco']/1e3, color = Zero1_pc_color, ls = '--')
#plt.plot(XSfactor*One_pc_V_482SZA['variable']['y'][:,One_pc_spec_482SZA.index('O3')] * 1e6, One_pc_V_482SZA['atm']['pco']/1e3, color = One_pc_color)
#plt.plot(One_pc_V_482SZA['variable']['y'][:,One_pc_spec_482SZA.index('O2')] * 1e6, One_pc_V_482SZA['atm']['pco']/1e3, color = One_pc_color, ls = '--')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e3, 1e-5)
plt.xlim(1e19, 1e24)

plt.figure()
plt.plot(XSfactor*PI_pc_V_482SZA['variable']['y'][:,PI_pc_spec_482SZA.index('O3')] * 1e6, PI_pc_V_482SZA['atm']['pco']/1e3, color = 'black')
plt.plot(PI_pc_V_482SZA['variable']['y'][:,PI_pc_spec_482SZA.index('O2')] * 1e6, PI_pc_V_482SZA['atm']['pco']/1e3, color = 'black', ls = '--')
plt.plot(XSfactor*Ten_pc_V_482SZA['variable']['y'][:,Ten_pc_spec_482SZA.index('O3')] * 1e6, Ten_pc_V_482SZA['atm']['pco']/1e3, color = Ten_pc_color)
plt.plot(Ten_pc_V_482SZA['variable']['y'][:,Ten_pc_spec_482SZA.index('O2')] * 1e6, Ten_pc_V_482SZA['atm']['pco']/1e3, color = Ten_pc_color, ls = '--')
plt.plot(XSfactor*One_pc_V_482SZA['variable']['y'][:,One_pc_spec_482SZA.index('O3')] * 1e6, One_pc_V_482SZA['atm']['pco']/1e3, color = One_pc_color)
plt.plot(One_pc_V_482SZA['variable']['y'][:,One_pc_spec_482SZA.index('O2')] * 1e6, One_pc_V_482SZA['atm']['pco']/1e3, color = One_pc_color, ls = '--')
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e3, 1e-5)
plt.xlim(1e15, 1e22)

#%% How does the flux change through the atmosphere?
level = 35
plt.figure(figsize = (10,5))
plt.plot(PCb_V_60SZA['variable']['bins'], PCb_V_60SZA['variable']['sflux'][level], color = 'black')
plt.plot(PCb_10pc_V_60SZA['variable']['bins'], PCb_10pc_V_60SZA['variable']['sflux'][level], color = Ten_pc_color)
plt.plot(PCb_1pc_V_60SZA['variable']['bins'], PCb_1pc_V_60SZA['variable']['sflux'][level], color = One_pc_color)
plt.plot(PCb_01pc_V_60SZA['variable']['bins'], PCb_01pc_V_60SZA['variable']['sflux'][level], color = Zero1_pc_color)
plt.xlim(150, 300)
plt.yscale('log')
plt.ylim(1e-9, 1e3)
plt.title(str(PCb_01pc_V_60SZA['atm']['pco'][level]/1e3) + ' hPa')

plt.figure(figsize = (10,5))
plt.plot(PI_pc_V_482SZA['variable']['bins'], PI_pc_V_482SZA['variable']['sflux'][level], color = 'black')
plt.plot(Ten_pc_V_482SZA['variable']['bins'], Ten_pc_V_482SZA['variable']['sflux'][level], color = Ten_pc_color)
plt.plot(One_pc_V_482SZA['variable']['bins'], One_pc_V_482SZA['variable']['sflux'][level], color = One_pc_color)
plt.plot(Z1_pc_V_482SZA['variable']['bins'], Z1_pc_V_482SZA['variable']['sflux'][level], color = Zero1_pc_color)
plt.xlim(150, 300)
plt.yscale('log')
plt.ylim(1e-9, 1e3)
plt.title(str(PI_pc_V_482SZA['atm']['pco'][level]/1e3) + ' hPa')


#%%



file_path = '/Users/gregcooke/photochem/examples/TOI1468c/atmosphere.txt'
output_path = '/Users/gregcooke/photochem/examples/TOI1468c/TOI1468c_atmosphere.txt'

# 1. Read the file
df = pd.read_csv(file_path, delim_whitespace=True)

# 2. Modify the specific columns
# We use .loc to ensure we are modifying the original dataframe
df['O'] = 0.00001
df['CH4'] = 1e-3
df['NH3'] = 1e-4
df['CO2'] = 1e-5
df['H2'] = 0.9

# 3. Write the file back out
# float_format='%.5E' ensures the scientific notation matches your original file
df.to_csv(output_path, sep=' ', index=False, float_format='%.5E')

print(f"Modified file saved to: {output_path}")


#%%
import numpy as np

# 1. Extract Temperature array
T = PI_V_482SZA['atm']['Tco']
P = PI_V_482SZA['atm']['pco']

# 2. Define rate parameters for the low-pressure limit
A_0 = 2.989029e-28
b_0 = -2.3

# 3. Calculate the low-pressure rate constant (k_0)
# Units: cm^6 / molecule^2 / s
k_0_photochem = A_0 * (T ** b_0)

# (Optional) If you want the exact termolecular falloff without multiplying by [M]:
A_inf = 2.8e-11
k_inf = A_inf * (T ** 0.0)
M_density = (P / (1.380649e-23 * T)) * 1e-6 
k_termolecular = k_0 / (1 + (k_0 * M_density / k_inf))
# Note: For this specific reaction in the atmosphere, k_0 and k_termolecular 
# will be virtually identical.

# 4. Plot the corrected low-pressure rate constant
plt.figure(figsize=(6, 8))
plt.plot(k_0, P / 1e3, label='New Low-P Rate (k_0)', color='blue')
plt.plot(k_termolecular, P / 1e3, label='New Low-P Rate (k_0)', color='blue')

# Overlay the original rate constant for a direct visual check
plt.plot(PI_V_482SZA['variable']['k'][1075], P / 1e3, 
         label='Original 1075', color='red', linestyle='--')

# Match your original plot formatting
plt.ylim(1e3, 1e-5)
plt.yscale('log')
plt.xlabel('Rate Constant (cm$^6$ molecule$^{-2}$ s$^{-1}$)')
plt.ylabel('P / 1e3')
plt.title('O + O$_2$ + M Termolecular Rate Constant')
plt.legend()
plt.show()


T = PI_V_482SZA['atm']['Tco']
P = PI_V_482SZA['atm']['pco']
k_atmos = 5.12e-34 * (298.0 / T)**2.6 * np.exp(40.5 / T)

A_0_tbdy = 5.70e-34
B_0_tbdy = -2.6
k_kasting = A_0_tbdy * (T / 300.0)**B_0_tbdy

#%%

plt.figure()
width = 15
plt.plot(k_0_photochem, P / 1e3, label='New Low-P Rate (k_0)', color='blue')
plt.plot(k_0, P / 1e3, label='New Low-P Rate (k_0)', color='blue')
plt.plot(k_atmos, P / 1e3, label="Atmos", color='darkorange')
plt.plot(k_kasting, P / 1e3, label="Atmos", color='teal')
plt.plot(LWAV_spec(Pre_h0.usr_O_O2, time = False, latitude = np.arange(47-width, 47+width)), Pre_pc_h0.lev, color = 'k');
plt.ylim(1e3,1e-5);  
plt.plot(PI_V_482SZA['variable']['k'][1075], PI_V_482SZA['atm']['pco']/1e3, color = 'm')
plt.ylim(1e3,1e-1)
plt.yscale('log');
plt.xlim(0.5e-33, 1.9e-33)

T = Ten_pc_V_482SZA['atm']['Tco']
P = Ten_pc_V_482SZA['atm']['pco']
k_atmos = 5.12e-34 * (298.0 / T)**2.6 * np.exp(40.5 / T)

A_0_tbdy = 5.70e-34
B_0_tbdy = -2.6
k_kasting = A_0_tbdy * (T / 300.0)**B_0_tbdy

# 2. Define rate parameters for the low-pressure limit
A_0 = 2.989029e-28
b_0 = -2.3

# 3. Calculate the low-pressure rate constant (k_0)
# Units: cm^6 / molecule^2 / s
k_0_photochem = A_0 * (T ** b_0)

plt.figure()
width = 15
plt.plot(k_0_photochem, P / 1e3, label='New Low-P Rate (k_0)', color='blue')
plt.plot(k_atmos, P / 1e3, label="Atmos", color='darkorange')
plt.plot(k_kasting, P / 1e3, label="Atmos", color='teal')
plt.plot(LWAV_spec(Pre_h0.usr_O_O2, time = False, latitude = np.arange(47-width, 47+width)), Pre_pc_h0.lev, color = 'k');
plt.ylim(1e3,1e-5);  
plt.plot(Ten_pc_V_482SZA['variable']['k'][1075], One_pc_V_482SZA['atm']['pco']/1e3, color = 'm')
plt.ylim(1e3,1e-1)
plt.yscale('log');
plt.xlim(1.3e-33, 3e-33)

T = One_pc_V_482SZA['atm']['Tco']
P = One_pc_V_482SZA['atm']['pco']
k_atmos = 5.12e-34 * (298.0 / T)**2.6 * np.exp(40.5 / T)

A_0_tbdy = 5.70e-34
B_0_tbdy = -2.6
k_kasting = A_0_tbdy * (T / 300.0)**B_0_tbdy

# 2. Define rate parameters for the low-pressure limit
A_0 = 2.989029e-28
b_0 = -2.3

# 3. Calculate the low-pressure rate constant (k_0)
# Units: cm^6 / molecule^2 / s
k_0_photochem = A_0 * (T ** b_0)

plt.figure()
width = 15
plt.plot(k_0_photochem, P / 1e3, label='New Low-P Rate (k_0)', color='blue')
plt.plot(k_atmos, P / 1e3, label="Atmos", color='darkorange')
plt.plot(k_kasting, P / 1e3, label="Atmos", color='teal')
plt.plot(LWAV_spec(Pre_h0.usr_O_O2, time = False, latitude = np.arange(47-width, 47+width)), Pre_pc_h0.lev, color = 'k');
plt.ylim(1e3,1e-5);  
plt.plot(One_pc_V_482SZA['variable']['k'][1075], One_pc_V_482SZA['atm']['pco']/1e3, color = 'm')
plt.ylim(1e3,1e-1)
plt.yscale('log');
plt.xlim(1.3e-33, 3e-33)

T = Z1_pc_V_482SZA['atm']['Tco']
P = Z1_pc_V_482SZA['atm']['pco']
k_atmos = 5.12e-34 * (298.0 / T)**2.6 * np.exp(40.5 / T)

A_0_tbdy = 5.70e-34
B_0_tbdy = -2.6
k_kasting = A_0_tbdy * (T / 300.0)**B_0_tbdy

# 2. Define rate parameters for the low-pressure limit
A_0 = 2.989029e-28
b_0 = -2.3

# 3. Calculate the low-pressure rate constant (k_0)
# Units: cm^6 / molecule^2 / s
k_0_photochem = A_0 * (T ** b_0)

plt.figure()
width = 15
plt.plot(k_0_photochem, P / 1e3, label='New Low-P Rate (k_0)', color='blue')
plt.plot(k_atmos, P / 1e3, label="Atmos", color='darkorange')
plt.plot(k_kasting, P / 1e3, label="Atmos", color='teal')
plt.plot(LWAV_spec(Pre_h0.usr_O_O2, time = False, latitude = np.arange(47-width, 47+width)), Pre_pc_h0.lev, color = 'k');
plt.ylim(1e3,1e-5);  
plt.plot(Z1_pc_V_482SZA['variable']['k'][1075], One_pc_V_482SZA['atm']['pco']/1e3, color = 'm')
plt.ylim(1e3,1e-1)
plt.yscale('log');
plt.xlim(1.3e-33, 3e-33)
#%% Ozone loss

xlim = (1e-2, 1e12)
ylim = (1e3, 1e-3)

plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)
gs.update(wspace = 0.2, hspace = 0.2)
plt.subplot(gs[0,0])
reaction = PI_V_482SZA['variable']['k'][707]* PI_V_482SZA['variable']['k'][707]*1e6
O =  PI_V_482SZA['variable']['y'][:,PI_spec_482SZA.index('O')] * 1e6
O3 = PI_V_482SZA['variable']['y'][:,PI_spec_482SZA.index('O3')] * 1e6
plt.plot((reaction * O * O3)*1e3, PI_V_482SZA['atm']['pco']/1e3, color = 'm', lw = 2, label = 'VULCAN')
o3_loss, pressure_hpa = compute_o3loss(
    mech_file=P_path+"100pc/zahnle_earth.yaml",
    settings_file=P_path+"100pc/settings_100pc.yaml",
    flux_file=P_path+"100pc/Sun_0.0Ga.txt",
    pt_file=P_path+"100pc/Earth_100pc_48.2.txt",
    atol=1e-23,
    verbose=0
)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.plot((o3_loss), pressure_hpa, color = 'b', lw = 2, label = 'Photochem')
plt.ylim(ylim); plt.xlim(xlim)
plt.xscale('log'); plt.yscale('log')
thick_axes(top = True)
plt.legend(loc = 0, fontsize = 15, frameon = False)

plt.subplot(gs[0,1])
reaction = Ten_pc_V_482SZA['variable']['k'][707]*Ten_pc_V_482SZA['variable']['k'][707]*1e6
O =  Ten_pc_V_482SZA['variable']['y'][:,Ten_pc_spec_482SZA.index('O')] * 1e6
O3 = Ten_pc_V_482SZA['variable']['y'][:,Ten_pc_spec_482SZA.index('O3')] * 1e6
plt.plot((reaction * O * O3)*1e3, Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = 2)
o3_loss, pressure_hpa = compute_o3loss(
    mech_file=P_path+"10pc/zahnle_earth.yaml",
    settings_file=P_path+"10pc/settings_100pc.yaml",
    flux_file=P_path+"10pc/Sun_0.0Ga.txt",
    pt_file=P_path+"10pc/Earth_10pc_48.2.txt",
    atol=1e-23,
    verbose=0
)
plt.plot((o3_loss), pressure_hpa, color = 'b', lw = 2)
plt.ylim(ylim); plt.xlim(xlim)
plt.xscale('log'); plt.yscale('log')
thick_axes(top = True)

plt.subplot(gs[1,0])
reaction = One_pc_V_482SZA['variable']['k'][707]* One_pc_V_482SZA['variable']['k'][707]*1e6
O =  One_pc_V_482SZA['variable']['y'][:,One_pc_spec_482SZA.index('O')] * 1e6
O3 = One_pc_V_482SZA['variable']['y'][:,One_pc_spec_482SZA.index('O3')] * 1e6
plt.plot((reaction * O * O3)*1e3, PI_V_482SZA['atm']['pco']/1e3, color = 'm', lw = 2)
o3_loss, pressure_hpa = compute_o3loss(
    mech_file=P_path+"1pc/zahnle_earth.yaml",
    settings_file=P_path+"1pc/settings_100pc.yaml",
    flux_file=P_path+"1pc/Sun_0.0Ga.txt",
    pt_file=P_path+"1pc/Earth_1pc_48.2.txt",
    atol=1e-23,
    verbose=0
)
plt.plot((o3_loss), pressure_hpa, color = 'b', lw = 2)
plt.ylim(ylim); plt.xlim(xlim)
plt.xscale('log'); plt.yscale('log')
thick_axes(top = True)
plt.xlabel('O'+sub(3)+' loss rate [molecules m'+sup(-3)+'s'+sup(-1)+']', fontsize = 15, weight = 'bold')
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[1,1])
reaction = Z1_pc_V_482SZA['variable']['k'][707]*Z1_pc_V_482SZA['variable']['k'][707]*1e6
O =  Z1_pc_V_482SZA['variable']['y'][:,Z1_pc_spec_482SZA.index('O')] * 1e6
O3 = Z1_pc_V_482SZA['variable']['y'][:,Z1_pc_spec_482SZA.index('O3')] * 1e6
plt.plot((reaction * O * O3)*1e3, Z1_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = 2)
o3_loss, pressure_hpa = compute_o3loss(
    mech_file=P_path+"0.1pc/zahnle_earth.yaml",
    settings_file=P_path+"0.1pc/settings_100pc.yaml",
    flux_file=P_path+"0.1pc/Sun_0.0Ga.txt",
    pt_file=P_path+"0.1pc/Earth_0.1pc_48.2.txt",
    atol=1e-23,
    verbose=0
)
plt.plot((o3_loss), pressure_hpa, color = 'b', lw = 2)
plt.ylim(ylim); plt.xlim(xlim)
plt.xscale('log'); plt.yscale('log')
thick_axes(top = True)
plt.xlabel('O'+sub(3)+' loss rate [molecules m'+sup(-3)+'s'+sup(-1)+']', fontsize = 15, weight = 'bold')
plt.savefig('/Users/gregcooke/python_output/O3_loss_all_models.png')


#%% VULCAN Proxima Centauri runs


plt.figure()

O3 =  PCb_V['variable']['ymix'][:,PCb_V_spec.index('O3')]
plt.plot(O3, PCb_V['atm']['pco']/1e3, color = 'm', lw = 2)
O3 =  PCb_10pc_V['variable']['ymix'][:,PCb_10pc_V_spec.index('O3')]
plt.plot(O3, PCb_1pc_V['atm']['pco']/1e3, color = 'm', ls = ':', lw = 2)
O3 =  PCb_1pc_V['variable']['ymix'][:,PCb_1pc_V_spec.index('O3')]
plt.plot(O3, PCb_1pc_V['atm']['pco']/1e3, color = 'm', ls = '--', lw = 2)

plt.yscale('log'); plt.xscale('log')
plt.ylim(1e3, 1e-5); plt.xlim(1e-8, 1e-5)




#%% TUV outputs


# ===================================================================
# Read multiple files
# ===================================================================

files = {
    "Z1": "/Users/gregcooke/Z1_pc.txt",
    "Z1V": "/Users/gregcooke/VULCAN.txt",
}

spectral = {}
dose = {}

for name, path in files.items():
    spectral[name], dose[name] = read_usrout(path)

# Optional aliases
df_Z1 = spectral["Z1"]
df_Z1V  = spectral["Z1V"]


# ===================================================================
# Example plots
# ===================================================================

DNA = pd.read_csv('/Users/gregcooke/V5.4/DATAS1/dna.setlow.new', delim_whitespace=True, 
                  names = ['Wav', 'Func'], skiprows = 10)


plt.figure()

for name, df in spectral.items():
    plt.plot(df["wavelength_nm"], df["sza_40"], label=name)
plt.plot(DNA['Wav'], DNA['Func'] )
plt.yscale('log')
plt.yscale('log'); plt.ylim(1e-4, 3)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Spectral irradiance")
plt.legend()

plt.figure()
for name, df in dose.items():
    plt.plot(df["sza_deg"], df["UV_index"], label=name)

plt.xlabel("Solar zenith angle")
plt.ylabel("UV index")
plt.legend()

plt.figure()
for name, df in dose.items():
    plt.plot(df["sza_deg"], df["DNA_damage"], label=name)
plt.xlabel("Solar zenith angle")
plt.ylabel("DNA damage")
plt.yscale('log')
plt.legend()

#%%

#%% O3 production and loss analysis (VULCAN, 48.2 deg SZA)

# var['y'] is number density in molecules cm^-3; var['k'] are rate coefficients.
# Total O3 production / loss sums every reaction (and reverse step) involving O3.

O3_482SZA_cases = [
    ('PI (100% O2)', PI_V_482SZA, PI_spec_482SZA, 'k'),
    ('10% O2', Ten_pc_V_482SZA, Ten_pc_spec_482SZA, Ten_pc_color),
    ('1% O2', One_pc_V_482SZA, One_pc_spec_482SZA, One_pc_color),
    ('0.1% O2', Z1_pc_V_482SZA, Z1_pc_spec_482SZA, Zero1_pc_color),
]

O3_482_budgets = {}
for label, vulcan_ds, spec_list, color in O3_482SZA_cases:
    O3_482_budgets[label] = vulcan_o3_production_loss(vulcan_ds, spec_list, to_m3=True)

# --- Profile plot: production, loss, net ---
xlim_o3 = (1e8, 1e15)
ylim_o3 = (1e3, 1e-1)

fig = plt.figure(figsize=(14, 10))
gs_o3 = gridspec.GridSpec(2, 2, wspace=0.2, hspace=0.2)

for ax_idx, (label, vulcan_ds, spec_list, color) in enumerate(O3_482SZA_cases):
    b = O3_482_budgets[label]
    ax = fig.add_subplot(gs_o3[ax_idx // 2, ax_idx % 2])
    ax.plot(b['loss'], b['pco'], color='m', lw=2, label='Loss')
    ax.plot(b['prod'], b['pco'], color='b', lw=2, label='Production')
    ax.plot(np.abs(b['net']), b['pco'], color=color, lw=2, ls='--', label='|Net|')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xlim_o3)
    ax.set_ylim(ylim_o3)
    ax.set_title(label, fontsize=15, weight='bold')
    thick_axes(top=True)
    if ax_idx in (0, 2):
        ax.set_ylabel('Pressure [hPa]', fontsize=15, weight='bold')
    if ax_idx >= 2:
        ax.set_xlabel(
            'O' + sub(3) + ' rate [molecules m' + sup(-3) + ' s' + sup(-1) + ']',
            fontsize=15,
            weight='bold',
        )
    if ax_idx == 0:
        ax.legend(loc=0, fontsize=12, frameon=False)

fig.suptitle('O' + sub(3) + ' production and loss (VULCAN, SZA = 48.2°)', fontsize=16, weight='bold')
plt.savefig('/Users/gregcooke/python_output/O3_prod_loss_482SZA.png', dpi=dpi, bbox_inches='tight')

# --- Chapman-cycle check: O + O2 + M -> O3 + M  vs  O + O3 -> O2 + O2 ---
fig2, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
chapman_prod = 'O + O2 + M -> O3 + M'
chapman_loss = 'O + O3 -> O2 + O2'

for ax, (label, _, _, color) in zip(axes, O3_482SZA_cases[:2]):
    b = O3_482_budgets[label]
    if chapman_prod in b['prod_by_rxn']:
        ax.plot(b['prod_by_rxn'][chapman_prod], b['pco'], color='b', lw=2, label='O + O2 + M')
    if chapman_loss in b['loss_by_rxn']:
        ax.plot(b['loss_by_rxn'][chapman_loss], b['pco'], color='m', lw=2, label='O + O3')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(xlim_o3)
    ax.set_ylim(ylim_o3)
    ax.set_title(label, fontsize=15, weight='bold')
    thick_axes(top=True)
    ax.legend(loc=0, fontsize=12, frameon=False)

axes[0].set_ylabel('Pressure [hPa]', fontsize=15, weight='bold')
for ax in axes:
    ax.set_xlabel(
        'O' + sub(3) + ' rate [molecules m' + sup(-3) + ' s' + sup(-1) + ']',
        fontsize=15,
        weight='bold',
    )
fig2.suptitle('Chapman reactions', fontsize=16, weight='bold')
plt.savefig('/Users/gregcooke/python_output/O3_chapman_482SZA.png', dpi=dpi, bbox_inches='tight')

# --- Print dominant channels near the O3 peak (1% O2 case) ---
ref_label = '1% O2'
ref_budget = O3_482_budgets[ref_label]
o3_prof = One_pc_V_482SZA['variable']['y'][:, One_pc_spec_482SZA.index('O3')]
iz_peak = int(np.argmax(o3_prof))
prod_top, loss_top = top_o3_channels(ref_budget, iz_peak, n=10)

print('\n=== O3 budget at O3 peak (' + ref_label + ', p = {:.3f} hPa) ==='.format(ref_budget['pco'][iz_peak]))
print('Production {:.3e}  Loss {:.3e}  Net {:.3e} molecules m^-3 s^-1'.format(
    ref_budget['prod'][iz_peak],
    ref_budget['loss'][iz_peak],
    ref_budget['net'][iz_peak],
))
print('\nTop production channels:')
for rxn, rate in prod_top:
    print('  {:.3e}  {}'.format(rate, rxn))
print('\nTop loss channels:')
for rxn, rate in loss_top:
    print('  {:.3e}  {}'.format(rate, rxn))

# --- NOx and HOx: individual reaction profiles (1% O2) ---
xlim_family = (1e5, 1e14)
ylim_family = (1000.0, 1e-2)
# All rates are instantaneous O3 loss/production: molecules m^-3 s^-1
o3_rate_label = 'O' + sub(3) + ' rate [molecules m' + sup(-3) + ' s' + sup(-1) + ']'
o3_loss_label = 'O' + sub(3) + ' loss rate [molecules m' + sup(-3) + ' s' + sup(-1) + ']'

for family, fname in (('nox', 'NO' + sub('x')), ('hox', 'HO' + sub('x'))):
    prod_fam, loss_fam = o3_channels_in_family(ref_budget, family=family)
    fig_fam, axes_fam = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for ax, channels, panel in zip(axes_fam, (loss_fam, prod_fam), ('Loss', 'Production')):
        for rxn, rate in sorted(channels.items(), key=lambda item: np.max(item[1]), reverse=True):
            ax.plot(rate, ref_budget['pco'], lw=1.5, label=rxn.replace('REV: ', '(rev) '))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(xlim_family)
        ax.set_ylim(ylim_family)
        ax.set_title(panel, fontsize=14, weight='bold')
        thick_axes(top=True)
        ax.legend(loc='best', fontsize=7, frameon=False)
    axes_fam[0].set_ylabel('Pressure [hPa]', fontsize=15, weight='bold')
    for ax in axes_fam:
        ax.set_xlabel(o3_rate_label, fontsize=14, weight='bold')
    fig_fam.suptitle(
        ref_label + ': O' + sub(3) + ' ' + fname + ' channels (SZA = 48.2°)',
        fontsize=16,
        weight='bold',
    )
    plt.tight_layout()
    plt.savefig(
        '/Users/gregcooke/python_output/O3_{}_channels_1pc_482SZA.png'.format(family),
        dpi=dpi,
        bbox_inches='tight',
    )

# --- Unit check (1% O2): channel sum = total; categories partition loss ---
_cat = o3_loss_by_category(ref_budget)
_ch_sum = sum(v[iz_peak] for v in ref_budget['loss_by_rxn'].values())
print('\n=== O3 rate units: {} ==='.format(ref_budget['rate_unit']))
print('Channel sum vs total loss at O3 peak: {:.6e} vs {:.6e}'.format(_ch_sum, ref_budget['loss'][iz_peak]))
print('Loss partition at O3 peak (molecules m^-3 s^-1):')
for _k in ('nox', 'hox_catalytic', 'photolysis', 'chapman', 'other', 'total'):
    print('  {:16s} {:.6e}'.format(_k, _cat[_k][iz_peak]))
print('  {:16s} {:.6e}  (NOx + HOx cat. + phot. + Chap.)'.format(
    'sum categories',
    _cat['nox'][iz_peak] + _cat['hox_catalytic'][iz_peak] + _cat['photolysis'][iz_peak]
    + _cat['chapman'][iz_peak] + _cat['other'][iz_peak],
))

# --- NOx vs HOx catalytic O3 loss (directly comparable; same units, no photolysis in HOx) ---
fig_tot, ax_tot = plt.subplots(figsize=(8, 6))

for label, _, _, color in O3_482SZA_cases:
    b = O3_482_budgets[label]
    cat = o3_loss_by_category(b)
    ax_tot.plot(cat['nox'], b['pco'], color=color, lw=2, ls='-', label='NO' + sub('x') + ' cat. ' + label)
    ax_tot.plot(
        cat['hox_catalytic'], b['pco'], color=color, lw=2, ls='--',
        label='HO' + sub('x') + ' cat. ' + label,
    )

ax_tot.set_xscale('log')
ax_tot.set_yscale('log')
ax_tot.set_xlim(xlim_family)
ax_tot.set_ylim(ylim_family)
ax_tot.set_xlabel(o3_loss_label, fontsize=15, weight='bold')
ax_tot.set_ylabel('Pressure [hPa]', fontsize=15, weight='bold')
ax_tot.set_title(
    'Catalytic O' + sub(3) + ' loss: NO' + sub('x') + ' vs HO' + sub('x') + ' (SZA = 48.2°)',
    fontsize=15, weight='bold',
)
thick_axes(top=True)
ax_tot.legend(loc='best', fontsize=8, frameon=False, ncol=2)
plt.tight_layout()
plt.savefig(
    '/Users/gregcooke/python_output/O3_NOx_vs_HOx_catalytic_loss_482SZA.png',
    dpi=dpi,
    bbox_inches='tight',
)

# --- Full O3 loss budget: catalytic NOx/HOx + photolysis + Chapman + other ---
fig_part, ax_part = plt.subplots(figsize=(9, 6))

for label, _, _, color in O3_482SZA_cases:
    b = O3_482_budgets[label]
    cat = o3_loss_by_category(b)
    if label == O3_482SZA_cases[0][0]:
        ax_part.plot(cat['nox'], b['pco'], color=color, lw=2, ls='-', label='NO' + sub('x') + ' cat.')
        ax_part.plot(cat['hox_catalytic'], b['pco'], color=color, lw=2, ls='--', label='HO' + sub('x') + ' cat.')
        ax_part.plot(cat['photolysis'], b['pco'], color='gray', lw=1.5, ls='-', label='Photolysis')
        ax_part.plot(cat['chapman'], b['pco'], color='gray', lw=1.5, ls='--', label='Chapman O+O' + sub(3))
        ax_part.plot(cat['other'], b['pco'], color='gray', lw=1.5, ls=':', label='Other')
    else:
        ax_part.plot(cat['nox'], b['pco'], color=color, lw=2, ls='-')
        ax_part.plot(cat['hox_catalytic'], b['pco'], color=color, lw=2, ls='--')

ax_part.set_xscale('log')
ax_part.set_yscale('log')
ax_part.set_xlim(xlim_family)
ax_part.set_ylim(ylim_family)
ax_part.set_xlabel(o3_loss_label, fontsize=15, weight='bold')
ax_part.set_ylabel('Pressure [hPa]', fontsize=15, weight='bold')
ax_part.set_title(
    'O' + sub(3) + ' loss budget by category (SZA = 48.2°)',
    fontsize=15, weight='bold',
)
thick_axes(top=True)
ax_part.legend(loc='best', fontsize=9, frameon=False)
plt.tight_layout()
plt.savefig(
    '/Users/gregcooke/python_output/O3_loss_partition_482SZA.png',
    dpi=dpi,
    bbox_inches='tight',
)

# --- Key O3 + NOx / HOx loss channels vs O2 (one panel per reaction; all include O3) ---
key_nox_loss = [
    'NO + O3 -> NO2 + O2',
    'NO2 + O3 -> NO3 + O2',
    'N + O3 -> O2 + NO',
]
key_hox_loss = [
    'O3 + H -> OH + O2',
    'OH + O3 -> HO2 + O2',
    'HO2 + O3 -> OH + O2 + O2',
]

for key_list, family_title, out_name in (
    (key_nox_loss, 'NO' + sub('x') + ' O' + sub(3) + ' loss vs O' + sub(2), 'NOx_loss_compare'),
    (key_hox_loss, 'HO' + sub('x') + ' O' + sub(3) + ' loss vs O' + sub(2), 'HOx_loss_compare'),
):
    fig_cmp, axes_cmp = plt.subplots(1, len(key_list), figsize=(4.5 * len(key_list), 5), sharey=True)
    if len(key_list) == 1:
        axes_cmp = [axes_cmp]
    for ax, rxn in zip(axes_cmp, key_list):
        for label, _, _, color in O3_482SZA_cases:
            b = O3_482_budgets[label]
            if rxn in b['loss_by_rxn']:
                ax.plot(b['loss_by_rxn'][rxn], b['pco'], color=color, lw=2, label=label)
            if label in O3_waccm_482_budgets:
                bw = O3_waccm_482_budgets[label]
                if rxn in bw['loss_by_rxn']:
                    ax.plot(
                        bw['loss_by_rxn'][rxn], bw['pco'],
                        color=color_WACCM, lw=2, ls=':',
                        label='WACCM6 ' + label,
                    )
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(xlim_family)
        ax.set_ylim(ylim_family)
        ax.set_title(rxn, fontsize=10, weight='bold')
        thick_axes(top=True)
        ax.legend(loc='best', fontsize=8, frameon=False)
    axes_cmp[0].set_ylabel('Pressure [hPa]', fontsize=14, weight='bold')
    for ax in axes_cmp:
        ax.set_xlabel(
            'O' + sub(3) + ' loss [molecules m' + sup(-3) + ' s' + sup(-1) + ']',
            fontsize=12,
            weight='bold',
        )
    fig_cmp.suptitle(family_title + ' (VULCAN SZA = 48.2°; WACCM6 zonal mean)', fontsize=15, weight='bold')
    plt.tight_layout()
    plt.savefig(
        '/Users/gregcooke/python_output/O3_{}_482SZA.png'.format(out_name),
        dpi=dpi,
        bbox_inches='tight',
    )

#%% O3 production and loss analysis (WACCM6 CAM h0, latitude-weighted mean)

# Requires h0 diagnostics <reactant>_O3 (e.g. NO_O3, OH_O3). Currently on Pre_h0;
# add matching cam.h0 files for other O2 levels when available.
O3_WACCM_482_DATASETS = [
    ('PI (100% O2)', Pre_h0),
    ('10% O2', Ten_pc_h0),
    ('1% O2', One_pc_h0),
    ('0.1% O2', Zero1_pc_h0),
]

print('\n=== Building WACCM6 O3 catalytic budgets (cam.h0 *_O3 diagnostics) ===')
O3_waccm_482_budgets = {}
for _wlabel, _wds in O3_WACCM_482_DATASETS:
    if waccm_has_o3_loss_diagnostics(_wds):
        O3_waccm_482_budgets[_wlabel] = waccm_o3_catalytic_budget(_wds)
        print('  {}: OK'.format(_wlabel))
    else:
        print('  {}: skipped (no *_O3 variables in this cam.h0 file)'.format(_wlabel))

#%% O3 production and loss analysis (Photochem, 48.2 deg SZA)

O3_PHOTO_482_CASES = [
    ('PI (100% O2)', '100pc', 'k'),
    ('10% O2', '10pc', Ten_pc_color),
    ('1% O2', '1pc', One_pc_color),
    ('0.1% O2', '0.1pc', Zero1_pc_color),
]

print('\n=== Building Photochem O3 budgets (48.2 deg SZA) ===')
O3_photo_482_budgets = {}
for label, folder, _color in O3_PHOTO_482_CASES:
    print('  {}'.format(label))
    O3_photo_482_budgets[label] = photochem_o3_budget_for_folder(folder, sza_tag='48.2')

# --- Photochem: catalytic NOx vs HOx O3 loss (all O2 levels) ---
fig_photo_nh, ax_photo_nh = plt.subplots(figsize=(8, 6))
for label, _folder, color in O3_PHOTO_482_CASES:
    b = O3_photo_482_budgets[label]
    cat = o3_loss_by_category(b)
    ax_photo_nh.plot(
        cat['nox'], b['pco'], color=color, lw=2, ls='-',
        label='NO' + sub('x') + ' cat. ' + label,
    )
    ax_photo_nh.plot(
        cat['hox_catalytic'], b['pco'], color=color, lw=2, ls='--',
        label='HO' + sub('x') + ' cat. ' + label,
    )

ax_photo_nh.set_xscale('log')
ax_photo_nh.set_yscale('log')
ax_photo_nh.set_xlim(xlim_family)
ax_photo_nh.set_ylim(ylim_family)
ax_photo_nh.set_xlabel(o3_loss_label, fontsize=15, weight='bold')
ax_photo_nh.set_ylabel('Pressure [hPa]', fontsize=15, weight='bold')
ax_photo_nh.set_title(
    'Photochem: catalytic O' + sub(3) + ' loss, NO' + sub('x') + ' vs HO' + sub('x') + ' (SZA = 48.2°)',
    fontsize=15, weight='bold',
)
thick_axes(top=True)
ax_photo_nh.legend(loc='best', fontsize=8, frameon=False, ncol=2)
plt.tight_layout()
plt.savefig(
    '/Users/gregcooke/python_output/O3_Photo_NOx_vs_HOx_catalytic_loss_482SZA.png',
    dpi=dpi,
    bbox_inches='tight',
)

# --- HOx catalytic O3 loss: VULCAN vs Photochem vs WACCM6 ---
fig_hox_cmp, ax_hox_cmp = plt.subplots(figsize=(9, 6))
for label, _folder, color in O3_PHOTO_482_CASES:
    bv = O3_482_budgets[label]
    bp = O3_photo_482_budgets[label]
    cat_v = o3_loss_by_category(bv)
    cat_p = o3_loss_by_category(bp)
    ax_hox_cmp.plot(
        cat_v['hox_catalytic'], bv['pco'], color=color, lw=2, ls='-',
        label='VULCAN ' + label,
    )
    ax_hox_cmp.plot(
        cat_p['hox_catalytic'], bp['pco'], color=color, lw=2, ls='--',
        label='Photochem ' + label,
    )
    if label in O3_waccm_482_budgets:
        bw = O3_waccm_482_budgets[label]
        cat_w = o3_loss_by_category(bw)
        ax_hox_cmp.plot(
            cat_w['hox_catalytic'], bw['pco'], color=color_WACCM, lw=2, ls=':',
            label='WACCM6 ' + label,
        )

ax_hox_cmp.set_xscale('log')
ax_hox_cmp.set_yscale('log')
ax_hox_cmp.set_xlim(xlim_family)
ax_hox_cmp.set_ylim(ylim_family)
ax_hox_cmp.set_xlabel(o3_loss_label, fontsize=15, weight='bold')
ax_hox_cmp.set_ylabel('Pressure [hPa]', fontsize=15, weight='bold')
ax_hox_cmp.set_title(
    'HO' + sub('x') + ' catalytic O' + sub(3) + ' loss: VULCAN vs Photochem vs WACCM6',
    fontsize=15, weight='bold',
)
thick_axes(top=True)
ax_hox_cmp.legend(loc='best', fontsize=7, frameon=False, ncol=2)
plt.tight_layout()
plt.savefig(
    '/Users/gregcooke/python_output/O3_HOx_VULCAN_Photo_WACCM_482SZA.png',
    dpi=dpi,
    bbox_inches='tight',
)

# --- NOx catalytic O3 loss: VULCAN vs Photochem vs WACCM6 ---
fig_nox_cmp, ax_nox_cmp = plt.subplots(figsize=(9, 6))
for label, _folder, color in O3_PHOTO_482_CASES:
    bv = O3_482_budgets[label]
    bp = O3_photo_482_budgets[label]
    cat_v = o3_loss_by_category(bv)
    cat_p = o3_loss_by_category(bp)
    ax_nox_cmp.plot(
        cat_v['nox'], bv['pco'], color=color, lw=2, ls='-',
        label='VULCAN ' + label,
    )
    ax_nox_cmp.plot(
        cat_p['nox'], bp['pco'], color=color, lw=2, ls='--',
        label='Photochem ' + label,
    )
    if label in O3_waccm_482_budgets:
        bw = O3_waccm_482_budgets[label]
        cat_w = o3_loss_by_category(bw)
        ax_nox_cmp.plot(
            cat_w['nox'], bw['pco'], color=color_WACCM, lw=2, ls=':',
            label='WACCM6 ' + label,
        )

ax_nox_cmp.set_xscale('log')
ax_nox_cmp.set_yscale('log')
ax_nox_cmp.set_xlim(xlim_family)
ax_nox_cmp.set_ylim(ylim_family)
ax_nox_cmp.set_xlabel(o3_loss_label, fontsize=15, weight='bold')
ax_nox_cmp.set_ylabel('Pressure [hPa]', fontsize=15, weight='bold')
ax_nox_cmp.set_title(
    'NO' + sub('x') + ' catalytic O' + sub(3) + ' loss: VULCAN vs Photochem vs WACCM6',
    fontsize=15, weight='bold',
)
thick_axes(top=True)
ax_nox_cmp.legend(loc='best', fontsize=7, frameon=False, ncol=2)
plt.tight_layout()
plt.savefig(
    '/Users/gregcooke/python_output/O3_NOx_VULCAN_Photo_WACCM_482SZA.png',
    dpi=dpi,
    bbox_inches='tight',
)

# --- Print Photochem vs VULCAN HOx/NOx at O3 peak (1% O2) ---
_ref = '1% O2'
_bv = O3_482_budgets[_ref]
_bp = O3_photo_482_budgets[_ref]
_o3_photo = Photo_1pc['O3'].values
_iz = int(np.argmax(_o3_photo))
_cv = o3_loss_by_category(_bv)
_cp = o3_loss_by_category(_bp)
print('\n=== Catalytic O3 loss at O3 peak (1% O2, Photochem O3 profile) ===')
print('  p = {:.3f} hPa (Photochem), {:.3f} hPa (VULCAN)'.format(
    _bp['pco'][_iz], _bv['pco'][_iz]))
print('  NOx loss   VULCAN {:.3e}   Photochem {:.3e} molecules m^-3 s^-1'.format(
    _cv['nox'][_iz], _cp['nox'][_iz]))
print('  HOx loss   VULCAN {:.3e}   Photochem {:.3e} molecules m^-3 s^-1'.format(
    _cv['hox_catalytic'][_iz], _cp['hox_catalytic'][_iz]))
if _ref in O3_waccm_482_budgets:
    _bw = O3_waccm_482_budgets[_ref]
    _cw = o3_loss_by_category(_bw)
    _iz_w = int(np.argmin(np.abs(_bw['pco'] - _bp['pco'][_iz])))
    print('  WACCM6 at {:.3f} hPa: NOx {:.3e}   HOx {:.3e} molecules m^-3 s^-1'.format(
        _bw['pco'][_iz_w], _cw['nox'][_iz_w], _cw['hox_catalytic'][_iz_w]))
elif 'PI (100% O2)' in O3_waccm_482_budgets:
    _bw = O3_waccm_482_budgets['PI (100% O2)']
    _cw = o3_loss_by_category(_bw)
    _iz_w = int(np.argmin(np.abs(_bw['pco'] - _bp['pco'][_iz])))
    print('  WACCM6 PI at {:.3f} hPa (nearest to 1% O2 peak): NOx {:.3e}   HOx {:.3e}'.format(
        _bw['pco'][_iz_w], _cw['nox'][_iz_w], _cw['hox_catalytic'][_iz_w]))

# --- Fractional O3 loss (NOx + HOx + Chapman); log fraction axis 1e-3 to 1 ---
_ylim_frac = ylim_family
_xlim_frac = (1e-3, 1.0)




_o3_loss_fraction_figure(
    O3_482SZA_cases, O3_482_budgets, 'VULCAN', 'O3_loss_fraction_VULCAN_482SZA',
    ylim_press=_ylim_frac, xlim_frac=_xlim_frac,
)
_o3_loss_fraction_figure(
    O3_PHOTO_482_CASES, O3_photo_482_budgets, 'Photochem',
    'O3_loss_fraction_Photochem_482SZA',
    ylim_press=_ylim_frac, xlim_frac=_xlim_frac,
)

# 1% O2: VULCAN vs Photochem vs WACCM6 (WACCM panel only if diagnostics exist for that O2 level)
_frac_rows = [
    ('VULCAN', O3_482_budgets),
    ('Photochem', O3_photo_482_budgets),
]
if _ref in O3_waccm_482_budgets:
    _frac_rows.append(('WACCM6', O3_waccm_482_budgets))
fig_1pc, axes_1pc = plt.subplots(len(_frac_rows), 1, figsize=(7, 3.5 * len(_frac_rows)), sharey=True, sharex=True)
if len(_frac_rows) == 1:
    axes_1pc = [axes_1pc]
for ax, (model, bdict) in zip(axes_1pc, _frac_rows):
    plot_o3_loss_fraction(
        bdict[_ref],
        model + ' (' + _ref + ')',
        ax,
        ylim_press=_ylim_frac,
        xlim_frac=_xlim_frac,
    )
    thick_axes(top=True)
axes_1pc[0].set_ylabel('Pressure [hPa]', fontsize=14, weight='bold')
handles, labels = axes_1pc[0].get_legend_handles_labels()
fig_1pc.legend(handles, labels, loc='lower center', ncol=3, fontsize=10, frameon=False)
fig_1pc.suptitle(
    'O' + sub(3) + ' loss (NO' + sub('x') + '+HO' + sub('x') + '+Chapman): VULCAN vs Photochem vs WACCM6 (1% O'
    + sub(2) + ')',
    fontsize=15, weight='bold',
)
plt.tight_layout()
plt.savefig(
    '/Users/gregcooke/python_output/O3_loss_fraction_VULCAN_Photo_WACCM_1pc_482SZA.png',
    dpi=dpi,
    bbox_inches='tight',
)

# --- Four-panel: VULCAN vs Photochem vs WACCM6 NOx and HOx fractions at each O2 level ---
fig_vxp, axes_vxp = plt.subplots(2, 2, figsize=(12, 10), sharex=True, sharey=True)
axes_vxp = axes_vxp.flatten()
for ax, case in zip(axes_vxp, O3_482SZA_cases):
    label = case[0]
    plot_o3_nox_hox_fraction_vulcan_vs_photo(
        O3_482_budgets[label],
        O3_photo_482_budgets[label],
        ax,
        label,
        xlim_frac=_xlim_frac,
        budget_waccm=O3_waccm_482_budgets.get(label),
    )
    ax.set_ylim(_ylim_frac)
    thick_axes(top=True)
for idx, ax in enumerate(axes_vxp):
    if idx in (0, 2):
        ax.set_ylabel('Pressure [hPa]', fontsize=13, weight='bold')
    else:
        ax.tick_params(labelleft=False)
handles, labels = axes_vxp[0].get_legend_handles_labels()
fig_vxp.legend(
    handles, labels, loc='lower center', ncol=3, fontsize=8,
    frameon=False, bbox_to_anchor=(0.5, -0.04),
)
fig_vxp.suptitle(
    'O' + sub(3) + ' loss fractions (NO' + sub('x') + ', HO' + sub('x') + '): VULCAN vs Photochem vs WACCM6',
    fontsize=16, weight='bold', y=1.02,
)
plt.tight_layout()
plt.savefig(
    '/Users/gregcooke/python_output/O3_NOx_HOx_fraction_VULCAN_Photo_WACCM_482SZA.png',
    dpi=dpi,
    bbox_inches='tight',
)

#%% Stellar spectra used by each photochemical framework
# UV only (100–400 nm), common unit: mW m^-2 nm^-1 at 1 AU.

HC_KASTING = 6.625e-27 * 3.0e10  # erg cm s^-1 (Planck constant * c), as in Kasting/Atmos photo.f
UV_WL_MIN, UV_WL_MAX = 100.0, 400.0


def kasting_phot_bin_to_mwm2_nm(wmin_a, wmax_a, flux_phot_cm2_s):
    """Kasting photos.pdat: photons cm^-2 s^-1 per bin -> mW m^-2 nm^-1 (CorrK4 CHEM/photo.f)."""
    wav_a = 0.5 * (wmin_a + wmax_a)
    delwav_a = wmax_a - wmin_a
    eflux_erg = (HC_KASTING * 1.0e8 / wav_a) * (1.0e4 / 1.0e7) * flux_phot_cm2_s / (0.1 * delwav_a)
    return eflux_erg * 1.0e3  # erg cm^-2 s^-1 nm^-1 -> mW m^-2 nm^-1


def kasting_faruv_bin_to_mwm2_nm(wl_a, sfx, lyman_alpha=False):
    """Kasting faruv_*.pdat: photons cm^-2 s^-1 per 50 A bin (Ly-a special case in photo.f)."""
    if lyman_alpha:
        eflux_erg = (HC_KASTING * 1.0e8 / wl_a) * (1.0e4 / 1.0e7) * sfx * 50.0
    else:
        eflux_erg = (HC_KASTING * 1.0e8 / wl_a) * (1.0e4 / 1.0e7) * sfx / 0.1
    return eflux_erg * 1.0e3


def vulcan_gueymard_to_mwm2_nm(wl_nm, flux_vulcan, wl_ref_mw, flux_ref_mw):
    """Map VULCAN Gueymard/VPL 'erg' files to mW m^-2 nm^-1 using UV match to Thuillier (Sun_now).

    VULCAN RT uses these values as erg cm^-2 s^-1 nm^-1 internally, but the archived
    Gueymard/VPL numbers are offset from Thuillier/WACCM mW m^-2 nm^-1 by ~4e4 in the UV.
  """
    wl_nm = np.asarray(wl_nm, dtype=float)
    flux_vulcan = np.asarray(flux_vulcan, dtype=float)
    ref_wl = np.linspace(UV_WL_MIN, UV_WL_MAX, 80)
    photo_uv = np.interp(ref_wl, wl_ref_mw, flux_ref_mw, left=np.nan, right=np.nan)
    vulcan_uv = np.interp(ref_wl, wl_nm, flux_vulcan, left=np.nan, right=np.nan)
    ok = np.isfinite(photo_uv) & np.isfinite(vulcan_uv) & (vulcan_uv > 0) & (photo_uv > 0)
    scale = np.median(photo_uv[ok] / vulcan_uv[ok])
    return flux_vulcan * scale


def read_two_column_spectrum(path, skiprows=0):
    """Read nm + flux text files (optional # comment or text header line)."""
    try:
        data = np.loadtxt(path, comments='#', skiprows=skiprows)
    except ValueError:
        data = np.loadtxt(path, skiprows=max(skiprows, 1))
    return data[:, 0], data[:, 1]


def read_kasting_faruv(path):
    """Parse faruv_sun.pdat (wavelength in Angstrom, photons cm^-2 s^-1 per bin)."""
    wl_a, sfx = [], []
    with open(path, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                wl_a.append(float(parts[0]))
                sfx.append(float(parts[1]))
            except ValueError:
                continue
    return np.array(wl_a), np.array(sfx)


# --- File paths ---
VULCAN_FLUX_PATH = '/Users/gregcooke/VIH_cases/atm/stellar_flux/Gueymard_solar.txt'
WACCM_FLUX_PATH = '/Users/gregcooke/stellar_files/SolarForcingCMIP6piControl_c160921.nc'
PHOTOCHEM_FLUX_PATH = '/Users/gregcooke/photochem/examples/ModernEarth/Sun_now.txt'
KASTING_PHOTOS_PATH = '/Users/gregcooke/Aoshuang_2_CorrK4_clean/DATA/photos.pdat'
KASTING_FARUV_PATH = '/Users/gregcooke/Aoshuang_2_CorrK4_clean/DATA/faruv_sun.pdat'
ATMOS_ATLAS_PATH = '/Users/gregcooke/atmos/PHOTOCHEM/DATA/FLUX/composite.atl1_1'



# Photochem Modern Earth: mW m^-2 nm^-1 (Thuillier ATLAS composite)
wl_photo, flux_photo = read_two_column_spectrum(PHOTOCHEM_FLUX_PATH)

# WACCM6 CMIP6 solar forcing: mW m^-2 nm^-1
with xr.open_dataset(WACCM_FLUX_PATH) as ds_waccm:
    wl_waccm = ds_waccm['wlen'].values
    ssi = ds_waccm['ssi'].values
    if ssi.ndim == 2:
        ssi = np.nanmean(ssi, axis=0)
    flux_waccm = np.asarray(ssi, dtype=float)

# VULCAN (Gueymard): file header says erg cm^-2 s^-1 nm^-1; scale to mW via UV Thuillier
wl_vulcan, flux_vulcan_raw = read_two_column_spectrum(VULCAN_FLUX_PATH, skiprows=1)
flux_vulcan = vulcan_gueymard_to_mwm2_nm(
    wl_vulcan, flux_vulcan_raw, wl_photo, flux_photo
)

# Kasting CorrK4: photos.pdat (photons cm^-2 s^-1 per bin) + faruv_sun.pdat (UV)
kasting_blocks = read_kasting_file(KASTING_PHOTOS_PATH)
df_kast_flux = kasting_blocks['flux_block']
df_kast_flux[['Wave_Min', 'Wave_Max']] = (
    df_kast_flux['Range_A'].str.split('-', expand=True).astype(float)
)
wl_kast_nm = (df_kast_flux['Wave_Min'].values + df_kast_flux['Wave_Max'].values) / 20.0
flux_kast = np.array([
    kasting_phot_bin_to_mwm2_nm(w0, w1, f)
    for w0, w1, f in zip(
        df_kast_flux['Wave_Min'].values,
        df_kast_flux['Wave_Max'].values,
        df_kast_flux['Flux'].values,
    )
])

wl_faruv_a, sfx_faruv = read_kasting_faruv(KASTING_FARUV_PATH)
flux_faruv = np.array([
    kasting_faruv_bin_to_mwm2_nm(wl, fx, lyman_alpha=(i == len(wl_faruv_a) - 1))
    for i, (wl, fx) in enumerate(zip(wl_faruv_a, sfx_faruv))
])
wl_faruv_nm = wl_faruv_a / 10.0

# Atmos photochem: same Thuillier ATLAS composite as readflux (pstar=sun)
wl_atmos, flux_atmos = read_two_column_spectrum(ATMOS_ATLAS_PATH, skiprows=0)

# --- Figure: UV 100–400 nm ---
_stellar_specs = [
    (wl_vulcan, flux_vulcan, 'VULCAN (Gueymard)', '#1f77b4', '-'),
    (wl_waccm, flux_waccm, 'WACCM6 (CMIP6 SSI)', '#ff7f0e', '-'),
    (wl_photo, flux_photo, 'Photochem (Sun_now)', '#2ca02c', '-'),
    (wl_kast_nm, flux_kast, 'Kasting (photos.pdat)', '#d62728', '-'),
    (wl_faruv_nm, flux_faruv, 'Kasting (faruv_sun.pdat)', '#d62728', '--'),
    (wl_atmos, flux_atmos, 'Atmos (ATLAS1 composite)', '#9467bd', '-'),
]

fig_stellar, ax_uv = plt.subplots(figsize=(9, 7))

for wl, flux, label, color, ls in _stellar_specs:
    m = (
        np.isfinite(wl) & np.isfinite(flux) & (flux > 0)
        & (wl >= UV_WL_MIN) & (wl <= UV_WL_MAX)
    )
    ax_uv.plot(wl[m], flux[m], ls=ls, lw=1.5, color=color, label=label)

ax_uv.set_xlim(UV_WL_MIN, UV_WL_MAX)
ax_uv.set_yscale('log')
ax_uv.set_xlabel('Wavelength [nm]', fontsize=12, weight='bold')
ax_uv.set_ylabel('Spectral irradiance [mW m$^{-2}$ nm$^{-1}$]', fontsize=12, weight='bold')
ax_uv.set_title('TOA solar spectra (UV, 100–400 nm)', fontsize=13, weight='bold')
ax_uv.legend(loc='upper right', fontsize=8, frameon=False)
thick_axes(top=True)
plt.ylim(1e0, 1e2)

fig_stellar.suptitle(
    'Stellar input spectra used in 1-D photochemical models',
    fontsize=14, weight='bold', y=1.02,
)
plt.tight_layout()
plt.savefig(
    '/Users/gregcooke/python_output/Stellar_spectra_model_comparison.png',
    dpi=dpi,
    bbox_inches='tight',
)
print('Saved stellar spectra figure to python_output/Stellar_spectra_model_comparison.png')

# --- Figure: ratio to WACCM on a common interpolated wavelength grid ---
_wl_common = np.linspace(UV_WL_MIN, UV_WL_MAX, 601)


def _interp_flux_onto_grid(wl, flux, wl_grid):
    """Linear interpolation onto wl_grid (nm); NaN outside source range."""
    wl = np.asarray(wl, dtype=float)
    flux = np.asarray(flux, dtype=float)
    m = np.isfinite(wl) & np.isfinite(flux) & (flux > 0)
    if m.sum() < 2:
        return np.full(wl_grid.shape, np.nan)
    order = np.argsort(wl[m])
    wl_s = wl[m][order]
    flux_s = flux[m][order]
    out = np.interp(wl_grid, wl_s, flux_s, left=np.nan, right=np.nan)
    out[(wl_grid < wl_s.min()) | (wl_grid > wl_s.max())] = np.nan
    return out


_flux_waccm_grid = _interp_flux_onto_grid(wl_waccm, flux_waccm, _wl_common)
_waccm_ok = np.isfinite(_flux_waccm_grid) & (_flux_waccm_grid > 0)

_ratio_specs = [
    (wl_vulcan, flux_vulcan, 'VULCAN / WACCM', '#1f77b4', '-'),
    (wl_photo, flux_photo, 'Photochem / WACCM', '#2ca02c', '-'),
    (wl_kast_nm, flux_kast, 'Kasting photos.pdat / WACCM', '#d62728', '-'),
    (wl_faruv_nm, flux_faruv, 'Kasting faruv_sun / WACCM', '#d62728', '--'),
    (wl_atmos, flux_atmos, 'Atmos / WACCM', '#9467bd', '-'),
]

fig_ratio, ax_ratio = plt.subplots(figsize=(9, 5.5))

for wl, flux, label, color, ls in _ratio_specs:
    flux_on_grid = _interp_flux_onto_grid(wl, flux, _wl_common)
    ratio = flux_on_grid / _flux_waccm_grid
    m_plot = _waccm_ok & np.isfinite(ratio) & (ratio > 0)
    ax_ratio.plot(_wl_common[m_plot], ratio[m_plot], ls=ls, lw=1.5, color=color, label=label)

ax_ratio.axhline(1.0, color='#ff7f0e', ls=':', lw=1.5, label='WACCM (reference)')
ax_ratio.set_xlim(UV_WL_MIN, UV_WL_MAX)
#ax_ratio.set_yscale('log')
ax_ratio.set_xlabel('Wavelength [nm]', fontsize=12, weight='bold')
ax_ratio.set_ylabel('Flux ratio (model / WACCM6 SSI)', fontsize=12, weight='bold')
ax_ratio.set_title('Spectral irradiance relative to WACCM (100–400 nm)', fontsize=13, weight='bold')
ax_ratio.legend(loc='upper right', fontsize=8, frameon=False)
thick_axes(top=True)
plt.ylim(0.5,1.5)
plt.xlim(110, 250)
#plt.axvline(121.6, color = 'm', lw = 3)

fig_ratio.suptitle(
    'Stellar spectra vs WACCM6 CMIP6 (interpolated to common grid)',
    fontsize=14, weight='bold', y=1.02,
)
plt.tight_layout()
plt.savefig(
    '/Users/gregcooke/python_output/Stellar_spectra_ratio_to_WACCM.png',
    dpi=dpi,
    bbox_inches='tight',
)
print('Saved ratio figure to python_output/Stellar_spectra_ratio_to_WACCM.png')

#%%

metrics = []

for wl, flux, label, *_ in _ratio_specs:

    flux_grid = _interp_flux_onto_grid(wl, flux, _wl_common)

    m = (
        (_wl_common >= 180) &
        (_wl_common <= 250) &
        np.isfinite(flux_grid) &
        np.isfinite(_flux_waccm_grid) &
        (_flux_waccm_grid > 0)
    )

    ratio = flux_grid[m] / _flux_waccm_grid[m]

    mad = np.mean(np.abs(ratio - 1.0))

    rms = np.sqrt(np.mean((ratio - 1.0)**2))

    log_rms = np.sqrt(
        np.mean(
            (np.log10(flux_grid[m]) -
             np.log10(_flux_waccm_grid[m]))**2
        )
    )

    int_bias = (
        np.trapz(flux_grid[m], _wl_common[m]) -
        np.trapz(_flux_waccm_grid[m], _wl_common[m])
    ) / np.trapz(_flux_waccm_grid[m], _wl_common[m])

    metrics.append(
        (label, mad, rms, log_rms, int_bias)
    )

for row in metrics:
    print(row)

