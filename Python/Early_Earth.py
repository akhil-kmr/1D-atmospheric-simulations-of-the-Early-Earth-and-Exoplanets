#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 19:52:47 2023

@author: Greg Cooke (gjc53@cam.ac.uk)
"""

#%% imports
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from matplotlib import rcParams
from photochem import EvoAtmosphere

#%% 
'''
Instructions for user to go here
'''

#%% Other functions

def compute_ox_production(
    mech_file="zahnle_earth.yaml",
    settings_file="input/settings_100pc.yaml",
    flux_file="input/Sun_0.0Ga.txt",
    pt_file="100pc_PT_profile/100%PAL_Sun_0.0Ga.txt",
    atol=1e-23,
    verbose=0
):
    """
    Compute Ox production from O2 photolysis.

    Returns
    -------
    ox_prod_m3 : np.ndarray
        Ox production rate (molecules/m^3/s)
    pressure_hpa : np.ndarray
        Pressure profile (hPa)
    """

    # Initialize model
    pc = EvoAtmosphere(
        mech_file,
        settings_file,
        flux_file,
        pt_file
    )
    
    pc.var.verbose = verbose
    pc.var.atol = atol

    # Get production/loss and solution
    pl = pc.production_and_loss('O2', pc.wrk.usol)
    sol = pc.mole_fraction_dict()
    pressure_hpa = sol['pressure'] / 1e2  # Pa → hPa

    # Find indices for O2 photolysis reactions
    idx_o_o = -1
    idx_o_o1d = -1

    for i, rx in enumerate(pl.loss_rx):
        if rx == 'O2 + hv => O + O':
            idx_o_o = i
        elif rx == 'O2 + hv => O + O1D':
            idx_o_o1d = i

    # Optional safety check
    if idx_o_o == -1 or idx_o_o1d == -1:
        raise ValueError("Could not find required O2 photolysis reactions in network.")

    # Compute Ox production (cm^-3 s^-1 → m^-3 s^-1)
    ox_prod_cm3 = (2 * pl.loss[:, idx_o_o]) + (2 * pl.loss[:, idx_o_o1d])
    ox_prod_m3 = ox_prod_cm3 * 1e6

    return ox_prod_m3, pressure_hpa

P_path = '/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/'
def compute_o3loss(
    mech_file="zahnle_earth.yaml",
    settings_file="input/settings_100pc.yaml",
    flux_file="input/Sun_0.0Ga.txt",
    pt_file="100pc_PT_profile/100%PAL_Sun_0.0Ga.txt",
    atol=1e-23,
    verbose=0
):
    """
    Compute Ox production from O2 photolysis.

    Returns
    -------
    ox_prod_m3 : np.ndarray
        Ox production rate (molecules/m^3/s)
    pressure_hpa : np.ndarray
        Pressure profile (hPa)
    """

    # Initialize model
    pc = EvoAtmosphere(
        mech_file,
        settings_file,
        flux_file,
        pt_file
    )
    
    pc.var.verbose = verbose
    pc.var.atol = atol

    # Get production/loss and solution
    pl = pc.production_and_loss('O3', pc.wrk.usol)
    sol = pc.mole_fraction_dict()
    pressure_hpa = sol['pressure'] / 1e3  # → hPa

    # Find indices for O2 photolysis reactions
    idx_o_o3 = -1
    for i, rx in enumerate(pl.loss_rx):
        if rx == 'O + O3 => O2 + O2':
            idx_o_o3 = i
            print('yay')

    # Optional safety check
    if idx_o_o3 == -1:
        raise ValueError("Could not find required O2 photolysis reactions in network.")

    # Compute Ox production (cm^-3 s^-1 → m^-3 s^-1)
    o3loss_cm3 = pl.loss[:, idx_o_o3]
    o3_loss = o3loss_cm3 * 1e6

    return o3_loss, pressure_hpa

rcParams['font.weight'] = 'bold' 

#%% return thick lw on axes

def thick_axes(top = False, labelleft = True, right = True, direction = 'in'):
    # Accessing the axes object and setting linewidth
    plt.gca().spines['top'].set_linewidth(2)  # Top axis
    plt.gca().spines['bottom'].set_linewidth(2)  # Bottom axis
    plt.gca().spines['left'].set_linewidth(2)  # Left axis
    plt.gca().spines['right'].set_linewidth(2)  # Right axis
    plt.grid(False)
    plt.tick_params(which = 'major', axis = 'both', direction = direction, labelsize = 15, length = 6, width = 2, labelleft = labelleft, right = right, top = top)
    plt.tick_params(which = 'minor',axis = 'both', direction = direction, labelsize = 15, length = 3, width = 1, labelleft = labelleft, right = right, top = top)

#%% return subscript number or text
def sub(num):
    return r'$_{'+str(num)+'}$'

#%% return superscript number or text

def sup(num):
    return r'$^{'+str(num)+'}$'

#%% return specific LWAV
def LWAV_spec(DS, time = False, latitude = np.arange(0,96)):
    File =  "/Users/gregcooke/H_escape/b.e21.BWma1850.f19_g17.PC_b.SSPO.016.cam.h0.0320-0320.nc" #file name

    Gw_file = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

    gw = Gw_file.gw.values
    
    DS = DS.isel(lat=latitude)
    gw = gw[latitude]
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

dpi = 200

#%% Colour bar definition

def cbar(model, label = True, clabel = 'O'+ sub(3) +' column [DU]', orientation='vertical', tick = False, ticks = [2,3], shrink = 1):
    cbar = plt.colorbar(model, orientation = orientation, shrink = shrink)
    if (label == True):
        clabel = clabel
        cbar.set_label(clabel,size=15, weight = 'bold')
    cbar.ax.tick_params(labelsize=15)
    if (tick == True):
        cbar.set_ticks(ticks)
        
#%% Read in VULCAN files

import pickle

def read_pickle_file(file_path):
    try:
        with open(file_path, 'rb') as file:
            data = pickle.load(file)
            return data
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
        return None
    except Exception as e:
        print(f"An error occurred while reading the pickle file: {e}")
        return None

figure = plt.figure(figsize = (8,5))
lw = 2

plt.title('Ozone profiles', fontsize = 15, weight = 'bold')

def Read_O3_Run(file_path=''):
    DS = read_pickle_file('/Users/gregcooke/VIH_cases/output/'+file_path)
    spec = DS['variable']['species']
    return DS, spec

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
Ten_pc_V_45SZA, Ten_pc_spec_45SZA = Read_O3_Run(file_path='Earth_10pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Five_pc_V_60SZA, Five_pc_spec_60SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_60SZA_WPT_1rtol.vul')
Five_pc_V_482SZA, Five_pc_spec_482SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_48.2SZA_WPT_1rtol.vul')
Five_pc_V_45SZA, Five_pc_spec_45SZA = Read_O3_Run(file_path='Earth_5pc_o2_1e12s_45SZA_WPT_1rtol.vul')

One_pc_V_60SZA, One_pc_spec_60SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_60SZA_WPT_1rtol.vul')
One_pc_V_482SZA, One_pc_spec_482SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_48.2SZA_WPT_1rtol.vul')
One_pc_V_45SZA, One_pc_spec_45SZA = Read_O3_Run(file_path='Earth_1pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Z5_pc_V_60SZA, Z5_pc_spec_60SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_60SZA_WPT_1rtol.vul')
Z5_pc_V_482SZA, Z5_pc_spec_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_WPT_1rtol.vul')
Z5_pc_V_45SZA, Z5_pc_spec_45SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_45SZA_WPT_1rtol.vul')

Z1_pc_V_60SZA, Z1_pc_spec_60SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_60SZA_WPT_1rtol.vul')
Z1_pc_V_482SZA, Z1_pc_spec_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_WPT_1rtol.vul')
Z1_pc_V_45SZA, Z1_pc_spec_45SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_45SZA_WPT_1rtol.vul')

#%% VULCAN cases with different methane
#One50_pc_V_10xCH4_482SZA, One50_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_150pc_o2_1e12s_48.2SZA_10xCH4_WPT_1rtol.vul')
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
'''
Z5_pc_V_10xCH4f_482SZA, Z5_pc_spec_10xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_10xCH4flux_WPT_1rtol.vul')
Z5_pc_V_5xCH4f_482SZA, Z5_pc_spec_5xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_5xCH4flux_WPT_1rtol.vul')
Z5_pc_V_1xCH4f_482SZA, Z5_pc_spec_1xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_1xCH4flux_WPT_1rtol.vul')
Z5_pc_V_05xCH4f_482SZA, Z5_pc_spec_05xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_0.5xCH4flux_WPT_1rtol.vul')
Z5_pc_V_01xCH4f_482SZA, Z5_pc_spec_01xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.5pc_o2_1e12s_48.2SZA_0.1xCH4flux_WPT_1rtol.vul')
'''
Z1_pc_V_10xCH4_482SZA, Z1_pc_spec_10xCH4_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_10xCH4_WPT_1rtol.vul')
Z1_pc_V_5xCH4_482SZA, Z1_pc_spec_5xCH4_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_5xCH4_WPT_1rtol.vul')
Z1_pc_V_05xCH4_482SZA, Z1_pc_spec_05xCH4_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_0.5xCH4_WPT_1rtol.vul')
Z1_pc_V_01xCH4_482SZA, Z1_pc_spec_01xCH4_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_0.1xCH4_WPT_1rtol.vul')
'''
Z1_pc_V_10xCH4f_482SZA, Z1_pc_spec_10xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_10xCH4flux_WPT_1rtol.vul')
Z1_pc_V_5xCH4f_482SZA, Z1_pc_spec_5xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_5xCH4flux_WPT_1rtol.vul')
Z1_pc_V_1xCH4f_482SZA, Z1_pc_spec_1xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_1xCH4flux_WPT_1rtol.vul')
Z1_pc_V_05xCH4f_482SZA, Z1_pc_spec_05xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_0.5xCH4flux_WPT_1rtol.vul')
Z1_pc_V_01xCH4f_482SZA, Z1_pc_spec_01xCH4f_482SZA = Read_O3_Run(file_path='Earth_0.1pc_o2_1e12s_48.2SZA_0.1xCH4flux_WPT_1rtol.vul')
'''
#%% Cross sections and checks

import pandas as pd
import io

file_path = '/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Kasting_1D_model/100pc/SZA_48.2/OUTPUT_PLOT.dat'

def read_kasting_file(path):
    with open(path, 'r') as f:
        lines = f.readlines()

    blocks = {}
    current_lines = []
    current_header = None

    for line in lines:
        # Detect the start of the Wavelength/Flux block
        if "Range (A)" in line:
            current_header = "flux_block"
            headers = ["Int", "Range_A", "Flux", "Ozone1", "Ozone2", "CH3CHO", "CH3COCH3"]
            continue
        
        # Detect the start of the Assorted Cross Sections block
        if "Assorted cross sections" in line:
            current_header = "cross_sections"
            # Move down to find the actual column names
            continue
        
        if current_header == "cross_sections" and "Int #" in line:
            headers = ["Int", "O2", "H2O", "CO2", "HO2", "N2O", "HCl", "CCl2F2", "CHClF2", "CH3Cl", "CH3CCL3"]
            continue

        # If we hit the separator "====", start collecting data
        if "====" in line and current_header:
            current_lines = []
            continue

        # Collect data lines (must start with a number)
        strip_line = line.strip()
        if strip_line and strip_line[0].isdigit():
            current_lines.append(strip_line)
        # If we hit a blank line or new text, save the block
        elif current_lines and not strip_line:
            df = pd.read_csv(io.StringIO("\n".join(current_lines)), sep='\s+', names=headers)
            blocks[current_header] = df
            current_header = None
            current_lines = []

    # Catch the last block if file doesn't end in blank line
    if current_lines and current_header:
        df = pd.read_csv(io.StringIO("\n".join(current_lines)), sep='\s+', names=headers)
        blocks[current_header] = df

    return blocks
'''
# Execute
data_blocks = read_kasting_file(file_path)

# Extract your two main tables
df_flux = data_blocks['flux_block']
df_cross = data_blocks['cross_sections']

# Optional: Split the 'Range_A' into min and max wavelengths
df_flux[['Wave_Min', 'Wave_Max']] = df_flux['Range_A'].str.split('-', expand=True).astype(float)

print("--- Wavelength & Flux Table ---")
print(df_flux[['Int', 'Wave_Min', 'Wave_Max', 'Flux']].head())

print("\n--- Cross Sections Table ---")
print(df_cross[['Int', 'O2', 'H2O', 'N2O']].head())

plt.figure()
plt.plot(PI_V_482SZA['variable']['bins'], PI_V_482SZA['variable']['cross_J']['N2O',1])
plt.plot((df_flux['Wave_Min']+df_flux['Wave_Max'])[:35]/20, df_cross['N2O'])
plt.xlim(1, 300)
plt.yscale('log')

plt.figure()
plt.plot(PI_V_482SZA['variable']['bins'], PI_V_482SZA['variable']['cross_J']['O2',1]+PI_V_482SZA['variable']['cross_J']['O2',2], color = 'k')
plt.plot(PI_V_482SZA['variable']['bins'], PI_V_482SZA['variable']['cross_J']['O3',1]+PI_V_482SZA['variable']['cross_J']['O3',2], color = 'm')
plt.ylim(1e-25, 1e-16)
#plt.plot((df_flux['Wave_Min']+df_flux['Wave_Max'])[:35]/20, df_cross['O2'], color = 'teal')
plt.xlim(1, 300)
plt.yscale('log')

plt.figure()
plt.plot(PI_V_482SZA['variable']['bins'], PI_V_482SZA['variable']['cross_J']['H2O',1], color = 'm')
plt.plot(PI_V_482SZA['variable']['bins'], PI_V_482SZA['variable']['cross_J']['H2O',2], color = 'm')
plt.plot((df_flux['Wave_Min']+df_flux['Wave_Max'])[:35]/20, df_cross['H2O'], color = 'teal')
plt.xlim(1, 300)
plt.yscale('log')

plt.figure()
plt.plot(PI_V_482SZA['variable']['bins'], PI_V_482SZA['variable']['cross_J']['CO2',1], color = 'm')
plt.plot(PI_V_482SZA['variable']['bins'], PI_V_482SZA['variable']['cross_J']['CO2',2], color = 'm')
plt.plot((df_flux['Wave_Min']+df_flux['Wave_Max'])[:35]/20, df_cross['CO2'], color = 'teal')
plt.xlim(1, 300)
plt.yscale('log')


plt.figure()
plt.plot(PI_V_482SZA['variable']['bins'], PI_V_482SZA['variable']['cross_J']['HO2',1], color = 'm')
plt.plot((df_flux['Wave_Min']+df_flux['Wave_Max'])[:35]/20, df_cross['HO2'], color = 'teal')
plt.xlim(1, 300)
plt.yscale('log')
'''
#%% Read in Kasting simulations

import re

K_150pc_J = pd.read_csv(
    "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Kasting_1D_model/150pc/SZA_48.2/OUTPUT_PLOT.dat",
    delim_whitespace=True, engine="python")
p_150pc = K_150pc_J['PRESS']/1000
z = K_150pc_J['Z']/1000

K_100pc_J = pd.read_csv(
    "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Kasting_1D_model/100pc/SZA_48.2/OUTPUT_PLOT.dat",
    delim_whitespace=True, engine="python")

p_100pc = K_100pc_J['PRESS']/1000
z = K_100pc_J['Z']/1000

K_10pc_J = pd.read_csv(
    "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Kasting_1D_model/10pc/SZA_48.2/OUTPUT_PLOT.dat",
    delim_whitespace=True, engine="python")
p_10pc = K_10pc_J['PRESS']/1000
z = K_10pc_J['Z']/1000

K_1pc_J = pd.read_csv(
    "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Kasting_1D_model/1pc/SZA_48.2/OUTPUT_PLOT.dat",
    delim_whitespace=True, engine="python")
p_1pc = K_1pc_J['PRESS']/1000
z = K_1pc_J['Z']/1000

K_01pc_J = pd.read_csv(
    "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Kasting_1D_model/0.1pc/SZA_48.2/OUTPUT_PLOT.dat",
    delim_whitespace=True, engine="python")
p_01pc = K_01pc_J['PRESS']/1000
z = K_01pc_J['Z']/1000

#%% J rates plot

def n_dens(ds, var = 'H2O'):
    k = 1.381e-23
    ndens = ds[var]*ds.lev*100/(ds.T*k)
    return ndens

def Kasting_HOX(DS):
    
    # Specify the columns involved in the calculation
    cols = ['FOH', 'FHO2', 'FH2O2', 'FH']
    
    # 1. Convert columns to numeric: malformed strings become NaN
    # 2. fillna(0): Replace those NaNs with 0.0
    temp_df = DS[cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    # Calculate the sum using the cleaned numeric values
    HOX = temp_df['FH'] + temp_df['FHO2'] + temp_df['FH2O2'] + temp_df['FOH']
    
    return HOX

def Kasting_NOX(DS):
    # Select the columns we need
    cols = ['FNO', 'FNO2']
    
    # 1. Convert to numeric: invalid strings (like 8.966-223) become NaN
    # 2. fillna(0): Replace those NaNs with 0
    temp_df = DS[cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    # Sum the cleaned numeric values
    return temp_df['FNO'] + temp_df['FNO2']


def J_rates_K(data):
    po2 = data['PO2']
    po2D = data['PO2D']
    for i in range(len(po2D)):
        try:
            po2D[i] = float(po2D[i])  
        except:
            po2D[i] = float(0)  
    JO2 = po2+po2D
    return 2*JO2 # factor for 2*O

def prox_ox_K(data):
    return J_rates_K(data) * data['NumO2'] * 1e6 # factor for density

def J_rates_W(DS):
    JO2 = (DS.jo2_a+DS.jo2_b)
    return 2*JO2

def prox_ox_W(DS):
    return J_rates_W(DS) * n_dens(DS, 'O2')

def J_rates_V(data):
    JO2 = (
        data['variable']['J_sp']['O2',0]
    )
    return 2*JO2

def prox_ox_V(data, data_spec):
    return J_rates_V(data) * data['variable']['y'][:,data_spec.index('O2')] * 1e6


#%% Read in Atmos files

import pandas as pd


def fix_exponents(line):
    # Fix cases like 5.004-232 → 5.004E-232
    return re.sub(r'(\d\.\d+)([-+]\d+)', r'\1E\2', line)

def Read_Atmos(file_path):

    # Read and fix file
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Fix all lines except header
    header = lines[0]
    data_lines = [fix_exponents(line) for line in lines[1:]]
    
    # Combine back
    fixed_text = header + ''.join(data_lines)
    
    # Read into pandas
    from io import StringIO
    df = pd.read_csv(StringIO(fixed_text), delim_whitespace=True)
    return df

def Atmos_dens(DS):
    k = 1.381e-23
    dens = DS["PRESS"]*1e5 / (k*DS["TEMP"])
    return dens

def Atmos_HOX(DS):
    # Specify the columns involved in the calculation
    cols = ['OH', 'HO2', 'H2O2', 'H']
    
    # 1. Convert columns to numeric: malformed strings become NaN
    # 2. fillna(0): Replace those NaNs with 0.0
    temp_df = DS[cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    # Calculate the sum using the cleaned numeric values
    HOX = temp_df['H'] + temp_df['HO2'] + temp_df['H2O2'] + temp_df['OH']
    
    return HOX
def Atmos_NOX(DS):
    # Select the columns we need
    cols = ['N', 'NO', 'NO2']
    
    # 1. Convert to numeric: invalid strings (like 8.966-223) become NaN
    # 2. fillna(0): Replace those NaNs with 0
    temp_df = DS[cols].apply(pd.to_numeric, errors='coerce').fillna(0)
    
    # Sum the cleaned numeric values
    return temp_df['N'] + temp_df['NO'] + temp_df['NO2']


'''
Now read in horribly formatted photolysis data
'''

import pandas as pd

file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/100pc/SZA_48.2/PTZ_mixingratios_out.dist"
Atmos_100pc = pd.read_csv(file_path, delim_whitespace=True)
file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/10pc/SZA_48.2/PTZ_mixingratios_out.dist"
Atmos_10pc = pd.read_csv(file_path, delim_whitespace=True)
file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/1pc/SZA_48.2/PTZ_mixingratios_out.dist"
Atmos_1pc = pd.read_csv(file_path, delim_whitespace=True)
file_path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/0.1pc/SZA_48.2/PTZ_mixingratios_out.dist"
Atmos_01pc = pd.read_csv(file_path, delim_whitespace=True)
'''
def Atmos_photo(file_path):
    
    def fix_exponents(s):
        # Fix cases like 5.44-173 → 5.44E-173
        return re.sub(r'(\d\.\d+)([-+]\d+)', r'\1E\2', s)
    
    def safe_float(x):
        try:
            return float(x)
        except ValueError:
            return 0.0
    
    with open(file_path) as f:
        lines = f.readlines()
    
    # --- find photolysis section ---
    start_idx = None
    for i, line in enumerate(lines):
        if "PHOTOLYSIS RATES" in line:
            start_idx = i
            break
    
    header = lines[start_idx + 2].split()
    
    # --- read data block ---
    data = []
    ncols = len(header)
    
    for line in lines[start_idx + 3:]:
        s = line.strip()
        
        if s == "":
            break
        
        parts = fix_exponents(s).split()
        
        # stop if structure changes
        if len(parts) != ncols:
            break
        
        # convert safely
        row = [safe_float(x) for x in parts]
        data.append(row)
    
    df = pd.DataFrame(data, columns=header)
    O2_Photo = df["PO2_O1D"]+df["PO2_O3P"]
    return df, df["Z"], O2_Photo
'''
def Atmos_photo(file_path):
    
    df = pd.read_csv(file_path, names = ['Z', 'PO2_1', "PO2_2"], delim_whitespace=True)
    JO2 = 2*(df["PO2_1"]+df["PO2_2"])
    return df["Z"], JO2

def Atmos_dens(DS):
    k = 1.381e-23
    dens = DS["PRESS"]*1e5/(k*DS["TEMP"])
    return dens

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
import pandas as pd
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
Photo_PI_01CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_0.1x_methane.txt", delim_whitespace=True)
Photo_PI_05CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_0.5x_methane.txt", delim_whitespace=True)
Photo_PI_1CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_1x_methane.txt", delim_whitespace=True)
Photo_PI_2CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_2x_methane.txt", delim_whitespace=True)
Photo_PI_3CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_3x_methane.txt", delim_whitespace=True)
Photo_PI_4CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_4x_methane.txt", delim_whitespace=True)
Photo_PI_5CH4 = pd.read_csv(path+"../Methane_Perturbations/Earth_100pc_48.2_5x_methane.txt", delim_whitespace=True)


path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/10pc/Earth_10pc_48.2.txt"
Photo_10pc = pd.read_csv(path, delim_whitespace=True)



path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/1pc/Earth_1pc_48.2.txt"
Photo_1pc = pd.read_csv(path, delim_whitespace=True)



path = "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/0.1pc/Earth_0.1pc_48.2.txt"
Photo_01pc = pd.read_csv(path, delim_whitespace=True)

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

PI_color = 'k'; 
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

#%% define ozone column calculation function

def O3_col_trapz(ds, time=False, O2_mr=0.21, lon=True, g=9.81):
    """
    Calculate ozone column in Dobson Units using explicit loops
    and trapezoidal integration over hybrid-pressure layers.
    """
    
    g = g  # gravity acceleration m/s^2
    to_DU = 1./2.6867e20  # convert mol/m^2 to DU
    N2_mr = 1-O2_mr
    m_avg = ((28*N2_mr) + (O2_mr*32)) * 1.661e-27  # kg per molecule? (used here for consistency)
    
    # O3 variable
    o3 = ds['O3'].values
    
    # Hybrid pressure variables
    ps = ds['PS'].values
    p0 = ds['P0'].values
    hyai = ds['hyai'].values
    hybi = ds['hybi'].values
    
    nt = ds['time'].size if 'time' in ds.dims else 0
    ny = ds['lat'].size
    nz = ds['lev'].size
    nx = ds['lon'].size
    
    # Initialize layer thickness array
    if nt > 0:
        dp = np.zeros((nz-1, ny, nx), dtype=np.float32)

    else:
        dp = np.zeros((nz-1, ny, nx), dtype=np.float32)

    
    # Compute layer thickness dp using trapezoid (average between interfaces)
    if nt > 0:
        for i in range(nt):
            for j in range(ny):
                for k in range(nx):
                    # Compute interface pressures
                    p_interfaces = hyai * p0 + hybi * ps[i, j, k]
                    for iz in range(nz-1):
                        dp[iz, j, k] = p_interfaces[iz+1] - p_interfaces[iz]

    else:
        for j in range(ny):
            for k in range(nx):
                p_interfaces = hyai * p0 + hybi * ps[j, k]
                for iz in range(nz-1):
                    dp[iz, j, k] = p_interfaces[iz+1] - p_interfaces[iz]

    
    # Apply trapezoidal integration: multiply layer thickness by mean of O3 at layer
    if nt > 0:
        o3_col = np.zeros((nt, ny, nx), dtype=np.float32)
        for i in range(nt):
            for j in range(ny):
                for k in range(nx):
                    # Trapezoid: sum over vertical layers
                    o3_col[i, j, k] = np.sum((o3[i, :-1, j, k] + o3[i, 1:, j, k])/2 * dp[i, :, j, k])
    else:
        o3_col = np.zeros((ny, nx), dtype=np.float32)
        for j in range(ny):
            for k in range(nx):
                o3_col[j, k] = np.sum((o3[:-1, j, k] + o3[1:, j, k])/2 * dp[:, j, k])
    
    # Divide by gravity to get column (mol/m^2)
    o3_col /= g
    
    # Convert to Dobson Units
    o3_col *= to_DU
    
    # Optional averaging over lon/time
    if not lon:
        o3_col = o3_col.mean(axis=-1)  # longitude
    if not time and nt > 0:
        o3_col = o3_col.mean(axis=0)  # time
    
    print("O3 columns calculated using trapezoid")
    return o3_col/m_avg 

def O3_col(ds, time = False, O2_mr = 0.21, lon = True, g = 9.81):
    
    g = g # gravity acceleration (cm/sec2)
    to_DU = 1./2.6867e20  # convert colunm density / cm^2 to Dobson units
    N2_mr = 1-O2_mr #nitrogen mixing ratio
    m_avg = ((28*N2_mr)+(O2_mr*32))*1.661e-27 # average weight of atmosphere
    
    #calcualte column integrals 
    o3 = ds['O3']
    
    # variables to calcualate hybrid pressure on model interfaces
    ps = ds['PS'].values # surface pressure
    p0 = ds['P0'].values
    hyai = ds['hyai'].values
    hybi = ds['hybi'].values
    
    #col = np.ndarray(ps.shape, np.float32)
    dp = np.ndarray(o3.shape, np.float32)
    
    nt = 1
    try:
        nt = ds['time'].size
    except:
        nt = 0
    ny = ds['lat'].size
    nz = ds['lev'].size
    nx = ds['lon'].size
    
    press = np.ndarray(nz, np.float32)
    if (nt > 0):
        for i in range(nt):
            for j in range(ny):
                for k in range(nx):
                    press = hyai * p0 + hybi * ps[i,j,k]
                    for iz in range(nz):
                        dp[i,iz,j,k] = (press[iz+1]-press[iz])
    else:
        for j in range(ny):
            for k in range(nx):
                press = hyai * p0 + hybi * ps[j,k]
                for iz in range(nz):
                    dp[iz,j,k] = (press[iz+1]-press[iz])
                    
    #calcualte column integrals 
    o3 = o3 * dp / m_avg
    
    if (nt > 0):
        for i in range(nt):
            for j in range(96):
                for k in range(144):
                    o3[i,:,j,k] = o3[i,:,j,k]/g    
    else:
        for j in range(96):
            for k in range(144):
                o3[:,j,k] = o3[:,j,k]/g
    
    if (nt > 0):
        o3col = np.sum(o3, axis = 1) * to_DU # convert to DU
    else:
        o3col = np.sum(o3, axis = 0) * to_DU # convert to DU
    if (lon == False):
        o3col  = o3col.mean(dim='lon')
    if (time == False):
        try:
            o3col  = o3col.mean(dim='time')
        except:
            o3col  = o3col
    print('O3 columns calculated')
    return o3col

#%%
'''
# Pre-industrial baseliine run (PI; 100% PAL of O2)
File =  "/Users/gregcooke/H_escape/b.e21.BWma1850.f19_g17.baseline.cam.h0.0013-0017.nc" #file name
Pre_pc_h0 = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

col = O3_col(Pre_pc_h0.mean(dim='time'), time = False, O2_mr = 0.21, lon = True, g = 9.81)
col_trapz = O3_col_trapz(Pre_pc_h0.mean(dim='time'), time = False, O2_mr = 0.21, lon = True, g = 9.81)

col_zm = col_trapz.mean(axis=1)   # shape (96,)
global_mean = np.sum(gw * col_zm) / 2.0

col_zm = O3_col.mean(axis=1)   # shape (96,)
global_mean = np.sum(gw * col_zm) / 2.0
'''
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
import cartopy.crs as ccrs  

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

#%%


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


def convert_photons_to_energy(wav, flux):
    
    # Physical constants
    h = 6.626e-34  # Planck constant (J·s)
    c = 2.998e8    # speed of light (m/s)

    # Convert wavelength to meters
    wav = wav * 1e-10

    # Conversion: photons/cm²/s → W/m²/nm
    irradiance = flux * (h * c / wav) * 1e4
    
    return irradiance

FarUV_SFX_3D_arr_corrected = convert_photons_to_energy(FarUV_3D_Wav_arr/10, FarUV_3D_SFX_arr)
FarUV_SFX_arr_corrected = convert_photons_to_energy(FarUV_Wav_arr/10, FarUV_SFX_arr)
Flux_arr_corrected = convert_photons_to_energy(wav/10, flux_arr)

import pandas as pd
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
plt.plot(wav/10, Flux_arr_corrected/Kasting_factor, color = 'k', lw = 2, label = 'Kasting 1D model')
#plt.plot(wav/10, 3*Flux_arr_corrected/Kasting_factor, color = 'k', lw = 2, label = 'Kasting 1D model')
#plt.plot(FarUV_Wav_arr/10, FarUV_SFX_arr_corrected/Kasting_factor, color = 'k', lw = 2)#, label = 'Kasting 1D model')
#plt.plot(FarUV_3D_Wav_arr/10, FarUV_SFX_3D_arr_corrected/Kasting_factor, color = 'k', lw = 2)#, label = 'Kasting 1D model')

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
plt.plot(Sun_0_0['Wav'], Sun_0_0['Flux'], color = 'k', lw = 2, label = 'Sun 0.0 Ga')
plt.plot(Sun_0_5['Wav'], Sun_0_5['Flux'], color = 'grey', lw = 2, label = 'Sun 0.5 Ga')
plt.plot(Sun_1_0['Wav'], Sun_1_0['Flux'], color = 'b', lw = 2, label = 'Sun 1.0 Ga')
plt.plot(Sun_2_0['Wav'], Sun_2_0['Flux'], color = 'cyan', lw = 2, label = 'Sun 2.0 Ga')
plt.plot(Sun_2_4['Wav'], Sun_2_4['Flux'], color = 'g', lw = 2, label = 'Sun 2.4 Ga')

#plt.yscale('log')
plt.xlim(10,1000)
#plt.ylim(1e-3)
#%% Cross sections

path = '/Users/gregcooke/K2_18_XS/'
v_path = '/Users/gregcooke/VULCAN/thermo/photo_cross/'
csv_a_H2O = pd.read_csv(path+'H2O.XS.dat', delim_whitespace=True, skiprows=4, names=['Wav', 'XS'])
csv_a_CO2 = pd.read_csv(path+'CO2.XS.dat', delim_whitespace=True, skiprows=4, names=['Wav', 'XS'])

csv_v_H2O = pd.read_csv(v_path+'H2O/H2O_cross.csv', skiprows=1, names=['Wav', 'XS', 'XS1', 'XS2'])
csv_v_CO2 = pd.read_csv(v_path+'CO2/CO2_cross.csv', skiprows=1, names=['Wav', 'XS', 'XS1', 'XS2'])

wav = (range_min+range_max)/2
plt.figure(figsize = (9,5))
plt.plot(wav/10, ozone1_arr)
plt.plot(wav/10, ozone2_arr)
plt.plot((wav/10)[:35], o2_arr, color = 'cyan', label = 'O2')
plt.plot((wav/10)[:35], n2o_arr, color = 'purple', label = 'N2O')
plt.plot((wav/10)[:35], co2_arr, color = 'g', lw = 2, ls = '-', label = 'CO2')
plt.plot(csv_a_CO2['Wav']/10, csv_a_CO2['XS'], color = 'g', lw = 2, ls = '--')
plt.plot(csv_v_CO2['Wav'], csv_v_CO2['XS'], color = 'g', lw = 2, ls = ':')
plt.plot((wav/10)[:35], h2o_arr, color = 'b', lw = 2, ls = '-', label = 'H2O')
plt.plot(csv_a_H2O['Wav']/10, csv_a_H2O['XS'], color = 'b', lw = 2, ls = '--')
plt.plot(csv_v_H2O['Wav'], csv_v_H2O['XS'], color = 'b', lw = 2, ls = ':')
thick_axes()
plt.ylim(1e-27,1e-16)
plt.legend(loc = 0, fontsize = 15, frameon = False)
plt.xlim(120,250)
plt.yscale('log')
plt.ylabel('Cross section', fontsize = 15, weight = 'bold')
plt.xlabel('Wavelength [nm]', fontsize = 15, weight = 'bold')

data = read_pickle_file(file_path); ls = '-'
ls1 = ':'
plt.figure(figsize = (9,5))
#species = data['variable']['species']
plt.plot(Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('O3')], Z1_pc_V_482SZA['atm']['pco']/1e6, lw = lw, color = Zero1_pc_color, ls = ls, label = 'O'+sub(3))
plt.plot(Z5_pc_V_482SZA['variable']['ymix'][:,Z5_pc_spec_482SZA.index('O3')], Z5_pc_V_482SZA['atm']['pco']/1e6, lw = lw, color = Zero5_pc_color, ls = ls, label = 'O'+sub(3))
plt.plot(One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('O3')], One_pc_V_482SZA['atm']['pco']/1e6, lw = lw, color = One_pc_color, ls = ls, label = 'O'+sub(3))
plt.plot(Five_pc_V_482SZA['variable']['ymix'][:,Five_pc_spec_482SZA.index('O3')], Five_pc_V_482SZA['atm']['pco']/1e6, lw = lw, color = Five_pc_color, ls = ls, label = 'O'+sub(3))
plt.plot(Ten_pc_V_482SZA['variable']['ymix'][:,Ten_pc_spec_482SZA.index('O3')], Ten_pc_V_482SZA['atm']['pco']/1e6, lw = lw, color = Ten_pc_color, ls = ls, label = 'O'+sub(3))
plt.plot(Fifty_pc_V_482SZA['variable']['ymix'][:,Fifty_pc_spec_482SZA.index('O3')], Fifty_pc_V_482SZA['atm']['pco']/1e6, lw = lw, color = Fifty_pc_color, ls = ls, label = 'O'+sub(3))
plt.plot(PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('O3')], PI_V_482SZA['atm']['pco']/1e6, lw = lw, color = 'k', ls = ls, label = 'O'+sub(3))
plt.plot(One50_pc_V_482SZA['variable']['ymix'][:,One50_pc_spec_482SZA.index('O3')], PI_V_482SZA['atm']['pco']/1e6, lw = lw, color = 'grey', ls = ls, label = 'O'+sub(3))
plt.ylim(1, 1e-8)
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-12, 1e-4)
thick_axes()


#%% 1D comparisons

def P_dens(DS):
    k = 1.381e-23
    dens = DS['press'] * 1e5 / (k*DS['temp'])
    return dens

def P_NOX(DS):
    N = DS['N']
    NO = DS['NO']
    NO2 = DS['NO2']
    NOX = N+NO+NO2
    return NOX

def P_HOX(DS):
    H = DS['H']
    OH = DS['OH']
    HO2 = DS['HO2']
    H2O2 = DS['H2O2']
    HOX = H+OH+HO2+H2O2
    return HOX

def calculate_column_height_method(df):
    # Convert altitude from km to cm for unit consistency
    # 1 km = 100,000 cm
    z_cm = df['alt'].values * 100000 
    
    # Calculate Ozone number density (molecules/cm^3)
    # Ensure mixing ratio is in absolute units (not ppm/ppb)
    o3_num_den = df['O3'].values * df['den'].values
    
    # Integrate using the trapezoidal rule
    # This sums the area under the ozone density vs. height curve
    total_molecules = np.trapz(o3_num_den, z_cm)
    
    # Convert to Dobson Units
    column_du = total_molecules / 2.687e16
    return column_du

def calculate_column_pressure_method(df):
    # Constants
    g = 9.80665          # m/s^2
    m_air = 4.81e-26     # kg (approx mass of one air molecule)
    
    # Convert pressure from bar to Pascals
    press_pa = df['press'].values * 1e5
    
    # Mixing ratio
    o3_mixing = df['O3'].values
    
    # Integrate mixing ratio over pressure
    # Note: we use negative np.trapz because pressure decreases with index/altitude
    integral = -np.trapz(o3_mixing, press_pa)
    
    # Calculate total molecules per square meter
    total_molecules_m2 = integral / (g * m_air)
    
    # Convert m^-2 to cm^-2 (divide by 10,000)
    total_molecules_cm2 = total_molecules_m2 / 10000
    
    # Convert to Dobson Units
    column_du = total_molecules_cm2 / 2.687e16
    return column_du

print(calculate_column_pressure_method(Photo_PI))
print(calculate_column_height_method(Photo_PI))

P_col_150pc = calculate_column_height_method(Photo_150pc)
P_col_PI = calculate_column_height_method(Photo_PI)
P_col_10pc = calculate_column_height_method(Photo_10pc)
P_col_1pc = calculate_column_height_method(Photo_1pc)
P_col_01pc = calculate_column_height_method(Photo_01pc)

#%%

def method_A(df):
    ozone_index = df['O3']
    #  no. density (molecules/cm^3)
    ozone_density = df['den']*ozone_index
    altitudes = df['alt']
    # density over alt to get column (molecules/cm^2)
    column_density = np.trapz(ozone_density, altitudes)
    dobson_units = column_density / 2.6867e16 # Convert to DU
    print(f"Ozone column: {dobson_units:.2f} Dobson Units")
    
    return 1e5*dobson_units

def method_W(df):
    #wogans method
    ozone_index = df['O3']
    ozone_density = df['den']*ozone_index
    dz = df['alt'][1] - df['alt'][0]
    column_density = np.sum(ozone_density * dz)
    dobson_units = column_density / 2.6867e16 # Convert to DU
    print(f"Ozone column Wogan: {dobson_units:.2f} Dobson Units")
    
    return 1e5*dobson_units


def o3_trapz(df):
    ozone_density = df['O3']*df['den']*1e6 #retrieving number density (molecules/cm3)
    altitudes = df['alt']*1e3
    column_density = np.trapz(ozone_density, altitudes) #trapezoidal rule approximation to integrate the number density
    dobson_units = column_density / 2.6867e20  #loschmidt constant
    return dobson_units


def o3_mass_conserving_trapz(df, g=9.81, MO2=0.21):
    # mean molecular mass of dry air [kg]
    m_air = ((1 - MO2) * 28 + MO2 * 32) * 1.66054e-27

    # sort so highest pressure first
    df = df.sort_values('press', ascending=False)

    chi = df['O3'].to_numpy()                 # volume mixing ratio
    p = df['press'].to_numpy() * 1e5          # pressure [Pa]

    # Trapezoidal integration in pressure coordinates
    # Integral chi dp ≈ sum 0.5 * (chi[i] + chi[i+1]) * (p[i] - p[i+1])
    integral = np.trapz(chi, p)

    # mass-conserving column
    N_O3 = integral / (m_air * g)

    return N_O3 / 2.6867e20  # Dobson units


def o3_mass_conserving(df, g=9.81, MO2=0.21):
    # mean molecular mass of dry air [kg]
    m_air = ((1-MO2)*28 + MO2*32)*1.66054e-27

    # sort so highest pressure first
    df = df.sort_values('press', ascending=False)

    chi = df['O3'].to_numpy()
    p = df['press'].to_numpy()*1e5

    dp = -np.diff(p)  # pressure difference per layer

    # mass-conserving column
    N_O3 = np.sum(chi[:-1] * dp) / (m_air * g)
    
    # See Salby 1996

    return N_O3 / 2.6867e20  # Dobson units

def compare_columns(df, tol=0.05):
    # compute columns both ways
    du_pressure = o3_mass_conserving(df)
    # height-based trapezoidal with correct units
    ozone_density = df['O3'] * df['den']*1e6
    heights = df['alt'] * 1e3  # km→cm
    du_height = np.trapz(ozone_density, heights) / 2.6867e20

    # percentage difference
    pct_diff = 100*(du_height - du_pressure)/du_pressure

    warning = abs(pct_diff) > 100*tol
    return {
        'pressure_DU': du_pressure,
        'height_DU': du_height,
        'pct_diff': pct_diff,
        'warning': warning
    }

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
def layerwise_truncation_error(df, g=9.81, MO2=0.21):
    m_air = ((1-MO2)*28 + MO2*32)*1.66054e-27  # kg

    df_sorted = df.sort_values('press', ascending=False)
    
    chi = df_sorted['O3'].to_numpy()
    n_m3 = chi* (df_sorted['den'].to_numpy() * 1e6 )   # cm3 → m3
    z_m = df_sorted['alt'].to_numpy() * 1e3    # km → m
    p = df_sorted['press'].to_numpy() * 1e5    # hPa → Pa

    dp = -np.diff(p)
    dz = np.diff(z_m)

    # Mass-conserving layer
    N_O3_layer_pressure = chi[:-1] * dp / (g * m_air)
    
    # Height-trapezoid layer
    N_O3_layer_height = 0.5*(n_m3[:-1]+n_m3[1:]) * dz

    # Convert to Dobson Units
    to_DU = 1/2.6867e20
    N_O3_layer_pressure_DU = N_O3_layer_pressure * to_DU
    print(np.sum(N_O3_layer_pressure_DU))
    N_O3_layer_height_DU = N_O3_layer_height * to_DU
    print(np.sum(N_O3_layer_height_DU))

    pct_error = 100*(N_O3_layer_height_DU - N_O3_layer_pressure_DU) / N_O3_layer_pressure_DU

    return {
        'layer_pressure_DU': N_O3_layer_pressure_DU,
        'layer_height_DU': N_O3_layer_height_DU,
        'pct_error': pct_error,
        'z_mid': 0.5*(z_m[:-1]+z_m[1:])/1e3  # km
    }

# Example usage
errors = layerwise_truncation_error(P_new_WACCM)

import matplotlib.pyplot as plt

plt.figure()
plt.plot(errors['pct_error'], errors['z_mid'])
plt.xlabel('Layerwise % error (height vs pressure)')
plt.ylabel('Altitude (km)')
plt.grid(True)
plt.show()

def total_pct_error(errors):
    """
    Compute the weighted total percentage error of height-based vs. pressure-based column.
    Each layer is weighted by its pressure-column contribution.
    """
    layer_pressure = errors['layer_pressure_DU']
    layer_height = errors['layer_height_DU']
    
    # Fractional weight of each layer in total column
    weights = layer_pressure / np.sum(layer_pressure)
    
    # Weighted average of percentage error
    total_error = np.sum(weights * 100*(layer_height - layer_pressure)/layer_pressure)
    
    return total_error

# Example usage
total_error = total_pct_error(errors)
print(f"Total weighted % error: {total_error:f}%")


def cumulative_column_and_error(errors):
    layer_pressure = errors['layer_pressure_DU']
    layer_height = errors['layer_height_DU']
    z_mid = errors['z_mid']

    # Cumulative sum from surface to top
    cum_pressure = np.cumsum(layer_pressure)
    cum_height = np.cumsum(layer_height)

    # Cumulative % error at each layer
    cum_pct_error = 100 * (cum_height - cum_pressure) / cum_pressure

    # Plot
    fig, ax1 = plt.subplots(figsize=(6,8))

    # Left axis: cumulative ozone column
    ax1.plot(cum_pressure, z_mid, label='Pressure-based', color='blue')
    ax1.plot(cum_height, z_mid, label='Height-trapezoid', color='red', linestyle='--')
    ax1.set_xlabel('Cumulative O$_3$ column [DU]')
    ax1.set_ylabel('Altitude [km]')
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)
    ax1.legend(loc='upper left')
    ax1.grid(True, linestyle=':')

    # Right axis: cumulative % error
    ax2 = ax1.twiny()
    ax2.plot(cum_pct_error, z_mid, color='green', label='Cumulative % error')
    ax2.set_xlabel('Cumulative % error')
    ax2.legend(loc='upper right')

    plt.title('Cumulative O$_3$ Column and Trapezoidal Error vs Altitude')
    plt.show()

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

def Photochem_O3_col_mid(df, g=9.81, MO2=0.21):
    # Physical constants
    m_air = ((1 - MO2) * 28 + MO2 * 32) * 1.66054e-27  # kg

    # Sort surface → top (pressure decreasing)
    df = df.sort_values('press', ascending=False)

    # Midpoint values
    chi = df['O3'].to_numpy()           # mixing ratio (dimensionless)
    p_mid = df['press'].to_numpy() * 1e5  # Pa

    # Reconstruct interface pressures (second-order)
    p_half = 0.5 * (p_mid[:-1] + p_mid[1:])

    # Add top and bottom boundaries by extrapolation
    p_top = p_mid[-1] + (p_mid[-1] - p_half[-1])
    p_surf = p_mid[0] + (p_mid[0] - p_half[0])

    p_interfaces = np.concatenate(([p_surf], p_half, [p_top]))

    # Layer pressure thicknesses
    dp = p_interfaces[:-1] - p_interfaces[1:]

    # Midpoint-rule column integral
    N_O3 = np.sum(chi * dp) / (g * m_air)

    # Convert to Dobson Units
    return N_O3 / 2.6867e20

def Photochem_O3_col(df, g=9.81, MO2 = 0.21):
    # Sort top-to-bottom or bottom-to-top safely
    # Physical constants
    g = 9.81             # m s^-2
    m_air = ((1-MO2)*28+(MO2*32)) * 1.66054e-27  # kg (mean molecular mass of air)

    # Sort by pressure (surface → top)
    df = df.sort_values('press', ascending=False)

    chi = df['O3'].to_numpy() # mixing ratio
    p = df['press'].to_numpy()*1e5

    # Pressure layer thicknesses
    dp = -np.diff(p)

    # Column density [molecules m^-2]
    N_O3 = np.sum(chi[:-1] * dp) / (g * m_air)

    # Convert to Dobson Units
    return N_O3 / 2.6867e20

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
def V_O3_col_z_sum(DS):
    species = DS['variable']['species']
    O3_dens = DS['variable']['y'][:,species.index('O3')] * 1e6
    
    z = np.array((DS['atm']['zco']/1e2))
    dz = np.empty(0)
    for i in range(len(DS['atm']['zco']/1e5)-1):
        dz_new = z[i+1] - z[i]
        dz = np.append(dz, dz_new)
        
    O3_col = (O3_dens * dz).sum() / 2.6867e20
    
    return O3_col

def V_O3_col_z_trapz(DS):
    species = DS['variable']['species']
    O3_dens = DS['variable']['y'][:,species.index('O3')] * 1e6
    
    z = np.array((DS['atm']['zmco']/1e2))
        
    O3_col = np.trapz(O3_dens, z) / 2.6867e20
    
    return O3_col

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
#One50_10xCH4_col = V_O3_col_z_trapz(One50_pc_V_10xCH4_482SZA)
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
One_10xCH4_col = V_O3_col_z_trapz(One_pc_V_10xCH4_482SZA)
One_5xCH4_col = V_O3_col_z_trapz(One_pc_V_5xCH4_482SZA)
One_05xCH4_col = V_O3_col_z_trapz(One_pc_V_05xCH4_482SZA)
One_01xCH4_col = V_O3_col_z_trapz(One_pc_V_01xCH4_482SZA)

Z5_col = V_O3_col_z_trapz(Z5_pc_V_482SZA)
Z5_10xCH4_col = V_O3_col_z_trapz(Z5_pc_V_10xCH4_482SZA)
Z5_5xCH4_col = V_O3_col_z_trapz(Z5_pc_V_5xCH4_482SZA)
Z5_05xCH4_col = V_O3_col_z_trapz(Z5_pc_V_05xCH4_482SZA)
Z5_01xCH4_col = V_O3_col_z_trapz(Z5_pc_V_01xCH4_482SZA)

Z1_col = V_O3_col_z_trapz(Z1_pc_V_482SZA)
Z1_10xCH4_col = V_O3_col_z_trapz(Z1_pc_V_10xCH4_482SZA)
Z1_5xCH4_col = V_O3_col_z_trapz(Z1_pc_V_5xCH4_482SZA)
Z1_05xCH4_col = V_O3_col_z_trapz(Z1_pc_V_05xCH4_482SZA)
Z1_01xCH4_col = V_O3_col_z_trapz(Z1_pc_V_01xCH4_482SZA)

# Data mapping to make looping possible
# Format: (Title, Methane_Data_List, Column_Data_List, has_fixed_flux)
plot_configs = [
    ('150% PAL', [0.1, 0.5, 1, 5], [One50_01xCH4_col, One50_05xCH4_col, One50_col, One50_5xCH4_col], True),
    ('PI', [0.1, 0.5, 1, 5, 10], [PI_01xCH4_col, PI_05xCH4_col, PI_col, PI_5xCH4_col, PI_10xCH4_col], True),
    ('50% PAL', [0.1, 0.5, 1, 5, 10], [Fifty_01xCH4_col, Fifty_05xCH4_col, Fifty_col, Fifty_5xCH4_col, Fifty_10xCH4_col], True),
    ('10% PAL', [0.1, 0.5, 1, 5, 10], [Ten_01xCH4_col, Ten_05xCH4_col, Ten_col, Ten_5xCH4_col, Ten_10xCH4_col], True),
    ('5% PAL', [0.1, 0.5, 1, 5, 10], [Five_01xCH4_col, Five_05xCH4_col, Five_col, Five_5xCH4_col, Five_10xCH4_col], True),
    ('1% PAL', [0.1, 0.5, 1, 5, 10], [One_01xCH4_col, One_05xCH4_col, One_col, One_5xCH4_col, One_10xCH4_col], False),
    ('0.5% PAL', [0.1, 0.5, 1, 5, 10], [Z1_01xCH4_col, Z1_05xCH4_col, Z5_col, Z5_5xCH4_col, Z5_10xCH4_col], False),
    ('0.1% PAL', [0.1, 0.5, 1, 5, 10], [Z1_01xCH4_col, Z1_05xCH4_col, Z1_col, Z1_5xCH4_col, Z1_10xCH4_col], False)
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

'''
Plot figure
'''

fig = plt.figure(figsize=(22, 12))
gs = gridspec.GridSpec(2, 4, hspace=0.15, wspace=0.25)

for i, (title, meth_vals, col_vals, has_flux) in enumerate(plot_configs):
    ax = fig.add_subplot(gs[i // 4, i % 4])
    
    # Plot standard Fixed MR
    methane_array = np.array(meth_vals) * 0.8e-6
    ax.plot(methane_array, col_vals, lw=lw, color='magenta', label='Fixed MR, VULCAN')
    
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
            ax.plot(flux_methane, flux_cols, lw=lw, color='m', ls='--', label='Fixed Flux, VULCAN')

        if (i == 1):
            flux_methane = [
                PI_V_01xCH4f_482SZA['variable']['ymix'][:,PI_pc_spec_01xCH4f_482SZA.index('CH4')][0],
                PI_V_05xCH4f_482SZA['variable']['ymix'][:,PI_pc_spec_05xCH4f_482SZA.index('CH4')][0],
                PI_V_1xCH4f_482SZA['variable']['ymix'][:,PI_pc_spec_1xCH4f_482SZA.index('CH4')][0],
                PI_V_5xCH4f_482SZA['variable']['ymix'][:,PI_pc_spec_5xCH4f_482SZA.index('CH4')][0],
                PI_V_10xCH4f_482SZA['variable']['ymix'][:,PI_pc_spec_10xCH4f_482SZA.index('CH4')][0]
            ]
            flux_cols = [PI_01xCH4f_col, PI_05xCH4f_col, PI_1xCH4f_col, PI_5xCH4f_col, PI_10xCH4f_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='magenta', ls='--', label='Fixed Flux, VULCAN')
            
            flux_methane = [
                Photo_PI_01CH4['CH4'][0],Photo_PI_05CH4['CH4'][0],Photo_PI_1CH4['CH4'][0],
                Photo_PI_2CH4['CH4'][0], Photo_PI_3CH4['CH4'][0], Photo_PI_4CH4['CH4'][0],
                Photo_PI_5CH4['CH4'][0]
            ]
            flux_cols = [Photo_PI_01CH4_col, Photo_PI_05CH4_col, Photo_PI_1CH4_col, Photo_PI_2CH4_col,
                         Photo_PI_3CH4_col, Photo_PI_4CH4_col, Photo_PI_5CH4_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='b', ls='-', label='Fixed MR, Photochem')

        if (i == 2):
            flux_methane = [
                Fifty_pc_V_01xCH4f_482SZA['variable']['ymix'][:,Fifty_pc_spec_01xCH4f_482SZA.index('CH4')][0],
                Fifty_pc_V_05xCH4f_482SZA['variable']['ymix'][:,Fifty_pc_spec_05xCH4f_482SZA.index('CH4')][0],
                Fifty_pc_V_1xCH4f_482SZA['variable']['ymix'][:,Fifty_pc_spec_1xCH4f_482SZA.index('CH4')][0],
                Fifty_pc_V_5xCH4f_482SZA['variable']['ymix'][:,Fifty_pc_spec_5xCH4f_482SZA.index('CH4')][0],
                Fifty_pc_V_10xCH4f_482SZA['variable']['ymix'][:,Fifty_pc_spec_10xCH4f_482SZA.index('CH4')][0]
                ]
            flux_cols = [Fifty_01xCH4f_col, Fifty_05xCH4f_col, Fifty_1xCH4f_col, Fifty_5xCH4f_col, Fifty_10xCH4f_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='magenta', ls='--', label='Fixed Flux, VULCAN')

        if (i == 3):
            flux_methane = [
                Ten_pc_V_01xCH4f_482SZA['variable']['ymix'][:,Ten_pc_spec_01xCH4f_482SZA.index('CH4')][0],
                Ten_pc_V_05xCH4f_482SZA['variable']['ymix'][:,Ten_pc_spec_05xCH4f_482SZA.index('CH4')][0],
                Ten_pc_V_1xCH4f_482SZA['variable']['ymix'][:,Ten_pc_spec_1xCH4f_482SZA.index('CH4')][0],
                Ten_pc_V_5xCH4f_482SZA['variable']['ymix'][:,Ten_pc_spec_5xCH4f_482SZA.index('CH4')][0],
                Ten_pc_V_10xCH4f_482SZA['variable']['ymix'][:,Ten_pc_spec_10xCH4f_482SZA.index('CH4')][0]
                ]
            flux_cols = [Ten_01xCH4f_col, Ten_05xCH4f_col, Ten_1xCH4f_col, Ten_5xCH4f_col, Ten_10xCH4f_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='magenta', ls='--', label='Fixed Flux, VULCAN')
        
        if (i == 4):
            flux_methane = [
                Five_pc_V_01xCH4f_482SZA['variable']['ymix'][:,Five_pc_spec_01xCH4f_482SZA.index('CH4')][0],
                Five_pc_V_05xCH4f_482SZA['variable']['ymix'][:,Five_pc_spec_05xCH4f_482SZA.index('CH4')][0],
                Five_pc_V_1xCH4f_482SZA['variable']['ymix'][:,Five_pc_spec_1xCH4f_482SZA.index('CH4')][0],
                Five_pc_V_5xCH4f_482SZA['variable']['ymix'][:,Five_pc_spec_5xCH4f_482SZA.index('CH4')][0],
                Five_pc_V_10xCH4f_482SZA['variable']['ymix'][:,Five_pc_spec_10xCH4f_482SZA.index('CH4')][0]
                ]
            flux_cols = [Five_01xCH4f_col, Five_05xCH4f_col, Five_1xCH4f_col, Five_5xCH4f_col, Five_10xCH4f_col]
            ax.plot(flux_methane, flux_cols, lw=lw, color='magenta', ls='--', label='Fixed Flux, VULCAN')
        

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
        ax.legend(loc='best', fontsize=12)

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
plt.plot(O2, cols, marker = 's', color = 'k')
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
plt.plot(TP1e_MS_60SZA['variable']['ymix'][:,TP1e_MS_60SZA_spec.index('O3')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'k', ls = '--', label = 'O'+sub(3))
plt.plot(TP1e_P19_60SZA['variable']['ymix'][:,TP1e_MS_60SZA_spec.index('O3')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'k', ls = '--', label = 'O'+sub(3))
plt.plot(LWAV(Pre_pc_h0.O3), Pre_pc_h0.lev/1e3, color = 'k', lw = 2)
plt.xscale('log')
plt.xlim(1e-9, 1e-3)
plt.ylim(1e0, 1e-7)
plt.yscale('log')


#%% Ozone column plot

alpha = 0.25
photochem = [P_col_01pc, P_col_1pc, P_col_10pc, P_col_PI, P_col_150pc]
photochem_45 = [61.82, 155.41, 278.6, 365.87, 361.38]
#photochem1 = [P_01pc_col*1.01, P_1pc_col*1.01, P_10pc_col*1.01, P_PI_col*1.01]

o2_conc_less = [0.001, 0.01, 0.1, 1, 1.5]
kasting_45sza = [60.5, 96.4, 244.1, 403.79, 422.4]
kasting_60sza = [32.96, 86.8, 202.3, 308.7, 316.6]
#WACCM_0obq = [LWAV(Zero1_0obq_col), LWAV(One_0obq_col), LWAV(Ten_0obq_col), LWAV(Pre_0obq_col)]
#WACCM_0obq_max = [Zero1_0obq_col.max(), One_0obq_col.max(), Ten_0obq_col.max(), Pre_0obq_col.max()]
#WACCM_0obq_min = [Zero1_0obq_col.min(), One_0obq_col.min(), Ten_0obq_col.min(), Pre_0obq_col.min()]

#Atmos
A_SZA_60 = [45.52, 145.28, 190.81, 207.26, 203.2]
A_SZA_45 = [54.21, 187.43, 251.14, 280.12, 277.57]

o2_conc = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 1.5]
WACCM = [LWAV(Zero1_col), LWAV(Zero5_col), LWAV(One_col), LWAV(Five_col), LWAV(Ten_col), LWAV(Fifty_col), LWAV(Pre_col), LWAV(One50_col)]
WACCM_max = [Zero1_col.max(), Zero5_col.max(), One_col.max(), Five_col.max(), Ten_col.max(), Fifty_col.max(), Pre_col.max(), One50_col.max()]
WACCM_min = [Zero1_col.min(), Zero5_col.min(), One_col.min(), Five_col.min(), Ten_col.min(), Fifty_col.min(), Pre_col.min(), One50_col.min()]
'''
SZA_58 = [E_WPT_V_01pc_60SZA, E_WPT_V_05pc_60SZA, E_WPT_V_1pc_60SZA,
          E_WPT_V_5pc_60SZA, E_WPT_V_10pc_60SZA, E_WPT_V_50pc_60SZA,
          E_WPT_V_60SZA, E_WPT_V_150pc_60SZA]
'''


SZA_58 = [Earth_V_58SZA_01pc_O3_col, Earth_V_58SZA_05pc_O3_col, Earth_V_58SZA_1pc_O3_col,
          Earth_V_58SZA_5pc_O3_col, Earth_V_58SZA_10pc_O3_col,  Earth_V_58SZA_50pc_O3_col,
          Earth_V_58SZA_O3_col, Earth_V_58SZA_150pc_O3_col]
SZA_45 = [Earth_V_45SZA_01pc_O3_col, Earth_V_45SZA_05pc_O3_col, Earth_V_45SZA_1pc_O3_col,
          Earth_V_45SZA_5pc_O3_col, Earth_V_45SZA_10pc_O3_col,  Earth_V_45SZA_50pc_O3_col,
          Earth_V_45SZA_O3_col, Earth_V_45SZA_150pc_O3_col]

SZA_58 = [E_WPT_V_01pc_60SZA_O3_col, E_WPT_V_05pc_60SZA_O3_col, E_WPT_V_1pc_60SZA_O3_col,
          E_WPT_V_5pc_60SZA_O3_col, E_WPT_V_10pc_60SZA_O3_col, E_WPT_V_50pc_60SZA_O3_col,
          E_WPT_V_60SZA_O3_col, E_WPT_V_150pc_60SZA_O3_col]

SZA_45 = [E_WPT_V_01pc_45SZA_O3_col, E_WPT_V_05pc_45SZA_O3_col, E_WPT_V_1pc_45SZA_O3_col,
          E_WPT_V_5pc_45SZA_O3_col, E_WPT_V_10pc_45SZA_O3_col, E_WPT_V_50pc_45SZA_O3_col,
          E_WPT_V_45SZA_O3_col, E_WPT_V_150pc_45SZA_O3_col]

plt.figure(figsize = (12,6))
#plt.fill_between(o2_conc_0obq, WACCM_0obq_min,WACCM_0obq_max, alpha = alpha, color = 'b', label = 'WACCM6; Cooke et al. (2022)')
#plt.plot(o2_conc_0obq, WACCM_0obq, color = 'b')
#plt.fill_between(o2_conc, photochem, photochem1, alpha = alpha, color = 'b', label = 'Photochem, ?'+r'$^\circ$'+' SZA')
#plt.plot(o2_conc_less, photochem, color = 'b', label = 'Photochem, 60'+r'$^\circ$'+' SZA')
plt.fill_between(o2_conc_less, photochem, photochem_45, alpha = alpha, color = 'b', label = 'Photochem, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, SZA_45, SZA_58, alpha = alpha, color = 'm', label = 'VULCAN, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, WACCM_min, WACCM_max, alpha = alpha, color = 'k', label = 'WACCM6; Cooke et al. (2022)', lw = 4)
plt.fill_between(o2_conc_less, A_SZA_60, A_SZA_45, alpha = alpha, color = 'darkorange', label = 'Atmos 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.plot(o2_conc, WACCM, color = 'k')
plt.xscale('log')
thick_axes(top = True)
#plt.title('Simulations with Cooke et al. (2022) WACCM6 boundary conditions', fontsize = 15, weight = 'bold')
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
plt.title('')
plt.xlim(1e-3, 1.5)
markersize = 7
plt.ylim(0,380)
plt.plot([0.001, 1], [18, 330], marker = 's', ls = '', color = 'teal', markersize = markersize, label = 'Kasting 1D model; Ji et al. (2024)')
plt.plot([1],[292], color = 'k', marker = 'o', markersize = markersize)
plt.legend(loc = 0, fontsize = 15, frameon = False)
plt.xlabel('Oxygen concentration [PAL]', fontsize = 15, weight = 'bold')
plt.ylabel('O'+sub(3)+' column [DU]', fontsize = 15, weight = 'bold')
#plt.savefig('/Users/gregcooke/python_output/Ozone_column_vs_o2_curve'+save+'.png', dpi = 200, bbox_inches = 'tight')
plt.savefig('/Users/gregcooke/python_output/Ozone_column_vs_o2_curve.png', dpi = 200, bbox_inches = 'tight')
#%%

# 1. Create a common grid for the 1D models to find the min/max envelope
o2_common = [1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 0.5, 1, 1.5] # Grid from 1e-3 to ~1.5 PAL

# Interpolate all 1D models onto the common grid
photo_interp = np.interp(o2_common, o2_conc_less, photochem)
photo45_interp = np.interp(o2_common, o2_conc_less, photochem_45)
vulcan45_interp = np.interp(o2_common, o2_conc, SZA_45)
vulcan58_interp = np.interp(o2_common, o2_conc, SZA_58)
atmos45_interp = np.interp(o2_common, o2_conc_less, A_SZA_45)
atmos60_interp = np.interp(o2_common, o2_conc_less, A_SZA_60)
kasting45_interp = np.interp(o2_common, o2_conc_less, kasting_45sza)
kasting60_interp = np.interp(o2_common, o2_conc_less, kasting_60sza)

# 2. Calculate the collective 1D envelope
all_1d_stack = np.vstack([photo_interp, photo45_interp, 
                          vulcan45_interp, vulcan58_interp, 
                          atmos45_interp, atmos60_interp, 
                          kasting45_interp, kasting60_interp])

one_d_min = np.min(all_1d_stack, axis=0)
one_d_max = np.max(all_1d_stack, axis=0)

# --- Updated Plotting Function ---
def plot_panel(ax, highlight=None, show_legend=False, show_ylabel=True):
    # WACCM Background (Always present)
    ax.fill_between(o2_conc, WACCM_min, WACCM_max, alpha=0.3, color='k', label='WACCM6 Range')
    ax.plot(o2_conc, WACCM, color='k', lw=2, label='WACCM6 Mean')

    # Top Panel: Unified 1D Range
    if highlight == 'unified_1d':
        ax.fill_between(o2_common, one_d_min, one_d_max, alpha=0.4, color='tab:blue', 
                        label='Total 1D Model Range', lw=0)
        
    # Bottom Panels: Specific Model Highlights
    elif highlight == 'kasting':
        ax.fill_between(o2_conc_less, kasting_60sza, kasting_45sza, color='teal',  alpha=0.3, label='Kasting 1D range')
    elif highlight == 'photochem':
        ax.fill_between(o2_conc_less, photochem, photochem_45, color='b', alpha=0.3, label='Photochem range')
    elif highlight == 'atmos':
        ax.fill_between(o2_conc_less, A_SZA_60, A_SZA_45, color='darkorange', alpha=0.3, label='Atmos range')
    elif highlight == 'vulcan':
        ax.fill_between(o2_conc, SZA_45, SZA_58, color='m', alpha=0.3, label='VULCAN range')

    # Formatting
    ax.set_xscale('log')
    ax.set_xlim(1e-3, 1.5)
    ax.set_ylim(0, 380)
    if show_ylabel: ax.set_ylabel('O$_3$ column [DU]', weight='bold', fontsize=15)
    if show_legend: ax.legend(loc='upper left', frameon=False, fontsize=10)
    thick_axes(top=True, labelleft = True)

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
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def plot_panel(ax, highlight=None, show_legend=False, show_ylabel=True):
    """
    highlight options: 'photochem', 'vulcan', 'atmos', 'kasting', 'all'
    """
    # 1. Always plot the WACCM background (The 3D benchmark)
    ax.fill_between(o2_conc, WACCM_min, WACCM_max, alpha=0.2, color='k', 
                    label='WACCM6 (Cooke 2022)', zorder=1)
    ax.plot(o2_conc, WACCM, color='k', lw=2, zorder=2)

    # 2. Plotting Logic
    # We define the 1D model plotting as separate blocks so we can call them individually or all at once
    
    def add_photochem():
        ax.fill_between(o2_conc_less, photochem, photochem_45, alpha=alpha, color='b', 
                        label='Photochem', lw=2)
    def add_vulcan():
        ax.fill_between(o2_conc, SZA_45, SZA_58, alpha=alpha, color='m', 
                        label='VULCAN', lw=2)
    def add_atmos():
        ax.fill_between(o2_conc_less, A_SZA_60, A_SZA_45, alpha=alpha, color='darkorange', 
                        label='Atmos', lw=2)
    def add_kasting():
        ax.plot([0.001, 1], [18, 330], marker='s', ls='', color='teal', 
                markersize=markersize, label='Kasting 1D (Ji 2024)')

    if highlight == 'all':
        add_photochem()
        add_vulcan()
        add_atmos()
        add_kasting()
    elif highlight == 'photochem': add_photochem()
    elif highlight == 'vulcan':    add_vulcan()
    elif highlight == 'atmos':     add_atmos()
    elif highlight == 'kasting':   add_kasting()
    
    # Global Styling
    ax.set_xscale('log')
    ax.set_xlim(1e-3, 1.5)
    ax.set_ylim(0, 380)
    if show_ylabel:
        ax.set_ylabel('O$_3$ column [DU]', fontsize=12, weight='bold')
    if show_legend:
        ax.legend(loc='upper left', fontsize=9, frameon=False, ncol=2)
    
    thick_axes(top=True)

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
plt.figure(figsize = (14,18))
gs = gridspec.GridSpec(3,4)

plt.subplot(gs[0,:])

plt.fill_between(o2_conc_less, photochem, photochem_45, alpha = alpha, color = 'b', label = 'Photochem, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, SZA_45, SZA_58, alpha = alpha, color = 'm', label = 'VULCAN, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, WACCM_min, WACCM_max, alpha = alpha, color = 'k', label = 'WACCM6; Cooke et al. (2022)', lw = 4)
plt.fill_between(o2_conc_less, A_SZA_60, A_SZA_45, alpha = alpha, color = 'darkorange', label = 'Atmos 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.plot(o2_conc, WACCM, color = 'k')
plt.xscale('log')
thick_axes(top = True)
#plt.title('Simulations with Cooke et al. (2022) WACCM6 boundary conditions', fontsize = 15, weight = 'bold')
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
plt.title('')
markersize = 7
plt.ylim(0,380)
plt.plot([0.001, 1], [18, 330], marker = 's', ls = '', color = 'teal', markersize = markersize, label = 'Kasting 1D model; Ji et al. (2024)')
plt.plot([1],[292], color = 'k', marker = 'o', markersize = markersize)
plt.legend(loc = 0, fontsize = 15, frameon = False)
plt.xlabel('Oxygen concentration [PAL]', fontsize = 15, weight = 'bold')
plt.ylabel('O'+sub(3)+' column [DU]', fontsize = 15, weight = 'bold')
plt.xlim(1e-3,1.5)

plt.subplot(gs[1,0:2])

plt.fill_between(o2_conc_less, photochem, photochem_45, alpha = alpha, color = 'b', label = 'Photochem, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, SZA_45, SZA_58, alpha = alpha, color = 'm', label = 'VULCAN, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, WACCM_min, WACCM_max, alpha = alpha, color = 'k', label = 'WACCM6; Cooke et al. (2022)', lw = 4)
plt.fill_between(o2_conc_less, A_SZA_60, A_SZA_45, alpha = alpha, color = 'darkorange', label = 'Atmos 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.plot(o2_conc, WACCM, color = 'k')
plt.xscale('log')
thick_axes(top = True)
#plt.title('Simulations with Cooke et al. (2022) WACCM6 boundary conditions', fontsize = 15, weight = 'bold')
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
plt.title('')
markersize = 7
plt.ylim(0,380)
plt.plot([0.001, 1], [18, 330], marker = 's', ls = '', color = 'teal', markersize = markersize, label = 'Kasting 1D model; Ji et al. (2024)')
plt.plot([1],[292], color = 'k', marker = 'o', markersize = markersize)
plt.xlabel('Oxygen concentration [PAL]', fontsize = 15, weight = 'bold')
plt.ylabel('O'+sub(3)+' column [DU]', fontsize = 15, weight = 'bold')
plt.xlim(1e-3,1.5)

plt.subplot(gs[1,2:4])

plt.fill_between(o2_conc_less, photochem, photochem_45, alpha = alpha, color = 'b', label = 'Photochem, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, SZA_45, SZA_58, alpha = alpha, color = 'm', label = 'VULCAN, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, WACCM_min, WACCM_max, alpha = alpha, color = 'k', label = 'WACCM6; Cooke et al. (2022)', lw = 4)
plt.fill_between(o2_conc_less, A_SZA_60, A_SZA_45, alpha = alpha, color = 'darkorange', label = 'Atmos 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.plot(o2_conc, WACCM, color = 'k')
plt.xscale('log')
thick_axes(top = True, labelleft = False)
#plt.title('Simulations with Cooke et al. (2022) WACCM6 boundary conditions', fontsize = 15, weight = 'bold')
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
plt.title('')
markersize = 7
plt.ylim(0,380)
plt.plot([0.001, 1], [18, 330], marker = 's', ls = '', color = 'teal', markersize = markersize, label = 'Kasting 1D model; Ji et al. (2024)')
plt.plot([1],[292], color = 'k', marker = 'o', markersize = markersize)
plt.xlabel('Oxygen concentration [PAL]', fontsize = 15, weight = 'bold')
plt.xlim(1e-3,1.5)

plt.subplot(gs[2,0:2])

plt.fill_between(o2_conc_less, photochem, photochem_45, alpha = alpha, color = 'b', label = 'Photochem, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, SZA_45, SZA_58, alpha = alpha, color = 'm', label = 'VULCAN, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, WACCM_min, WACCM_max, alpha = alpha, color = 'k', label = 'WACCM6; Cooke et al. (2022)', lw = 4)
plt.fill_between(o2_conc_less, A_SZA_60, A_SZA_45, alpha = alpha, color = 'darkorange', label = 'Atmos 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.plot(o2_conc, WACCM, color = 'k')
plt.xscale('log')
thick_axes(top = True)
#plt.title('Simulations with Cooke et al. (2022) WACCM6 boundary conditions', fontsize = 15, weight = 'bold')
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
plt.title('')
markersize = 7
plt.ylim(0,380)
plt.plot([0.001, 1], [18, 330], marker = 's', ls = '', color = 'teal', markersize = markersize, label = 'Kasting 1D model; Ji et al. (2024)')
plt.plot([1],[292], color = 'k', marker = 'o', markersize = markersize)
plt.xlabel('Oxygen concentration [PAL]', fontsize = 15, weight = 'bold')
plt.ylabel('O'+sub(3)+' column [DU]', fontsize = 15, weight = 'bold')
plt.xlim(1e-3,1.5)

plt.subplot(gs[2,2:4])

plt.fill_between(o2_conc_less, photochem, photochem_45, alpha = alpha, color = 'b', label = 'Photochem, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, SZA_45, SZA_58, alpha = alpha, color = 'm', label = 'VULCAN, 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.fill_between(o2_conc, WACCM_min, WACCM_max, alpha = alpha, color = 'k', label = 'WACCM6; Cooke et al. (2022)', lw = 4)
plt.fill_between(o2_conc_less, A_SZA_60, A_SZA_45, alpha = alpha, color = 'darkorange', label = 'Atmos 45'+r'$^\circ$'+'- 60'+r'$^\circ$'+' SZA', lw = 4)
plt.plot(o2_conc, WACCM, color = 'k')
plt.xscale('log')
thick_axes(top = True, labelleft = False)
#plt.title('Simulations with Cooke et al. (2022) WACCM6 boundary conditions', fontsize = 15, weight = 'bold')
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
plt.title('')
markersize = 7
plt.ylim(0,380)
plt.plot([0.001, 1], [18, 330], marker = 's', ls = '', color = 'teal', markersize = markersize, label = 'Kasting 1D model; Ji et al. (2024)')
plt.plot([1],[292], color = 'k', marker = 'o', markersize = markersize)
plt.xlabel('Oxygen concentration [PAL]', fontsize = 15, weight = 'bold')
plt.xlim(1e-3,1.5)

plt.savefig('/Users/gregcooke/python_output/Ozone_column_vs_o2_curve.png', dpi = 200, bbox_inches = 'tight')

#%%
plt.figure(figsize = (11,6))
plt.plot(o2_conc, np.array(SZA_58) / np.array(WACCM), color = 'm', ls = '-', lw = 2, label = 'VULCAN 60'+r'$^\circ$'+' SZA')
plt.plot(o2_conc, np.array(SZA_45) / np.array(WACCM), color = 'm', ls = '--', lw = 2,  label = 'VULCAN 45'+r'$^\circ$'+' SZA')
WACCM = [LWAV(Zero1_col), LWAV(One_col), LWAV(Ten_col),  LWAV(Pre_col), LWAV(One50_col)]
plt.plot(o2_conc_less, np.array(A_SZA_60) / np.array(WACCM), color = 'darkorange', ls = '-', lw = 2,  label = 'Atmos 60'+r'$^\circ$'+' SZA')
plt.plot(o2_conc_less, np.array(A_SZA_45) / np.array(WACCM), color = 'darkorange', ls = '--', lw = 2,  label = 'Atmos 45'+r'$^\circ$'+' SZA')
plt.plot(o2_conc_less, np.array(photochem) / np.array(WACCM), color = 'b', ls = '-', lw = 2,  label = 'Photochem 60'+r'$^\circ$'+' SZA')
plt.plot(o2_conc_less, np.array(photochem_45) / np.array(WACCM), color = 'b', ls = '--', lw = 2,  label = 'Photochem 45'+r'$^\circ$'+' SZA')
plt.legend(loc = 0, fontsize = 15, frameon = False)
plt.plot([1],[292/LWAV(Pre_col)], color = 'k', marker = 'o', markersize = markersize)
plt.xscale('log')
thick_axes(top = True)
plt.xlabel('Oxygen concentration [PAL]', fontsize = 15, weight = 'bold')
plt.axhline(1, color = 'k', lw = 2, ls = ':')
plt.ylabel('O'+sub(3)+' column ratio to WACCM6 mean', fontsize = 15, weight = 'bold')

#%%
xlim = 1e-5

plt.figure(figsize = (10,6))
gs = gridspec.GridSpec(2, 2)
plt.subplot(gs[0,0])
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('O3')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'k', ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V_45SZA['variable']['ymix'][:,Earth_V_spec.index('O3')], Earth_V_45SZA['atm']['pco']/1e6, lw = lw, color = 'k', ls = ':', label = 'O'+sub(3))
plt.plot(LWAV(Pre_pc_h0.O3), Pre_pc_h0.lev/1e3, color = 'k', lw = 2)
plt.xscale('log')
plt.xlim(1e-9, xlim)
plt.ylim(1e0, 1e-7)
plt.yscale('log')

plt.subplot(gs[0,1])
plt.plot(Earth_V_10pc['variable']['ymix'][:,Earth_V_spec_10pc.index('O3')], Earth_V_10pc['atm']['pco']/1e6, lw = lw, color = Ten_pc_color, ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V_10pc_45SZA['variable']['ymix'][:,Earth_V_spec_10pc_45SZA.index('O3')], Earth_V_10pc_45SZA['atm']['pco']/1e6, lw = lw, color = Ten_pc_color, ls = ':', label = 'O'+sub(3))
plt.plot(LWAV(Ten_pc_h0.O3), Ten_pc_h0.lev/1e3, color = Ten_pc_color, lw = 2)
plt.xscale('log')
plt.xlim(1e-9, xlim)
plt.ylim(1e0, 1e-7)
plt.yscale('log')

plt.subplot(gs[1,0])
plt.plot(Earth_V_1pc['variable']['ymix'][:,Earth_V_spec_1pc.index('O3')], Earth_V_1pc['atm']['pco']/1e6, lw = lw, color = One_pc_color, ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V_1pc_45SZA['variable']['ymix'][:,Earth_V_spec_1pc_45SZA.index('O3')], Earth_V_1pc_45SZA['atm']['pco']/1e6, lw = lw, color = One_pc_color, ls = ':', label = 'O'+sub(3))
plt.plot(LWAV(One_pc_h0.O3), One_pc_h0.lev/1e3, color = One_pc_color, lw = 2)
plt.xscale('log')
plt.xlim(1e-9, xlim)
plt.ylim(1e0, 1e-7)
plt.yscale('log')

plt.subplot(gs[1,1])
plt.plot(Earth_V_01pc['variable']['ymix'][:,Earth_V_spec_01pc.index('O3')], Earth_V_01pc['atm']['pco']/1e6, lw = lw, color = Zero1_pc_color, ls = '--', label = 'O'+sub(3))
#plt.plot(Earth_V_01pc_45SZA['variable']['ymix'][:,Earth_V_spec_01pc_45SZA.index('O3')], Earth_V_01pc_45SZA['atm']['pco']/1e6, lw = lw, color = Zero1_pc_color, ls = ':', label = 'O'+sub(3))
plt.plot(LWAV(Zero1_pc_h0.O3), Zero1_pc_h0.lev/1e3, color = Zero1_pc_color, lw = 2)
plt.xscale('log')
plt.xlim(1e-9, xlim)
plt.ylim(1e0, 1e-7)
plt.yscale('log')

xlim = 1e-4
xlim1 = 1e-8

plt.figure(figsize = (10,6))
gs = gridspec.GridSpec(2, 2)
plt.subplot(gs[0,0])
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('H2O')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'k', ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V_45SZA['variable']['ymix'][:,Earth_V_spec.index('H2O')], Earth_V_45SZA['atm']['pco']/1e6, lw = lw, color = 'k', ls = ':', label = 'O'+sub(3))
plt.plot(LWAV(Pre_pc_h0.H2O), Pre_pc_h0.lev/1e3, color = 'k', lw = 2)
plt.xscale('log')
plt.xlim(xlim1, xlim)
plt.ylim(1e0, 1e-7)
plt.yscale('log')

plt.subplot(gs[0,1])
plt.plot(Earth_V_10pc['variable']['ymix'][:,Earth_V_spec_10pc.index('H2O')], Earth_V_10pc['atm']['pco']/1e6, lw = lw, color = Ten_pc_color, ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V_10pc_45SZA['variable']['ymix'][:,Earth_V_spec_10pc_45SZA.index('H2O')], Earth_V_10pc_45SZA['atm']['pco']/1e6, lw = lw, color = Ten_pc_color, ls = ':', label = 'O'+sub(3))
plt.plot(LWAV(Ten_pc_h0.H2O), Ten_pc_h0.lev/1e3, color = Ten_pc_color, lw = 2)
plt.xscale('log')
plt.xlim(xlim1, xlim)
plt.ylim(1e0, 1e-7)
plt.yscale('log')

plt.subplot(gs[1,0])
plt.plot(Earth_V_1pc['variable']['ymix'][:,Earth_V_spec_1pc.index('H2O')], Earth_V_1pc['atm']['pco']/1e6, lw = lw, color = One_pc_color, ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V_1pc_45SZA['variable']['ymix'][:,Earth_V_spec_1pc_45SZA.index('H2O')], Earth_V_1pc_45SZA['atm']['pco']/1e6, lw = lw, color = One_pc_color, ls = ':', label = 'O'+sub(3))
plt.plot(LWAV(One_pc_h0.H2O), One_pc_h0.lev/1e3, color = One_pc_color, lw = 2)
plt.xscale('log')
plt.xlim(xlim1, xlim)
plt.ylim(1e0, 1e-7)
plt.yscale('log')

plt.subplot(gs[1,1])
plt.plot(Earth_V_01pc['variable']['ymix'][:,Earth_V_spec_01pc.index('H2O')], Earth_V_01pc['atm']['pco']/1e6, lw = lw, color = Zero1_pc_color, ls = '--', label = 'O'+sub(3))
#plt.plot(Earth_V_01pc_45SZA['variable']['ymix'][:,Earth_V_spec_01pc_45SZA.index('H2O')], Earth_V_01pc_45SZA['atm']['pco']/1e6, lw = lw, color = Zero1_pc_color, ls = ':', label = 'O'+sub(3))
plt.plot(LWAV(Zero1_pc_h0.H2O), Zero1_pc_h0.lev/1e3, color = Zero1_pc_color, lw = 2)
plt.xscale('log')
plt.xlim(xlim1, xlim)
plt.ylim(1e0, 1e-7)
plt.yscale('log')


xlim = 1e-6
xlim1 = 1e-9
mol = 'CH4'

plt.figure(figsize = (10,6))
gs = gridspec.GridSpec(2, 2)
plt.subplot(gs[0,0])
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index(mol)], Earth_V['atm']['pco']/1e6, lw = lw, color = 'k', ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V_45SZA['variable']['ymix'][:,Earth_V_spec.index(mol)], Earth_V_45SZA['atm']['pco']/1e6, lw = lw, color = 'k', ls = ':', label = 'O'+sub(3))
plt.plot(LWAV(Pre_pc_h0[mol]), Pre_pc_h0.lev/1e3, color = 'k', lw = 2)
plt.xscale('log')
plt.xlim(xlim1, xlim)
plt.ylim(1e0, 1e-7)
plt.yscale('log')

plt.subplot(gs[0,1])
plt.plot(Earth_V_10pc['variable']['ymix'][:,Earth_V_spec_10pc.index(mol)], Earth_V_10pc['atm']['pco']/1e6, lw = lw, color = Ten_pc_color, ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V_10pc_45SZA['variable']['ymix'][:,Earth_V_spec_10pc_45SZA.index(mol)], Earth_V_10pc_45SZA['atm']['pco']/1e6, lw = lw, color = Ten_pc_color, ls = ':', label = 'O'+sub(3))
plt.plot(LWAV(Ten_pc_h0[mol]), Ten_pc_h0.lev/1e3, color = Ten_pc_color, lw = 2)
plt.xscale('log')
plt.xlim(xlim1, xlim)
plt.ylim(1e0, 1e-7)
plt.yscale('log')

plt.subplot(gs[1,0])
plt.plot(Earth_V_1pc['variable']['ymix'][:,Earth_V_spec_1pc.index(mol)], Earth_V_1pc['atm']['pco']/1e6, lw = lw, color = One_pc_color, ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V_1pc_45SZA['variable']['ymix'][:,Earth_V_spec_1pc_45SZA.index(mol)], Earth_V_1pc_45SZA['atm']['pco']/1e6, lw = lw, color = One_pc_color, ls = ':', label = 'O'+sub(3))
plt.plot(LWAV(One_pc_h0[mol]), One_pc_h0.lev/1e3, color = One_pc_color, lw = 2)
plt.xscale('log')
plt.xlim(xlim1, xlim)
plt.ylim(1e0, 1e-7)
plt.yscale('log')

plt.subplot(gs[1,1])
plt.plot(Earth_V_01pc['variable']['ymix'][:,Earth_V_spec_01pc.index(mol)], Earth_V_01pc['atm']['pco']/1e6, lw = lw, color = Zero1_pc_color, ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V_01pc_45SZA['variable']['ymix'][:,Earth_V_spec_01pc_45SZA.index(mol)], Earth_V_01pc_45SZA['atm']['pco']/1e6, lw = lw, color = Zero1_pc_color, ls = ':', label = 'O'+sub(3))
plt.plot(LWAV(Zero1_pc_h0[mol]), Zero1_pc_h0.lev/1e3, color = Zero1_pc_color, lw = 2)
plt.xscale('log')
plt.xlim(xlim1, xlim)
plt.ylim(1e0, 1e-7)
plt.yscale('log')

#%%


plt.figure()
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('O3')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'k', ls = '-', label = 'O'+sub(3))
plt.plot(np.flip(LWAV(Pre_pc_h0.O3)), np.flip(Pre_pc_h0.lev/1e3), lw = lw, color = 'k', ls = '--', label = 'O'+sub(3))
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('H2O')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'blue', ls = '-', label = 'O'+sub(3))
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('CO2')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'green', ls = '-', label = 'O'+sub(3))
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('O2')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'grey', ls = '-', label = 'O'+sub(3))
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('N2O')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'purple', ls = '-', label = 'O'+sub(3))
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('CH4')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'orange', ls = '-', label = 'O'+sub(3))
plt.xscale('log'); plt.xlim(1e-9,1)
plt.text(1e-5, 1e-4, str(int(Earth_V_58SZA_O3_col)) + ' DU')
plt.yscale('log'); plt.ylim(1,1e-7)

plt.figure()
plt.plot(Earth_V['variable']['ymix'][:,Earth_V_spec.index('O3')], Earth_V['atm']['pco']/1e6, lw = lw, color = 'k', ls = '-', label = 'O'+sub(3))
plt.plot(Earth_V_45SZA['variable']['ymix'][:,Earth_V_spec.index('O3')], Earth_V_45SZA['atm']['pco']/1e6, lw = lw, color = 'k', ls = ':', label = 'O'+sub(3))
plt.plot(np.flip(LWAV(Pre_pc_h0.O3)), np.flip(Pre_pc_h0.lev/1e3), lw = lw, color = 'k', ls = '--', label = 'O'+sub(3))
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
One_pc_V_482SZA_O3_col = V_O3_col_z_trapz(One_pc_V_482SZA)
Five_pc_V_482SZA_O3_col = V_O3_col_z_trapz(Five_pc_V_482SZA)
Z5_pc_V_482SZA_O3_col = V_O3_col_z_trapz(Z5_pc_V_482SZA)
Z1_pc_V_482SZA_O3_col = V_O3_col_z_trapz(Z1_pc_V_482SZA)

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
plt.plot([1.5], One50_pc_V_482SZA_O3_col, marker = 'o', color = 'k', markersize = markersize, lw = 0, label = 'VULCAN SZA = 48.2'+r'$^\circ$')
plt.plot([1], PI_V_482SZA_O3_col, marker = 's', color = 'k', markersize = markersize, lw = 2)
plt.plot([0.5], Fifty_pc_V_482SZA_O3_col, marker = 'o', color = 'k', markersize = markersize, lw = 2)
plt.plot([0.1], Ten_pc_V_482SZA_O3_col, marker = 'o', color = 'k', markersize = markersize, lw = 2)
plt.plot([0.05], Five_pc_V_482SZA_O3_col, marker = 'o', color = 'k', markersize = markersize, lw = 2)
plt.plot([0.01], One_pc_V_482SZA_O3_col, marker = 'o', color = 'k', markersize = markersize, lw = 2)
plt.plot([0.005], Z5_pc_V_482SZA_O3_col, marker = 'o', color = 'k', markersize = markersize, lw = 2)
plt.plot([0.001], Z1_pc_V_482SZA_O3_col, marker = 'o', color = 'k', markersize = markersize, lw = 2)

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

#%% Create new VULCAN files

DS = PCb
# --------------------------------------------------
# Input / output files
# --------------------------------------------------
infile = "/Users/gregcooke/VIH_cases/atm/atm_Earth_Jan_Kzz.txt"
outfile = "/Users/gregcooke/VIH_cases/atm/atm_Earth_PCb_0.1pc_o2_WACCM_Kzz.txt"
#outfile = "/Users/gregcooke/Downloads/atm_Earth_PCb_PI_WACCM_Kzz.txt"
#outfile = "/Users/gregcooke/VIH_cases/atm/atm_Earth_0.1pc_o2_WACCM_Kzz.txt"

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
    "# (bar) (K)     (cm2/s)\n"
    "# Pressure        Temp     Kzz\n"
)

with open(outfile, "w") as f:
    f.write(header)
    for p, t, k in zip(P_new, T_new, Kzz_new):
        f.write(f"{p:12.5E} {t:8.2f} {k:10.3E}\n")

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
plt.plot(PCb_V_45SZA['variable']['ymix'][:,PCb_spec_V_45SZA.index('O3')], PCb_1pc_V_45SZA['atm']['pco']/1e6, lw = lw, color = 'k', ls = ls, label = 'O'+sub(3))
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
plt.errorbar([0.001, 0.01, 0.1, 1], W_cols, yerr=[W_cols-W_min_cols,W_max_cols-W_cols], capsize = capsize,  lw = lw, marker = 's', markersize = markersize, color = 'k')
V_cols = [PCb_01pc_V_45SZA_O3_col, PCb_1pc_V_45SZA_O3_col, PCb_10pc_V_45SZA_O3_col, PCb_V_45SZA_O3_col]
plt.plot([0.001, 0.01, 0.1, 1], V_cols, color = 'magenta', lw = lw, marker = 'o', markersize = markersize)
V_cols = [PCb_01pc_V_60SZA_O3_col, PCb_1pc_V_60SZA_O3_col, PCb_10pc_V_60SZA_O3_col, PCb_V_60SZA_O3_col]
plt.plot([0.001, 0.01, 0.1, 1], V_cols, color = 'pink', lw = lw, marker = 'o', markersize = markersize)
plt.ylim(5,2000)
plt.yscale('log')
plt.xscale('log')
thick_axes(labelleft=True, top = True, right = True)

#%% 
import pandas as pd
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
plt.plot(GJ551['Wav'], flux, color = 'k')

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

import pandas as pd
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
plt.plot(GJ551['Wav'], flux, color = 'k')
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
plt.plot(WACCM_wavelength, WACCM_Flux, color = 'k')
plt.plot(WACCM_wavelength, P19_Flux, color = 'darkorange')
solar_file = 'TRAPPIST1_flux_at_e_W21.nc'#read in file
ds = xr.open_dataset(solar_dir+solar_file)#attach file to dataset
ssi = ds['ssi'].isel(time=0) #define dataset from file
W21_Flux = ssi.values
WACCM_wavelength = ssi.wavelength.values
plt.plot(WACCM_wavelength, W21_Flux, color = 'navy')

flux = GJ551['Flux']*((r_star*r_sun)/(au*orbit_radius))**2
plt.plot(GJ551['Wav'], flux, color = 'k')
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
plt.plot(GJ551['Wav'], flux, color = 'k')
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

#%% Include Her

plt.figure(figsize=(14, 7))
gs = gridspec.GridSpec(1, 2)

# Subplot 1: GCSE
plt.subplot(gs[0, 0])
plt.title('Ages 14 - 16 (GCSE / NQ5)', fontsize=15, weight = 'bold')
bars1 = plt.bar(['Women', 'Men'], [1, 76], color=['#1c579e', '#1c9e79'], alpha = 0.9)
plt.ylabel('Number of scientists', fontsize=15, weight = 'bold')
plt.xticks(fontsize=15)  # Set x-axis fontsize
plt.ylim(0, 80)          # Increased slightly to fit labels
plt.grid(False)


# Add labels to GCSE bars
counter = 0
for bar in bars1:
    counter = counter +1
    yval = bar.get_height()
    if (counter == 1):
        plt.text(bar.get_x() + bar.get_width()/2, 3 , yval, ha='center', va='bottom', fontsize=15)
    else:
        plt.text(bar.get_x() + bar.get_width()/2, yval/2, yval, ha='center', va='bottom', fontsize=15)
thick_axes(top = True)
# Subplot 2: A-Level
plt.subplot(gs[0, 1])
plt.title('Ages 16 - 18 (A-Level / Scottish Highers)', fontsize=15, weight = 'bold')
bars2 = plt.bar(['Women', 'Men'], [3, 162], color=['#1c579e', '#1c9e79'], alpha = 0.9)
plt.xticks(fontsize=15)  # Set x-axis fontsize
plt.ylim(0, 170)         # Increased slightly to fit labels
plt.grid(False)
# Add labels to A-Level bars
counter = 0
for bar in bars2:
    counter = counter +1
    yval = bar.get_height()
    if (counter == 1):
        plt.text(bar.get_x() + bar.get_width()/2, 3 , yval, ha='center', va='bottom', fontsize=15)
    else:
        plt.text(bar.get_x() + bar.get_width()/2, yval/2, yval, ha='center', va='bottom', fontsize=15)
thick_axes(top = True)
# Save the complete figure
plt.savefig('/Users/gregcooke/python_output/A-Level_Number_of_Scientists.png', bbox_inches = 'tight')

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

def cumulative_col(ds, time = False, var = 'O2', O2_mr = 0.21, lon = True, g = 9.80616):
    
    g = g # gravity acceleration (cm/sec2)
    to_DU = 1./2.6867e20  # convert colunm density / cm^2 to Dobson units
    N2_mr = 1-O2_mr #nitrogen mixing ratio
    m_avg = ((28*N2_mr)+(O2_mr*32))*1.661e-27 # average weight of atmosphere
    
    #calcualte column integrals 
    var = ds[var]
    
    # variables to calcualate hybrid pressure on model interfaces
    ps = ds['PS'].values # surface pressure
    p0 = ds['P0'].values
    hyai = ds['hyai'].values
    hybi = ds['hybi'].values
    

    #col = np.ndarray(ps.shape, np.float32)
    dp = np.ndarray(var.shape, np.float32)
    
    nt = 1
    try:
        nt = ds['time'].size
    except:
        nt = 0
    ny = ds['lat'].size
    nz = ds['lev'].size
    nx = ds['lon'].size
    
    press = np.ndarray(nz, np.float32)
    if (nt > 0):
        for i in range(nt):
            for j in range(ny):
                for k in range(nx):
                    press = hyai * p0 + hybi * ps[i,j,k]
                    for iz in range(nz):
                        dp[i,iz,j,k] = (press[iz+1]-press[iz])
    else:
        for j in range(ny):
            for k in range(nx):
                press = hyai * p0 + hybi * ps[j,k]
                for iz in range(nz):
                    dp[iz,j,k] = (press[iz+1]-press[iz])
                    
    #calcualte column integrals 
    var_col = var * dp / (m_avg * g)
    
    #o3col = np.sum(o3, axis = 0)  # convert to DU
    if (lon == False):
        var_col  = var_col.mean(dim='lon')
    if (time == False):
        try:
            var_col  = var_col.mean(dim='time')
        except:
            var_col  = var_col
    print('Columns calculated')
    return var_col

def int_O2_photo(DS, O2 = 0.21, g = 9.1454):
    O2_num_dens = cumulative_col(DS, time = False, var = 'O2', O2_mr = O2, lon = True, g = g)
    int_JO2 = O2_num_dens * 2*(DS.jo2_a+DS.jo2_b)
    int_JO2 = int_JO2.sum(dim='lev').mean(dim = 'lon')
    return int_JO2

def cum_int_O2_photo(DS, O2 = 0.21, g = 9.1454):
    O2_num_dens = cumulative_col(DS, time = False, var = 'O2', O2_mr = O2, lon = True, g = g)
    int_JO2 = O2_num_dens * 2*(DS.jo2_a+DS.jo2_b)
    int_JO2 = LWAV(int_JO2.cumsum(dim='lev'))
    return int_JO2

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

ylim = (1e3, 1e-3)
xlim = (1e11, 6e12)

P_path = '/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/'

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

plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)
gs.update(wspace = 0.1, hspace = 0.1)
plt.subplot(gs[0,0])
plt.title('100% PAL  ', fontsize = 15, weight = 'bold', y = 0.9, loc = 'right')
plt.plot(LWAV(prox_ox_W(base)), base.lev, color = 'k', label = 'WACCM6', lw = 2)
plt.plot(prox_ox_K(K_100pc_J), p_100pc, color = 'teal', label = 'Kasting\n1D model', lw = 2)
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
plt.plot(LWAV( prox_ox_W(Ten_pc_new)), base.lev, color = 'k', lw = 2)
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
plt.plot(LWAV( prox_ox_W(One_pc_new)), base.lev, color = 'k', lw = 2)
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
plt.plot(LWAV( prox_ox_W(Zero1_pc_new)), base.lev, color = 'k', lw = 2)
plt.plot(prox_ox_K(K_01pc_J), p_01pc, color = 'teal', lw = 2)
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
plt.plot(LWAV(n_dens(Pre_pc_h0, "O3")), Pre_pc_h0.lev, color = 'k', label = 'WACCM6')
plt.plot(PI_V_482SZA['variable']['y'][:,PI_spec_482SZA.index('O3')]*1e6, PI_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN')
plt.plot(Photo_PI["O3"]*P_dens(Photo_PI), Photo_PI["press"]*1e3, color = 'b', label = 'Photochem', lw = 2)
plt.plot(Atmos_100pc["O3"]*Atmos_dens(Atmos_100pc), Atmos_100pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = 2)
plt.plot(K_100pc_J["FO3"]*K_100pc_J["DEN"]*1e6, K_100pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting\n1D model', lw = 2)
plt.yscale('log')
plt.xlim(xlim)
plt.ylim(ylim)
plt.xscale('log');  thick_axes(top = True)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[0,1])
plt.title('10% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(n_dens(Ten_pc_h0, "O3")), Pre_pc_h0.lev, color = 'k', lw = lw)
plt.plot(Ten_pc_V_482SZA['variable']['y'][:,Ten_pc_spec_482SZA.index('O3')]*1e6, Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(Photo_10pc["O3"]*P_dens(Photo_10pc), Photo_10pc["press"]*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_10pc["O3"]*Atmos_dens(Atmos_10pc), Atmos_10pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_10pc_J["FO3"]*K_10pc_J["DEN"]*1e6, K_10pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting\n1D model', lw = 2)
plt.yscale('log')
plt.xlim(xlim)
plt.ylim(ylim)
plt.xscale('log');  thick_axes(top = True, labelleft = False)

plt.subplot(gs[1,0])
plt.title('1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(n_dens(One_pc_h0, "O3")), Pre_pc_h0.lev, color = 'k', lw = lw, label = 'WACCM6')
plt.plot(One_pc_V_482SZA['variable']['y'][:,One_pc_spec_482SZA.index('O3')]*1e6, Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN', lw = lw)
plt.plot(Photo_1pc["O3"]*P_dens(Photo_1pc), Photo_1pc["press"]*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_1pc["O3"]*Atmos_dens(Atmos_1pc), Atmos_1pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_1pc_J["FO3"]*K_1pc_J["DEN"]*1e6, K_1pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting\n1D model', lw = 2)
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
plt.plot(LWAV(n_dens(Zero1_pc_h0, "O3")), Pre_pc_h0.lev, color = 'k', lw = lw)
plt.plot(Z1_pc_V_482SZA['variable']['y'][:,Z1_pc_spec_482SZA.index('O3')]*1e6, Z1_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(Photo_01pc["O3"]*P_dens(Photo_01pc), Photo_01pc["press"]*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_01pc["O3"]*Atmos_dens(Atmos_01pc), Atmos_01pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_01pc_J["FO3"]*K_01pc_J["DEN"]*1e6, K_01pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting\n1D model', lw = 2)
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
plt.plot(LWAV(Pre_pc_h2.NOX), Pre_pc_h0.lev, color = 'k', label = 'WACCM6', lw = lw)
N = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('N')]
NO = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('NO')]
NO2 = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('NO2')]
NOX = N+NO+NO2
plt.plot(NOX, PI_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN', lw = lw)
plt.plot(P_NOX(Photo_PI), Photo_PI['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_NOX(Atmos_100pc), Atmos_100pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(Kasting_NOX(K_100pc_J), K_100pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting 1D model', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-5)
plt.xscale('log');  thick_axes(top = True)
plt.xlim(1e-12, 1e-6)
plt.legend(loc = (0.03, 0.4), handlelength = 1,
           fontsize = 15, frameon = False)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[0,1])
plt.title('10% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Ten_pc_h2.NOX), Pre_pc_h0.lev, color = 'k', lw = lw)
N = Ten_pc_V_482SZA['variable']['ymix'][:,Ten_pc_spec_482SZA.index('N')]
NO = Ten_pc_V_482SZA['variable']['ymix'][:,Ten_pc_spec_482SZA.index('NO')]
NO2 = Ten_pc_V_482SZA['variable']['ymix'][:,Ten_pc_spec_482SZA.index('NO2')]
NOX = N+NO+NO2
plt.plot(Atmos_NOX(Atmos_10pc), Atmos_10pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(NOX, Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(P_NOX(Photo_10pc), Photo_10pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Kasting_NOX(K_10pc_J), K_10pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting 1D model', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-5)
plt.xscale('log');  thick_axes(top = True, labelleft = False)
plt.xlim(1e-12, 1e-6)

plt.subplot(gs[1,0])
plt.title('1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(One_pc_h2.NOX), Pre_pc_h0.lev, color = 'k', lw = lw)
N = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('N')]
NO = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('NO')]
NO2 = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('NO2')]
NOX = N+NO+NO2
plt.plot(NOX, One_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(Atmos_NOX(Atmos_1pc), Atmos_1pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(P_NOX(Photo_1pc), Photo_1pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Kasting_NOX(K_1pc_J), K_1pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting 1D model', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-5)
plt.xscale('log');  thick_axes(top = True)
plt.xlim(1e-12, 1e-6)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.xlabel('NO'+sub('\\rm x')+' mixing ratio', fontsize = 15, weight = 'bold')

plt.subplot(gs[1,1])
plt.title('0.1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Zero1_pc_h2.NOX), Pre_pc_h0.lev, color = 'k', lw = lw)
N = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('N')]
NO = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('NO')]
NO2 = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('NO2')]
NOX = N+NO+NO2
plt.plot(NOX, Z1_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(P_NOX(Photo_01pc), Photo_01pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_NOX(Atmos_01pc), Atmos_01pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(Kasting_NOX(K_01pc_J), K_01pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting 1D model', lw = lw)
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
plt.plot(LWAV(Pre_pc_h2.HOX), Pre_pc_h0.lev, color = 'k', label = 'WACCM6', lw = lw)
H = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('H')]
OH = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('OH')]
HO2 = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('HO2')]
H2O2 = PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('H2O2')]
HOX = H+OH+HO2+H2O2
plt.plot(HOX, PI_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN', lw = lw)
plt.plot(P_HOX(Photo_PI), Photo_PI['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_HOX(Atmos_100pc), Atmos_100pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(Kasting_HOX(K_100pc_J), K_100pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting\n1D model', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-3)
plt.xscale('log');  thick_axes(top = True)
plt.xlim(xlim)
plt.legend(loc = (0.03, 0.4), handlelength = 1,
           fontsize = 15, frameon = False)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[0,1])
plt.title('10% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Ten_pc_h2.HOX), Pre_pc_h0.lev, color = 'k', lw = lw)
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
plt.plot(LWAV(One_pc_h2.HOX), Pre_pc_h0.lev, color = 'k', lw = lw)
H = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('H')]
OH = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('OH')]
HO2 = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('HO2')]
H2O2 = One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('H2O2')]
HOX = H+OH+HO2+H2O2
plt.plot(HOX, One_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(P_HOX(Photo_1pc), Photo_1pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_HOX(Atmos_1pc), Atmos_1pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(Kasting_HOX(K_1pc_J), K_1pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting 1D model', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-3)
plt.xscale('log');  thick_axes(top = True)
plt.xlim(xlim)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.xlabel('HO'+sub('\\rm x')+' mixing ratio', fontsize = 15, weight = 'bold')

plt.subplot(gs[1,1])
plt.title('0.1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Zero1_pc_h2.HOX), Pre_pc_h0.lev, color = 'k', lw = lw)
H = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('H')]
OH = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('OH')]
HO2 = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('HO2')]
H2O2 = Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('H2O2')]
HOX = H+OH+HO2+H2O2
plt.plot(HOX, Z1_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(P_HOX(Photo_01pc), Photo_01pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_HOX(Atmos_01pc), Atmos_01pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(Kasting_HOX(K_01pc_J), K_01pc_J["PRESS"]/1e3, color = color_Kasting, label = 'Kasting 1D model', lw = lw)
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
plt.plot(LWAV(Pre_pc_h0.H2O), Pre_pc_h0.lev, color = 'k', label = 'WACCM6', lw = lw)
plt.plot(PI_V_482SZA['variable']['ymix'][:,PI_spec_482SZA.index('H2O')], PI_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN')
plt.plot(Photo_PI["H2O"], Photo_PI['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_100pc["H2O"], Atmos_100pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_100pc_J["FH2O"], K_100pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting\n1D model', lw = lw)
plt.yscale('log')
plt.ylim(1e3, ylim)
plt.xscale('log');  thick_axes(top = True)
plt.xlim(1e-7, 1e-4)
plt.legend(loc = (0.01, 0.01), handlelength = 1,
           fontsize = 15, frameon = False)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[0,1])
plt.title('10% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Ten_pc_h0.H2O), Ten_pc_h0.lev, color = 'k', lw = lw)
plt.plot(Ten_pc_V_482SZA['variable']['ymix'][:,Ten_pc_spec_482SZA.index('H2O')], Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN', lw = lw)
plt.plot(Photo_10pc["H2O"], Photo_10pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_10pc["H2O"], Atmos_10pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_10pc_J["FH2O"], K_10pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting\n1D model', lw = lw)
plt.yscale('log')
plt.ylim(1e3, ylim)
plt.xscale('log');  thick_axes(top = True, labelleft = False)
plt.xlim(1e-7, 1e-4)

plt.subplot(gs[1,0])
plt.title('1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(One_pc_h0.H2O), One_pc_h0.lev, color = 'k', lw = lw)
plt.plot(One_pc_V_482SZA['variable']['ymix'][:,One_pc_spec_482SZA.index('H2O')], One_pc_V_482SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN', lw = lw)
plt.plot(Photo_1pc["H2O"], Photo_1pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_1pc["H2O"], Atmos_1pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_1pc_J["FH2O"], K_1pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting\n1D model', lw = lw)
plt.yscale('log')
plt.ylim(1e3, 1e-3)
plt.xscale('log');  thick_axes(top = True)
plt.xlim(1e-7, ylim)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.xlabel('H'+sub(2)+'O mixing ratio', fontsize = 15, weight = 'bold')

plt.subplot(gs[1,1])
plt.title('0.1% PAL', fontsize = 15, weight = 'bold', y = 0.9)
plt.plot(LWAV(Zero1_pc_h0.H2O), Zero1_pc_h0.lev, color = 'k', lw = lw)
plt.plot(Z1_pc_V_482SZA['variable']['ymix'][:,Z1_pc_spec_482SZA.index('H2O')], Z1_pc_V_60SZA['atm']['pco']/1e3, color = 'm', label = 'VULCAN', lw = lw)
plt.plot(Photo_01pc["H2O"], Photo_01pc['press']*1e3, color = 'b', label = 'Photochem', lw = lw)
plt.plot(Atmos_01pc["H2O"], Atmos_01pc["PRESS"]*1e3, color = 'darkorange', label = 'Atmos', lw = lw)
plt.plot(K_01pc_J["FH2O"], K_01pc_J["PRESS"]/1e3, color = 'teal', label = 'Kasting\n1D model', lw = lw)
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

def photochem_integrated_JO2(O2_photo):
    
    # local production per layer (cm^-2 s^-1)
    layer_prod = O2_photo * 1000
    
    # cumulative integral from top of atmosphere downward
    int_JO2 = layer_prod[::-1].cumsum()[::-1]
    
    return int_JO2


def vulcan_integrated_JO2(data, Vulcan_factor=Vulcan_factor):
    
    # species list
    species = data['variable']['species']
    O2_ind = species.index('O2')
    
    # O2 number density (cm^-3)
    O2_dens = data['variable']['y'][:, O2_ind] * 1e6
    
    # total O2 photolysis rate
    JO2 = (
        data['variable']['J_sp']['O2',0] +
        data['variable']['J_sp']['O2',1] +
        data['variable']['J_sp']['O2',2]
    )
    
    # layer thickness (cm)
    dz = data['atm']['dz']/1000
    
    # local production per layer (cm^-2 s^-1)
    layer_prod = Vulcan_factor * JO2 * O2_dens * dz
    
    # cumulative integral from top of atmosphere downward
    int_JO2 = layer_prod[::-1].cumsum()[::-1]
    
    return int_JO2


# --- Compute VULCAN integrated JO2 ---
int_JO2_PI = vulcan_integrated_JO2(PI_pc_V_482SZA)
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

# --- Plot ---
plt.figure(figsize=(14,10))
gs = gridspec.GridSpec(2,2)

xlim = (5e15, 2e17)

plt.subplot(gs[0,0])
# WACCM6
plt.plot(cum_int_base, base.lev, color='k', label='WACCM6', lw = lw)
# VULCAN case
plt.plot(int_JO2_PI, PI_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw, label='VULCAN')
plt.plot(int_ox_prod_100pc_photo, Photo_PI['press']*1e3, color = 'b', lw = lw, label='Photochem')
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.legend(loc = 0, fontsize = 15, frameon = False)
thick_axes(top = True)
plt.yscale('log'); plt.ylim(1e3, 1e-5)
plt.xscale('log'); plt.xlim(xlim)

plt.subplot(gs[0,1])
plt.plot(int_JO2_Ten, Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(int_ox_prod_10pc_photo, Photo_10pc['press']*1e3, color = 'b', lw = lw, label='Photochem')
plt.plot(cum_int_ten, base.lev, color = 'k', lw = lw)
thick_axes(top = True)
plt.yscale('log'); plt.ylim(1e3, 1e-5)
plt.xscale('log');  plt.xlim(xlim)

plt.subplot(gs[1,0])
plt.plot(int_JO2_One, One_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(cum_int_one, base.lev, color = 'k', lw = lw)
plt.plot(int_ox_prod_1pc_photo, Photo_1pc['press']*1e3, color = 'b', lw = lw, label='Photochem')
thick_axes(top = True)
plt.yscale('log'); plt.ylim(1e3, 1e-5)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')
plt.xscale('log');  plt.xlim(xlim)
plt.xlabel('Integrated O$_2$ photolysis', fontsize = 15, weight = 'bold')

plt.subplot(gs[1,1])
plt.plot(int_JO2_Zero1, Z1_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = lw)
plt.plot(int_ox_prod_01pc_photo, Photo_01pc['press']*1e3, color = 'b', lw = lw, label='Photochem')
plt.plot(cum_int_zero1, base.lev, color = 'k', lw = lw)

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

plt.plot(PI_photo, base.lev, lw = lw, color = 'k')

lw = 2
V_PI_photo = (
    PI_pc_V_60SZA['variable']['J_sp']['O2',0] +
    PI_pc_V_60SZA['variable']['J_sp']['O2',1] +
    PI_pc_V_60SZA['variable']['J_sp']['O2',2]
)
species = PI_pc_V_60SZA['variable']['species']
O2_dens_PI = PI_pc_V_60SZA['variable']['y'][:,species.index('O2')] * 1e6
V_PI_prod_ox = Vulcan_factor*V_PI_photo*O2_dens_PI
plt.plot(V_PI_prod_ox, PI_pc_V_60SZA['atm']['pco']/1e3, lw = lw, ls = ':', color = 'k')
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

from scipy.interpolate import interp1d

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
plt.plot(LWAV(jo2), base.lev, lw = lw, color = 'k')
thick_axes()
plt.ylim(1e3,1e-5)
plt.xscale('log')
plt.yscale('log')
#plt.xlim(1e10,1e13)

plt.figure(figsize = (7,5))

plt.plot(PI_PCb_photo , base.lev, lw = lw, color = 'k')
#plt.plot(Ten_photo, base.lev, lw = lw, color = Ten_pc_color)
#plt.plot(One_photo, base.lev, lw = lw, color = One_pc_color)
plt.plot(Zero1_photo, base.lev, lw = lw, color = Zero1_pc_color)

plt.title('PCb simulations', fontsize = 15, weight = 'bold')
lw = 2
V_PCb_PI_photo = PCb_V_60SZA['variable']['J_sp']['O2', 0] + PCb_V_60SZA['variable']['J_sp']['O2', 1] + PCb_V_60SZA['variable']['J_sp']['O2', 2]
species = PCb_V_60SZA['variable']['species']
O2_dens = PCb_V_60SZA['variable']['y'][:,species.index('O2')] * 1e5
plt.plot(Vulcan_factor*V_PCb_PI_photo*O2_dens, PCb_V_60SZA['atm']['pco']/1e3, lw = lw, ls = ':', color = 'k')

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

plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)

plt.subplot(gs[0, 0])
PCb_W_JO2 = prox_ox_W(PCb)
PCb_V_JO2 = prox_ox_V(PCb_V_60SZA, PCb_spec_V_60SZA)
plt.plot(LWAV(PCb_W_JO2), PCb.lev, lw = lw, color = 'k', label = 'WACCM6')
plt.plot(PCb_V_JO2, PCb_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)


plt.subplot(gs[0, 1])
PCb_W_JO2 = prox_ox_W(PCb_10pc)
PCb_V_JO2 = prox_ox_V(PCb_10pc_V_60SZA, PCb_10pc_spec_V_60SZA)
plt.plot(LWAV(PCb_W_JO2), PCb.lev, lw = lw, color = 'k', label = 'WACCM6')
plt.plot(PCb_V_JO2, PCb_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)


plt.subplot(gs[1, 0])
PCb_W_JO2 = prox_ox_W(PCb_1pc)
PCb_V_JO2 = prox_ox_V(PCb_1pc_V_60SZA, PCb_1pc_spec_V_60SZA)
plt.plot(LWAV(PCb_W_JO2), PCb.lev, lw = lw, color = 'k', label = 'WACCM6')
plt.plot(PCb_V_JO2, PCb_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)

plt.subplot(gs[1, 1])
PCb_W_JO2 = prox_ox_W(PCb_01pc)
PCb_V_JO2 = prox_ox_V(PCb_01pc_V_60SZA, PCb_01pc_spec_V_60SZA)
plt.plot(LWAV(PCb_W_JO2), PCb.lev, lw = lw, color = 'k', label = 'WACCM6')
plt.plot(PCb_V_JO2, PCb_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)


'''
Now plot integrated photolysis rate
'''

int_JO2_PI = vulcan_integrated_JO2(PCb_V_60SZA)
int_JO2_Ten = vulcan_integrated_JO2(PCb_10pc_V_60SZA)
int_JO2_One = vulcan_integrated_JO2(PCb_1pc_V_60SZA)
int_JO2_Zero1 = vulcan_integrated_JO2(PCb_01pc_V_60SZA)
# --- Existing cumulative calculations ---
cum_int_base = cum_int_O2_photo(PCb, O2=0.21, g=9.81)
cum_int_ten = cum_int_O2_photo(PCb_10pc, O2=0.021, g=9.81)
cum_int_one = cum_int_O2_photo(PCb_1pc, O2=0.0021, g=9.81)
cum_int_zero1 = cum_int_O2_photo(PCb_01pc, O2=0.00021, g=9.81)

xlim = 1e14, 1e17
plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)

plt.subplot(gs[0, 0])
plt.plot(int_JO2_PI, PCb_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.plot(cum_int_base, PCb.lev, color = 'k', lw =lw, label= 'WACCM6')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')


plt.subplot(gs[0, 1])
plt.plot(int_JO2_Ten, PCb_10pc_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.plot(cum_int_ten, PCb.lev, color = 'k', lw =lw, label= 'WACCM6')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)


plt.subplot(gs[1, 0])
plt.plot(int_JO2_One, PCb_1pc_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.plot(cum_int_one, PCb.lev, color = 'k', lw =lw, label= 'WACCM6')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)
plt.xlabel('Integrated O2 photolysis', fontsize = 15, weight = 'bold')
plt.ylabel('Pressure [hPa]', fontsize = 15, weight = 'bold')

plt.subplot(gs[1, 1])
plt.plot(int_JO2_Zero1, PCb_01pc_V_60SZA['atm']['pco']/1e3, lw = lw, color = 'm', label = 'VULCAN')
plt.plot(cum_int_zero1, PCb.lev, color = 'k', lw =lw, label= 'WACCM6')
plt.ylim(ylim); plt.yscale('log'); plt.xlim(xlim)
plt.xscale('log'); thick_axes(top = True)
plt.xlabel('Integrated O2 photolysis', fontsize = 15, weight = 'bold')

#%% J rates calculation check with aflux


def calculate_J_rates(data):
    sigma1 = data['variable']['cross_J']['O2', 1] # Shape: (2610,)
    sigma2 = data['variable']['cross_J']['O2', 2] # Shape: (2610,)
    flux = data['variable']['aflux']             # Shape: (120, 2610)
    
    # Multiply sigma across the second dimension of flux
    # Resulting j_raw will be (120, 2610)
    j_raw = flux * (sigma1+sigma2)
    
    # Sum over the wavelength axis (axis 1) 
    # This leaves you with one value for each of the 120 layers
    J_rates = np.sum(j_raw, axis=1)
    
    return J_rates

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
plt.plot(WACCM_wavelength, WACCM_Flux/1000, lw = 2, color = 'k', label = 'Sun at Earth', zorder = 3, alpha = alpha)
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
plt.plot(PI_pc_V_482SZA['variable']['bins'], XSfactor * (PI_pc_V_482SZA['variable']['cross_J']['O2',1]+PI_pc_V_482SZA['variable']['cross_J']['O2',2]), color = 'k')
plt.plot(PI_pc_V_482SZA['variable']['bins'], PI_pc_V_482SZA['variable']['cross_J']['O3',1]+PI_pc_V_482SZA['variable']['cross_J']['O3',2], color = 'm')
plt.ylim(1e-21, 1e-16)
#plt.plot((df_flux['Wave_Min']+df_flux['Wave_Max'])[:35]/20, df_cross['O2'], color = 'teal')
plt.xlim(150, 250)
plt.yscale('log')


plt.figure()
plt.plot(PCb_V_60SZA['variable']['y'][:,PCb_spec_V_60SZA.index('O3')] * 1e6, PCb_V_60SZA['atm']['pco']/1e3, color = 'k')
plt.plot(LWAV(n_dens(PCb, "O3")), PCb.lev, color = 'k', ls = '--')
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
plt.plot(PCb_V_60SZA['variable']['y'][:,PCb_spec_V_60SZA.index('O3')] * 1e6, PCb_V_60SZA['atm']['pco']/1e3, color = 'k')
O3 = LWAV_spec(n_dens(PCb, "O3"), time = False, latitude = np.arange(47-5, 47+5))
plt.plot(O3, PCb.lev, color = 'k', ls = '--')
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
plt.plot(XSfactor*PCb_V_60SZA['variable']['y'][:,PCb_spec_V_60SZA.index('O3')] * 1e6, PCb_V_60SZA['atm']['pco']/1e3, color = 'k')
plt.plot(PCb_V_60SZA['variable']['y'][:,PCb_spec_V_60SZA.index('O2')] * 1e6, PCb_V_60SZA['atm']['pco']/1e3, color = 'k', ls = '--')
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
plt.plot(XSfactor*PI_pc_V_482SZA['variable']['y'][:,PI_pc_spec_482SZA.index('O3')] * 1e6, PI_pc_V_482SZA['atm']['pco']/1e3, color = 'k')
plt.plot(PI_pc_V_482SZA['variable']['y'][:,PI_pc_spec_482SZA.index('O2')] * 1e6, PI_pc_V_482SZA['atm']['pco']/1e3, color = 'k', ls = '--')
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
plt.plot(PCb_V_60SZA['variable']['bins'], PCb_V_60SZA['variable']['sflux'][level], color = 'k')
plt.plot(PCb_10pc_V_60SZA['variable']['bins'], PCb_10pc_V_60SZA['variable']['sflux'][level], color = Ten_pc_color)
plt.plot(PCb_1pc_V_60SZA['variable']['bins'], PCb_1pc_V_60SZA['variable']['sflux'][level], color = One_pc_color)
plt.plot(PCb_01pc_V_60SZA['variable']['bins'], PCb_01pc_V_60SZA['variable']['sflux'][level], color = Zero1_pc_color)
plt.xlim(150, 300)
plt.yscale('log')
plt.ylim(1e-9, 1e3)
plt.title(str(PCb_01pc_V_60SZA['atm']['pco'][level]/1e3) + ' hPa')

plt.figure(figsize = (10,5))
plt.plot(PI_pc_V_482SZA['variable']['bins'], PI_pc_V_482SZA['variable']['sflux'][level], color = 'k')
plt.plot(Ten_pc_V_482SZA['variable']['bins'], Ten_pc_V_482SZA['variable']['sflux'][level], color = Ten_pc_color)
plt.plot(One_pc_V_482SZA['variable']['bins'], One_pc_V_482SZA['variable']['sflux'][level], color = One_pc_color)
plt.plot(Z1_pc_V_482SZA['variable']['bins'], Z1_pc_V_482SZA['variable']['sflux'][level], color = Zero1_pc_color)
plt.xlim(150, 300)
plt.yscale('log')
plt.ylim(1e-9, 1e3)
plt.title(str(PI_pc_V_482SZA['atm']['pco'][level]/1e3) + ' hPa')


#%%

import pandas as pd

import pandas as pd

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



#%% Ozone loss

xlim = (1e-2, 1e9)
ylim = (1e3, 1e-3)

plt.figure(figsize = (14,10))
gs = gridspec.GridSpec(2,2)
gs.update(wspace = 0.2, hspace = 0.2)
plt.subplot(gs[0,0])
reaction = PI_V_482SZA['variable']['k'][707]* PI_V_482SZA['variable']['k'][707]*1e6
O =  PI_V_482SZA['variable']['y'][:,PI_spec_482SZA.index('O')] * 1e6
O3 = PI_V_482SZA['variable']['y'][:,PI_spec_482SZA.index('O3')] * 1e6
plt.plot((reaction * O * O3), PI_V_482SZA['atm']['pco']/1e3, color = 'm', lw = 2, label = 'VULCAN')
o3_loss, pressure_hpa = compute_o3loss(
    mech_file=P_path+"zahnle_earth.yaml",
    settings_file=P_path+"input/settings_100pc.yaml",
    flux_file=P_path+"input/Sun_0.0Ga.txt",
    pt_file=P_path+"100pc_PT_profile/100%PAL_Sun_0.0Ga.txt",
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
plt.plot((reaction * O * O3), Ten_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = 2)
o3_loss, pressure_hpa = compute_o3loss(
    mech_file=P_path+"zahnle_earth.yaml",
    settings_file=P_path+"input/settings_100pc.yaml",
    flux_file=P_path+"input/Sun_0.0Ga.txt",
    pt_file=P_path+"10pc_PT_profile/10%PAL_Sun_0.0Ga.txt",
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
plt.plot((reaction * O * O3), PI_V_482SZA['atm']['pco']/1e3, color = 'm', lw = 2)
o3_loss, pressure_hpa = compute_o3loss(
    mech_file=P_path+"zahnle_earth.yaml",
    settings_file=P_path+"input/settings_100pc.yaml",
    flux_file=P_path+"input/Sun_0.0Ga.txt",
    pt_file=P_path+"1pc_PT_profile/1%PAL_Sun_0.0Ga.txt",
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
plt.plot((reaction * O * O3), Z1_pc_V_482SZA['atm']['pco']/1e3, color = 'm', lw = 2)
o3_loss, pressure_hpa = compute_o3loss(
    mech_file=P_path+"zahnle_earth.yaml",
    settings_file=P_path+"input/settings_100pc.yaml",
    flux_file=P_path+"input/Sun_0.0Ga.txt",
    pt_file=P_path+"0.1pc_PT_profile/0.1%PAL_Sun_0.0Ga.txt",
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
PCb_V, PCb_V_spec = Read_O3_Run(file_path='PCb_1e12s_482SZA_WBC_WPT_1rtol.vul')
PCb_50pc_V, PCb_50pc_V_spec = Read_O3_Run(file_path='PCb_50pc_o2_1e12s_482SZA_WBC_WPT_1rtol.vul')
PCb_10pc_V, PCb_10pc_V_spec = Read_O3_Run(file_path='PCb_10pc_o2_1e12s_482SZA_WBC_WPT_1rtol.vul')
PCb_5pc_V, PCb_5pc_V_spec = Read_O3_Run(file_path='PCb_5pc_o2_1e12s_482SZA_WBC_WPT_1rtol.vul')
PCb_1pc_V, PCb_1pc_V_spec = Read_O3_Run(file_path='PCb_1pc_o2_1e12s_482SZA_WBC_WPT_1rtol.vul')

PCb_V_col = V_O3_col_z_trapz(PCb_V)
PCb_10pc_V_col = V_O3_col_z_trapz(PCb_10pc_V)
PCb_1pc_V_col = V_O3_col_z_trapz(PCb_1pc_V)

plt.figure()

O3 =  PCb_V['variable']['ymix'][:,PCb_V_spec.index('O3')]
plt.plot(O3, PCb_V['atm']['pco']/1e3, color = 'm', lw = 2)
O3 =  PCb_10pc_V['variable']['ymix'][:,PCb_10pc_V_spec.index('O3')]
plt.plot(O3, PCb_1pc_V['atm']['pco']/1e3, color = 'm', ls = ':', lw = 2)
O3 =  PCb_1pc_V['variable']['ymix'][:,PCb_1pc_V_spec.index('O3')]
plt.plot(O3, PCb_1pc_V['atm']['pco']/1e3, color = 'm', ls = '--', lw = 2)

plt.yscale('log'); plt.xscale('log')
plt.ylim(1e3, 1e-5); plt.xlim(1e-8, 1e-5)