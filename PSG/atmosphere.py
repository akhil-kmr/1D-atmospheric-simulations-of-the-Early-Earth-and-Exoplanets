#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
process_atmosphere.py

Convert raw atmosphere data (VULCAN, WACCM, Photochem) into PSG-ready dictionaries.
"""

import numpy as np
import pickle

# ----------------------------------------------------------
# Utility functions
# ----------------------------------------------------------

def read_pickle_file(file_path):
    try:
        with open(file_path, 'rb') as file:
            return pickle.load(file)
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
        return None
    except Exception as e:
        print(f"Error reading pickle file: {e}")
        return None

def save_array(data, species_name='H'):
    species = data['variable']['species']
    return data['variable']['ymix'][:, species.index(species_name)]

# ----------------------------------------------------------
# Main processing function
# ----------------------------------------------------------

def process_atmosphere(file_path, 
                       clouds=False, 
                       sulfur=False, 
                       model='VULCAN'):
    """
    Convert raw atmosphere data into a PSG-ready atmosphere dictionary.

    Parameters
    ----------
    file_path : str
        Path to the file.
    clouds : bool
        Include cloud profiles (optional).
    sulfur : bool
        Include sulfur species (SO2, CS2, H2S).
    model : str
        Atmosphere model ('VULCAN' and 'Photochem' currently supported).

    Returns
    -------
    atmosphere : dict
        PSG-ready dictionary with molecules, layers, and header info.
    """

    if model != 'VULCAN' and model != 'Photochem' and model != 'Atmos':
        raise NotImplementedError("Only Photochem and VULCAN models are currently supported")
     
    if (model == 'Atmos'):
        import pandas as pd
        
        data = pd.read_csv(file_path, 
             delim_whitespace=True,   # split on arbitrary whitespace
             engine="python"          # more robust for irregular spacing
        )

        #P-T profile
        temp = data['TEMP']
        press = data['PRESS']
        # --- Base species ---
        H2 = data['H2']
        H2O = data['H2O']
        CH4 = data['CH4']
        CO2 = data['CO2']
        NH3 = data['NH3']
        O2 = data['O2']
        O3 = data['O3']
        N2O = data['N2O']
        #He = data['HE']
        N2 = 1-H2-H2O-O2-CO2-CH4-O3
        
    if (model == 'Photochem'):
        import pandas as pd
        data = pd.read_csv(file_path, delim_whitespace=True)
        #P-T profile
        temp = data['temp']
        press = data['press']
        # --- Base species ---
        H2 = data['H2']
        H2O = data['H2O']
        CH4 = data['CH4']
        CO2 = data['CO2']
        N2 = data['N2']
        NH3 = data['NH3']
        O2 = data['O2']
        O3 = data['O3']
        N2O = data['N2O']
        He = data['He']
    
    if (model == 'VULCAN'):
        # --- Load data ---
        data = read_pickle_file(file_path)
        if data is None:
            return None
        #P-T profile
        temp = data['atm']['Tco']
        press = data['atm']['pco'] / 1e6  # convert Pa -> bar
    
        # --- Base species ---
        H2 = save_array(data, 'H2')
        H2O = save_array(data, 'H2O')
        CH4 = save_array(data, 'CH4')
        CO2 = save_array(data, 'CO2')
        N2 = save_array(data, 'N2')
        NH3 = save_array(data, 'NH3')
        O2 = save_array(data, 'O2')
        O3 = save_array(data, 'O3')
        N2O = save_array(data, 'N2O')
        He = save_array(data, 'He')
        

    molecules = {
        'P': press,
        'T': temp,
        'H2': H2,
        'H2O': H2O,
        'CH4': CH4,
        'CO2': CO2,
        'N2': N2,
        'NH3': NH3,
        'O2': O2,
        'O3': O3,
        'N2O': N2O,
    }

    # --- Optional sulfur species ---
    if sulfur:
        SO2 = save_array(data, 'SO2')
        CS2 = save_array(data, 'CS2')
        H2S = save_array(data, 'H2S')
        molecules.update({'SO2': SO2, 'CS2': CS2, 'H2S': H2S})

    # --- Define molecule strings for PSG headers ---
    gas_keys = [k for k in molecules.keys() if k not in ['P', 'T']]
    GASES = ",".join(gas_keys)
    NGAS = len(gas_keys)

    # HIT, ABUN, ATM_UNIT strings (must match GASES order)
    HIT_map = {
        'H2': 45, 'H2O': 2, 'CH4': 6, 'CO2': 2, 'N2': 22,
        'NH3': 11, 'O2': 7, 'O3': 3, 'N2O': 4,
        'SO2': 9, 'CS2': 53, 'H2S': 22
    }
    HIT = ",".join([f"HIT[{HIT_map[g]}]" for g in gas_keys])
    ABUN = ",".join(["1"]*NGAS)
    ATM_UNIT = ",".join(["scl"]*NGAS)

    # --- Compute mean molecular weight ---
    mmw_map = {
        'He': 4, 'H2': 2, 'H2O': 18, 'CH4': 16, 'CO2': 44,
        'N2': 28, 'NH3': 17, 'O2': 32, 'O3': 48, 'N2O': 44,
        'SO2': 64, 'CS2': 76, 'H2S': 34
    }
    mmw = sum([molecules[g][0]*mmw_map[g] for g in molecules if g in mmw_map])
    print(f"Mean molecular weight = {mmw:.3f}")

    layers = len(press)

    # --- Build PSG-ready dictionary ---
    atmosphere_dict = {
        'molecules': molecules,
        'GASES': GASES,
        'HIT': HIT,
        'ABUN': ABUN,
        'ATM_UNIT': ATM_UNIT,
        'layers': layers,
        'NGAS': NGAS,
        'MMW': mmw
    }

    return atmosphere_dict

# ----------------------------------------------------------
# PSG writer
# ----------------------------------------------------------

def write_psg_atmosphere(atmosphere, output_file):
    """
    Write a PSG-ready atmosphere file from the dictionary.
    """

    mol = atmosphere['molecules']
    gas_keys = [k for k in mol.keys() if k not in ['P', 'T']]
    layer_header = ['P', 'T'] + gas_keys
    layers = atmosphere['layers']

    header = [
        "<ATMOSPHERE-DESCRIPTION>Custom Atmosphere",
        "<ATMOSPHERE-STRUCTURE>Equilibrium",
        f"<ATMOSPHERE-LAYERS>{layers}",
        "<ATMOSPHERE-PUNIT>bar",
        f"<ATMOSPHERE-WEIGHT>{atmosphere['MMW']}",
        "<ATMOSPHERE-CONTINUUM>Rayleigh,Refraction,CIA_all,UV_all",
        f"<ATMOSPHERE-NGAS>{atmosphere['NGAS']}",
        f"<ATMOSPHERE-GAS>{atmosphere['GASES']}",
        f"<ATMOSPHERE-TYPE>{atmosphere['HIT']}",
        f"<ATMOSPHERE-ABUN>{atmosphere['ABUN']}",
        f"<ATMOSPHERE-UNIT>{atmosphere['ATM_UNIT']}",
        "<ATMOSPHERE-LAYERS-MOLECULES>" + ",".join(layer_header)
    ]

    with open(output_file, "w") as fw:
        for line in header:
            fw.write(line + "\n")

        # Write layers
        for i in range(layers):
            row = [f"{mol['P'][i]:.5E}", f"{mol['T'][i]:.5E}"]
            row += [f"{mol[g][i]:.5E}" for g in gas_keys]
            fw.write(f"<ATMOSPHERE-LAYER-{i+1}>" + ",".join(row) + "\n")
        
        # --- Add extra fixed atmosphere/surface tags ---
        fw.write("\n<ATMOSPHERE-TEMPERATURE>320\n")
        fw.write("<ATMOSPHERE-TAU>0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07\n")
        fw.write("<ATMOSPHERE-NMAX>3\n")
        fw.write("<ATMOSPHERE-LMAX>40\n")
    
        fw.write("<SURFACE-TEMPERATURE>320\n")
        fw.write("<SURFACE-ALBEDO>0.06\n")
        fw.write("<SURFACE-EMISSIVITY>0.95\n")
        fw.write("<SURFACE-GAS-RATIO>1.0\n")
        fw.write("<SURFACE-GAS-UNIT>ratio\n")
        fw.write("<SURFACE-NSURF>0\n")
        fw.write("<SURFACE-MODEL>Lambert\n")
