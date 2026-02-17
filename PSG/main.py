#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config import CONFIG
from atmosphere import process_atmosphere
from observation import build_observation
from instruments import INSTRUMENTS
from psg_writer import build_psg_config
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from plot_spectra import plot_spectra
from plot_spectra import sub
from plot_spectra import filter_molecules_by_wavelength

#%% Instructions

'''
    First make sure the config.py file is set up correctly for
    the planet you are interested in modelling
    Note that noise sources are not correctlt implemented at the moment

    Choose whether to auto upload your file to PSG or not
    auto_upload = False or auto_upload = True
    If you want to play around with plot, chose auto_upload = False
    Choose wavelength range of plot. Usually between 0.2-20 microns 
    has features for Earth, but that may not be the best way to visualise the data
'''

auto_upload = True
wavelength_range = (0.2, 3)

#%% Start main programme

def main():
    cfg = CONFIG
    atm_cfg = cfg["atmosphere"]
    metadata = cfg["metadata"]
    misc = cfg["misc"]
    output_dir = Path(cfg["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    for atm_file in atm_cfg["files"]:
        print(f"Loading atmosphere file: {atm_file}")
        atmosphere = process_atmosphere(
            atm_file,
            clouds=atm_cfg["clouds"],
            sulfur=atm_cfg["sulfur"],
            model=(atm_cfg["model"]),
        )

        planet = cfg["planet"]
        star = cfg["star"]
        observation = build_observation(**cfg["observation"])
        instrument = INSTRUMENTS[cfg["telescope"]]


        # Determine PSG input and output file paths
        atm_stem = atm_file.stem if isinstance(atm_file, Path) else Path(atm_file).stem
        psg_input = output_dir / cfg["output_template"].format(atm_stem=atm_stem)
        psg_output = output_dir / f"{atm_stem}_psg_output.txt"
        
        # Build PSG input lines
        build_psg_config(
            planet=planet,
            star=star,
            geometry=cfg["geometry"],
            atmosphere=atmosphere,
            observation=observation,
            instrument=instrument,
            config=cfg,
            output_file=psg_input  # <-- this writes the file
        )

        print(psg_input)
        print(psg_output)
        api = "https://psg.gsfc.nasa.gov//api.php"
        
        # Run curl and save the output properly
        import os

        if (auto_upload == True):
            api = ' https://psg.gsfc.nasa.gov//api.php >'
            #command = 'curl -s -d app=globes --data-urlencode file@'
            command = 'curl --data-urlencode file@'
            
            #Calls PSG from the command line and outputs the file
            print('Upload to PSG started')
            #os.system(command+psg_input+api+' '+psg_output)
            os.system(f"{command}{psg_input}{api} {psg_output}")

            print('Upload to PSG finished')
        else:
            print('no upload')
        
        print(f"Upload to PSG finished, saved to {psg_output}")
        
        # Load the data, skipping header lines
        
        molecules = {
            metadata['name']: [
                # H2O
                {'label':'H'+sub(2)+'O', 'x':1.42, 'y':82, 'color':'blue',
                 'gradient': {'xrange':(1.41,1.45), 'cmap':plt.cm.Blues, 'alpha':0.3}},
                {'label':'H'+sub(2)+'O', 'x':6, 'y':82, 'color':'blue',
                 'gradient': {'xrange':(5.3,7.4), 'cmap':plt.cm.Blues, 'alpha':0.3}},
                # O2-X
                {'label':'O'+sub(2)+'-X', 'x':6.6, 'y':82, 'color':'red',
                 'gradient': {'xrange':(6,7.1), 'cmap':plt.cm.Reds, 'alpha':0.3}},
                # O2
                {'label':'O'+sub(2), 'x':0.76, 'y':82, 'color':'red',
                 'gradient': {'xrange':(0.75,0.77), 'cmap':plt.cm.Reds, 'alpha':0.3}},
                {'label':'O'+sub(2), 'x':0.68, 'y':85, 'color':'red',
                 'gradient': {'xrange':(0.68,0.7), 'cmap':plt.cm.Reds, 'alpha':0.3}},
                {'label':'O'+sub(2), 'x':1.27, 'y':82, 'color':'red',
                 'gradient': {'xrange':(1.25,1.29), 'cmap':plt.cm.Reds, 'alpha':0.3}},
                # CH4
                {'label':'CH'+sub(4), 'x':3.3, 'y':82, 'color':'orange',
                 'gradient': {'xrange':(3.2, 3.4), 'cmap':plt.cm.Oranges, 'alpha':0.3}},
                {'label':'CH'+sub(4), 'x':7.5, 'y':82, 'color':'orange',
                 'gradient': {'xrange':(7.1, 8.1), 'cmap':plt.cm.Oranges, 'alpha':0.3}},
                # O3
                {'label':'O'+sub(3), 'x':0.6, 'y':82, 'color':'green',
                 'gradient': {'xrange':(0.48, 1), 'cmap':plt.cm.Greens, 'alpha':0.3}},
                {'label':'O'+sub(3), 'x':0.3, 'y':82, 'color':'green',
                 'gradient': {'xrange':(0.2, 0.33), 'cmap':plt.cm.Greens, 'alpha':0.3}},
                {'label':'O'+sub(3), 'x':4.8, 'y':82, 'color':'green',
                 'gradient': {'xrange':(4.65, 5), 'cmap':plt.cm.Greens, 'alpha':0.3}},
                {'label':'O'+sub(3), 'x':9.6, 'y':82, 'color':'green',
                 'gradient': {'xrange':(9, 10), 'cmap':plt.cm.Greens, 'alpha':0.3}},
                # CO2
                {'label':'CO'+sub(2), 'x':15, 'y':82, 'color':'black',
                 'gradient': {'xrange':(13, 18), 'cmap':plt.cm.binary, 'alpha':0.3}},
                {'label':'CO'+sub(2), 'x':4.2, 'y':82, 'color':'black',
                 'gradient': {'xrange':(4, 4.65), 'cmap':plt.cm.binary, 'alpha':0.3}},
                {'label':'CO'+sub(2), 'x':2.7, 'y':82, 'color':'black',
                 'gradient': {'xrange':(2.6, 2.95), 'cmap':plt.cm.binary, 'alpha':0.3}},
                {'label':'CO'+sub(2), 'x':2, 'y':82, 'color':'black',
                 'gradient': {'xrange':(1.9, 2.1), 'cmap':plt.cm.binary, 'alpha':0.3}},
                # N2O
                {'label':'N'+sub(2)+'O', 'x':17, 'y':82, 'color':'purple',
                 'gradient': {'xrange':(16.8, 17.2), 'cmap':plt.cm.Purples, 'alpha':0.3}},
            ]
        }
        
        filtered_molecules = filter_molecules_by_wavelength(
            molecules,
            wavelength_range
        )
        
        plot_spectra(
            files={metadata['name']: [psg_output]},
            figsize=(12, 7),
            transit=True,
            wavelength_range=wavelength_range,
            molecule_options=filtered_molecules,
            output_file=True,
            xscale='linear',
        )

if __name__ == "__main__":
    main()
