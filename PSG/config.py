#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 15:30:09 2024

@author: gregcooke

Central configuration for PSG file generation.

This file should contain *no logic* — only parameters.
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Central configuration for PSG file generation.

This file should contain *no logic* — only parameters.
"""

from pathlib import Path
from stars import STARS

star_name = "Sun"
star = STARS[star_name]

CONFIG = {

    # ------------------------------------------------------------------
    # Simulation metadata
    # ------------------------------------------------------------------
    
    "metadata": {
        "name":'Earth 100% PAL VULCAN',
        "model":'VULCAN',
        "boundary":'O2 = 0.21',
        
    },

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------
    "output_dir": Path("/Users/gregcooke/psg_outputs"),
    "output_template": "{atm_stem}_psg.txt",

    # ------------------------------------------------------------------
    # Target selection
    # ------------------------------------------------------------------
    "star": star,

    
    # ------------------------------------------------------------------
    # Target (planet)
    # ------------------------------------------------------------------
    "planet": {
        "name": "Earth",
        "type": "Exoplanet",
        "diameter": 12742,
        "gravity": 9.81,
        "gravity_unit": "g",
        "period": 365.0,
        "eccentricity": 0.0167,
        "inclination": 90.0,
        "season": 180.0,
        "phase": 180.0,
        "obs_velocity": 29.0,
        "periapsis": 0.0,
        "solar_longitude": 0.0,
        "solar_latitude": 0.0,
        "obs_longitude": 0.0,
        "obs_latitude": 0.0,
    },

    # ------------------------------------------------------------------
    # Atmosphere
    # ------------------------------------------------------------------
    "atmosphere": {
        "model": "VULCAN",
        #"model": "Photochem",
        #"model": "Atmos",
        "files": [
            
            Path(
                "/Users/gregcooke/VIH_cases/output/"
                "Earth_1e12s_60SZA_WPT_1rtol.vul"
            ),
            
            #Path(
            #    "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos/100pc_PT_profile/"
            #    "100%PAL_Sun_0.0Ga.dist"
            #),
            
            #Path(
            #    "/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/100pc_PT_profile/"
            #    "100%PAL_Sun_0.0Ga.txt"
            #),
        ],
        "clouds": False, #option to be included later
        "sulfur": False, #option to be included later
    },

    # ------------------------------------------------------------------
    # Observation geometry
    # ------------------------------------------------------------------
    "observation": {
        "transit": True,
        "reflection": False,
        "emission": False,
        "exposure_time": 3600.0,
        "n_exposures": 10,
        "phase": 0.0,
    },

    # ------------------------------------------------------------------
    # Telescope / instrument
    # ------------------------------------------------------------------
    "telescope": "Ideal", #pick telescope - see options in instruments.py file
    
    "geometry": {
        "type": "Observatory",
        "offset_ns": 0.0,
        "offset_ew": 0.0,
        "offset_unit": "arcsec",
        "obs_altitude": 10.0,
        "altitude_unit": "pc",
        "user_param": 0.0,
        "solar_angle": 90.0,
        "obs_angle": 48.121,
        "planet_fraction": 1,
        "star_distance": 1,
        "star_fraction": 0.09959069827155864,
        "phase": 180.0,
        "season": 180.0,
        "ref": "User",
        "disk_angles": 1,
        "rotation": "-0.00,0.00",
        "azimuth": 0.0,
        "brdfscaler": 1.0,
    },
    
    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------
    
    "plot": {
        "transit": True,
        "wavelength_range": (0.2, 10),
        "molecules": ["O3", "H2O"],
        #save name below
        "output_file": "/Users/gregcooke/python_output/planet_spectra.png"
    },

    # ------------------------------------------------------------------
    # PSG / numerical options
    # ------------------------------------------------------------------
    "misc": {
        "units": "um",
        "resolution": 250,
        "noise": True,
        "rad_units": "rkm",
    },
}
