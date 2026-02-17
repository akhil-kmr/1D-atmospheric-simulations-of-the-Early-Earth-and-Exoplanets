#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 15:30:09 2024

@author: gregcooke
"""

'''
Instructions to user:
    
    filename is the config file that is produced
    file_path is the file path of the atmospheric file you want get spectra from
    filout is the output config file
    changing these three to the correct inputs and running the script should work
    provided your conda environment is set up correctly
    
    Star options include the Sun, K2-18, and TOI-270
    
    Other changes can be made
    
    
    Need a way to get the correct mean molecular weight
    
'''
from config import CONFIG
cfg = CONFIG
misc = cfg["misc"]


# instruments.py

INSTRUMENTS = {
    
    "Ideal": {
        "instrument_name": "Ideal",
        "description": (
            "User specified instrument"
        ),

        "range": (0.2, 20),
        "range_unit": "um",

        "resolution": misc['resolution'],
        "resolution_unit": "RP",

        "telescope_type": "SINGLE",
        "diameter": 15.0,

        "beam": 1,
        "beam_unit": "diffrac",

        "noise": {
            "model": "CCD",
            "read": "0.2@0.2,0.2@1,2.5@1.01,2.5@2.5",
            "dark": "3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.5",
            "opt_temp": 270,
            "opt_emis": 0.1,
            "pixels": 8,
        },
    },

    "LUVOIR HDI": {
        "instrument_name": "LUVOIR_HDI",
        "description": (
            "The HDI design provides a 2 x 3 arcminute field-of-view "
            "with UVIS (0.2–1.0 µm) and NIR (0.8–2.5 µm) channels."
        ),

        "range": (0.2, 2.5),
        "range_unit": "um",

        "resolution": 125,
        "resolution_unit": "RP",

        "telescope_type": "SINGLE",
        "diameter": 15.0,

        "beam": 1,
        "beam_unit": "diffrac",

        "noise": {
            "model": "CCD",
            "read": "0.2@0.2,0.2@1,2.5@1.01,2.5@2.5",
            "dark": "3e-5@0.2,3e-5@1,2e-3@1.01,2e-3@2.5",
            "opt_temp": 270,
            "opt_emis": 0.1,
            "pixels": 8,
        },
    },

    "LUVOIR HDI 6m": {
        "instrument_name": "LUVOIR_HDI",
        "range": (0.2, 2.5),
        "resolution": 500,
        "diameter": 6.0,
        "telescope_type": "SINGLE",
    },

    # Add more instruments here
}





