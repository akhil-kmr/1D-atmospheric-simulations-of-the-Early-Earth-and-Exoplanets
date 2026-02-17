"""
PSG input file writer.

This module converts structured Python dictionaries
(target, geometry, atmosphere, observation, instrument, config)
into PSG-formatted text lines.

No science logic should live here — only formatting.
"""
import numpy as np

def build_psg_config(
    planet,
    star,
    geometry,
    atmosphere,
    observation,
    instrument,
    config,
    output_file=None,  # <--- new argument
):
    """
    Build full PSG configuration file.

    Parameters
    ----------
    output_file : Path or str, optional
        If given, write the PSG input to this file.

    Returns
    -------
    lines : list of str
        PSG input lines
    """
    lines = []

    write_object(lines, planet, star)
    write_geometry(lines, geometry)
    write_atmosphere(lines, atmosphere)
    write_instrument(lines, instrument, observation, config)
    write_psg_options(lines, config)

    # If a file path is given, write lines to disk
    if output_file is not None:
        with open(output_file, "w") as f:
            for line in lines:
                f.write(line + "\n")

    return lines


# ---------------------------------------------------------------------
# Object (planet + star)
# ---------------------------------------------------------------------

def write_object(lines, planet, star):
    """
    Write OBJECT and OBJECT-STAR fields.
    """

    # -----------------------------
    # Object type
    # -----------------------------
    lines.append("<OBJECT>Exoplanet")

    # -----------------------------
    # Planet fields
    # -----------------------------
    planet_map = {
        "name": "NAME",
        "date": "DATE",
        "diameter": "DIAMETER",
        "gravity": "GRAVITY",
        "gravity_unit": "GRAVITY-UNIT",
        "period": "PERIOD",
        "periapsis": "PERIAPSIS",
        "inclination": "INCLINATION",
        "obs_velocity": "OBS-VELOCITY",
        "season": "SEASON",
        "solar_longitude": "SOLAR-LONGITUDE",
        "solar_latitude": "SOLAR-LATITUDE",
        "obs_longitude": "OBS-LONGITUDE",
        "obs_latitude": "OBS-LATITUDE",
        #"eccentricity": "ECCENTRICITY",
        #"phase": "PHASE",
    }
    
    

    for key, tag in planet_map.items():
        if key in planet:
            lines.append(f"<OBJECT-{tag}>{planet[key]}")

    # -----------------------------
    # Star fields (THIS is what you need)
    # -----------------------------
    star_map = {
        "type": "STAR-TYPE",
        "temperature": "STAR-TEMPERATURE",
        "radius": "STAR-RADIUS",
        "metallicity": "STAR-METALLICITY",
        "distance": "STAR-DISTANCE",      # ✅ REQUIRED
        "velocity": "STAR-VELOCITY",      # ✅ REQUIRED
    }

    for key, tag in star_map.items():
        if key in star:
            lines.append(f"<OBJECT-{tag}>{star[key]}")



# ---------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------

def write_geometry(lines, geometry):
    """
    Write all <GEOMETRY-*> fields
    """
    defaults = {
        "type": "Observatory",
        "offset_ns": 0.0,
        "offset_ew": 0.0,
        "offset_unit": "arcsec",
        "obs_altitude": 10.0,
        "altitude_unit": "pc",
        "user_param": 0.0,
        "stellar_type": "G",
        "stellar_temperature": 5777,
        "stellar_magnitude": 0.0,
        "solar_angle": 90.0,
        "obs_angle": 90.0,
        "planet_fraction": 1.0,
        "star_distance": 0.09959069827155864,
        "star_fraction": 0,
        "phase": 90.0,
        "ref": "User",
        "disk_angles": 1,
        "rotation": "-0.00,0.00",
        "azimuth": 0.0,
        "brdfscaler": 1.0,
    }

    for key, default in defaults.items():
        value = geometry.get(key, default)
        if (key == 'type'):
            lines.append("<GEOMETRY>")
        else:
            lines.append(f"<GEOMETRY-{key.replace('_', '-').upper()}>{value}")

    # Generic <GEOMETRY> type
    lines.append(f"<GEOMETRY>{geometry.get('type', defaults['type'])}")


# ---------------------------------------------------------------------
# Atmosphere
# ---------------------------------------------------------------------

def write_atmosphere(lines, atmosphere):

    mol = atmosphere["molecules"]
    layers = atmosphere["layers"]

    lines.append("<ATMOSPHERE-DESCRIPTION>Custom Atmosphere")
    lines.append("<ATMOSPHERE-STRUCTURE>Equilibrium")
    lines.append(f"<ATMOSPHERE-LAYERS>{layers}")
    lines.append("<ATMOSPHERE-PUNIT>bar")
    lines.append(f"<ATMOSPHERE-WEIGHT>{np.mean(atmosphere['MMW'])}")

    lines.append("<ATMOSPHERE-CONTINUUM>Rayleigh,Refraction,CIA_all,UV_all")
    lines.append(f"<ATMOSPHERE-NGAS>{atmosphere['NGAS']}")
    lines.append(f"<ATMOSPHERE-GAS>{atmosphere['GASES']}")
    lines.append(f"<ATMOSPHERE-TYPE>{atmosphere['HIT']}")
    lines.append(f"<ATMOSPHERE-ABUN>{atmosphere['ABUN']}")
    lines.append(f"<ATMOSPHERE-UNIT>{atmosphere['ATM_UNIT']}")

    gas_keys = list(mol.keys())
    mol_keys = list(mol.keys())[2:]
    lines.append(
        "<ATMOSPHERE-LAYERS-MOLECULES>" + ",".join(mol_keys)
    )

    for i in range(layers):
        vals = []
        for g in gas_keys:
            vals.append(f"{mol[g][i]:.5E}")
        lines.append(f"<ATMOSPHERE-LAYER-{i+1}>" + ",".join(vals))
    
    lines.append("<ATMOSPHERE-PRESSURE>1")
    lines.append("<ATMOSPHERE-TEMPERATURE>320")
    lines.append("<ATMOSPHERE-TAU>0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07")
    lines.append("<ATMOSPHERE-NMAX>3")
    lines.append("<ATMOSPHERE-LMAX>40")
    lines.append("<SURFACE-TEMPERATURE>320")
    lines.append("<SURFACE-ALBEDO>0.3")
    lines.append("<SURFACE-EMISSIVITY>0.95")
    lines.append("<SURFACE-GAS-RATIO>1.0")
    lines.append("<SURFACE-GAS-UNIT>ratio")
    lines.append("<SURFACE-NSURF>0")
    lines.append("<SURFACE-MODEL>Lambert")



# ---------------------------------------------------------------------
# Instrument / telescope
# ---------------------------------------------------------------------

def write_instrument(lines, instrument, observation, config):
    lines.append(f"<GENERATOR-INSTRUMENT>{instrument.get('instrument_name','user')}")
    r1, r2 = instrument.get("range",(0,0))
    lines.append(f"<GENERATOR-RANGE1>{r1}")
    lines.append(f"<GENERATOR-RANGE2>{r2}")
    lines.append(f"<GENERATOR-RANGEUNIT>{instrument.get('range_unit','um')}")
    lines.append(f"<GENERATOR-RESOLUTION>{instrument.get('resolution',0)}")
    lines.append(f"<GENERATOR-RESOLUTIONUNIT>{instrument.get('resolution_unit','RP')}")
    lines.append(f"<GENERATOR-TELESCOPE>{instrument.get('telescope_type','Generic')}")
    lines.append(f"<GENERATOR-DIAMTELE>{instrument.get('diameter',0)}")
    beam = instrument.get("beam")
    if beam is not None:
        lines.append(f"<GENERATOR-BEAM>{beam}")
        lines.append(f"<GENERATOR-BEAM-UNIT>arcsec")
    # Optional noise
    if config["misc"].get("noise", False):
        noise = instrument.get("noise", {})
        lines.append(f"<GENERATOR-NOISE>{noise.get('model','CCD')}")
        lines.append(f"<GENERATOR-NOISETIME>{observation.get('exposure_time',0.0)}")
        lines.append(f"<GENERATOR-NOISEFRAMES>{observation.get('n_exposures',1)}")
        for field, tag in {
            "dark_current":"NOISE1",
            "read_noise":"NOISE2",
            "optical_temp":"NOISEOTEMP",
            "optical_efficiency":"NOISEOEFF",
            "optical_emissivity":"NOISEOEMIS",
            "pixels":"NOISEPIXELS"
        }.items():
            if field in noise:
                lines.append(f"<GENERATOR-{tag}>{noise[field]}")


# ---------------------------------------------------------------------
# PSG physics & numerical options
# ---------------------------------------------------------------------
def write_psg_options(lines, config):

    psg = config["misc"]

    # Radiative transfer
    if "gas_model" in psg:
        lines.append(f"<GENERATOR-GAS-MODEL>{psg['gas_model']}")
    if "continuum_model" in psg:
        lines.append(f"<GENERATOR-CONT-MODEL>{psg['continuum_model']}")
    if "stellar_continuum" in psg:
        lines.append(
            f"<GENERATOR-CONT-STELLAR>{psg['stellar_continuum']}"
        )
    if "rad_units" in psg:
        lines.append(f"<GENERATOR-RADUNITS>{psg['rad_units']}")
    if "lograd" in psg:
        lines.append(f"<GENERATOR-LOGRAD>{psg['lograd']}")
    if "resolution_kernel" in psg:
        lines.append(
            f"<GENERATOR-RESOLUTIONKERNEL>{psg['resolution_kernel']}"
        )
    # Transmission presets
    if "trans" in psg:
        lines.append(f"<GENERATOR-TRANS>{psg['trans']}")

    if "trans_apply" in psg:
        lines.append(
            f"<GENERATOR-TRANS-APPLY>{psg['trans_apply']}"
        )
    if "trans_show" in psg:
        lines.append(
            f"<GENERATOR-TRANS-SHOW>{psg['trans_show']}"
        )
