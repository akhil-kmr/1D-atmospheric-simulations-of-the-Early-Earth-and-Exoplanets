# observation.py

def build_observation(transit=True, reflection=False, emission=False,
                      exposure_time=3600.0, n_exposures=1, phase=0.0):
    return {
        "transit": transit,
        "reflection": reflection,
        "emission": emission,
        "exposure_time": exposure_time,
        "n_exposures": n_exposures,
        "phase": phase,
    }


