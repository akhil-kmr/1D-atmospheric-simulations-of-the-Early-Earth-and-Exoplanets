# io.py

from pathlib import Path
import pickle

def load_atmosphere_file(file_path, model="vulcan"):
    """
    Load atmosphere data from different models.
    Returns a dict compatible with process_atmosphere.
    """
    if model.lower() == "vulcan":
        with open(file_path, "rb") as f:
            raw = pickle.load(f)
        # Wrap in atm dict to match process_atmosphere expectation
        atm_dict = {
            "pco": raw.get("P", raw.get("press")),
            "Tco": raw.get("T", raw.get("temp")),
            "H2": raw.get("H2"),
            "H2O": raw.get("H2O"),
            "CH4": raw.get("CH4"),
            "CO2": raw.get("CO2"),
            "N2": raw.get("N2"),
            "NH3": raw.get("NH3"),
            "O2": raw.get("O2"),
            "O3": raw.get("O3"),
            "N2O": raw.get("N2O"),
        }
        return {"model": "vulcan", "atm": atm_dict}

    elif model.lower() == "waccm":
        # WACCM files could be netCDF or another object
        import netCDF4
        nc = netCDF4.Dataset(file_path)
        return nc  # process_atmosphere knows how to handle WACCM objects

    elif model.lower() == "photochem":
        # Implement Photochem loader here
        import h5py
        with h5py.File(file_path, "r") as f:
            return dict(f)  # wrap as needed

    else:
        raise ValueError(f"Unknown atmosphere model: {model}")



def save_psg_file(path, lines):
    """
    Write PSG input file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        for line in lines:
            f.write(line + "\n")

