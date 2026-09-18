#!/usr/bin/env python3
"""Build with make and run the TUV executable.

Install atmosphere columns from WACCM, VULCAN, Photochem, Atmos, or Kasting 1D,
run TUV over a solar-zenith-angle sweep, and write comparison outputs.

Atmosphere files use 48.2° SZA equilibrium outputs (100% PAL Atmos uses SZA_48.5).
TUV computes surface UV over SZA = 0–85° every 5°.

Examples
--------
- ``python run_TUV.py --compare`` — all six models × four PAL cases (incl. WACCM equator).
- ``python run_TUV.py --compare --model photochem,atmos`` — 1D models only.
- ``python run_TUV.py --compare --cases pal0p1,pal1`` — selected cases only.
- ``python run_TUV.py --compare --wstart 100`` — spectral window from 100–420 nm.
- ``python run_TUV.py --model vulcan --case pal0p1`` — one VULCAN run.
- ``python run_TUV.py --cam-only --case pal10`` — write ``ussa.*`` only (WACCM).
- ``python run_TUV.py --list-cases`` — registered WACCM NetCDF paths.

Outputs land in ``OUTPUT/compare/`` as ``.txt`` (TUV native), ``.csv`` (parsed),
and ``comparison_summary.csv`` / ``comparison_spectral.csv`` for side-by-side use.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import pickle
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import xarray as xr

# ---------------------------------------------------------------------------
# USER: paths and toggles
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent

USE_CAM_ATM = True
CESM_DATA = Path("/Users/gregcooke/CESM_data")
VIH_OUTPUT_PATH = Path("/Users/gregcooke/VIH_cases/output")
ATMOS_BASE = Path("/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos")
PHOTOCHEM_BASE = Path("/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem")
KASTING_BASE = Path("/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Kasting_1D_model")
OUTPUT_COMPARE_DIR = ROOT / "OUTPUT" / "compare"

ALL_COMPARISON_MODELS = (
    "waccm", "waccm_eq", "vulcan", "photochem", "atmos", "kasting",
)
WACCM_EQ_ILAT = 47  # nearest model latitude to equator (~-0.95 deg on f19 grid)

# VULCAN .vul filenames encode the photochemical equilibrium SZA (48.2°).
VULCAN_ATM_SZA_DEG = 48.2

# TUV illumination grid (lzenit = T): independent of atmosphere source.
SZA_SWEEP_START = 0.0
SZA_SWEEP_STOP = 85.0
SZA_SWEEP_STEP = 5.0

# TUV spectral window (rdinp enforces wstart >= 100 nm). Default matches legacy runs.
DEFAULT_WSTART_NM = 280.0
DEFAULT_WSTOP_NM = 420.0
CURRENT_WSTART_NM = DEFAULT_WSTART_NM
CURRENT_WSTOP_NM = DEFAULT_WSTOP_NM

DU_PER_COL = 2.687e16


def sza_sweep_count() -> int:
    n = int(round((SZA_SWEEP_STOP - SZA_SWEEP_START) / SZA_SWEEP_STEP)) + 1
    if n < 1:
        raise ValueError("SZA sweep settings must produce at least one angle.")
    return n


def set_wavelength_range(wstart_nm: float, wstop_nm: float) -> None:
    """Set the TUV spectral window used for subsequent usrinp / archive naming."""
    global CURRENT_WSTART_NM, CURRENT_WSTOP_NM
    if wstart_nm < 100.0:
        raise ValueError(f"TUV requires wstart >= 100 nm (got {wstart_nm})")
    if wstop_nm <= wstart_nm:
        raise ValueError(f"wstop ({wstop_nm}) must be > wstart ({wstart_nm})")
    CURRENT_WSTART_NM = float(wstart_nm)
    CURRENT_WSTOP_NM = float(wstop_nm)


def wavelength_nwint(wstart_nm: float | None = None, wstop_nm: float | None = None) -> int:
    """Number of wavelength intervals for usrinp ``nwint``.

    Equal spacing is only valid for windows entirely at/above 205.8 nm. For any
    window that extends into the Lyman-alpha / SRB range, TUV requires the
    pre-specified atmospheric grid ``nwint = -156`` (``DATAE1/GRIDS/combined.grid``,
    which spans ~120–735 nm).
    """
    ws = CURRENT_WSTART_NM if wstart_nm is None else wstart_nm
    we = CURRENT_WSTOP_NM if wstop_nm is None else wstop_nm
    if we <= ws:
        raise ValueError(f"Invalid wavelength window {ws}–{we} nm")
    if ws < 205.8:
        return -156
    n = int(round(we - ws))
    if n < 1:
        raise ValueError(f"Invalid wavelength window {ws}–{we} nm")
    if n > 999:
        raise ValueError(f"nwint={n} exceeds TUV kw-1=999; narrow the wavelength window")
    return n


def wavelength_tag(wstart_nm: float | None = None, wstop_nm: float | None = None) -> str:
    ws = CURRENT_WSTART_NM if wstart_nm is None else wstart_nm
    we = CURRENT_WSTOP_NM if wstop_nm is None else wstop_nm
    return f"w{int(round(ws))}-{int(round(we))}"


def comparison_summary_path(wstart_nm: float | None = None, wstop_nm: float | None = None) -> Path:
    """Primary summary CSV for the active (or given) wavelength window."""
    tag = wavelength_tag(wstart_nm, wstop_nm)
    return OUTPUT_COMPARE_DIR / f"comparison_summary_{tag}.csv"


def meta_matches_wavelength(meta: dict, wstart_nm: float | None = None, wstop_nm: float | None = None) -> bool:
    """Match archived run meta to a wavelength window (legacy metas → 280–420)."""
    ws = CURRENT_WSTART_NM if wstart_nm is None else wstart_nm
    we = CURRENT_WSTOP_NM if wstop_nm is None else wstop_nm
    meta_ws = float(meta.get("wstart_nm", DEFAULT_WSTART_NM))
    meta_we = float(meta.get("wstop_nm", DEFAULT_WSTOP_NM))
    return abs(meta_ws - ws) < 0.05 and abs(meta_we - we) < 0.05

DATASET_CASES: dict[str, Path] = {
    # Newer CESM_data WACCM6 runs from Early_Earth.py (NOT the older H_escape *_pc_h0 set).
    # Comparison cases use: Pre_h0, Ten_h0, One_h0, Zero1_h0.
    "pal150": CESM_DATA / "Earth_150pc_o2.cam.h0.0034-0037.nc",  # One50_h0
    "baseline": CESM_DATA / "Earth_100pc_o2.cam.h0.0009-0012.nc",  # Pre_h0  (pal100)
    "pal50": CESM_DATA / "Earth_50pc_o2.cam.h0.0040-0043.nc",  # Fifty_h0
    "pal10": CESM_DATA / "Earth_10pc_o2.cam.h0.0037-0040.nc",  # Ten_h0
    "pal5": CESM_DATA / "Earth_5pc_o2.cam.h0.0048-0051.nc",  # Five_h0
    "pal1": CESM_DATA / "Earth_1pc_o2.cam.h0.0045-0048.nc",  # One_h0
    "pal0p5": CESM_DATA / "Earth_0.5pc_o2.cam.h0.0055-0058.nc",  # Zero5_h0
    "pal0p5_ys": CESM_DATA / "Earth_0.5pc_o2_YS.cam.h0.0057-0060.nc",  # Zero5_h0_YS
    "pal0p1": CESM_DATA / "Earth_0.1pc_o2.cam.h0.0033-0036.nc",  # Zero1_h0
    "pal0p1_ys": CESM_DATA / "Earth_0.1pc_o2_YS.cam.h0.0261-0264.nc",  # Zero1_h0_YS
}

DEFAULT_CAM_CASE = "pal0p1"
CAM_FILE = DATASET_CASES[DEFAULT_CAM_CASE]

# PAL comparison cases (see Early_Earth.py for file naming).
COMPARISON_CASES: dict[str, dict] = {
    "pal0p1": {
        "label": "0.1% PAL",
        "waccm_case": "pal0p1",
        "vulcan_file": "Earth_0.1pc_o2_1e12s_48.2SZA_WPT_1rtol.vul",
        "photochem_file": "0.1pc/Earth_0.1pc_48.2.txt",
        "atmos_ptz": "0.1pc/SZA_48.2/PTZ_mixingratios_out.dist",
        "kasting_plot": "0.1pc/SZA_48.2/OUTPUT_PLOT.dat",
        "outfil_waccm": "Wp01W",
        "outfil_waccm_eq": "Wp01E",
        "outfil_vulcan": "Vp01V",
        "outfil_photochem": "Pp01P",
        "outfil_atmos": "Ap01A",
        "outfil_kasting": "Kp01K",
    },
    "pal1": {
        "label": "1% PAL",
        "waccm_case": "pal1",
        "vulcan_file": "Earth_1pc_o2_1e12s_48.2SZA_WPT_1rtol.vul",
        "photochem_file": "1pc/Earth_1pc_48.2.txt",
        "atmos_ptz": "1pc/SZA_48.2/PTZ_mixingratios_out.dist",
        "kasting_plot": "1pc/SZA_48.2/OUTPUT_PLOT.dat",
        "outfil_waccm": "W01PW",
        "outfil_waccm_eq": "W01PE",
        "outfil_vulcan": "V01PV",
        "outfil_photochem": "P01PP",
        "outfil_atmos": "A01PA",
        "outfil_kasting": "K01PK",
    },
    "pal10": {
        "label": "10% PAL",
        "waccm_case": "pal10",
        "vulcan_file": "Earth_10pc_o2_1e12s_48.2SZA_WPT_1rtol.vul",
        "photochem_file": "10pc/Earth_10pc_48.2.txt",
        "atmos_ptz": "10pc/SZA_48.2/PTZ_mixingratios_out.dist",
        "kasting_plot": "10pc/SZA_48.2/OUTPUT_PLOT.dat",
        "outfil_waccm": "W10PW",
        "outfil_waccm_eq": "W10PE",
        "outfil_vulcan": "V10PV",
        "outfil_photochem": "P10PP",
        "outfil_atmos": "A10PA",
        "outfil_kasting": "K10PK",
    },
    "pal100": {
        "label": "100% PAL",
        "waccm_case": "baseline",
        "vulcan_file": "Earth_1e12s_48.2SZA_WPT_1rtol.vul",
        "photochem_file": "100pc/Earth_100pc_48.2.txt",
        "atmos_ptz": "100pc/SZA_48.5/PTZ_mixingratios_out.dist",
        "kasting_plot": "100pc/SZA_48.2/OUTPUT_PLOT.dat",
        "outfil_waccm": "W100W",
        "outfil_waccm_eq": "W100E",
        "outfil_vulcan": "V100V",
        "outfil_photochem": "P100P",
        "outfil_atmos": "A100A",
        "outfil_kasting": "K100K",
    },
}

HORIZONTAL_MODE = "column"
TARGET_LAT_DEG = 45.0
LWAV_AVERAGE_TIME = True
GW_REFERENCE_FILE: Path | None = None

VAR_T = "T"
VAR_O3 = "O3"
VAR_Z3 = "Z3"
VAR_P = "lev"

O3_IS_PPMV = False
PRESSURE_IN_HPA = True

ITIME = 0
ILAT = 0
ILON = 0

DIM_TIME = "time"
DIM_LAT = "lat"
DIM_LON = "lon"

REQUIRE_EXPLICIT_PRESSURE = True

_OZONE_FILE_Z_KM = (
    0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0,
    26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0, 46.0, 48.0, 50.0,
    52.0, 54.0, 56.0, 58.0, 60.0, 62.0, 64.0, 66.0, 68.0, 70.0, 72.0, 74.0,
)

K_B = 1.380649e-23
TEMP_PROFILE_TOP_KM = 120.0

DOSE_COLUMNS = [
    "UVB_280_315",
    "UVBstar_280_320",
    "UVA_315_400",
    "DNA_damage",
    "UV_index",
    "P_Dam_C_1971",
    "P_Dam_FC_2003",
    "P_Dam_FC_2003_ext390",
]


def pressure_from_dataset(ds) -> object:
    raise NotImplementedError("Define pressure_from_dataset or set VAR_P / flags.")


def vmr_to_o3_number_density(o3_vmr: object, p_pa: object, t_k: object) -> np.ndarray:
    vmr = np.asarray(o3_vmr, dtype=np.float64)
    if O3_IS_PPMV:
        vmr = vmr / 1.0e6
    p = np.asarray(p_pa, dtype=np.float64)
    tk = np.asarray(t_k, dtype=np.float64)
    n_air_m3 = p / (K_B * tk)
    n_air_cm3 = n_air_m3 / 1.0e6
    return vmr * n_air_cm3


def _load_gw_lat(ds_cam: xr.Dataset, cam_path: Path) -> np.ndarray:
    path = GW_REFERENCE_FILE if GW_REFERENCE_FILE is not None else cam_path
    with xr.open_dataset(path, decode_times=False) as gds:
        if "gw" not in gds:
            raise KeyError(
                f"No variable 'gw' in {path}. Set GW_REFERENCE_FILE to a grid "
                "dataset that contains ``gw`` on the same latitude dimension as CAM_FILE."
            )
        gw = np.asarray(gds["gw"].values, dtype=np.float64).ravel()
    nlat = int(ds_cam.sizes.get(DIM_LAT, -1))
    if nlat < 0:
        raise KeyError(f"Dataset has no dimension {DIM_LAT!r} (DIM_LAT).")
    if gw.size != nlat:
        raise ValueError(
            f"gw length {gw.size} does not match CAM {DIM_LAT} size {nlat}. "
            "Use a GW_REFERENCE_FILE on the same grid as CAM_FILE."
        )
    return gw


def _gaussian_lat_mean(da: xr.DataArray, gw: np.ndarray, lat_dim: str) -> xr.DataArray:
    if lat_dim not in da.dims:
        return da
    gw_da = xr.DataArray(
        gw,
        dims=[lat_dim],
        coords={lat_dim: da.coords[lat_dim]},
    )
    return (da * gw_da).sum(dim=lat_dim) / float(np.sum(gw))


def lwav(
    da: xr.DataArray,
    gw: np.ndarray,
    lat_dim: str,
    lon_dim: str,
    time_dim: str,
    average_time: bool,
) -> xr.DataArray:
    out = da
    if lon_dim in out.dims:
        out = out.mean(dim=lon_dim)
    if average_time and time_dim in out.dims:
        out = out.mean(dim=time_dim)
    elif time_dim in out.dims:
        out = out.isel({time_dim: ITIME})
    return _gaussian_lat_mean(out, gw, lat_dim)


def _reduce_to_lev_column(
    ds: xr.Dataset,
    var_name: str,
    cam_path: Path,
    *,
    target_lat_deg: float | None = None,
    ilat: int | None = None,
    horizontal_mode: str | None = None,
) -> xr.DataArray:
    if var_name not in ds:
        raise KeyError(f"Missing variable {var_name!r} in dataset.")

    da = ds[var_name]
    has_lat = DIM_LAT in da.dims
    has_lon = DIM_LON in da.dims
    mode = horizontal_mode if horizontal_mode is not None else HORIZONTAL_MODE

    if not has_lat and not has_lon:
        out = da
        if DIM_TIME in out.dims:
            out = out.isel({DIM_TIME: ITIME})
        return out.squeeze(drop=True)

    if ilat is not None:
        isel_kw = {}
        for dim_name, idx in (
            (DIM_TIME, ITIME),
            (DIM_LAT, ilat),
            (DIM_LON, ILON),
        ):
            if dim_name in da.dims:
                isel_kw[dim_name] = idx
        out = da.isel(isel_kw)

    elif target_lat_deg is not None:
        isel_kw = {}
        for dim_name, idx in (
            (DIM_TIME, ITIME),
            (DIM_LON, ILON),
        ):
            if dim_name in da.dims:
                isel_kw[dim_name] = idx
        out = da.isel(isel_kw)
        if DIM_LAT in out.dims:
            out = out.sel({DIM_LAT: target_lat_deg}, method="nearest")

    elif mode == "column":
        isel_kw = {}
        for dim_name, idx in (
            (DIM_TIME, ITIME),
            (DIM_LAT, ILAT),
            (DIM_LON, ILON),
        ):
            if dim_name in da.dims:
                isel_kw[dim_name] = idx
        out = da.isel(isel_kw)

    elif mode == "zonal_mean":
        out = da
        if DIM_TIME in out.dims:
            out = out.isel({DIM_TIME: ITIME})
        if DIM_LON in out.dims:
            out = out.mean(dim=DIM_LON)
        if DIM_LAT not in out.dims:
            raise ValueError("zonal_mean requires a latitude dimension after lon mean.")
        out = out.sel({DIM_LAT: TARGET_LAT_DEG}, method="nearest")

    elif mode == "gaussian_lat_mean":
        gw = _load_gw_lat(ds, cam_path)
        out = lwav(da, gw, DIM_LAT, DIM_LON, DIM_TIME, LWAV_AVERAGE_TIME)
    else:
        raise ValueError(
            f"Unknown horizontal_mode={mode!r}; "
            'use "column", "zonal_mean", or "gaussian_lat_mean".'
        )

    remaining = [d for d in out.dims if out.sizes[d] > 1]
    if len(remaining) != 1:
        raise ValueError(
            f"After horizontal reduction, {var_name!r} has dims {out.dims!r} "
            f"sizes {dict(out.sizes)}; expected one vertical dimension."
        )
    return out.squeeze(drop=True)


def read_cam_vertical_profile(
    cam_path: Path,
    *,
    target_lat_deg: float | None = None,
    ilat: int | None = None,
    horizontal_mode: str | None = None,
):
    with xr.open_dataset(cam_path, decode_times=False) as ds:
        reduce_kw = {
            "target_lat_deg": target_lat_deg,
            "ilat": ilat,
            "horizontal_mode": horizontal_mode,
        }
        t_da = _reduce_to_lev_column(ds, VAR_T, cam_path, **reduce_kw)
        z_da = _reduce_to_lev_column(ds, VAR_Z3, cam_path, **reduce_kw)
        o3_da = _reduce_to_lev_column(ds, VAR_O3, cam_path, **reduce_kw)

        if VAR_P in ds:
            p_da = _reduce_to_lev_column(ds, VAR_P, cam_path, **reduce_kw)
        elif not REQUIRE_EXPLICIT_PRESSURE:
            p_da = pressure_from_dataset(ds)
        else:
            raise KeyError(
                f"Pressure variable {VAR_P!r} not found. Set VAR_P, set "
                "REQUIRE_EXPLICIT_PRESSURE False, or implement pressure_from_dataset."
            )

    def vertical_values(da):
        return np.asarray(da.values, dtype=np.float64).ravel()

    t_k = vertical_values(t_da)
    z_m = vertical_values(z_da)
    o3 = vertical_values(o3_da)
    p_raw = vertical_values(p_da)

    p_pa = p_raw * 100.0 if PRESSURE_IN_HPA else p_raw
    z_km = z_m / 1000.0

    order = np.argsort(z_km)
    z_km = z_km[order]
    t_k = t_k[order]
    p_pa = p_pa[order]
    o3 = o3[order]

    if np.any(np.diff(z_km) <= 0):
        uz, inv, cnt = np.unique(z_km, return_inverse=True, return_counts=True)
        acc_t = np.zeros_like(uz)
        acc_p = np.zeros_like(uz)
        acc_o = np.zeros_like(uz)
        np.add.at(acc_t, inv, t_k)
        np.add.at(acc_p, inv, p_pa)
        np.add.at(acc_o, inv, o3)
        cnt = cnt.astype(np.float64)
        z_km, t_k, p_pa, o3 = uz, acc_t / cnt, acc_p / cnt, acc_o / cnt

    return z_km, t_k, p_pa, o3


def read_vulcan_vertical_profile(vul_path: Path):
    """Return sorted z_km, t_k, n_o3_cm3 from a VULCAN .vul pickle."""
    with open(vul_path, "rb") as f:
        ds = pickle.load(f)

    spec = ds["variable"]["species"]
    o3_idx = spec.index("O3")
    n_o3 = np.asarray(ds["variable"]["y"][:, o3_idx], dtype=np.float64)
    z_km = np.asarray(ds["atm"]["zmco"], dtype=np.float64) / 1e5
    t_k = np.asarray(ds["atm"]["Tco"], dtype=np.float64)

    order = np.argsort(z_km)
    z_km = z_km[order]
    t_k = t_k[order]
    n_o3 = n_o3[order]

    if np.any(np.diff(z_km) <= 0):
        uz, inv, cnt = np.unique(z_km, return_inverse=True, return_counts=True)
        acc_t = np.zeros_like(uz)
        acc_o = np.zeros_like(uz)
        np.add.at(acc_t, inv, t_k)
        np.add.at(acc_o, inv, n_o3)
        cnt = cnt.astype(np.float64)
        z_km, t_k, n_o3 = uz, acc_t / cnt, acc_o / cnt

    return z_km, t_k, n_o3


def _fix_atmos_exponents(line: str) -> str:
    """Fix Atmos exponent typos like 5.004-232 -> 5.004E-232."""
    return re.sub(r"(\d\.\d+)([-+]\d+)", r"\1E\2", line)


def read_atmos_vertical_profile(atmos_path: Path):
    """Return z_km, t_k, p_bar, o3_vmr, n_o3_cm3 from Atmos PTZ_mixingratios_out.dist."""
    with open(atmos_path, encoding="ascii", errors="replace") as f:
        lines = f.readlines()
    header = lines[0]
    data_lines = [_fix_atmos_exponents(ln) for ln in lines[1:]]
    from io import StringIO
    import pandas as pd

    df = pd.read_csv(StringIO(header + "".join(data_lines)), sep=r"\s+")
    # Atmos PTZ_mixingratios_out.dist uses ALT units of centimeters in the files
    # (e.g. 2.500E+04 => 0.25 km). Convert cm -> km.
    z_km = np.asarray(df["ALT"], dtype=np.float64) / 1e5
    t_k = np.asarray(df["TEMP"], dtype=np.float64)
    p_bar = np.asarray(df["PRESS"], dtype=np.float64)
    o3_vmr = np.asarray(df["O3"], dtype=np.float64)
    p_pa = p_bar * 1e5
    n_air = p_pa / (K_B * t_k) / 1e6
    n_o3 = o3_vmr * n_air

    order = np.argsort(z_km)
    z_km, t_k, p_bar, o3_vmr, n_o3 = (
        z_km[order], t_k[order], p_bar[order], o3_vmr[order], n_o3[order]
    )
    return z_km, t_k, p_bar, o3_vmr, n_o3


def read_photochem_vertical_profile(photo_path: Path):
    """Return z_km, t_k, p_bar, o3_vmr, n_o3_cm3 from Photochem Earth_*_48.2.txt."""
    import pandas as pd

    df = pd.read_csv(photo_path, sep=r"\s+")
    z_km = np.asarray(df["alt"], dtype=np.float64)
    t_k = np.asarray(df["temp"], dtype=np.float64)
    p_bar = np.asarray(df["press"], dtype=np.float64)
    o3_vmr = np.asarray(df["O3"], dtype=np.float64)
    n_air = np.asarray(df["den"], dtype=np.float64)
    n_o3 = o3_vmr * n_air

    order = np.argsort(z_km)
    z_km, t_k, p_bar, o3_vmr, n_o3 = (
        z_km[order], t_k[order], p_bar[order], o3_vmr[order], n_o3[order]
    )
    return z_km, t_k, p_bar, o3_vmr, n_o3


def read_kasting_vertical_profile(kasting_path: Path):
    """Return z_km, t_k, p_bar, o3_vmr, n_o3_cm3 from Kasting OUTPUT_PLOT.dat."""
    import pandas as pd

    df = pd.read_csv(kasting_path, delim_whitespace=True, engine="python")
    # Z uses centimeters in Kasting OUTPUT_PLOT.dat (same convention as Atmos ALT).
    z_km = np.asarray(df["Z"], dtype=np.float64) / 1e5
    # PRESS is dyne cm^-2; divide by 1e6 for bar.
    p_bar = np.asarray(df["PRESS"], dtype=np.float64) / 1e6
    den_cm3 = np.asarray(df["DEN"], dtype=np.float64)
    o3_vmr = np.asarray(df["FO3"], dtype=np.float64)
    n_o3 = o3_vmr * den_cm3
    p_pa = p_bar * 1e5
    n_air_m3 = den_cm3 * 1e6
    t_k = np.where(n_air_m3 > 0, p_pa / (n_air_m3 * K_B), np.nan)

    order = np.argsort(z_km)
    z_km, t_k, p_bar, o3_vmr, n_o3 = (
        z_km[order], t_k[order], p_bar[order], o3_vmr[order], n_o3[order]
    )
    return z_km, t_k, p_bar, o3_vmr, n_o3


def read_waccm_vertical_profile(
    cam_path: Path,
    *,
    target_lat_deg: float | None = None,
    ilat: int | None = None,
    horizontal_mode: str | None = None,
):
    """Return z_km, t_k, p_bar, o3_vmr, n_o3_cm3 from WACCM/CAM NetCDF."""
    z_km, t_k, p_pa, o3_vmr = read_cam_vertical_profile(
        cam_path,
        target_lat_deg=target_lat_deg,
        ilat=ilat,
        horizontal_mode=horizontal_mode,
    )
    n_o3 = vmr_to_o3_number_density(o3_vmr, p_pa, t_k)
    return z_km, t_k, p_pa / 1e5, o3_vmr, n_o3


def read_vulcan_vertical_profile_full(vul_path: Path):
    """Return z_km, t_k, p_bar, o3_vmr, n_o3_cm3 from VULCAN .vul."""
    with open(vul_path, "rb") as f:
        ds = pickle.load(f)
    spec = ds["variable"]["species"]
    o3_idx = spec.index("O3")
    n_o3 = np.asarray(ds["variable"]["y"][:, o3_idx], dtype=np.float64)
    z_km = np.asarray(ds["atm"]["zmco"], dtype=np.float64) / 1e5
    t_k = np.asarray(ds["atm"]["Tco"], dtype=np.float64)
    p_bar = np.asarray(ds["atm"]["pco"], dtype=np.float64) / 1e6

    order = np.argsort(z_km)
    z_km, t_k, p_bar, n_o3 = z_km[order], t_k[order], p_bar[order], n_o3[order]
    p_pa = p_bar * 1e5
    n_air = p_pa / (K_B * t_k) / 1e6
    o3_vmr = np.where(n_air > 0, n_o3 / n_air, 0.0)
    return z_km, t_k, p_bar, o3_vmr, n_o3


def compute_o3_column_du(z_km: np.ndarray, n_o3_cm3: np.ndarray) -> float:
    """Total O3 column (Dobson Units) from molec cm^-3 on geometric km grid."""
    z_cm = np.asarray(z_km, dtype=np.float64) * 1e5
    n = np.asarray(n_o3_cm3, dtype=np.float64)
    col_molec_cm2 = float(np.trapz(n, z_cm))
    return col_molec_cm2 / DU_PER_COL


def _interp_clamped(xq: object, x: object, y: object) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    xq = np.asarray(xq, dtype=np.float64)
    return np.interp(xq, x, y, left=y[0], right=y[-1])


def extend_temperature_profile(z_km, t_k, top_km: float = TEMP_PROFILE_TOP_KM):
    z_km = np.asarray(z_km, dtype=np.float64)
    t_k = np.asarray(t_k, dtype=np.float64)
    zmax = float(z_km[-1])
    if zmax >= top_km - 1e-6:
        return z_km, t_k
    z_extra = np.arange(np.nextafter(zmax, np.inf), top_km + 1.0e-6, 2.0)
    if z_extra.size == 0:
        z_extra = np.array([top_km], dtype=np.float64)
    elif z_extra[-1] < top_km:
        z_extra = np.append(z_extra, top_km)
    t_extra = np.full_like(z_extra, t_k[-1])
    return np.concatenate([z_km, z_extra]), np.concatenate([t_k, t_extra])


def write_ussa_temp(out_path: Path, z_km, t_k, source: str) -> None:
    # TUV/vptmp.f uses fixed-size buffers zd(kdata=150), td(kdata=150) and
    # reads until end-of-file with no bounds checking. Keep the number of
    # (z, T) pairs <= 150 to avoid memory corruption / segfaults.
    MAX_VPTMP_POINTS = 150

    lines = [
        f"# Temperature profile from {source} (via run_TUV.py)",
        "# Column 1: geometric altitude (km), column 2: temperature (K)",
        "#",
    ]

    z_km = np.asarray(z_km, dtype=np.float64).ravel()
    t_k = np.asarray(t_k, dtype=np.float64).ravel()

    if z_km.size != t_k.size:
        raise ValueError(f"write_ussa_temp expects z and t same length; got {z_km.size} vs {t_k.size}")

    if z_km.size > MAX_VPTMP_POINTS:
        z_target = np.linspace(float(z_km[0]), float(z_km[-1]), MAX_VPTMP_POINTS)
        t_k = np.interp(z_target, z_km, t_k)
        z_km = z_target

    for z, t in zip(z_km.astype(float), t_k.astype(float)):
        lines.append(f"{z:8.3f} {t:10.4f}")

    out_path.write_text("\n".join(lines) + "\n", encoding="ascii")


def write_ussa_ozone(out_path: Path, z_grid_km, n_o3_cm3, source: str) -> None:
    lines = [
        f"# Ozone profile from {source} (via run_TUV.py)",
        "# Column 1: geometric altitude (km), column 2: n(O3) in molec cm-3",
        "# values at 0 and 1 km: interpolated from model; 2-74 km knots match USSA grid",
        "# (total column scaling in TUV still uses o3col if set > 0 in the input deck)",
        "#",
        "#",
        "#",
    ]
    if len(z_grid_km) != 39 or len(n_o3_cm3) != 39:
        raise ValueError("Ozone template requires exactly 39 levels.")
    for z, n in zip(z_grid_km, n_o3_cm3):
        lines.append(f"{z:4.0f} {n:8.2E}")
    out_path.write_text("\n".join(lines) + "\n", encoding="ascii")


def install_profile_atm_files(
    root: Path,
    z_raw: np.ndarray,
    t_raw: np.ndarray,
    n_o3_raw: np.ndarray,
    source_label: str,
) -> float:
    """Write ussa.temp / ussa.ozone; return integrated O3 column (DU)."""
    z_t, t_t = extend_temperature_profile(z_raw, t_raw)
    write_ussa_temp(root / "DATAE1" / "ATM" / "ussa.temp", z_t, t_t, source_label)

    z_tgt = np.array(_OZONE_FILE_Z_KM, dtype=np.float64)
    n_o3_on_grid = _interp_clamped(z_tgt, z_raw, n_o3_raw)
    write_ussa_ozone(root / "DATAE1" / "ATM" / "ussa.ozone", z_tgt, n_o3_on_grid, source_label)

    o3col = compute_o3_column_du(z_raw, n_o3_raw)
    print(
        f"Installed {source_label} profiles -> DATAE1/ATM/ussa.* "
        f"(O3 column = {o3col:.3f} DU)",
        file=sys.stderr,
    )
    return o3col


def install_photochem_atm_files(root: Path, photo_path: Path) -> float:
    z_raw, t_raw, _p, _vmr, n_o3_raw = read_photochem_vertical_profile(photo_path)
    return install_profile_atm_files(
        root, z_raw, t_raw, n_o3_raw, f"Photochem ({photo_path.name})"
    )


def install_atmos_atm_files(root: Path, atmos_path: Path) -> float:
    z_raw, t_raw, _p, _vmr, n_o3_raw = read_atmos_vertical_profile(atmos_path)
    return install_profile_atm_files(
        root, z_raw, t_raw, n_o3_raw, f"Atmos ({atmos_path.name})"
    )


def install_kasting_atm_files(root: Path, kasting_path: Path) -> float:
    z_raw, t_raw, _p, _vmr, n_o3_raw = read_kasting_vertical_profile(kasting_path)
    return install_profile_atm_files(
        root, z_raw, t_raw, n_o3_raw, f"Kasting ({kasting_path.name})"
    )


def install_cam_atm_files(root: Path, cam_path: Path) -> float:
    z_raw, t_raw, p_raw, o3_raw = read_cam_vertical_profile(cam_path)
    n_o3_raw = vmr_to_o3_number_density(o3_raw, p_raw, t_raw)
    return install_profile_atm_files(
        root, z_raw, t_raw, n_o3_raw, f"CAM/WACCM ({cam_path.name})"
    )


def install_vulcan_atm_files(root: Path, vul_path: Path) -> float:
    z_raw, t_raw, _p, _vmr, n_o3_raw = read_vulcan_vertical_profile_full(vul_path)
    return install_profile_atm_files(
        root, z_raw, t_raw, n_o3_raw, f"VULCAN ({vul_path.name})"
    )


def archive_o3_profile(
    case_key: str,
    model: str,
    label: str,
    source_file: str,
    z_km: np.ndarray,
    t_k: np.ndarray,
    p_bar: np.ndarray,
    o3_vmr: np.ndarray,
    n_o3_cm3: np.ndarray,
    o3col_du: float,
) -> Path:
    """Write per-model O3 profile CSV and return path."""
    OUTPUT_COMPARE_DIR.mkdir(parents=True, exist_ok=True)
    stem = _friendly_stem(case_key, model)
    out_path = OUTPUT_COMPARE_DIR / f"{stem}_o3_profile.csv"
    fieldnames = [
        "case", "label", "model", "source_file", "o3col_du",
        "z_km", "t_k", "p_bar", "o3_vmr", "n_o3_cm3",
    ]
    with out_path.open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i in range(len(z_km)):
            writer.writerow({
                "case": case_key,
                "label": label,
                "model": model,
                "source_file": source_file,
                "o3col_du": round(o3col_du, 4),
                "z_km": z_km[i],
                "t_k": t_k[i],
                "p_bar": p_bar[i],
                "o3_vmr": o3_vmr[i],
                "n_o3_cm3": n_o3_cm3[i],
            })
    print(f"Wrote O3 profile -> {out_path}", file=sys.stderr)
    return out_path


def _backup_atm_files(root: Path) -> None:
    atm_dir = root / "DATAE1" / "ATM"
    for name in ("ussa.temp", "ussa.ozone"):
        src = atm_dir / name
        bak = atm_dir / f"{name}.orig_backup"
        if src.is_file() and not bak.is_file():
            shutil.copy2(src, bak)
            print(f"Backed up {src.name} -> {bak.name}", file=sys.stderr)


def write_batch_usrinp(
    root: Path,
    outfil: str,
    o3col: float,
    sza_start: float = SZA_SWEEP_START,
    sza_stop: float = SZA_SWEEP_STOP,
    nt: int | None = None,
    wstart_nm: float | None = None,
    wstop_nm: float | None = None,
    nwint: int | None = None,
) -> Path:
    """Write INPUTS/usrinp for batch TUV (SZA sweep, UV damage actions enabled)."""
    if len(outfil) > 6:
        raise ValueError(f"outfil must be <= 6 characters for TUV, got {outfil!r}")

    if nt is None:
        nt = sza_sweep_count()

    ws = CURRENT_WSTART_NM if wstart_nm is None else float(wstart_nm)
    we = CURRENT_WSTOP_NM if wstop_nm is None else float(wstop_nm)
    if nwint is None:
        nwint = wavelength_nwint(ws, we)

    template_path = root / "INPUTS" / "Z1_pc"
    lines = template_path.read_text(encoding="ascii").splitlines()

    # Patch only the fields that change; keep Fortran column alignment from Z1_pc.
    lines[2] = (
        f"inpfil =      usrinp   outfil =      {outfil:<6}   nstr =            -2"
    )
    # Fortran FORMAT 125: wstart = F11.3, wstop = 1X+F11.3, nwint = 1X+I11
    lines[6] = (
        f"wstart = {ws:11.3f}"
        f"   wstop =  {we:11.3f}"
        f"   nwint =  {nwint:11d}"
    )
    lines[7] = (
        "tstart = "
        + f"{sza_start:11.3f}"
        + "   tstop = "
        + f"{sza_stop:11.3f}"
        + "   nt = "
        + f"{nt:11d}"
    )
    lines[9] = (
        "o3col = "
        + f"{o3col:11.3f}"
        + "   so2col = "
        + f"{0.0:11.3f}"
        + "   no2col = "
        + f"{0.0:11.3f}"
    )

    out_path = root / "INPUTS" / "usrinp"
    out_path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return out_path


def tuv_raw_output_path(outfil: str) -> Path:
    """TUV writes ../{outfil}.txt relative to the V5.4 directory."""
    return ROOT.parent / f"{outfil}.txt"


def run_tuv(root: Path, outfil: str | None = None) -> None:
    """Compile (if needed) and run TUV; fail loudly if the Fortran binary STOPs."""
    if outfil:
        stale = tuv_raw_output_path(outfil)
        if stale.is_file():
            stale.unlink()
    subprocess.run(["make"], cwd=root, check=True)
    result = subprocess.run(
        ["./tuv"],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    combined = (result.stdout or "") + (result.stderr or "")
    if combined.strip():
        print(combined, end="" if combined.endswith("\n") else "\n", file=sys.stderr)
    if result.returncode != 0 or "STOP" in combined:
        raise RuntimeError(
            f"TUV failed (exit={result.returncode}). "
            "Check wavelength grid (nwint=-156 required below 205.8 nm)."
        )
    if outfil and not tuv_raw_output_path(outfil).is_file():
        raise RuntimeError(f"TUV finished but output missing: {tuv_raw_output_path(outfil)}")


def _sza_column_name(sza: float) -> str:
    if abs(sza - round(sza)) < 0.01:
        return f"sza_{int(round(sza))}"
    return f"sza_{sza:.1f}".replace(".", "p")


def read_tuv_output(filepath: Path) -> tuple[list[dict], list[dict]]:
    """
    Parse TUV text output (compatible with early_earth_lib.read_usrout).

    Returns (spectral_rows, dose_rows) as lists of dicts.
    """
    text = filepath.read_text(encoding="ascii", errors="replace")

    spectral_header = re.search(r"^\s*wc, nm\s+(.*)$", text, re.M)
    spectral_match = re.search(r"^\s*wc, nm.*?\n((?:\s*\d.*\n)+)-", text, re.M)
    if not spectral_match or not spectral_header:
        raise ValueError(f"Could not find spectral table in {filepath}")

    sza_vals = [float(x) for x in spectral_header.group(1).split()]
    spec_col_names = [_sza_column_name(v) for v in sza_vals]

    spec_lines = [ln for ln in spectral_match.group(1).splitlines() if ln.strip()]
    spectral_rows: list[dict] = []
    for line in spec_lines:
        parts = [float(x) for x in line.split()]
        row = {"wavelength_nm": parts[0]}
        for name, val in zip(spec_col_names, parts[1:]):
            row[name] = val
        spectral_rows.append(row)

    dose_match = re.search(r"sza, deg\..*?\n(.*?)\n-+", text, re.S)
    if not dose_match:
        raise ValueError(f"Could not find dose table in {filepath}")

    dose_lines = [ln for ln in dose_match.group(1).splitlines() if ln.strip()]
    dose_rows: list[dict] = []
    for line in dose_lines:
        parts = [float(x) for x in line.split()]
        row = {"sza_deg": parts[0]}
        for name, val in zip(DOSE_COLUMNS, parts[1:]):
            row[name] = val
        dose_rows.append(row)

    return spectral_rows, dose_rows


def _friendly_stem(case_key: str, model: str) -> str:
    atm_label = str(VULCAN_ATM_SZA_DEG).replace(".", "p")
    sweep = f"sza{int(SZA_SWEEP_START)}-{int(SZA_SWEEP_STOP)}"
    return f"{case_key}_{model}_atm{atm_label}_{sweep}_{wavelength_tag()}"


def archive_tuv_output(
    outfil: str,
    case_key: str,
    model: str,
    meta: dict,
) -> Path:
    """Copy raw TUV output and write parsed CSV + metadata JSON."""
    raw_src = tuv_raw_output_path(outfil)
    if not raw_src.is_file():
        raise FileNotFoundError(f"TUV output not found: {raw_src}")

    OUTPUT_COMPARE_DIR.mkdir(parents=True, exist_ok=True)
    stem = _friendly_stem(case_key, model)

    txt_dst = OUTPUT_COMPARE_DIR / f"{stem}.txt"
    shutil.copy2(raw_src, txt_dst)

    spectral, dose = read_tuv_output(raw_src)

    spec_csv = OUTPUT_COMPARE_DIR / f"{stem}_spectral.csv"
    if spectral:
        with spec_csv.open("w", newline="", encoding="ascii") as f:
            writer = csv.DictWriter(f, fieldnames=list(spectral[0].keys()))
            writer.writeheader()
            writer.writerows(spectral)

    dose_csv = OUTPUT_COMPARE_DIR / f"{stem}_dose.csv"
    if dose:
        with dose_csv.open("w", newline="", encoding="ascii") as f:
            writer = csv.DictWriter(f, fieldnames=list(dose[0].keys()))
            writer.writeheader()
            writer.writerows(dose)

    meta_path = OUTPUT_COMPARE_DIR / f"{stem}_meta.json"
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="ascii")

    print(f"Archived -> {txt_dst}", file=sys.stderr)
    return txt_dst


def _parse_output_stem(stem: str) -> tuple[str, str]:
    """Extract (case_key, model) from archive filename stem."""
    atm_label = str(VULCAN_ATM_SZA_DEG).replace(".", "p")
    sweep = f"sza{int(SZA_SWEEP_START)}-{int(SZA_SWEEP_STOP)}"
    # Prefer wavelength-tagged stems; also accept legacy untagged archives.
    for suffix in (
        f"_atm{atm_label}_{sweep}_{wavelength_tag()}",
        f"_atm{atm_label}_{sweep}",
    ):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    case_key, model = stem.rsplit("_", 1)
    return case_key, model


def run_single_comparison(
    root: Path,
    case_key: str,
    model: str,
) -> list[dict]:
    """Run one atmosphere model through TUV; return one summary dict per SZA."""
    if case_key not in COMPARISON_CASES:
        raise KeyError(f"Unknown comparison case {case_key!r}")

    cfg = COMPARISON_CASES[case_key]
    _backup_atm_files(root)

    waccm_lat_deg: float | None = None
    waccm_ilat: int | None = None
    waccm_horizontal_mode: str | None = None
    if model == "waccm":
        cam_path = DATASET_CASES[cfg["waccm_case"]].expanduser().resolve()
        if not cam_path.is_file():
            raise FileNotFoundError(f"WACCM file not found: {cam_path}")
        waccm_horizontal_mode = "gaussian_lat_mean"
        z_km, t_k, p_bar, o3_vmr, n_o3 = read_waccm_vertical_profile(
            cam_path, horizontal_mode=waccm_horizontal_mode
        )
        o3col = install_profile_atm_files(
            root,
            z_km,
            t_k,
            n_o3,
            f"CAM/WACCM LWAV ({cam_path.name})",
        )
        outfil = cfg["outfil_waccm"]
        source_path = str(cam_path)
        source_file = cam_path.name
    elif model == "waccm_eq":
        cam_path = DATASET_CASES[cfg["waccm_case"]].expanduser().resolve()
        if not cam_path.is_file():
            raise FileNotFoundError(f"WACCM file not found: {cam_path}")
        waccm_ilat = WACCM_EQ_ILAT
        z_km, t_k, p_bar, o3_vmr, n_o3 = read_waccm_vertical_profile(
            cam_path, ilat=waccm_ilat
        )
        with xr.open_dataset(cam_path, decode_times=False) as ds:
            waccm_lat_deg = float(ds[DIM_LAT].isel({DIM_LAT: waccm_ilat}).values)
        o3col = install_profile_atm_files(
            root,
            z_km,
            t_k,
            n_o3,
            f"CAM/WACCM equator ILAT={waccm_ilat} ({cam_path.name})",
        )
        outfil = cfg["outfil_waccm_eq"]
        source_path = str(cam_path)
        source_file = cam_path.name
    elif model == "vulcan":
        vul_path = (VIH_OUTPUT_PATH / cfg["vulcan_file"]).resolve()
        if not vul_path.is_file():
            raise FileNotFoundError(f"VULCAN file not found: {vul_path}")
        z_km, t_k, p_bar, o3_vmr, n_o3 = read_vulcan_vertical_profile_full(vul_path)
        o3col = install_profile_atm_files(
            root, z_km, t_k, n_o3, f"VULCAN ({vul_path.name})"
        )
        outfil = cfg["outfil_vulcan"]
        source_path = str(vul_path)
        source_file = vul_path.name
    elif model == "photochem":
        photo_path = (PHOTOCHEM_BASE / cfg["photochem_file"]).resolve()
        if not photo_path.is_file():
            raise FileNotFoundError(f"Photochem file not found: {photo_path}")
        z_km, t_k, p_bar, o3_vmr, n_o3 = read_photochem_vertical_profile(photo_path)
        o3col = install_profile_atm_files(
            root, z_km, t_k, n_o3, f"Photochem ({photo_path.name})"
        )
        outfil = cfg["outfil_photochem"]
        source_path = str(photo_path)
        source_file = photo_path.name
    elif model == "atmos":
        atmos_path = (ATMOS_BASE / cfg["atmos_ptz"]).resolve()
        if not atmos_path.is_file():
            raise FileNotFoundError(f"Atmos file not found: {atmos_path}")
        z_km, t_k, p_bar, o3_vmr, n_o3 = read_atmos_vertical_profile(atmos_path)
        o3col = install_profile_atm_files(
            root, z_km, t_k, n_o3, f"Atmos ({atmos_path.name})"
        )
        outfil = cfg["outfil_atmos"]
        source_path = str(atmos_path)
        source_file = atmos_path.name
    elif model == "kasting":
        kasting_path = (KASTING_BASE / cfg["kasting_plot"]).resolve()
        if not kasting_path.is_file():
            raise FileNotFoundError(f"Kasting file not found: {kasting_path}")
        z_km, t_k, p_bar, o3_vmr, n_o3 = read_kasting_vertical_profile(kasting_path)
        o3col = install_profile_atm_files(
            root, z_km, t_k, n_o3, f"Kasting ({kasting_path.name})"
        )
        outfil = cfg["outfil_kasting"]
        source_path = str(kasting_path)
        source_file = kasting_path.name
    else:
        raise ValueError(
            f"model must be one of {ALL_COMPARISON_MODELS}, got {model!r}"
        )

    archive_o3_profile(
        case_key, model, cfg["label"], source_file,
        z_km, t_k, p_bar, o3_vmr, n_o3, o3col,
    )

    write_batch_usrinp(root, outfil=outfil, o3col=o3col)
    run_tuv(root, outfil=outfil)

    meta = {
        "case": case_key,
        "label": cfg["label"],
        "model": model,
        "atmosphere_equilibrium_sza_deg": VULCAN_ATM_SZA_DEG,
        "tuv_sza_start_deg": SZA_SWEEP_START,
        "tuv_sza_stop_deg": SZA_SWEEP_STOP,
        "tuv_sza_step_deg": SZA_SWEEP_STEP,
        "tuv_sza_count": sza_sweep_count(),
        "wstart_nm": CURRENT_WSTART_NM,
        "wstop_nm": CURRENT_WSTOP_NM,
        "nwint": wavelength_nwint(),
        "o3col_du": round(o3col, 4),
        "outfil": outfil,
        "source_file": source_file,
        "source_path": source_path,
    }
    if waccm_lat_deg is not None:
        meta["waccm_latitude_deg"] = waccm_lat_deg
    if waccm_ilat is not None:
        meta["waccm_ilat"] = waccm_ilat
    if waccm_horizontal_mode is not None:
        meta["waccm_horizontal_mode"] = waccm_horizontal_mode
    archive_tuv_output(outfil, case_key, model, meta)

    _, dose = read_tuv_output(tuv_raw_output_path(outfil))
    if not dose:
        raise RuntimeError(f"No dose rows parsed for {case_key}/{model}")

    summaries: list[dict] = []
    for drow in dose:
        summaries.append({
            "case": case_key,
            "label": cfg["label"],
            "model": model,
            "sza_deg": drow["sza_deg"],
            "o3col_du": round(o3col, 4),
            "source_file": source_file,
            **{k: drow[k] for k in DOSE_COLUMNS if k in drow},
        })
    return summaries


def rebuild_comparison_from_archives() -> None:
    """Rebuild combined comparison CSVs from archived runs for the active wavelength window."""
    OUTPUT_COMPARE_DIR.mkdir(parents=True, exist_ok=True)

    summary_fieldnames = [
        "case", "label", "model", "sza_deg", "o3col_du", "source_file",
        *DOSE_COLUMNS,
    ]
    summary_rows: list[dict] = []
    spectral_rows: list[dict] = []
    spectral_fieldnames = ["case", "model"]
    tag = wavelength_tag()

    for meta_path in sorted(OUTPUT_COMPARE_DIR.glob("*_meta.json")):
        meta = json.loads(meta_path.read_text(encoding="ascii"))
        if not meta_matches_wavelength(meta):
            continue
        stem = meta_path.name[: -len("_meta.json")]
        dose_path = OUTPUT_COMPARE_DIR / f"{stem}_dose.csv"
        if dose_path.is_file():
            with dose_path.open(encoding="ascii") as f:
                for drow in csv.DictReader(f):
                    summary_rows.append({
                        "case": meta["case"],
                        "label": meta["label"],
                        "model": meta["model"],
                        "sza_deg": drow["sza_deg"],
                        "o3col_du": meta["o3col_du"],
                        "source_file": meta["source_file"],
                        **{k: drow.get(k, "") for k in DOSE_COLUMNS},
                    })

        spec_path = OUTPUT_COMPARE_DIR / f"{stem}_spectral.csv"
        if spec_path.is_file():
            with spec_path.open(encoding="ascii") as f:
                reader = csv.DictReader(f)
                for name in reader.fieldnames or []:
                    if name not in spectral_fieldnames:
                        spectral_fieldnames.append(name)
                for srow in reader:
                    spectral_rows.append({
                        "case": meta["case"],
                        "model": meta["model"],
                        **srow,
                    })

    summary_rows.sort(
        key=lambda row: (row["case"], row["model"], float(row["sza_deg"]))
    )
    tagged_summary = comparison_summary_path()
    with tagged_summary.open("w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(
            f, fieldnames=summary_fieldnames, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Wrote {tagged_summary} ({len(summary_rows)} rows)", file=sys.stderr)

    # Keep legacy untagged name in sync for the default 280–420 window only.
    if (
        abs(CURRENT_WSTART_NM - DEFAULT_WSTART_NM) < 0.05
        and abs(CURRENT_WSTOP_NM - DEFAULT_WSTOP_NM) < 0.05
    ):
        legacy_summary = OUTPUT_COMPARE_DIR / "comparison_summary.csv"
        shutil.copy2(tagged_summary, legacy_summary)
        print(f"Wrote {legacy_summary} (default-window copy)", file=sys.stderr)

    if spectral_rows:
        tagged_spec = OUTPUT_COMPARE_DIR / f"comparison_spectral_{tag}.csv"
        with tagged_spec.open("w", newline="", encoding="ascii") as out_f:
            writer = csv.DictWriter(
                out_f, fieldnames=spectral_fieldnames, extrasaction="ignore"
            )
            writer.writeheader()
            writer.writerows(spectral_rows)
        print(f"Wrote {tagged_spec} ({len(spectral_rows)} rows)", file=sys.stderr)
        if (
            abs(CURRENT_WSTART_NM - DEFAULT_WSTART_NM) < 0.05
            and abs(CURRENT_WSTOP_NM - DEFAULT_WSTOP_NM) < 0.05
        ):
            legacy_spec = OUTPUT_COMPARE_DIR / "comparison_spectral.csv"
            shutil.copy2(tagged_spec, legacy_spec)

    prof_fieldnames = [
        "case", "label", "model", "source_file", "o3col_du",
        "z_km", "t_k", "p_bar", "o3_vmr", "n_o3_cm3",
    ]
    out_path = OUTPUT_COMPARE_DIR / "ozone_profiles_comparison.csv"
    prof_count = 0
    with out_path.open("w", newline="", encoding="ascii") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=prof_fieldnames)
        writer.writeheader()
        for prof_path in sorted(OUTPUT_COMPARE_DIR.glob("*_o3_profile.csv")):
            stem = prof_path.name[: -len("_o3_profile.csv")]
            is_tagged = stem.endswith(f"_{tag}")
            is_legacy_default = (
                abs(CURRENT_WSTART_NM - DEFAULT_WSTART_NM) < 0.05
                and abs(CURRENT_WSTOP_NM - DEFAULT_WSTOP_NM) < 0.05
                and re.search(r"_w\d+-\d+$", stem) is None
            )
            if not (is_tagged or is_legacy_default):
                continue
            with prof_path.open(encoding="ascii") as in_f:
                for row in csv.DictReader(in_f):
                    writer.writerow({k: row[k] for k in prof_fieldnames if k in row})
                    prof_count += 1
    if prof_count:
        print(f"Wrote {out_path} ({prof_count} rows)", file=sys.stderr)


def write_comparison_tables(summaries: list[dict], spectral_stems: list[str]) -> None:
    del summaries, spectral_stems
    rebuild_comparison_from_archives()


def write_ozone_profiles_combined(case_keys: list[str], models: list[str]) -> None:
    """Merge per-run *_o3_profile.csv into one comparison table."""
    OUTPUT_COMPARE_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUTPUT_COMPARE_DIR / "ozone_profiles_comparison.csv"
    fieldnames = [
        "case", "label", "model", "source_file", "o3col_du",
        "z_km", "t_k", "p_bar", "o3_vmr", "n_o3_cm3",
    ]
    wrote_header = False
    with out_path.open("w", newline="", encoding="ascii") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames)
        for case_key in case_keys:
            for model in models:
                stem = _friendly_stem(case_key, model)
                prof_path = OUTPUT_COMPARE_DIR / f"{stem}_o3_profile.csv"
                if not prof_path.is_file():
                    continue
                with prof_path.open(encoding="ascii") as in_f:
                    reader = csv.DictReader(in_f)
                    if not wrote_header:
                        writer.writeheader()
                        wrote_header = True
                    for row in reader:
                        writer.writerow({k: row[k] for k in fieldnames if k in row})
    if wrote_header:
        print(f"Wrote {out_path}", file=sys.stderr)


def run_comparison(
    root: Path,
    case_keys: list[str],
    models: list[str],
) -> None:
    summaries: list[dict] = []
    spectral_stems: list[str] = []

    for case_key in case_keys:
        for model in models:
            print(f"\n=== {case_key} / {model} ===", file=sys.stderr)
            summaries.extend(run_single_comparison(root, case_key, model))
            spectral_stems.append(_friendly_stem(case_key, model))

    write_comparison_tables(summaries, spectral_stems)


def _build_arg_parser() -> argparse.ArgumentParser:
    waccm_keys = ", ".join(sorted(DATASET_CASES))
    cmp_keys = ", ".join(sorted(COMPARISON_CASES))
    p = argparse.ArgumentParser(
        description="Install WACCM/VULCAN atmospheres, run TUV, archive comparison outputs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            f"WACCM --case keys: {waccm_keys}\n"
            f"Comparison --cases keys: {cmp_keys}\n"
            f"VULCAN files under: {VIH_OUTPUT_PATH}\n"
            f"Outputs under: {OUTPUT_COMPARE_DIR}"
        ),
    )
    p.add_argument("--case", metavar="NAME", help="WACCM case from DATASET_CASES.")
    p.add_argument("--cases", metavar="LIST", help="Comma-separated comparison cases (for --compare).")
    p.add_argument("--cam-file", type=Path, metavar="PATH", help="Explicit WACCM NetCDF path.")
    p.add_argument("--vulcan-file", type=Path, metavar="PATH", help="Explicit VULCAN .vul path.")
    p.add_argument(
        "--model",
        choices=(
            "waccm", "waccm_eq", "vulcan", "photochem", "atmos", "kasting",
            "both", "all",
        ),
        default="all",
        help=(
            "Atmosphere source for --compare (default: all). "
            "'both' = waccm+vulcan only."
        ),
    )
    p.add_argument(
        "--compare",
        action="store_true",
        help=(
            "Run all requested models for each PAL case. "
            "VULCAN/Photochem/Atmos/Kasting use 48.2-deg-equilibrium outputs; "
            f"TUV SZA sweep {SZA_SWEEP_START:g}-{SZA_SWEEP_STOP:g} "
            f"every {SZA_SWEEP_STEP:g} deg. Writes ozone_profiles_comparison.csv."
        ),
    )
    p.add_argument("--list-cases", action="store_true", help="List WACCM and comparison cases.")
    p.add_argument("--cam-only", action="store_true", help="Only write ussa.* from WACCM; skip TUV.")
    p.add_argument("--no-cam", action="store_true", help="Skip WACCM install on a plain run.")
    p.add_argument(
        "--wstart",
        type=float,
        default=DEFAULT_WSTART_NM,
        metavar="NM",
        help=(
            f"TUV spectral start wavelength in nm (default: {DEFAULT_WSTART_NM:g}; "
            "minimum 100). Use 100 to include UVC in DNA-damage weighting."
        ),
    )
    p.add_argument(
        "--wstop",
        type=float,
        default=DEFAULT_WSTOP_NM,
        metavar="NM",
        help=f"TUV spectral stop wavelength in nm (default: {DEFAULT_WSTOP_NM:g}).",
    )
    return p


def resolve_cam_path(args: argparse.Namespace) -> Path:
    if args.cam_file is not None:
        p = Path(args.cam_file).expanduser().resolve()
        if not p.is_file():
            raise SystemExit(f"--cam-file not found: {p}")
        return p

    case = (args.case or os.environ.get("TUV_CAM_CASE") or "").strip()
    if case:
        if case not in DATASET_CASES:
            keys = ", ".join(sorted(DATASET_CASES))
            raise SystemExit(f"Unknown case {case!r}. Registered keys: {keys}")
        p = DATASET_CASES[case].expanduser().resolve()
        if not p.is_file():
            raise SystemExit(f"Case {case!r} file not found: {p}")
        return p

    p = DATASET_CASES[DEFAULT_CAM_CASE].expanduser().resolve()
    if not p.is_file():
        raise SystemExit(
            f"Default case {DEFAULT_CAM_CASE!r} file not found: {p}. "
            "Set --case, --cam-file, or fix DATASET_CASES / CESM_DATA."
        )
    return p


def _models_from_arg(model: str) -> list[str]:
    if model == "all":
        return list(ALL_COMPARISON_MODELS)
    if model == "both":
        return ["waccm", "vulcan"]
    return [model]


def _comparison_case_list(args: argparse.Namespace) -> list[str]:
    if args.cases:
        keys = [k.strip() for k in args.cases.split(",") if k.strip()]
    else:
        keys = list(COMPARISON_CASES)
    unknown = [k for k in keys if k not in COMPARISON_CASES]
    if unknown:
        valid = ", ".join(sorted(COMPARISON_CASES))
        raise SystemExit(f"Unknown comparison case(s) {unknown}. Valid: {valid}")
    return keys


def main(argv: list[str] | None = None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    args = _build_arg_parser().parse_args(argv)
    root = ROOT

    try:
        set_wavelength_range(args.wstart, args.wstop)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    print(
        f"TUV spectral window: {CURRENT_WSTART_NM:g}–{CURRENT_WSTOP_NM:g} nm "
        f"(nwint={wavelength_nwint()}, tag={wavelength_tag()})",
        file=sys.stderr,
    )

    if args.list_cases:
        print("WACCM cases:")
        for name in sorted(DATASET_CASES):
            print(f"  {name:12}  {DATASET_CASES[name]}")
        print("\nComparison cases (48.2 deg SZA atmospheres):")
        for name in sorted(COMPARISON_CASES):
            cfg = COMPARISON_CASES[name]
            print(f"  {name:12}  {cfg['label']}")
            print(f"               WACCM:     {DATASET_CASES[cfg['waccm_case']].name}")
            print(f"               VULCAN:    {cfg['vulcan_file']}")
            print(f"               Photochem: {cfg['photochem_file']}")
            print(f"               Atmos:     {cfg['atmos_ptz']}")
            print(f"               Kasting:   {cfg['kasting_plot']}")
        return

    if args.compare:
        run_comparison(root, _comparison_case_list(args), _models_from_arg(args.model))
        return

    if args.case and args.case in COMPARISON_CASES:
        models = _models_from_arg(args.model)
        for model in models:
            run_single_comparison(root, args.case, model)
        return

    if args.model == "vulcan" and (args.vulcan_file or args.case):
        _backup_atm_files(root)
        if args.vulcan_file:
            vul_path = Path(args.vulcan_file).expanduser().resolve()
        elif args.case and args.case in COMPARISON_CASES:
            vul_path = (VIH_OUTPUT_PATH / COMPARISON_CASES[args.case]["vulcan_file"]).resolve()
        else:
            raise SystemExit("VULCAN run requires --vulcan-file or --case from COMPARISON_CASES.")
        if not vul_path.is_file():
            raise SystemExit(f"VULCAN file not found: {vul_path}")
        o3col = install_vulcan_atm_files(root, vul_path)
        outfil = "VulRun"
        if args.case in COMPARISON_CASES:
            outfil = COMPARISON_CASES[args.case]["outfil_vulcan"]
        write_batch_usrinp(root, outfil=outfil, o3col=o3col)
        if args.cam_only:
            return
        run_tuv(root)
        return

    want_cam = not args.no_cam and (
        USE_CAM_ATM
        or args.cam_only
        or bool(args.case or args.cam_file or os.environ.get("TUV_CAM_CASE"))
    )

    if want_cam:
        cam_path = resolve_cam_path(args)
        _backup_atm_files(root)
        install_cam_atm_files(root, cam_path)

    if args.cam_only:
        return

    run_tuv(root)


if __name__ == "__main__":
    main()
