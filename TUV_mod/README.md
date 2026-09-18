# TUV 5.4 — Early Earth UV comparison workflow

This directory contains **NCAR TUV v5.4** (Tropospheric Ultraviolet–Visible radiation model), extended with a Python driver and input decks so that **surface UV and biological dose rates** can be computed consistently from:

- **WACCM6 / CAM** three-dimensional atmospheres (CESM NetCDF output)
- **VULCAN** photochemical equilibrium columns (`.vul` pickle files)
- **Photochem** 1D model outputs
- **Atmos** 1D model outputs
- **Kasting** 1D model outputs

The scientific goal is to compare how different atmospheric-chemistry models (and different oxygen levels / PAL cases) affect **spectral surface irradiance** and **biologically weighted UV** (DNA damage, UV index, plant damage, etc.) over a **solar zenith angle (SZA) sweep**, using the same radiative-transfer solver for all inputs.

Upstream TUV reference: Madronich et al., NCAR TUV 5.4 (November 2018). See also `README.txt` and `VERSION_NOTES.txt` for the stock model documentation.

---

## Table of contents

1. [Conceptual workflow](#conceptual-workflow)
2. [Atmosphere sources and how they connect](#atmosphere-sources-and-how-they-connect)
3. [Quick start](#quick-start)
4. [Using the model](#using-the-model)
5. [Interpreting output](#interpreting-output)
6. [Spectral vs weighted (biological) UV](#spectral-vs-weighted-biological-uv)
7. [Output files and directory layout](#output-files-and-directory-layout)
8. [Changes compared with stock TUV 5.4](#changes-compared-with-stock-tuv-54)
9. [Where to find modifications in the code](#where-to-find-modifications-in-the-code)
10. [Dependencies and external data paths](#dependencies-and-external-data-paths)
11. [Troubleshooting](#troubleshooting)

The section [Atmosphere sources and how they connect](#atmosphere-sources-and-how-they-connect) includes the Early_Earth.py ↔ TUV case crosswalk (`Pre_h0`, `Ten_h0`, `One_h0`, `Zero1_h0`) and the full 3D + 1D PAL matrix.

---

## Conceptual workflow

Each comparison run follows the same pipeline:

```
  Chemistry model output                TUV radiative transfer           Archived products
 ┌─────────────────────────┐          ┌──────────────────────┐         ┌─────────────────────────┐
 │ WACCM6 NetCDF           │          │                      │         │ *_spectral.csv          │
 │ VULCAN .vul             │  ──►     │  Install T(z),       │  ──►    │ *_dose.csv              │
 │ Photochem / Atmos /     │  write   │  n(O₃)(z) into      │  parse  │ *_o3_profile.csv        │
 │ Kasting 1D files        │  ussa.*  │  ussa.temp /         │         │ *_meta.json             │
 └─────────────────────────┘          │  ussa.ozone          │         │ comparison_summary*.csv │
                                      │                      │         └─────────────────────────┘
                                      │  SZA sweep 0–85°     │
                                      │  (lzenit = T)        │
                                      └──────────────────────┘
```

Important design choice:

- **Atmosphere chemistry** is taken from each model at its **photochemical equilibrium SZA** (48.2° for most cases; 100% PAL Atmos uses 48.5°).
- **Illumination geometry** in TUV is a **separate SZA sweep** from 0° to 85° in 5° steps.
- That means TUV answers: *“Given this ozone/temperature profile from model X at equilibrium SZA ≈ 48°, what is surface UV at noon, morning, evening, etc.?”*

This decoupling is intentional: the 1D and VULCAN outputs are single-column equilibrium states, while WACCM provides a global mean or equatorial column; TUV then explores diurnal illumination on that fixed vertical structure.

---

## Atmosphere sources and how they connect

This extension feeds **the same PAL oxygen cases** through **one 3D model (WACCM6)** and **four 1D chemistry models**, then runs **identical TUV radiative transfer** on each column. The chemistry models supply O₃, T, and altitude; TUV supplies surface UV and biological dose rates.

| Model key (`--model`) | Source data | Horizontal treatment | Typical use |
|----------------------|-------------|----------------------|-------------|
| `waccm` | CESM/CAM `.nc` from `CESM_data/` | Gaussian latitude mean (LWAV over `gw` weights) | Global-mean WACCM6 column for each PAL case |
| `waccm_eq` | Same NetCDF | Single column at `ILAT = 47` (~ equator) | Equatorial WACCM6 column |
| `vulcan` | `.vul` under `VIH_cases/output/` | 1D column | VULCAN photochemical equilibrium |
| `photochem` | `Photochem/.../Earth_*_48.2.txt` | 1D column | Photochem 1D model |
| `atmos` | `Atmos/.../PTZ_mixingratios_out.dist` | 1D column | Atmos 1D model |
| `kasting` | `Kasting_1D_model/.../OUTPUT_PLOT.dat` | 1D column | Kasting 1D model |

### Link to `Early_Earth.py` WACCM variables

`Early_Earth.py` loads **two** WACCM families. TUV uses only the **newer `CESM_data` set**:

| PAL | TUV case | `DATASET_CASES` key | NetCDF under `CESM_data/` | `Early_Earth.py` variable |
|-----|----------|---------------------|---------------------------|---------------------------|
| 100% | `pal100` | `baseline` | `Earth_100pc_o2.cam.h0.0009-0012.nc` | **`Pre_h0`** |
| 10% | `pal10` | `pal10` | `Earth_10pc_o2.cam.h0.0037-0040.nc` | **`Ten_h0`** |
| 1% | `pal1` | `pal1` | `Earth_1pc_o2.cam.h0.0045-0048.nc` | **`One_h0`** |
| 0.1% | `pal0p1` | `pal0p1` | `Earth_0.1pc_o2.cam.h0.0033-0036.nc` | **`Zero1_h0`** |


### PAL comparison matrix (3D + 1D)

Four standard oxygen levels are registered in `run_TUV.py` (`COMPARISON_CASES`). For each PAL, every chemistry model contributes one atmosphere column:

| Case | Label | WACCM (`Pre_h0` family) | VULCAN | Photochem | Atmos | Kasting |
|------|-------|-------------------------|--------|-----------|-------|---------|
| `pal100` | 100% PAL | `Earth_100pc_o2.cam.h0.0009-0012.nc` (`Pre_h0`) | `Earth_1e12s_48.2SZA_WPT_1rtol.vul` | `100pc/Earth_100pc_48.2.txt` | `100pc/SZA_48.5/PTZ_…` | `100pc/SZA_48.2/OUTPUT_PLOT.dat` |
| `pal10` | 10% PAL | `Earth_10pc_o2.cam.h0.0037-0040.nc` (`Ten_h0`) | `Earth_10pc_o2_…48.2SZA….vul` | `10pc/Earth_10pc_48.2.txt` | `10pc/SZA_48.2/PTZ_…` | `10pc/SZA_48.2/OUTPUT_PLOT.dat` |
| `pal1` | 1% PAL | `Earth_1pc_o2.cam.h0.0045-0048.nc` (`One_h0`) | `Earth_1pc_o2_…48.2SZA….vul` | `1pc/Earth_1pc_48.2.txt` | `1pc/SZA_48.2/PTZ_…` | `1pc/SZA_48.2/OUTPUT_PLOT.dat` |
| `pal0p1` | 0.1% PAL | `Earth_0.1pc_o2.cam.h0.0033-0036.nc` (`Zero1_h0`) | `Earth_0.1pc_o2_…48.2SZA….vul` | `0.1pc/Earth_0.1pc_48.2.txt` | `0.1pc/SZA_48.2/PTZ_…` | `0.1pc/SZA_48.2/OUTPUT_PLOT.dat` |

Notes:

- **WACCM** is the only 3D source. `--model waccm` uses a Gaussian global mean of `Pre_h0`/`Ten_h0`/`One_h0`/`Zero1_h0`; `--model waccm_eq` uses the equatorial latitude index (`ILAT = 47`).
- **1D models** (VULCAN, Photochem, Atmos, Kasting) use their own equilibrium columns at ~48.2° SZA (Atmos 100% PAL uses SZA 48.5°).
- All six models for a given PAL are then illuminated with the **same** TUV SZA sweep (0–85°).

Run `python run_TUV.py --list-cases` to print the full registry.

### What gets transferred into TUV

For every source, `run_TUV.py`:

1. Reads **geometric altitude** `z` (km), **temperature** `T` (K), and **O₃ number density** `n(O₃)` (molec cm⁻³).
2. Computes the **total O₃ column** (Dobson Units) by integrating `n(O₃)` over height.
3. Writes:
   - `DATAE1/ATM/ussa.temp` — full temperature profile (extended to 120 km if needed, capped at 150 levels for TUV safety)
   - `DATAE1/ATM/ussa.ozone` — O₃ on the fixed 39-level USSA altitude grid (0, 1, 2, 4, …, 74 km)
4. Sets `o3col` in the TUV input deck to the integrated column (DU), so TUV uses the correct total ozone while preserving the profile shape.

Original USSA files are backed up once as `ussa.temp.orig_backup` and `ussa.ozone.orig_backup`.

---

## Quick start

### Prerequisites

- `gfortran` (see `Makefile`)
- Python 3 with `numpy`, `xarray`, `pandas`, `matplotlib` (plotting only)
- External atmosphere files at the paths configured in `run_TUV.py` (see [Dependencies](#dependencies-and-external-data-paths))

### Build TUV

```bash
cd /Users/gregcooke/V5.4
make
```

### Run the full comparison matrix (recommended)

All six atmosphere models × four PAL cases, SZA 0–85°, default spectral window 280–420 nm:

```bash
python run_TUV.py --compare
```

Include UVC (100–420 nm) for DNA-damage calculations that weight short wavelengths:

```bash
python run_TUV.py --compare --wstart 100
```

Run only 1D models for selected cases:

```bash
python run_TUV.py --compare --model photochem,atmos,kasting,vulcan --cases pal0p1,pal1
```

### Run a single model / case

```bash
python run_TUV.py --model vulcan --case pal0p1
python run_TUV.py --model waccm_eq --case pal10
python run_TUV.py --cam-only --case pal10   # install WACCM profiles only, skip TUV
```

### Plot archived results

```bash
python OUTPUT/compare/plot_tuv_comparison.py
python OUTPUT/compare/plot_tuv_comparison.py --summary OUTPUT/compare/comparison_summary_w100-420.csv --sza 0 45
```

---

## Using the model

There are two ways to run TUV in this installation: **batch mode** (default, used by `run_TUV.py`) and **interactive mode** (stock TUV behaviour).

### Batch mode (automated)

In `TUV.f`, batch mode is enabled:

```fortran
intrct = .FALSE.
IF ( .NOT. intrct) inpfil = 'usrinp'
```

`run_TUV.py` writes `INPUTS/usrinp` from the template `INPUTS/Z1_pc`, patching:

| Field | Batch value | Meaning |
|-------|-------------|---------|
| `lzenit` | `T` | `tstart` / `tstop` are **solar zenith angles** (degrees), not local time |
| `tstart`, `tstop` | `0`, `85` | SZA sweep endpoints |
| `nt` | `18` | Number of SZA steps (0, 5, …, 85) |
| `o3col` | computed | Total O₃ column (DU) from source atmosphere |
| `wstart`, `wstop`, `nwint` | CLI / defaults | Spectral window; `nwint = -156` if `wstart < 205.8` nm |
| `lrates` | `T` | Enable biologically weighted dose rates |
| `lirrad` | `T` | Enable spectral irradiance table |
| `nms` | `8` | Eight action spectra enabled (see below) |

Then TUV is executed:

```bash
./tuv
```

Raw Fortran output is written to **`../{outfil}.txt`** (one directory **above** `V5.4/`, i.e. `/Users/gregcooke/{outfil}.txt`). `run_TUV.py` copies and parses these into `OUTPUT/compare/`.

### Interactive mode (manual)

To restore stock interactive behaviour, edit `TUV.f`:

```fortran
intrct = .TRUE.
c     intrct = .FALSE.
```

Recompile (`make`), then run `./tuv` and follow the menu. Useful variables:

| Variable | Purpose |
|----------|---------|
| `lzenit` | Toggle: `T` = SZA mode, `F` = local solar time mode |
| `tstart`, `tstop`, `nt` | SZA or time grid |
| `o3col` | Total O₃ column (DU); also scales profile if > 0 |
| `wstart`, `wstop`, `nwint` | Wavelength grid |
| `nms` | Open action-spectrum menu |
| `outfil` | Output basename (max 6 characters) |

When prompted `write new value for tstart`, enter **only the number** (e.g. `48.2`), not `tstart = 48.2`.

Stock default input decks remain in `INPUTS/defin1`–`defin5`. Custom Early-Earth decks: `PI`, `One_pc`, `Ten_pc`, `Z1_pc`, `Z1_V`.

---

## Interpreting output

### Raw TUV text (`*.txt`)

Two main tables:

1. **Spectral irradiance** — header `Spectral Irradiance, W m-2 nm-1`
   - Rows: wavelength (nm)
   - Columns: SZA values
   - Units: W m⁻² nm⁻¹ (unweighted)

2. **Dose rates** — header `Dose rates, W m-2`
   - Rows: SZA (when `lzenit = T`)
   - Columns: biologically weighted irradiances

With the `Z1_pc` / batch template, the dose-rate legend maps to:

| Index | Name in TUV | Parsed CSV column |
|-------|-------------|-------------------|
| 1 | UV-B (280–315 nm) | `UVB_280_315` |
| 2 | UV-B* (280–320 nm) | `UVBstar_280_320` |
| 3 | UV-A (315–400 nm) | `UVA_315_400` |
| 4 | DNA damage (Setlow 1974) | `DNA_damage` |
| 5 | UV index | `UV_index` |
| 6 | Plant damage (Caldwell 1971) | `P_Dam_C_1971` |
| 7 | Plant damage (Flint & Caldwell 2003) | `P_Dam_FC_2003` |
| 8 | Plant damage (FC 2003, ext. 390 nm) | `P_Dam_FC_2003_ext390` |

**DNA damage** values are **weighted dose rates** in W m⁻², not a separate exported action-spectrum table. To obtain a dose over time, multiply by exposure duration (e.g. 1 h → J m⁻² = W m⁻² × 3600 s).

Example row from `pal0p1_waccm` at SZA = 0° (100–420 nm window):

```
DNA_damage = 0.659 W m-2
UV_index   = 148.8
```

### Parsed CSV outputs

| File | Contents |
|------|----------|
| `*_spectral.csv` | Wavelength × SZA irradiance matrix |
| `*_dose.csv` | SZA × biological dose rates |
| `*_o3_profile.csv` | Input O₃, T, pressure profile used for that run |
| `*_meta.json` | Run metadata (case, model, O₃ column, wavelength window, source file) |
| `comparison_summary_{wstart-wstop}.csv` | All models/cases/SZAs — primary analysis table |
| `comparison_spectral_{wstart-wstop}.csv` | Combined spectral data |
| `ozone_profiles_comparison.csv` | All input O₃ profiles side by side |

---

## Spectral vs weighted (biological) UV

TUV produces **both**, and they answer different questions:

| Quantity | What it is | Where in output |
|----------|------------|-----------------|
| **Spectral irradiance** | Physical UV flux vs wavelength | Top table / `*_spectral.csv` |
| **Action spectrum** | Biological sensitivity vs wavelength (input data) | `DATAS1/dna.setlow.new` (Setlow DNA), other files in `DATAS1/` |
| **Dose rate (weighted)** | ∫ irradiance(λ) × action(λ) dλ | `Dose rates` table / `*_dose.csv` |

The **DNA action spectrum** (Setlow 1974, normalized at 300 nm, per quantum, converted to energy basis inside TUV) lives in:

```
DATAS1/dna.setlow.new
```

It is read by `swbiol.f` and applied internally; it is **not** written as a separate column in `usrout.txt`. Papers referring to “the DNA action spectrum in TUV” mean this weighting function applied to the computed spectral irradiance.

---

## Output files and directory layout

```
V5.4/
├── TUV.f, *.f              # Fortran source (see modifications below)
├── Makefile
├── run_TUV.py              # Main driver: install atmospheres, run TUV, archive
├── INPUTS/
│   ├── defin1–defin5       # Stock NCAR defaults
│   ├── Z1_pc               # Batch template (0.1% PAL settings, SZA sweep, bio weights)
│   ├── PI, One_pc, Ten_pc  # Other PAL input decks
│   └── usrinp              # Generated each batch run
├── DATAE1/ATM/
│   ├── ussa.temp           # Overwritten by run_TUV.py each run
│   ├── ussa.ozone          # Overwritten by run_TUV.py each run
│   └── *.orig_backup       # One-time backup of stock USSA files
├── DATAS1/                 # Biological action spectra (DNA, erythema, plant, …)
├── OUTPUT/compare/         # Archived runs, summary CSVs, plots
│   └── plot_tuv_comparison.py
├── README.md               # This file
└── README.txt              # Stock NCAR TUV documentation

/Users/gregcooke/{outfil}.txt   # Raw TUV output (parent of V5.4/)
/Users/gregcooke/tuvlog.txt     # Input log from last TUV run
```

Archive filenames follow:

```
{case}_{model}_atm48p2_sza0-85_{wstart-wstop}.{txt,csv,json}
```

Example: `pal0p1_waccm_atm48p2_sza0-85_w100-420_dose.csv`

---

## Changes compared with stock TUV 5.4

### Summary

| Area | Stock TUV 5.4 | This installation |
|------|---------------|-------------------|
| Execution | Interactive by default | **Batch mode** default (`intrct = .FALSE.`) |
| Atmosphere | Fixed USSA O₃ and T | **Dynamic** `ussa.temp` / `ussa.ozone` from WACCM6 & 1D models |
| Driver | Manual menu | **`run_TUV.py`** orchestrates multi-model comparison |
| Input decks | `defin1`–`defin5` | Added **`PI`, `One_pc`, `Ten_pc`, `Z1_pc`, `Z1_V`** for PAL cases |
| SZA control | Often via local time | **`lzenit = T`** with explicit 0–85° sweep in batch template |
| Wavelength start | Typically ≥ 280 nm in examples | **`wstart` allowed down to 100 nm** (for UVC / DNA short-wave tail) |
| Outputs | Single run text file | Parsed **CSV + JSON + comparison tables** |
| Analysis | External | **`plot_tuv_comparison.py`** for PAL vs metric plots |

### Detailed change list

#### 1. Batch / non-interactive execution

- **File:** `TUV.f` (lines ~237–239)
- **Change:** `intrct = .FALSE.` so `./tuv` reads `INPUTS/usrinp` without prompts.
- **Revert:** Set `intrct = .TRUE.` and recompile.

#### 2. Extended minimum wavelength

- **File:** `rdinp.f` (`chkval`, ~line 756)
- **Change:** Validation allows `wstart >= 100` nm (stock examples use 280 nm; exact stock lower limit may differ).
- **Purpose:** Include UVC when using `--wstart 100`; requires `nwint = -156` (predefined 120–735 nm grid) for windows extending below 205.8 nm.

#### 3. Interactive menu entries for custom input decks

- **File:** `rdinp.f` (~lines 123–127)
- **Change:** Added menu options for `PI`, `One_pc`, `Ten_pc`, `Z1_pc`, `Z1_pc_VULCAN`.
- **Purpose:** Quick manual access to Early-Earth PAL configurations.

#### 4. Custom input decks (`INPUTS/`)

New files (not in stock TUV distribution):

- `Z1_pc` — 0.1% PAL template: SZA sweep, eight biological weights enabled
- `One_pc`, `Ten_pc`, `PI` — other PAL levels with pre-set `o3col` examples
- `Z1_V` — VULCAN-oriented variant

These enable DNA damage (Setlow), UV index, and three plant-damage spectra by default (`nms = 8`).

#### 5. Python driver `run_TUV.py` (new file)

This is the main integration layer. It:

- Reads WACCM6/CAM NetCDF (`T`, `Z3`, `O3`, `lev`) with optional Gaussian lat mean or equatorial column
- Reads VULCAN `.vul`, Photochem text, Atmos `.dist`, Kasting `OUTPUT_PLOT.dat`
- Normalises altitude units (cm → km where needed) and fixes Atmos exponent typos
- Writes `ussa.temp` / `ussa.ozone`, computes O₃ column (DU)
- Generates `INPUTS/usrinp` from `Z1_pc` template
- Runs `make` + `./tuv`, parses output, archives to `OUTPUT/compare/`
- Builds `comparison_summary_*.csv` across all runs

#### 6. Runtime atmosphere files (`DATAE1/ATM/`)

- **Not Fortran changes**, but operational: `ussa.temp` and `ussa.ozone` are overwritten each run.
- Read by unchanged stock subroutines `vptmp.f` and `vpo3.f`.
- Backups preserved in `*.orig_backup`.

#### 7. Post-processing (`OUTPUT/compare/plot_tuv_comparison.py`, new)

Plots dose metrics vs PAL for selected SZAs from `comparison_summary*.csv`.

#### 8. Fortran code **not** modified for chemistry coupling

The following remain stock TUV 5.4 behaviour:

- Radiative transfer (`rtrans.f`, DISORT / 2-stream)
- Biological weighting (`swbiol.f`) — uses standard `DATAS1/` action spectra
- O₃ cross sections, Rayleigh, aerosol, cloud optics
- Schumann–Runge / Lyman-α bands (`la_srb.f`)

No changes were made to the DNA action spectrum data file itself (`DATAS1/dna.setlow.new`); it is the standard Setlow (1974) spectrum shipped with TUV 5.4.

---

## Where to find modifications in the code

Use this checklist when diffing against a pristine TUV 5.4 tarball:

| Location | Modification |
|----------|--------------|
| `TUV.f` | Batch mode default |
| `rdinp.f` | `wstart >= 100` validation; extra input-file menu entries |
| `INPUTS/Z1_pc`, `One_pc`, `Ten_pc`, `PI`, `Z1_V` | New input decks |
| `run_TUV.py` | **New** — entire WACCM6 / 1D integration |
| `OUTPUT/compare/plot_tuv_comparison.py` | **New** — plotting |
| `DATAE1/ATM/ussa.temp`, `ussa.ozone` | **Runtime data** (not source); restored from `*.orig_backup` if needed |
| `Makefile` | Unchanged except local compiler (`gfortran`) |

All other `.f` files match stock TUV 5.4 unless you have made additional local edits.

---

## Dependencies and external data paths

Configured at the top of `run_TUV.py` (edit if your paths differ):

| Variable | Default path | Contents |
|----------|--------------|----------|
| `CESM_DATA` | `/Users/gregcooke/CESM_data` | WACCM6/CAM `Earth_*pc_o2.cam.h0.*.nc` |
| `VIH_OUTPUT_PATH` | `/Users/gregcooke/VIH_cases/output` | VULCAN `.vul` files |
| `PHOTOCHEM_BASE` | `/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem` | Photochem 1D outputs |
| `ATMOS_BASE` | `/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Atmos` | Atmos 1D outputs |
| `KASTING_BASE` | `/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Kasting_1D_model` | Kasting 1D outputs |

WACCM variables expected in NetCDF:

- `T` — temperature (K)
- `Z3` — geopotential height (m)
- `O3` — ozone volume mixing ratio
- `lev` — pressure (hPa)
- `gw` — Gaussian latitude weights (for `--model waccm` global mean)

WACCM files for the four comparison cases are the **`CESM_data`** NetCDFs loaded in `Early_Earth.py` as `Pre_h0`, `Ten_h0`, `One_h0`, and `Zero1_h0` (see [Atmosphere sources](#atmosphere-sources-and-how-they-connect)). Do not point TUV at the older `H_escape/` `*_pc_h0` files.

---

## Troubleshooting

### `incorrect value for variable: tstart`

You are in **time mode** (`lzenit = F`) but entered an SZA. Set `lzenit = T`, ensure `tstop >= tstart`, and use 0–180°.

### Output files not in `V5.4/`

Stock TUV writes to the **parent directory**: `/Users/gregcooke/{outfil}.txt`. `run_TUV.py` copies them to `OUTPUT/compare/`.

### TUV STOP / segfault with short wavelengths

If `wstart < 205.8` nm, you **must** use `nwint = -156`. `run_TUV.py` sets this automatically. Also keep ≤ 150 levels in `ussa.temp` (handled in `write_ussa_temp`).

### `wstart` validation error below 100 nm

This build enforces `wstart >= 100` nm in `rdinp.f`.

### Restoring stock USSA atmosphere

```bash
cp DATAE1/ATM/ussa.temp.orig_backup DATAE1/ATM/ussa.temp
cp DATAE1/ATM/ussa.ozone.orig_backup DATAE1/ATM/ussa.ozone
```

---

## References

- Madronich, S., et al. — NCAR TUV model v5.4 (`README.txt`)
- Setlow, R. B. (1974) — DNA damage action spectrum (`DATAS1/dna.setlow.new`)
- WACCM6 / CESM — three-dimensional chemistry–climate columns
- VULCAN, Photochem, Atmos, Kasting 1D — equilibrium column atmospheres at ~48° SZA

For stock TUV variable definitions and FAQ, see `README.txt` sections 1–5.
