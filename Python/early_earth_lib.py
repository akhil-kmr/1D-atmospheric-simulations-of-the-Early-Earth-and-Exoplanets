"""
Early Earth analysis — shared functions, constants, and utilities.

Imported by Early_Earth.py. Do not run this file directly for full analysis.
"""

#%% imports
import io
import pickle
import re

import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib import rcParams
from photochem import EvoAtmosphere

try:
    import cartopy.crs as ccrs
except ImportError:
    ccrs = None

try:
    from scipy.interpolate import interp1d
except ImportError:
    interp1d = None
    
File =  "/Users/gregcooke/H_escape/b.e21.BWma1850.f19_g17.PC_b.SSPO.016.cam.h0.0320-0320.nc" #file name

Gw_file = xr.open_dataset(File,decode_times=False) #open the file and decode time as false

gw = Gw_file.gw.values

#%% paths and global style
P_path = '/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/'
VIH_OUTPUT_PATH = '/Users/gregcooke/VIH_cases/output/'

rcParams['font.weight'] = 'bold'
dpi = 200
Vulcan_factor = 2

#%% O3 budget and Photochem helpers

def compute_ox_production(
    mech_file="zahnle_earth.yaml",
    settings_file="input/settings_100pc.yaml",
    flux_file="input/Sun_0.0Ga.txt",
    pt_file="100pc_PT_profile/100%PAL_Sun_0.0Ga.txt",
    atol=1e-23,
    verbose=0
):
    """
    Run Photochem EvoAtmosphere and return odd-oxygen (Ox) production from O2 photolysis.

    Ox here means O + O3 (two O atoms per O2 photolysis event). The rate sums
    O2 + hv -> O + O and O2 + hv -> O + O1D branches from production_and_loss.

    Returns
    -------
    ox_prod_m3 : np.ndarray
        Ox production rate (molecules m^-3 s^-1).
    pressure_hpa : np.ndarray
        Pressure profile (hPa).
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


def compute_o3loss(
    mech_file="zahnle_earth.yaml",
    settings_file="input/settings_100pc.yaml",
    flux_file="input/Sun_0.0Ga.txt",
    pt_file="100pc_PT_profile/100%PAL_Sun_0.0Ga.txt",
    atol=1e-23,
    verbose=0
):
    """
    Run the Photochem EvoAtmosphere model and return the O3 loss rate from the
    Chapman reaction O + O3 -> O2 + O2 only.

    This is a separate one-dimensional check against the Photochem code, not
    the full VULCAN reaction-by-reaction budget built by vulcan_o3_production_loss.

    Returns
    -------
    o3_loss : np.ndarray
        O3 loss rate from O + O3 -> O2 + O2 (molecules m^-3 s^-1).
    pressure_hpa : np.ndarray
        Pressure profile (hPa).
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


def _parse_reaction_side(side):
    """
    Parse one side of a VULCAN reaction string into stoichiometric coefficients.

    Examples
    --------
    'O + O3'        -> {'O': 1, 'O3': 1}
    '2 O2 + M'      -> {'O2': 2}
    'O2 + O'        -> {'O2': 1, 'O': 1}

    The third-body symbol M is ignored here; it is handled separately when
    building the reaction rate.
    """
    from collections import defaultdict
    import re

    stoich = defaultdict(int)
    for part in side.split('+'):
        part = part.strip()
        if not part:
            continue
        match = re.match(r'^(\d+)\s+(\S+)$', part)
        if match:
            stoich[match.group(2)] += int(match.group(1))
        else:
            stoich[part] += 1
    return dict(stoich)


def _vulcan_reaction_rate_cm3(var, species, rid, reverse=False):
    """
    Evaluate the volumetric rate of one VULCAN reaction at every model level.

    The rate is computed the same way VULCAN does in its rate laws:
        rate = k * [reactant1]^nu1 * [reactant2]^nu2 * ... * M (if M appears)

    Parameters
    ----------
    var : dict
        vulcan_ds['variable'] from a saved .vul file.
    species : list
        Species names (same order as var['y'] columns).
    rid : int
        Reaction index in var['k'] and (for forward steps) var['Rf'].
    reverse : bool
        If False, use reactants on the left of var['Rf'][rid].
        If True, use products on the right of var['Rf'][rid-1] with k[rid]
        (VULCAN stores reverse rate coefficients on even indices).

    Returns
    -------
    rate : np.ndarray
        Reaction rate in molecules cm^-3 s^-1 at each vertical level, or None
        if a required species is missing from the network.

    Notes
    -----
    - var['y'] is number density in molecules cm^-3.
    - var['k'][rid] is already the temperature/pressure-dependent rate coefficient
      from the run (including J-values for photolysis, in s^-1 or cm^3 mol^-1 s^-1).
    - Total density sum(var['y'], axis=1) is used when M appears in the reaction.
    """
    if reverse:
        if rid not in var['k'] or (rid - 1) not in var['Rf']:
            return None
        rf = var['Rf'][rid - 1]
        side = rf.split('->', 1)[1].strip()
    else:
        if rid not in var['Rf']:
            return None
        rf = var['Rf'][rid]
        side = rf.split('->', 1)[0].strip()

    k = np.asarray(var['k'][rid], dtype=float)
    reactants = _parse_reaction_side(side)
    rate = k.copy()
    number_density = var['y'].sum(axis=1)

    for sp, coeff in reactants.items():
        if sp == 'M':
            rate *= number_density
        elif sp not in species:
            return None
        else:
            rate *= var['y'][:, species.index(sp)] ** coeff
    return rate


def vulcan_o3_production_loss(vulcan_ds, species_list, to_m3=True):
    """
    Build a complete instantaneous O3 production and loss budget from a VULCAN run.

    This loops over every reaction in var['Rf'], computes its rate from var['k']
    and var['y'], and accumulates only the part that creates or destroys O3.
    Reverse steps (even reaction indices in k, not listed in Rf) are included.

    For each reaction, if O3 appears on the left (reactants), the O3 loss term is:
        nu_O3 * rate   (e.g. O + O3 -> ... has nu_O3 = 1)
    If O3 appears on the right (products), the O3 production term is:
        nu_O3 * rate

    Parameters
    ----------
    vulcan_ds : dict
        Loaded VULCAN output (e.g. One_pc_V_482SZA from Read_O3_Run).
    species_list : list
        Species names in column order, typically vulcan_ds['variable']['species'].
    to_m3 : bool
        If True (default), multiply all rates by 1e6 to convert
        molecules cm^-3 s^-1 -> molecules m^-3 s^-1.

    Returns
    -------
    dict with keys:
        prod : np.ndarray
            Total O3 production rate vs pressure (sum of all producing channels).
        loss : np.ndarray
            Total O3 loss rate vs pressure (sum of all destroying channels).
        net : np.ndarray
            prod - loss (positive means net O3 production at that level).
        prod_by_rxn : dict
            {reaction_string: rate_profile} for each O3-producing channel.
        loss_by_rxn : dict
            {reaction_string: rate_profile} for each O3-destroying channel.
            Reverse reactions are labelled 'REV: <forward reaction>'.
        pco : np.ndarray
            Pressure in hPa (vulcan_ds['atm']['pco'] / 1e3).
        rate_unit : str
            Label for plot axes when to_m3=True.

    Example
    -------
        budget = vulcan_o3_production_loss(One_pc_V_482SZA, One_pc_spec_482SZA)
        budget['loss_by_rxn']['O + O3 -> O2 + O2']  # Chapman loss profile
    """
    var = vulcan_ds['variable']
    nz = var['y'].shape[0]
    prod = np.zeros(nz)
    loss = np.zeros(nz)
    prod_by_rxn = {}
    loss_by_rxn = {}

    for rid, rf in var['Rf'].items():
        if '->' not in rf:
            continue
        left, right = [part.strip() for part in rf.split('->', 1)]
        nu_left = _parse_reaction_side(left)
        nu_right = _parse_reaction_side(right)
        o3_left = nu_left.get('O3', 0)
        o3_right = nu_right.get('O3', 0)
        if o3_left == 0 and o3_right == 0:
            continue

        rate = _vulcan_reaction_rate_cm3(var, species_list, rid, reverse=False)
        if rate is None:
            continue
        if o3_left > 0:
            term = o3_left * rate
            loss += term
            loss_by_rxn[rf] = term
        if o3_right > 0:
            term = o3_right * rate
            prod += term
            prod_by_rxn[rf] = term

    # Reverse rates: k on even indices, reactants = forward products
    for rid in var['k'].keys():
        if rid % 2 != 0 or rid in var['Rf'] or (rid - 1) not in var['Rf']:
            continue
        rf_fwd = var['Rf'][rid - 1]
        if '->' not in rf_fwd:
            continue
        left, right = [part.strip() for part in rf_fwd.split('->', 1)]
        nu_left = _parse_reaction_side(left)
        nu_right = _parse_reaction_side(right)
        o3_on_products = nu_right.get('O3', 0)
        o3_on_reactants = nu_left.get('O3', 0)

        rate = _vulcan_reaction_rate_cm3(var, species_list, rid, reverse=True)
        if rate is None:
            continue
        label = 'REV: ' + rf_fwd
        if o3_on_products > 0:
            term = o3_on_products * rate
            loss += term
            loss_by_rxn[label] = term
        if o3_on_reactants > 0:
            term = o3_on_reactants * rate
            prod += term
            prod_by_rxn[label] = term

    scale = 1e6 if to_m3 else 1.0
    return {
        'prod': prod * scale,
        'loss': loss * scale,
        'net': (prod - loss) * scale,
        'prod_by_rxn': {k: v * scale for k, v in prod_by_rxn.items()},
        'loss_by_rxn': {k: v * scale for k, v in loss_by_rxn.items()},
        'pco': vulcan_ds['atm']['pco'] / 1e3,
        # Instantaneous O3 chemical rate; y in molecules cm^-3, k in VULCAN cgs units
        'rate_unit': 'molecules m^-3 s^-1 (O3)',
        'rate_unit_cm3': 'molecules cm^-3 s^-1 (O3)',
    }


def _normalize_reaction_label(reaction_string):
    """Unify VULCAN (->) and Photochem (=>) reaction labels."""
    return reaction_string.replace('REV: ', '').replace('=>', ' -> ').strip()


def _reaction_species_tokens(reaction_string):
    """
    Extract the set of chemical species from a reaction label.

    Used to test whether a reaction belongs to the NOx or HOx family. Stoichiometric
    coefficients are not returned (unlike _parse_reaction_side); only species names.
    The third-body symbol M and photolysis placeholder hv are excluded.
    """
    label = _normalize_reaction_label(reaction_string)
    if '->' not in label:
        return set()
    left, right = [part.strip() for part in label.split('->', 1)]
    skip = {'M', 'hv'}
    tokens = set()
    for side in (left, right):
        for part in side.split('+'):
            part = part.strip()
            if part and part not in skip:
                tokens.add(part)
    return tokens


# Species used to tag O3 channels for NOx vs HOx family plots (must also contain O3).
O3_NOX_SPECIES = {'N', 'NO', 'NO2'}
O3_HOX_SPECIES = {'H', 'OH', 'HO2', 'H2O2'}


def o3_channels_in_family(budget, family='nox'):
    """
    Split the O3 budget into channels that involve NOx or HOx chemistry.

    A reaction is included only if:
      1. It already appears in budget['prod_by_rxn'] or budget['loss_by_rxn']
         (so it must change O3), and
      2. O3 is among the species on either side of the arrow, and
      3. At least one species from O3_NOX_SPECIES or O3_HOX_SPECIES is present.

    Parameters
    ----------
    budget : dict
        Output of vulcan_o3_production_loss.
    family : str
        'nox' (N, NO, NO2) or 'hox' (H, OH, HO2, H2O2).

    Returns
    -------
    prod : dict
        Subset of budget['prod_by_rxn'] matching the family.
    loss : dict
        Subset of budget['loss_by_rxn'] matching the family.

    Notes
    -----
    Reactions that destroy O3 but only involve e.g. CH3 or S on the other side
    are not in either family. Photolysis O3 -> O2 + O is not HOx by this definition
    (no H, OH, HO2, or H2O2 on the arrow).
    """
    if family.lower() == 'nox':
        species = O3_NOX_SPECIES
    elif family.lower() == 'hox':
        species = O3_HOX_SPECIES
    else:
        raise ValueError("family must be 'nox' or 'hox'")

    def matches(reaction_string):
        tokens = _reaction_species_tokens(reaction_string)
        return ('O3' in tokens) and bool(tokens & species)

    prod = {k: v for k, v in budget['prod_by_rxn'].items() if matches(k)}
    loss = {k: v for k, v in budget['loss_by_rxn'].items() if matches(k)}
    return prod, loss


def o3_family_total_rates(budget, family='nox'):
    """
    Sum all O3 production and loss rates in a NOx or HOx family.

    This adds every matching channel from o3_channels_in_family vertically,
    giving two profiles: total NOx (or HOx) O3 production and total loss.

    Parameters
    ----------
    budget : dict
        Output of vulcan_o3_production_loss.
    family : str
        'nox' or 'hox'.

    Returns
    -------
    prod_total, loss_total : np.ndarray
        1D arrays vs pressure, same units as budget['prod'] and budget['loss'].
    """
    prod_ch, loss_ch = o3_channels_in_family(budget, family=family)
    prod_total = np.zeros_like(budget['prod'])
    loss_total = np.zeros_like(budget['loss'])
    for rate in prod_ch.values():
        prod_total += rate
    for rate in loss_ch.values():
        loss_total += rate
    return prod_total, loss_total


def _is_o3_photolysis_loss(reaction_string):
    """
    Return True if the reaction is O3 photolysis (O3 alone on the reactant side).

    Examples that return True:  O3 -> O2 + O,  O3 + hv -> O + O2
    Examples that return False: OH + O3 -> HO2 + O2,  O + O3 -> O2 + O2
    """
    label = _normalize_reaction_label(reaction_string)
    if '->' not in label:
        return False
    left = label.split('->', 1)[0].strip()
    stoich = _parse_reaction_side(left)
    stoich = {sp: nu for sp, nu in stoich.items() if sp not in ('M', 'hv')}
    return stoich == {'O3': 1}


def o3_photolysis_total_loss(budget):
    """
    Sum O3 loss from photolysis branches only (O3 is the sole reactant).

    These are reactions of the form O3 -> products with no other reactant,
    where the rate is k * [O3] and k is a photolysis frequency (s^-1) from
    the radiative transfer in VULCAN.

    Returns
    -------
    loss_total : np.ndarray
        Total photolytic O3 loss rate vs pressure (molecules m^-3 s^-1 if
        budget was built with to_m3=True).
    """
    loss_total = np.zeros_like(budget['loss'])
    for rxn, rate in budget['loss_by_rxn'].items():
        if _is_o3_photolysis_loss(rxn):
            loss_total += rate
    return loss_total


def o3_hox_catalytic_total_loss(budget):
    """
    Sum HOx-family O3 loss, excluding photolysis.

    Starts from o3_channels_in_family(..., 'hox') and drops any channel where
    O3 alone is photolysed. The result is comparable to NOx catalytic loss:
    both are bimolecular (or termolecular) chemistry, not light-driven O3 -> O2 + ...

    In the current Earth network, HOx catalytic loss is usually from
    O3 + H, OH + O3, and HO2 + O3 (H2O2 does not appear in any O3 reaction).

    Returns
    -------
    loss_total : np.ndarray
        HOx catalytic O3 loss rate vs pressure.
    """
    _, hox_loss = o3_channels_in_family(budget, family='hox')
    loss_total = np.zeros_like(budget['loss'])
    for rxn, rate in hox_loss.items():
        if not _is_o3_photolysis_loss(rxn):
            loss_total += rate
    return loss_total


def o3_chapman_total_loss(budget, reaction='O + O3 -> O2 + O2'):
    """
    Extract the Chapman O3 loss channel O + O3 -> O2 + O2.

    Parameters
    ----------
    budget : dict
        Output of vulcan_o3_production_loss.
    reaction : str
        Reaction label exactly as stored in budget['loss_by_rxn'].

    Returns
    -------
    loss_total : np.ndarray
        Rate profile for that reaction, or zeros if the reaction is absent.
    """
    loss_total = np.zeros_like(budget['loss'])
    if reaction in budget['loss_by_rxn']:
        loss_total += budget['loss_by_rxn'][reaction]
    return loss_total


def o3_loss_by_category(budget):
    """
    Split total O3 loss into non-overlapping categories for comparison plots.

    Categories (all molecules m^-3 s^-1, same as budget['loss']):
        nox          - loss via N, NO, or NO2 on the arrow (see O3_NOX_SPECIES)
        hox_catalytic - loss via H, OH, HO2, or H2O2, excluding photolysis
        photolysis   - O3 -> ... with O3 the only reactant
        chapman      - O + O3 -> O2 + O2 only
        other        - remainder (CH3+O3, sulfur chemistry, etc.)
        total        - full budget['loss']

    The categories are mutually exclusive and sum to total loss (except for
    tiny numerical round-off, clamped in 'other').

    Returns
    -------
    dict with keys: total, nox, hox_catalytic, photolysis, chapman, other,
    and nox_prod (NOx O3 production, useful if you need the NOx source term).
    """
    nox_prod, nox_loss = o3_family_total_rates(budget, 'nox')
    hox_catalytic = o3_hox_catalytic_total_loss(budget)
    photolysis = o3_photolysis_total_loss(budget)
    chapman = o3_chapman_total_loss(budget)
    total = budget['loss']
    other = total - nox_loss - hox_catalytic - photolysis - chapman
    # Guard against small negative other from numerical edge cases
    other = np.maximum(other, 0.0)
    return {
        'total': total,
        'nox': nox_loss,
        'hox_catalytic': hox_catalytic,
        'photolysis': photolysis,
        'chapman': chapman,
        'other': other,
        'nox_prod': nox_prod,
    }


# WACCM6 CAM h0: bimolecular rate coefficients named <reactant>_O3 (units cm^3 molecules^-1 s^-1).
WACCM_O3_LOSS_CHANNEL_TABLE = (
    ('NO_O3', 'NO', 'NO + O3 -> NO2 + O2'),
    ('NO2_O3', 'NO2', 'NO2 + O3 -> NO3 + O2'),
    ('H_O3', 'H', 'O3 + H -> OH + O2'),
    ('OH_O3', 'OH', 'OH + O3 -> HO2 + O2'),
    ('HO2_O3', 'HO2', 'HO2 + O3 -> OH + O2 + O2'),
    ('O_O3', 'O', 'O + O3 -> O2 + O2'),
)
WACCM_O3_LOSS_REQUIRED_VARS = tuple(row[0] for row in WACCM_O3_LOSS_CHANNEL_TABLE)


def waccm_has_o3_loss_diagnostics(cam_ds):
    """True if the CAM history file includes *_O3 bimolecular rate coefficients."""
    return all(v in cam_ds for v in WACCM_O3_LOSS_REQUIRED_VARS)


def _waccm_number_density_m3(cam_ds, species, latitude=np.arange(0, 96)):
    """Number density [m^-3] from mol/mol, hybrid mid-level pressure [hPa], and T [K]."""
    k_b = 1.381e-23
    n = cam_ds[species] * cam_ds.lev * 100.0 / (cam_ds.T * k_b)
    return np.asarray(LWAV_spec(n, latitude=latitude), dtype=float)


def _waccm_bimolecular_loss_m3(k_cm3_molec_s, n_a_m3, n_b_m3):
    """
    O3 loss rate [molecules m^-3 s^-1] from k [cm^3 molecules^-1 s^-1] and densities [m^-3].
    """
    return np.asarray(k_cm3_molec_s, dtype=float) * np.asarray(n_a_m3, dtype=float) * np.asarray(n_b_m3, dtype=float) / 1e6


def waccm_o3_catalytic_budget(cam_ds, latitude=np.arange(0, 96)):
    """
    Build an O3 catalytic-loss budget from WACCM6 CAM h0 diagnostics.

    Uses latitude-weighted zonal means (LWAV_spec) of <reactant>_O3 rate coefficients
    multiplied by number densities. Photolysis and other loss are not in standard h0
    *_O3 outputs; returned budget['loss'] is NOx + HOx + Chapman only.

    Returns None if required variables are missing (older cam.h0 files).
    Otherwise returns a dict compatible with o3_loss_by_category and the O3
    comparison plot helpers (same keys as vulcan_o3_production_loss).
    """
    if not waccm_has_o3_loss_diagnostics(cam_ds):
        return None

    pco = np.asarray(cam_ds.lev.values, dtype=float)
    n_o3 = _waccm_number_density_m3(cam_ds, 'O3', latitude=latitude)
    loss_by_rxn = {}
    nox = np.zeros_like(pco, dtype=float)
    hox = np.zeros_like(pco, dtype=float)
    chapman = np.zeros_like(pco, dtype=float)

    for kvar, species, label in WACCM_O3_LOSS_CHANNEL_TABLE:
        k = np.asarray(LWAV_spec(cam_ds[kvar], latitude=latitude), dtype=float)
        n_sp = _waccm_number_density_m3(cam_ds, species, latitude=latitude)
        loss = _waccm_bimolecular_loss_m3(k, n_sp, n_o3)
        loss_by_rxn[label] = loss
        if species in ('NO', 'NO2'):
            nox += loss
        elif species in ('H', 'OH', 'HO2'):
            hox += loss
        elif species == 'O':
            chapman += loss

    total = nox + hox + chapman
    return {
        'prod': np.zeros_like(pco, dtype=float),
        'loss': total,
        'net': -total,
        'prod_by_rxn': {},
        'loss_by_rxn': loss_by_rxn,
        'pco': pco,
        'rate_unit': 'molecules m^-3 s^-1 (O3, WACCM bimolecular)',
        'rate_unit_cm3': 'molecules cm^-3 s^-1 (O3, WACCM bimolecular)',
        'source': 'WACCM6',
    }


# Fractional O3 loss plots: NOx + HOx catalytic + Chapman only (excludes photolysis & other).
O3_LOSS_FRAC_NHC_CATEGORIES = (
    ('nox', 'NO' + 'x', 'k'),
    ('hox_catalytic', 'HO' + 'x', 'b'),
    ('chapman', 'Chapman', 'r'),
)


def o3_loss_fractional_nhc(budget):
    """
    Fractional O3 loss among NOx, HOx catalytic, and Chapman channels only.

    At each level, the three fractions sum to 1 where their combined loss > 0.
    Photolysis and other loss are excluded from numerator and denominator.
    """
    cat = o3_loss_by_category(budget)
    subset_total = cat['nox'] + cat['hox_catalytic'] + cat['chapman']
    frac = {}
    with np.errstate(divide='ignore', invalid='ignore'):
        for key, _label, _color in O3_LOSS_FRAC_NHC_CATEGORIES:
            frac[key] = np.where(subset_total > 0, cat[key] / subset_total, 0.0)
    frac['subset_total'] = subset_total
    frac['total'] = cat['total']
    return frac


def plot_o3_loss_fraction(
    budget,
    title,
    ax,
    ylim_press=None,
    xlim_frac=(1e-3, 1.0),
):
    """
    Draw NOx / HOx / Chapman fractional loss (sums to 1) vs pressure on a log fraction axis.
    """
    frac = o3_loss_fractional_nhc(budget)
    p = budget['pco']
    xlo, xhi = xlim_frac

    for key, label, color in O3_LOSS_FRAC_NHC_CATEGORIES:
        f_plot = np.asarray(frac[key], dtype=float)
        # Log axis cannot display exactly zero; clip only for rendering, not for masking data
        f_plot = np.where(f_plot > 0, f_plot, np.nan)
        ax.plot(f_plot, p, color=color, lw=2, label=label)

    ax.set_xscale('log')
    ax.set_xlim(xlo, xhi)
    ax.set_xlabel(
        'Fraction of NO' + sub('x') + '+HO' + sub('x') + '+Chapman loss',
        fontsize=11, weight='bold',
    )
    ax.set_title(title, fontsize=12, weight='bold')
    if ylim_press is not None:
        ax.set_ylim(ylim_press)
    ax.set_yscale('log')


def plot_o3_nox_hox_fraction_vulcan_vs_photo(
    budget_vulcan,
    budget_photo,
    ax,
    title,
    xlim_frac=(1e-3, 1.0),
    budget_waccm=None,
):
    """Plot NOx (black) and HOx (blue) fractions: VULCAN solid, Photochem dashed, WACCM dotted."""
    frac_v = o3_loss_fractional_nhc(budget_vulcan)
    frac_p = o3_loss_fractional_nhc(budget_photo)
    p_v = budget_vulcan['pco']
    p_p = budget_photo['pco']
    xlo, xhi = xlim_frac

    def _line(f, p, color, ls, label):
        y = np.asarray(f, dtype=float)
        y = np.where(y > 0, y, np.nan)
        ax.plot(y, p, color=color, ls=ls, lw=2, label=label)

    _line(frac_v['nox'], p_v, 'k', '-', 'NO' + sub('x') + ' VULCAN')
    _line(frac_v['hox_catalytic'], p_v, 'b', '-', 'HO' + sub('x') + ' VULCAN')
    _line(frac_p['nox'], p_p, 'k', '--', 'NO' + sub('x') + ' Photochem')
    _line(frac_p['hox_catalytic'], p_p, 'b', '--', 'HO' + sub('x') + ' Photochem')

    if budget_waccm is not None:
        frac_w = o3_loss_fractional_nhc(budget_waccm)
        p_w = budget_waccm['pco']
        _line(frac_w['nox'], p_w, 'k', ':', 'NO' + sub('x') + ' WACCM6')
        _line(frac_w['hox_catalytic'], p_w, 'b', ':', 'HO' + sub('x') + ' WACCM6')

    ax.set_xscale('log')
    ax.set_xlim(xlo, xhi)
    ax.set_yscale('log')
    ax.set_title(title, fontsize=12, weight='bold')
    ax.set_xlabel(
        'Fraction of NO' + sub('x') + '+HO' + sub('x') + '+Chapman loss',
        fontsize=11, weight='bold',
    )


def top_o3_channels(budget, level_index, n=8):
    """
    Find the largest O3 production and loss channels at one vertical level.

    Useful for printed diagnostics (e.g. at the O3 mixing-ratio peak).

    Parameters
    ----------
    budget : dict
        Output of vulcan_o3_production_loss.
    level_index : int
        Vertical index (0 = top of model, same ordering as var['y']).
    n : int
        Number of channels to return in each list.

    Returns
    -------
    prod_top, loss_top : list of (reaction_string, rate) tuples
        Rates are scalars at level_index in molecules m^-3 s^-1, sorted
        largest first.
    """
    def at_level(channel_dict):
        return sorted(
            ((name, float(profile[level_index])) for name, profile in channel_dict.items()),
            key=lambda item: item[1],
            reverse=True,
        )[:n]

    return at_level(budget['prod_by_rxn']), at_level(budget['loss_by_rxn'])


def photochem_o3_production_loss(
    mech_file,
    settings_file,
    flux_file,
    pt_file,
    atol=1e-23,
    verbose=0,
    to_m3=True,
):
    """
    Build an O3 production/loss budget from Photochem EvoAtmosphere.

    Uses production_and_loss('O3', usol), which already returns per-reaction O3
    rates in molecules cm^-3 s^-1. Reaction labels use Photochem syntax (=>);
    they are normalised to ' -> ' for compatibility with the VULCAN budget helpers.

    Parameters
    ----------
    mech_file, settings_file, flux_file, pt_file : str
        Paths passed to EvoAtmosphere (pt_file is the saved atmosphere, e.g.
        Earth_100pc_48.2.txt, the same file used for Photo_PI etc.).
    atol, verbose : passed to pc.var.
    to_m3 : bool
        If True, multiply rates by 1e6 (molecules m^-3 s^-1).

    Returns
    -------
    dict with keys prod, loss, net, prod_by_rxn, loss_by_rxn, pco, rate_unit
    (same structure as vulcan_o3_production_loss).
    """
    pc = EvoAtmosphere(mech_file, settings_file, flux_file, pt_file)
    pc.var.atol = atol
    pc.var.verbose = verbose
    pl = pc.production_and_loss('O3', pc.wrk.usol)
    scale = 1e6 if to_m3 else 1.0

    prod_by_rxn = {
        _normalize_reaction_label(rx): pl.production[:, i] * scale
        for i, rx in enumerate(pl.production_rx)
    }
    loss_by_rxn = {
        _normalize_reaction_label(rx): pl.loss[:, i] * scale
        for i, rx in enumerate(pl.loss_rx)
    }
    prod = pl.production.sum(axis=1) * scale
    loss = pl.loss.sum(axis=1) * scale
    sol = pc.mole_fraction_dict()
    # sol['pressure'] in dyne cm^-2; /1e3 gives hPa (matches Photo_* press * 1e3)
    pco = sol['pressure'] / 1e3

    return {
        'prod': prod,
        'loss': loss,
        'net': (prod - loss),
        'prod_by_rxn': prod_by_rxn,
        'loss_by_rxn': loss_by_rxn,
        'pco': pco,
        'rate_unit': 'molecules m^-3 s^-1 (O3)',
        'rate_unit_cm3': 'molecules cm^-3 s^-1 (O3)',
    }


def photochem_o3_budget_for_folder(folder, sza_tag='48.2', atol=1e-23, verbose=0, to_m3=True):
    """
    Convenience wrapper: O3 budget for a standard Early-Earth Photochem case folder.

    folder examples: '100pc', '10pc', '1pc', '0.1pc'
    """
    pt_names = {
        '100pc': 'Earth_100pc_{}.txt',
        '10pc': 'Earth_10pc_{}.txt',
        '1pc': 'Earth_1pc_{}.txt',
        '0.1pc': 'Earth_0.1pc_{}.txt',
    }
    if folder not in pt_names:
        raise ValueError('folder must be one of {}'.format(list(pt_names.keys())))
    base = P_path + folder + '/'
    pt_file = base + pt_names[folder].format(sza_tag)
    return photochem_o3_production_loss(
        mech_file=base + 'zahnle_earth.yaml',
        settings_file=base + 'settings_100pc.yaml',
        flux_file=base + 'Sun_0.0Ga.txt',
        pt_file=pt_file,
        atol=atol,
        verbose=verbose,
        to_m3=to_m3,
    )


#%% plotting helpers (thick_axes, sub, sup, cbar, LWAV)

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

def cbar(model, label = True, clabel = 'O'+ sub(3) +' column [DU]', orientation='vertical', tick = False, ticks = [2,3], shrink = 1):
    cbar = plt.colorbar(model, orientation = orientation, shrink = shrink)
    if (label == True):
        clabel = clabel
        cbar.set_label(clabel,size=15, weight = 'bold')
    cbar.ax.tick_params(labelsize=15)
    if (tick == True):
        cbar.set_ticks(ticks)


#%% VULCAN I/O

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


def Read_O3_Run(file_path=''):
    """
    Load a saved VULCAN .vul pickle and return the dataset plus species list.

    Parameters
    ----------
    file_path : str
        Filename under VIH_cases/output/ (e.g. 'Earth_1pc_o2_1e12s_48.2SZA_WPT_1rtol.vul').

    Returns
    -------
    DS : dict
        Full VULCAN output (use DS['variable'], DS['atm'], etc.).
    spec : list
        Species names in the same order as DS['variable']['y'] columns.
    """
    DS = read_pickle_file(VIH_OUTPUT_PATH + file_path)
    spec = DS['variable']['species']
    return DS, spec


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
        if current_header and strip_line and strip_line[0].isdigit():
            current_lines.append(strip_line)
        # If we hit a blank line or new text, save the block
        elif current_header and current_lines and not strip_line:
            df = pd.read_csv(io.StringIO("\n".join(current_lines)), sep='\s+', names=headers)
            blocks[current_header] = df
            current_header = None
            current_lines = []

    # Catch the last block if file doesn't end in blank line
    if current_lines and current_header:
        df = pd.read_csv(io.StringIO("\n".join(current_lines)), sep='\s+', names=headers)
        blocks[current_header] = df

    return blocks


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


def Atmos_photo(file_path):
    
    df = pd.read_csv(file_path, names = ['Z', 'PO2_1', "PO2_2"], delim_whitespace=True)
    JO2 = 2*(df["PO2_1"]+df["PO2_2"])
    return df["Z"], JO2


def Atmos_dens(DS):
    k = 1.381e-23
    dens = DS["PRESS"]*1e5/(k*DS["TEMP"])
    return dens


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


def rayleigh_sigma_nm(wavelength_nm: float) -> float:
    """Same closure as mo_jshort.F90 rayleigh_sigma_nm (cm^2 per molecule)."""
    w_um = max(wavelength_nm * 1e-3, 1e-6)
    w2 = w_um * w_um
    w4 = w2 * w2
    return 4.006e-28 * (1.0 + 0.0113 / w2 + 0.00013 / w4) / w4


def air_density_molec_cm3(P_pa: np.ndarray, T_k: np.ndarray) -> np.ndarray:
    """Ideal gas: n = P/(k_B T), return molecules cm^-3."""
    return (P_pa / (k_B * np.maximum(T_k, 1.0))) / 1e6


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


def convert_photons_to_energy(wav, flux):
    
    # Physical constants
    h = 6.626e-34  # Planck constant (J·s)
    c = 2.998e8    # speed of light (m/s)

    # Convert wavelength to meters
    wav = wav * 1e-10

    # Conversion: photons/cm²/s → W/m²/nm
    irradiance = flux * (h * c / wav) * 1e4
    
    return irradiance


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
        data['variable']['J_sp']['O2',0]# +
        #data['variable']['J_sp']['O2',1] +
        #data['variable']['J_sp']['O2',2]
    )
    
    # layer thickness (cm)
    dz = data['atm']['dz']/100
    
    # local production per layer (cm^-2 s^-1)
    layer_prod = Vulcan_factor * JO2 * O2_dens * dz * 3/8
    
    # cumulative integral from top of atmosphere downward
    int_JO2 = layer_prod[::-1].cumsum()[::-1]
    
    return int_JO2


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


def read_usrout(filepath):
    """
    Read UVSPEC-style output file.

    Returns
    -------
    spectral_df : pandas.DataFrame
    dose_df : pandas.DataFrame
    """

    with open(filepath, "r") as f:
        text = f.read()

    # ================================================================
    # Spectral irradiance table
    # ================================================================

    spectral_pattern = r"wc, nm.*?\n(.*?)\n-+"

    match = re.search(spectral_pattern, text, re.S)

    if not match:
        raise ValueError(f"Could not find spectral table in {filepath}")

    lines = [
        line for line in match.group(1).splitlines()
        if line.strip()
    ]

    spectral_data = [
        [float(x) for x in line.split()]
        for line in lines
    ]

    spectral_columns = [
        "wavelength_nm",
        "sza_0",
        "sza_10",
        "sza_20",
        "sza_30",
        "sza_40",
        "sza_50",
        "sza_60",
        "sza_70",
        "sza_80",
        "sza_90",
    ]

    spectral_df = pd.DataFrame(
        spectral_data,
        columns=spectral_columns
    )

    # ================================================================
    # Dose rate table
    # ================================================================

    dose_pattern = r"sza, deg\..*?\n(.*?)\n-+"

    match = re.search(dose_pattern, text, re.S)

    if not match:
        raise ValueError(f"Could not find dose table in {filepath}")

    lines = [
        line for line in match.group(1).splitlines()
        if line.strip()
    ]

    dose_data = [
        [float(x) for x in line.split()]
        for line in lines
    ]

    dose_columns = [
        "sza_deg",
        "UVB_280_315",
        "UVBstar_280_320",
        "UVA_315_400",
        "DNA_damage",
        "UV_index",
        "P_Dam_C_1971",
        "P_Dam_FC_2003",
        "P_Dam_FC_2003_ext390",
    ]

    dose_df = pd.DataFrame(
        dose_data,
        columns=dose_columns
    )

    return spectral_df, dose_df


def _o3_loss_fraction_figure(
    case_list,
    budgets_dict,
    model_name,
    out_stem,
    ylim_press=(1000.0, 1e-2),
    xlim_frac=(1e-3, 1.0),
):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharey=True)
    axes = axes.flatten()
    for idx, case in enumerate(case_list):
        label = case[0]
        ax = axes[idx]
        plot_o3_loss_fraction(
            budgets_dict[label],
            label,
            ax,
            ylim_press=ylim_press,
            xlim_frac=xlim_frac,
        )
        thick_axes(top=True)
        if idx in (0, 2):
            ax.set_ylabel('Pressure [hPa]', fontsize=13, weight='bold')
        else:
            ax.tick_params(labelleft=False)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles, labels, loc='lower center', ncol=3, fontsize=10,
        frameon=False, bbox_to_anchor=(0.5, -0.02),
    )
    fig.suptitle(
        'O' + sub(3) + ' loss: NO' + sub('x') + ', HO' + sub('x') + ', Chapman fractions ('
        + model_name + ', SZA = 48.2°)',
        fontsize=16, weight='bold', y=1.02,
    )
    plt.tight_layout()
    plt.savefig(
        '/Users/gregcooke/python_output/{}.png'.format(out_stem),
        dpi=dpi,
        bbox_inches='tight',
    )

