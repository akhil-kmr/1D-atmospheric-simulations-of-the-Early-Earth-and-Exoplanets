"""
Convert saved Photochem steady-state atmosphere files to Chempath HDF5 input.

Uses your existing Photochem install (no chempath fork required). Reaction rates
come from pc.wrk.rx_rates; non-chemical steady-state terms (transport, rainout,
deposition, etc.) are folded into per-species transport pseudo-reactions so the
rate budget closes at each level.

Example:
    python steadystate_to_chempath.py --sza 48.2 --layers 20 30 40
    python steadystate_to_chempath.py --sza 48.2 --all-layers
"""

from __future__ import annotations

import argparse
import pathlib

import h5py
import numpy as np
from photochem import EvoAtmosphere

F128 = np.longdouble

PHOTOCHEM_ROOT = pathlib.Path('/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem')
ANALYSIS_ROOT = PHOTOCHEM_ROOT / 'chempath_analysis'

MECH_FILE = PHOTOCHEM_ROOT / 'Old sims' / 'zahnle_earth.yaml'
FLUX_FILE = PHOTOCHEM_ROOT / 'Old sims' / 'input' / 'Sun_0.0Ga.txt'

SETTINGS_BY_PAL_SZA = {
    '100pc': {
        '45': ANALYSIS_ROOT / 'settings' / 'settings_100pc_sza45.yaml',
        '48.2': PHOTOCHEM_ROOT / 'Proxima_Centauri' / 'input' / 'settings_100pc.yaml',
        '60': PHOTOCHEM_ROOT / 'Old sims' / 'input' / 'settings_100pc.yaml',
    },
    '1pc': {
        '48.2': PHOTOCHEM_ROOT / 'Proxima_Centauri' / 'input' / 'settings_1pc.yaml',
    },
}

DEFAULT_OUTPUT = ANALYSIS_ROOT / 'chempath_input'


def get_reactant_indices(species_names, reaction_equations):
    """Return 0-based reactant species indices for each reaction.

    Indices refer to ``pc.wrk.densities`` rows (includes hv and M pseudo-species).
    """
    r1 = []
    r2 = []
    for equation in reaction_equations:
        reaction = equation.replace(' ', '').replace('=>', '=')
        reactants = reaction.split('=')[0].split('+')
        r1.append(species_names.index(reactants[0]))
        if len(reactants) == 1:
            r2.append(species_names.index('hv'))
        else:
            r2.append(species_names.index(reactants[1]))
    return np.array(r1), np.array(r2)


def get_sij(species_names, reaction_equations):
    """Stoichiometry matrix: positive = produced, negative = consumed."""
    ni = len(species_names)
    nj = len(reaction_equations)
    sij = np.zeros((ni, nj), dtype=np.float64)
    for i, species in enumerate(species_names):
        for j, equation in enumerate(reaction_equations):
            reaction = equation.replace(' ', '').replace('=>', '=')
            reactants = reaction.split('=')[0].split('+')
            products = reaction.split('=')[1].split('+')
            if species in reactants:
                sij[i, j] -= reactants.count(species)
            if species in products:
                sij[i, j] += products.count(species)
    return sij


def load_photochem_state(pt_file: pathlib.Path, settings_file: pathlib.Path):
    """Load a saved atmosphere and compute reaction rates at steady state."""
    pc = EvoAtmosphere(str(MECH_FILE), str(settings_file), str(FLUX_FILE), str(pt_file))
    pc.prep_atmosphere(pc.wrk.usol)

    all_species_names = [str(x) for x in pc.dat.species_names]
    species_names = all_species_names[:-2]
    reactions = [str(x).replace('=>', '=') for x in pc.dat.reaction_equations]
    r1, r2 = get_reactant_indices(all_species_names, reactions)

    densities_ext = np.vstack(
        [pc.wrk.densities.astype(F128), pc.wrk.densities[-1, :]]
    )
    if 'M' in all_species_names:
        m_idx = all_species_names.index('M')
        if m_idx >= densities_ext.shape[0]:
            air_density = pc.wrk.density.astype(F128)
            densities_ext = np.vstack([densities_ext, air_density])
        else:
            densities_ext[m_idx, :] = pc.wrk.density.astype(F128)
    densities = densities_ext[: pc.dat.nq, :]

    rx_rates = pc.wrk.rx_rates.astype(F128)
    reaction_rates = np.multiply(rx_rates.T, densities_ext[r1, :])
    reaction_rates = np.multiply(reaction_rates, densities_ext[r2, :])

    return species_names, reactions, densities, reaction_rates


def write_layer_hdf5(
    output_dir: pathlib.Path,
    layer_idx: int,
    species_names: list[str],
    reactions: list[str],
    densities: np.ndarray,
    reaction_rates: np.ndarray,
):
    """Write one altitude layer to Chempath HDF5 format."""
    sij = get_sij(species_names, reactions)
    layer_rates = reaction_rates[:, layer_idx].astype(F128)
    layer_dens = densities[:, layer_idx].astype(F128)

    chemprod = np.dot(sij, layer_rates)
    transport_rates = (-chemprod).astype(F128)
    rainout_rates = np.zeros(len(species_names), dtype=F128)
    dist_flux_rates = np.zeros(len(species_names), dtype=F128)
    error_rates = np.zeros(len(species_names), dtype=F128)

    rainout_reactions = [f'{species}_rainout = {species}' for species in species_names]
    transport_reactions = [f'{species}_transport = {species}' for species in species_names]
    dist_flux_reactions = [f'{species}_distflx = {species}' for species in species_names]
    err_reactions = [f'{species}_err = {species}' for species in species_names]

    all_reactions = (
        reactions
        + rainout_reactions
        + transport_reactions
        + dist_flux_reactions
        + err_reactions
    )
    all_rates = np.concatenate(
        [
            layer_rates,
            rainout_rates,
            transport_rates,
            dist_flux_rates,
            error_rates,
        ]
    )

    nd = np.array([layer_dens, layer_dens], dtype=F128)
    model_time = np.array([0.0, 1.0], dtype=F128)

    layer_dir = output_dir / str(layer_idx)
    layer_dir.mkdir(parents=True, exist_ok=True)
    h5_path = layer_dir / 'photochem_output_0.hdf5'

    with h5py.File(h5_path, 'w') as datafile:
        datafile.create_dataset('rates', all_rates.shape, dtype='f8', data=all_rates)
        datafile.create_dataset('model_time', model_time.shape, dtype='f8', data=model_time)
        datafile.create_dataset('num_densities', nd.shape, dtype='f8', data=nd)
        datafile.create_dataset(
            'reaction_equations',
            (len(all_reactions),),
            dtype=h5py.string_dtype(),
            data=all_reactions,
        )
        datafile.create_dataset(
            'species_names',
            (len(species_names),),
            dtype=h5py.string_dtype(),
            data=species_names,
        )

    return h5_path


def convert_case(
    sza: str,
    pal: str = '100pc',
    output_dir: pathlib.Path | None = None,
    layers: list[int] | None = None,
    all_layers: bool = False,
):
    """Convert one SZA case to Chempath HDF5 files."""
    if pal not in SETTINGS_BY_PAL_SZA:
        raise ValueError(f'PAL must be one of {list(SETTINGS_BY_PAL_SZA)}')
    settings_by_sza = SETTINGS_BY_PAL_SZA[pal]
    if sza not in settings_by_sza:
        raise ValueError(f'SZA must be one of {list(settings_by_sza)} for {pal}')

    pt_file = PHOTOCHEM_ROOT / pal / f'Earth_{pal}_{sza}.txt'
    if not pt_file.exists():
        raise FileNotFoundError(pt_file)

    settings_file = settings_by_sza[sza]
    out = output_dir or (DEFAULT_OUTPUT / f'{pal}_sza_{sza}')
    out.mkdir(parents=True, exist_ok=True)

    species_names, reactions, densities, reaction_rates = load_photochem_state(
        pt_file, settings_file
    )
    nz = densities.shape[1]
    if all_layers:
        layer_list = list(range(nz))
    elif layers is not None:
        layer_list = layers
    else:
        layer_list = list(range(0, nz, 5))

    written = []
    for layer_idx in layer_list:
        if layer_idx < 0 or layer_idx >= nz:
            raise ValueError(f'Layer index {layer_idx} out of range 0..{nz - 1}')
        written.append(
            write_layer_hdf5(
                out, layer_idx, species_names, reactions, densities, reaction_rates
            )
        )

    alt_km = 0.5 + np.arange(nz)
    print(f'Converted {pt_file.name} -> {out}')
    print(f'  layers: {layer_list} (altitudes km: {[round(float(alt_km[i]), 1) for i in layer_list]})')
    print(f'  wrote {len(written)} HDF5 files')
    return written


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--pal',
        default='100pc',
        choices=sorted(SETTINGS_BY_PAL_SZA.keys()),
        help='Present atmospheric level case (default: 100pc)',
    )
    parser.add_argument(
        '--sza',
        required=True,
        help='Solar zenith angle tag matching Earth_<pal>_<sza>.txt',
    )
    parser.add_argument(
        '--output',
        type=pathlib.Path,
        default=None,
        help='Output directory (default: chempath_analysis/chempath_input/sza_<sza>)',
    )
    parser.add_argument(
        '--layers',
        type=int,
        nargs='+',
        default=None,
        help='Altitude layer indices to export (0 = lowest layer)',
    )
    parser.add_argument(
        '--all-layers',
        action='store_true',
        help='Export all 100 layers',
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None):
    args = parse_args(argv)
    convert_case(
        sza=args.sza,
        pal=args.pal,
        output_dir=args.output,
        layers=args.layers,
        all_layers=args.all_layers,
    )


if __name__ == '__main__':
    main()
