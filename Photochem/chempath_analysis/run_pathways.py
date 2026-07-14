"""
Run Chempath pathway analysis on converted Photochem steady-state input.

Example:
    python run_pathways.py --sza 48.2 --layers 20 30
    python run_pathways.py --sza 48.2 --species O3 --timeout 600
"""

from __future__ import annotations

import argparse
import pathlib

import numpy as np
import pandas as pd
from chempath import Chempath
ANALYSIS_ROOT = pathlib.Path('/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/chempath_analysis')
DEFAULT_INPUT = ANALYSIS_ROOT / 'chempath_input'
DEFAULT_OUTPUT = ANALYSIS_ROOT / 'pathways'

DEFAULT_IGNORED_SB = [
    'O2', 'H2O', 'O3', 'CH4', 'HV', 'N2', 'CO2', 'CO', 'H2', 'M',
]
DEFAULT_SPECIES = ['O3', 'O2', 'CH4', 'CO', 'H2']


def get_fmin(chempath):
    species_idxs = np.array(
        [chempath.species_list.index(sp) for sp in DEFAULT_SPECIES if sp in chempath.species_list]
    )
    possij = np.multiply(chempath.sij, chempath.sij > 0)
    production_rates = np.dot(possij, chempath.fk)[species_idxs]
    return min(production_rates) / 1e4


def run_layer(
    input_root: pathlib.Path,
    output_root: pathlib.Path,
    layer_idx: int,
    species_list: list[str],
    timeout: int,
):
    h5_path = input_root / str(layer_idx) / 'photochem_output_0.hdf5'
    if not h5_path.exists():
        raise FileNotFoundError(h5_path)

    chempath = Chempath(
        h5py_path=str(h5_path),
        transport_species=True,
        delete_error_reactions=False,
    )
    chempath.ignored_sb += DEFAULT_IGNORED_SB
    chempath.f_min = get_fmin(chempath)

    status = chempath.find_all_pathways(timeout=timeout)
    chempath.del_error_reactions()
    if status == 'timed out':
        print(f'  layer {layer_idx}: timed out')
        return None, None

    alt_km = 0.5 + layer_idx
    prod_dfs = []
    loss_dfs = []
    for species in species_list:
        if species not in chempath.species_list:
            continue
        prod = chempath.get_pathways_contributions(species, on='production')
        loss = chempath.get_pathways_contributions(species, on='loss')
        prod['species'] = species
        loss['species'] = species
        prod['alt_km'] = alt_km
        loss['alt_km'] = alt_km
        prod_dfs.append(prod)
        loss_dfs.append(loss)

    if not prod_dfs:
        return None, None

    prod_df = pd.concat(prod_dfs, ignore_index=True)
    loss_df = pd.concat(loss_dfs, ignore_index=True)

    out_dir = output_root / str(layer_idx)
    out_dir.mkdir(parents=True, exist_ok=True)
    prod_df.to_csv(out_dir / 'prod_pathways.csv', index=False)
    loss_df.to_csv(out_dir / 'loss_pathways.csv', index=False)
    print(f'  layer {layer_idx} ({alt_km:.1f} km): {len(prod_df)} prod rows, {len(loss_df)} loss rows')
    return prod_df, loss_df


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pal', default='100pc', choices=['100pc', '1pc'])
    parser.add_argument('--sza', required=True)
    parser.add_argument('--layers', type=int, nargs='+', default=[20, 30, 40])
    parser.add_argument('--species', nargs='+', default=['O3'])
    parser.add_argument('--timeout', type=int, default=3600, help='Seconds per layer')
    parser.add_argument('--input', type=pathlib.Path, default=None)
    parser.add_argument('--output', type=pathlib.Path, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None):
    args = parse_args(argv)
    input_root = args.input or (DEFAULT_INPUT / f'{args.pal}_sza_{args.sza}')
    output_root = args.output or (DEFAULT_OUTPUT / f'{args.pal}_sza_{args.sza}')

    if not input_root.exists():
        raise FileNotFoundError(
            f'No Chempath input at {input_root}. Run steadystate_to_chempath.py first.'
        )

    print(f'Chempath pathways: sza={args.sza}, layers={args.layers}, species={args.species}')
    all_prod = []
    all_loss = []
    for layer_idx in args.layers:
        prod_df, loss_df = run_layer(
            input_root, output_root, layer_idx, args.species, args.timeout
        )
        if prod_df is not None:
            all_prod.append(prod_df)
        if loss_df is not None:
            all_loss.append(loss_df)

    if all_prod:
        pd.concat(all_prod, ignore_index=True).to_csv(
            output_root / 'prod_pathways_summary.csv', index=False
        )
    if all_loss:
        pd.concat(all_loss, ignore_index=True).to_csv(
            output_root / 'loss_pathways_summary.csv', index=False
        )
    print(f'Done. Results in {output_root}')


if __name__ == '__main__':
    main()
