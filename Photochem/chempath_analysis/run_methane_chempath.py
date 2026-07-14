# -*- coding: utf-8 -*-
"""
Convert Methane Perturbation atmospheres and run Chempath O3 pathway analysis.

Default: 1pc at 0.1x, 1x, 10x methane, troposphere + stratosphere layers.

Example:
    python run_methane_chempath.py
    python run_methane_chempath.py --pal 100pc --methane 0.1x 1x 10x --layers 14 23
    python run_methane_chempath.py --plot
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

PHOTOCHEM_ROOT = pathlib.Path('/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem')
ANALYSIS_ROOT = PHOTOCHEM_ROOT / 'chempath_analysis'
METHANE_DIR = PHOTOCHEM_ROOT / 'Methane_Perturbations'
SZA = '48.2'

SETTINGS_BY_PAL = {
    '100pc': PHOTOCHEM_ROOT / 'Proxima_Centauri' / 'input' / 'settings_100pc.yaml',
    '10pc': PHOTOCHEM_ROOT / 'Proxima_Centauri' / 'input' / 'settings_10pc.yaml',
    '1pc': PHOTOCHEM_ROOT / 'Proxima_Centauri' / 'input' / 'settings_1pc.yaml',
    '0.1pc': PHOTOCHEM_ROOT / 'Proxima_Centauri' / 'input' / 'settings_0.1pc.yaml',
}

DEFAULT_METHANE = ['0.1x', '1x', '10x']
DEFAULT_LAYERS = [14, 23]
DEFAULT_PAL = ['1pc']


def convert_methane_case(
    pal: str,
    methane: str,
    layers: list[int],
) -> pathlib.Path:
    pt_file = METHANE_DIR / f'Earth_{pal}_{SZA}_{methane}_methane.txt'
    out = ANALYSIS_ROOT / 'chempath_input' / f'methane_{pal}_{methane}'
    settings = SETTINGS_BY_PAL[pal]

    sys.path.insert(0, str(ANALYSIS_ROOT))
    from steadystate_to_chempath import load_photochem_state, write_layer_hdf5

    if not pt_file.exists():
        raise FileNotFoundError(pt_file)

    out.mkdir(parents=True, exist_ok=True)
    species_names, reactions, densities, reaction_rates = load_photochem_state(
        pt_file, settings
    )
    nz = densities.shape[1]
    for layer_idx in layers:
        if layer_idx < 0 or layer_idx >= nz:
            raise ValueError(f'Layer {layer_idx} out of range')
        write_layer_hdf5(
            out, layer_idx, species_names, reactions, densities, reaction_rates
        )
    print(f'Converted {pt_file.name} -> {out}')
    return out


def run_pathways(
    input_root: pathlib.Path,
    output_root: pathlib.Path,
    layers: list[int],
    timeout: int,
    python: str,
):
    cmd = [
        python,
        str(ANALYSIS_ROOT / 'run_pathways.py'),
        '--pal', '100pc',
        '--sza', SZA,
        '--input', str(input_root),
        '--output', str(output_root),
        '--layers', *[str(x) for x in layers],
        '--species', 'O3',
        '--timeout', str(timeout),
    ]
    subprocess.run(cmd, check=True)


def plot_pathways(output_root: pathlib.Path, layers: list[int], label: str, python: str):
    cmd = [
        python,
        str(ANALYSIS_ROOT / 'plot_chempath_pathways.py'),
        '--input', str(output_root),
        '--output', str(ANALYSIS_ROOT / 'plots' / 'chempath_pathways' / output_root.name),
        '--label', label,
        '--layers', *[str(x) for x in layers],
        '--process', 'loss',
        '--top', '10',
    ]
    subprocess.run(cmd, check=True)


def main(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pal', nargs='+', default=DEFAULT_PAL, choices=sorted(SETTINGS_BY_PAL))
    parser.add_argument('--methane', nargs='+', default=DEFAULT_METHANE)
    parser.add_argument('--layers', type=int, nargs='+', default=DEFAULT_LAYERS)
    parser.add_argument('--timeout', type=int, default=600)
    parser.add_argument('--plot', action='store_true', help='Also generate pathway bar charts')
    parser.add_argument('--python', default=sys.executable)
    args = parser.parse_args(argv)

    for pal in args.pal:
        for methane in args.methane:
            tag = f'methane_{pal}_{methane}'
            print(f'=== {tag} ===')
            input_root = convert_methane_case(pal, methane, args.layers)
            output_root = ANALYSIS_ROOT / 'pathways' / tag
            run_pathways(input_root, output_root, args.layers, args.timeout, args.python)
            if args.plot:
                plot_pathways(output_root, args.layers, f'{pal} {methane} methane', args.python)


if __name__ == '__main__':
    main()
