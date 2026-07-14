# -*- coding: utf-8 -*-
"""
O3 reaction budget and Chempath pathway layer selection for 48.2 deg SZA.

Produces per PAL case (100pc / 1pc):
  - o3_reactions_by_altitude.csv   (every O3-affecting reaction at each level)
  - o3_top_reactions_trop.csv      (ranked reactions in troposphere)
  - o3_top_reactions_strat.csv     (ranked reactions in stratosphere)
  - o3_summary.md                  (human-readable summary)
  - pathway_layers.csv             (representative layers for Chempath)

Troposphere: altitude < 15 km
Stratosphere: 15 km <= altitude < 55 km
"""

from __future__ import annotations

import argparse
import pathlib
import re

import numpy as np
import pandas as pd
from photochem import EvoAtmosphere

PHOTOCHEM_ROOT = pathlib.Path('/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem')
ANALYSIS_ROOT = PHOTOCHEM_ROOT / 'chempath_analysis'

MECH = PHOTOCHEM_ROOT / 'Old sims' / 'zahnle_earth.yaml'
FLUX = PHOTOCHEM_ROOT / 'Old sims' / 'input' / 'Sun_0.0Ga.txt'

SETTINGS_BY_PAL = {
    '100pc': PHOTOCHEM_ROOT / 'Proxima_Centauri' / 'input' / 'settings_100pc.yaml',
    '1pc': PHOTOCHEM_ROOT / 'Proxima_Centauri' / 'input' / 'settings_1pc.yaml',
}

TROP_TOP_KM = 15.0
STRAT_TOP_KM = 55.0


def layer_to_alt_km(layer_idx: int) -> float:
    return 0.5 + layer_idx


def normalize_rx(label: str) -> str:
    return re.sub(r'\s+', ' ', label.replace('=>', '->').replace('=', '->')).strip()


def region(alt_km: float) -> str:
    if alt_km < TROP_TOP_KM:
        return 'troposphere'
    if alt_km < STRAT_TOP_KM:
        return 'stratosphere'
    return 'mesosphere'


def build_reaction_budget(pc: EvoAtmosphere) -> pd.DataFrame:
    """Per-altitude O3 production and loss for every reaction."""
    pl = pc.production_and_loss('O3', pc.wrk.usol)
    nz = pl.production.shape[0]
    rows = []

    prod_rx = [normalize_rx(r) for r in pl.production_rx]
    loss_rx = [normalize_rx(r) for r in pl.loss_rx]

    for j in range(nz):
        alt = layer_to_alt_km(j)
        press_hpa = pc.wrk.pressure[j] / 1e3
        for i, rx in enumerate(prod_rx):
            rate = pl.production[j, i]
            if rate > 0:
                rows.append({
                    'alt_km': alt,
                    'layer': j,
                    'pressure_hPa': press_hpa,
                    'region': region(alt),
                    'process': 'production',
                    'reaction': rx,
                    'rate_cm3_s': rate,
                })
        for i, rx in enumerate(loss_rx):
            rate = pl.loss[j, i]
            if rate > 0:
                rows.append({
                    'alt_km': alt,
                    'layer': j,
                    'pressure_hPa': press_hpa,
                    'region': region(alt),
                    'process': 'loss',
                    'reaction': rx,
                    'rate_cm3_s': rate,
                })

    return pd.DataFrame(rows)


def rank_regional_reactions(df: pd.DataFrame, reg: str, top_n: int = 25) -> pd.DataFrame:
    """Sum reaction rates over altitudes in a region, return top contributors."""
    sub = df[df['region'] == reg].copy()
    grouped = (
        sub.groupby(['process', 'reaction'], as_index=False)['rate_cm3_s']
        .sum()
        .sort_values('rate_cm3_s', ascending=False)
    )
    total_prod = grouped.loc[grouped['process'] == 'production', 'rate_cm3_s'].sum()
    total_loss = grouped.loc[grouped['process'] == 'loss', 'rate_cm3_s'].sum()

    grouped['fraction_of_regional_prod_or_loss'] = grouped.apply(
        lambda r: r['rate_cm3_s'] / total_prod if r['process'] == 'production' else r['rate_cm3_s'] / total_loss,
        axis=1,
    )
    grouped['region'] = reg
    return grouped.head(top_n)


def representative_altitudes(df: pd.DataFrame, reg: str, n: int = 3) -> list[int]:
    """Pick layers with highest total O3 turnover in a region."""
    sub = df[df['region'] == reg]
    turnover = (
        sub.groupby('layer')['rate_cm3_s']
        .sum()
        .sort_values(ascending=False)
    )
    return [int(x) for x in turnover.head(n).index.tolist()]


def run_case(pal: str) -> list[int]:
    if pal not in SETTINGS_BY_PAL:
        raise ValueError(f'PAL must be one of {list(SETTINGS_BY_PAL)}')

    pt_file = PHOTOCHEM_ROOT / pal / f'Earth_{pal}_48.2.txt'
    if not pt_file.exists():
        raise FileNotFoundError(pt_file)

    out = ANALYSIS_ROOT / f'o3_48.2_{pal}'
    out.mkdir(parents=True, exist_ok=True)

    print(f'Loading {pt_file.name} and computing O3 production/loss ...')
    pc = EvoAtmosphere(str(MECH), str(SETTINGS_BY_PAL[pal]), str(FLUX), str(pt_file))
    pc.prep_atmosphere(pc.wrk.usol)

    budget = build_reaction_budget(pc)
    budget.to_csv(out / 'o3_reactions_by_altitude.csv', index=False)

    trop_rank = rank_regional_reactions(budget, 'troposphere')
    strat_rank = rank_regional_reactions(budget, 'stratosphere')
    trop_rank.to_csv(out / 'o3_top_reactions_trop.csv', index=False)
    strat_rank.to_csv(out / 'o3_top_reactions_strat.csv', index=False)

    trop_layers = representative_altitudes(budget, 'troposphere', 3)
    strat_layers = representative_altitudes(budget, 'stratosphere', 3)

    summary_lines = [
        f'# O3 analysis - 48.2 deg SZA ({pal} Early Earth)',
        '',
        '## Regional definition',
        f'- Troposphere: altitude < {TROP_TOP_KM} km',
        f'- Stratosphere: {TROP_TOP_KM}-{STRAT_TOP_KM} km',
        '',
        '## Representative altitudes for pathway analysis',
        f'- Troposphere layers {trop_layers}: alt km {[layer_to_alt_km(x) for x in trop_layers]}',
        f'- Stratosphere layers {strat_layers}: alt km {[layer_to_alt_km(x) for x in strat_layers]}',
        '',
        '## Top O3 production reactions (troposphere, column-integrated)',
    ]
    for _, row in trop_rank[trop_rank['process'] == 'production'].head(10).iterrows():
        summary_lines.append(
            f"- {row['fraction_of_regional_prod_or_loss']*100:5.1f}%  {row['reaction']}"
        )
    summary_lines.append('')
    summary_lines.append('## Top O3 loss reactions (troposphere, column-integrated)')
    for _, row in trop_rank[trop_rank['process'] == 'loss'].head(10).iterrows():
        summary_lines.append(
            f"- {row['fraction_of_regional_prod_or_loss']*100:5.1f}%  {row['reaction']}"
        )
    summary_lines.append('')
    summary_lines.append('## Top O3 production reactions (stratosphere, column-integrated)')
    for _, row in strat_rank[strat_rank['process'] == 'production'].head(10).iterrows():
        summary_lines.append(
            f"- {row['fraction_of_regional_prod_or_loss']*100:5.1f}%  {row['reaction']}"
        )
    summary_lines.append('')
    summary_lines.append('## Top O3 loss reactions (stratosphere, column-integrated)')
    for _, row in strat_rank[strat_rank['process'] == 'loss'].head(10).iterrows():
        summary_lines.append(
            f"- {row['fraction_of_regional_prod_or_loss']*100:5.1f}%  {row['reaction']}"
        )

    (out / 'o3_summary.md').write_text('\n'.join(summary_lines) + '\n')

    all_layers = sorted(set(trop_layers + strat_layers))
    pd.Series(all_layers, name='layer').to_csv(out / 'pathway_layers.csv', index=False)

    print(f'Wrote results to {out}')
    print('Top troposphere O3 production:')
    print(trop_rank[trop_rank['process'] == 'production'].head(5).to_string(index=False))
    print('Top stratosphere O3 production:')
    print(strat_rank[strat_rank['process'] == 'production'].head(5).to_string(index=False))
    return all_layers


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--pal',
        nargs='+',
        default=['100pc', '1pc'],
        choices=sorted(SETTINGS_BY_PAL.keys()),
        help='PAL cases to analyse (default: 100pc 1pc)',
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None):
    args = parse_args(argv)
    for pal in args.pal:
        run_case(pal)


if __name__ == '__main__':
    main()
