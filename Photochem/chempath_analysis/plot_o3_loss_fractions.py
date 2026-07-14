# -*- coding: utf-8 -*-
"""
Plot top O3 loss reactions: fractional contribution vs altitude/pressure.

For each PAL case, selects the top-N loss reactions by column-integrated rate,
bins rates in altitude, applies a light rolling median, then plots each case in
its own panel (filtered by PAL).

Example:
    python plot_o3_loss_fractions.py
    python plot_o3_loss_fractions.py --pal 100pc 10pc 1pc 0.1pc --top 7
"""

from __future__ import annotations

import argparse
import pathlib
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from photochem import EvoAtmosphere

PHOTOCHEM_ROOT = pathlib.Path('/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem')
ANALYSIS_ROOT = PHOTOCHEM_ROOT / 'chempath_analysis'
PLOT_DIR = ANALYSIS_ROOT / 'plots'

MECH = PHOTOCHEM_ROOT / 'Old sims' / 'zahnle_earth.yaml'
FLUX = PHOTOCHEM_ROOT / 'Old sims' / 'input' / 'Sun_0.0Ga.txt'
SZA = '48.2'

SETTINGS_BY_PAL = {
    '100pc': PHOTOCHEM_ROOT / 'Proxima_Centauri' / 'input' / 'settings_100pc.yaml',
    '10pc': PHOTOCHEM_ROOT / 'Proxima_Centauri' / 'input' / 'settings_10pc.yaml',
    '1pc': PHOTOCHEM_ROOT / 'Proxima_Centauri' / 'input' / 'settings_1pc.yaml',
    '0.1pc': PHOTOCHEM_ROOT / 'Proxima_Centauri' / 'input' / 'settings_0.1pc.yaml',
}

REACTION_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
    '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
]

MIN_FRAC = 1e-5
LOSS_FLOOR_FRAC = 1e-3
MAX_ALT_KM = 55.0


def normalize_rx(label: str) -> str:
    return re.sub(r'\s+', ' ', label.replace('=>', '->').replace('=', '->')).strip()


def short_label(reaction: str, max_len: int = 40) -> str:
    text = reaction.replace(' + M', '+M').replace(' -> ', '->')
    if len(text) <= max_len:
        return text
    return text[: max_len - 3] + '...'


def load_loss_budget_from_file(
    pt_file: pathlib.Path,
    scenario: str,
    pal: str | None = None,
) -> pd.DataFrame:
    """Return per-level O3 loss rates from an arbitrary steady-state file."""
    if pal is None:
        pal = scenario
    if not pt_file.exists():
        raise FileNotFoundError(pt_file)

    settings = SETTINGS_BY_PAL.get(pal)
    if settings is None:
        raise ValueError(f'No settings for PAL {pal}')

    pc = EvoAtmosphere(str(MECH), str(settings), str(FLUX), str(pt_file))
    pc.prep_atmosphere(pc.wrk.usol)

    pl = pc.production_and_loss('O3', pc.wrk.usol)
    nz = pl.loss.shape[0]
    loss_rx = [normalize_rx(r) for r in pl.loss_rx]

    rows = []
    for j in range(nz):
        alt_km = 0.5 + j
        press_hpa = pc.wrk.pressure[j] / 1e3
        for i, rx in enumerate(loss_rx):
            rate = float(pl.loss[j, i])
            if rate <= 0:
                continue
            rows.append({
                'pal': pal,
                'scenario': scenario,
                'layer': j,
                'alt_km': alt_km,
                'pressure_hPa': press_hpa,
                'reaction': rx,
                'rate_cm3_s': rate,
            })

    return pd.DataFrame(rows)


def load_loss_budget(pal: str) -> pd.DataFrame:
    """Return per-level O3 loss rates."""
    pt_file = PHOTOCHEM_ROOT / pal / f'Earth_{pal}_{SZA}.txt'
    return load_loss_budget_from_file(pt_file, scenario=pal, pal=pal)


def bin_and_fraction(df: pd.DataFrame, bin_km: float, group_col: str = 'pal') -> pd.DataFrame:
    """Bin rates in altitude, then compute per-bin fractions."""
    work = df.copy()
    work['alt_bin'] = np.floor(work['alt_km'] / bin_km) * bin_km + 0.5 * bin_km
    work['press_bin'] = work.groupby([group_col, 'alt_bin'])['pressure_hPa'].transform('mean')

    binned = (
        work.groupby([group_col, 'alt_bin', 'press_bin', 'reaction'], as_index=False)['rate_cm3_s']
        .sum()
    )
    totals = (
        binned.groupby([group_col, 'alt_bin'], as_index=False)['rate_cm3_s']
        .sum()
        .rename(columns={'rate_cm3_s': 'total_loss_cm3_s'})
    )
    binned = binned.merge(totals, on=[group_col, 'alt_bin'])
    binned['fraction'] = binned['rate_cm3_s'] / binned['total_loss_cm3_s']
    return binned.rename(columns={'alt_bin': 'alt_km'})


def mask_insignificant(
    df: pd.DataFrame,
    floor_frac: float,
    max_alt_km: float,
    group_col: str = 'pal',
) -> pd.DataFrame:
    """Drop bins with negligible column loss or above the plotting ceiling."""
    max_loss = df.groupby(group_col)['total_loss_cm3_s'].transform('max')
    keep = (df['total_loss_cm3_s'] >= floor_frac * max_loss) & (df['alt_km'] <= max_alt_km)
    return df.loc[keep].copy()


def smooth_fractions(df: pd.DataFrame, window: int, group_col: str = 'pal') -> pd.DataFrame:
    """Rolling median of fractions along altitude for each group/reaction."""
    if window <= 1:
        return df.copy()

    parts = []
    for (group, reaction), group_df in df.groupby([group_col, 'reaction']):
        g = group_df.sort_values('alt_km').copy()
        g['fraction'] = (
            g['fraction']
            .rolling(window=window, center=True, min_periods=1)
            .median()
        )
        parts.append(g)
    return pd.concat(parts, ignore_index=True)


def top_reactions(df: pd.DataFrame, n: int) -> list[str]:
    ranked = (
        df.groupby('reaction', as_index=False)['rate_cm3_s']
        .sum()
        .sort_values('rate_cm3_s', ascending=False)
    )
    return ranked.head(n)['reaction'].tolist()


def build_color_map(all_reactions: list[str]) -> dict[str, str]:
    return {rx: REACTION_COLORS[i % len(REACTION_COLORS)] for i, rx in enumerate(all_reactions)}


def frac_axis_formatter(value: float, _pos: int) -> str:
    """Tick labels for log fractional axis."""
    if value >= 0.01:
        return f'{value * 100:g}%'
    return f'{value:.0e}'


def plot_scenario(
    df: pd.DataFrame,
    scenario: str,
    top_rx: list[str],
    color_map: dict[str, str],
    ax_alt,
    ax_press,
    title: str,
    group_col: str = 'pal',
):
    sub = df[(df[group_col] == scenario) & (df['reaction'].isin(top_rx))].copy()

    for rx in top_rx:
        rx_df = sub[sub['reaction'] == rx].sort_values('alt_km')
        if rx_df.empty:
            continue
        frac = rx_df['fraction'].values
        valid = np.isfinite(frac) & (frac >= MIN_FRAC)
        if not np.any(valid):
            continue
        label = short_label(rx)
        color = color_map[rx]
        ax_alt.plot(
            frac[valid],
            rx_df['alt_km'].values[valid],
            color=color,
            lw=2.2,
            label=label,
        )
        ax_press.plot(
            frac[valid],
            rx_df['press_bin'].values[valid],
            color=color,
            lw=2.2,
            label=label,
        )

    for ax in (ax_alt, ax_press):
        ax.set_xscale('log')
        ax.set_xlim(MIN_FRAC, 1.0)
        ax.set_xlabel('Fractional O3 loss')
        ax.xaxis.set_major_formatter(plt.FuncFormatter(frac_axis_formatter))
        ax.grid(True, which='both', alpha=0.3)

    ax_alt.set_title(title)
    ax_alt.set_ylabel('Altitude (km)')
    ax_alt.set_ylim(0, MAX_ALT_KM)

    ax_press.set_title(title)
    ax_press.set_ylabel('Pressure (hPa)')
    ax_press.set_yscale('log')
    ax_press.invert_yaxis()


def plot_case(
    df: pd.DataFrame,
    pal: str,
    top_rx: list[str],
    color_map: dict[str, str],
    ax_alt,
    ax_press,
):
    plot_scenario(
        df,
        pal,
        top_rx,
        color_map,
        ax_alt,
        ax_press,
        title=f'{pal} - 48.2 deg SZA',
        group_col='pal',
    )


def make_plots(
    pals: list[str],
    top_n: int,
    output_dir: pathlib.Path,
    bin_km: float = 3.0,
    smooth_window: int = 3,
    loss_floor_frac: float = LOSS_FLOOR_FRAC,
    max_alt_km: float = MAX_ALT_KM,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    all_data = []
    top_by_pal: dict[str, list[str]] = {}
    for pal in pals:
        print(f'Loading {pal} ...')
        raw = load_loss_budget(pal)
        binned = bin_and_fraction(raw, bin_km=bin_km)
        binned = mask_insignificant(binned, loss_floor_frac, max_alt_km)
        binned = smooth_fractions(binned, smooth_window)
        all_data.append(binned)
        top_by_pal[pal] = top_reactions(binned, top_n)
        print(f'  top {top_n}:')
        for rx in top_by_pal[pal]:
            print(f'    - {rx}')

    combined = pd.concat(all_data, ignore_index=True)
    combined['fraction_pct'] = combined['fraction'] * 100
    combined.to_csv(output_dir / 'o3_loss_fraction_by_altitude.csv', index=False)

    all_top_rx = sorted({rx for rxs in top_by_pal.values() for rx in rxs})
    color_map = build_color_map(all_top_rx)

    n = len(pals)
    fig_alt, axes_alt = plt.subplots(1, n, figsize=(4.8 * n, 8), sharey=True)
    fig_press, axes_press = plt.subplots(1, n, figsize=(4.8 * n, 8), sharey=True)
    if n == 1:
        axes_alt = [axes_alt]
        axes_press = [axes_press]

    for ax_alt, ax_press, pal in zip(axes_alt, axes_press, pals):
        plot_case(combined, pal, top_by_pal[pal], color_map, ax_alt, ax_press)

    handles, labels = axes_alt[0].get_legend_handles_labels()
    for ax in axes_alt[1:]:
        h, l = ax.get_legend_handles_labels()
        for handle, label in zip(h, l):
            if label not in labels:
                handles.append(handle)
                labels.append(label)

    fig_alt.legend(
        handles,
        labels,
        loc='lower center',
        bbox_to_anchor=(0.5, -0.05),
        ncol=3,
        fontsize=7.5,
        frameon=False,
    )
    fig_press.legend(
        handles,
        labels,
        loc='lower center',
        bbox_to_anchor=(0.5, -0.05),
        ncol=3,
        fontsize=7.5,
        frameon=False,
    )

    note = (
        f'Top {top_n} O3 loss reactions | {bin_km:g} km bins | '
        f'{smooth_window}-bin median smooth | loss > {loss_floor_frac:.0e} of peak | '
        f'z <= {max_alt_km:g} km'
    )
    fig_alt.suptitle(f'{note}\nvs altitude', fontsize=11, y=1.02)
    fig_press.suptitle(f'{note}\nvs pressure', fontsize=11, y=1.02)

    fig_alt.tight_layout()
    fig_press.tight_layout()

    alt_path = output_dir / f'o3_loss_top{top_n}_vs_altitude.png'
    press_path = output_dir / f'o3_loss_top{top_n}_vs_pressure.png'
    fig_alt.savefig(alt_path, dpi=180, bbox_inches='tight')
    fig_press.savefig(press_path, dpi=180, bbox_inches='tight')
    plt.close(fig_alt)
    plt.close(fig_press)

    print(f'Wrote {alt_path}')
    print(f'Wrote {press_path}')


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--pal',
        nargs='+',
        default=['100pc', '10pc', '1pc', '0.1pc'],
        choices=sorted(SETTINGS_BY_PAL.keys()),
    )
    parser.add_argument('--top', type=int, default=7, help='Number of top loss reactions')
    parser.add_argument('--bin-km', type=float, default=3.0, help='Altitude bin width (km)')
    parser.add_argument(
        '--smooth-window',
        type=int,
        default=3,
        help='Rolling median window in altitude bins (1 = no smoothing)',
    )
    parser.add_argument(
        '--loss-floor-frac',
        type=float,
        default=LOSS_FLOOR_FRAC,
        help='Mask bins with loss below this fraction of the column peak',
    )
    parser.add_argument(
        '--max-alt-km',
        type=float,
        default=MAX_ALT_KM,
        help='Maximum altitude shown (km)',
    )
    parser.add_argument('--output', type=pathlib.Path, default=PLOT_DIR)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None):
    args = parse_args(argv)
    make_plots(
        args.pal,
        args.top,
        args.output,
        bin_km=args.bin_km,
        smooth_window=args.smooth_window,
        loss_floor_frac=args.loss_floor_frac,
        max_alt_km=args.max_alt_km,
    )


if __name__ == '__main__':
    main()
