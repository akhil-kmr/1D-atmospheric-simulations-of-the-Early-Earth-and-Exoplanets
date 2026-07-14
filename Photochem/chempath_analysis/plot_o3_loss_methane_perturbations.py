# -*- coding: utf-8 -*-
"""
O3 loss pathway plots for Methane Perturbation simulations.

One three-panel figure per oxygen level (100pc, 10pc, 1pc, 0.1pc), comparing:
  - minimum methane (0.1x)
  - original methane (1x)
  - maximum methane (10x)

Example:
    python plot_o3_loss_methane_perturbations.py
"""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.pyplot as plt
import pandas as pd

from plot_o3_loss_fractions import (
    MAX_ALT_KM,
    MIN_FRAC,
    PHOTOCHEM_ROOT,
    SETTINGS_BY_PAL,
    bin_and_fraction,
    build_color_map,
    load_loss_budget_from_file,
    mask_insignificant,
    plot_scenario,
    smooth_fractions,
    top_reactions,
)

METHANE_DIR = PHOTOCHEM_ROOT / 'Methane_Perturbations'
OUTPUT_DIR = PHOTOCHEM_ROOT / 'chempath_analysis' / 'plots' / 'methane_perturbations'
SZA = '48.2'

METHANE_CASES = [
    ('0.1x', 'Min methane (0.1x)'),
    ('1x', 'Original (1x)'),
    ('10x', 'Max methane (10x)'),
]

PALS = ['100pc', '10pc', '1pc', '0.1pc']


def methane_pt_file(pal: str, methane_tag: str) -> pathlib.Path:
    return METHANE_DIR / f'Earth_{pal}_{SZA}_{methane_tag}_methane.txt'


def load_pal_methane_cases(pal: str) -> pd.DataFrame:
    frames = []
    for methane_tag, _label in METHANE_CASES:
        pt_file = methane_pt_file(pal, methane_tag)
        scenario = f'{pal}_{methane_tag}'
        frames.append(load_loss_budget_from_file(pt_file, scenario=scenario, pal=pal))
    out = pd.concat(frames, ignore_index=True)
    out['methane'] = out['scenario'].str.replace(f'{pal}_', '', regex=False)
    return out


def prepare_binned_data(
    raw: pd.DataFrame,
    bin_km: float,
    smooth_window: int,
    loss_floor_frac: float,
    max_alt_km: float,
) -> pd.DataFrame:
    binned = bin_and_fraction(raw, bin_km=bin_km, group_col='scenario')
    binned = mask_insignificant(
        binned, loss_floor_frac, max_alt_km, group_col='scenario'
    )
    binned = smooth_fractions(binned, smooth_window, group_col='scenario')
    binned['fraction_pct'] = binned['fraction'] * 100
    return binned


def column_ranked_loss(df: pd.DataFrame, scenario: str) -> pd.DataFrame:
  sub = df[df['scenario'] == scenario]
  ranked = (
      sub.groupby('reaction', as_index=False)['rate_cm3_s']
      .sum()
      .sort_values('rate_cm3_s', ascending=False)
  )
  total = ranked['rate_cm3_s'].sum()
  ranked['fraction_of_total_loss'] = ranked['rate_cm3_s'] / total
  ranked['scenario'] = scenario
  return ranked


def write_summary(pal: str, binned: pd.DataFrame, top_n: int, output_dir: pathlib.Path):
    rows = []
    for methane_tag, label in METHANE_CASES:
        scenario = f'{pal}_{methane_tag}'
        ranked = column_ranked_loss(binned, scenario)
        for rank, (_, row) in enumerate(ranked.head(top_n).iterrows(), start=1):
            rows.append({
                'pal': pal,
                'methane': methane_tag,
                'methane_label': label,
                'rank': rank,
                'reaction': row['reaction'],
                'column_loss_rate': row['rate_cm3_s'],
                'fraction_of_total_loss': row['fraction_of_total_loss'],
            })
    summary = pd.DataFrame(rows)
    summary.to_csv(output_dir / f'o3_loss_ranking_{pal}.csv', index=False)

    lines = [f'# O3 loss pathways - {pal} methane perturbations (48.2 deg SZA)', '']
    for methane_tag, label in METHANE_CASES:
        scenario = f'{pal}_{methane_tag}'
        lines.append(f'## {label}')
        ranked = column_ranked_loss(binned, scenario)
        for _, row in ranked.head(top_n).iterrows():
            pct = row['fraction_of_total_loss'] * 100
            lines.append(f'- {pct:5.1f}%  {row["reaction"]}')
        lines.append('')

    min_s = f'{pal}_0.1x'
    max_s = f'{pal}_10x'
    min_rank = column_ranked_loss(binned, min_s).set_index('reaction')['fraction_of_total_loss']
    max_rank = column_ranked_loss(binned, max_s).set_index('reaction')['fraction_of_total_loss']
    all_rx = min_rank.index.union(max_rank.index)
    diff = (max_rank.reindex(all_rx, fill_value=0) - min_rank.reindex(all_rx, fill_value=0)).sort_values(
        ascending=False
    )
    lines.append('## Biggest gain in importance (10x minus 0.1x methane, column-integrated)')
    for rx, delta in diff.head(8).items():
        if delta <= 0:
            break
        lines.append(f'- +{delta*100:4.1f} pp  {rx}')
    lines.append('')
    lines.append('## Biggest loss in importance (10x minus 0.1x methane)')
    for rx, delta in diff.tail(8).sort_values().items():
        if delta >= 0:
            break
        lines.append(f'- {delta*100:4.1f} pp  {rx}')
    lines.append('')

    (output_dir / f'o3_loss_summary_{pal}.md').write_text('\n'.join(lines) + '\n')


def plot_pal_figure(
    pal: str,
    binned: pd.DataFrame,
    top_n: int,
    output_dir: pathlib.Path,
):
    top_by_case = {}
    for methane_tag, label in METHANE_CASES:
        scenario = f'{pal}_{methane_tag}'
        case_df = binned[binned['scenario'] == scenario]
        top_by_case[methane_tag] = top_reactions(case_df, top_n)

    all_top_rx = sorted({rx for rxs in top_by_case.values() for rx in rxs})
    color_map = build_color_map(all_top_rx)

    fig_alt, axes_alt = plt.subplots(1, 3, figsize=(14, 8), sharey=True)
    fig_press, axes_press = plt.subplots(1, 3, figsize=(14, 8), sharey=True)

    for ax_alt, ax_press, (methane_tag, label) in zip(axes_alt, axes_press, METHANE_CASES):
        scenario = f'{pal}_{methane_tag}'
        plot_scenario(
            binned,
            scenario,
            top_by_case[methane_tag],
            color_map,
            ax_alt,
            ax_press,
            title=label,
            group_col='scenario',
        )

    handles, labels = axes_alt[0].get_legend_handles_labels()
    for ax in axes_alt[1:]:
        h, lab = ax.get_legend_handles_labels()
        for handle, label in zip(h, lab):
            if label not in labels:
                handles.append(handle)
                labels.append(label)

    for fig, axes, suffix in (
        (fig_alt, axes_alt, 'altitude'),
        (fig_press, axes_press, 'pressure'),
    ):
        fig.legend(
            handles,
            labels,
            loc='lower center',
            bbox_to_anchor=(0.5, -0.04),
            ncol=3,
            fontsize=7.5,
            frameon=False,
        )
        fig.suptitle(
            f'{pal} O3 loss pathways | methane perturbations | 48.2 deg SZA | top {top_n} reactions',
            fontsize=12,
            y=1.01,
        )
        fig.tight_layout()
        out = output_dir / f'o3_loss_{pal}_methane_vs_{suffix}.png'
        fig.savefig(out, dpi=180, bbox_inches='tight')
        plt.close(fig)
        print(f'Wrote {out}')


def make_plots(
    pals: list[str],
    top_n: int,
    output_dir: pathlib.Path,
    bin_km: float = 3.0,
    smooth_window: int = 3,
    loss_floor_frac: float = 1e-4,
    max_alt_km: float = MAX_ALT_KM,
):
    output_dir.mkdir(parents=True, exist_ok=True)
    all_binned = []

    for pal in pals:
        print(f'Processing {pal} ...')
        raw = load_pal_methane_cases(pal)
        binned = prepare_binned_data(
            raw, bin_km, smooth_window, loss_floor_frac, max_alt_km
        )
        all_binned.append(binned)
        write_summary(pal, binned, top_n, output_dir)
        plot_pal_figure(pal, binned, top_n, output_dir)

    combined = pd.concat(all_binned, ignore_index=True)
    combined.to_csv(output_dir / 'o3_loss_fraction_by_altitude.csv', index=False)
    print(f'Wrote {output_dir / "o3_loss_fraction_by_altitude.csv"}')


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--pal',
        nargs='+',
        default=PALS,
        choices=sorted(SETTINGS_BY_PAL.keys()),
    )
    parser.add_argument('--top', type=int, default=7)
    parser.add_argument('--bin-km', type=float, default=3.0)
    parser.add_argument('--smooth-window', type=int, default=3)
    parser.add_argument('--loss-floor-frac', type=float, default=1e-4)
    parser.add_argument('--max-alt-km', type=float, default=MAX_ALT_KM)
    parser.add_argument('--output', type=pathlib.Path, default=OUTPUT_DIR)
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
