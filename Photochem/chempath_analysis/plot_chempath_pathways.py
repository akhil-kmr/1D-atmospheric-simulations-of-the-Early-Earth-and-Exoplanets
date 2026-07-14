# -*- coding: utf-8 -*-
"""
Plot Chempath O3 pathway contributions from run_pathways.py CSV output.

Chempath itself only writes tables (CSV). This script makes horizontal bar
charts of the top chemical pathways, excluding transport budget-closure terms.

Example:
    python plot_chempath_pathways.py --input pathways/100pc_sza_48.2 --layers 14 23
    python plot_chempath_pathways.py --input pathways/methane_1pc_1x --process loss
"""

from __future__ import annotations

import argparse
import pathlib
import textwrap

import matplotlib.pyplot as plt
import pandas as pd

ANALYSIS_ROOT = pathlib.Path('/Users/gregcooke/1D-Simulations-of-the-Early-Earth/Photochem/chempath_analysis')
DEFAULT_OUTPUT = ANALYSIS_ROOT / 'plots' / 'chempath_pathways'


def pathway_net_label(pathway: str, width: int = 55) -> str:
    """Short label: first O3 chemistry step, else Net line."""
    steps = [
        ln.strip()
        for ln in pathway.split('\n')
        if ln.strip() and not ln.startswith('Net:')
    ]
    chem = [s for s in steps if 'transport' not in s.lower() and 'O3' in s]
    if chem:
        label = chem[0]
    else:
        net = next(
            (ln.replace('Net:', '').strip() for ln in pathway.split('\n') if ln.startswith('Net:')),
            None,
        )
        label = net or (steps[0] if steps else pathway)
    return textwrap.fill(label, width=width)


def pathway_steps(pathway: str) -> str:
    """Full pathway as arrow-separated steps for text output."""
    steps = [ln.strip() for ln in pathway.split('\n') if ln.strip() and not ln.startswith('Net:')]
    net = next((ln.replace('Net:', '').strip() for ln in pathway.split('\n') if ln.startswith('Net:')), '')
    body = ' -> '.join(steps)
    if net:
        return f'{body}  ||  Net: {net}'
    return body


def is_budget_pathway(row: pd.Series) -> bool:
    """Drop transport closure and bookkeeping rows."""
    pid = str(row.get('pathway_id', ''))
    pathway = str(row.get('pathway', ''))
    if pid in ('del', 'err'):
        return True
    if 'O3_transport -> O3' in pathway.replace(' ', ''):
        return True
    if row.get('contribution', 0) > 0.99 and 'transport' in pathway.lower():
        return True
    return False


def load_pathways(csv_path: pathlib.Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if df.empty:
        return df
    df = df[~df.apply(is_budget_pathway, axis=1)].copy()
    df = df.sort_values('contribution', ascending=False).reset_index(drop=True)
    return df


def plot_top_pathways(
    df: pd.DataFrame,
    title: str,
    output_path: pathlib.Path,
    top_n: int = 10,
    min_contribution: float = 1e-8,
):
    sub = df[df['contribution'] >= min_contribution].head(top_n).copy()
    if sub.empty:
        print(f'  no chemical pathways above threshold for {output_path.name}')
        return

    sub = sub.sort_values('contribution', ascending=True)
    labels = [pathway_net_label(p) for p in sub['pathway']]
    pct = sub['contribution'].values * 100

    fig_h = max(4, 0.45 * len(sub) + 1.5)
    fig, ax = plt.subplots(figsize=(10, fig_h))
    ax.barh(range(len(sub)), pct, color='#4472C4', edgecolor='white')
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel('Fraction of O3 loss (%)')
    ax.set_xscale('log')
    ax.set_xlim(max(min_contribution * 100, pct.min() * 0.5), 100)
    ax.set_title(title, fontsize=11)
    ax.grid(True, axis='x', alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f'Wrote {output_path}')


def write_text_summary(
    df: pd.DataFrame,
    output_path: pathlib.Path,
    title: str,
    top_n: int = 15,
):
    lines = [f'# {title}', '']
    for i, row in df.head(top_n).iterrows():
        pct = row['contribution'] * 100
        lines.append(f'## {i + 1}. {pct:.4g}%')
        for line in str(row['pathway']).split('\n'):
            lines.append(f'    {line}')
        lines.append('')
    output_path.write_text('\n'.join(lines) + '\n')
    print(f'Wrote {output_path}')


def plot_case(
    input_root: pathlib.Path,
    output_root: pathlib.Path,
    layers: list[int],
    process: str,
    top_n: int,
    case_label: str,
):
    output_root.mkdir(parents=True, exist_ok=True)
    csv_name = 'loss_pathways.csv' if process == 'loss' else 'prod_pathways.csv'

    for layer in layers:
        csv_path = input_root / str(layer) / csv_name
        if not csv_path.exists():
            print(f'Missing {csv_path}')
            continue
        alt = 0.5 + layer
        df = load_pathways(csv_path)
        tag = f'{case_label}_layer{layer}_{alt:.1f}km_{process}'
        plot_top_pathways(
            df,
            title=f'{case_label} | O3 {process} pathways | {alt:.1f} km',
            output_path=output_root / f'{tag}.png',
            top_n=top_n,
        )
        write_text_summary(
            df,
            output_path=output_root / f'{tag}.md',
            title=f'{case_label} | O3 {process} | {alt:.1f} km',
            top_n=top_n,
        )


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=pathlib.Path, required=True)
    parser.add_argument('--output', type=pathlib.Path, default=None)
    parser.add_argument('--label', default=None, help='Title prefix for plots')
    parser.add_argument('--layers', type=int, nargs='+', default=[14, 23])
    parser.add_argument('--process', choices=['loss', 'production', 'both'], default='loss')
    parser.add_argument('--top', type=int, default=10)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None):
    args = parse_args(argv)
    label = args.label or args.input.name
    out = args.output or (DEFAULT_OUTPUT / args.input.name)
    processes = ['loss', 'production'] if args.process == 'both' else [args.process]
    for proc in processes:
        plot_case(args.input, out, args.layers, proc, args.top, label)


if __name__ == '__main__':
    main()
