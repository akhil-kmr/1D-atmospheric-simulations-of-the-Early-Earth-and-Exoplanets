#!/usr/bin/env python3
"""Plot TUV comparison metrics vs oxygen concentration (PAL).

Reads ``comparison_summary.csv`` in this directory and writes one figure per
metric at solar zenith angles 0 and 45 degrees.

Examples
--------
    python plot_tuv_comparison.py
    python plot_tuv_comparison.py --summary /path/to/comparison_summary.csv
    python plot_tuv_comparison.py --sza 0 45 60
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def thick_axes(top = False, labelleft = True, right = True, direction = 'in'):
    # Accessing the axes object and setting linewidth
    ax = plt.gca()
    ax.spines['top'].set_linewidth(2)  # Top axis
    ax.spines['bottom'].set_linewidth(2)  # Bottom axis
    ax.spines['left'].set_linewidth(2)  # Left axis
    ax.spines['right'].set_linewidth(2)  # Right axis
    plt.grid(False)
    ax.tick_params(which = 'major', axis = 'both', direction = direction, labelsize = 15, length = 6, width = 2, labelleft = labelleft, right = right, top = top)
    ax.tick_params(which = 'minor',axis = 'both', direction = direction, labelsize = 15, length = 3, width = 1, labelleft = labelleft, right = right, top = top)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')


def _set_bold_tick_labels(ax) -> None:
    """Force bold tick labels, including log-scale mathtext.

    Matplotlib log ticks use ``\\mathdefault{...}``, which ignores fontweight.
    Use a permanent ``\\mathbf`` formatter on log axes, and pin linear
    tick labels with ``set_*ticklabels(..., fontweight='bold')``.
    """
    import math
    from matplotlib.ticker import FuncFormatter

    def _bold_log_label(x, _pos=None):
        if x <= 0:
            return ""
        exp = int(round(math.log10(x)))
        if not math.isclose(x, 10 ** exp, rel_tol=1e-8, abs_tol=0.0):
            return ""
        return rf"$\mathbf{{10^{{{exp}}}}}$"

    if ax.get_xscale() == "log":
        ax.xaxis.set_major_formatter(FuncFormatter(_bold_log_label))
    else:
        xticks = ax.get_xticks()
        ax.set_xticks(xticks)
        ax.set_xticklabels(
            [lab.get_text() for lab in ax.get_xticklabels()],
            fontweight="bold",
            fontsize=15,
        )

    if ax.get_yscale() == "log":
        ax.yaxis.set_major_formatter(FuncFormatter(_bold_log_label))
    else:
        yticks = list(ax.get_yticks())
        ylim = ax.get_ylim()
        yticks = [y for y in yticks if ylim[0] <= y <= ylim[1]]
        ax.set_yticks(yticks)
        ax.set_yticklabels(
            [f"{int(y)}" if abs(y - round(y)) < 1e-8 else f"{y:g}" for y in yticks],
            fontweight="bold",
            fontsize=15,
        )


def _boldify_legend_texts(legend) -> None:
    for text in legend.get_texts():
        text.set_fontweight("bold")
        content = text.get_text()
        if "$" in content and "\\mathbf" not in content:
            # Bold any remaining math fragments (e.g. degree symbols).
            text.set_text(content.replace("$^\\circ$", r"$\mathbf{^\circ}$"))


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SUMMARY = SCRIPT_DIR / "comparison_summary.csv"
DEFAULT_OUTDIR = SCRIPT_DIR / "plots"

# PAL case key -> O2 mixing ratio relative to present atmospheric level.
CASE_O2_PAL: dict[str, float] = {
    "pal0p1": 0.001,
    "pal1": 0.01,
    "pal10": 0.1,
    "pal100": 1.0,
}

CASE_ORDER = ["pal0p1", "pal1", "pal10", "pal100"]

MODEL_ORDER = [
    "waccm",
    "waccm_eq",
    "vulcan",
    "photochem",
    "atmos",
    "kasting",
]

REFERENCE_WACCM_MODEL = "waccm"

RATIO_MODEL_ORDER = [m for m in MODEL_ORDER if m != REFERENCE_WACCM_MODEL]

MODEL_STYLE: dict[str, dict] = {
    "waccm": {
        "label": "WACCM6 (LWAV)",
        "color": "black",
        "marker": "s",
        "linestyle": "-",
        "zorder": 5,
    },
    "waccm_eq": {
        "label": "WACCM6 (equator)",
        "color": "0.35",
        "marker": "D",
        "linestyle": "--",
        "zorder": 4,
    },
    "vulcan": {
        "label": "VULCAN",
        "color": "magenta",
        "marker": "o",
        "linestyle": "-",
        "zorder": 3,
    },
    "photochem": {
        "label": "Photochem",
        "color": "blue",
        "marker": "^",
        "linestyle": "-",
        "zorder": 3,
    },
    "atmos": {
        "label": "Atmos",
        "color": "darkorange",
        "marker": "v",
        "linestyle": "-",
        "zorder": 3,
    },
    "kasting": {
        "label": "Kasting 1D",
        "color": "teal",
        "marker": "P",
        "linestyle": "-",
        "zorder": 3,
    },
}

PLOT_SPECS: list[tuple[str, str, str]] = [
    ("DNA_damage", "DNA damage rate", "DNA_damage_vs_O2"),
    ("UV_index", "UV index", "UV_index_vs_O2"),
    ("P_Dam_C_1971", "Plant damage (Caldwell 1971)", "plant_damage_C1971_vs_O2"),
    ("P_Dam_FC_2003", "Plant damage (Flint & Caldwell 2003)", "plant_damage_FC2003_vs_O2"),
    (
        "P_Dam_FC_2003_ext390",
        "Plant damage (FC 2003, ext. 390 nm)",
        "plant_damage_FC2003_ext390_vs_O2",
    ),
    (
        "UVB_280_315",
        r"UVB irradiance (280--315 nm) [mW m$^{-2}$ nm$^{-1}$]",
        "UVB_280_315_vs_O2",
    ),
    (
        "UVBstar_280_320",
        r"UVB* irradiance (280--320 nm) [mW m$^{-2}$ nm$^{-1}$]",
        "UVBstar_280_320_vs_O2",
    ),
    (
        "UVA_315_400",
        r"UVA irradiance (315--400 nm) [mW m$^{-2}$ nm$^{-1}$]",
        "UVA_315_400_vs_O2",
    ),
]


def load_summary(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Comparison summary not found: {path}")
    df = pd.read_csv(path)
    required = {"case", "model", "sza_deg", *{spec[0] for spec in PLOT_SPECS}}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"comparison_summary.csv missing columns: {sorted(missing)}")
    df = df.copy()
    df["o2_pal"] = df["case"].map(CASE_O2_PAL)
    if df["o2_pal"].isna().any():
        unknown = sorted(df.loc[df["o2_pal"].isna(), "case"].unique())
        raise ValueError(f"Unknown comparison case(s) for O2 mapping: {unknown}")
    df["sza_deg"] = df["sza_deg"].astype(float)
    return df


def subset_at_sza(df: pd.DataFrame, sza_deg: float) -> pd.DataFrame:
    tol = 0.25
    sub = df[np.isclose(df["sza_deg"], sza_deg, atol=tol)].copy()
    if sub.empty:
        available = sorted(df["sza_deg"].unique())
        raise ValueError(
            f"No rows near SZA = {sza_deg:g} deg (available: {available})"
        )
    return sub


def model_points(sub: pd.DataFrame, model: str) -> pd.DataFrame:
    rows = []
    for case in CASE_ORDER:
        match = sub[(sub["case"] == case) & (sub["model"] == model)]
        if match.empty:
            continue
        rows.append(match.iloc[0])
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    return out.sort_values("o2_pal")


def model_dna_ratio_points(sub: pd.DataFrame, model: str) -> pd.DataFrame:
    """Return O2 PAL and WACCM LWAV / model DNA_damage for one model."""
    ref = model_points(sub, REFERENCE_WACCM_MODEL)
    pts = model_points(sub, model)
    if ref.empty or pts.empty:
        return pd.DataFrame()

    ref_by_case = ref.set_index("case")["DNA_damage"]
    rows = []
    for _, row in pts.iterrows():
        ref_val = ref_by_case.get(row["case"])
        if ref_val is None or ref_val <= 0 or row["DNA_damage"] <= 0:
            continue
        rows.append({
            "case": row["case"],
            "o2_pal": row["o2_pal"],
            "dna_ratio": ref_val / row["DNA_damage"],
        })
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values("o2_pal")


DNA_RATIO_YLABEL = (
    "WACCM6 global mean / model DNA damage"#"\n"
    #"(>1: WACCM higher; <1: model higher)"
)
DNA_RATIO_TITLE = "WACCM6 global-mean DNA damage vs other models"


def _annotate_dna_ratio_plot(ax, fontsize: int = 13) -> None:
    """Arrow and note: WACCM excess damage grows as O2 falls (vs 1D models)."""
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    log_x0, log_x1 = np.log10(xlim[0]), np.log10(xlim[1])
    # Arrow from higher O2 (right) toward lower O2 (left), ratio increasing upward.
    x_arrow_tail = 10 ** (0.35 * log_x0 + 0.65 * log_x1)
    x_arrow_head = 10 ** (0.62 * log_x0 + 0.22 * log_x1)
    y_arrow_tail = ylim[0] + 0.30 * (ylim[1] - ylim[0])
    y_arrow_head = ylim[0] + 0.78 * (ylim[1] - ylim[0])

    ax.annotate(
        "",
        xy=(x_arrow_head, y_arrow_head),
        xytext=(x_arrow_tail, y_arrow_tail),
        arrowprops=dict(arrowstyle="-|>", lw=2.2, color="0.25"),
        zorder=6,
    )
    ax.text(
        x_arrow_head,
        y_arrow_head + 0.06 * (ylim[1] - ylim[0]),
        "WACCM6 excess DNA damage\nincreases at lower O₂",
        fontsize=fontsize - 1,
        ha="center",
        va="bottom",
        color="0.15",
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.75", alpha=0.93),
        zorder=6,
    )
    '''
    ax.text(
        0.98,
        1.02,
        "Ratio = WACCM6 global mean / model\n(reference line: equal to WACCM6)",
        transform=ax.transAxes,
        fontsize=fontsize - 2,
        ha="right",
        va="bottom",
        color="0.2",
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.75", alpha=0.93),
        zorder=6,
    )
    '''

def plot_metric_vs_o2(
    df: pd.DataFrame,
    metric: str,
    ylabel: str,
    sza_deg: float,
    out_path: Path,
) -> None:
    sub = subset_at_sza(df, sza_deg)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for model in MODEL_ORDER:
        pts = model_points(sub, model)
        if pts.empty:
            continue
        style = MODEL_STYLE[model]
        ax.plot(
            pts["o2_pal"],
            pts[metric],
            color=style["color"],
            linestyle=style["linestyle"],
            marker=style["marker"],
            markersize=8,
            linewidth=2,
            label=style["label"],
            zorder=style["zorder"],
        )

    ax.set_xscale("log")
    ax.set_xlabel("Oxygen concentration [PAL]", fontsize=13, weight="bold")
    ax.set_ylabel(ylabel, fontsize=13, weight="bold")
    ax.set_title(
        rf"{ylabel} vs O$_2$  (SZA = {sza_deg:g}$^\circ$)",
        fontsize=14,
        weight="bold",
    )
    ax.grid(True, which="both", linestyle=":", alpha=0.4)
    ax.legend(loc="best", fontsize=9, framealpha=0.92)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    #plt.close(fig)
    print(f"Wrote {out_path}")


def plot_dna_damage_ratio_vs_o2(
    df: pd.DataFrame,
    sza_deg: float,
    out_path: Path,
) -> None:
    """Plot WACCM6 global-mean / model DNA damage vs O2 PAL."""
    sub = subset_at_sza(df, sza_deg)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    for model in RATIO_MODEL_ORDER:
        pts = model_dna_ratio_points(sub, model)
        if pts.empty:
            continue
        style = MODEL_STYLE[model]
        ax.plot(
            pts["o2_pal"],
            pts["dna_ratio"],
            color=style["color"],
            linestyle=style["linestyle"],
            marker=style["marker"],
            markersize=8,
            linewidth=2,
            label=style["label"],
            zorder=style["zorder"],
        )

    ax.axhline(
        1.0,
        color=MODEL_STYLE[REFERENCE_WACCM_MODEL]["color"],
        linestyle=":",
        linewidth=1.5,
        label="Equal to WACCM6\nglobal mean",
        zorder=1,
    )
    ax.set_xscale("log")
    ax.set_xlabel("Oxygen concentration [PAL]", fontsize=13, weight="bold")
    ax.set_ylabel(DNA_RATIO_YLABEL, fontsize=13, weight="bold")
    ax.set_title(
        rf"{DNA_RATIO_TITLE}  (SZA = {sza_deg:g}$^\circ$)",
        fontsize=14,
        weight="bold",
    )
    ax.grid(True, which="both", linestyle=":", alpha=0.4)
    ax.legend(loc="best", fontsize=9, framealpha=0.92)
    _annotate_dna_ratio_plot(ax, fontsize=13)
    fig.tight_layout()
    thick_axes(top=True, right=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    #plt.close(fig)
    print(f"Wrote {out_path}")


def plot_dna_damage_ratio_vs_o2_dual_sza(
    df: pd.DataFrame,
    sza_thick: float,
    sza_dashed: float,
    out_path: Path,
) -> None:
    """Plot DNA damage ratios at two SZAs on one axes (thick = low SZA, dashed = high SZA)."""
    sub_thick = subset_at_sza(df, sza_thick)
    sub_dashed = subset_at_sza(df, sza_dashed)

    # Bold regular text + mathtext (log ticks use mathtext; fontweight alone is ignored).
    with plt.rc_context({
        "font.weight": "bold",
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "mathtext.default": "bf",
        "mathtext.fontset": "dejavusans",
    }):
        fig, ax = plt.subplots(figsize=(11, 6))
        for model in RATIO_MODEL_ORDER:
            style = MODEL_STYLE[model]
            pts0 = model_dna_ratio_points(sub_thick, model)
            if not pts0.empty:
                ax.plot(
                    pts0["o2_pal"],
                    pts0["dna_ratio"],
                    color=style["color"],
                    linestyle=style["linestyle"],
                    marker=style["marker"],
                    markersize=7,
                    linewidth=3,
                    label=style["label"],
                    zorder=style["zorder"],
                )
            pts45 = model_dna_ratio_points(sub_dashed, model)
            if not pts45.empty:
                ax.plot(
                    pts45["o2_pal"],
                    pts45["dna_ratio"],
                    color=style["color"],
                    linestyle="--",
                    marker=style["marker"],
                    markersize=7,
                    linewidth=2,
                    zorder=style["zorder"] - 0.1,
                )

        ax.axhline(
            1.0,
            color=MODEL_STYLE[REFERENCE_WACCM_MODEL]["color"],
            linestyle=":",
            linewidth=1.5,
            label="Equal to WACCM6\nglobal mean",
            zorder=1,
        )
        ax.set_xscale("log")
        ax.set_xlabel("Oxygen concentration [PAL]", fontsize=15, weight="bold")
        ax.set_ylabel(DNA_RATIO_YLABEL, fontsize=15, weight="bold")
        ax.set_title(
            rf"{DNA_RATIO_TITLE}  (SZA = {sza_thick:g}$^\circ$ thick, "
            rf"{sza_dashed:g}$^\circ$ dashed)",
            fontsize=15,
            weight="bold",
        )
        thick_axes(top=True, right=True)

        from matplotlib.lines import Line2D

        legend_prop = {"weight": "bold", "size": 15}
        model_leg = ax.legend(loc=0, prop=legend_prop, frameon=False)
        ax.add_artist(model_leg)
        _boldify_legend_texts(model_leg)
        sza_handles = [
            Line2D([0], [0], color="0.2", lw=2, ls="-", label=rf"SZA = {sza_thick:g}$^\circ$"),
            Line2D([0], [0], color="0.2", lw=2, ls="--", label=rf"SZA = {sza_dashed:g}$^\circ$"),
        ]
        sza_leg = ax.legend(
            handles=sza_handles, loc=(0.75, 0.28), prop=legend_prop, frameon=False
        )
        _boldify_legend_texts(sza_leg)
        _annotate_dna_ratio_plot(ax, fontsize=15)

        fig.tight_layout()
        fig.canvas.draw()
        _set_bold_tick_labels(ax)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
    print(f"Wrote {out_path}")
def make_all_plots(
    summary_path: Path,
    outdir: Path,
    sza_values: list[float],
    *,
    dna_only: bool = False,
) -> None:
    df = load_summary(summary_path)
    for sza_deg in sza_values:
        if abs(sza_deg - round(sza_deg)) < 0.01:
            sza_tag = f"sza{int(round(sza_deg))}"
        else:
            sza_tag = f"sza{sza_deg:g}".replace(".", "p")
        if dna_only:
            plot_metric_vs_o2(
                df,
                "DNA_damage",
                "DNA damage rate",
                sza_deg,
                outdir / f"DNA_damage_vs_O2_{sza_tag}.png",
            )
        else:
            for metric, ylabel, stem in PLOT_SPECS:
                out_path = outdir / f"{stem}_{sza_tag}.png"
                plot_metric_vs_o2(df, metric, ylabel, sza_deg, out_path)
        ratio_path = outdir / f"DNA_damage_ratio_to_WACCM_{sza_tag}.png"
        plot_dna_damage_ratio_vs_o2(df, sza_deg, ratio_path)

    dual_ratio_path = outdir / "DNA_damage_ratio_to_WACCM_sza0_sza45.png"
    plot_dna_damage_ratio_vs_o2_dual_sza(df, 0.0, 45.0, dual_ratio_path)


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Plot TUV comparison dose metrics vs O2 concentration."
    )
    p.add_argument(
        "--summary",
        type=Path,
        default=DEFAULT_SUMMARY,
        help=f"Path to comparison_summary.csv (default: {DEFAULT_SUMMARY})",
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Directory for PNG output (default: plots/ next to summary, or plots_<tag>/).",
    )
    p.add_argument(
        "--sza",
        type=float,
        nargs="+",
        default=[0.0, 45.0],
        help="Solar zenith angle(s) to plot (default: 0 45).",
    )
    p.add_argument(
        "--dna-only",
        action="store_true",
        help="Only write DNA damage absolute and ratio figures.",
    )
    p.add_argument(
        "--wtag",
        type=str,
        default=None,
        metavar="TAG",
        help="Wavelength tag such as w100-420; loads comparison_summary_<tag>.csv if --summary omitted.",
    )
    return p


def main() -> None:
    args = _build_arg_parser().parse_args()
    summary = args.summary.resolve()
    if args.wtag:
        tag = args.wtag if args.wtag.startswith("w") else f"w{args.wtag}"
        tagged = SCRIPT_DIR / f"comparison_summary_{tag}.csv"
        if args.summary == DEFAULT_SUMMARY:
            summary = tagged.resolve()
        if args.outdir is None:
            args.outdir = SCRIPT_DIR / f"plots_{tag}"
    outdir = (args.outdir or DEFAULT_OUTDIR).resolve()
    make_all_plots(summary, outdir, args.sza, dna_only=args.dna_only)


if __name__ == "__main__":
    main()
