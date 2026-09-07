#!/usr/bin/env python3
"""Regenerate Supplementary Figure S9 from frozen GSE284005 result tables."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures" / "supplementary"


def main() -> None:
    coverage = pd.read_csv(RESULTS / "GSE284005_module_coverage.csv")
    summary = pd.read_csv(RESULTS / "GSE284005_paired_summary.csv")
    genes = pd.read_csv(RESULTS / "GSE284005_donor_gene_differences.csv")
    FIGURES.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans"})
    fig = plt.figure(figsize=(8.2, 7.0))
    fig.subplots_adjust(left=.11, right=.98, top=.89, bottom=.09, hspace=.52, wspace=.34)
    grid = fig.add_gridspec(2, 2, height_ratios=[.88, 1.12], width_ratios=[.82, 1.18])
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, :])]

    ax = axes[0]
    shown = coverage.copy()
    shown["short"] = shown["module"].replace({
        "Endothelial immune activation": "Immune activation",
        "Adhesion and trafficking": "Adhesion/trafficking",
        "IFN and antigen presentation": "IFN/antigen presentation",
    })
    y = np.arange(len(shown))
    ax.barh(y, shown.coverage_fraction * 100, color=["#B33B3B", "#C56A4F", "#8F4A79"], height=.62)
    ax.set_yticks(y, shown.short)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Prespecified genes measured (%)")
    ax.set_title("A  Targeted-panel module coverage", loc="left", fontweight="bold")
    for i, row in enumerate(shown.itertuples()):
        ax.text(row.coverage_fraction * 100 + 2, i, f"{int(row.genes_measured)}/{int(row.genes_expected)}", va="center", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    directions = summary[summary.interpretation.str.startswith("single-gene")].copy()
    directions = directions.sort_values(["donors_higher_in_Vas_Imm", "mean_paired_difference"])
    y = np.arange(len(directions))
    colors = np.where(directions.donors_higher_in_Vas_Imm >= 3, "#B33B3B", "#4676A9")
    ax.barh(y, directions.donors_higher_in_Vas_Imm, color=colors, alpha=.86)
    ax.axvline(2.5, color="#555555", lw=.8, ls="--")
    ax.set_yticks(y, directions.feature)
    ax.set_xlim(0, 5.35)
    ax.set_xticks(range(6))
    ax.set_xlabel("Donors higher in vascular–immune region (of 5)", fontsize=8.5)
    ax.set_title("B  Gene-level directional consistency", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    gene_summary = (
        genes.groupby("gene")["log_CPM_difference_Vas_Imm_minus_DMWM"]
        .agg(["mean", "median"])
        .sort_values("mean")
    )
    shown_genes = gene_summary.iloc[np.r_[0:min(6, len(gene_summary)), max(6, len(gene_summary)-7):len(gene_summary)]].drop_duplicates()
    ax = axes[2]
    y = np.arange(len(shown_genes))
    ax.barh(y, shown_genes["mean"], color=np.where(shown_genes["mean"] >= 0, "#B33B3B", "#4676A9"), alpha=.85)
    ax.axvline(0, color="#555555", lw=.8)
    ax.set_yticks(y, shown_genes.index)
    ax.set_xlabel("Mean paired log(CPM+1) difference")
    ax.set_title("C  Covered genes (selected extremes)", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "GSE284005 donor-level spatial audit: vascular–immune regions versus demyelinated white matter",
        fontsize=11, y=.985,
    )
    fig.savefig(FIGURES / "Figure_S9_GSE284005_spatial_validation.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
