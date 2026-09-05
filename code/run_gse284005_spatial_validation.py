#!/usr/bin/env python3
"""Independent spatial audit of endothelial inflammatory programs in GSE284005.

The analysis is deliberately donor based.  Each MERFISH section represents one
donor in the deposited series (14 MS and three non-neurological controls).  The
primary paired contrast is restricted to MS donors with at least 50 annotated
endothelial cells in both the vascular-immune (Vas_Imm) and demyelinated white
matter (DMWM) regions.  Only prespecified modules with at least four measured
genes and at least 40% panel coverage are scored.  CLDN5 is shown separately;
the targeted 500-gene panel does not support a composite BBB-specialization
score.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "public_data" / "GSE284005"
RESULTS = ROOT / "results"
SOURCE = ROOT / "supplementary_tables"
FIGURES = ROOT / "figures" / "supplementary"

ENDOTHELIAL_SUBTYPES = {
    "Homeo.Endo",
    "Stress.Endo",
    "Repair.Endo",
    "inflammatory.Endo",
}

MODULES = {
    "Endothelial immune activation": [
        "STAT1", "STAT3", "IRF1", "CXCL9", "CXCL10", "CCL2", "B2M",
        "TAP1", "HLA-DRA", "HLA-DRB1", "H2-AA", "H2-AB1", "ICAM1",
        "VCAM1", "SELE", "SELP", "IFIT1", "IFIT2", "IFIT3", "ISG15",
    ],
    "Adhesion and trafficking": [
        "ICAM1", "VCAM1", "SELE", "SELP", "CCL2", "CXCL10",
    ],
    "IFN and antigen presentation": [
        "STAT1", "IRF1", "CXCL9", "CXCL10", "B2M", "TAP1", "HLA-DRA",
        "HLA-DRB1", "IFIT1", "IFIT2", "IFIT3", "ISG15",
    ],
}

REGIONS = ("DMWM", "Vas_Imm")
MIN_EC_PER_REGION = 50


def exact_sign_flip_p(values: np.ndarray) -> float:
    """Two-sided exact sign-flip P value for a paired mean difference."""
    values = np.asarray(values, dtype=float)
    obs = abs(values.mean())
    perm = []
    for signs in itertools.product((-1.0, 1.0), repeat=len(values)):
        perm.append(abs(np.mean(values * np.asarray(signs))))
    return float(np.mean(np.asarray(perm) >= obs - 1e-12))


def zscore_rows(frame: pd.DataFrame) -> pd.DataFrame:
    mean = frame.mean(axis=1)
    sd = frame.std(axis=1, ddof=1).replace(0, 1.0)
    return frame.sub(mean, axis=0).div(sd, axis=0)


def load_sample(meta_path: Path) -> tuple[str, pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(meta_path, sep="\t")
    sample = str(meta.loc[0, "cells"]).split("_")[0]
    count_path = Path(str(meta_path).replace("_celltypes.tsv.gz", "_count.tsv.gz"))
    counts = pd.read_csv(count_path, sep="\t", index_col=0)
    counts.index = counts.index.astype(str).str.upper()
    return sample, meta, counts


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    SOURCE.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    aggregated: dict[tuple[str, str], pd.Series] = {}
    cell_rows: list[dict[str, object]] = []

    for meta_path in sorted(INPUT.glob("GSM*_celltypes.tsv.gz")):
        sample, meta, counts = load_sample(meta_path)
        ec = meta[meta["clean_sub"].isin(ENDOTHELIAL_SUBTYPES)].copy()
        for region in REGIONS:
            cells = ec.loc[ec["Region_banksy"].eq(region), "cells"].tolist()
            cells = [c for c in cells if c in counts.columns]
            sub = ec[ec["cells"].isin(cells) & ec["Region_banksy"].eq(region)]
            cell_rows.append(
                {
                    "donor": sample,
                    "region": region,
                    "endothelial_cells": len(cells),
                    "homeostatic_endothelial_cells": int(sub["clean_sub"].eq("Homeo.Endo").sum()),
                    "stress_endothelial_cells": int(sub["clean_sub"].eq("Stress.Endo").sum()),
                    "repair_endothelial_cells": int(sub["clean_sub"].eq("Repair.Endo").sum()),
                    "inflammatory_endothelial_cells": int(sub["clean_sub"].eq("inflammatory.Endo").sum()),
                }
            )
            if cells:
                aggregated[(sample, region)] = counts.loc[:, cells].sum(axis=1)

    cell_counts = pd.DataFrame(cell_rows)
    pivot = cell_counts.pivot(index="donor", columns="region", values="endothelial_cells").fillna(0)
    eligible = sorted(
        pivot.index[
            (pivot["DMWM"] >= MIN_EC_PER_REGION)
            & (pivot["Vas_Imm"] >= MIN_EC_PER_REGION)
        ].tolist()
    )

    profiles = []
    profile_names = []
    for donor in eligible:
        for region in REGIONS:
            profiles.append(aggregated[(donor, region)])
            profile_names.append(f"{donor}|{region}")
    raw = pd.concat(profiles, axis=1)
    raw.columns = profile_names
    cpm = raw.div(raw.sum(axis=0).clip(lower=1), axis=1) * 1e6
    log_cpm = np.log1p(cpm)
    z = zscore_rows(log_cpm)

    coverage_rows = []
    score_rows = []
    measured = set(raw.index)
    retained_modules = {}
    for module, expected in MODULES.items():
        present = [g for g in expected if g in measured]
        retained = len(present) >= 4 and len(present) / len(expected) >= 0.40
        coverage_rows.append(
            {
                "module": module,
                "genes_expected": len(expected),
                "genes_measured": len(present),
                "coverage_fraction": len(present) / len(expected),
                "measured_genes": "; ".join(present),
                "retained_for_scoring": retained,
            }
        )
        if retained:
            retained_modules[module] = present
            scores = z.loc[present].mean(axis=0)
            for profile, value in scores.items():
                donor, region = profile.split("|")
                score_rows.append(
                    {"donor": donor, "region": region, "module": module, "score": value}
                )

    scores = pd.DataFrame(score_rows)
    summary_rows = []
    for module in retained_modules:
        wide = scores[scores["module"].eq(module)].pivot(index="donor", columns="region", values="score")
        diff = (wide["Vas_Imm"] - wide["DMWM"]).dropna()
        summary_rows.append(
            {
                "feature": module,
                "n_donors": len(diff),
                "mean_paired_difference": diff.mean(),
                "median_paired_difference": diff.median(),
                "donors_higher_in_Vas_Imm": int((diff > 0).sum()),
                "exact_two_sided_sign_flip_p": exact_sign_flip_p(diff.to_numpy()),
                "interpretation": "prespecified module score",
            }
        )

    gene_rows = []
    focus_genes = sorted(set(sum(retained_modules.values(), [])) | {"CLDN5"})
    for gene in focus_genes:
        if gene not in log_cpm.index:
            continue
        diffs = []
        for donor in eligible:
            d = float(log_cpm.loc[gene, f"{donor}|Vas_Imm"] - log_cpm.loc[gene, f"{donor}|DMWM"])
            diffs.append(d)
            gene_rows.append(
                {
                    "donor": donor,
                    "gene": gene,
                    "log_CPM_difference_Vas_Imm_minus_DMWM": d,
                }
            )
        summary_rows.append(
            {
                "feature": gene,
                "n_donors": len(diffs),
                "mean_paired_difference": float(np.mean(diffs)),
                "median_paired_difference": float(np.median(diffs)),
                "donors_higher_in_Vas_Imm": int(np.sum(np.asarray(diffs) > 0)),
                "exact_two_sided_sign_flip_p": exact_sign_flip_p(np.asarray(diffs)),
                "interpretation": "single-gene log(CPM+1) difference",
            }
        )

    # State proportions are descriptive because the subtype labels were learned
    # in the source study and are not independent molecular measurements.
    eligible_cells = cell_counts[cell_counts["donor"].isin(eligible)].copy()
    eligible_cells["stress_or_inflammatory_fraction"] = (
        eligible_cells["stress_endothelial_cells"]
        + eligible_cells["inflammatory_endothelial_cells"]
    ) / eligible_cells["endothelial_cells"]

    fraction_wide = eligible_cells.pivot(
        index="donor", columns="region", values="stress_or_inflammatory_fraction"
    )
    fraction_diff = (fraction_wide["Vas_Imm"] - fraction_wide["DMWM"]).dropna()
    summary_rows.append(
        {
            "feature": "Stress/inflammatory endothelial fraction",
            "n_donors": len(fraction_diff),
            "mean_paired_difference": fraction_diff.mean(),
            "median_paired_difference": fraction_diff.median(),
            "donors_higher_in_Vas_Imm": int((fraction_diff > 0).sum()),
            "exact_two_sided_sign_flip_p": exact_sign_flip_p(fraction_diff.to_numpy()),
            "interpretation": "source-study endothelial-state annotation; descriptive",
        }
    )

    coverage = pd.DataFrame(coverage_rows)
    summary = pd.DataFrame(summary_rows)
    genes = pd.DataFrame(gene_rows)

    coverage.to_csv(RESULTS / "GSE284005_module_coverage.csv", index=False)
    scores.to_csv(RESULTS / "GSE284005_donor_region_module_scores.csv", index=False)
    genes.to_csv(RESULTS / "GSE284005_donor_gene_differences.csv", index=False)
    summary.to_csv(RESULTS / "GSE284005_paired_summary.csv", index=False)
    eligible_cells.to_csv(RESULTS / "GSE284005_endothelial_cell_counts.csv", index=False)

    source = summary.merge(
        coverage[["module", "genes_expected", "genes_measured", "coverage_fraction", "measured_genes"]],
        how="left", left_on="feature", right_on="module",
    ).drop(columns="module")
    source.insert(0, "dataset", "GSE284005")
    source.insert(1, "contrast", "endothelial Vas_Imm versus DMWM within donor")
    source.to_csv(SOURCE / "S16_GSE284005_spatial_endothelial_validation.csv", index=False)

    audit = pd.DataFrame(
        [
            {
                "dataset": "GSE168202",
                "proposed_role": "fourth EAE endothelial meta-analysis cohort",
                "eligibility_decision": "excluded",
                "reason": "brain parenchyma has one library at each of naive, onset and peak stages; no biological replication for a brain endothelial disease effect",
            },
            {
                "dataset": "GSE284005",
                "proposed_role": "independent spatial MS validation",
                "eligibility_decision": "included with restricted scope",
                "reason": "paired donor-level Vas_Imm versus DMWM analysis was feasible in five MS donors with at least 50 endothelial cells per region; the targeted panel supports inflammatory modules but not a composite BBB module",
            },
        ]
    )
    audit.to_csv(SOURCE / "S17_public_dataset_eligibility_audit.csv", index=False)

    # Supplementary figure: panel coverage and expanded gene-level audit.  The
    # paired module and cell-state plots now appear in main Figure 3.
    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans"})
    fig = plt.figure(figsize=(8.2, 7.0))
    fig.subplots_adjust(left=.11, right=.98, top=.89, bottom=.09, hspace=.52, wspace=.34)
    gs = fig.add_gridspec(2, 2, height_ratios=[.88, 1.12], width_ratios=[.82, 1.18])
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])]

    ax = axes[0]
    coverage_plot = coverage.copy()
    coverage_plot["short"] = coverage_plot["module"].replace({
        "Endothelial immune activation": "Immune activation",
        "Adhesion and trafficking": "Adhesion/trafficking",
        "IFN and antigen presentation": "IFN/antigen presentation",
    })
    y = np.arange(len(coverage_plot))
    ax.barh(y, coverage_plot.coverage_fraction * 100, color=["#B33B3B", "#C56A4F", "#8F4A79"], height=.62)
    ax.set_yticks(y, coverage_plot.short); ax.invert_yaxis(); ax.set_xlim(0, 100)
    ax.set_xlabel("Prespecified genes measured (%)")
    ax.set_title("A  Targeted-panel module coverage", loc="left", fontweight="bold")
    for i,row in enumerate(coverage_plot.itertuples()):
        ax.text(row.coverage_fraction*100+2,i,f"{int(row.genes_measured)}/{int(row.genes_expected)}",va='center',fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    gene_direction = summary[summary.interpretation.str.startswith("single-gene")].copy()
    gene_direction = gene_direction.sort_values(["donors_higher_in_Vas_Imm", "mean_paired_difference"])
    y = np.arange(len(gene_direction))
    direction_colors = np.where(gene_direction.donors_higher_in_Vas_Imm >= 3, "#B33B3B", "#4676A9")
    ax.barh(y, gene_direction.donors_higher_in_Vas_Imm, color=direction_colors, alpha=.86)
    ax.axvline(2.5, color="#555555", lw=.8, ls="--")
    ax.set_yticks(y, gene_direction.feature); ax.set_xlim(0,5.35); ax.set_xticks(range(6))
    ax.set_xlabel("Donors higher in vascular–immune region (of 5)", fontsize=8.5)
    ax.set_title("B  Gene-level directional consistency", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    gene_summary = (
        genes.groupby("gene")["log_CPM_difference_Vas_Imm_minus_DMWM"]
        .agg(["mean", "median"])
        .sort_values("mean")
    )
    show = gene_summary.iloc[np.r_[0:min(6, len(gene_summary)), max(6, len(gene_summary)-7):len(gene_summary)]].drop_duplicates()
    ax = axes[2]
    y = np.arange(len(show))
    ax.barh(y, show["mean"], color=np.where(show["mean"] >= 0, "#B33B3B", "#4676A9"), alpha=0.85)
    ax.axvline(0, color="#555555", lw=0.8)
    ax.set_yticks(y, show.index)
    ax.set_xlabel("Mean paired log(CPM+1) difference")
    ax.set_title("C  Covered genes (selected extremes)", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "GSE284005 donor-level spatial audit: vascular–immune regions versus demyelinated white matter",
        fontsize=11, y=.985,
    )
    for ext, dpi in [("png", 300), ("tiff", 600)]:
        fig.savefig(FIGURES / f"Figure_S9_GSE284005_spatial_validation.{ext}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    print("Eligible MS donors:", ", ".join(eligible))
    print("\nModule coverage:")
    print(coverage.to_string(index=False))
    print("\nPaired summary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
