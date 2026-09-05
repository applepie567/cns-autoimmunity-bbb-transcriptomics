"""Create the detailed supplementary figure for published functional/protein data."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = ROOT / "figures" / "supplementary"
OUT.mkdir(parents=True, exist_ok=True)

COLORS = {"IgG": "#7D8790", "r84": "#2A9D8F", "barrier": "#2F6F9F", "lesion": "#C84C35", "nawm": "#D98BA5", "control": "#3A9D6F"}


def panel(ax: plt.Axes, letter: str, title: str) -> None:
    ax.set_title(f"{letter}  {title}", loc="left", fontsize=9, fontweight="bold", pad=8)


def main() -> None:
    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.2, "axes.spines.top": False, "axes.spines.right": False})
    mouse = pd.read_csv(RESULTS / "published_source_mouse_level_barrier_and_vascular_data.csv")
    summary = pd.read_csv(RESULTS / "published_source_barrier_and_vascular_summary.csv")
    teer = pd.read_csv(RESULTS / "published_mBEC_TEER_well_level_source_data.csv")
    ihc = pd.read_csv(RESULTS / "published_human_MS_IHC_image_level_source_data.csv")

    fig, axs = plt.subplots(2, 2, figsize=(7.6, 6.7))
    fig.subplots_adjust(left=.10, right=.98, top=.91, bottom=.11, hspace=.55, wspace=.40)
    rng = np.random.default_rng(20260901)
    fig.legend(
        handles=[
            Line2D([0], [0], marker="o", ls="", color=COLORS["IgG"], label="IgG"),
            Line2D([0], [0], marker="o", ls="", color=COLORS["r84"], label="r84"),
        ],
        frameon=False, loc="upper center", bbox_to_anchor=(.5, .985), ncol=2, fontsize=6.6,
    )

    ax = axs[0, 0]
    panel(ax, "A", "Chronic EAE: vascular target engagement")
    endpoints = ["Venous endothelial proliferation", "Venous coverage", "Endomucin coverage", "Venous width", "Ki67-positive endothelial cells"]
    for i, endpoint in enumerate(endpoints):
        block = mouse[mouse.endpoint == endpoint]
        pooled = block.mouse_mean.to_numpy()
        pooled_mean = np.mean(pooled)
        pooled_sd = np.std(pooled, ddof=1) or 1
        for j, group in enumerate(["IgG", "r84"]):
            vals = block[block.group.str.lower() == group.lower()].mouse_mean.to_numpy()
            x = i + (-.16 if group == "IgG" else .16) + rng.uniform(-.035, .035, len(vals))
            z = (vals - pooled_mean) / pooled_sd
            ax.scatter(x, z, s=16, color=COLORS[group], alpha=.85, edgecolor="white", lw=.3)
            ax.plot([i + (-.23 if group == "IgG" else .09), i + (-.09 if group == "IgG" else .23)], [z.mean(), z.mean()], color="#26343D", lw=1)
    ax.axhline(0, color="#CDD4D8", lw=.7)
    ax.set_xticks(range(len(endpoints)), ["Proliferation", "Venous\ncoverage", "Endomucin", "Venous\nwidth", "Ki67+ EC"], rotation=25, ha="right")
    ax.set_ylabel("Within-endpoint standardized mouse mean")
    ax.text(.99, .02, "n=9 per group; fields averaged within mouse", transform=ax.transAxes, ha="right", fontsize=6, color="#66737C")

    ax = axs[0, 1]
    panel(ax, "B", "Chronic EAE: direct leakage measurements")
    leakage_max = mouse[mouse.endpoint.isin(["IgG leakage", "Fibrinogen leakage"])].mouse_mean.max()
    ax.set_ylim(-.5, leakage_max * 1.18)
    for i, endpoint in enumerate(["IgG leakage", "Fibrinogen leakage"]):
        block = mouse[mouse.endpoint == endpoint]
        for j, group in enumerate(["IgG", "r84"]):
            vals = block[block.group.str.lower() == group.lower()].mouse_mean.to_numpy()
            x = i + (-.14 if group == "IgG" else .14) + rng.uniform(-.045, .045, len(vals))
            ax.scatter(x, vals, s=23, color=COLORS[group], edgecolor="white", lw=.4)
            ax.plot([i + (-.22 if group == "IgG" else .06), i + (-.06 if group == "IgG" else .22)], [vals.mean(), vals.mean()], color="#26343D", lw=1.2)
        row = summary[summary.endpoint == endpoint].iloc[0]
        ax.text(i, .95, f"exact P={row.p_exact_two_sided:.3f}", transform=ax.get_xaxis_transform(), ha="center", fontsize=6.3)
    ax.set_xticks([0, 1], ["IgG leakage", "Fibrinogen leakage"])
    ax.set_ylabel("Mouse-level mean")

    ax = axs[1, 0]
    panel(ax, "C", "Primary brain endothelial TEER")
    phases = ["acute drop", "chronic drop"]
    conditions = ["Ctrl", "VEGFa + IgG", "VEGFa + r84"]
    cols = ["#7D8790", "#C84C35", "#2A9D8F"]
    for i, phase in enumerate(phases):
        for j, cond in enumerate(conditions):
            vals = teer[(teer.phase == phase) & (teer.condition == cond)].TEER_change.to_numpy()
            x0 = i + (j - 1) * .23
            ax.scatter(x0 + rng.uniform(-.035, .035, len(vals)), vals, s=18, color=cols[j], edgecolor="white", lw=.3)
            ax.plot([x0-.07, x0+.07], [vals.mean(), vals.mean()], color="#26343D", lw=1.1)
    ax.set_xticks([0, 1], ["Acute drop", "Chronic drop"]); ax.set_ylabel("Change in resistance")
    for cond, col in zip(conditions, cols): ax.scatter([], [], color=col, label=cond)
    ax.legend(frameon=False, fontsize=5.8, ncol=3, loc="lower center", bbox_to_anchor=(.5, -.22))
    ax.text(.99, .97, "7–8 wells/condition; descriptive only", transform=ax.transAxes, ha="right", va="top", fontsize=6, color="#66737C")

    ax = axs[1, 1]
    panel(ax, "D", "Human MS vascular immunohistochemistry")
    markers = ["CD31", "EGFL7", "MCAM"]; regions = ["Lesion MS", "NAWM MS", "NAWM HC"]
    rcols = [COLORS["lesion"], COLORS["nawm"], COLORS["control"]]
    for i, marker in enumerate(markers):
        for j, region in enumerate(regions):
            vals = ihc[(ihc.marker == marker) & (ihc.region == region)].image_level_value.to_numpy()
            x0 = i + (j - 1) * .23
            ax.scatter(x0 + rng.uniform(-.035, .035, len(vals)), vals, s=17, color=rcols[j], alpha=.8, edgecolor="white", lw=.3)
            ax.plot([x0-.07, x0+.07], [vals.mean(), vals.mean()], color="#26343D", lw=1.1)
    ax.set_xticks(range(3), markers); ax.set_ylabel("Positive vessel fragments")
    ax.set_ylim(0, max(24, ihc.image_level_value.max() * 1.15))
    for region, col in zip(regions, rcols): ax.scatter([], [], color=col, label=region)
    ax.legend(frameon=False, fontsize=5.6, ncol=3, loc="upper right")
    ax.text(.99, .02, "Image-level source values; donor IDs unavailable", transform=ax.transAxes, ha="right", fontsize=6, color="#66737C")

    stem = "Figure_S8_published_functional_protein_source_data"
    png = OUT / f"{stem}.png"; tif = OUT / f"{stem}.tiff"
    png_tmp = OUT / f".{stem}.tmp.png"; tif_tmp = OUT / f".{stem}.tmp.tiff"
    fig.savefig(png_tmp, dpi=350, bbox_inches="tight", format="png")
    fig.savefig(tif_tmp, dpi=600, bbox_inches="tight", format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    os.replace(png_tmp, png); os.replace(tif_tmp, tif)
    plt.close(fig)
    print(png)


if __name__ == "__main__":
    main()
