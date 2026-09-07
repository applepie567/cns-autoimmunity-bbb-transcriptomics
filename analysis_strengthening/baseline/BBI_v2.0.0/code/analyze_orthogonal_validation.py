"""Reanalyse published source data used as orthogonal checks in the BBI manuscript.

The script treats each mouse—not each microscopic field—as the independent unit.
Human immunohistochemistry values are summarized descriptively because the Prism
source file does not identify which image fields came from the same donor.
"""

from __future__ import annotations

import itertools
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
FIG7 = ROOT / "public_data/published_validation/Figure 7_b,c,e,f,g.pzfx"
EXT7 = ROOT / "public_data/published_validation/ED_Figure7_c-h.pzfx"
EXT6 = ROOT / "public_data/published_validation/ED_Figure6_a.pzfx"
FIG5 = ROOT / "public_data/published_validation/Figure 5_d-f.pzfx"


def lname(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def children(el: ET.Element, name: str) -> list[ET.Element]:
    return [node for node in el if lname(node.tag) == name]


def child(el: ET.Element, name: str) -> ET.Element | None:
    return next((node for node in el if lname(node.tag) == name), None)


def text(el: ET.Element | None) -> str:
    return "" if el is None else "".join(el.itertext()).strip()


def tables(path: Path) -> list[dict]:
    root = ET.parse(path).getroot()
    parsed = []
    for table in [node for node in root if lname(node.tag) == "Table"]:
        row_node = child(table, "RowTitlesColumn")
        row_titles: list[str] = []
        if row_node is not None:
            sub = children(row_node, "Subcolumn")
            if sub:
                row_titles = [text(d) for d in children(sub[0], "d")]
        cols = []
        for ycol in children(table, "YColumn"):
            values = []
            for sub in children(ycol, "Subcolumn"):
                values.append(
                    [float(text(d)) if text(d) not in {"", "nan"} else np.nan for d in children(sub, "d")]
                )
            cols.append({"title": text(child(ycol, "Title")), "subcolumns": values})
        parsed.append({"id": table.attrib.get("ID", ""), "title": text(child(table, "Title")), "rows": row_titles, "columns": cols})
    return parsed


def mouse_means(table: dict, endpoint: str, source: str) -> pd.DataFrame:
    rows = table["rows"]
    records = []
    for column in table["columns"]:
        if not column["title"].strip():
            continue
        mat = np.full((len(column["subcolumns"]), len(rows)), np.nan)
        for i, sub in enumerate(column["subcolumns"]):
            mat[i, : min(len(sub), len(rows))] = sub[: len(rows)]
        for j, mouse in enumerate(rows):
            vals = mat[:, j]
            vals = vals[np.isfinite(vals)]
            records.append(
                {
                    "source": source,
                    "endpoint": endpoint,
                    "group": "r84" if column["title"].strip().lower() == "r84" else column["title"].strip(),
                    "mouse": mouse.title(),
                    "n_fields": len(vals),
                    "mouse_mean": float(np.mean(vals)),
                }
            )
    return pd.DataFrame.from_records(records)


def exact_permutation_p(a: np.ndarray, b: np.ndarray) -> float:
    """Two-sided exact randomization P value for a difference in means."""
    pooled = np.concatenate([a, b])
    n_a = len(a)
    observed = abs(float(np.mean(a) - np.mean(b)))
    extreme = 0
    total = 0
    all_idx = set(range(len(pooled)))
    for idx_a in itertools.combinations(range(len(pooled)), n_a):
        idx_a = set(idx_a)
        idx_b = list(all_idx - idx_a)
        diff = abs(float(np.mean(pooled[list(idx_a)]) - np.mean(pooled[idx_b])))
        extreme += diff >= observed - 1e-12
        total += 1
    return extreme / total


def hedges_g(treat: np.ndarray, control: np.ndarray) -> tuple[float, float, float]:
    n1, n0 = len(treat), len(control)
    df = n1 + n0 - 2
    pooled_sd = math.sqrt(((n1 - 1) * np.var(treat, ddof=1) + (n0 - 1) * np.var(control, ddof=1)) / df)
    if pooled_sd == 0:
        return float("nan"), float("nan"), float("nan")
    d = (np.mean(treat) - np.mean(control)) / pooled_sd
    correction = 1 - 3 / (4 * (n1 + n0) - 9)
    g = correction * d
    var_g = (n1 + n0) / (n1 * n0) + g * g / (2 * df)
    se = math.sqrt(var_g)
    return g, g - 1.96 * se, g + 1.96 * se


def bh(values: pd.Series) -> pd.Series:
    x = values.to_numpy(float)
    order = np.argsort(x)
    ranked = x[order] * len(x) / np.arange(1, len(x) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.minimum(ranked, 1)
    return pd.Series(out, index=values.index)


def summarize_mouse_data(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for endpoint, block in data.groupby("endpoint", sort=False):
        control = block.loc[block.group.str.lower() == "igg", "mouse_mean"].to_numpy(float)
        treat = block.loc[block.group.str.lower() == "r84", "mouse_mean"].to_numpy(float)
        g, lo, hi = hedges_g(treat, control)
        rows.append(
            {
                "endpoint": endpoint,
                "n_IgG": len(control),
                "n_r84": len(treat),
                "IgG_mean": np.mean(control),
                "IgG_sd": np.std(control, ddof=1),
                "r84_mean": np.mean(treat),
                "r84_sd": np.std(treat, ddof=1),
                "r84_minus_IgG": np.mean(treat) - np.mean(control),
                "hedges_g": g,
                "g_ci_low": lo,
                "g_ci_high": hi,
                "p_exact_two_sided": exact_permutation_p(treat, control),
            }
        )
    out = pd.DataFrame(rows)
    permeability = out.endpoint.isin(["IgG leakage", "Fibrinogen leakage"])
    target = ~permeability
    out["FDR_within_family"] = np.nan
    out.loc[permeability, "FDR_within_family"] = bh(out.loc[permeability, "p_exact_two_sided"])
    out.loc[target, "FDR_within_family"] = bh(out.loc[target, "p_exact_two_sided"])
    out["family"] = np.where(permeability, "permeability", "vascular target engagement")
    return out


def protein_data(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = []
    for table in tables(path):
        marker = table["title"].replace("EGLF7", "EGFL7")
        for col in table["columns"]:
            values = [v for sub in col["subcolumns"] for v in sub if np.isfinite(v)]
            for i, value in enumerate(values, 1):
                records.append({"marker": marker, "region": col["title"], "observation": i, "image_level_value": value})
    raw = pd.DataFrame(records)
    summary = (
        raw.groupby(["marker", "region"], sort=False).image_level_value
        .agg(n_image_observations="size", mean="mean", sd="std", median="median", q1=lambda x: x.quantile(.25), q3=lambda x: x.quantile(.75))
        .reset_index()
    )
    summary["inference_note"] = "Descriptive only; source file has no donor identifiers for image-level observations"
    return raw, summary


def teer_data(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    table = next(t for t in tables(path) if t["title"] == "FinalFigure_VEGFaR84ResistanceDrop_Acute_Chronic")
    cols = [c for c in table["columns"] if c["title"].strip()]
    records = []
    for i, col in enumerate(cols):
        phase = "acute drop" if i < 3 else "chronic drop"
        condition = " ".join(col["title"].split()).replace("R84", "r84")
        values = [v for sub in col["subcolumns"] for v in sub if np.isfinite(v)]
        for well, value in enumerate(values, 1):
            records.append({"phase": phase, "condition": condition, "well": well, "TEER_change": value})
    raw = pd.DataFrame(records)
    summary = raw.groupby(["phase", "condition"], sort=False).TEER_change.agg(n_wells="size", mean="mean", sd="std").reset_index()
    summary["inference_note"] = "Technical-well summary; experiment-level replicate values were unavailable"
    return raw, summary


def main() -> None:
    fig7_names = {
        "Figure 7b-Chronic_IgG vs r84_Venous proliferation": "Venous endothelial proliferation",
        "Figure 7c-Chronic_IgG vs r84_Venous coverage": "Venous coverage",
        "Figure 7e-Chronic_IgG vs r84_Endomucin coverage": "Endomucin coverage",
        "Figure 7f-Chronic_IgG vs r84_IgG Leakage": "IgG leakage",
        "Figure 7g-Chronic_IgG vs r84_Venous width": "Venous width",
    }
    ext7_names = {
        "Number of ki67 pod+": "Ki67-positive endothelial cells",
        "precentage of fibrinogen area": "Fibrinogen leakage",
    }
    frames = []
    for table in tables(FIG7):
        if table["title"] in fig7_names:
            frames.append(mouse_means(table, fig7_names[table["title"]], "Shahriar et al. 2024, Figure 7"))
    for table in tables(EXT7):
        if table["title"] in ext7_names:
            frames.append(mouse_means(table, ext7_names[table["title"]], "Shahriar et al. 2024, Extended Data Figure 7"))
    mouse = pd.concat(frames, ignore_index=True)
    mouse.to_csv(OUT / "published_source_mouse_level_barrier_and_vascular_data.csv", index=False)
    summary = summarize_mouse_data(mouse)
    summary.to_csv(OUT / "published_source_barrier_and_vascular_summary.csv", index=False)

    ihc, ihc_summary = protein_data(FIG5)
    ihc.to_csv(OUT / "published_human_MS_IHC_image_level_source_data.csv", index=False)
    ihc_summary.to_csv(OUT / "published_human_MS_IHC_descriptive_summary.csv", index=False)

    teer, teer_summary = teer_data(EXT6)
    teer.to_csv(OUT / "published_mBEC_TEER_well_level_source_data.csv", index=False)
    teer_summary.to_csv(OUT / "published_mBEC_TEER_descriptive_summary.csv", index=False)

    print(summary.to_string(index=False))
    print("\nHuman IHC descriptive summary\n", ihc_summary.to_string(index=False))
    print("\nTEER descriptive summary\n", teer_summary.to_string(index=False))


if __name__ == "__main__":
    main()
