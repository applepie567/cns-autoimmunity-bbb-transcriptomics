#!/usr/bin/env python3
"""Reproduce Table S15 from frozen effects, or repeat the recovered HCOP filter.

The analysis functions are restored from public_data_strengthen_bbi.py in the
4 September 2026 public-data analysis package. The original focused-set
bootstrap seed is retained so the frozen confidence intervals are reproduced.
Document-editing and obsolete figure-generation code are intentionally omitted.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
BOOTSTRAP_SEED = 20260904


def fisher_ci(rho: float, n: int) -> tuple[float, float]:
    z = np.arctanh(rho)
    se = 1 / np.sqrt(n - 3)
    crit = stats.norm.ppf(0.975)
    return tuple(np.tanh([z - crit * se, z + crit * se]))


def bootstrap_spearman(x: np.ndarray, y: np.ndarray,
                       seed: int = BOOTSTRAP_SEED, n_boot: int = 2000) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(x)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.unique(x[idx]).size < 2 or np.unique(y[idx]).size < 2:
            continue
        vals.append(stats.spearmanr(x[idx], y[idx]).statistic)
    return tuple(np.quantile(vals, [0.025, 0.975]))


def hcop_symbols(path: Path) -> set[str]:
    """Apply the original reciprocal uniqueness, support and symbol filters.

    Uniqueness is evaluated before filtering support. The final equal-symbol
    restriction makes this a sensitivity subset of the primary symbol-matched
    comparison. It does not add pairs with different mouse and human symbols.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        hcop = pd.read_csv(handle, sep="\t", dtype=str)
    hcop.columns = [c.strip().lower().replace(" ", "_") for c in hcop.columns]
    human = next(c for c in hcop.columns if c in {"human_symbol", "human_gene_symbol"})
    mouse = next(c for c in hcop.columns if c in {"mouse_symbol", "mouse_gene_symbol"})
    support = next(c for c in hcop.columns if "support" in c)
    hcop = hcop.dropna(subset=[human, mouse, support]).copy()
    hcop[human] = hcop[human].str.upper()
    hcop[mouse] = hcop[mouse].str.upper()
    hcop = hcop[hcop[human].ne("-") & hcop[mouse].ne("-")
                & hcop[human].str.len().gt(0) & hcop[mouse].str.len().gt(0)].copy()
    human_counts = hcop.groupby(human)[mouse].nunique()
    mouse_counts = hcop.groupby(mouse)[human].nunique()
    pairs = hcop[hcop[human].map(human_counts).eq(1)
                 & hcop[mouse].map(mouse_counts).eq(1)
                 & hcop[support].str.contains("Ensembl", case=False, na=False)
                 & hcop[support].str.contains("NCBI", case=False, na=False)].copy()
    return set(pairs.loc[pairs[human].eq(pairs[mouse]), human])


def summarize(df: pd.DataFrame, mapping: str, scope: str) -> dict:
    x = df.mouse_pooled_g.to_numpy(float)
    y = df.human_hedges_g.to_numpy(float)
    rho, p = stats.spearmanr(x, y)
    if len(df) >= 50:
        lo, hi = fisher_ci(float(rho), len(df))
        method = "Fisher z"
    else:
        lo, hi = bootstrap_spearman(x, y)
        method = "bootstrap (2,000 resamples)"
    return dict(mapping=mapping, scope=scope, n_genes=len(df),
                spearman_rho=float(rho), ci_low=float(lo), ci_high=float(hi),
                p_value=float(p), same_direction_fraction=float(np.mean(np.sign(x) == np.sign(y))),
                ci_method=method)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hcop", type=Path,
                        help="Optional original HCOP export dated 2026-09-04. Without it, use the frozen selected gene set.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "reconstruction_results" / "strict_orthology")
    args = parser.parse_args()
    effects = pd.read_csv(ROOT / "results" / "cross_species_gene_effects.csv")
    archived = pd.read_csv(ROOT / "results" / "strict_orthology_gene_effects.csv")
    frozen_symbols = set(archived.gene.astype(str).str.upper())
    symbols = hcop_symbols(args.hcop) if args.hcop else frozen_symbols
    strict = effects[effects.gene.astype(str).str.upper().isin(symbols)].copy()
    if set(strict.gene.astype(str).str.upper()) != frozen_symbols:
        raise ValueError("The supplied HCOP export changes the frozen selected gene set. Use the original dated export.")
    pd.testing.assert_frame_equal(strict.reset_index(drop=True), archived.reset_index(drop=True),
                                  check_exact=False, rtol=1e-12, atol=1e-14)
    focused = effects.prespecified_state_gene.astype(str).str.lower().isin({"true", "1", "yes"})
    strict_focused = strict.prespecified_state_gene.astype(str).str.lower().isin({"true", "1", "yes"})
    summary = pd.DataFrame([
        summarize(effects, "Uppercase-symbol match", "All shared genes"),
        summarize(strict, "HCOP reciprocal 1:1; Ensembl+NCBI support", "All shared genes"),
        summarize(effects[focused], "Uppercase-symbol match", "Prespecified focused genes"),
        summarize(strict[strict_focused], "HCOP reciprocal 1:1; Ensembl+NCBI support", "Prespecified focused genes"),
    ])
    expected = pd.read_csv(ROOT / "results" / "strict_orthology_sensitivity_summary.csv")
    pd.testing.assert_frame_equal(summary, expected, check_exact=False, rtol=1e-12, atol=1e-14)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_dir / "strict_orthology_sensitivity_summary.csv", index=False)
    strict.to_csv(args.output_dir / "strict_orthology_gene_effects.csv", index=False)
    print(summary.to_string(index=False))
    print("All Table S15 values match the frozen summary. Bootstrap seed:", BOOTSTRAP_SEED)
    print("Selection source:", str(args.hcop) if args.hcop else "archived selected gene effects; original HCOP filtering not rerun")


if __name__ == "__main__":
    main()
