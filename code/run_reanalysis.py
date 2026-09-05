#!/usr/bin/env python3
"""Deterministic validation reanalysis for the MS/BBB manuscript.

The script re-analyzes uploaded count matrices where matrix-readable inputs are
available, audits uploaded Space Ranger archives, and exports tidy CSV files for
figures and the source-data workbook. All stochastic operations use SEED.
"""

from __future__ import annotations

import gc
import gzip
import itertools
import json
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import mannwhitneyu
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
UPLOAD = ROOT / "public_data"
OUT = ROOT / "reconstruction_results"
SEED = 20260901
RNG = np.random.default_rng(SEED)


BBB_MODULES = {
    "Tight_junction": ["CLDN5", "OCLN", "TJP1", "LSR"],
    "BBB_transport_identity": ["MFSD2A", "SLC2A1", "ABCG2", "ABCB1"],
    "Caveolae_structural": ["CAV1", "CAV2", "EHD2"],
    "PLVAP_permeability": ["PLVAP", "ESM1", "APLN", "ANGPT2"],
    "Adhesion_trafficking": ["ICAM1", "VCAM1", "SELE", "SELP", "CCL2", "CXCL10"],
    "IFN_antigen_presentation": [
        "STAT1", "IRF1", "CXCL9", "CXCL10", "B2M", "TAP1",
        "HLA-DRA", "HLA-DRB1", "IFIT1", "IFIT2", "IFIT3", "ISG15",
    ],
    "Wnt_BBB": ["LEF1", "AXIN2", "APCDD1", "ADGRA2", "RECK", "TCF7L2"],
    "TGF_response": ["TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7", "SERPINE1", "CTGF"],
    "VEGF_response": ["KDR", "FLT1", "PGF", "ANGPT2", "ESM1", "APLN"],
    "ROS_Src": ["SRC", "RAC1", "NOX4", "CYBB", "NCF1", "NCF2", "HIF1A"],
    "ECM_protease": ["MMP2", "MMP9", "TIMP1", "TIMP2"],
}

EC_MARKERS = ["PECAM1", "VWF", "EMCN", "CDH5", "RAMP2", "RGCC", "CA4", "EPAS1", "TEK", "ENG", "ESAM", "KLF2", "KLF4"]

IMMUNE_MODULES = {
    "IFN_antiviral": ["IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "OAS2", "STAT1", "IRF7"],
    "TNF_NFkB": ["TNF", "NFKB1", "NFKBIA", "RELA", "TNFAIP3", "ICAM1"],
    "IL1_IL6": ["IL1B", "IL1R1", "IL6", "IL6R", "STAT3", "SOCS3"],
    "TGF_remodeling": ["TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SERPINE1", "CTGF"],
    "Myeloid_proteolytic": ["S100A8", "S100A9", "CTSS", "CTSB", "LST1", "FCN1", "MMP9"],
    "Chemokine": ["CCL2", "CCL3", "CCL4", "CXCL8", "CXCL9", "CXCL10"],
    "VEGF_angiogenic": ["VEGFA", "FLT1", "KDR", "PGF", "ANGPT2", "HIF1A"],
}

CELL_MARKERS = {
    "Neutrophil": ["S100A8", "S100A9", "CSF3R", "FCGR3B", "CXCR2"],
    "Monocyte": ["LST1", "FCN1", "CTSS", "LYZ", "S100A10"],
    "Macrophage": ["C1QA", "C1QB", "C1QC", "APOE", "CTSD"],
    "mDC": ["FCER1A", "CD1C", "CST3", "CLEC10A"],
    "pDC": ["GZMB", "JCHAIN", "TCF4", "IL3RA"],
    "T": ["CD3D", "CD3E", "TRAC", "IL7R"],
    "NK": ["NKG7", "GNLY", "PRF1", "KLRD1"],
    "B": ["CD79A", "MS4A1", "CD74", "CD37"],
    "Plasma": ["MZB1", "JCHAIN", "SDC1", "IGHG1"],
}


def ensure_dirs() -> None:
    for p in [OUT, OUT / "source_tables", OUT / "gse163005", OUT / "gse163194", OUT / "gse279183", OUT / "gse208747", OUT / "representative_images"]:
        p.mkdir(parents=True, exist_ok=True)


def bh_fdr(pvalues: list[float] | np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    inv = np.empty_like(order)
    inv[order] = np.arange(len(order))
    out[ok] = np.minimum(adj[inv], 1.0)
    return out


def exact_median_permutation(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Two-sided exact label permutation for difference in medians."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    obs = float(np.median(x) - np.median(y))
    z = np.r_[x, y]
    n1 = len(x)
    total = 0
    extreme = 0
    for idx in itertools.combinations(range(len(z)), n1):
        mask = np.zeros(len(z), dtype=bool)
        mask[list(idx)] = True
        stat = np.median(z[mask]) - np.median(z[~mask])
        total += 1
        if abs(stat) >= abs(obs) - 1e-12:
            extreme += 1
    return obs, extreme / total


def zscore_rows(a: np.ndarray) -> np.ndarray:
    mu = np.nanmean(a, axis=1, keepdims=True)
    sd = np.nanstd(a, axis=1, ddof=1, keepdims=True)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    return (a - mu) / sd


def module_scores_from_gene_by_sample(logexpr: np.ndarray, genes: list[str], modules: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes_u = np.array([str(g).upper() for g in genes])
    lookup = {g: i for i, g in enumerate(genes_u)}
    z = zscore_rows(logexpr)
    scores = {}
    coverage = []
    for module, members in modules.items():
        idx = [lookup[g.upper()] for g in members if g.upper() in lookup]
        scores[module] = np.nanmean(z[idx, :], axis=0) if idx else np.full(logexpr.shape[1], np.nan)
        coverage.append({"module": module, "genes_expected": len(members), "genes_found": len(idx), "genes_used": "; ".join([genes_u[i] for i in idx])})
    return pd.DataFrame(scores), pd.DataFrame(coverage)


def export_source_tables() -> None:
    import openpyxl
    import shutil

    workbook = UPLOAD / "MS_BBB_complete_source_data_tables.xlsx"
    destination = OUT / "source_tables"
    destination.mkdir(parents=True, exist_ok=True)
    existing = sorted((OUT / "source_tables").glob("S*.csv"))
    if not workbook.exists():
        if existing:
            print(
                "Source tables are already present in reconstruction_results/source_tables; "
                "skipping workbook extraction."
            )
            return
        bundled = sorted((ROOT / "supplementary_tables").glob("S*.csv"))
        if bundled:
            for source in bundled:
                shutil.copyfile(source, destination / source.name)
            print("Copied the bundled supplementary CSV tables to reconstruction_results/source_tables.")
            return
        raise FileNotFoundError(
            "Source tables are unavailable. Keep the repository-provided "
            "supplementary_tables directory or place "
            "MS_BBB_complete_source_data_tables.xlsx in public_data/."
        )

    wb = openpyxl.load_workbook(workbook, read_only=True, data_only=True)
    for sheet in wb.sheetnames:
        rows = list(wb[sheet].values)
        if not rows:
            continue
        pd.DataFrame(rows[1:], columns=rows[0]).to_csv(OUT / "source_tables" / f"{sheet}.csv", index=False)


def analyze_gse163005() -> None:
    ddir = OUT / "gse163005"
    barcodes = pd.read_csv(UPLOAD / "GSE163005_barcodes.tsv.gz", sep="\t", header=None)[0].astype(str)
    features = pd.read_csv(UPLOAD / "GSE163005_features.tsv.gz", sep="\t", header=None, names=["id", "gene", "type"])
    dx = pd.read_csv(UPLOAD / "GSE163005_annotation_dx.csv.gz", index_col=0).iloc[:, 0].astype(str)
    patient = pd.read_csv(UPLOAD / "GSE163005_annotation_patients.csv.gz", index_col=0).iloc[:, 0].astype(str)
    cluster = pd.read_csv(UPLOAD / "GSE163005_annotation_cluster.csv.gz", index_col=0).iloc[:, 0].astype(str)

    with gzip.open(UPLOAD / "GSE163005_matrix.mtx.gz", "rb") as fh:
        mat = mmread(fh).tocsc()
    if mat.shape != (len(features), len(barcodes)):
        raise ValueError(f"GSE163005 matrix mismatch: {mat.shape}, features={len(features)}, barcodes={len(barcodes)}")

    pos = pd.Series(np.arange(len(barcodes)), index=barcodes.values)
    common = dx.index.intersection(patient.index).intersection(cluster.index).intersection(pos.index)
    idx = pos.loc[common].to_numpy()
    sub = mat[:, idx]
    totals = np.asarray(sub.sum(axis=0)).ravel()
    genes_detected = np.asarray(sub.getnnz(axis=0)).ravel()
    mt_idx = np.where(features["gene"].str.upper().str.startswith("MT-").to_numpy())[0]
    mt = np.asarray(sub[mt_idx, :].sum(axis=0)).ravel()
    qc = pd.DataFrame({
        "barcode": common,
        "diagnosis": dx.loc[common].values,
        "patient": patient.loc[common].values,
        "cluster": cluster.loc[common].values,
        "total_UMI": totals,
        "genes_detected": genes_detected,
        "mitochondrial_fraction": np.divide(mt, totals, out=np.zeros_like(mt, dtype=float), where=totals > 0),
    })
    qc.to_csv(ddir / "cell_qc.csv.gz", index=False, compression="gzip")

    selected = qc[qc["diagnosis"].isin(["MS", "IIH", "VE"])].copy()
    patients = selected[["patient", "diagnosis"]].drop_duplicates().sort_values(["diagnosis", "patient"]).reset_index(drop=True)
    pmap = {p: i for i, p in enumerate(patients["patient"])}
    sel_cols = pos.loc[selected["barcode"]].to_numpy()
    pcols = selected["patient"].map(pmap).to_numpy()
    indicator = sparse.csr_matrix((np.ones(len(sel_cols)), (np.arange(len(sel_cols)), pcols)), shape=(len(sel_cols), len(patients)))
    pseudobulk = mat[:, sel_cols] @ indicator
    libsize = np.asarray(pseudobulk.sum(axis=0)).ravel()
    logcpm = np.log1p(np.asarray(pseudobulk.toarray(), dtype=np.float64) / np.maximum(libsize, 1)[None, :] * 1e6)
    scores, coverage = module_scores_from_gene_by_sample(logcpm, features["gene"].tolist(), IMMUNE_MODULES)
    scores.insert(0, "diagnosis", patients["diagnosis"].values)
    scores.insert(0, "patient", patients["patient"].values)
    scores.to_csv(ddir / "patient_module_scores.csv", index=False)
    coverage.to_csv(ddir / "module_gene_coverage.csv", index=False)

    contrasts = []
    for label, g1, g2 in [("MS vs IIH", "MS", "IIH"), ("VE vs IIH", "VE", "IIH"), ("VE vs MS", "VE", "MS")]:
        block = []
        for module in IMMUNE_MODULES:
            x = scores.loc[scores["diagnosis"] == g1, module].to_numpy(float)
            y = scores.loc[scores["diagnosis"] == g2, module].to_numpy(float)
            effect, p = exact_median_permutation(x, y)
            block.append({"module": module, "contrast": label, "delta_median": effect, "p_exact": p, "n1": len(x), "n2": len(y)})
        fdr = bh_fdr([r["p_exact"] for r in block])
        for row, q in zip(block, fdr):
            row["FDR_BH"] = q
        contrasts.extend(block)
    pd.DataFrame(contrasts).to_csv(ddir / "raw_reconstructed_module_contrasts.csv", index=False)

    broad = cluster.replace({
        "CD4": "T", "CD8": "T", "Treg": "T", "cycling": "T",
        "mono1": "Monocyte", "mono2": "Monocyte", "mono3": "Monocyte",
        "mDC1": "mDC", "mDC2": "mDC", "matDC": "mDC",
        "granulo1": "Granulocyte", "granulo2": "Granulocyte",
        "naiveBc": "B", "plasma": "Plasma", "pDC": "pDC", "NK": "NK",
    })
    comp = pd.DataFrame({"patient": patient.loc[common], "diagnosis": dx.loc[common], "celltype": broad.loc[common]}).reset_index(drop=True)
    comp = comp[comp["diagnosis"].isin(["MS", "IIH", "VE"])]
    counts = comp.groupby(["patient", "diagnosis", "celltype"]).size().rename("cells").reset_index()
    counts["fraction"] = counts["cells"] / counts.groupby("patient")["cells"].transform("sum")
    counts.to_csv(ddir / "patient_cell_composition.csv", index=False)

    summary = qc.groupby("diagnosis").agg(cells=("barcode", "size"), patients=("patient", "nunique"), median_UMI=("total_UMI", "median"), median_genes=("genes_detected", "median"), median_mito_fraction=("mitochondrial_fraction", "median")).reset_index()
    summary.to_csv(ddir / "qc_summary.csv", index=False)
    del logcpm, pseudobulk, indicator, sub, mat
    gc.collect()


GSE163194_STAGE = {
    "C47": "S6", "C55": "S7", "C56": "S6", "C57": "S7", "C58": "S6", "C59": "S6",
    "C60": "S7", "C61": "S7", "C62": "S8", "C65": "S3", "C66": "S7", "C67": "S7",
    "C69": "S8", "C71": "S7", "C72": "S7", "C73": "S8", "C77": "S4", "C80": "S8",
    "C96": "S3", "C100": "S1", "C101": "S3", "C102": "S4", "C103": "S1", "C106": "S3",
    "C107": "S3", "C114": "S7", "C119": "S4", "C121": "S1", "C122": "S3", "C126": "S3",
    "C129": "S8", "C131": "S1", "C139": "S4",
}


def analyze_gse163194() -> None:
    ddir = OUT / "gse163194"
    targets = sorted(set(sum(IMMUNE_MODULES.values(), []) + sum(CELL_MARKERS.values(), [])))
    target_set = set(targets)
    patient_gene_counts: dict[str, pd.Series] = {}
    patient_composition: dict[str, dict[str, int]] = {}
    library_qc = []

    with tarfile.open(UPLOAD / "GSE163194_RAW.tar", "r") as tar:
        for member in tar.getmembers():
            if not member.name.endswith(".csv.gz"):
                continue
            sample_match = re.search(r"_(C\d+)(?:-\d+)?-S", member.name)
            if not sample_match:
                raise ValueError(member.name)
            patient = sample_match.group(1)
            raw = tar.extractfile(member)
            assert raw is not None
            with gzip.GzipFile(fileobj=raw, mode="rb") as gz:
                chunks = pd.read_csv(gz, chunksize=2000)
                totals = None
                genes_detected = None
                target_by_cell = None
                cell_names = None
                for chunk in chunks:
                    gene_col = chunk.columns[0]
                    if cell_names is None:
                        cell_names = list(chunk.columns[1:])
                        totals = np.zeros(len(cell_names), dtype=np.float64)
                        genes_detected = np.zeros(len(cell_names), dtype=np.int32)
                        target_by_cell = pd.DataFrame(0.0, index=targets, columns=cell_names)
                    arr = chunk.iloc[:, 1:].to_numpy(dtype=np.float32, copy=False)
                    totals += arr.sum(axis=0)
                    genes_detected += (arr > 0).sum(axis=0)
                    genes = chunk[gene_col].astype(str).str.upper().to_numpy()
                    keep = np.where(np.isin(genes, list(target_set)))[0]
                    for i in keep:
                        target_by_cell.loc[genes[i], :] += arr[i, :]
            assert target_by_cell is not None and totals is not None and genes_detected is not None
            gene_sum = target_by_cell.sum(axis=1)
            patient_gene_counts[patient] = patient_gene_counts.get(patient, pd.Series(0.0, index=targets)).add(gene_sum, fill_value=0)

            logt = np.log1p(target_by_cell)
            score_mat = []
            labels = []
            for ct, genes in CELL_MARKERS.items():
                present = [g for g in genes if g in logt.index]
                score_mat.append(logt.loc[present].mean(axis=0).to_numpy() if present else np.zeros(logt.shape[1]))
                labels.append(ct)
            score_mat = np.vstack(score_mat)
            max_score = score_mat.max(axis=0)
            assignment = np.array(labels, dtype=object)[score_mat.argmax(axis=0)]
            assignment[max_score <= 0] = "Other"
            counts = pd.Series(assignment).value_counts().to_dict()
            pc = patient_composition.setdefault(patient, {})
            for ct, n in counts.items():
                pc[ct] = pc.get(ct, 0) + int(n)
            library_qc.append({
                "library": member.name.split("/")[-1].replace(".csv.gz", ""),
                "patient": patient,
                "stage": GSE163194_STAGE[patient],
                "cells": len(totals),
                "median_UMI": float(np.median(totals)),
                "median_genes": float(np.median(genes_detected)),
            })
            del target_by_cell, logt, score_mat
            gc.collect()

    patients = sorted(patient_gene_counts)
    genes = targets
    counts = np.column_stack([patient_gene_counts[p].reindex(genes, fill_value=0).to_numpy(float) for p in patients])
    libsize = counts.sum(axis=0)
    logcpm = np.log1p(counts / np.maximum(libsize, 1)[None, :] * 1e6)
    scores, coverage = module_scores_from_gene_by_sample(logcpm, genes, IMMUNE_MODULES)
    meta = pd.DataFrame({"patient": patients, "stage": [GSE163194_STAGE[p] for p in patients]})
    meta["group"] = np.where(meta["stage"].eq("S1"), "Acute", np.where(meta["stage"].isin(["S4", "S8"]), "Recovery", "Other"))
    scores.insert(0, "group", meta["group"])
    scores.insert(0, "stage", meta["stage"])
    scores.insert(0, "patient", meta["patient"])
    scores.to_csv(ddir / "patient_module_scores.csv", index=False)
    coverage.to_csv(ddir / "module_gene_coverage.csv", index=False)

    contrast = []
    for module in IMMUNE_MODULES:
        x = scores.loc[scores["group"] == "Acute", module].to_numpy(float)
        y = scores.loc[scores["group"] == "Recovery", module].to_numpy(float)
        effect = float(np.median(x) - np.median(y))
        p = float(mannwhitneyu(x, y, alternative="two-sided", method="exact").pvalue)
        cliff = float((np.sum(x[:, None] > y[None, :]) - np.sum(x[:, None] < y[None, :])) / (len(x) * len(y)))
        contrast.append({"module": module, "contrast": "Acute_vs_All_recovery", "delta_median": effect, "p_MWU": p, "n1": len(x), "n2": len(y), "cliffs_delta": cliff})
    fdr = bh_fdr([r["p_MWU"] for r in contrast])
    for row, q in zip(contrast, fdr):
        row["FDR_BH"] = q
    pd.DataFrame(contrast).to_csv(ddir / "raw_reconstructed_module_contrasts.csv", index=False)

    comp_rows = []
    for patient, counts_d in patient_composition.items():
        total = sum(counts_d.values())
        for ct, n in counts_d.items():
            comp_rows.append({"patient": patient, "stage": GSE163194_STAGE[patient], "group": "Acute" if GSE163194_STAGE[patient] == "S1" else ("Recovery" if GSE163194_STAGE[patient] in ["S4", "S8"] else "Other"), "celltype": ct, "cells": n, "fraction": n / total})
    pd.DataFrame(comp_rows).to_csv(ddir / "marker_reconstructed_composition.csv", index=False)
    pd.DataFrame(library_qc).to_csv(ddir / "library_qc.csv", index=False)
    meta.to_csv(ddir / "sample_stage_map.csv", index=False)


def read_mtx_triplet(prefix: str) -> tuple[sparse.csc_matrix, pd.DataFrame, pd.Series]:
    with gzip.open(prefix + "matrix", "rb") as fh:
        mat = mmread(fh).tocsc()
    features = pd.read_csv(prefix + "features", sep="\t", header=None, names=["id", "gene", "type"])
    barcodes = pd.read_csv(prefix + "barcodes", sep="\t", header=None)[0].astype(str)
    return mat, features, barcodes


def analyze_gse279183() -> None:
    ddir = OUT / "gse279183"
    coverage_rows = []

    def uploaded(prefix: str, stem: str) -> Path:
        """Resolve GEO filenames with or without browser-added '(1)' suffixes."""
        exact = UPLOAD / f"{prefix}_{stem}.gz"
        if exact.exists():
            return exact
        candidates = sorted(UPLOAD.glob(f"{prefix}_{stem}*.gz"))
        if len(candidates) == 1:
            return candidates[0]
        if not candidates:
            raise FileNotFoundError(
                f"Missing {prefix}_{stem}.gz in {UPLOAD}. Download the associated "
                "GEO supplementary file and keep the GSM/sample prefix unchanged."
            )
        raise RuntimeError(
            f"Ambiguous files for {prefix}_{stem}: "
            + ", ".join(path.name for path in candidates)
        )

    specs = [
        (
            "MS497T",
            uploaded("GSM8563712_MS497T", "matrix.mtx"),
            uploaded("GSM8563712_MS497T", "features.tsv"),
            uploaded("GSM8563712_MS497T", "barcodes.tsv"),
            uploaded("GSM8563712_MS497T", "tissue_positions_list.csv"),
            uploaded("GSM8563712_MS497T", "scalefactors_json.json"),
            uploaded("GSM8563712_MS497T", "tissue_lowres_image.png"),
        ),
        (
            "MS549H",
            uploaded("GSM8563713_MS549H", "matrix.mtx"),
            uploaded("GSM8563713_MS549H", "features.tsv"),
            uploaded("GSM8563713_MS549H", "barcodes.tsv"),
            uploaded("GSM8563713_MS549H", "tissue_positions_list.csv"),
            uploaded("GSM8563713_MS549H", "scalefactors_json.json"),
            uploaded("GSM8563713_MS549H", "tissue_lowres_image.png"),
        ),
        (
            "MS549T",
            uploaded("GSM8563714_MS549T", "matrix.mtx"),
            uploaded("GSM8563714_MS549T", "features.tsv"),
            uploaded("GSM8563714_MS549T", "barcodes.tsv"),
            uploaded("GSM8563714_MS549T", "tissue_positions_list.csv"),
            uploaded("GSM8563714_MS549T", "scalefactors_json.json"),
            uploaded("GSM8563714_MS549T", "tissue_lowres_image.png"),
        ),
    ]
    all_spots = []
    summaries = []
    for sample, mtxf, featf, barf, posf, scalef, imgf in specs:
        with gzip.open(mtxf, "rb") as fh:
            mat = mmread(fh).tocsc()
        features = pd.read_csv(featf, sep="\t", header=None, names=["id", "gene", "type"])
        barcodes = pd.read_csv(barf, sep="\t", header=None)[0].astype(str)
        pos = pd.read_csv(posf)
        if mat.shape != (len(features), len(barcodes)):
            raise ValueError(f"{sample} matrix mismatch {mat.shape}")
        pos = pos.set_index("barcode").reindex(barcodes)
        tissue = pos["in_tissue"].fillna(0).astype(int).to_numpy() == 1
        totals = np.asarray(mat.sum(axis=0)).ravel()
        detected = np.asarray(mat.getnnz(axis=0)).ravel()
        genes_u = features["gene"].astype(str).str.upper().to_numpy()
        lookup = {g: i for i, g in enumerate(genes_u)}
        target_genes = sorted(set(EC_MARKERS + sum(BBB_MODULES.values(), [])))
        target_idx = [lookup[g] for g in target_genes if g in lookup]
        target_names = [genes_u[i] for i in target_idx]
        norm = np.log1p(np.asarray(mat[target_idx, :].toarray(), dtype=float) / np.maximum(totals, 1)[None, :] * 1e4)
        norm_df = pd.DataFrame(norm, index=target_names)
        ec_present = [g for g in EC_MARKERS if g in norm_df.index]
        spot = pd.DataFrame({
            "sample": sample,
            "barcode": barcodes,
            "in_tissue": tissue,
            "array_row": pos["array_row"].to_numpy(),
            "array_col": pos["array_col"].to_numpy(),
            "pxl_row_in_fullres": pos["pxl_row_in_fullres"].to_numpy(),
            "pxl_col_in_fullres": pos["pxl_col_in_fullres"].to_numpy(),
            "total_UMI": totals,
            "genes_detected": detected,
            "EC_score": norm_df.loc[ec_present].mean(axis=0).to_numpy(),
        })
        for module, members in BBB_MODULES.items():
            present = [g for g in members if g in norm_df.index]
            spot[module] = norm_df.loc[present].mean(axis=0).to_numpy() if present else np.nan
            coverage_rows.append({"sample": sample, "module": module, "genes_expected": len(members), "genes_found": len(present), "genes_used": "; ".join(present)})
        all_spots.append(spot[spot["in_tissue"]].copy())
        summaries.append({"sample": sample, "matrix_genes": mat.shape[0], "matrix_barcodes": mat.shape[1], "spots_under_tissue": int(tissue.sum()), "median_UMI": float(np.median(totals[tissue])), "median_genes": float(np.median(detected[tissue]))})
        with gzip.open(scalef, "rt") as fh:
            scale = json.load(fh)
        with gzip.open(imgf, "rb") as fh:
            img = Image.open(fh).convert("RGB")
            img.save(OUT / "representative_images" / f"{sample}_tissue_lowres.png")
        with open(ddir / f"{sample}_scalefactors.json", "w", encoding="utf-8") as fh:
            json.dump(scale, fh, indent=2)
        del mat, norm, norm_df
        gc.collect()
    pd.concat(all_spots, ignore_index=True).to_csv(ddir / "spot_qc_and_module_scores.csv.gz", index=False, compression="gzip")
    pd.DataFrame(summaries).to_csv(ddir / "sample_qc_summary.csv", index=False)
    pd.DataFrame(coverage_rows).to_csv(ddir / "module_gene_coverage.csv", index=False)


def parse_metric(html: str, name: str) -> float:
    matches = re.findall(r'\["' + re.escape(name) + r'",\s*"([0-9,.]+)%?"\]', html)
    if not matches:
        return np.nan
    return float(matches[-1].replace(",", ""))


def analyze_gse208747() -> None:
    ddir = OUT / "gse208747"
    expected = ["A1", "A2", "A3", "A4", "C1", "C2", "M1", "M2", "M3", "M4", "M5", "M6", "N1", "N2", "N3"]
    files = sorted(UPLOAD.glob("GSM636*.tar.gz"))
    rows = []
    for f in files:
        sample = re.search(r"_([ACMN]\d+)\.tar\.gz$", f.name).group(1)
        group = {"A": "Active lesion", "C": "Control white matter", "M": "Mixed active/inactive lesion", "N": "NAWM"}[sample[0]]
        with tarfile.open(f, "r:gz") as tar:
            web = tar.extractfile(f"{sample}/web_summary.html")
            assert web is not None
            html = web.read().decode("utf-8", errors="replace")
            names = set(tar.getnames())
            required = [f"{sample}/filtered_feature_bc_matrix.h5", f"{sample}/spatial/tissue_positions_list.csv", f"{sample}/spatial/scalefactors_json.json", f"{sample}/spatial/tissue_lowres_image.png", f"{sample}/spatial/tissue_hires_image.png"]
            posf = tar.extractfile(f"{sample}/spatial/tissue_positions_list.csv")
            assert posf is not None
            positions = pd.read_csv(posf, header=None)
            spots_under_tissue = int(pd.to_numeric(positions.iloc[:, 1], errors="coerce").fillna(0).sum())
            rows.append({
                "sample": sample,
                "group": group,
                "archive": f.name,
                "core_files_complete": all(x in names for x in required),
                "spots_under_tissue": spots_under_tissue,
                "mean_reads_per_spot": parse_metric(html, "Mean Reads per Spot"),
                "median_genes_per_spot": parse_metric(html, "Median Genes per Spot"),
                "median_UMI_per_spot": parse_metric(html, "Median UMI Counts per Spot"),
                "sequencing_saturation_pct": parse_metric(html, "Sequencing Saturation"),
                "fraction_reads_under_tissue_pct": parse_metric(html, "Fraction Reads in Spots Under Tissue"),
            })
            if sample in {"A1", "C1", "M1", "N1"}:
                imf = tar.extractfile(f"{sample}/spatial/tissue_lowres_image.png")
                assert imf is not None
                Image.open(imf).convert("RGB").save(OUT / "representative_images" / f"GSE208747_{sample}_tissue_lowres.png")
    pd.DataFrame(rows).sort_values("sample").to_csv(ddir / "uploaded_archive_qc.csv", index=False)
    manifest = pd.DataFrame({"sample": expected})
    manifest["uploaded"] = manifest["sample"].isin([r["sample"] for r in rows])
    manifest["note"] = np.where(manifest["uploaded"], "available", "not uploaded; obtain from GEO before full raw rerun")
    manifest.to_csv(ddir / "expected_archive_manifest.csv", index=False)


def compare_reconstructed_with_frozen() -> None:
    rows = []
    checks = [
        ("GSE163005", OUT / "gse163005/raw_reconstructed_module_contrasts.csv", OUT / "source_tables/S11_GSE163005.csv", "contrast", "delta_median"),
        ("GSE163194", OUT / "gse163194/raw_reconstructed_module_contrasts.csv", OUT / "source_tables/S12_GSE163194.csv", "contrast", "delta_median"),
    ]
    for dataset, reconf, frozenf, contrast_col, effect_col in checks:
        recon = pd.read_csv(reconf)
        frozen = pd.read_csv(frozenf)
        if dataset == "GSE163194":
            frozen = frozen[(frozen["celltype"] == "WholeCSF") & (frozen["contrast"] == "Acute_vs_All_recovery")]
        merged = recon.merge(frozen[["module", contrast_col, effect_col]], on=["module", contrast_col], suffixes=("_reconstructed", "_frozen"))
        for _, r in merged.iterrows():
            a = float(r[f"{effect_col}_reconstructed"])
            b = float(r[f"{effect_col}_frozen"])
            rows.append({"dataset": dataset, "module": r["module"], "contrast": r[contrast_col], "effect_reconstructed": a, "effect_frozen": b, "direction_concordant": bool(np.sign(a) == np.sign(b))})
    pd.DataFrame(rows).to_csv(OUT / "raw_vs_frozen_direction_audit.csv", index=False)


def main() -> None:
    ensure_dirs()
    with open(OUT / "analysis_parameters.json", "w", encoding="utf-8") as fh:
        json.dump({
            "random_seed": SEED,
            "exact_permutations": "fully enumerated and therefore deterministic",
            "gse208747_spot_filters": {"minimum_UMI": 500, "minimum_genes": 200},
            "gse208747_EC_rich_primary_quantile": 0.30,
            "gse208747_EC_rich_sensitivity_quantiles": [0.20, 0.40],
            "notes": "No stochastic embedding is used in the final figure set. Canonical marker reconstruction for GSE163194 is an independent validation, not an exact reproduction of the historical cell labels.",
        }, fh, indent=2)
    export_source_tables()
    analyze_gse163005()
    analyze_gse163194()
    analyze_gse279183()
    analyze_gse208747()
    compare_reconstructed_with_frozen()
    print(json.dumps({"status": "completed", "seed": SEED, "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
