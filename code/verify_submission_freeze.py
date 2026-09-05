#!/usr/bin/env python3
"""Verify package completeness and agreement with the frozen manuscript."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd
from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[1]
TITLE = "Brain endothelial responses in acute EAE partially overlap with vascular changes in multiple sclerosis"
SEED = 20260901


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def docx_text_and_media_hashes(path: Path) -> tuple[str, set[str]]:
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
        text = "\n".join(node.text or "" for node in root.iter() if node.tag.endswith("}t"))
        hashes = {
            hashlib.sha256(archive.read(name)).hexdigest()
            for name in archive.namelist()
            if name.startswith("word/media/")
        }
    return text, hashes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manuscript", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    required = [
        "README.md",
        "CITATION.cff",
        "environment.yml",
        "requirements.txt",
        "software_versions.csv",
        "analysis_parameters.json",
        "results/analysis_summary.json",
        "supplementary_tables/Supplementary_Tables_S1-S18.xlsx",
        "figures/BBI_Supplementary_Figures_S1-S9.pdf",
        "figure_source_data/figure_source_index.csv",
    ]
    missing = [item for item in required if not (ROOT / item).is_file()]
    check("required files", not missing, "missing: " + ", ".join(missing) if missing else "all required files present")

    main_figures = sorted((ROOT / "figures" / "main").glob("Figure_*.png"))
    supplementary_figures = sorted((ROOT / "figures" / "supplementary").glob("Figure_S*.png"))
    check("main figure count", len(main_figures) == 5, f"found {len(main_figures)}")
    check("supplementary figure count", len(supplementary_figures) == 9, f"found {len(supplementary_figures)}")

    table_files = list((ROOT / "supplementary_tables").glob("S*.csv"))
    represented = set()
    for path in table_files:
        token = path.name.split("_", 1)[0]
        if token.startswith("S") and token[1:].isdigit():
            represented.add(int(token[1:]))
    check("supplementary tables S1 to S18", represented == set(range(1, 19)), f"represented tables: {sorted(represented)}")

    workbook_path = ROOT / "supplementary_tables" / "Supplementary_Tables_S1-S18.xlsx"
    workbook = load_workbook(workbook_path, read_only=True, data_only=False)
    expected_sheets = {
        "README", "S1_Datasets", "S2_Signatures", "S3_EAE_modules", "S4_EAE_genes",
        "S5_Mouse_sender", "S6_EAE_time", "S7_VEGFA_venous", "S8_MS_snRNA",
        "S9_MS_spatial", "S10_Lesion_stage", "S11_Viral_context", "S12_Bacterial_context",
        "S13_Infection_context", "S14_Result_audit", "S15_Orthology", "S16_MERFISH",
        "S17_Eligibility", "S18_Gene_sets", "S18_Coverage",
    }
    check("combined workbook sheets", set(workbook.sheetnames) == expected_sheets, f"found {len(workbook.sheetnames)} sheets")
    formula_errors = []
    for sheet in workbook.worksheets:
        for row in sheet.iter_rows():
            for cell in row:
                if isinstance(cell.value, str) and cell.value in {"#REF!", "#DIV/0!", "#VALUE!", "#NAME?", "#N/A"}:
                    formula_errors.append(f"{sheet.title}!{cell.coordinate}={cell.value}")
    check("workbook formula errors", not formula_errors, ", ".join(formula_errors) if formula_errors else "none")

    parameters = json.loads((ROOT / "analysis_parameters.json").read_text(encoding="utf-8"))
    check("random seed", parameters.get("random_seed") == SEED, str(parameters.get("random_seed")))
    figure_script = (ROOT / "code" / "make_bbi_figures.py").read_text(encoding="utf-8")
    check("selected Figure 4 background", "background_color: str = '#BFC4C8'" in figure_script and "background_alpha: float = .48" in figure_script, "light gray settings present")

    stale_names = [
        str(path.relative_to(ROOT))
        for path in ROOT.rglob("*")
        if path.is_file() and any(token in path.name.lower() for token in ["candidate", "figure_s10", "ms_gwas"])
    ]
    check("stale artifact filenames", not stale_names, ", ".join(stale_names) if stale_names else "none")

    stale_seed_files = []
    text_suffixes = {".py", ".md", ".csv", ".json", ".yml", ".txt", ".cff"}
    for path in ROOT.rglob("*"):
        if path.is_file() and path.suffix.lower() in text_suffixes:
            try:
                obsolete_seed = "2026" + "0820"
                if obsolete_seed in path.read_text(encoding="utf-8"):
                    stale_seed_files.append(str(path.relative_to(ROOT)))
            except UnicodeDecodeError:
                pass
    check("obsolete random seed", not stale_seed_files, ", ".join(stale_seed_files) if stale_seed_files else "none")

    meta = pd.read_csv(ROOT / "results" / "acute_EAE_random_effects_meta.csv").set_index("module")
    check("acute EAE pooled immune effect", math.isclose(meta.loc["Endothelial_immune_activation", "pooled_g"], 2.596911, abs_tol=1e-6), f"{meta.loc['Endothelial_immune_activation', 'pooled_g']:.6f}")
    check("acute EAE pooled BBB effect", math.isclose(meta.loc["BBB_specialization", "pooled_g"], -1.790610, abs_tol=1e-6), f"{meta.loc['BBB_specialization', 'pooled_g']:.6f}")
    check("modified Knapp Hartung intervals cross zero", meta.loc["Endothelial_immune_activation", "mKH_ci_low"] < 0 < meta.loc["Endothelial_immune_activation", "mKH_ci_high"] and meta.loc["BBB_specialization", "mKH_ci_low"] < 0 < meta.loc["BBB_specialization", "mKH_ci_high"], "both principal intervals cross zero")

    cross = pd.read_csv(ROOT / "results" / "cross_species_concordance_summary.csv").set_index("gene_set")
    row = cross.loc["All_common_genes_sensitivity"]
    check("cross species all gene result", int(row.n_genes) == 11633 and math.isclose(row.spearman_rho, 0.1972107329, abs_tol=1e-9), f"n={int(row.n_genes)}, rho={row.spearman_rho:.6f}")
    row = cross.loc["All_prespecified_state_genes"]
    check("cross species focused result", int(row.n_genes) == 43 and math.isclose(row.spearman_rho, 0.243431, abs_tol=1e-6), f"n={int(row.n_genes)}, rho={row.spearman_rho:.6f}")

    strict = pd.read_csv(ROOT / "results" / "strict_orthology_sensitivity_summary.csv")
    row = strict[(strict.mapping.str.startswith("HCOP")) & (strict.scope == "All shared genes")].iloc[0]
    check("strict orthology all gene result", int(row.n_genes) == 7705 and math.isclose(row.spearman_rho, 0.1868056122, abs_tol=1e-9), f"n={int(row.n_genes)}, rho={row.spearman_rho:.6f}")

    merfish = pd.read_csv(ROOT / "results" / "GSE284005_paired_summary.csv")
    row = merfish[merfish.feature == "Stress/inflammatory endothelial fraction"].iloc[0]
    check("MERFISH paired composition", int(row.donors_higher_in_Vas_Imm) == 5 and math.isclose(row.exact_two_sided_sign_flip_p, 0.0625, abs_tol=1e-12), f"{int(row.donors_higher_in_Vas_Imm)} of 5, P={row.exact_two_sided_sign_flip_p}")

    functional = pd.read_csv(ROOT / "results" / "published_source_barrier_and_vascular_summary.csv").set_index("endpoint")
    check("IgG leakage result", math.isclose(functional.loc["IgG leakage", "p_exact_two_sided"], 0.373385, abs_tol=1e-6), f"P={functional.loc['IgG leakage', 'p_exact_two_sided']:.6f}")
    check("fibrinogen leakage result", math.isclose(functional.loc["Fibrinogen leakage", "p_exact_two_sided"], 0.913945, abs_tol=1e-6), f"P={functional.loc['Fibrinogen leakage', 'p_exact_two_sided']:.6f}")

    source_index = pd.read_csv(ROOT / "figure_source_data" / "figure_source_index.csv")
    bad_source_rows = []
    for row in source_index.itertuples():
        source = ROOT / row.source_file
        if not source.is_file() or sha256(source) != row.sha256:
            bad_source_rows.append(f"{row.figure}: {row.source_file}")
    check("figure source index", not bad_source_rows, ", ".join(bad_source_rows) if bad_source_rows else f"{len(source_index)} mappings verified")

    if args.manuscript:
        manuscript_text, media_hashes = docx_text_and_media_hashes(args.manuscript)
        check("manuscript title", TITLE in manuscript_text, TITLE)
        missing_embedded = [path.name for path in main_figures if sha256(path) not in media_hashes]
        check("main figures embedded in manuscript", not missing_embedded, ", ".join(missing_embedded) if missing_embedded else "all five PNG hashes found in the DOCX media")

    passed = all(item["passed"] for item in checks)
    report = {"package": str(ROOT), "passed": passed, "checks": checks}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
