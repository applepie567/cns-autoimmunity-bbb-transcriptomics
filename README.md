# Endothelial responses in EAE and multiple sclerosis

## Current analysis

**Endothelial responses in acute EAE show limited overlap with those in chronic multiple sclerosis**

Analysis release **v3.0.0**, prepared on 2026-09-07.
Reserved version DOI 10.5281/zenodo.22641087. Registration is completed when the corresponding new Zenodo version is published.

The current analysis is under `analysis_strengthening/fbcns_revision`.
It compares all ten pairs of five cohorts on the same 6,611 genes. It includes
animal and donor resampling, an additional human endothelial cohort, capillary
comparisons in the same donors, and focused barrier functions. The original
acute EAE and spatial results are retained as supporting analyses.

Mean correlation was 0.389 for mouse cohort pairs and 0.060 for all six mouse
and human pairs. These values refer to exactly the same genes and comparisons.
Eight of eleven focused gene sets met the strict coverage threshold. None of
their 40 tests was significant after correction for the false discovery rate.
The resampling ranges describe uncertainty within the selected small cohorts.

## Find the current files

| Material | Location |
| --- | --- |
| Five current main figures | `analysis_strengthening/fbcns_revision/figures/Figure_1.png` through `Figure_5.png`, with PDF copies |
| Supplementary Figures S1 to S9 | `analysis_strengthening/fbcns_revision/figures/supplementary` |
| Supplementary Figure S10 | `analysis_strengthening/fbcns_revision/figures/Figure_S10.png` and PDF |
| Supplementary methods and legends | `analysis_strengthening/fbcns_revision/documentation/Additional_file_1_Supplementary_methods_and_figures.docx` |
| All tables S1 to S22 | `SUPPLEMENTARY_TABLE_INDEX.csv` maps each table to its file |
| Current figure source data | `FIGURE_SOURCE_INDEX.csv` maps panels to included files |
| Animal and donor inputs | `analysis_strengthening/derived` |
| Original input filenames and URLs | `analysis_strengthening/output/INPUT_FILES.csv` |
| Five cohort and focused function code | `analysis_strengthening/fbcns_revision/code` |
| Human validation and reconstruction code | `analysis_strengthening/code` |
| Reproduction instructions | `REPRODUCIBILITY.md` |

## Reproduce the revised analysis

Use Python 3.12 in a separate environment. From the repository root run

```bash
python -m pip install -r requirements.txt
python code/verify_release.py
python analysis_strengthening/fbcns_revision/code/analyze_fbcns.py
python analysis_strengthening/fbcns_revision/code/make_figures.py
```

Check the manifest before regenerating files. New plot exports and numerical
formatting may change binary checksums after regeneration. See
`REPRODUCIBILITY.md` for upstream reconstructions and all raw download steps.
The included derived expression profiles are sufficient to rerun the revised
five cohort comparison. They are not a complete collection of raw sequencing
matrices. Source annotations and the retained strict mapping define the
eligibility and matching rules for this version.

## Previous files and version history

The root folders `code`, `results`, `figure_source_data`, `figures` and
`supplementary_tables` retain v2 source material so previous links continue
to work. Their README files point to the current outputs. Use the current
figure paths above for the revised manuscript. `code/verify_release.py` is
the v3 manifest verifier. Other root code belongs to the earlier analysis.
The separate baseline copy under `analysis_strengthening/baseline/BBI_v2.0.0`
provides the exact relative input paths used by the revised scripts.
Historical infection related tables are retained for provenance and are not
part of the current EAE and MS endothelial comparison.

Previous v2.0.0 archive: https://doi.org/10.5281/zenodo.22340814

Version series: https://doi.org/10.5281/zenodo.22031404

Code repository: https://github.com/applepie567/cns-autoimmunity-bbb-transcriptomics

The software archive retains Yuan Feng as creator, consistent with the existing
release. The manuscript has its own author list. This archive supplies the
analysis materials and the manuscript checksum. The full manuscript and cover
letter are supplied separately for journal submission.

## Licenses

Code is MIT licensed. Newly generated derived tables and metadata are licensed
under CC BY 4.0. Source study materials and third party annotations retain their
original terms. See `LICENSE` and `LICENSE-DATA.md`.
