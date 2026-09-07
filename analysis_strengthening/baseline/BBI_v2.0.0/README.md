# BBI submission data freeze

## Manuscript

**Brain endothelial responses in acute EAE partially overlap with vascular changes in multiple sclerosis**

Version 2.0.0, frozen on 5 September 2026.

This archive contains the analysis scripts, version-pinned software specifications, derived intermediate results, figure source data, final figures, Supplementary Tables S1 to S18 and an audit that links reported results to their source tables. The original GEO matrices are not redistributed.

## Main finding and interpretation boundary

Across three acute EAE cohorts, endothelial immune activation increased and BBB specialization decreased in the same direction. The modified Knapp Hartung intervals crossed zero, so the pooled effect magnitudes remain uncertain. Agreement with chronic active MS was weak at the gene level and the clearest human changes occurred in local inflammatory vascular regions. Published VEGF A blockade data showed effects on selected transcripts and vascular proliferation without a detectable reduction in IgG or fibrinogen leakage. The archive therefore supports transcriptomic associations and independent functional comparisons. It does not establish a new causal endothelial mechanism or functional BBB recovery.

## Archive layout

* `code` contains the analysis and figure scripts used for the final manuscript.
* `results` contains biological sample level scores, effect estimates, meta analysis results, cross species comparisons and published functional data summaries.
* `figure_source_data` contains the exact derived tables read by each final figure and an index that maps panels to files.
* `supplementary_tables` contains CSV files for Tables S1 to S18 and the combined workbook `Supplementary_Tables_S1-S18.xlsx`.
* `figures/main` contains the five final main figures. Figure 4 uses the selected light gray background.
* `figures/supplementary` contains Supplementary Figures S1 to S9.
* `figures/BBI_Supplementary_Figures_S1-S9.pdf` contains the merged supplementary figure set.
* `MANUSCRIPT_ALIGNMENT.md` records the final manuscript checksum and figure checksums.
* `FREEZE_AUDIT.md` records the manuscript, figure, table and numerical consistency checks.
* `public_data/README.md` lists input locations expected by the reconstruction scripts.
* `REPRODUCIBILITY.md` describes the verified workflow and the limits of reconstruction from original inputs.
* `MANIFEST.sha256` records every public file except the manifest itself.

## Reproduce figures from frozen tables

```bash
conda env create -f environment.yml
conda activate bbi-endothelial-freeze
python code/verify_submission_freeze.py --check-manifest
python code/make_submission_figures.py
python code/run_strict_orthology.py
python code/verify_submission_freeze.py
```

Check the manifest on the downloaded archive before regenerating outputs. Regeneration also creates TIFF exports and rebuilds the supplementary PDF, whose timestamp metadata can change its binary checksum. The figure workflow reads only files included in this archive. Run reconstruction from original inputs in a separate working copy, following `public_data/README.md` and `REPRODUCIBILITY.md`.

## Statistical units

Animals, donors or independent tissue samples are used as the units for statistical inference. Cells, nuclei, microscopic fields and spatial spots remain nested within their biological samples. Human image level immunohistochemistry and well level electrical resistance data are summarized descriptively when donor or experiment identifiers were unavailable.

## Repository and version DOI

* Code repository: https://github.com/applepie567/cns-autoimmunity-bbb-transcriptomics
* Version 2.0.0 DOI: [doi:10.5281/zenodo.22340814](https://doi.org/10.5281/zenodo.22340814)
* Previous version 1.0.0 DOI: [doi:10.5281/zenodo.22031405](https://doi.org/10.5281/zenodo.22031405)
* DOI for the version series: [doi:10.5281/zenodo.22031404](https://doi.org/10.5281/zenodo.22031404)

The version 2.0.0 DOI was reserved in an existing Zenodo new-version draft during package preparation on 5 September 2026. Reservation does not constitute publication. Zenodo registers the DOI when that draft is published. Cite the version DOI to identify these particular archived files.

The project seed is 20260901. The focused-set bootstrap intervals in Table S15 use the recovered original seed 20260904 with 2000 resamples, recorded in analysis_parameters.json.

The software archive retains Yuan Feng as its creator, consistent with the preceding software release. The manuscript has its own author list.

## Licenses

Code is distributed under the MIT License. Newly generated derived tables and metadata are distributed under CC BY 4.0. Original public data and source study files retain their original terms.
