# Reproducibility scope for version 2.0.0

## Verified from the deposited package

1. The supplied manuscript SHA256 matches MANUSCRIPT_ALIGNMENT.md, and the five main PNG figures are the exact images embedded in that manuscript.
2. All frozen figure-source tables, numerical result tables, the combined S1–S18 workbook and final PNG/PDF files are retained without numerical or binary changes in this packaging update.
3. The figure workflow was executed in a separate copy using Python 3.12.13 and the available package versions recorded in software_versions.csv. All five main and nine supplementary PNG figures matched the deposited PNG files byte for byte.
4. The verification script checks the specified results, figure-source file hashes, table coverage, citation DOI and absence of the two legacy directories. Its --check-manifest option checks all files and detects missing, changed or unlisted archive files.
5. The four cross-species comparison rows are independently checked against archived gene effects for gene count, Spearman correlation and agreement in direction. The restored code/run_strict_orthology.py additionally reproduces the confidence intervals and P values from the archived gene effects.

## Run the verified workflow

Use a separate copy of the downloaded archive for regenerated outputs.

```bash
conda env create -f environment.yml
conda activate bbi-endothelial-freeze
python code/verify_submission_freeze.py --check-manifest
python code/make_submission_figures.py
python code/run_strict_orthology.py
python code/verify_submission_freeze.py
```

Run --check-manifest before regeneration. The plotting scripts also create TIFF exports, and the rebuilt supplementary PDF can have different timestamp metadata. Those new outputs are not part of the original immutable checksum set. The archived supplementary PDF has nine pages. No changed figure or PDF from the verification copy was substituted into the release.

## Reconstruction from original inputs

The main reconstruction scripts require separately downloaded GEO and source-study files at the locations documented in public_data/README.md. The raw matrices were not available in this packaging session, so reconstruction from raw inputs was not repeated. h5py is pinned for that workflow but was not installed or executed in the figure-verification runtime.

The package preserves scripts for the main EAE/MS reconstruction, the targeted MERFISH analysis and the published functional comparisons. The original orthology-analysis functions were recovered and adapted into code/run_strict_orthology.py. Its frozen-input mode reproduces all four Table S15 rows, including confidence limits and P values. Its --hcop mode implements the original filter, but the original dated HCOP export is not included and that mode was not rerun against the original export. The package supports reproduction of the final figures and frozen S15 calculations, but should not be described as a tested single-command reproduction of every analysis from original downloads.

Infection-context source details remain in supplementary tables S11 and S12, with their summary in S13. They are outside the main EAE/MS vascular conclusions. Retaining these supplemental CSV files is separate from removing the obsolete analysis/ and analysis_results/ directories from the release tree.

## Citation and publication state

The new version DOI is 10.5281/zenodo.22340814. It was reserved during preparation and becomes registered when the existing Zenodo version draft is published. The DOI printed in this package does not itself establish that publication has occurred.

## Recovered S15 provenance and seed exception

The strict-orthology functions were recovered from public_data_strengthen_bbi.py in BBI_公共数据补强分析包.zip, created on 4 September 2026. The source script hash is 7ee97a14f83f8dde88a56a64c9a4b623ba88ba7f463a9a421cd91cdf12a96f32. The restored script changes local paths and separates analysis from manuscript editing and obsolete figure generation. It writes only to reconstruction_results/strict_orthology by default, so the published results and Figure S9 are preserved.

The recovered focused-set bootstrap seed is 20260904, with 2000 resamples. Both frozen S15 confidence intervals match that seed exactly. A seed of 20260901 produces different intervals. The general project seed remains 20260901, and the exception is explicitly recorded in analysis_parameters.json and software_versions.csv. The supplied manuscript should have this exception added to Section 2.10 when its Data Availability statement is updated for submission.
