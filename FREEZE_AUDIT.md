# Freeze and archive preparation audit

Status: **PASS for package and frozen-output validation**.

The supplied manuscript and its five main images match the frozen archive. All scientific result files, figure-source tables, the combined workbook and the deposited figures remain byte-identical. Fourteen PNG figures were regenerated in a separate copy and matched byte for byte. The original deposited PDF is retained.

Version DOI **10.5281/zenodo.22340814** was reserved in the existing Zenodo new-version draft. Publication has not been performed during package preparation.

The restored strict-orthology script reproduces all four S15 summary rows from archived gene effects, including the original 2000-resample confidence intervals with seed **20260904**. The general project seed remains **20260901**. The seed exception must be added to the manuscript's Section 2.10 when Data Availability is updated for submission.

The original raw GEO matrices and dated HCOP export were not available for re-execution in this packaging session. See REPRODUCIBILITY.md for the verified scope and remaining external inputs.

The obsolete analysis/ and analysis_results/ directories are excluded from the current release tree. Input instructions are restored, source-table fallback is repaired and the verifier can now detect missing files, unexpected archive files and obsolete directories.

## Checks

* PASS — required files: all required files present
* PASS — legacy directories removed: none
* PASS — version DOI metadata: 10.5281/zenodo.22340814 (reserved during package preparation)
* PASS — main figure count: found 5
* PASS — supplementary figure count: found 9
* PASS — supplementary tables S1 to S18: represented tables: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
* PASS — combined workbook sheets: found 20 sheets
* PASS — workbook formula errors: none
* PASS — random seed: 20260901
* PASS — orthology bootstrap seed: 20260904; 2000 resamples, as used by the original S15 analysis
* PASS — selected Figure 4 background: light gray settings present
* PASS — stale artifact filenames: none
* PASS — obsolete random seed: none
* PASS — acute EAE pooled immune effect: 2.596911
* PASS — acute EAE pooled BBB effect: -1.790610
* PASS — modified Knapp Hartung intervals cross zero: both principal intervals cross zero
* PASS — cross species all gene result: n=11633, rho=0.197211
* PASS — cross species focused result: n=43, rho=0.243431
* PASS — strict orthology all gene result: n=7705, rho=0.186806
* PASS — concordance recomputed from frozen gene effects: all four counts, correlations and direction fractions match
* PASS — MERFISH paired composition: 5 of 5, P=0.0625
* PASS — IgG leakage result: P=0.373385
* PASS — fibrinogen leakage result: P=0.913945
* PASS — figure source index: 39 mappings verified
* PASS — manuscript title: Brain endothelial responses in acute EAE partially overlap with vascular changes in multiple sclerosis
* PASS — main figures embedded in manuscript: all five PNG hashes found in the DOCX media
* PASS — scientific files preserved: 96 result, figure-source, figure and supplementary files are byte-identical to the input package
* PASS — PNG figure reproduction: 14 of 14 regenerated PNG files are byte-identical
* PASS — S15 full numerical reproduction: All four rows, including confidence intervals and P values, match using bootstrap seed 20260904
* PASS — source-table fallback: 19 bundled CSV tables copied with identical hashes and reused on a second call

## Archive integrity

MANIFEST.sha256 is generated after this audit and contains every public file except itself. Run `python code/verify_submission_freeze.py --check-manifest` on the original downloaded archive, before regenerating outputs.
