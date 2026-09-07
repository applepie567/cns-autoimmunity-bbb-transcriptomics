# Reproduce the v3 analysis

Use Python 3.12 and the packages in requirements.txt. Work on a copy of the
archive when recalculating results. The frozen files remain the reference.

## Primary five cohort analysis

From the repository root run

```bash
python code/verify_release.py
python analysis_strengthening/fbcns_revision/code/analyze_fbcns.py
python analysis_strengthening/fbcns_revision/code/make_figures.py
```

The first analysis command reads the included biological sample NPZ profiles,
the eligible validation genes, and the original gene sets. It writes S21 and
S22 files under analysis_strengthening/fbcns_revision/results. It attempts
2,000 stratified resamples with seed 20260910. The focused function calculation
uses seed 20260911. The supplied results contain 1,841 estimable draws for
the primary comparison. Undefined standardized effects are removed jointly
across cohort pairs in each draw. Ranges are descriptive conditional percentiles.

The figure command regenerates Figures 1, 2 and 4, Supplementary Figures S1
to S10, and copies the supplied current Figures 3 and 5 to the final directory.
Figure export metadata may change PDF checksums without changing the plotted
values. Figures 3 and 5 can be regenerated separately as described below.

## Human validation and earlier four cohort reference

Run these commands in order if recomputing the inputs to the primary analysis

```bash
python analysis_strengthening/code/benchmark_concordance.py
python analysis_strengthening/code/analyze_macnair.py
python analysis_strengthening/code/joint_validation_uncertainty.py
python analysis_strengthening/fbcns_revision/code/analyze_fbcns.py
```

The earlier four cohort reference uses 7,705 genes. The primary five cohort
analysis uses 6,611 genes. The capillary sensitivity comparison uses 5,811
genes shared by its three donor and subtype subsets. These are distinct
comparisons. The S19 file labels ci_low and ci_high denote descriptive
percentiles, consistent with the current manuscript.

The human validation reads the included white matter endothelial count sums
and donor grouping records. These are derived from published source nuclei.
No original raw sequencing download is needed for this downstream step.

To regenerate Figure 3 and Figure 5 before collecting the final figures run

```bash
python analysis_strengthening/code/make_figure3_label_revision.py
python -c "import sys; sys.path.insert(0, 'analysis_strengthening/code'); from make_strengthening_figures import validation; validation()"
python analysis_strengthening/fbcns_revision/code/make_figures.py
```

Only the validation function is used from make_strengthening_figures.py.
Its other figures depict the earlier four cohort analysis rather than the
current primary comparison.

## Original input reconstruction

The input manifest lists original filenames, URLs, sizes and available MD5
values. The source download command accepts a dataset filter

```bash
python analysis_strengthening/code/download_inputs.py --raw-dir /absolute/path/source_data --dataset Macnair
```

The Macnair count matrix is approximately 4.56 GB compressed and is not
bundled. After downloading its row, column and matrix files, set BBI_RAW_DIR
to the same source directory and run extract_macnair.py followed by
analyze_macnair.py. The other reconstruction scripts use the same variable.
Source study documentation supplies original collection and processing details.

The archived strict orthologue table supplies the fixed mapping used here.
The full downloaded HCOP source snapshot is not bundled. The original v2
REPRODUCIBILITY.md and public_data/README.md under the baseline directory
describe the original reconstruction scope. No new assay data were collected.
