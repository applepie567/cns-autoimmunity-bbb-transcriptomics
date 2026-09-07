# Inputs downloaded separately

Downloaded GEO and source-study files belong here. Original data are not distributed in this archive. Download each source from its GEO accession or the original publication listed in DATA_SOURCES.md. Keep original filenames, including sample identifiers. The version-pinned requirements include h5py for raw Visium HDF5 input.

## Paths used by the main reconstruction scripts

All paths below are relative to public_data/.

| Input | Expected location or naming pattern | Script |
| --- | --- | --- |
| GSE210776 | GSE210776/*.tar.gz; filenames include _DS001_processed through the deposited DS samples, or _BEVAC001_processed through BEVAC003 | code/run_bbi_extension.py |
| GSE199460 annotations | GSE199460/GSE199460_cell_annotation.meta_data.cd31_selection.csv.gz | code/run_bbi_extension.py |
| GSE199460 normalized expression | GSE199460/GSE199460_normalized_expr.cd31_selection.sctransform.csv.gz | code/run_bbi_extension.py |
| GSE95401 source workbook | GSE95401/GSE95401_MunjiSoungDaneman_RNAseq_BBB_Health_Disease.xlsx; worksheet Replicates and Averages | code/run_bbi_extension.py |
| GSE279183 single-nucleus files | GSE279183_snRNA/GSM*_matrix.mtx.gz, *_features.tsv.gz and *_barcodes.tsv.gz; sample prefixes are listed in analyze_gse279183() | code/run_bbi_extension.py |
| GSE279183 deposited source tables | GSE279183_Supplementary_Tables_1_10.xlsx; worksheet STable_6_Cell_type_DEG | code/run_bbi_extension.py |
| GSE208747 Visium archives | GSE208747/*.tar.gz with original _A1.tar.gz, _C1.tar.gz, _M1.tar.gz and _N1.tar.gz sample suffixes, etc. | code/run_bbi_extension.py |
| GSE284005 MERFISH | GSE284005/GSM*_celltypes.tsv.gz and corresponding *_count.tsv.gz files | code/run_gse284005_spatial_validation.py |
| Published r84, TEER and human protein source files | published_validation/Figure 7_b,c,e,f,g.pzfx; published_validation/ED_Figure7_c-h.pzfx; published_validation/ED_Figure6_a.pzfx; published_validation/Figure 5_d-f.pzfx | code/analyze_orthogonal_validation.py |

## Earlier infection-context reconstruction helper

code/run_reanalysis.py is a supplementary reconstruction and audit helper. It expects GSE163005_*.gz files, GSE163194_RAW.tar, the spatial GSE279183 files listed inside that script, and GSM636*.tar.gz Visium files directly in public_data/. Its outputs are written to reconstruction_results/, not to the frozen results/ directory.

The helper can initialize reconstruction_results/source_tables from the bundled supplementary_tables/*.csv files. If supplied, the older workbook MS_BBB_complete_source_data_tables.xlsx in this directory remains an alternative input.

## Strict orthology

The manuscript reports an HCOP download date of 4 September 2026. The final selected gene effects and sensitivity summary are archived in results/strict_orthology_gene_effects.csv and results/strict_orthology_sensitivity_summary.csv. The original analytical functions were recovered from the 4 September 2026 public-data package and are supplied as code/run_strict_orthology.py. By default, that script uses the archived selected gene effects and reproduces all four S15 rows, including the 2000-resample bootstrap intervals with the original seed 20260904. The original dated HCOP export is not included. With that original file available, pass --hcop public_data/HCOP/human_mouse_hcop_fifteen_column.txt.gz to repeat the recovered filter. The filter retains reciprocal unique pairs supported by Ensembl and NCBI, restricted to equal uppercase symbols within the primary comparison. A current HCOP export may select a different set, which the script detects.

See REPRODUCIBILITY.md before running reconstruction in a separate working copy. Do not overwrite the frozen release with newly downloaded or regenerated analysis results.
