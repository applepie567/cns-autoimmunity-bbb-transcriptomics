# Manuscript alignment

Final title: **Brain endothelial responses in acute EAE partially overlap with vascular changes in multiple sclerosis**

Final manuscript file: `BBI_修订稿_投稿冻结版.docx`

Final manuscript SHA256: `14624b7c54e43292970310b8a763e872581d9d9edf7f3ca0b20c704a1da8083e`

The manuscript itself is not included in the public reproducibility archive. Its checksum fixes the precise version against which this package was audited.

## Fixed scope

* Main figures: Figure 1 to Figure 5
* Supplementary figures: Figure S1 to Figure S9
* Supplementary tables: Table S1 to Table S18
* Project seed: 20260901
* S15 focused-set bootstrap seed: 20260904 (2000 resamples); recovered original setting
* Figure 4: selected light gray background
* Infection context: source detail tables S11 and S12 and summary table S13; excluded from the main vascular conclusions
* GWAS analysis: not part of the final results or figure set
* Functional statement: transcript changes and reduced vascular proliferation do not establish BBB permeability recovery

## Final figure checksums

* `figures/main/Figure_1_study_design.png`  `adc69e327b740d08ed4c48e1d1625c15fb4ccae51f7d4b944a0b6d5d6bcb77fb`
* `figures/main/Figure_2_acute_EAE_meta.png`  `3fefdc5d9e1071e7c9690c4a324e45bed84cfc77beac8d68f0695c623f6e03e3`
* `figures/main/Figure_3_human_MS_microenvironment_stage.png`  `1ec24e102a38dc013654db7d25621e73a091e7444b945a954f88cb5d771febea`
* `figures/main/Figure_4_cross_species_concordance.png`  `38017578abc641d5724eb11628e25ff48dbeb6332117b50f34e369293e063db8`
* `figures/main/Figure_5_VEGFA_perturbation.png`  `fd4dd2e1ae0797fa4d15059fefdfda5fd35a5a0f6120ffd47517889851f6ac22`
* `figures/supplementary/Figure_S1_GSE210776_QC.png`  `d88573b37e69e77bea7a48a3db6f7a4d567b617ece631a3595a4bc23534925d6`
* `figures/supplementary/Figure_S2_acute_EAE_sample_scores.png`  `e1bc7d0a4bc92f1d3e9a3f777ac962c40d966937f7a6f1863bbbd624e139c0e5`
* `figures/supplementary/Figure_S3_GSE199460_subtype_sensitivity.png`  `b09dbb8e6a287ab8da42525efe72666c01303831ea3a13f1162734b3c488466f`
* `figures/supplementary/Figure_S4_GSE95401_time_course.png`  `18d98484f88705528e7ca7614d98dc16c783a1ac8149543304a9d4e50acc9476`
* `figures/supplementary/Figure_S5_GSE279183_EC_QC.png`  `ea410fae354c7e388f1ed9739cb3d91ed68d6cce6ef0a90ad59076d065ac2c90`
* `figures/supplementary/Figure_S6_GSE208747_threshold_sensitivity.png`  `39df662111d7fc95dc0ecaf17b0d5eb6caca84502b1627254da92bb196315e6b`
* `figures/supplementary/Figure_S7_VEGFA_all_modules.png`  `b0fe6f131d992fb05c18316c07815f274b63797986d02b0c837962084c90a3fd`
* `figures/supplementary/Figure_S8_published_functional_protein_source_data.png`  `73ac3554646d4fd04839afd861d4092f7fbd09fa7de0268cdeeaae26e325a200`
* `figures/supplementary/Figure_S9_GSE284005_spatial_validation.png`  `b4939314e1f8b9f8e25e71c239f106ce603dfe060a58dd1e183c3c110d88bac4`

## Archive preparation

The version DOI 10.5281/zenodo.22340814 was reserved on 5 September 2026. This update changes packaging, source-table fallback paths and archive verification. The frozen numerical tables, combined workbook and figure files remain byte-identical to the manuscript-aligned input package. The checksum above identifies the supplied manuscript before its Data Availability paragraph is updated with the published archive DOI.

The supplied manuscript describes the project seed but omits the S15 bootstrap seed exception. Before submission, its statistical reproducibility sentence should specify that seed 20260904 was used for the focused-set bootstrap intervals in the orthology sensitivity analysis. The frozen confidence limits and Figure 4 are preserved.
