# Public input data and dataset roles

The original matrices and archives are hosted by NCBI GEO and are not copied
into this repository. The repository contains newly computed outputs and frozen
source-to-figure tables needed to audit the manuscript.

| GEO accession | Organism / modality | Role in this study |
|---|---|---|
| [GSE210776](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210776) | Mouse endothelial single-cell transcriptomics | Acute EAE reference state and VEGF-A perturbation |
| [GSE199460](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199460) | Mouse single-cell transcriptomics | Independent endothelial support and exploratory sender analysis |
| [GSE95401](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95401) | Mouse CNS endothelial bulk RNA-seq | Temporal directionality support |
| [GSE279183](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279183) | Human MS single-nucleus and spatial transcriptomics | Human vascular-niche validation and partial raw spatial audit |
| [GSE208747](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208747) | Human MS Visium spatial transcriptomics | Independent lesion-stage and threshold-sensitivity analysis |
| [GSE163005](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163005) | Human CSF single-cell transcriptomics | Viral-encephalitis specificity comparator |
| [GSE163194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163194) | Human CSF single-cell transcriptomics | Bacterial-meningitis specificity comparator |

## Input policy

- Keep public GEO inputs in the ignored local directory `upload/`.
- Do not commit human-level raw matrices, archives, histology images, or GEO
  supplementary files to this repository.
- The complete filename manifest and dataset-specific sample mapping are recorded
  in the manuscript, supplementary methods, and `analysis_results/`.
- Cells and spatial spots are not treated as independent biological replicates.

## Derived outputs

The reconstructed QC, sample-level summaries, module scores, statistical
contrasts, and direction-audit tables were generated for this study. Their
licensing and the boundary for third-party source-provided values are described
in `LICENSE-DATA.md`.
