# Cross-species BBB remodeling in CNS autoimmunity

This repository contains the analysis code, prespecified signatures, derived
result tables, software environment, and reproducibility metadata for:

> Cross-species single-cell and spatial transcriptomic analysis reveals
> context-dependent blood–brain barrier remodeling in central nervous system
> autoimmunity

The study is an original, question-driven computational analysis of accession-
level count matrices, processed expression objects, archived files, and sample
metadata. It is not a literature review or a pooling of published effect
estimates.

## Public datasets

The analysis uses public data from GEO accessions
[GSE210776](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210776),
[GSE199460](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199460),
[GSE95401](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95401),
[GSE279183](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279183),
[GSE208747](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208747),
[GSE163005](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163005), and
[GSE163194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163194).

Original GEO files are not duplicated in this repository. Download the relevant
GEO supplementary files into a local `upload/` directory before a complete raw
rerun. The scripts accept the original GEO filenames and tolerate browser-added
suffixes such as `(1)`.

## Repository contents

- `analysis/run_reanalysis.py`: deterministic raw-matrix reconstruction, sample
  aggregation, exact permutation testing, Mann–Whitney testing, and QC exports.
- `analysis/make_submission_figures.py`: generation of the main and supplementary
  figures from the included source and reconstructed result tables.
- `analysis_results/`: aggregate numerical outputs and the source-to-figure
  tables used for the manuscript.
- `requirements.txt` and `environment.yml`: pinned Python environment.
- `CITATION.cff`: author and software citation metadata.
- `LICENSE` and `LICENSE-DATA.md`: code and data licensing terms.

## Reproducibility settings

- Fixed random seed: `20260820`
- Primary inferential unit: animal, donor, or patient
- Cells and Visium spots: used only for within-sample aggregation
- Exact label permutations: fully enumerated and deterministic
- Multiple testing: Benjamini–Hochberg FDR within prespecified families
- Signed-effect palette: blue–white–red, centered at zero
- No stochastic embedding is used in the submitted figure set

The fixed values are also recorded in
`analysis_results/analysis_parameters.json`.

## Public-data boundary

This public repository intentionally excludes patient-, donor-, animal-, and
slide-level derived score tables. Those tables can be reconstructed locally
from the cited public GEO inputs with `analysis/run_reanalysis.py`. The deposited
tables are limited to aggregate contrasts, prespecified gene-set coverage,
software metadata, and source-to-figure summaries.

## Installation

Using conda:

```bash
conda env create -f environment.yml
conda activate bbb-cns-autoimmunity
```

Or using Python 3.12 and pip:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Reproduce the deposited figures

The deposited result tables are sufficient to regenerate the figure set:

```bash
python analysis/make_submission_figures.py
```

Figures are written under `outputs/bmc_submission/figures/`.

## Complete raw-matrix rerun

1. Create `upload/` in the repository root.
2. Download the accession-associated matrices and archives from GEO.
3. Keep the GEO/GSM prefixes in the filenames unchanged.
4. Run:

```bash
python analysis/run_reanalysis.py
python analysis/make_submission_figures.py
```

Because some public cohorts distribute processed objects or atlas-level outputs
rather than a single uniform raw format, the role of each dataset and the
raw-versus-frozen audit are documented in the manuscript and in
`analysis_results/raw_vs_frozen_direction_audit.csv`.

## Interpretation boundary

The GSE163005 and GSE163194 reconstructions are independent canonical-marker
audits and are not claimed to recreate every historical atlas preprocessing or
cell-labeling step. GSE279183 raw spatial reconstruction covers the supplied
slides and is used for QC and directional robustness; full-cohort inference
remains anchored to donor-level processed source tables. Transcriptomic changes
are not direct measurements of blood–brain barrier permeability.

## Citation

Please cite the Zenodo version 1.0.0 archive and the associated manuscript:

- DOI: [10.5281/zenodo.22031405](https://doi.org/10.5281/zenodo.22031405)
- GitHub repository:
  [applepie567/cns-autoimmunity-bbb-transcriptomics](https://github.com/applepie567/cns-autoimmunity-bbb-transcriptomics)

This DOI identifies the version 1.0.0 Zenodo record and becomes resolvable
when that record is published. Citation metadata are provided in `CITATION.cff`.

## Contact

Yuan Feng ([ORCID 0000-0002-1839-1316](https://orcid.org/0000-0002-1839-1316))

## Licensing

The analysis code is released under the MIT License. Newly generated derived
tables and metadata are released under CC BY 4.0. Public source datasets and
source-provided values remain subject to their original repository and study
terms; see `LICENSE-DATA.md`.
