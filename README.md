# Shanghai food-animal *E. coli* AMR reproducibility package

This repository reproduces the manuscript analyses from public **processed
inputs**. Raw sequencing reads are intentionally outside the workflow boundary.
The starting inputs are isolate metadata, acquired-ARG presence/absence calls,
core-SNP distances and PCG assignments, phenotype-linked WGS calibration
weights, and predicted plasmid-type/ARG assignments.

## One-command workflow

### Option 1: Conda

```bash
conda env create -f environmentRevised.yml
conda activate shanghai-ecoli-amr
python run_allRevised.py
```

### Option 2: existing Python and R installations

```bash
pip install -r requirementsRevised.txt
Rscript install_R_dependencies.R
python run_allRevised.py
```

The complete workflow normally finishes in several minutes on a standard
desktop computer. A successful run writes:

- `outputs/workflow/run_allRevised.log`
- `outputs/workflow/workflow_status.json` with `"status": "success"`
- recomputed metadata-derived source tables;
- WGS-representativeness, breakpoint, ST-ARG and farm-adjusted analyses;
- the 494-PCG threshold reconciliation;
- the second-round PCG sampling-intensity analysis;
- PCoA/PERMANOVA, supplied LEfSe-score validation and predicted-PT summaries;
- synchronized copies and SHA-256 hashes for Supplementary Tables 4/9/11 and Source
  Data Fig. 1-6;
- unit-test results.

For a faster smoke test of the PCoA/PERMANOVA stage, use:

```bash
python run_allRevised.py --pcoa-permutations 99
```

The default reproducibility run uses 999 permutations.

## Workflow modules

| Script | Purpose |
|---|---|
| `run_pipelineRevised.py` | Validates 5,299 phenotype records, 2,622 genomes, ARG counts, AMR rates, MDR and antibiogram length. |
| `run_pcg_reconciledRevised.py` | Applies R `hclust(method="ward.D2")` at 5/10/20/25, audits pairwise distances and verifies exact recovery of the primary 494-PCG partition. |
| `run_farm_clusteredRevised.py` | Reproduces farm-clustered robust GLMs used in Supplementary Table 2. |
| `run_extendedRevised.py` | Runs WGS selection/weighting, breakpoint and ST-ARG analyses; its optional `pcg` target delegates to the authoritative R Ward.D2 workflow. |
| `run_pcg_sampling_intensityRevised.py` | Adds period denominators, normalized PCG measures, calibration weighting, farm bootstrap and sampling-density/farm-adjusted models. |
| `run_derived_analysesRevised.py` | Recomputes PCoA/PERMANOVA and predicted-PT summaries; validates supplied LEfSe scores. |
| `update_fig2c_supp_table4Revised.R` | Rebuilds the MDR-only Fig. 2c accumulation curves, Supplementary Table 4 and synchronized Fig. 2 source data while preserving fixed profile codes. |
| `run_publication_syncRevised.py` | Validates and exports canonical Supplementary Tables 4/9/11 and Source Data Fig. 1-6 without overwriting publication-order files. |
| `run_allRevised.py` | Runs all modules and tests in dependency order. |

## Why raw reads are not required

Read assembly, ARG calling, MLST assignment and plasmid prediction are upstream
bioinformatics processes. This repository reproduces the statistical analyses
reported in the manuscript from their processed outputs. The exact LEfSe scores
and short-read predicted plasmid assignments are included as derived inputs,
because the original external LEfSe environment and raw plasmid-prediction files
were not available. The workflow validates these inputs and regenerates all
downstream summaries used here; it does not claim to reproduce read assembly or
plasmid prediction from FASTQ files.

## Input layout

- `data/raw/phenotype_metadata.xlsx`: 5,299 phenotype records.
- `data/raw/genomic_metadata.xlsx`: 2,622 genomes and acquired-ARG matrix.
- `data/raw/core_snp_distances.csv`: core-SNP distance matrix.
- `data/raw/pcg_original_assignments.xlsx`: historical PCG labels retained for audit only.
- `data/raw/pcg_metadata.xlsx`: metadata aligned to the SNP matrix.
- Ward.D2 isolate-level PCG assignments used by the sampling-intensity models
  are generated under `outputs/extended/pcg_reconciled/`.
- `data/derived/fig3_resistome_derived_input.xlsx`: exact LEfSe/figure-derived resistome inputs.
- `data/derived/fig6_plasmid_derived_input.xlsx`: predicted PT-ARG, PT burden and PT-ST inputs.
- `data/publication_inputs/`: canonical Supplementary Tables 4/9/11 and Source Data
  Fig. 1-6 used by the synchronization and validation module.

## Key output locations

- `outputs/extended/pcg_sampling_intensity/`: new second-round PCG analysis.
- `outputs/derived_reproducibility/`: PCoA/PERMANOVA, LEfSe validation and PT summaries.
- `outputs/fig2_profile_update/`: reproducible Fig. 2c, Supplementary Table 4 and synchronized Fig. 2 source-data outputs.
- `outputs/publication_sync/`: synchronized Tables 4/9/11 and Source Data Fig. 1-6 plus checksums.
- `outputs/extended/`: all other reviewer-response analyses.

## Archived analysis modules

The executable entry points are the repository-relative `*Revised.py` scripts
listed above. A small set of earlier Python analysis modules is retained under
`archive/legacy_python/` because the revised wrappers import their statistical
functions after overriding all input and output paths. They are not intended to
be run directly. Old path-dependent R working scripts are not included. R is
required only for the explicitly generated `glmmTMB` models in the PCG
sampling-intensity module.

## Synchronized publication data

To keep the public repository focused on code, processed inputs and
reproducibility outputs, manuscript files and publication figure artwork are
not included. Canonical, workflow-validated copies of Supplementary Tables 4/9/11
and Source Data Fig. 1-6 are provided under `outputs/publication_sync/`,
together with SHA-256 checksums and a validation report.

## Licensing

Code is released under the MIT License. Processed metadata, derived tables and
source data are released under CC BY 4.0 unless otherwise stated; see
`DATA_LICENSE.md`.
