# Reproducibility report — 20 November 2026

## Status

`python run_allRevised.py` completed successfully from the repository root.
All six unit tests passed. See `outputs/workflow/workflow_status.json` and
`outputs/workflow/run_allRevised.log`.

## Verified inputs

- 5,299 unique phenotype records: 2,305 poultry and 2,994 swine.
- 2,622 unique genomic records.
- 2,547 phenotype-linked WGS records.
- Zero ARG-count/matrix discrepancies.
- Zero corrected antibiogram-length or MDR discrepancies.

## New PCG sampling-intensity analysis

- Primary threshold: 10 core-genome SNPs.
- R `hclust(method="ward.D2")` with `cutree(h=10)` exactly reproduced the
  historical primary partition; all 494 PCGs had maximum pairwise distances
  of 10 SNPs or fewer.
- 2,610 PCG-analyzable isolates.
- 494 PCGs: 308 localized and 186 non-localized.
- Localized-PCG isolate proportions:
  - poultry: 258/655 (39.39%) pre-2017 and 155/473 (32.77%) post-2017;
  - swine: 236/774 (30.49%) pre-2017 and 181/708 (25.56%) post-2017.
- Calibration-weighted localized-PCG proportions:
  - poultry: 34.26% pre-2017 and 31.76% post-2017;
  - swine: 28.99% pre-2017 and 25.84% post-2017.
- Farm-cluster bootstrap confidence intervals for weighted post-minus-pre
  differences included zero in both hosts.
- Joint calibration-weighted, sampling-density/fraction-adjusted models with
  farm-clustered standard errors were non-significant in both hosts.

Full output:
`outputs/extended/pcg_sampling_intensity/Supplementary_PCG_sampling_intensity_normalized_analysis.xlsx`.

## Resistome and predicted-plasmid reproducibility

- PCoA/PERMANOVA was recomputed from the 2,622-isolate acquired-ARG matrix:
  poultry R2 = 0.005638; swine R2 = 0.038540; both permutation P = 0.001
  with 999 permutations.
- Supplied LEfSe score tables contain 41 poultry and 44 swine period-associated
  ARG rows and are now included as explicit derived inputs.
- Predicted PT richness was reproduced as 67 in poultry, 108 in swine and 50
  shared between hosts.

## Publication-file synchronization

Verified copies of Supplementary Tables 9/11 and Source Data Fig. 1-6 are under
`outputs/publication_sync/`. The SHA-256 manifest proves that these copies are
identical to the canonical publication-order files in this repository.

## Workflow boundary

The executable workflow starts from processed phenotype/genomic metadata,
acquired-ARG calls, core-SNP distances/PCG labels and predicted plasmid-type
assignments. Raw FASTQ processing, genome assembly and de novo plasmid
prediction are not included and are not claimed to be reproduced here.
