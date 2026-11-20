# File Manifest

## Public repository scope

- Manuscript files and publication figure artwork: not included.
- Full supplementary-table collection: not included.
- Synchronized Supplementary Table 9: 1 workbook.
- Synchronized source-data workbooks: Source Data Fig. 1-6 (6 workbooks).

The synchronized publication files, checksum manifest and validation report
are stored in `outputs/publication_sync/`.

## Executable workflow additions (20 November 2026)

- `run_allRevised.py`: one-command entry point.
- `run_pcg_sampling_intensityRevised.py`: period-normalized PCG analysis.
- `run_derived_analysesRevised.py`: resistome/plasmid derived-input analysis.
- `run_publication_syncRevised.py`: Table 9 and Source Data synchronization.
- `environmentRevised.yml` and `install_R_dependencies.R`: environment setup.
- `data/derived/`: supplied derived inputs needed in place of unavailable raw
  LEfSe/plasmid-prediction intermediates.
- `outputs/publication_sync/`: verified copies of Supplementary Table 9 and
  Source Data Fig. 1-6 with SHA-256 manifest.
