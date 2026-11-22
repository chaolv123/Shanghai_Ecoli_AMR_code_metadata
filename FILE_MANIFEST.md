# File Manifest

## Public repository scope

- Manuscript files and publication figure artwork: not included.
- Full supplementary-table collection: not included.
- Synchronized Supplementary Tables 4, 9 and 11: 3 workbooks.
- Synchronized source-data workbooks: Source Data Fig. 1-6 (6 workbooks).

The synchronized publication files, checksum manifest and validation report
are stored in `outputs/publication_sync/`.

## Executable workflow additions (20 November 2026)

- `run_allRevised.py`: one-command entry point.
- `run_pcg_sampling_intensityRevised.py`: period-normalized PCG analysis.
- `run_derived_analysesRevised.py`: resistome/plasmid derived-input analysis.
- `update_fig2c_supp_table4Revised.R`: MDR-only Fig. 2c and Supplementary Table 4 rebuild.
- `run_fig2_profile_updateRevised.py`: cross-platform launcher for the R rebuild, including non-ASCII Windows paths.
- `run_publication_syncRevised.py`: Tables 4/9/11 and Source Data synchronization.
- `environmentRevised.yml` and `install_R_dependencies.R`: environment setup.
- `data/derived/`: supplied derived inputs needed in place of unavailable raw
  LEfSe/plasmid-prediction intermediates.
- `outputs/publication_sync/`: verified copies of Supplementary Tables 4/9/11 and
  Source Data Fig. 1-6 with SHA-256 manifest.
