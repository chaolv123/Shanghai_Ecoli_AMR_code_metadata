# Shanghai food-animal E. coli AMR analysis code and metadata

This repository contains the public reproducibility package for the revised
food-animal *Escherichia coli* AMR surveillance analysis. It includes cleaned
metadata, reviewer-response analysis code and lightweight tests.

Figures, supplementary figures, supplementary tables, manuscript files and
pre-generated output workbooks are intentionally not included in this public
GitHub package.

## Repository Structure

- `06_Metadata_and_Input_Data/`: public input metadata and SNP-distance data.
- `code/`: Python scripts for regenerating the reviewer-response analyses.
- `code/archive/legacy_python/`: original reviewer-analysis modules called by
  the wrapper scripts.
- `data/supporting/`: small supporting files for ST69 virulence/tree analyses.
- `tests/`: basic validation tests for regenerated outputs.
- `LICENSE`: MIT License for code.
- `DATA_LICENSE.md`: CC BY 4.0 license for processed metadata and derived data.

## Input Data

The canonical public input files are:

- `phenotype_metadata.xlsx`: 5,299 phenotyped isolates with corrected
  `antibiogram_length` and `MDR_binary` values.
- `genomic_metadata.xlsx`: 2,622 WGS isolates with ARG profiles.
- `pcg_metadata.xlsx`: metadata for isolates used in PCG/SNP analyses.
- `pcg_original_assignments.xlsx`: audited author PCG assignments.
- `core_snp_distances.csv.zip`: zipped core-genome SNP distance matrix, kept
  below GitHub's browser-upload single-file limit.
- `metadata_farm_analysis_original.xlsx`: farm-adjusted model input metadata.
- `ARG_metadata.xlsx`: ARG metadata supplied for sensitivity analyses.

Farm and district fields are intentionally retained as coded metadata to support
reproducible farm- and district-adjusted analyses.

## Reproducing Analyses

Install Python dependencies:

```bash
pip install -r requirements.txt
```

Run the full public workflow:

```bash
python code/run_all.py
```

Or run individual components:

```bash
python code/run_pipeline.py --check
python code/run_pcg_reconciled.py
python code/run_farm_clustered.py
python code/run_extended.py --analyses wgs weighted_arg breakpoint st_arg
python code/export_antibiogram_mdr_discrepancies.py
```

Generated files are written to `outputs/`. This folder is ignored by Git and can
be regenerated from the public inputs.

After running the scripts, validate the regenerated outputs:

```bash
python -m unittest discover -s tests -v
```

## Notes

Raw sequencing reads are not stored in this repository and should be linked via
the final sequence-read accession or repository record when available.
