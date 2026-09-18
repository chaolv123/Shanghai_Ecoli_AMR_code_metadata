# Code

Run these scripts from the repository root or directly with Python paths.

```bash
pip install -r requirements.txt
python code/run_pipeline.py --check
python code/run_pcg_reconciled.py
python code/run_farm_clustered.py
python code/run_extended.py --analyses wgs weighted_arg breakpoint st_arg
python code/export_antibiogram_mdr_discrepancies.py
```

To run the full public workflow:

```bash
python code/run_all.py
```

Inputs are read from `06_Metadata_and_Input_Data/`. Generated files are written
to `outputs/`, which is intentionally excluded from version control.
