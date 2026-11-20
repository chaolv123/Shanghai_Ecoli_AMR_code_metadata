# Public Release Checklist

Before uploading or making the GitHub repository public:

- Farm and district metadata are intentionally included as coded public
  metadata for reproducibility.
- Supplementary figures included in the supplied revision package are
  Supplementary Fig. 1-4.
- `LICENSE` has been added for code under MIT.
- `DATA_LICENSE.md` has been added for processed metadata, tables, source data
  and figure materials under CC BY 4.0.
- Raw reads are outside the declared workflow boundary; do not describe this
  repository as a FASTQ-to-figure pipeline.
- Run `python run_allRevised.py` from the repository root and confirm that
  `outputs/workflow/workflow_status.json` reports `status: success`.
- Confirm that the manuscript and figure files included here are intended for
  public release, not only journal submission.
