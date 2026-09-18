"""Reproduce the supplied Supplementary Table 2 provenance analysis.

This runs the archived Python farm-clustered GLM analysis (sandwich SEs), not
the distinct R GLMM. It writes only under outputs/extended/farm_clustered.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path


CODE_DIR = Path(__file__).resolve().parent
ROOT = CODE_DIR.parents[0]
RAW = ROOT / "06_Metadata_and_Input_Data"
OUT = ROOT / "outputs" / "extended" / "farm_clustered"


def main():
    path = CODE_DIR / "archive" / "legacy_python" / "reanalyze_updated_metadata.py"
    spec = importlib.util.spec_from_file_location("reanalyze_updated_metadata", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    (OUT / "new_add_1").mkdir(parents=True, exist_ok=True)
    (OUT / "new_add_4").mkdir(parents=True, exist_ok=True)
    module.BASE = OUT
    module.PHENO = RAW / "phenotype_metadata.xlsx"
    module.WGS = module.PHENO
    module.ARG = RAW / "genomic_metadata.xlsx"
    module.GENOME = module.ARG
    module.main()


if __name__ == "__main__":
    main()
