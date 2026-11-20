"""Run supplied reviewer-sensitivity scripts against repository-local inputs.

Examples: python run_extendedRevised.py --analyses wgs breakpoint st_arg pcg
The scripts are retained verbatim in archive/legacy_python; this wrapper
overrides only their input/output paths, not their statistical algorithms.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parent
LEGACY = ROOT / "archive" / "legacy_python"
RAW = ROOT / "data" / "raw"
OUT = ROOT / "outputs" / "extended"


def load(name: str):
    path = LEGACY / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def run_wgs():
    module = load("analyze_wgs_sampling_bias")
    module.INPUT = RAW / "phenotype_metadata.xlsx"
    module.OUTDIR = OUT / "wgs"
    module.main()


def run_weighted_arg():
    module = load("analyze_weighted_arg_updated")
    module.WEIGHTS = OUT / "wgs" / "WGS_calibration_weights.csv"
    module.PHENO = RAW / "phenotype_metadata.xlsx"
    module.GENOME = RAW / "genomic_metadata.xlsx"
    module.OUTPUT = OUT / "wgs" / "Supplementary_weighted_ARG_burden_updated.xlsx"
    if not module.WEIGHTS.exists():
        raise FileNotFoundError("Run --analyses wgs before weighted_arg")
    module.main()


def run_breakpoint():
    module = load("analyze_breakpoint_sensitivity_review4")
    module.PHENO_PATH = RAW / "phenotype_metadata.xlsx"
    module.ARG_PATH = RAW / "genomic_metadata.xlsx"
    module.OUT_DIR = OUT / "breakpoint"
    module.OUT_XLSX = module.OUT_DIR / "Supplementary_breakpoint_panel_sensitivity.xlsx"
    module.OUT_SUMMARY = module.OUT_DIR / "summary.txt"
    module.main()


def run_st_arg():
    module = load("analyze_st_arg_association_review6")
    module.IN_PATH = RAW / "genomic_metadata.xlsx"
    module.OUT_DIR = OUT / "st_arg"
    module.OUT_XLSX = module.OUT_DIR / "Supplementary_ST_ARG_association_sensitivity.xlsx"
    module.OUT_SUMMARY = module.OUT_DIR / "summary.txt"
    module.main()


def run_pcg():
    module = load("analyze_pcg_threshold_sensitivity_review5")
    module.META_PATH = RAW / "pcg_metadata.xlsx"
    module.DIST_PATH = RAW / "core_snp_distances.csv"
    module.OUT_DIR = OUT / "pcg"
    module.OUT_XLSX = module.OUT_DIR / "Supplementary_PCG_threshold_sensitivity.xlsx"
    module.OUT_SUMMARY = module.OUT_DIR / "summary.txt"
    module.main()


ANALYSES = {"wgs": run_wgs, "weighted_arg": run_weighted_arg, "breakpoint": run_breakpoint,
            "st_arg": run_st_arg, "pcg": run_pcg}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analyses", choices=list(ANALYSES), nargs="+", required=True)
    args = parser.parse_args()
    for analysis in args.analyses:
        print(f"Running {analysis}...", flush=True)
        ANALYSES[analysis]()
        print(f"Done: {analysis}", flush=True)


if __name__ == "__main__":
    main()
