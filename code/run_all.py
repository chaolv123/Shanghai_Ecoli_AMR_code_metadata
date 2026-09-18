"""Run the public reproducibility workflow.

This script regenerates analysis outputs under ../outputs from the public
metadata files in ../06_Metadata_and_Input_Data.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


CODE_DIR = Path(__file__).resolve().parent


def run(args: list[str]) -> None:
    print("Running:", " ".join(args), flush=True)
    subprocess.run(args, cwd=CODE_DIR.parent, check=True)


def main() -> None:
    python = sys.executable
    run([python, str(CODE_DIR / "run_pipeline.py"), "--check"])
    run([python, str(CODE_DIR / "run_pcg_reconciled.py")])
    run([python, str(CODE_DIR / "run_farm_clustered.py")])
    run([python, str(CODE_DIR / "run_extended.py"), "--analyses", "wgs", "weighted_arg", "breakpoint", "st_arg"])
    run([python, str(CODE_DIR / "export_antibiogram_mdr_discrepancies.py")])


if __name__ == "__main__":
    main()
