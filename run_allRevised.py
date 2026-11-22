"""One-command processed-data-to-results workflow for the revision package.

Raw sequencing reads are not required. The workflow starts from the public
processed metadata, ARG calls, PCG/SNP inputs and predicted plasmid-type calls.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
LOG_DIR = ROOT / "outputs" / "workflow"


def run(label: str, command: list[str], log) -> None:
    print(f"\n[{label}] {' '.join(command)}", flush=True)
    log.write(f"\n[{label}] {' '.join(command)}\n")
    process = subprocess.Popen(
        command,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    assert process.stdout is not None
    for line in process.stdout:
        console_encoding = sys.stdout.encoding or "utf-8"
        safe_line = line.encode(console_encoding, errors="replace").decode(
            console_encoding, errors="replace"
        )
        print(safe_line, end="")
        log.write(line)
    code = process.wait()
    if code:
        raise RuntimeError(f"{label} failed with exit code {code}")


def check_r() -> None:
    rscript = shutil.which("Rscript")
    if not rscript:
        raise RuntimeError("Rscript was not found. Install R and run install_R_dependencies.R")
    command = [
        rscript,
        "-e",
        "stopifnot(requireNamespace('glmmTMB', quietly=TRUE), requireNamespace('sandwich', quietly=TRUE), requireNamespace('readxl', quietly=TRUE), requireNamespace('openxlsx', quietly=TRUE), requireNamespace('dplyr', quietly=TRUE), requireNamespace('tidyr', quietly=TRUE), requireNamespace('ggplot2', quietly=TRUE), requireNamespace('svglite', quietly=TRUE), requireNamespace('ragg', quietly=TRUE))",
    ]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(
            "Required R packages glmmTMB and sandwich are unavailable. "
            "Run: Rscript install_R_dependencies.R\n" + result.stderr
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pcoa-permutations",
        type=int,
        default=999,
        help="permutations for the reproducible PCoA/PERMANOVA check",
    )
    args = parser.parse_args()
    check_r()
    rscript = shutil.which("Rscript")
    assert rscript is not None
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / "run_allRevised.log"
    started = datetime.now(timezone.utc).isoformat()
    py = sys.executable
    commands = [
        ("Core metadata/source-data rebuild", [py, "run_pipelineRevised.py", "--check"]),
        ("PCG 494 reconciliation", [py, "run_pcg_reconciledRevised.py"]),
        ("Farm-clustered sensitivity models", [py, "run_farm_clusteredRevised.py"]),
        (
            "Reviewer sensitivity analyses",
            [py, "run_extendedRevised.py", "--analyses", "wgs", "weighted_arg", "breakpoint", "st_arg"],
        ),
        ("PCG sampling-intensity normalization", [py, "run_pcg_sampling_intensityRevised.py"]),
        (
            "Resistome/plasmid derived-input reproducibility",
            [py, "run_derived_analysesRevised.py", "--permutations", str(args.pcoa_permutations)],
        ),
        ("Antibiogram/MDR discrepancy export", [py, "export_antibiogram_mdr_discrepanciesRevised.py"]),
        ("Fig. 2c and Supplementary Table 4 rebuild", [rscript, "update_fig2c_supp_table4Revised.R"]),
        ("Publication-file synchronization", [py, "run_publication_syncRevised.py"]),
        ("Unit tests", [py, "-m", "unittest", "discover", "-s", "tests", "-v"]),
    ]
    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"Started UTC: {started}\n")
        for label, command in commands:
            run(label, command, log)
        completed = datetime.now(timezone.utc).isoformat()
        log.write(f"\nCompleted UTC: {completed}\n")

    summary = {
        "status": "success",
        "started_utc": started,
        "completed_utc": completed,
        "workflow_boundary": "processed metadata/ARG/PCG/predicted-PT inputs; raw reads excluded",
        "steps": [label for label, _ in commands],
        "log": str(log_path.relative_to(ROOT)),
    }
    (LOG_DIR / "workflow_status.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
