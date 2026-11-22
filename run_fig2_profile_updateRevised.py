"""Run the R Fig. 2c/Table 4 rebuild, including on non-ASCII Windows paths."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCRIPT = "update_fig2c_supp_table4Revised.R"
INPUTS = ("Source_Data_Fig_02.xlsx", "Supplementary_Table_04.xlsx")


def is_ascii_path(path: Path) -> bool:
    try:
        str(path).encode("ascii")
        return True
    except UnicodeEncodeError:
        return False


def run_r(root: Path) -> None:
    rscript = shutil.which("Rscript")
    if not rscript:
        raise RuntimeError("Rscript was not found; run install_R_dependencies.R first")
    subprocess.run([rscript, SCRIPT], cwd=root, check=True)


def main() -> None:
    if is_ascii_path(ROOT):
        run_r(ROOT)
        return

    # Some older Windows R builds cannot open xlsx files below a non-ASCII path.
    # Stage only the required files in an ASCII-safe temporary directory, run the
    # unchanged R analysis there, and copy its outputs back to the repository.
    with tempfile.TemporaryDirectory(prefix="shanghai_ecoli_fig2_") as tmp_name:
        tmp = Path(tmp_name)
        input_dir = tmp / "data" / "publication_inputs"
        input_dir.mkdir(parents=True)
        shutil.copy2(ROOT / SCRIPT, tmp / SCRIPT)
        for name in INPUTS:
            shutil.copy2(ROOT / "data" / "publication_inputs" / name, input_dir / name)
        run_r(tmp)
        source_output = tmp / "outputs" / "fig2_profile_update"
        target_output = ROOT / "outputs" / "fig2_profile_update"
        target_output.mkdir(parents=True, exist_ok=True)
        for source in source_output.iterdir():
            shutil.copy2(source, target_output / source.name)
    print(f"Fig. 2c/Table 4 outputs written to {target_output}")


if __name__ == "__main__":
    main()
