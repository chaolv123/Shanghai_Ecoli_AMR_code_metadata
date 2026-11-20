"""Validate and export canonical Supplementary Table 9 and Source Data 1-6.

The script never overwrites the publication-order files. It writes verified
copies and a SHA-256 manifest under outputs/publication_sync.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "publication_sync"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_table9(path: Path) -> dict[str, object]:
    book = pd.ExcelFile(path)
    required = {
        "Annual_ST_ARG_correlations",
        "Isolate_level_ST_ARG_models",
        "ST69_blaCTXM55_models",
        "ST10_ARG_burden_models",
    }
    missing = required - set(book.sheet_names)
    if missing:
        raise ValueError(f"Supplementary Table 9 missing sheets: {sorted(missing)}")
    annual = pd.read_excel(path, sheet_name="Annual_ST_ARG_correlations")
    annual_positive_fdr = annual[
        annual["Spearman_rho"].gt(0) & annual["FDR"].lt(0.05)
    ]
    st69 = pd.read_excel(path, sheet_name="ST69_blaCTXM55_models")
    st10 = pd.read_excel(path, sheet_name="ST10_ARG_burden_models")
    adjusted69 = st69.loc[st69["analysis"].eq("adjusted_host_period_district")].iloc[0]
    adjusted10 = st10.loc[st10["analysis"].eq("negative_binomial_period_host_district")].iloc[0]
    if len(annual_positive_fdr) != 8:
        raise ValueError("Expected eight positive annual ST-ARG rows with FDR < 0.05")
    return {
        "annual_positive_FDR_rows": int(len(annual_positive_fdr)),
        "annual_host_specific_rows": int(annual_positive_fdr["host"].isin(["Poultry", "Swine"]).sum()),
        "ST69_adjusted_OR": float(adjusted69["OR"]),
        "ST69_adjusted_p": float(adjusted69["p_value"]),
        "ST10_adjusted_ratio": float(adjusted10["ratio_or_difference"]),
        "ST10_adjusted_p": float(adjusted10["p_value"]),
    }


def validate_source_data(path: Path, figure: int) -> dict[str, object]:
    book = pd.ExcelFile(path)
    result: dict[str, object] = {"sheets": book.sheet_names}
    if figure == 2:
        panel = pd.read_excel(path, sheet_name="b")
        result["panel_b_rows"] = int(len(panel))
        if len(panel) != 5299:
            raise ValueError("Source Data Fig. 2b must contain 5,299 isolates")
    elif figure == 3:
        panel = pd.read_excel(path, sheet_name="a")
        result["panel_a_rows"] = int(len(panel))
        if len(panel) != 2622:
            raise ValueError("Source Data Fig. 3a must contain 2,622 genomes")
    elif figure == 5:
        localized = pd.read_excel(path, sheet_name="a")
        nonlocalized = pd.read_excel(path, sheet_name="d&f")
        result["localized_isolate_rows"] = int(len(localized))
        result["nonlocalized_isolate_rows"] = int(len(nonlocalized))
    elif figure == 6:
        panel = pd.read_excel(path, sheet_name="a")
        result["PT_ARG_rows"] = int(len(panel))
        result["predicted_PT_count"] = int(panel["PTs"].nunique())
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    publication_inputs = ROOT / "data" / "publication_inputs"
    table9 = publication_inputs / "Supplementary_Table_09.xlsx"
    validation: dict[str, object] = {"Supplementary_Table_09": validate_table9(table9)}
    files = [table9]
    for figure in range(1, 7):
        path = publication_inputs / f"Source_Data_Fig_{figure:02d}.xlsx"
        validation[f"Source_Data_Fig_{figure:02d}"] = validate_source_data(path, figure)
        files.append(path)

    manifest_rows = []
    for source in files:
        target = OUT / source.name
        shutil.copy2(source, target)
        manifest_rows.append(
            {
                "file": source.name,
                "source_relative_path": str(source.relative_to(ROOT)),
                "output_relative_path": str(target.relative_to(ROOT)),
                "sha256": sha256(target),
                "bytes": target.stat().st_size,
            }
        )
    pd.DataFrame(manifest_rows).to_csv(OUT / "publication_sync_manifest.csv", index=False)
    (OUT / "publication_sync_validation.json").write_text(
        json.dumps(validation, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(validation, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
