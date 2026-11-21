"""Validate and export canonical Supplementary Tables 9/11 and Source Data 1-6.

The script never overwrites the publication-order files. It writes verified
copies and a SHA-256 manifest under outputs/publication_sync.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs" / "publication_sync"


def norm_id(value: object) -> str:
    text = str(value).upper()
    text = re.split(r"(?:[._-]R1|_HQ|-8(?:\.|$)|\.FILTERED)", text, maxsplit=1)[0]
    return re.sub(r"[^A-Z0-9]", "", text)


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
    pooled_rows = int(annual["host"].astype(str).str.casefold().eq("all").sum())
    if pooled_rows:
        raise ValueError("Supplementary Table 9 must contain host-stratified analyses only")
    if len(annual_positive_fdr) != 4:
        raise ValueError("Expected four positive host-stratified annual ST-ARG rows with FDR < 0.05")
    return {
        "annual_positive_FDR_rows": int(len(annual_positive_fdr)),
        "annual_host_specific_rows": int(annual_positive_fdr["host"].isin(["Poultry", "Swine"]).sum()),
        "annual_pooled_rows": pooled_rows,
        "ST69_adjusted_OR": float(adjusted69["OR"]),
        "ST69_adjusted_p": float(adjusted69["p_value"]),
        "ST10_adjusted_ratio": float(adjusted10["ratio_or_difference"]),
        "ST10_adjusted_p": float(adjusted10["p_value"]),
    }


def validate_table11(path: Path) -> dict[str, object]:
    table = pd.read_excel(path)
    required = {
        "Host",
        "Pre-2017, n/N (%)",
        "Post-2017, n/N (%)",
        "Pre-2017 weighted proportion (%)",
        "Post-2017 weighted proportion (%)",
        "Adjusted OR (95% CI)",
        "P value",
    }
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Supplementary Table 11 missing columns: {sorted(missing)}")
    if table["Host"].tolist() != ["Poultry", "Swine"]:
        raise ValueError("Supplementary Table 11 must contain Poultry and Swine rows")
    expected = {
        "Poultry": {"pre": 34.3, "post": 31.8, "p": 0.094},
        "Swine": {"pre": 29.0, "post": 25.8, "p": 0.409},
    }
    for _, row in table.iterrows():
        e = expected[row["Host"]]
        if abs(float(row["Pre-2017 weighted proportion (%)"]) - e["pre"]) > 0.05:
            raise ValueError(f"Unexpected pre-2017 weighted proportion for {row['Host']}")
        if abs(float(row["Post-2017 weighted proportion (%)"]) - e["post"]) > 0.05:
            raise ValueError(f"Unexpected post-2017 weighted proportion for {row['Host']}")
        if abs(float(row["P value"]) - e["p"]) > 0.0005:
            raise ValueError(f"Unexpected adjusted P value for {row['Host']}")
    return {"rows": int(len(table)), "hosts": table["Host"].tolist(), "adjusted_P_values": dict(zip(table["Host"], table["P value"].astype(float)))}


def validate_source_data(path: Path, figure: int) -> dict[str, object]:
    book = pd.ExcelFile(path)
    result: dict[str, object] = {"sheets": book.sheet_names}
    if figure == 1:
        phenotype = pd.read_excel(ROOT / "data" / "raw" / "phenotype_metadata.xlsx")
        panel = pd.read_excel(path, sheet_name="b")
        expected = (
            phenotype.pivot_table(index="year", columns="host", values="isolate_id", aggfunc="count", fill_value=0)
            .reset_index()
            .rename(columns={"year": "Year"})
        )
        expected.columns.name = None
        if not panel.equals(expected):
            raise ValueError("Source Data Fig. 1b does not match phenotype metadata")
        result["panel_b_total"] = int(panel[["Poultry", "Swine"]].to_numpy().sum())
    elif figure == 2:
        panel = pd.read_excel(path, sheet_name="b")
        result["panel_b_rows"] = int(len(panel))
        if len(panel) != 5299:
            raise ValueError("Source Data Fig. 2b must contain 5,299 isolates")
        phenotype_ids = set(pd.read_excel(ROOT / "data" / "raw" / "phenotype_metadata.xlsx")["isolate_id"].astype(str))
        if set(panel["Identity_Code"].astype(str)) != phenotype_ids or panel["Identity_Code"].duplicated().any():
            raise ValueError("Source Data Fig. 2b isolate IDs do not uniquely match phenotype metadata")
    elif figure == 3:
        panel = pd.read_excel(path, sheet_name="a")
        result["panel_a_rows"] = int(len(panel))
        if len(panel) != 2622:
            raise ValueError("Source Data Fig. 3a must contain 2,622 genomes")
        genomic = pd.read_excel(ROOT / "data" / "raw" / "genomic_metadata.xlsx")
        if set(panel["ID"].astype(str)) != set(genomic["isolate_id"].astype(str)) or panel["ID"].duplicated().any():
            raise ValueError("Source Data Fig. 3a isolate IDs do not uniquely match genomic metadata")
        merged = panel.merge(genomic[["isolate_id", "year", "ARG_count", "host"]], left_on="ID", right_on="isolate_id", validate="1:1")
        if not (merged["Year"].eq(merged["year"]) & merged["ARGs"].eq(merged["ARG_count"]) & merged["Host"].eq(merged["host"])).all():
            raise ValueError("Source Data Fig. 3a values do not match genomic metadata")
    elif figure == 5:
        localized = pd.read_excel(path, sheet_name="a")
        nonlocalized = pd.read_excel(path, sheet_name="d&f")
        result["localized_isolate_rows"] = int(len(localized))
        result["nonlocalized_isolate_rows"] = int(len(nonlocalized))
        if localized["PCGs"].nunique() != 308 or nonlocalized["PCGs"].nunique() != 186:
            raise ValueError("Source Data Fig. 5 does not contain 308 localized and 186 non-localized PCGs")
        if set(localized["PCGs"]) & set(nonlocalized["PCGs"]):
            raise ValueError("Source Data Fig. 5 localized/non-localized PCG labels overlap")
        pcg_ids = set(pd.read_excel(ROOT / "data" / "raw" / "pcg_metadata.xlsx")["isolate_id"].map(norm_id))
        source_raw_ids = pd.concat([localized["ID"], nonlocalized["ID"]]).astype(str)
        source_ids = set(source_raw_ids.map(norm_id))
        if not source_ids.issubset(pcg_ids) or len(source_raw_ids) != 1459 or source_raw_ids.duplicated().any():
            raise ValueError("Source Data Fig. 5 isolate IDs do not match the 1,459 PCG-associated isolates")
    elif figure == 6:
        panel = pd.read_excel(path, sheet_name="a")
        burden = pd.read_excel(path, sheet_name="b")
        result["PT_ARG_rows"] = int(len(panel))
        result["predicted_PT_count"] = int(panel["PTs"].nunique())
        if set(zip(panel["ID"], panel["PTs"])) != set(zip(burden["ID"], burden["PTs"])):
            raise ValueError("Source Data Fig. 6a and 6b ID-PT pairs differ")
        matrix_counts = panel.drop(columns=["ID", "PTs"]).apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
        expected = panel[["ID", "PTs"]].assign(expected_ARG_counts=matrix_counts)
        checked = burden.merge(expected, on=["ID", "PTs"], validate="1:1")
        if not pd.to_numeric(checked["ARG_counts"], errors="coerce").eq(checked["expected_ARG_counts"]).all():
            raise ValueError("Source Data Fig. 6b ARG counts do not equal Fig. 6a matrix row sums")
        genomic_ids = set(pd.read_excel(ROOT / "data" / "raw" / "genomic_metadata.xlsx")["isolate_id"].map(norm_id))
        if not set(panel["ID"].map(norm_id)).issubset(genomic_ids):
            raise ValueError("Source Data Fig. 6 contains unmatched isolate IDs")
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    publication_inputs = ROOT / "data" / "publication_inputs"
    table9 = publication_inputs / "Supplementary_Table_09.xlsx"
    table11 = publication_inputs / "Supplementary_Table_11.xlsx"
    validation: dict[str, object] = {
        "Supplementary_Table_09": validate_table9(table9),
        "Supplementary_Table_11": validate_table11(table11),
    }
    files = [table9, table11]
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
