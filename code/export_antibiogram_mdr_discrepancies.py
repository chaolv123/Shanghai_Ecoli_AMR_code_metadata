"""Export isolate-level antibiogram-length and MDR discrepancies.

The rule checked here is the manuscript class rule: one resistance class is
counted when any antimicrobial in that class is R or I. MDR is recalculated as
resistance to at least three antimicrobial classes.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "06_Metadata_and_Input_Data"
OUT = ROOT / "outputs"

DRUGS = ["AMP", "AMC", "TIO", "CAZ", "MEM", "SPT", "GEN", "DOX", "TCY", "FF", "SFX", "SXT", "ENR", "OFX", "COL"]
CLASSES = {
    "penicillin": ["AMP"],
    "inhibitor_combination": ["AMC"],
    "cephalosporin": ["TIO", "CAZ"],
    "carbapenem": ["MEM"],
    "aminoglycoside": ["SPT", "GEN"],
    "tetracycline": ["DOX", "TCY"],
    "phenicol": ["FF"],
    "sulfonamide": ["SFX", "SXT"],
    "fluoroquinolone": ["ENR", "OFX"],
    "colistin": ["COL"],
}


def class_reason(row: pd.Series) -> str:
    reasons = []
    for name, drugs in CLASSES.items():
        calls = {drug: row.get(drug) for drug in drugs}
        if any(value in {"R", "I"} for value in calls.values()):
            reasons.append(f"{name}=" + "/".join(f"{drug}:{calls[drug]}" for drug in drugs))
    return "; ".join(reasons)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    phenotype = pd.read_excel(RAW / "phenotype_metadata.xlsx").rename(columns={"Year": "year"})

    class_matrix = pd.DataFrame({
        f"class_{name}": phenotype[drugs].isin(["R", "I"]).any(axis=1).astype(int)
        for name, drugs in CLASSES.items()
    })
    recomputed_length = class_matrix.sum(axis=1)
    recomputed_mdr = recomputed_length.ge(3).astype(int)

    base_columns = ["isolate_id", "year", "host", "period", "districts", "farm_id"]
    optional = [col for col in ["WGS", "WGSisolateid"] if col in phenotype.columns]
    result = phenotype[base_columns + optional + DRUGS + ["antibiogram_length", "MDR_binary"]].copy()
    result = pd.concat([result, class_matrix], axis=1)
    result["recomputed_antibiogram_length"] = recomputed_length
    result["antibiogram_length_difference"] = recomputed_length - phenotype["antibiogram_length"]
    result["recomputed_MDR_binary"] = recomputed_mdr
    result["MDR_difference"] = recomputed_mdr - phenotype["MDR_binary"]
    result["resistant_or_intermediate_classes"] = result.apply(class_reason, axis=1)

    antibiogram = result[result["antibiogram_length_difference"].ne(0)].copy()
    mdr = result[result["MDR_difference"].ne(0)].copy()
    both = result[result["antibiogram_length_difference"].ne(0) | result["MDR_difference"].ne(0)].copy()
    summary = pd.DataFrame([
        {"item": "antibiogram_length_discrepancies", "n": len(antibiogram)},
        {"item": "MDR_binary_discrepancies", "n": len(mdr)},
        {"item": "any_discrepancy", "n": len(both)},
        {"item": "rule", "n": "Class positive if any drug in the class is R or I; MDR if >=3 positive classes."},
    ])

    target = OUT / "antibiogram_MDR_discrepancies_Revised.xlsx"
    with pd.ExcelWriter(target, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="Summary", index=False)
        antibiogram.to_excel(writer, sheet_name="Antibiogram_329", index=False)
        mdr.to_excel(writer, sheet_name="MDR_5", index=False)
        both.to_excel(writer, sheet_name="All_discrepancies", index=False)
    antibiogram.to_csv(OUT / "antibiogram_329_discrepancies_Revised.csv", index=False, encoding="utf-8-sig")
    mdr.to_csv(OUT / "MDR_5_discrepancies_Revised.csv", index=False, encoding="utf-8-sig")
    print(target)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
