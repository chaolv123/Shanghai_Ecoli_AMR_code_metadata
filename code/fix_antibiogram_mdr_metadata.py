"""Correct antibiogram length and MDR fields in metadata workbooks.

The correction rule is the manuscript class rule:
- one class is positive if any antimicrobial in that class is R or I;
- MDR_binary is 1 if at least three classes are positive.

Only workbooks containing isolate_id, drug columns, antibiogram_length and
MDR_binary are modified.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"

DRUG_CLASSES = {
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

TARGETS = [
    ROOT / "06_Metadata_and_Input_Data" / "phenotype_metadata.xlsx",
]


def recompute(df: pd.DataFrame) -> tuple[pd.Series, pd.Series, pd.DataFrame]:
    class_matrix = pd.DataFrame({
        f"class_{name}": df[drugs].isin(["R", "I"]).any(axis=1).astype(int)
        for name, drugs in DRUG_CLASSES.items()
    })
    length = class_matrix.sum(axis=1)
    mdr = length.ge(3).astype(int)
    return length, mdr, class_matrix


def fix_workbook(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_excel(path)
    if "Year" in df.columns and "year" not in df.columns:
        df = df.rename(columns={"Year": "year"})
    required = {"isolate_id", "antibiogram_length", "MDR_binary", *sum(DRUG_CLASSES.values(), [])}
    if not required.issubset(df.columns):
        return pd.DataFrame()

    new_length, new_mdr, _ = recompute(df)
    mask = df["antibiogram_length"].ne(new_length) | df["MDR_binary"].ne(new_mdr)
    context_columns = [col for col in ["isolate_id", "year", "host", "period", "farm_id", "antibiogram_length", "MDR_binary"] if col in df.columns]
    changes = df.loc[mask, context_columns].copy()
    changes.insert(0, "file", str(path))
    changes["corrected_antibiogram_length"] = new_length.loc[mask].to_numpy()
    changes["corrected_MDR_binary"] = new_mdr.loc[mask].to_numpy()
    changes["antibiogram_length_difference"] = changes["corrected_antibiogram_length"] - changes["antibiogram_length"]
    changes["MDR_difference"] = changes["corrected_MDR_binary"] - changes["MDR_binary"]

    df.loc[:, "antibiogram_length"] = new_length
    df.loc[:, "MDR_binary"] = new_mdr
    df.to_excel(path, index=False)
    return changes


def main() -> None:
    OUT.mkdir(exist_ok=True)
    all_changes = []
    for target in TARGETS:
        changes = fix_workbook(target)
        if not changes.empty:
            all_changes.append(changes)
            print(f"Corrected {len(changes)} rows in {target}")
        else:
            print(f"No correction applied to {target}")
    if all_changes:
        pd.concat(all_changes, ignore_index=True).to_excel(OUT / "metadata_antibiogram_MDR_corrections_Revised.xlsx", index=False)
    else:
        pd.DataFrame().to_excel(OUT / "metadata_antibiogram_MDR_corrections_Revised.xlsx", index=False)


if __name__ == "__main__":
    main()
