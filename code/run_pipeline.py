"""Rebuild auditable metadata-derived tables and source-data panels.

Usage: python run_pipelineRevised.py [--check]
Input files live under data/raw. This script never modifies those files.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "06_Metadata_and_Input_Data"
OUT = ROOT / "outputs"
DRUGS = ["AMP", "AMC", "TIO", "CAZ", "MEM", "SPT", "GEN", "DOX", "TCY", "FF", "SFX", "SXT", "ENR", "OFX", "COL"]
LIMITED = {"CAZ", "MEM", "DOX"}


def norm_id(value: object) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(value).upper())


def bh(p: pd.Series) -> pd.Series:
    valid = p.dropna().sort_values()
    result = pd.Series(np.nan, index=p.index)
    if len(valid):
        adjusted = valid.to_numpy() * len(valid) / np.arange(1, len(valid) + 1)
        result.loc[valid.index] = np.clip(np.minimum.accumulate(adjusted[::-1])[::-1], 0, 1)
    return result


def test_binary(pre_r: int, pre_n: int, post_r: int, post_n: int) -> tuple[str, float]:
    if not pre_n or not post_n:
        return "Not estimable", np.nan
    table = np.array([[pre_r, pre_n - pre_r], [post_r, post_n - post_r]])
    if (table.sum(axis=0) == 0).any():
        return "Not estimable", np.nan
    expected = np.outer(table.sum(axis=1), table.sum(axis=0)) / table.sum()
    if (expected < 5).any():
        return "Fisher exact", float(stats.fisher_exact(table).pvalue)
    return "Chi-square", float(stats.chi2_contingency(table, correction=False).pvalue)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="validate inputs and outputs")
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    phenotype = pd.read_excel(RAW / "phenotype_metadata.xlsx").rename(columns={"Year": "year"})
    genomic = pd.read_excel(RAW / "genomic_metadata.xlsx").rename(columns={"Year": "year"})
    if len(phenotype) != 5299 or not phenotype.isolate_id.is_unique:
        raise ValueError("Expected 5,299 unique phenotyped isolates")
    if len(genomic) != 2622 or not genomic.isolate_id.is_unique:
        raise ValueError("Expected 2,622 unique sequenced isolates")
    if not phenotype.period.eq(np.where(phenotype.year < 2017, "pre-2017", "post-2017")).all():
        raise ValueError("Phenotype period labels disagree with year")

    rows = []
    annual = []
    for host in ["Poultry", "Swine"]:
        h = phenotype[phenotype.host.eq(host)]
        for drug in DRUGS:
            observed = h[h[drug].isin(["R", "I", "S"])].copy()
            observed["non_susceptible"] = observed[drug].isin(["R", "I"])
            pre = observed[observed.year.lt(2017)]
            post = observed[observed.year.ge(2017)]
            pre_r, post_r = int(pre.non_susceptible.sum()), int(post.non_susceptible.sum())
            test, p = test_binary(pre_r, len(pre), post_r, len(post)) if drug not in LIMITED else ("Not tested (limited panel)", np.nan)
            rows.append({"Host": host, "Antimicrobial": drug, "n_pre": len(pre), "resistant_pre": pre_r,
                         "rate_pre_pct": 100 * pre_r / len(pre) if len(pre) else np.nan,
                         "n_post": len(post), "resistant_post": post_r,
                         "rate_post_pct": 100 * post_r / len(post) if len(post) else np.nan,
                         "test": test, "P_value": p,
                         "years_pre": ",".join(map(str, sorted(pre.year.unique()))),
                         "years_post": ",".join(map(str, sorted(post.year.unique())))})
            for year, year_data in observed.groupby("year"):
                annual.append({"Year": int(year), "Host": host, "Percentage": 100 * year_data.non_susceptible.mean(),
                               "Antimicroials": drug, "n_tested": len(year_data)})
    table = pd.DataFrame(rows)
    table["q_value_BH"] = table.groupby("Host")["P_value"].transform(bh)
    table.to_excel(OUT / "supplementary_table_1_recomputed.xlsx", index=False)
    pd.DataFrame(annual).to_excel(OUT / "fig1d_annual_resistance_recomputed.xlsx", index=False)

    fig1b = phenotype.groupby(["year", "host"]).size().unstack(fill_value=0).reset_index()
    fig1b.columns.name = None
    fig1b = fig1b.rename(columns={"year": "Year"})[["Year", "Poultry", "Swine"]]
    fig1b.to_excel(OUT / "fig1b_host_year_counts_recomputed.xlsx", index=False)
    table.rename(columns={"Antimicrobial": "Agents", "rate_pre_pct": "Pre2017", "rate_post_pct": "Post2017"})[["Agents", "Pre2017", "Post2017", "Host"]].to_excel(OUT / "fig1c_period_resistance_recomputed.xlsx", index=False)

    fig2b = phenotype.reset_index(drop=True).reset_index().rename(columns={"index": "Code"})
    fig2b["Code"] += 1
    fig2b = fig2b.rename(columns={"year": "Year", "districts": "Districts", "host": "Host", "farm_id": "Farms",
                                 "isolate_id": "Identity_Code", "antibiogram_length": "Antibiogram length", "MDR_binary": "MDR"})
    fig2b[["Code", "Year", "Districts", "Host", "Farms", "Identity_Code", "Antibiogram length", "MDR"]].to_excel(OUT / "fig2b_isolate_metadata_recomputed.xlsx", index=False)

    genomic[["year", "isolate_id", "ARG_count", "host"]].rename(columns={"year": "Year", "isolate_id": "ID", "ARG_count": "ARGs", "host": "Host"}).to_excel(OUT / "fig3a_ARG_burden_recomputed.xlsx", index=False)
    classes = {
        "penicillin": ["AMP"], "inhibitor": ["AMC"], "cephalosporin": ["TIO", "CAZ"],
        "carbapenem": ["MEM"], "aminoglycoside": ["SPT", "GEN"],
        "tetracycline": ["DOX", "TCY"], "phenicol": ["FF"],
        "sulfonamide": ["SFX", "SXT"], "fluoroquinolone": ["ENR", "OFX"],
        "colistin": ["COL"],
    }
    class_matrix = pd.DataFrame({name: phenotype[drugs].isin(["R", "I"]).any(axis=1).astype(int)
                                 for name, drugs in classes.items()})
    recomputed_length = class_matrix.sum(axis=1)
    class_check = phenotype[["isolate_id", "year", "host", "period", "ENR", "OFX", "antibiogram_length", "MDR_binary"]].copy()
    class_check["length_any_R_or_I_in_class"] = recomputed_length
    class_check["length_difference"] = recomputed_length - phenotype.antibiogram_length
    class_check["MDR_any_R_or_I_in_class"] = recomputed_length.ge(3).astype(int)
    class_check = class_check[class_check.length_difference.ne(0) |
                              class_check.MDR_any_R_or_I_in_class.ne(class_check.MDR_binary)]
    class_check.to_csv(OUT / "antibiogram_class_definition_discrepancies.csv", index=False, encoding="utf-8-sig")
    check = {
        "phenotype_n": len(phenotype), "poultry_n": int(phenotype.host.eq("Poultry").sum()),
        "swine_n": int(phenotype.host.eq("Swine").sum()), "genomic_n": len(genomic),
        "wgs_linked_n": int(phenotype.WGS.eq(1).sum()),
        "ARG_matrix_count_mismatch": int((genomic.iloc[:, 8:].apply(pd.to_numeric, errors="coerce").sum(axis=1) != genomic.ARG_count).sum()),
        "tested_comparisons_per_host": table.groupby("Host").P_value.count().to_dict(),
        "significant_q_lt_0_05": table.assign(sig=table.q_value_BH.lt(.05)).groupby("Host").sig.sum().astype(int).to_dict(),
        "class_length_rule_discrepancies": len(class_check),
        "MDR_rule_discrepancies": int(class_check.MDR_any_R_or_I_in_class.ne(class_check.MDR_binary).sum()),
    }
    (OUT / "validation.json").write_text(json.dumps(check, indent=2), encoding="utf-8")
    if args.check and (check["wgs_linked_n"] != 2547 or check["ARG_matrix_count_mismatch"]):
        raise ValueError(check)
    print(json.dumps(check, indent=2))


if __name__ == "__main__":
    main()
