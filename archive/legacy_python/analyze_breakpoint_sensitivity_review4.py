from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


BASE = Path("New_analysis") / "new_add_4"
PHENO_PATH = BASE.parent / "原始数据整理和匹配.xlsx"
ARG_PATH = BASE.parent / "new_add_6" / "metadata6.xlsx"
OUT_DIR = BASE / "review4_breakpoint_sensitivity"
OUT_XLSX = OUT_DIR / "Supplementary_breakpoint_panel_sensitivity.xlsx"
OUT_SUMMARY = OUT_DIR / "review4_breakpoint_panel_sensitivity_summary.txt"

DRUG_CLASSES = {
    "AMP": "Penicillin",
    "AMC": "Penicillin-inhibitor",
    "TIO": "Cephalosporin",
    "CAZ": "Cephalosporin",
    "MEM": "Carbapenem",
    "SPT": "Aminoglycoside",
    "GEN": "Aminoglycoside",
    "DOX": "Tetracycline",
    "TCY": "Tetracycline",
    "FF": "Phenicol",
    "SFX": "Sulfonamide",
    "SXT": "Sulfonamide",
    "ENR": "Fluoroquinolone",
    "OFX": "Fluoroquinolone",
    "COL": "Colistin",
}

DRUGS = list(DRUG_CLASSES)
CUTOFFS = [2016, 2017, 2018]


def bh_fdr(p_values: pd.Series) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.dropna()
    if valid.empty:
        return out
    order = valid.sort_values().index
    ranked = valid.loc[order].to_numpy()
    m = len(ranked)
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out.loc[order] = adj
    return out


def clean_result(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip().upper()
    # Treat intermediate as non-susceptible for binary resistance analyses.
    if s in {"R", "I", "RESISTANT", "INTERMEDIATE", "1", "TRUE"}:
        return 1
    if s in {"S", "SUSCEPTIBLE", "0", "FALSE"}:
        return 0
    return np.nan


def binary_prepost(df: pd.DataFrame, outcome: str, cutoff: int, host: str) -> dict:
    sub = df if host == "All" else df[df["host"] == host]
    sub = sub[["year", outcome]].dropna()
    pre = sub[sub["year"] < cutoff][outcome].astype(float)
    post = sub[sub["year"] >= cutoff][outcome].astype(float)
    row = {
        "host": host,
        "outcome": outcome,
        "transition_year": cutoff,
        "pre_years": f"{int(sub['year'].min())}-{cutoff-1}" if len(pre) else None,
        "post_years": f"{cutoff}-{int(sub['year'].max())}" if len(post) else None,
        "n_pre": int(len(pre)),
        "n_post": int(len(post)),
        "pre_value": pre.mean() if len(pre) else np.nan,
        "post_value": post.mean() if len(post) else np.nan,
        "change_post_minus_pre": (post.mean() - pre.mean()) if len(pre) and len(post) else np.nan,
        "relative_change": ((post.mean() - pre.mean()) / pre.mean()) if len(pre) and len(post) and pre.mean() != 0 else np.nan,
        "test": "chi_square_or_fisher",
        "p_value": np.nan,
    }
    if len(pre) and len(post):
        table = np.array([
            [int(pre.sum()), int(len(pre) - pre.sum())],
            [int(post.sum()), int(len(post) - post.sum())],
        ])
        if (table < 5).any():
            row["test"] = "fisher_exact"
            row["p_value"] = stats.fisher_exact(table).pvalue
        else:
            row["p_value"] = stats.chi2_contingency(table, correction=False).pvalue
    return row


def continuous_prepost(df: pd.DataFrame, outcome: str, cutoff: int, host: str) -> dict:
    sub = df if host == "All" else df[df["host"] == host]
    sub = sub[["year", outcome]].dropna()
    pre = sub[sub["year"] < cutoff][outcome].astype(float)
    post = sub[sub["year"] >= cutoff][outcome].astype(float)
    row = {
        "host": host,
        "outcome": outcome,
        "transition_year": cutoff,
        "pre_years": f"{int(sub['year'].min())}-{cutoff-1}" if len(pre) else None,
        "post_years": f"{cutoff}-{int(sub['year'].max())}" if len(post) else None,
        "n_pre": int(len(pre)),
        "n_post": int(len(post)),
        "pre_value": pre.mean() if len(pre) else np.nan,
        "post_value": post.mean() if len(post) else np.nan,
        "change_post_minus_pre": (post.mean() - pre.mean()) if len(pre) and len(post) else np.nan,
        "relative_change": ((post.mean() - pre.mean()) / pre.mean()) if len(pre) and len(post) and pre.mean() != 0 else np.nan,
        "test": "mann_whitney_u",
        "p_value": np.nan,
    }
    if len(pre) and len(post):
        row["p_value"] = stats.mannwhitneyu(pre, post, alternative="two-sided").pvalue
    return row


def make_class_variables(df: pd.DataFrame, drug_list: list[str], prefix: str) -> tuple[pd.DataFrame, list[str]]:
    out = df.copy()
    classes = []
    for cls in sorted({DRUG_CLASSES[d] for d in drug_list}):
        ds = [d for d in drug_list if DRUG_CLASSES[d] == cls]
        tested = out[ds].notna().any(axis=1)
        resistant = out[ds].eq(1).any(axis=1)
        col = f"{prefix}_class_{cls}"
        out[col] = np.where(tested, resistant.astype(int), np.nan)
        classes.append(col)
    out[f"{prefix}_resistant_class_count"] = out[classes].sum(axis=1, min_count=1)
    out[f"{prefix}_MDR_by_class"] = np.where(
        out[f"{prefix}_resistant_class_count"].notna(),
        (out[f"{prefix}_resistant_class_count"] >= 3).astype(int),
        np.nan,
    )
    out[f"{prefix}_drug_resistance_count"] = out[drug_list].sum(axis=1, min_count=1)
    return out, classes


def summarize_directions(df: pd.DataFrame, outcome_col: str = "change_post_minus_pre") -> pd.DataFrame:
    rows = []
    for (host, outcome), g in df.groupby(["host", "outcome"], dropna=False):
        g = g.dropna(subset=[outcome_col])
        if g.empty:
            continue
        ref = g[g["transition_year"] == 2017]
        if ref.empty or ref[outcome_col].iloc[0] == 0:
            ref_sign = np.nan
        else:
            ref_sign = math.copysign(1, ref[outcome_col].iloc[0])
        same = []
        for _, r in g.iterrows():
            sign = np.nan if r[outcome_col] == 0 else math.copysign(1, r[outcome_col])
            same.append((not pd.isna(ref_sign)) and sign == ref_sign)
        rows.append({
            "host": host,
            "outcome": outcome,
            "n_cutoffs_evaluable": len(g),
            "change_2016": g.loc[g["transition_year"] == 2016, outcome_col].iloc[0] if (g["transition_year"] == 2016).any() else np.nan,
            "change_2017": g.loc[g["transition_year"] == 2017, outcome_col].iloc[0] if (g["transition_year"] == 2017).any() else np.nan,
            "change_2018": g.loc[g["transition_year"] == 2018, outcome_col].iloc[0] if (g["transition_year"] == 2018).any() else np.nan,
            "same_direction_as_2017_all_available": all(same) if same else np.nan,
        })
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pheno = pd.read_excel(PHENO_PATH).rename(columns={"Year": "year"})
    arg = pd.read_excel(ARG_PATH)[["isolate_id", "host", "year", "period", "districts", "farm_id", "ST", "ARG_count"]]
    for d in DRUGS:
        pheno[d] = pheno[d].map(clean_result)
    pheno["MDR_binary"] = pd.to_numeric(pheno["MDR_binary"], errors="coerce")
    pheno["antibiogram_length"] = pd.to_numeric(pheno["antibiogram_length"], errors="coerce")
    pheno["year"] = pd.to_numeric(pheno["year"], errors="coerce").astype(int)
    arg["year"] = pd.to_numeric(arg["year"], errors="coerce").astype(int)
    arg["ARG_count"] = pd.to_numeric(arg["ARG_count"], errors="coerce")

    years_all = set(range(int(pheno["year"].min()), int(pheno["year"].max()) + 1))
    panel_rows = []
    for d in DRUGS:
        tested = pheno.dropna(subset=[d])
        years = sorted(tested["year"].unique())
        panel_rows.append({
            "antimicrobial": d,
            "class": DRUG_CLASSES[d],
            "n_tested": int(len(tested)),
            "first_year_tested": int(min(years)) if years else np.nan,
            "last_year_tested": int(max(years)) if years else np.nan,
            "n_years_tested": len(years),
            "years_tested": ",".join(map(str, years)),
            "available_in_every_2008_2023_year": set(years) == years_all,
        })
    panel = pd.DataFrame(panel_rows)
    core_drugs = panel.loc[panel["available_in_every_2008_2023_year"], "antimicrobial"].tolist()

    pheno, all_class_cols = make_class_variables(pheno, DRUGS, "all_available")
    pheno, core_class_cols = make_class_variables(pheno, core_drugs, "core_panel")

    hosts = ["All"] + sorted(pheno["host"].dropna().unique().tolist())

    rows = []
    # Individual antimicrobial trends; later-introduced drugs become not evaluable where pre or post is absent.
    for cutoff in CUTOFFS:
        for host in hosts:
            for d in DRUGS:
                r = binary_prepost(pheno, d, cutoff, host)
                r["analysis_family"] = "individual_antimicrobial"
                r["panel"] = "observed_available"
                rows.append(r)
            for outcome in ["MDR_binary", "all_available_MDR_by_class"]:
                r = binary_prepost(pheno, outcome, cutoff, host)
                r["analysis_family"] = "MDR"
                r["panel"] = "all_available"
                rows.append(r)
            for outcome in ["antibiogram_length", "all_available_drug_resistance_count", "all_available_resistant_class_count"]:
                r = continuous_prepost(pheno, outcome, cutoff, host)
                r["analysis_family"] = "burden"
                r["panel"] = "all_available"
                rows.append(r)
            for outcome in ["core_panel_MDR_by_class"]:
                r = binary_prepost(pheno, outcome, cutoff, host)
                r["analysis_family"] = "MDR"
                r["panel"] = "fixed_core_panel"
                rows.append(r)
            for outcome in ["core_panel_drug_resistance_count", "core_panel_resistant_class_count"]:
                r = continuous_prepost(pheno, outcome, cutoff, host)
                r["analysis_family"] = "burden"
                r["panel"] = "fixed_core_panel"
                rows.append(r)
            for outcome in all_class_cols:
                r = binary_prepost(pheno, outcome, cutoff, host)
                r["analysis_family"] = "phenotypic_class"
                r["panel"] = "all_available"
                rows.append(r)
            for outcome in core_class_cols:
                r = binary_prepost(pheno, outcome, cutoff, host)
                r["analysis_family"] = "phenotypic_class"
                r["panel"] = "fixed_core_panel"
                rows.append(r)

    pheno_results = pd.DataFrame(rows)
    pheno_results["FDR"] = pheno_results.groupby(["analysis_family", "panel", "transition_year"])["p_value"].transform(bh_fdr)

    arg_rows = []
    for cutoff in CUTOFFS:
        for host in ["All"] + sorted(arg["host"].dropna().unique().tolist()):
            r = continuous_prepost(arg, "ARG_count", cutoff, host)
            r["analysis_family"] = "ARG_burden"
            r["panel"] = "WGS_ARG_count"
            arg_rows.append(r)
    arg_results = pd.DataFrame(arg_rows)
    arg_results["FDR"] = arg_results.groupby(["analysis_family", "transition_year"])["p_value"].transform(bh_fdr)

    pheno_direction = summarize_directions(
        pheno_results[
            pheno_results["outcome"].isin([
                "MDR_binary",
                "antibiogram_length",
                "all_available_MDR_by_class",
                "all_available_drug_resistance_count",
                "all_available_resistant_class_count",
                "core_panel_MDR_by_class",
                "core_panel_drug_resistance_count",
                "core_panel_resistant_class_count",
            ])
        ]
    )
    arg_direction = summarize_directions(arg_results)

    key = pheno_results[
        pheno_results["outcome"].isin([
            "MDR_binary",
            "antibiogram_length",
            "all_available_MDR_by_class",
            "all_available_drug_resistance_count",
            "all_available_resistant_class_count",
            "core_panel_MDR_by_class",
            "core_panel_drug_resistance_count",
            "core_panel_resistant_class_count",
        ])
    ].copy()

    annual_rows = []
    annual_outcomes = ["MDR_binary", "antibiogram_length", "core_panel_MDR_by_class", "core_panel_resistant_class_count"]
    for host in hosts:
        sub_host = pheno if host == "All" else pheno[pheno["host"] == host]
        for outcome in annual_outcomes:
            for year, g in sub_host.groupby("year"):
                vals = g[outcome].dropna().astype(float)
                annual_rows.append({
                    "dataset": "phenotype",
                    "host": host,
                    "outcome": outcome,
                    "year": int(year),
                    "n": int(len(vals)),
                    "value": vals.mean() if len(vals) else np.nan,
                })
    for host in ["All"] + sorted(arg["host"].dropna().unique().tolist()):
        sub_host = arg if host == "All" else arg[arg["host"] == host]
        for year, g in sub_host.groupby("year"):
            vals = g["ARG_count"].dropna().astype(float)
            annual_rows.append({
                "dataset": "ARG",
                "host": host,
                "outcome": "ARG_count",
                "year": int(year),
                "n": int(len(vals)),
                "value": vals.mean() if len(vals) else np.nan,
            })
    annual = pd.DataFrame(annual_rows)

    validation = pd.DataFrame([
        {"check": "phenotype_rows", "observed": len(pheno), "expected": 5299, "status": "PASS" if len(pheno) == 5299 else "CHECK"},
        {"check": "ARG_rows", "observed": len(arg), "expected": 2622, "status": "PASS" if len(arg) == 2622 else "CHECK"},
        {"check": "ARG_columns", "observed": ",".join(arg.columns), "expected": "ARG_count extracted from current metadata6.xlsx", "status": "PASS"},
        {"check": "core_drugs", "observed": ",".join(core_drugs), "expected": "drugs with data in every year 2008-2023", "status": "PASS"},
    ])

    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
        validation.to_excel(writer, sheet_name="Data_validation", index=False)
        panel.to_excel(writer, sheet_name="Panel_availability", index=False)
        pd.DataFrame({"core_drug": core_drugs, "class": [DRUG_CLASSES[d] for d in core_drugs]}).to_excel(
            writer, sheet_name="Core_panel_definition", index=False
        )
        key.to_excel(writer, sheet_name="Key_phenotype_sensitivity", index=False)
        pheno_results.to_excel(writer, sheet_name="All_phenotype_sensitivity", index=False)
        pheno_direction.to_excel(writer, sheet_name="Phenotype_direction_summary", index=False)
        arg_results.to_excel(writer, sheet_name="ARG_sensitivity", index=False)
        arg_direction.to_excel(writer, sheet_name="ARG_direction_summary", index=False)
        annual.to_excel(writer, sheet_name="Annual_observed_values", index=False)

    # Compact text summary for manuscript drafting.
    def get_change(table, host, outcome, cutoff):
        row = table[(table["host"] == host) & (table["outcome"] == outcome) & (table["transition_year"] == cutoff)]
        if row.empty:
            return np.nan
        return row["change_post_minus_pre"].iloc[0]

    lines = []
    lines.append("Review 4 breakpoint and panel sensitivity analysis")
    lines.append("=" * 58)
    lines.append(f"Phenotype table: n={len(pheno)}; years={pheno.year.min()}-{pheno.year.max()}")
    lines.append(f"ARG table: n={len(arg)}; years={arg.year.min()}-{arg.year.max()}; columns={list(arg.columns)}")
    lines.append("")
    lines.append("Panel availability")
    lines.append("------------------")
    for _, r in panel.iterrows():
        lines.append(
            f"{r['antimicrobial']} ({r['class']}): n={r['n_tested']}, years={r['years_tested']}, "
            f"core={r['available_in_every_2008_2023_year']}"
        )
    lines.append("")
    lines.append(f"Fixed core panel drugs: {', '.join(core_drugs)}")
    lines.append("Core panel excludes CAZ and MEM because they were introduced in 2016; excludes DOX because it is available only from 2010 to 2017; excludes SPT and OFX because one early year is absent.")
    lines.append("")
    lines.append("Key 2017 comparisons")
    lines.append("--------------------")
    for host in ["Poultry", "Swine"]:
        for outcome in ["MDR_binary", "antibiogram_length", "core_panel_MDR_by_class", "core_panel_resistant_class_count"]:
            row = key[(key["host"] == host) & (key["outcome"] == outcome) & (key["transition_year"] == 2017)]
            if not row.empty:
                rr = row.iloc[0]
                lines.append(
                    f"{host} {outcome}: pre={rr['pre_value']:.4g}, post={rr['post_value']:.4g}, "
                    f"change={rr['change_post_minus_pre']:.4g}, P={rr['p_value']:.3g}"
                )
        row = arg_results[(arg_results["host"] == host) & (arg_results["transition_year"] == 2017)].iloc[0]
        lines.append(
            f"{host} ARG_count: pre={row['pre_value']:.4g}, post={row['post_value']:.4g}, "
            f"change={row['change_post_minus_pre']:.4g}, P={row['p_value']:.3g}"
        )
    lines.append("")
    lines.append("Direction robustness across 2016/2017/2018 cutoffs")
    lines.append("---------------------------------------------------")
    for host in ["Poultry", "Swine"]:
        for outcome in ["MDR_binary", "antibiogram_length", "core_panel_MDR_by_class", "core_panel_resistant_class_count"]:
            ds = pheno_direction[(pheno_direction["host"] == host) & (pheno_direction["outcome"] == outcome)]
            if not ds.empty:
                r = ds.iloc[0]
                lines.append(
                    f"{host} {outcome}: changes 2016/2017/2018 = "
                    f"{r['change_2016']:.4g}/{r['change_2017']:.4g}/{r['change_2018']:.4g}; "
                    f"same direction={r['same_direction_as_2017_all_available']}"
                )
        ds = arg_direction[(arg_direction["host"] == host) & (arg_direction["outcome"] == "ARG_count")]
        if not ds.empty:
            r = ds.iloc[0]
            lines.append(
                f"{host} ARG_count: changes 2016/2017/2018 = "
                f"{r['change_2016']:.4g}/{r['change_2017']:.4g}/{r['change_2018']:.4g}; "
                f"same direction={r['same_direction_as_2017_all_available']}"
            )
    lines.append("")
    lines.append("Limitation")
    lines.append("----------")
    lines.append("This workbook analyzes total ARG_count; class-level ARG sensitivity is not included here although metadata6.xlsx contains individual ARG columns.")
    OUT_SUMMARY.write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote {OUT_XLSX}")
    print(f"Wrote {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
