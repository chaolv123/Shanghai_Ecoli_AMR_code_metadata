from __future__ import annotations

import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "New_analysis"
PHENO = BASE / "原始数据整理和匹配.xlsx"
ARG = BASE / "new_add_6" / "metadata6.xlsx"
GENOME = BASE / "new_add_6" / "metadata6.xlsx"
WGS = BASE / "原始数据整理和匹配.xlsx"
DRUGS = ["AMP", "AMC", "TIO", "CAZ", "MEM", "SPT", "GEN", "DOX", "TCY", "FF", "SFX", "SXT", "ENR", "OFX", "COL"]


def fdr_bh(series: pd.Series) -> pd.Series:
    p = pd.to_numeric(series, errors="coerce")
    result = pd.Series(np.nan, index=series.index)
    valid = p.dropna().sort_values()
    if valid.empty:
        return result
    m = len(valid)
    adjusted = valid.to_numpy() * m / np.arange(1, m + 1)
    adjusted = np.clip(np.minimum.accumulate(adjusted[::-1])[::-1], 0, 1)
    result.loc[valid.index] = adjusted
    return result


def read_tables():
    pheno = pd.read_excel(PHENO).rename(columns={"Year": "year"})
    arg = pd.read_excel(ARG)[["isolate_id", "host", "year", "period", "districts", "farm_id", "ST", "ARG_count"]]
    genome = pd.read_excel(GENOME)
    wgs = pd.read_excel(WGS).rename(columns={"Year": "year"})
    for table in (pheno, arg, genome, wgs):
        table["year"] = pd.to_numeric(table["year"], errors="coerce")
    if not pheno["isolate_id"].is_unique or not genome["isolate_id"].is_unique:
        raise ValueError("Duplicate isolate_id in updated phenotype or genome data")
    return pheno, arg, genome, wgs


def comparison_table(pheno):
    rows = []
    for host in ["Overall", "Poultry", "Swine"]:
        source = pheno if host == "Overall" else pheno[pheno["host"] == host]
        for drug in DRUGS:
            pre = source.loc[source["year"] < 2017, drug].astype(str).str.upper()
            post = source.loc[source["year"] >= 2017, drug].astype(str).str.upper()
            pre = pre[pre.isin(["R", "S", "I"])]
            post = post[post.isin(["R", "S", "I"])]
            a, b = int(pre.isin(["R", "I"]).sum()), int(pre.eq("S").sum())
            c, d = int(post.isin(["R", "I"]).sum()), int(post.eq("S").sum())
            n_pre, n_post = a + b, c + d
            p_value, test = np.nan, "not_evaluable"
            if n_pre and n_post:
                table = np.array([[a, b], [c, d]])
                expected = np.outer(table.sum(axis=1), table.sum(axis=0)) / table.sum()
                if (expected < 5).any():
                    test, p_value = "fisher_exact", stats.fisher_exact(table).pvalue
                else:
                    test, p_value = "chi_square", stats.chi2_contingency(table, correction=False).pvalue
            rows.append({
                "group": host, "drug": drug, "n_pre": n_pre, "resistant_pre": a,
                "resistance_pre_%": 100 * a / n_pre if n_pre else np.nan,
                "n_post": n_post, "resistant_post": c,
                "resistance_post_%": 100 * c / n_post if n_post else np.nan,
                "post_minus_pre_pp": 100 * (c / n_post - a / n_pre) if n_pre and n_post else np.nan,
                "test": test, "p_value": p_value,
            })
    out = pd.DataFrame(rows)
    out["q_value_BH"] = out.groupby("group")["p_value"].transform(fdr_bh)
    return out


def model_row(source, host, outcome, subset, is_binary):
    data = source.copy()
    data["value"] = data[outcome].map({"R": 1, "I": 1, "S": 0}) if outcome in DRUGS else pd.to_numeric(data[outcome], errors="coerce")
    data = data.dropna(subset=["value", "year", "farm_id", "districts", "host"])
    data["farm_id"] = data["farm_id"].astype(str)
    data["period_post"] = (data["year"] >= 2017).astype(int)
    data["year_centered"] = data["year"] - 2017
    pre, post = data.loc[data.period_post.eq(0), "value"], data.loc[data.period_post.eq(1), "value"]
    row = {
        "subset": subset, "host": host, "outcome": outcome,
        "outcome_type": "binary" if is_binary else "count",
        "effect_measure": "adjusted OR" if is_binary else "adjusted rate ratio",
        "n": len(data), "n_pre": len(pre), "n_post": len(post),
        "n_farms": data.farm_id.nunique(), "n_districts": data.districts.nunique(),
        "event_or_mean_pre": pre.mean(), "event_or_mean_post": post.mean(),
        "model_family": "binomial_logit" if is_binary else "poisson_cluster_robust",
        "key_term": "period_post", "estimate": np.nan,
        "conf_low": np.nan, "conf_high": np.nan, "p_value": np.nan, "note": "",
    }
    if len(pre) < 10 or len(post) < 10 or data.farm_id.nunique() < 3 or data.value.nunique() < 2:
        row["note"] = "Insufficient observations or variation"
        return row
    covariates = "period_post + year_centered + C(districts)"
    if host == "All":
        covariates += " + C(host)"
    family = sm.families.Binomial() if is_binary else sm.families.Poisson()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = smf.glm(f"value ~ {covariates}", data=data, family=family).fit(
                cov_type="cluster", cov_kwds={"groups": data["farm_id"]}, maxiter=200
            )
        beta, se = model.params["period_post"], model.bse["period_post"]
        if not np.isfinite(beta) or not np.isfinite(se):
            raise ValueError("nonfinite coefficient")
        row.update({"estimate": math.exp(beta), "conf_low": math.exp(beta - 1.96 * se),
                    "conf_high": math.exp(beta + 1.96 * se), "p_value": model.pvalues["period_post"]})
    except (ValueError, OverflowError, np.linalg.LinAlgError) as exc:
        row["note"] = f"Model not estimable: {type(exc).__name__}"
    return row


def farm_models(pheno, arg):
    rows = []
    for subset in ["All farms", "Farms sampled >5 years"]:
        for dataset, outcomes in [(pheno, DRUGS + ["MDR_binary", "antibiogram_length"]), (arg, ["ARG_count"])]:
            source = dataset.copy()
            if subset != "All farms":
                eligible = source.groupby("farm_id")["year"].nunique()
                source = source[source.farm_id.isin(eligible[eligible > 5].index)]
            for host in ["All", "Poultry", "Swine"]:
                part = source if host == "All" else source[source.host.eq(host)]
                for outcome in outcomes:
                    rows.append(model_row(part, host, outcome, subset, outcome in DRUGS or outcome == "MDR_binary"))
    return pd.DataFrame(rows)


def validation(pheno, arg, genome, wgs):
    wgs_ids = set(wgs.loc[wgs.WGS.eq(1), "WGSisolateid"].dropna().astype(str))
    genome_ids = set(genome.isolate_id.astype(str))
    arg_ids = set(arg.isolate_id.astype(str))
    linked = wgs[wgs.WGS.eq(1)][["WGSisolateid", "year", "host", "districts", "farm_id"]].merge(
        genome[["isolate_id", "year", "host", "districts", "farm_id"]], left_on="WGSisolateid", right_on="isolate_id", suffixes=("_pheno", "_genome")
    )
    facts = {
        "phenotype_rows": len(pheno), "WGS_phenotype_rows": int(wgs.WGS.eq(1).sum()),
        "unique_WGS_ids": len(wgs_ids), "genome_rows": len(genome), "ARG_rows": len(arg),
        "WGS_to_genome_matched": len(wgs_ids & genome_ids),
        "genome_without_phenotype": len(genome_ids - wgs_ids),
        "phenotype_WGS_without_genome": len(wgs_ids - genome_ids),
        "ARG_vs_genome_id_difference": len(arg_ids ^ genome_ids),
        "linked_host_disagreements": int(linked.host_pheno.ne(linked.host_genome).sum()),
        "linked_year_disagreements": int(linked.year_pheno.ne(linked.year_genome).sum()),
        "linked_farm_disagreements": int(linked.farm_id_pheno.ne(linked.farm_id_genome).sum()),
        "linked_district_disagreements": int(linked.districts_pheno.ne(linked.districts_genome).sum()),
        "intermediate_test_results": int(pheno[DRUGS].eq("I").sum().sum()),
    }
    discrepancies = linked[linked.year_pheno.ne(linked.year_genome) | linked.farm_id_pheno.ne(linked.farm_id_genome)].copy()
    return pd.DataFrame([{"check": k, "observed": v} for k, v in facts.items()]), discrepancies


def arg_farm_sensitivity(arg, wgs):
    matched_farms = wgs.loc[wgs.WGS.eq(1), ["WGSisolateid", "farm_id"]].rename(
        columns={"WGSisolateid": "isolate_id", "farm_id": "phenotype_farm_id"}
    )
    alternate = arg.merge(matched_farms, on="isolate_id", how="left", validate="one_to_one")
    alternate["farm_id"] = alternate["phenotype_farm_id"]
    rows = []
    for label, source in [("genome_farm_id", arg), ("phenotype_farm_id", alternate)]:
        for subset in ["All farms", "Farms sampled >5 years"]:
            data = source.copy()
            if subset != "All farms":
                eligible = data.groupby("farm_id")["year"].nunique()
                data = data[data.farm_id.isin(eligible[eligible > 5].index)]
            for host in ["All", "Poultry", "Swine"]:
                part = data if host == "All" else data[data.host.eq(host)]
                row = model_row(part, host, "ARG_count", subset, False)
                row["farm_id_source"] = label
                rows.append(row)
    return pd.DataFrame(rows)


def main():
    pheno, arg, genome, wgs = read_tables()
    comparisons = comparison_table(pheno)
    models = farm_models(pheno, arg)
    checks, discrepancies = validation(pheno, arg, genome, wgs)
    farm_sensitivity = arg_farm_sensitivity(arg, wgs)
    with pd.ExcelWriter(BASE / "new_add_4" / "pre_post_drug_resistance_tests.xlsx") as writer:
        comparisons[comparisons.group.eq("Overall")].to_excel(writer, "Overall", index=False)
        comparisons[~comparisons.group.eq("Overall")].to_excel(writer, "By_host", index=False)
    with pd.ExcelWriter(BASE / "new_add_1" / "Supplementary_Table_7_farm_adjusted_sensitivity_models.xlsx") as writer:
        models.to_excel(writer, "Supplementary Table 7", index=False)
        checks.to_excel(writer, "Data_validation", index=False)
        farm_sensitivity.to_excel(writer, "ARG_farm_ID_sensitivity", index=False)
        discrepancies.to_excel(writer, "ID_discrepancies", index=False)
    print(checks.to_string(index=False))
    print("Farm models estimable:", models.estimate.notna().sum(), "of", len(models))
    print(farm_sensitivity[["farm_id_source", "subset", "host", "estimate", "p_value"]].round(4).to_string(index=False))
    print("Drug comparisons q<0.05 by group:", comparisons.assign(sig=comparisons.q_value_BH.lt(.05)).groupby("group").sig.sum().to_dict())


if __name__ == "__main__":
    main()
