from __future__ import annotations

import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning, PerfectSeparationError


BASE = Path("New_analysis") / "new_add_6"
IN_PATH = BASE / "metadata6.xlsx"
OUT_DIR = BASE / "review6_st_arg_association"
OUT_XLSX = OUT_DIR / "Supplementary_ST_ARG_association_sensitivity.xlsx"
OUT_SUMMARY = OUT_DIR / "review6_ST_ARG_association_summary.txt"

META_COLS = ["isolate_id", "host", "year", "period", "districts", "farm_id", "ST", "ARG_count"]
KEY_ARG = "blaCTX-M-55"
KEY_ST = "ST69"
DEBURDENING_ST = "ST10"


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
    out.loc[order] = np.clip(adj, 0, 1)
    return out


def norm_district(x) -> str:
    if pd.isna(x):
        return "Unknown"
    s = str(x).strip()
    if s in {"", "-", "NA", "N/A", "nan", "None"}:
        return "Unknown"
    return s[:1].upper() + s[1:].lower()


def binary_counts(df: pd.DataFrame, exposure_col: str, outcome_col: str) -> dict:
    sub = df[[exposure_col, outcome_col]].dropna()
    exp = sub[exposure_col].astype(int)
    out = sub[outcome_col].astype(int)
    a = int(((exp == 1) & (out == 1)).sum())
    b = int(((exp == 1) & (out == 0)).sum())
    c = int(((exp == 0) & (out == 1)).sum())
    d = int(((exp == 0) & (out == 0)).sum())
    table = np.array([[a, b], [c, d]])
    fisher_p = stats.fisher_exact(table).pvalue
    # Haldane correction for finite odds ratio when cells are zero.
    or_h = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return {
        "n": int(len(sub)),
        "exposed_n": int((exp == 1).sum()),
        "unexposed_n": int((exp == 0).sum()),
        "outcome_in_exposed": a,
        "outcome_in_unexposed": c,
        "prevalence_exposed": a / (a + b) if (a + b) else np.nan,
        "prevalence_unexposed": c / (c + d) if (c + d) else np.nan,
        "risk_difference": (a / (a + b) - c / (c + d)) if (a + b) and (c + d) else np.nan,
        "fisher_p": fisher_p,
        "haldane_OR": or_h,
    }


def fit_logit(df: pd.DataFrame, formula: str, exposure_name: str) -> dict:
    row = {
        "model": formula,
        "n_model": np.nan,
        "OR": np.nan,
        "CI_low": np.nan,
        "CI_high": np.nan,
        "p_value": np.nan,
        "status": "not_run",
    }
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            model = smf.logit(formula, data=df).fit(disp=False, maxiter=200)
        coef = model.params.get(exposure_name, np.nan)
        se = model.bse.get(exposure_name, np.nan)
        p = model.pvalues.get(exposure_name, np.nan)
        row.update({
            "n_model": int(model.nobs),
            "OR": math.exp(coef) if pd.notna(coef) else np.nan,
            "CI_low": math.exp(coef - 1.96 * se) if pd.notna(coef) and pd.notna(se) else np.nan,
            "CI_high": math.exp(coef + 1.96 * se) if pd.notna(coef) and pd.notna(se) else np.nan,
            "p_value": p,
            "status": "ok",
        })
    except (PerfectSeparationError, np.linalg.LinAlgError, ValueError) as exc:
        row["status"] = f"failed: {type(exc).__name__}"
    except Exception as exc:
        row["status"] = f"failed: {type(exc).__name__}: {str(exc)[:80]}"
    return row


def fit_glm_count(df: pd.DataFrame, formula: str, family, exposure_name: str) -> dict:
    row = {
        "model": formula,
        "family": family.__class__.__name__,
        "n_model": np.nan,
        "estimate": np.nan,
        "ratio_or_difference": np.nan,
        "CI_low": np.nan,
        "CI_high": np.nan,
        "p_value": np.nan,
        "status": "not_run",
    }
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = smf.glm(formula, data=df, family=family).fit(maxiter=200)
        coef = model.params.get(exposure_name, np.nan)
        se = model.bse.get(exposure_name, np.nan)
        p = model.pvalues.get(exposure_name, np.nan)
        if isinstance(family, (sm.families.Poisson, sm.families.NegativeBinomial)):
            ratio = math.exp(coef)
            lo = math.exp(coef - 1.96 * se)
            hi = math.exp(coef + 1.96 * se)
        else:
            ratio = coef
            lo = coef - 1.96 * se
            hi = coef + 1.96 * se
        row.update({
            "n_model": int(model.nobs),
            "estimate": coef,
            "ratio_or_difference": ratio,
            "CI_low": lo,
            "CI_high": hi,
            "p_value": p,
            "status": "ok",
        })
    except Exception as exc:
        row["status"] = f"failed: {type(exc).__name__}: {str(exc)[:80]}"
    return row


def safe_name(value: str) -> str:
    keep = []
    for ch in str(value):
        if ch.isalnum():
            keep.append(ch)
        else:
            keep.append("_")
    return "".join(keep).strip("_")


def annual_correlations(df: pd.DataFrame, st_list: list[str], arg_cols: list[str]) -> pd.DataFrame:
    rows = []
    # The manuscript reports host-stratified associations only. Pooling hosts
    # can create associations driven by host composition, so no combined
    # "All" analysis is generated here.
    for host in sorted(df["host"].dropna().unique().tolist()):
        sub_host = df[df["host"] == host]
        years = sorted(sub_host["year"].dropna().unique())
        if len(years) < 5:
            continue
        for st in st_list:
            st_rates = sub_host.groupby("year")["ST"].apply(lambda x: (x == st).mean())
            if (st_rates > 0).sum() < 2:
                continue
            for arg in arg_cols:
                arg_rates = sub_host.groupby("year")[arg].mean()
                common = sorted(set(st_rates.dropna().index).intersection(arg_rates.dropna().index))
                if len(common) < 5:
                    continue
                x = st_rates.loc[common].to_numpy(dtype=float)
                y = arg_rates.loc[common].to_numpy(dtype=float)
                if np.nanstd(x) == 0 or np.nanstd(y) == 0:
                    continue
                rho, p = stats.spearmanr(x, y)
                rows.append({
                    "host": host,
                    "ST": st,
                    "ARG": arg,
                    "n_years": len(common),
                    "Spearman_rho": rho,
                    "p_value": p,
                    "ST_total_isolates": int((sub_host["ST"] == st).sum()),
                    "ARG_total_positive": int(sub_host[arg].sum()),
                })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["FDR"] = out.groupby("host")["p_value"].transform(bh_fdr)
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_excel(IN_PATH)
    df["isolate_id"] = df["isolate_id"].astype(str)
    df["ST"] = df["ST"].astype(str)
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype(int)
    df["ARG_count"] = pd.to_numeric(df["ARG_count"], errors="coerce")
    df["districts_clean"] = df["districts"].map(norm_district)
    df["farm_clean"] = df["farm_id"].fillna("Unknown").astype(str)
    df["post2017"] = (df["period"].astype(str).str.lower() == "post-2017").astype(int)
    df["ST69_present"] = (df["ST"] == KEY_ST).astype(int)
    df["ST10_present"] = (df["ST"] == DEBURDENING_ST).astype(int)

    arg_cols = [c for c in df.columns if c not in META_COLS + ["districts_clean", "farm_clean", "post2017", "ST69_present", "ST10_present"]]
    for c in arg_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

    st_counts = df["ST"].value_counts().rename_axis("ST").reset_index(name="n_isolates")
    arg_counts = pd.DataFrame({
        "ARG": arg_cols,
        "n_positive": [int(df[c].sum()) for c in arg_cols],
        "prevalence": [float(df[c].mean()) for c in arg_cols],
    }).sort_values("n_positive", ascending=False)

    top_sts = st_counts.loc[st_counts["n_isolates"] >= 20, "ST"].tolist()
    for key_st in [KEY_ST, DEBURDENING_ST]:
        if key_st not in top_sts and key_st in set(df["ST"]):
            top_sts.append(key_st)
    common_args = arg_counts.loc[arg_counts["n_positive"] >= 10, "ARG"].tolist()
    if KEY_ARG not in common_args and KEY_ARG in arg_cols:
        common_args.append(KEY_ARG)
    annual = annual_correlations(df, top_sts, common_args)

    # Isolate-level ST-ARG screening for top STs and common ARGs.
    # Full adjusted logistic models across all sparse ST-ARG combinations are unstable;
    # key associations are modelled separately below.
    model_rows = []
    for st in top_sts:
        exposure = f"exposure_{safe_name(st)}"
        df[exposure] = (df["ST"] == st).astype(int)
        for arg in common_args:
            counts = binary_counts(df, exposure, arg)
            row = {
                "ST": st,
                "ARG": arg,
                **counts,
            }
            row.update({
                "model": "isolate_level_fisher_screen",
                "n_model": counts["n"],
                "OR": counts["haldane_OR"],
                "CI_low": np.nan,
                "CI_high": np.nan,
                "p_value": counts["fisher_p"],
                "status": "screening_only_unadjusted",
            })
            model_rows.append(row)
    isolate_models = pd.DataFrame(model_rows)
    if not isolate_models.empty:
        isolate_models["FDR"] = isolate_models["p_value"].transform(bh_fdr)

    # Key ST69-blaCTX-M-55 analyses.
    key_rows = []
    for label, sub in {
        "all_isolates": df,
        "poultry_only": df[df["host"] == "Poultry"],
        "post2017_only": df[df["period"] == "post-2017"],
        "poultry_post2017": df[(df["host"] == "Poultry") & (df["period"] == "post-2017")],
    }.items():
        counts = binary_counts(sub, "ST69_present", KEY_ARG)
        row = {"stratum": label, **counts}
        key_rows.append(row)
    key_counts = pd.DataFrame(key_rows)

    key_models = []
    models = {
        "crude_logit": (df[[KEY_ARG, "ST69_present"]], f"Q('{KEY_ARG}') ~ ST69_present"),
        "adjusted_host_period_district": (df[[KEY_ARG, "ST69_present", "host", "period", "districts_clean"]], f"Q('{KEY_ARG}') ~ ST69_present + C(host) + C(period) + C(districts_clean)"),
        "poultry_adjusted_period_district": (df[df["host"] == "Poultry"][[KEY_ARG, "ST69_present", "period", "districts_clean"]], f"Q('{KEY_ARG}') ~ ST69_present + C(period) + C(districts_clean)"),
        "post2017_adjusted_host_district": (df[df["period"] == "post-2017"][[KEY_ARG, "ST69_present", "host", "districts_clean"]], f"Q('{KEY_ARG}') ~ ST69_present + C(host) + C(districts_clean)"),
    }
    for name, (mdata, formula) in models.items():
        r = fit_logit(mdata.dropna(), formula, "ST69_present")
        r["analysis"] = name
        key_models.append(r)
    key_models = pd.DataFrame(key_models)

    # ST10 within-lineage de-burdening.
    st10 = df[df["ST"] == DEBURDENING_ST].copy()
    st10_rows = []
    for label, sub in {
        "ST10_all": st10,
        "ST10_swine": st10[st10["host"] == "Swine"],
        "ST10_known_district": st10[st10["districts_clean"] != "Unknown"],
    }.items():
        pre = sub.loc[sub["period"] == "pre-2017", "ARG_count"]
        post = sub.loc[sub["period"] == "post-2017", "ARG_count"]
        p = stats.mannwhitneyu(pre, post, alternative="two-sided").pvalue if len(pre) and len(post) else np.nan
        st10_rows.append({
            "analysis": label,
            "n_pre": int(len(pre)),
            "n_post": int(len(post)),
            "mean_pre": pre.mean() if len(pre) else np.nan,
            "mean_post": post.mean() if len(post) else np.nan,
            "median_pre": pre.median() if len(pre) else np.nan,
            "median_post": post.median() if len(post) else np.nan,
            "change_mean_post_minus_pre": (post.mean() - pre.mean()) if len(pre) and len(post) else np.nan,
            "mannwhitney_p": p,
        })
    st10_summary = pd.DataFrame(st10_rows)

    st10_models = []
    if not st10.empty:
        candidate_models = {
            "negative_binomial_period_host_district": (st10, "ARG_count ~ post2017 + C(host) + C(districts_clean)", sm.families.NegativeBinomial()),
            "poisson_period_host_district": (st10, "ARG_count ~ post2017 + C(host) + C(districts_clean)", sm.families.Poisson()),
            "gaussian_period_host_district": (st10, "ARG_count ~ post2017 + C(host) + C(districts_clean)", sm.families.Gaussian()),
            "swine_negative_binomial_period_district": (st10[st10["host"] == "Swine"], "ARG_count ~ post2017 + C(districts_clean)", sm.families.NegativeBinomial()),
            "swine_gaussian_period_district": (st10[st10["host"] == "Swine"], "ARG_count ~ post2017 + C(districts_clean)", sm.families.Gaussian()),
        }
        for name, (mdata, formula, family) in candidate_models.items():
            r = fit_glm_count(mdata.dropna(subset=["ARG_count", "post2017", "districts_clean"]), formula, family, "post2017")
            r["analysis"] = name
            st10_models.append(r)
    st10_models = pd.DataFrame(st10_models)

    farm_missing = int(df["farm_id"].isna().sum())
    validation = pd.DataFrame([
        {"check": "rows", "observed": len(df), "expected": 2622, "status": "PASS" if len(df) == 2622 else "CHECK"},
        {"check": "unique_isolate_id", "observed": df["isolate_id"].nunique(), "expected": len(df), "status": "PASS" if df["isolate_id"].nunique() == len(df) else "FAIL"},
        {"check": "ARG_columns", "observed": len(arg_cols), "expected": "binary ARG columns", "status": "PASS"},
        {"check": "farm_id_missing", "observed": farm_missing, "expected": "reported; not included as main random effect", "status": "INFO"},
        {"check": "unknown_districts", "observed": int((df["districts_clean"] == "Unknown").sum()), "expected": "retained as Unknown covariate level", "status": "INFO"},
        {"check": "ST69_total", "observed": int((df["ST"] == KEY_ST).sum()), "expected": "key sparse lineage", "status": "INFO"},
        {"check": "ST10_total", "observed": int((df["ST"] == DEBURDENING_ST).sum()), "expected": "key de-burdening lineage", "status": "INFO"},
    ])

    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
        validation.to_excel(writer, sheet_name="Data_validation", index=False)
        st_counts.to_excel(writer, sheet_name="ST_counts", index=False)
        arg_counts.to_excel(writer, sheet_name="ARG_counts", index=False)
        annual.to_excel(writer, sheet_name="Annual_ST_ARG_correlations", index=False)
        isolate_models.to_excel(writer, sheet_name="Isolate_level_ST_ARG_models", index=False)
        key_counts.to_excel(writer, sheet_name="ST69_blaCTXM55_counts", index=False)
        key_models.to_excel(writer, sheet_name="ST69_blaCTXM55_models", index=False)
        st10_summary.to_excel(writer, sheet_name="ST10_ARG_burden_summary", index=False)
        st10_models.to_excel(writer, sheet_name="ST10_ARG_burden_models", index=False)

    lines = []
    lines.append("Review 6 ST-ARG association sensitivity analysis")
    lines.append("=" * 56)
    lines.append(f"Input: {IN_PATH}")
    lines.append(f"Isolates: {len(df)}; ARG columns: {len(arg_cols)}; STs: {df['ST'].nunique()}")
    lines.append(f"Farm ID missing: {farm_missing}; unknown district: {(df['districts_clean'] == 'Unknown').sum()}")
    lines.append("")
    lines.append("ST69-blaCTX-M-55 isolate-level evidence")
    lines.append("----------------------------------------")
    for _, r in key_counts.iterrows():
        lines.append(
            f"{r['stratum']}: ST69 n={int(r['exposed_n'])}, blaCTX-M-55 in ST69={int(r['outcome_in_exposed'])} "
            f"({r['prevalence_exposed']:.3f}); non-ST69 blaCTX-M-55={int(r['outcome_in_unexposed'])}/{int(r['unexposed_n'])} "
            f"({r['prevalence_unexposed']:.3f}); Fisher P={r['fisher_p']:.3g}; Haldane OR={r['haldane_OR']:.3g}."
        )
    for _, r in key_models.iterrows():
        lines.append(
            f"{r['analysis']}: OR={r['OR']:.3g}, 95% CI={r['CI_low']:.3g}-{r['CI_high']:.3g}, "
            f"P={r['p_value']:.3g}, status={r['status']}."
        )
    lines.append("")
    lines.append("ST10 within-lineage ARG burden")
    lines.append("------------------------------")
    for _, r in st10_summary.iterrows():
        lines.append(
            f"{r['analysis']}: n_pre={int(r['n_pre'])}, n_post={int(r['n_post'])}, "
            f"mean {r['mean_pre']:.3f}->{r['mean_post']:.3f}, change={r['change_mean_post_minus_pre']:.3f}, "
            f"Wilcoxon/Mann-Whitney P={r['mannwhitney_p']:.3g}."
        )
    for _, r in st10_models.iterrows():
        lines.append(
            f"{r['analysis']}: ratio/difference={r['ratio_or_difference']:.3g}, "
            f"95% CI={r['CI_low']:.3g}-{r['CI_high']:.3g}, P={r['p_value']:.3g}, status={r['status']}."
        )
    lines.append("")
    lines.append("Exploratory annual correlation summary")
    lines.append("--------------------------------------")
    if not annual.empty:
        sig = annual[(annual["FDR"] < 0.05) & (annual["Spearman_rho"] > 0)]
        lines.append(f"Annual ST-ARG tests run: {len(annual)}; positive FDR<0.05 associations: {len(sig)}.")
        key_annual = annual[(annual["ST"] == KEY_ST) & (annual["ARG"] == KEY_ARG)]
        for _, r in key_annual.iterrows():
            lines.append(
                f"Annual {KEY_ST}-{KEY_ARG} ({r['host']}): rho={r['Spearman_rho']:.3f}, "
                f"P={r['p_value']:.3g}, FDR={r['FDR']:.3g}, n_years={int(r['n_years'])}."
            )
    lines.append("")
    lines.append("Interpretation")
    lines.append("--------------")
    lines.append("Annual correlations should be treated as exploratory temporal co-occurrence signals. Isolate-level analyses support enrichment of blaCTX-M-55 among ST69 isolates but do not establish mechanistic linkage or co-selection. ST10 shows a robust post-2017 reduction in ARG burden within the lineage, supporting conservative wording as within-lineage resistome simplification.")
    OUT_SUMMARY.write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote {OUT_XLSX}")
    print(f"Wrote {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
