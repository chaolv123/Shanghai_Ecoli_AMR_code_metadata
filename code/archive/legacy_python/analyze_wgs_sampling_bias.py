from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import chi2_contingency, spearmanr


INPUT = Path("metadata2.xlsx")
OUTDIR = Path("WGS_sampling_bias_analysis")
SEED = 20260904
N_BOOT = 1000

DRUGS = [
    "AMP", "AMC", "TIO", "CAZ", "MEM", "SPT", "GEN", "DOX",
    "TCY", "FF", "SFX", "SXT", "ENR", "OFX", "COL",
]
CORE_PROFILE_DRUGS = ["AMP", "AMC", "TIO", "GEN", "TCY", "FF", "SXT", "COL"]


def binary_smd(p1: float, p0: float) -> float:
    pooled_variance = (p1 * (1 - p1) + p0 * (1 - p0)) / 2
    denom = math.sqrt(max(float(pooled_variance), 0.0))
    return (p1 - p0) / denom if denom > 0 else 0.0


def safe_chi2(tab: pd.DataFrame) -> tuple[float, float]:
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return np.nan, np.nan
    stat, p, _, _ = chi2_contingency(tab)
    return float(stat), float(p)


def resistance_indicator(series: pd.Series) -> pd.Series:
    return series.map({"R": 1.0, "S": 0.0})


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    ok = values.notna() & weights.notna()
    if not ok.any() or weights[ok].sum() <= 0:
        return np.nan
    return float(np.average(values[ok], weights=weights[ok]))


def raking_weights(
    sample: pd.DataFrame,
    target: pd.DataFrame,
    margins: list[str],
    max_iter: int = 500,
    tol: float = 1e-8,
) -> tuple[pd.Series, int, float]:
    weights = pd.Series(len(target) / len(sample), index=sample.index, dtype=float)
    final_error = np.inf
    for iteration in range(1, max_iter + 1):
        for col in margins:
            target_counts = target[col].value_counts(dropna=False)
            current = sample.groupby(col, dropna=False).apply(
                lambda x: weights.loc[x.index].sum()
            )
            for level, target_total in target_counts.items():
                current_total = current.get(level, 0.0)
                if current_total <= 0:
                    raise ValueError(f"No WGS support for raking margin {col}={level!r}")
                if pd.isna(level):
                    idx = sample[col].isna()
                else:
                    idx = sample[col].eq(level)
                weights.loc[idx] *= target_total / current_total

        errors = []
        for col in margins:
            target_prop = target[col].value_counts(normalize=True, dropna=False)
            weighted_totals = sample.groupby(col, dropna=False).apply(
                lambda x: weights.loc[x.index].sum()
            )
            weighted_prop = weighted_totals / weighted_totals.sum()
            errors.extend(
                abs(weighted_prop.get(level, 0.0) - prop)
                for level, prop in target_prop.items()
            )
        final_error = max(errors) if errors else 0.0
        if final_error < tol:
            break
    weights *= len(target) / weights.sum()
    return weights, iteration, final_error


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    df = pd.read_excel(INPUT).rename(columns={"Year": "year"})
    required = {
        "isolate_id", "year", "period", "host", "districts", "farm_id",
        "MDR_binary", "antibiogram_length", "WGS", "WGSisolateid", *DRUGS,
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if not df["isolate_id"].is_unique:
        raise ValueError("isolate_id must be unique")
    if not set(df["WGS"].dropna().unique()).issubset({0, 1}):
        raise ValueError("WGS must be coded 0/1")

    df["resistance_class_count"] = df["antibiogram_length"]
    df["resistance_count_category"] = df["resistance_class_count"].clip(upper=9).astype(str)
    for drug in DRUGS:
        df[f"{drug}_bin"] = resistance_indicator(df[drug])
        df[f"{drug}_status"] = df[drug].fillna("Missing")

    df["profile"] = df[DRUGS].fillna("NA").astype(str).agg("|".join, axis=1)
    df["core8_profile"] = df[CORE_PROFILE_DRUGS].fillna("NA").astype(str).agg("|".join, axis=1)

    expected = df.loc[df["WGSisolateid"].notna(), "host"].value_counts().to_dict()
    observed = df.loc[df.WGS.eq(1), "host"].value_counts().to_dict()
    validation = [
        {"check": "phenotypic_rows", "observed": len(df), "expected": 5299,
         "status": "PASS" if len(df) == 5299 else "FAIL"},
        {"check": "WGS_rows_total", "observed": int(df.WGS.sum()),
         "expected": int(df["WGSisolateid"].notna().sum()),
         "status": "PASS" if int(df.WGS.sum()) == int(df["WGSisolateid"].notna().sum()) else "FAIL"},
    ]
    for host, exp in expected.items():
        obs = int(observed.get(host, 0))
        validation.append({"check": f"WGS_rows_{host}", "observed": obs,
                           "expected": exp, "status": "PASS" if obs == exp else "FAIL"})
    validation_df = pd.DataFrame(validation)

    # Descriptive selection fractions requested by the reviewer.
    selection_rows = []
    for variable in ["year", "host", "districts", "farm_id", "period"]:
        grouped = df.groupby(variable, dropna=False)
        for level, g in grouped:
            selection_rows.append({
                "variable": variable,
                "level": level,
                "n_full": len(g),
                "n_WGS": int(g.WGS.sum()),
                "selection_fraction": float(g.WGS.mean()),
                "pct_of_full": len(g) / len(df) * 100,
                "pct_of_WGS": g.WGS.sum() / df.WGS.sum() * 100,
            })
    selection_distribution = pd.DataFrame(selection_rows)

    selection_tests = []
    for variable in ["year", "host", "districts", "farm_id", "period"]:
        tab = pd.crosstab(df[variable], df.WGS)
        stat, p = safe_chi2(tab)
        selection_tests.append({
            "variable": variable, "chi_square": stat, "p_value": p,
            "n_levels": df[variable].nunique(),
        })
    selection_tests = pd.DataFrame(selection_tests)

    # Full dataset versus WGS subset; inferential comparison uses WGS vs non-WGS.
    phenotype_rows = []
    outcomes = [(d, f"{d}_bin") for d in DRUGS] + [
        ("MDR", "MDR_binary"),
    ]
    for label, col in outcomes:
        valid = df[col].notna()
        full = df.loc[valid, col]
        wgs = df.loc[valid & df.WGS.eq(1), col]
        non = df.loc[valid & df.WGS.eq(0), col]
        tab = pd.crosstab(df.loc[valid, col], df.loc[valid, "WGS"])
        stat, p = safe_chi2(tab)
        phenotype_rows.append({
            "outcome": label,
            "n_full": len(full), "n_WGS": len(wgs), "n_nonWGS": len(non),
            "prevalence_full": full.mean(), "prevalence_WGS": wgs.mean(),
            "prevalence_nonWGS": non.mean(),
            "WGS_minus_full": wgs.mean() - full.mean(),
            "SMD_WGS_vs_nonWGS": binary_smd(wgs.mean(), non.mean()),
            "chi_square_WGS_vs_nonWGS": stat, "p_value": p,
        })
    phenotype_comparison = pd.DataFrame(phenotype_rows)

    count_comparison = pd.DataFrame([
        {
            "group": group,
            "n": len(g),
            "mean_resistance_class_count": g.resistance_class_count.mean(),
            "median_resistance_class_count": g.resistance_class_count.median(),
            "IQR_low": g.resistance_class_count.quantile(0.25),
            "IQR_high": g.resistance_class_count.quantile(0.75),
            "MDR_prevalence": g.MDR_binary.mean(),
            "unique_profiles": g.profile.nunique(),
        }
        for group, g in [("Full", df), ("WGS", df[df.WGS.eq(1)]),
                         ("Non-WGS", df[df.WGS.eq(0)])]
    ])

    # Profile coverage and divergence by host and period.
    profile_rows = []
    for profile_definition, profile_col in [
        ("all_15_agents_with_missing_state", "profile"),
        ("core_8_agents", "core8_profile"),
    ]:
        for host in ["All", *sorted(df.host.unique())]:
            for period in ["All", *sorted(df.period.unique())]:
                full = df.copy()
                if host != "All":
                    full = full[full.host.eq(host)]
                if period != "All":
                    full = full[full.period.eq(period)]
                wgs = full[full.WGS.eq(1)]
                levels = sorted(set(full[profile_col]) | set(wgs[profile_col]))
                p = full[profile_col].value_counts(normalize=True).reindex(levels, fill_value=0)
                q = wgs[profile_col].value_counts(normalize=True).reindex(levels, fill_value=0)
                profile_rows.append({
                    "profile_definition": profile_definition,
                    "host": host, "period": period,
                    "n_full": len(full), "n_WGS": len(wgs),
                    "profiles_full": full[profile_col].nunique(),
                    "profiles_WGS": wgs[profile_col].nunique(),
                    "profile_richness_coverage": wgs[profile_col].nunique() / full[profile_col].nunique(),
                    "isolate_coverage_by_WGS_profiles": full[profile_col].isin(set(wgs[profile_col])).mean(),
                    "Jensen_Shannon_divergence": float(jensenshannon(p, q, base=2) ** 2),
                })
    profile_comparison = pd.DataFrame(profile_rows)

    top_profiles = (
        df.groupby(["host", "period", "profile"])
        .agg(n_full=("isolate_id", "size"), n_WGS=("WGS", "sum"))
        .reset_index()
    )
    top_profiles["pct_full"] = top_profiles.groupby(["host", "period"])["n_full"].transform(lambda x: x / x.sum())
    top_profiles["pct_WGS"] = top_profiles.groupby(["host", "period"])["n_WGS"].transform(lambda x: x / x.sum())
    top_profiles = (
        top_profiles.sort_values(["host", "period", "n_full"], ascending=[True, True, False])
        .groupby(["host", "period"]).head(20).reset_index(drop=True)
    )

    # Does the WGS subset reproduce annual phenotypic trends? Use common years.
    trend_rows = []
    annual_rows = []
    trend_outcomes = outcomes + [("Resistance_class_count", "resistance_class_count")]
    for host in sorted(df.host.unique()):
        for label, col in trend_outcomes:
            temp = df[df.host.eq(host) & df[col].notna()]
            full_annual = temp.groupby("year")[col].mean()
            wgs_annual = temp[temp.WGS.eq(1)].groupby("year")[col].mean()
            common = full_annual.index.intersection(wgs_annual.index)
            if len(common) >= 3:
                rho, p = spearmanr(full_annual.loc[common], wgs_annual.loc[common])
                mae = np.mean(abs(full_annual.loc[common] - wgs_annual.loc[common]))
            else:
                rho, p, mae = np.nan, np.nan, np.nan
            pre = temp[temp.year.between(2010, 2016)]
            post = temp[temp.year.between(2017, 2023)]
            pre_wgs = pre[pre.WGS.eq(1)]
            post_wgs = post[post.WGS.eq(1)]
            d_full = post[col].mean() - pre[col].mean()
            d_wgs = post_wgs[col].mean() - pre_wgs[col].mean()
            trend_rows.append({
                "host": host, "outcome": label, "common_years": len(common),
                "first_common_year": min(common) if len(common) else np.nan,
                "last_common_year": max(common) if len(common) else np.nan,
                "Spearman_rho": rho, "Spearman_p": p, "annual_MAE": mae,
                "prevalence_or_mean_pre_full_2010_2016": pre[col].mean(),
                "prevalence_or_mean_post_full_2017_2023": post[col].mean(),
                "change_full": d_full,
                "prevalence_or_mean_pre_WGS_2010_2016": pre_wgs[col].mean(),
                "prevalence_or_mean_post_WGS_2017_2023": post_wgs[col].mean(),
                "change_WGS": d_wgs,
                "same_change_direction": np.sign(d_full) == np.sign(d_wgs),
            })
            for year in sorted(common):
                annual_rows.append({
                    "host": host, "outcome": label, "year": year,
                    "n_full": int((temp.year.eq(year)).sum()),
                    "n_WGS": int((temp.year.eq(year) & temp.WGS.eq(1)).sum()),
                    "full_rate_or_mean": full_annual.loc[year],
                    "WGS_rate_or_mean": wgs_annual.loc[year],
                })
    trend_concordance = pd.DataFrame(trend_rows)
    annual_values = pd.DataFrame(annual_rows)

    # Weighting target is restricted to observed support: 2010-2023, non-Jiading/
    # non-Minhang districts, and farms with at least one sequenced isolate.
    supported_farms = set(df.loc[df.WGS.eq(1), "farm_id"])
    overlap = df[
        df.year.between(2010, 2023)
        & ~df.districts.isin(["Jiading", "Minhang"])
        & df.farm_id.isin(supported_farms)
    ].copy()
    wgs_overlap = overlap[overlap.WGS.eq(1)].copy()
    margins = ["year", "host", "districts", "farm_id", "MDR_binary", "resistance_count_category"]
    for drug in DRUGS:
        # Only include margins for which each target category has WGS support.
        if set(overlap[f"{drug}_status"].unique()).issubset(set(wgs_overlap[f"{drug}_status"].unique())):
            margins.append(f"{drug}_status")
    weights, iterations, rake_error = raking_weights(wgs_overlap, overlap, margins)
    wgs_overlap["calibration_weight"] = weights
    lower, upper = weights.quantile([0.01, 0.99])
    wgs_overlap["trimmed_weight"] = weights.clip(lower=lower, upper=upper)
    wgs_overlap["trimmed_weight"] *= len(overlap) / wgs_overlap.trimmed_weight.sum()

    effective_n = weights.sum() ** 2 / np.square(weights).sum()
    weight_diagnostics = pd.DataFrame([{
        "target_definition": "2010-2023; Jiading/Minhang excluded; farms with >=1 WGS isolate",
        "n_target": len(overlap), "n_WGS": len(wgs_overlap),
        "n_target_farms": overlap.farm_id.nunique(),
        "raking_iterations": iterations, "maximum_margin_error": rake_error,
        "weight_min": weights.min(), "weight_p01": lower,
        "weight_median": weights.median(), "weight_mean": weights.mean(),
        "weight_p99": upper, "weight_max": weights.max(),
        "effective_sample_size": effective_n,
        "number_of_raking_margins": len(margins),
    }])

    weighted_balance_rows = []
    for label, col in outcomes + [("Resistance_class_count", "resistance_class_count")]:
        target_mean = overlap[col].mean()
        raw_mean = wgs_overlap[col].mean()
        weighted = weighted_mean(wgs_overlap[col], wgs_overlap.calibration_weight)
        trimmed = weighted_mean(wgs_overlap[col], wgs_overlap.trimmed_weight)
        weighted_balance_rows.append({
            "outcome": label, "target_full_overlap": target_mean,
            "WGS_unweighted": raw_mean, "WGS_calibrated": weighted,
            "WGS_trimmed": trimmed,
            "absolute_bias_unweighted": abs(raw_mean - target_mean),
            "absolute_bias_calibrated": abs(weighted - target_mean),
            "absolute_bias_trimmed": abs(trimmed - target_mean),
        })
    weighted_balance = pd.DataFrame(weighted_balance_rows)

    # Weighted estimates and bootstrap directional stability for period contrasts.
    rng = np.random.default_rng(SEED)
    weighted_trend_rows = []
    bootstrap_rows = []
    for host in sorted(overlap.host.unique()):
        target_h = overlap[overlap.host.eq(host)]
        sample_h = wgs_overlap[wgs_overlap.host.eq(host)]
        probs = sample_h.calibration_weight.to_numpy()
        probs = probs / probs.sum()
        for label, col in trend_outcomes:
            target_pre = target_h[target_h.year.between(2010, 2016)][col].mean()
            target_post = target_h[target_h.year.between(2017, 2023)][col].mean()
            raw_pre = sample_h[sample_h.year.between(2010, 2016)][col].mean()
            raw_post = sample_h[sample_h.year.between(2017, 2023)][col].mean()
            cal_pre = weighted_mean(
                sample_h.loc[sample_h.year.between(2010, 2016), col],
                sample_h.loc[sample_h.year.between(2010, 2016), "calibration_weight"],
            )
            cal_post = weighted_mean(
                sample_h.loc[sample_h.year.between(2017, 2023), col],
                sample_h.loc[sample_h.year.between(2017, 2023), "calibration_weight"],
            )
            weighted_trend_rows.append({
                "host": host, "outcome": label,
                "target_change": target_post - target_pre,
                "WGS_raw_change": raw_post - raw_pre,
                "WGS_calibrated_change": cal_post - cal_pre,
                "raw_same_direction": np.sign(raw_post - raw_pre) == np.sign(target_post - target_pre),
                "calibrated_same_direction": np.sign(cal_post - cal_pre) == np.sign(target_post - target_pre),
            })

            boot_changes = []
            for _ in range(N_BOOT):
                sampled_idx = rng.choice(sample_h.index.to_numpy(), size=len(sample_h), replace=True, p=probs)
                b = sample_h.loc[sampled_idx]
                pre_b = b[b.year.between(2010, 2016)][col].mean()
                post_b = b[b.year.between(2017, 2023)][col].mean()
                boot_changes.append(post_b - pre_b)
            boot_changes = np.asarray(boot_changes, dtype=float)
            target_direction = np.sign(target_post - target_pre)
            bootstrap_rows.append({
                "host": host, "outcome": label,
                "target_change": target_post - target_pre,
                "bootstrap_mean_change": np.nanmean(boot_changes),
                "bootstrap_CI_low": np.nanquantile(boot_changes, 0.025),
                "bootstrap_CI_high": np.nanquantile(boot_changes, 0.975),
                "proportion_same_direction_as_target": np.nanmean(np.sign(boot_changes) == target_direction),
                "n_bootstrap": N_BOOT,
            })
    weighted_trends = pd.DataFrame(weighted_trend_rows)
    bootstrap_stability = pd.DataFrame(bootstrap_rows)

    weight_export_cols = [
        "isolate_id", "year", "period", "host", "districts", "farm_id",
        "WGS", "MDR_binary", "resistance_class_count", "calibration_weight",
        "trimmed_weight",
    ]
    wgs_overlap[weight_export_cols].to_csv(
        OUTDIR / "WGS_calibration_weights.csv", index=False, encoding="utf-8-sig"
    )

    workbook = OUTDIR / "Supplementary_WGS_sampling_bias_analysis.xlsx"
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        validation_df.to_excel(writer, "Data_validation", index=False)
        selection_distribution.to_excel(writer, "Selection_distribution", index=False)
        selection_tests.to_excel(writer, "Selection_tests", index=False)
        phenotype_comparison.to_excel(writer, "Phenotype_comparison", index=False)
        count_comparison.to_excel(writer, "Burden_profile_summary", index=False)
        profile_comparison.to_excel(writer, "Profile_comparison", index=False)
        top_profiles.to_excel(writer, "Top_profiles", index=False)
        trend_concordance.to_excel(writer, "Trend_concordance", index=False)
        annual_values.to_excel(writer, "Annual_values", index=False)
        weight_diagnostics.to_excel(writer, "Weight_diagnostics", index=False)
        weighted_balance.to_excel(writer, "Weighted_balance", index=False)
        weighted_trends.to_excel(writer, "Weighted_period_trends", index=False)
        bootstrap_stability.to_excel(writer, "Bootstrap_stability", index=False)

    significant_bias = phenotype_comparison.loc[
        phenotype_comparison.p_value.lt(0.05), "outcome"
    ].tolist()
    direction_rate = trend_concordance.same_change_direction.mean()
    high_corr = trend_concordance.Spearman_rho.ge(0.7).mean()
    report = f"""WGS sampling-bias analysis summary
==================================
Input: {INPUT}
Phenotypic isolates: {len(df)}
WGS-coded isolates: {int(df.WGS.sum())}
Matched to phenotypic metadata: {int(df['WGSisolateid'].notna().sum())} (Poultry {expected.get('Poultry', 0)}; Swine {expected.get('Swine', 0)})
Observed: Poultry {observed.get('Poultry', 0)}; Swine {observed.get('Swine', 0)}

Derived variables
-----------------
MDR_binary was analysed as the supplied 0/1 indicator. Antibiogram_length (0-10) was
also copied to resistance_class_count for clear labelling in output tables.

Main descriptive findings
-------------------------
Overall WGS selection fraction: {df.WGS.mean():.3%}.
No WGS-coded isolates were present in 2008-2009.
Jiading had {int(df.loc[df.districts.eq('Jiading'), 'WGS'].sum())} WGS-coded isolates;
Minhang had {int(df.loc[df.districts.eq('Minhang'), 'WGS'].sum())}.
Farms represented in WGS: {df.groupby('farm_id').WGS.sum().gt(0).sum()} of {df.farm_id.nunique()}.
Phenotypes with WGS-vs-non-WGS chi-square P<0.05: {', '.join(significant_bias)}.

Phenotypic trend reproduction
-----------------------------
Across host-outcome combinations, {direction_rate:.1%} had the same 2010-2016 versus
2017-2023 direction in the full dataset and WGS subset. {high_corr:.1%} had annual
Spearman rho >= 0.70. See Trend_concordance for drug-specific results.

Weighting/resampling
--------------------
Because there is no WGS support in 2008-2009 and several districts/farms, weights cannot
identify genomic outcomes for all 5299 isolates. Calibration was therefore restricted to
the overlap population: 2010-2023, excluding Jiading/Minhang, and farms with >=1 WGS isolate.
Target n={len(overlap)}; WGS n={len(wgs_overlap)}; effective weighted n={effective_n:.1f}.
Raking stopped after {iterations} iterations (maximum margin error={rake_error:.3g}; target tolerance=1e-8).

Unresolved requirement
----------------------
ARG_count, plasmid-associated ARG burden/count, PT, ST and other genome-derived fields are
not present in metadata2.xlsx. The exported isolate-level calibration weights can be merged
with genomic data using WGSisolateid (phenotype table) and isolate_id (genome table); 75
additional poultry genomes from 2008-2009 have no phenotype-table match and are outside
the weighting target. Weighted ARG-burden results are in the companion workbook
Supplementary_weighted_ARG_burden_updated.xlsx; plasmid burden and PT were unavailable.
"""
    (OUTDIR / "analysis_summary.txt").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
