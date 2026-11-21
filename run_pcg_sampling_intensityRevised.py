from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint


ANALYSIS_ROOT = Path(__file__).resolve().parent
PCG_XLSX = (
    ANALYSIS_ROOT
    / "outputs"
    / "extended"
    / "pcg_reconciled"
    / "Supplementary_PCG_threshold_reconciled.xlsx"
)
PHENO_XLSX = ANALYSIS_ROOT / "data" / "raw" / "phenotype_metadata.xlsx"
WEIGHTS_CSV = ANALYSIS_ROOT / "outputs" / "extended" / "wgs" / "WGS_calibration_weights.csv"
OUTDIR = ANALYSIS_ROOT / "outputs" / "extended" / "pcg_sampling_intensity"
OUT_XLSX = OUTDIR / "Supplementary_PCG_sampling_intensity_normalized_analysis.xlsx"


def strip_suffix(value: object) -> str | None:
    if pd.isna(value):
        return None
    x = str(value).strip()
    x = x.split(".R1")[0].split("_HQ")[0]
    x = x.replace("-8.filtered", "")
    if x.endswith("-8"):
        x = x[:-2]
    return x


def normalize_id(value: object) -> str | None:
    x = strip_suffix(value)
    return None if x is None else x.replace("_", "")


def wilson(x: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    return tuple(proportion_confint(x, n, alpha=0.05, method="wilson"))


def weighted_prop(d: pd.DataFrame, outcome: str, weight: str) -> tuple[float, float, float]:
    w = d[weight].to_numpy(float)
    y = d[outcome].to_numpy(float)
    p = float(np.sum(w * y) / np.sum(w))
    neff = float(np.sum(w) ** 2 / np.sum(w**2))
    se = math.sqrt(max(p * (1 - p) / neff, 0))
    return p, neff, se


def farm_cluster_bootstrap(
    d: pd.DataFrame,
    outcome: str,
    weight: str,
    reps: int = 2000,
    seed: int = 20261116,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for host, h in d.groupby("host", observed=True):
        farms = np.asarray(h["farm_id"].dropna().unique())
        # Pre-aggregate by farm and period so that cluster resampling is fast.
        g = (
            h.assign(weighted_event=h[weight] * h[outcome])
            .groupby(["farm_id", "period"], observed=True)
            .agg(weight_sum=(weight, "sum"), weighted_event_sum=("weighted_event", "sum"))
            .reset_index()
        )
        periods = ["pre-2017", "post-2017"]
        wmat = np.zeros((len(farms), 2), dtype=float)
        emat = np.zeros((len(farms), 2), dtype=float)
        farm_index = {f: i for i, f in enumerate(farms)}
        period_index = {p: i for i, p in enumerate(periods)}
        for row in g.itertuples(index=False):
            if row.period in period_index:
                i, j = farm_index[row.farm_id], period_index[row.period]
                wmat[i, j] = row.weight_sum
                emat[i, j] = row.weighted_event_sum
        vals = []
        for _ in range(reps):
            idx = rng.integers(0, len(farms), size=len(farms))
            ws = wmat[idx].sum(axis=0)
            es = emat[idx].sum(axis=0)
            if np.all(ws > 0):
                ps = es / ws
                vals.append(ps[1] - ps[0])
        vals = np.asarray(vals, float)
        rows.append(
            {
                "host": host,
                "outcome": outcome,
                "weight": weight,
                "bootstrap_reps_requested": reps,
                "bootstrap_reps_valid": len(vals),
                "post_minus_pre": float(np.mean(vals)),
                "ci_low": float(np.quantile(vals, 0.025)),
                "ci_high": float(np.quantile(vals, 0.975)),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    assn = pd.read_excel(PCG_XLSX, sheet_name="Cluster_assignments")
    assn = assn.loc[assn["threshold"].eq(10)].copy()
    # The compact public workbook records a cluster identifier for every
    # isolate. Only clusters with at least two isolates are PCGs.
    assn["PCG_id"] = pd.Series(pd.NA, index=assn.index, dtype="object")
    in_pcg = assn["cluster_size"].ge(2)
    assn.loc[in_pcg, "PCG_id"] = (
        "PCG10_Ward_"
        + assn.loc[in_pcg, "cluster_id"].astype(int).astype(str).str.zfill(4)
    )
    assn["localized_same_year_host_district"] = assn["PCG_type"].eq("localized")
    assn["nonlocalized_or_unknown"] = assn["PCG_type"].isin(
        ["nonlocalized", "non-localized"]
    )
    pcs = (
        assn.loc[assn["PCG_id"].notna(), [
            "PCG_id", "localized_same_year_host_district", "nonlocalized_or_unknown"
        ]]
        .drop_duplicates("PCG_id")
        .copy()
    )
    assn["localized_pcg"] = assn["localized_same_year_host_district"].fillna(False).astype(int)
    assn["any_pcg"] = assn["PCG_id"].notna().astype(int)
    assn["nonlocalized_pcg"] = (
        assn["nonlocalized_or_unknown"].fillna(False) & assn["PCG_id"].notna()
    ).astype(int)

    # Period-specific denominators and normalized outcomes for all PCG-analyzable isolates.
    denom_rows = []
    for (host, period), d in assn.groupby(["host", "period"], observed=True):
        n = len(d)
        n_pcg = int(d["any_pcg"].sum())
        n_local = int(d["localized_pcg"].sum())
        n_nonlocal = int(d["nonlocalized_pcg"].sum())
        local_clusters = int(d.loc[d["localized_pcg"].eq(1), "PCG_id"].nunique())
        lo, hi = wilson(n_local, n)
        denom_rows.append(
            {
                "host": host,
                "period": period,
                "n_PCG_analyzable_WGS_isolates": n,
                "n_isolates_in_any_PCG": n_pcg,
                "proportion_in_any_PCG": n_pcg / n,
                "n_isolates_in_localized_PCG": n_local,
                "proportion_in_localized_PCG": n_local / n,
                "localized_prop_CI_low": lo,
                "localized_prop_CI_high": hi,
                "n_isolates_in_nonlocalized_PCG": n_nonlocal,
                "n_distinct_localized_PCGs": local_clusters,
                "localized_PCGs_per_100_WGS_isolates": local_clusters / n * 100,
            }
        )
    denominators = pd.DataFrame(denom_rows)

    test_rows = []
    for host, d in assn.groupby("host", observed=True):
        for outcome in ["localized_pcg", "any_pcg"]:
            tab = pd.crosstab(d["period"], d[outcome]).reindex(
                index=["pre-2017", "post-2017"], columns=[0, 1], fill_value=0
            )
            chi2, p, _, _ = chi2_contingency(tab.to_numpy(), correction=False)
            pre_event, pre_n = int(tab.loc["pre-2017", 1]), int(tab.loc["pre-2017"].sum())
            post_event, post_n = int(tab.loc["post-2017", 1]), int(tab.loc["post-2017"].sum())
            # Haldane-Anscombe correction is only used if a cell is zero.
            a, b = post_event, post_n - post_event
            c, e = pre_event, pre_n - pre_event
            if min(a, b, c, e) == 0:
                a, b, c, e = [x + 0.5 for x in (a, b, c, e)]
            odds_ratio = (a / b) / (c / e)
            se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / e)
            test_rows.append(
                {
                    "host": host,
                    "outcome": outcome,
                    "n_pre": pre_n,
                    "events_pre": pre_event,
                    "proportion_pre": pre_event / pre_n,
                    "n_post": post_n,
                    "events_post": post_event,
                    "proportion_post": post_event / post_n,
                    "post_minus_pre": post_event / post_n - pre_event / pre_n,
                    "unadjusted_OR_post_vs_pre": odds_ratio,
                    "OR_CI_low": math.exp(math.log(odds_ratio) - 1.96 * se),
                    "OR_CI_high": math.exp(math.log(odds_ratio) + 1.96 * se),
                    "chi_square": chi2,
                    "p_value": p,
                }
            )
    tests = pd.DataFrame(test_rows)
    tests["p_FDR_BH_within_outcome"] = np.nan
    for outcome, idx in tests.groupby("outcome").groups.items():
        tests.loc[idx, "p_FDR_BH_within_outcome"] = multipletests(
            tests.loc[idx, "p_value"], method="fdr_bh"
        )[1]

    # Link PCG assignments to phenotype-linked WGS isolates and existing calibration weights.
    pheno = pd.read_excel(PHENO_XLSX)
    if "Year" not in pheno.columns and "year" in pheno.columns:
        pheno = pheno.rename(columns={"year": "Year"})
    selected = pheno.loc[
        pheno["WGS"].eq(1),
        ["isolate_id", "WGSisolateid", "Year", "period", "host", "districts", "farm_id"],
    ].copy()
    selected = selected.rename(columns={"isolate_id": "phenotype_isolate_id", "Year": "year_pheno"})
    weights = pd.read_csv(WEIGHTS_CSV).rename(columns={"isolate_id": "phenotype_isolate_id"})
    selected = selected.merge(
        weights[["phenotype_isolate_id", "calibration_weight", "trimmed_weight"]],
        on="phenotype_isolate_id",
        how="left",
        validate="1:1",
    )
    selected["match_exact"] = selected["WGSisolateid"].map(strip_suffix)
    selected["match_norm"] = selected["WGSisolateid"].map(normalize_id)

    aa = assn.copy()
    aa["match_exact"] = aa["isolate_id"].map(strip_suffix)
    aa["match_norm"] = aa["isolate_id"].map(normalize_id)
    # Prefer direct ID matches; only use underscore-insensitive matching when unique.
    exact = aa.merge(selected, on="match_exact", how="inner", suffixes=("_pcg", "_pheno"))
    exact_keys = set(exact["match_norm_pcg"] if "match_norm_pcg" in exact else exact["match_norm"])
    aa_left = aa.loc[~aa["match_norm"].isin(exact_keys)].copy()
    sel_left = selected.loc[~selected["match_norm"].isin(exact_keys)].copy()
    aa_unique = aa_left.loc[~aa_left["match_norm"].duplicated(keep=False)]
    fallback = aa_unique.merge(sel_left, on="match_norm", how="inner", suffixes=("_pcg", "_pheno"))
    linked = pd.concat([exact, fallback], ignore_index=True, sort=False)

    # Harmonize columns after the two match routes.
    def pick(frame: pd.DataFrame, candidates: list[str]) -> pd.Series:
        for c in candidates:
            if c in frame.columns:
                return frame[c]
        raise KeyError(candidates)

    model = pd.DataFrame(
        {
            "pcg_isolate_id": pick(linked, ["isolate_id", "isolate_id_pcg"]),
            "phenotype_isolate_id": linked["phenotype_isolate_id"],
            "year": pick(linked, ["year", "year_pcg"]),
            "period": pick(linked, ["period_pcg", "period"]),
            "host": pick(linked, ["host_pcg", "host"]),
            "districts": pick(linked, ["districts_pheno", "districts_y", "districts"]),
            "farm_id": linked["farm_id"],
            "PCG_id": linked["PCG_id"],
            "localized_pcg": linked["localized_pcg"].astype(int),
            "any_pcg": linked["any_pcg"].astype(int),
            "calibration_weight": linked["calibration_weight"],
            "trimmed_weight": linked["trimmed_weight"],
        }
    )
    model = model.drop_duplicates("phenotype_isolate_id", keep="first").copy()
    model["districts"] = model["districts"].fillna("Unknown").astype(str)
    model["farm_id"] = model["farm_id"].fillna("Unknown").astype(str)
    model["post"] = model["period"].eq("post-2017").astype(int)
    model["wgs_stratum_n"] = model.groupby(["host", "year", "districts"])["host"].transform("size")

    # Phenotype denominator and sequencing fraction in the same host-year-district stratum.
    ph = pheno.copy()
    ph["districts"] = ph["districts"].fillna("Unknown").astype(str)
    ph_counts = (
        ph.groupby(["host", "Year", "districts"], observed=True)
        .size()
        .rename("phenotype_stratum_n")
        .reset_index()
        .rename(columns={"Year": "year"})
    )
    model = model.merge(ph_counts, on=["host", "year", "districts"], how="left", validate="m:1")
    model["wgs_sampling_fraction"] = model["wgs_stratum_n"] / model["phenotype_stratum_n"]

    density = (
        model.groupby(["host", "year", "period", "districts"], observed=True)
        .agg(
            n_WGS=("pcg_isolate_id", "size"),
            n_phenotype=("phenotype_stratum_n", "first"),
            n_localized=("localized_pcg", "sum"),
            n_any_PCG=("any_pcg", "sum"),
        )
        .reset_index()
    )
    density["WGS_sampling_fraction"] = density["n_WGS"] / density["n_phenotype"]
    density["localized_proportion"] = density["n_localized"] / density["n_WGS"]

    weighted_rows = []
    for (host, period), d in model.groupby(["host", "period"], observed=True):
        for outcome in ["localized_pcg", "any_pcg"]:
            for weight in ["calibration_weight", "trimmed_weight"]:
                p, neff, se = weighted_prop(d, outcome, weight)
                weighted_rows.append(
                    {
                        "host": host,
                        "period": period,
                        "outcome": outcome,
                        "weight": weight,
                        "n_unweighted": len(d),
                        "sum_weights": d[weight].sum(),
                        "effective_sample_size": neff,
                        "weighted_proportion": p,
                        "approx_CI_low": max(0, p - 1.96 * se),
                        "approx_CI_high": min(1, p + 1.96 * se),
                    }
                )
    weighted_props = pd.DataFrame(weighted_rows)

    boots = []
    for outcome in ["localized_pcg", "any_pcg"]:
        for weight in ["calibration_weight", "trimmed_weight"]:
            boots.append(farm_cluster_bootstrap(model, outcome, weight))
    bootstrap = pd.concat(boots, ignore_index=True)

    # Farm random-intercept logistic models in R/glmmTMB.
    with tempfile.TemporaryDirectory(prefix="pcg_review2_") as td:
        td = Path(td)
        data_csv = td / "model_input.csv"
        result_csv = td / "glmm_results.csv"
        r_script = td / "fit_models.R"
        model.to_csv(data_csv, index=False)
        r_script.write_text(
            r'''
args <- commandArgs(trailingOnly=TRUE)
d <- read.csv(args[1], stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(glmmTMB))
suppressPackageStartupMessages(library(sandwich))
out <- list(); k <- 1
for (h in c("Poultry", "Swine")) {
  z <- d[d$host == h & !is.na(d$farm_id) & !is.na(d$districts), ]
  z$districts <- factor(z$districts)
  z$farm_id <- factor(z$farm_id)
  z$log_density_z <- as.numeric(scale(log1p(z$wgs_stratum_n)))
  z$sampling_fraction_z <- as.numeric(scale(z$wgs_sampling_fraction))
  for (y in c("localized_pcg", "any_pcg")) {
    for (model_name in c("density_adjusted", "density_and_fraction_adjusted")) {
      if (model_name == "density_adjusted") {
        form <- as.formula(paste0(y, " ~ post + log_density_z + districts + (1|farm_id)"))
      } else {
        form <- as.formula(paste0(y, " ~ post + log_density_z + sampling_fraction_z + districts + (1|farm_id)"))
      }
      fit <- try(glmmTMB(form, data=z, family=binomial()), silent=TRUE)
      if (inherits(fit, "try-error")) {
        out[[k]] <- data.frame(host=h, outcome=y, model=model_name, n=nrow(z), n_farms=nlevels(z$farm_id),
                               OR=NA, CI_low=NA, CI_high=NA, p_value=NA,
                               random_intercept_variance=NA, convergence="failed")
      } else {
        co <- summary(fit)$coefficients$cond
        b <- co["post", "Estimate"]; se <- co["post", "Std. Error"]
        rv <- NA_real_
        vc <- VarCorr(fit)$cond
        if ("farm_id" %in% names(vc)) rv <- as.numeric(vc$farm_id[1,1])
        pdh <- fit$sdr$pdHess
        conv <- if (is.null(pdh)) "not_reported" else if (isTRUE(pdh)) "converged" else "non-positive-definite Hessian"
        out[[k]] <- data.frame(host=h, outcome=y, model=model_name, n=nrow(z), n_farms=nlevels(z$farm_id),
                               OR=exp(b), CI_low=exp(b-1.96*se), CI_high=exp(b+1.96*se),
                               p_value=co["post", "Pr(>|z|)"], random_intercept_variance=rv,
                               convergence=conv)
      }
      k <- k + 1
    }
    # Joint sensitivity analysis: calibration weighting + sequencing-density
    # adjustment + farm-clustered robust standard errors.
    for (wname in c("calibration_weight", "trimmed_weight")) {
      formw <- as.formula(paste0(y, " ~ post + log_density_z + sampling_fraction_z + districts"))
      fitw <- glm(formw, data=z, family=quasibinomial(), weights=z[[wname]])
      V <- sandwich::vcovCL(fitw, cluster=z$farm_id, type="HC0")
      b <- coef(fitw)["post"]; se <- sqrt(diag(V))["post"]
      out[[k]] <- data.frame(host=h, outcome=y,
                             model=paste0(wname, "_density_fraction_farm_clustered"),
                             n=nrow(z), n_farms=nlevels(z$farm_id), OR=exp(b),
                             CI_low=exp(b-1.96*se), CI_high=exp(b+1.96*se),
                             p_value=2*pnorm(-abs(b/se)), random_intercept_variance=NA,
                             convergence="robust_clustered")
      k <- k + 1
    }
  }
}
write.csv(do.call(rbind, out), args[2], row.names=FALSE)
''',
            encoding="utf-8",
        )
        cmd = [shutil.which("Rscript") or "Rscript", str(r_script), str(data_csv), str(result_csv)]
        run = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if run.returncode != 0 or not result_csv.exists():
            raise RuntimeError(f"R model failed\nSTDOUT:\n{run.stdout}\nSTDERR:\n{run.stderr}")
        glmm = pd.read_csv(result_csv)

    validation = pd.DataFrame(
        [
            {"check": "PCG threshold", "value": 10, "note": "Primary Ward.D2 height cut"},
            {"check": "Total PCG-analyzable isolates", "value": len(assn), "note": "SNP/PCG assignment denominator"},
            {"check": "Distinct PCGs", "value": assn["PCG_id"].nunique(), "note": "Expected 494"},
            {"check": "Localized PCGs", "value": pcs["localized_same_year_host_district"].sum(), "note": "Expected 308"},
            {"check": "Phenotype-linked WGS isolates with weights", "value": len(selected), "note": "2010-2023 overlap population"},
            {"check": "Uniquely linked isolates used in adjusted models", "value": len(model), "note": "After resolving aliases and excluding unmatched PCG IDs"},
            {"check": "Farms in adjusted model", "value": model["farm_id"].nunique(), "note": "Farm random intercept / cluster bootstrap"},
            {"check": "Missing calibration weights", "value": model["calibration_weight"].isna().sum(), "note": "Expected 0"},
        ]
    )

    core_rows = []
    for host in ["Poultry", "Swine"]:
        pre = denominators[(denominators.host.eq(host)) & (denominators.period.eq("pre-2017"))].iloc[0]
        post = denominators[(denominators.host.eq(host)) & (denominators.period.eq("post-2017"))].iloc[0]
        tst = tests[(tests.host.eq(host)) & (tests.outcome.eq("localized_pcg"))].iloc[0]
        wp_pre = weighted_props[
            weighted_props.host.eq(host)
            & weighted_props.period.eq("pre-2017")
            & weighted_props.outcome.eq("localized_pcg")
            & weighted_props.weight.eq("calibration_weight")
        ].iloc[0]
        wp_post = weighted_props[
            weighted_props.host.eq(host)
            & weighted_props.period.eq("post-2017")
            & weighted_props.outcome.eq("localized_pcg")
            & weighted_props.weight.eq("calibration_weight")
        ].iloc[0]
        bt = bootstrap[
            bootstrap.host.eq(host)
            & bootstrap.outcome.eq("localized_pcg")
            & bootstrap.weight.eq("calibration_weight")
        ].iloc[0]
        jm = glmm[
            glmm.host.eq(host)
            & glmm.outcome.eq("localized_pcg")
            & glmm.model.eq("calibration_weight_density_fraction_farm_clustered")
        ].iloc[0]
        core_rows.append(
            {
                "host": host,
                "pre_localized_isolates_over_total": f"{int(pre.n_isolates_in_localized_PCG)}/{int(pre.n_PCG_analyzable_WGS_isolates)}",
                "pre_unadjusted_proportion": pre.proportion_in_localized_PCG,
                "post_localized_isolates_over_total": f"{int(post.n_isolates_in_localized_PCG)}/{int(post.n_PCG_analyzable_WGS_isolates)}",
                "post_unadjusted_proportion": post.proportion_in_localized_PCG,
                "unadjusted_FDR_P": tst.p_FDR_BH_within_outcome,
                "pre_calibration_weighted_proportion": wp_pre.weighted_proportion,
                "post_calibration_weighted_proportion": wp_post.weighted_proportion,
                "weighted_post_minus_pre": wp_post.weighted_proportion - wp_pre.weighted_proportion,
                "farm_bootstrap_difference_CI_low": bt.ci_low,
                "farm_bootstrap_difference_CI_high": bt.ci_high,
                "joint_adjusted_OR_post_vs_pre": jm.OR,
                "joint_adjusted_OR_CI_low": jm.CI_low,
                "joint_adjusted_OR_CI_high": jm.CI_high,
                "joint_adjusted_P": jm.p_value,
            }
        )
    core_table = pd.DataFrame(core_rows)

    readme = pd.DataFrame(
        {
            "item": [
                "Primary outcome",
                "Denominator",
                "Unadjusted comparison",
                "Sampling-density adjustment",
                "Phenotype-guided selection adjustment",
                "Repeated sampling",
                "Interpretation",
            ],
            "description": [
                "Proportion of PCG-analyzable sequenced isolates assigned to a localized PCG.",
                "All isolates with a valid 10-SNP PCG/singleton assignment in each host-period stratum.",
                "Two-sided Pearson chi-square test; BH-FDR is applied across the two host comparisons for each outcome.",
                "Host-stratified binomial GLMM adjusted for log host-year-district WGS count and district; a sensitivity model additionally adjusts for the WGS/phenotype sampling fraction.",
                "Existing raking calibration weights balance year, host, district, farm, MDR, resistance-class count and available resistance phenotypes in the phenotype-WGS overlap population.",
                "Farm is included as a random intercept in GLMMs; weighted differences use a farm-cluster bootstrap (2,000 replicates).",
                "Raw PCG counts are descriptive only. A lower count is not interpreted as direct evidence of reduced dissemination.",
            ],
        }
    )

    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
        readme.to_excel(writer, sheet_name="README", index=False)
        core_table.to_excel(writer, sheet_name="Publication_core_table", index=False)
        denominators.to_excel(writer, sheet_name="Period_denominators", index=False)
        tests.to_excel(writer, sheet_name="Unadjusted_tests", index=False)
        glmm.to_excel(writer, sheet_name="Adjusted_GLMM", index=False)
        weighted_props.to_excel(writer, sheet_name="Weighted_proportions", index=False)
        bootstrap.to_excel(writer, sheet_name="Farm_bootstrap", index=False)
        density.to_excel(writer, sheet_name="Sampling_density", index=False)
        validation.to_excel(writer, sheet_name="Validation", index=False)
        model.to_excel(writer, sheet_name="Analysis_dataset", index=False)

    summary_lines = [
        "Second-round PCG sampling-intensity analysis",
        f"Output: {OUT_XLSX}",
        "",
        "Normalized denominators:",
        denominators.to_string(index=False),
        "",
        "Unadjusted tests:",
        tests.to_string(index=False),
        "",
        "Adjusted GLMM:",
        glmm.to_string(index=False),
        "",
        "Weighted proportions:",
        weighted_props.to_string(index=False),
        "",
        "Farm-cluster bootstrap:",
        bootstrap.to_string(index=False),
        "",
        "Validation:",
        validation.to_string(index=False),
    ]
    (OUTDIR / "analysis_summary.txt").write_text("\n".join(summary_lines), encoding="utf-8")
    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()
