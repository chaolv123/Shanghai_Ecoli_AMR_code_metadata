from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "New_analysis" / "new_add_2"
WEIGHTS = BASE / "WGS_sampling_bias_analysis" / "WGS_calibration_weights.csv"
PHENO = BASE / "metadata2.xlsx"
GENOME = ROOT / "New_analysis" / "new_add_6" / "metadata6.xlsx"
OUTPUT = BASE / "WGS_sampling_bias_analysis" / "Supplementary_weighted_ARG_burden_updated.xlsx"
SEED = 20260917


def weighted_mean(x, w):
    return float(np.average(x, weights=w))


def main():
    weights = pd.read_csv(WEIGHTS)
    pheno = pd.read_excel(PHENO)[["isolate_id", "WGSisolateid", "WGS"]]
    genome = pd.read_excel(GENOME)[["isolate_id", "year", "host", "period", "ARG_count"]]
    selected = pheno[pheno.WGS.eq(1)].merge(weights, on="isolate_id", validate="one_to_one")
    linked = selected.merge(genome, left_on="WGSisolateid", right_on="isolate_id", validate="one_to_one", suffixes=("_pheno", "_genome"))
    if len(linked) != len(selected):
        raise ValueError("Some WGS selected phenotype records did not match genome IDs")
    discordant_years = int(linked.year_pheno.ne(linked.year_genome).sum())
    linked = linked[linked.year_pheno.eq(linked.year_genome)].copy()
    linked["year"] = linked["year_genome"]
    linked["host"] = linked["host_genome"]
    linked["period"] = np.where(linked["year"] >= 2017, "post-2017", "pre-2017")
    rng = np.random.default_rng(SEED)
    rows = []
    for host in ["Poultry", "Swine"]:
        source = linked[linked.host.eq(host)]
        per = {}
        for period in ["pre-2017", "post-2017"]:
            part = source[source.period.eq(period)]
            per[period] = part
            rows.append({"host": host, "period": period, "n": len(part),
                         "unweighted_mean_ARG_count": part.ARG_count.mean(),
                         "calibrated_mean_ARG_count": weighted_mean(part.ARG_count, part.calibration_weight),
                         "trimmed_mean_ARG_count": weighted_mean(part.ARG_count, part.trimmed_weight)})
        pre, post = per["pre-2017"], per["post-2017"]
        boot = []
        for _ in range(1000):
            bpre = pre.iloc[rng.integers(len(pre), size=len(pre))]
            bpost = post.iloc[rng.integers(len(post), size=len(post))]
            boot.append(weighted_mean(bpost.ARG_count, bpost.calibration_weight) - weighted_mean(bpre.ARG_count, bpre.calibration_weight))
        rows.append({"host": host, "period": "post_minus_pre", "n": len(source),
                     "unweighted_mean_ARG_count": post.ARG_count.mean() - pre.ARG_count.mean(),
                     "calibrated_mean_ARG_count": weighted_mean(post.ARG_count, post.calibration_weight) - weighted_mean(pre.ARG_count, pre.calibration_weight),
                     "trimmed_mean_ARG_count": weighted_mean(post.ARG_count, post.trimmed_weight) - weighted_mean(pre.ARG_count, pre.trimmed_weight),
                     "bootstrap_CI_low": np.quantile(boot, .025), "bootstrap_CI_high": np.quantile(boot, .975)})
    result = pd.DataFrame(rows)
    qc = pd.DataFrame([
        {"check": "phenotype_WGS_linked_to_genome", "value": len(selected)},
        {"check": "excluded_discordant_year", "value": discordant_years},
        {"check": "genomes_without_phenotype", "value": len(genome) - len(selected)},
        {"check": "weighting_target", "value": "2010-2023 phenotype overlap: farms and districts with WGS support"},
        {"check": "bootstrap", "value": "1000 stratified isolate resamples within host and period; fixed calibration weights"},
    ])
    readme = pd.DataFrame([
        {"Field": "Title", "Description": "Supplementary Table. Calibration-weighted sensitivity analysis of acquired ARG burden in the WGS subset"},
        {"Field": "Cohort", "Description": "Phenotypically linked WGS isolates from 2010-2023 in farms and districts with WGS support; one isolate with discrepant sampling year was excluded."},
        {"Field": "Weights", "Description": "Raking weights calibrated to the phenotype overlap population by year, host, district, farm, MDR, resistance class count and available antimicrobial phenotypes."},
        {"Field": "Outcome", "Description": "Mean ARG_count by host and period; post_minus_pre is the post-2017 minus pre-2017 difference."},
        {"Field": "Uncertainty", "Description": "Percentile intervals from 1000 isolate bootstrap resamples within host and period, retaining calibration weights fixed; these do not account for uncertainty in the weights or farm clustering."},
        {"Field": "Limit", "Description": "Weights cannot recover genomic outcomes in farms or years with no WGS support, and no plasmid burden or PT information was provided."},
    ])
    with pd.ExcelWriter(OUTPUT) as writer:
        readme.to_excel(writer, "README", index=False)
        result.to_excel(writer, "ARG_burden_weighted", index=False)
        qc.to_excel(writer, "Data_validation", index=False)
    print(result.round(4).to_string(index=False))
    print(OUTPUT)


if __name__ == "__main__":
    main()
