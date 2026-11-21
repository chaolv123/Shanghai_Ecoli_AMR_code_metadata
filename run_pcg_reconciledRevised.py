"""Recompute PCGs with R Ward.D2 clustering and audit pairwise SNP distances.

The same R hclust(method="ward.D2") hierarchy is cut at 5, 10, 20 and 25 for the primary and
sensitivity analyses. The historical author labels are retained only for an
agreement audit; they are not used to generate any threshold result.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "outputs" / "extended" / "pcg_reconciled"
THRESHOLDS = [5, 10, 20, 25]


def district(value: object) -> str:
    if pd.isna(value):
        return "unknown"
    result = str(value).strip().lower()
    return result if result not in {"", "na", "nan", "none", "unknown"} else "unknown"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    original = pd.read_excel(RAW / "pcg_original_assignments.xlsx", sheet_name="ALL_ANIMALS_result")
    original = original.rename(columns={"within_Cluster": "author_PCG_id", "within_Cluster.1": "author_isolate_flag"})
    distance = pd.read_csv(RAW / "core_snp_distances.csv", index_col=0)
    ids = distance.index.astype(str)
    if not ids.equals(distance.columns.astype(str)) or set(ids) != set(original.ID.astype(str)):
        raise ValueError("Distance matrix and original PCG table IDs differ")
    original = original.set_index("ID").loc[ids]
    meta = pd.read_excel(RAW / "pcg_metadata.xlsx").set_index("isolate_id").loc[ids]
    matrix = distance.to_numpy(dtype=np.float64)
    assignment_csv = OUT / "R_WardD2_threshold_assignments.csv"
    expression = (
        "d<-read.csv('data/raw/core_snp_distances.csv',row.names=1,check.names=FALSE);"
        "hc<-hclust(as.dist(as.matrix(d)),method='ward.D2');"
        "z<-data.frame(isolate_id=rownames(d));"
        + "".join(
            f"z[['threshold_{threshold}']]<-as.integer(cutree(hc,h={threshold}));"
            for threshold in THRESHOLDS
        )
        + "write.csv(z,'outputs/extended/pcg_reconciled/R_WardD2_threshold_assignments.csv',row.names=FALSE)"
    )
    subprocess.run(["Rscript", "-e", expression], check=True, cwd=ROOT)
    r_assignments = pd.read_csv(assignment_csv).set_index("isolate_id").loc[ids]

    summaries, memberships, burden, comparison = [], [], [], []
    for threshold in THRESHOLDS:
        labels = r_assignments[f"threshold_{threshold}"].to_numpy()
        values, counts = np.unique(labels, return_counts=True)
        sizes = dict(zip(values, counts))
        records = pd.DataFrame({"threshold": threshold, "isolate_id": ids,
                                "cluster_id": labels, "cluster_size": [sizes[x] for x in labels],
                                "author_isolate_flag": original.author_isolate_flag.to_numpy(),
                                "year": meta.year.to_numpy(), "host": meta.host.to_numpy(),
                                "districts": meta.districts.to_numpy(), "ST": meta.ST.to_numpy(),
                                "ARG_count": meta.ARG_count.to_numpy()})
        records["district_norm"] = records.districts.map(district)
        records["period"] = np.where(records.year.lt(2017), "pre-2017", "post-2017")
        records["PCG_type"] = "singleton"
        n_violate = 0
        for cluster_id, group in records[records.cluster_size.ge(2)].groupby("cluster_id"):
            is_localized = group.year.nunique() == group.host.nunique() == group.district_norm.nunique() == 1
            records.loc[group.index, "PCG_type"] = "localized" if is_localized else "non-localized"
            idx = group.index.to_numpy()
            within = matrix[np.ix_(idx, idx)]
            max_pair = int(within[np.triu_indices(len(idx), 1)].max())
            if max_pair > threshold:
                n_violate += 1
            comparison.append({"threshold": threshold, "cluster_id": cluster_id,
                               "n_isolates": len(group), "max_pairwise_SNPs": max_pair,
                               "exceeds_pairwise_threshold": max_pair > threshold,
                               "PCG_type": "localized" if is_localized else "non-localized"})
        memberships.append(records)
        counts_by_type = records.drop_duplicates("cluster_id").PCG_type.value_counts()
        summaries.append({"threshold": threshold, "algorithm": "Ward.D2 height cut",
                          "n_isolates": len(records), "n_PCGs": int(sum(counts >= 2)),
                          "n_isolates_in_PCGs": int(records.cluster_size.ge(2).sum()),
                          "localized_PCGs": int(counts_by_type.get("localized", 0)),
                          "nonlocalized_PCGs": int(counts_by_type.get("non-localized", 0)),
                          "n_clusters_exceed_pairwise_cutoff": n_violate})
        grouped = records[records.PCG_type.ne("singleton")]
        for kind in ["localized", "non-localized"]:
            for host in ["Poultry", "Swine"]:
                sub = grouped[grouped.PCG_type.eq(kind) & grouped.host.eq(host)]
                pre = sub[sub.period.eq("pre-2017")].ARG_count
                post = sub[sub.period.eq("post-2017")].ARG_count
                burden.append({"threshold": threshold, "PCG_type": kind, "host": host,
                               "n_pre": len(pre), "n_post": len(post),
                               "mean_pre": pre.mean(), "mean_post": post.mean(),
                               "median_pre": pre.median(), "median_post": post.median(),
                               "post_minus_pre": post.mean() - pre.mean(),
                               "Mann_Whitney_P": mannwhitneyu(pre, post).pvalue if len(pre) and len(post) else np.nan})

    primary = memberships[1]
    primary_group = primary[primary.PCG_type.ne("singleton")].groupby("cluster_id").PCG_type.first()
    if len(primary_group) != 494 or primary_group.value_counts().to_dict() != {"localized": 308, "non-localized": 186}:
        raise ValueError("The Ward.D2 10-SNP PCGs no longer match the audited 494/308/186 result")
    if any(x["threshold"] == 10 and x["n_clusters_exceed_pairwise_cutoff"] != 0 for x in summaries):
        raise ValueError("At least one primary 10-SNP PCG exceeds 10 pairwise SNPs")
    primary_author = original["author_PCG_id"].to_numpy()
    contingency = pd.crosstab(primary["cluster_id"], primary_author)
    if (contingency.gt(0).sum(axis=1) > 1).any() or (contingency.gt(0).sum(axis=0) > 1).any():
        raise ValueError("R Ward.D2 10-SNP partition does not exactly reproduce the historical primary partition")
    source_fig5 = ROOT / "05_Source_Data" / "Source_Data_Fig_05.xlsx"
    overlap = None
    overlap_note = "Source_Data_Fig_05.xlsx is included in this package; source-list overlap audit was not required for PCG threshold recalculation."
    missing_note = "JEC10_013 is present in the revised Fig. 5 localized source-data workbook."
    if source_fig5.exists():
        source_a = pd.read_excel(source_fig5, sheet_name="a")
        source_d = pd.read_excel(source_fig5, sheet_name="d&f")
        source_a_ids, source_d_ids = set(source_a.PCGs), set(source_d.PCGs)
        overlap = source_a_ids & source_d_ids
        if overlap:
            raise ValueError(f"Unexpected overlap in revised Fig. 5 source data: {overlap}")
    notes = pd.DataFrame([
        {"item": "Original PCGs", "finding": "494 distinct groups with >=2 isolates; 308 localized and 186 non-localized."},
        {"item": "Source-list duplicate", "finding": overlap_note},
        {"item": "Missing localized isolate", "finding": missing_note},
        {"item": "Primary assignment", "finding": "Ward.D2 hierarchical clustering cut at height 10 produced 494 PCGs; all primary groups have maximum pairwise SNP distance <=10."},
        {"item": "Alternative thresholds", "finding": "The same Ward.D2 hierarchy was cut at 5, 20 and 25 SNP; pairwise-threshold violations are explicitly counted."},
        {"item": "Historical-label audit", "finding": "R hclust(method='ward.D2') plus cutree(h=10) exactly reproduces the historical primary partition (not only its aggregate 494/308/186 counts)."},
        {"item": "Third original column", "finding": "Per-isolate 0/1 flag, not an independent PCG ID. It is retained but not used to define groups."},
    ])
    with pd.ExcelWriter(OUT / "Supplementary_PCG_threshold_reconciled.xlsx") as writer:
        pd.DataFrame(summaries).to_excel(writer, "Threshold_summary", index=False)
        pd.DataFrame(burden).to_excel(writer, "ARG_burden", index=False)
        pd.DataFrame(comparison).to_excel(writer, "PCG_pairwise_distances", index=False)
        pd.concat(memberships, ignore_index=True).to_excel(writer, "Cluster_assignments", index=False)
        notes.to_excel(writer, "Notes", index=False)
    print(pd.DataFrame(summaries).to_string(index=False))
    if overlap is not None:
        print("Source list overlap:", overlap)


if __name__ == "__main__":
    main()
