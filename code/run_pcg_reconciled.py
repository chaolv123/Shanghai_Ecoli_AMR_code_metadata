"""Reconcile the author's primary PCG assignments with SNP sensitivity cuts.

Primary <=10-SNP groups come from pcg_original_assignments.xlsx, column 2.
Column 3 is a per-isolate flag and is retained for audit, not a group ID.
Other cutoffs use Ward.D2-equivalent hierarchical clustering and are explicitly
labelled as sensitivity analyses. Pairwise threshold violations are reported.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "06_Metadata_and_Input_Data"
OUT = ROOT / "outputs" / "extended" / "pcg_reconciled"
THRESHOLDS = [5, 10, 20, 25]


def snp_distance_path() -> Path:
    zipped = RAW / "core_snp_distances.csv.zip"
    plain = RAW / "core_snp_distances.csv"
    return zipped if zipped.exists() else plain


def district(value: object) -> str:
    if pd.isna(value):
        return "unknown"
    result = str(value).strip().lower()
    return result if result not in {"", "na", "nan", "none", "unknown"} else "unknown"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    original = pd.read_excel(RAW / "pcg_original_assignments.xlsx", sheet_name="ALL_ANIMALS_result")
    original = original.rename(columns={"within_Cluster": "author_PCG_id", "within_Cluster.1": "author_isolate_flag"})
    distance = pd.read_csv(snp_distance_path(), index_col=0)
    ids = distance.index.astype(str)
    if not ids.equals(distance.columns.astype(str)) or set(ids) != set(original.ID.astype(str)):
        raise ValueError("Distance matrix and original PCG table IDs differ")
    original = original.set_index("ID").loc[ids]
    meta = pd.read_excel(RAW / "pcg_metadata.xlsx").set_index("isolate_id").loc[ids]
    matrix = distance.to_numpy(dtype=np.float64)
    dendrogram = linkage(squareform(matrix, checks=False), method="ward")

    summaries, memberships, burden, comparison = [], [], [], []
    for threshold in THRESHOLDS:
        labels = original.author_PCG_id.to_numpy() if threshold == 10 else fcluster(dendrogram, threshold, criterion="distance")
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
        summaries.append({"threshold": threshold, "algorithm": "author assignment" if threshold == 10 else "Ward.D2 height cut",
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
    author_group = primary[primary.PCG_type.ne("singleton")].groupby("cluster_id").PCG_type.first()
    if len(author_group) != 494 or author_group.value_counts().to_dict() != {"localized": 308, "non-localized": 186}:
        raise ValueError("The author-assigned primary PCGs no longer match the audited 494/308/186 result")
    overlap_note = "Figure source-data workbooks are not included in this public code/data release."
    missing_note = "Figure-source validation was performed during revision; figure source-data files are not part of this public upload package."
    notes = pd.DataFrame([
        {"item": "Original PCGs", "finding": "494 distinct groups with >=2 isolates; 308 localized and 186 non-localized."},
        {"item": "Source-list duplicate", "finding": overlap_note},
        {"item": "Missing localized isolate", "finding": missing_note},
        {"item": "Primary assignment", "finding": "Author's original column-2 group labels are preserved at <=10 SNP. All primary groups have maximum pairwise SNP distance <=10."},
        {"item": "Alternative thresholds", "finding": "Ward.D2 hierarchical height cuts at 5, 20 and 25; pairwise-threshold violations are explicitly counted."},
        {"item": "Third original column", "finding": "Per-isolate 0/1 flag, not an independent PCG ID. It is retained but not used to define groups."},
    ])
    with pd.ExcelWriter(OUT / "Supplementary_PCG_threshold_reconciled.xlsx") as writer:
        pd.DataFrame(summaries).to_excel(writer, "Threshold_summary", index=False)
        pd.DataFrame(burden).to_excel(writer, "ARG_burden", index=False)
        pd.DataFrame(comparison).to_excel(writer, "PCG_pairwise_distances", index=False)
        pd.concat(memberships, ignore_index=True).to_excel(writer, "Cluster_assignments", index=False)
        notes.to_excel(writer, "Notes", index=False)
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
