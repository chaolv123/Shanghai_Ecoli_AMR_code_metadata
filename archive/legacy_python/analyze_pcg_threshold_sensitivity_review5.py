from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


BASE = Path("New_analysis") / "new_add_5"
DIST_PATH = BASE / "all_animals_dists.csv"
META_PATH = BASE / "metadata3.xlsx"
OUT_DIR = BASE / "review5_pcg_threshold_sensitivity"
OUT_XLSX = OUT_DIR / "Supplementary_PCG_threshold_sensitivity.xlsx"
OUT_SUMMARY = OUT_DIR / "review5_pcg_threshold_sensitivity_summary.txt"

THRESHOLDS = [5, 10, 20, 25]
UNKNOWN_DISTRICTS = {"", "-", "NA", "N/A", "UNKNOWN", "NAN", "NONE"}


class UnionFind:
    def __init__(self, n: int) -> None:
        self.parent = np.arange(n, dtype=np.int32)
        self.rank = np.zeros(n, dtype=np.int8)

    def find(self, x: int) -> int:
        parent = self.parent
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return int(x)

    def union(self, a: int, b: int) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1


def norm_district(x) -> str:
    if pd.isna(x):
        return "Unknown"
    s = str(x).strip()
    if s.upper() in UNKNOWN_DISTRICTS:
        return "Unknown"
    return s[:1].upper() + s[1:].lower()


def classify_pcg(g: pd.DataFrame) -> dict:
    known_districts = sorted([d for d in g["districts_clean"].dropna().unique() if d != "Unknown"])
    all_districts = sorted(g["districts_clean"].dropna().unique())
    years = sorted(g["year"].dropna().astype(int).unique())
    hosts = sorted(g["host"].dropna().unique())
    periods = sorted(g["period"].dropna().unique())
    sts = sorted(g["ST"].dropna().astype(str).unique())
    has_unknown_district = "Unknown" in set(all_districts)
    localized_known = len(years) == 1 and len(known_districts) == 1 and not has_unknown_district
    return {
        "n_isolates": len(g),
        "n_years": len(years),
        "first_year": min(years) if years else np.nan,
        "last_year": max(years) if years else np.nan,
        "year_span": (max(years) - min(years)) if years else np.nan,
        "n_hosts": len(hosts),
        "hosts": ";".join(hosts),
        "n_known_districts": len(known_districts),
        "known_districts": ";".join(known_districts),
        "has_unknown_district": has_unknown_district,
        "n_periods": len(periods),
        "periods": ";".join(periods),
        "n_STs": len(sts),
        "STs": ";".join(sts[:20]),
        "mean_ARG_count": g["ARG_count"].mean(),
        "median_ARG_count": g["ARG_count"].median(),
        "localized_known_same_year_district": localized_known,
        "nonlocalized_or_unknown": not localized_known,
        "multi_year": len(years) > 1,
        "multi_period": len(periods) > 1,
        "cross_host": len(hosts) > 1,
        "cross_known_district": len(known_districts) > 1,
        "single_known_district_multiyear": len(known_districts) == 1 and len(years) > 1 and not has_unknown_district,
    }


def make_components(dist: np.ndarray, ids: list[str], threshold: int) -> pd.DataFrame:
    n = len(ids)
    uf = UnionFind(n)
    mask = np.triu(dist <= threshold, k=1)
    rows, cols = np.where(mask)
    for a, b in zip(rows, cols):
        uf.union(int(a), int(b))
    roots = [uf.find(i) for i in range(n)]
    root_to_members: dict[int, list[int]] = defaultdict(list)
    for i, r in enumerate(roots):
        root_to_members[r].append(i)
    records = []
    pcg_number = 0
    for _, members in sorted(root_to_members.items(), key=lambda kv: (-len(kv[1]), min(kv[1]))):
        if len(members) < 2:
            for i in members:
                records.append({"threshold": threshold, "isolate_id": ids[i], "PCG_id": pd.NA, "component_size": 1})
            continue
        pcg_number += 1
        pcg_id = f"PCG{threshold}_{pcg_number:04d}"
        for i in members:
            records.append({"threshold": threshold, "isolate_id": ids[i], "PCG_id": pcg_id, "component_size": len(members)})
    return pd.DataFrame(records)


def mannwhitney_summary(values_pre: pd.Series, values_post: pd.Series) -> dict:
    values_pre = pd.to_numeric(values_pre, errors="coerce").dropna()
    values_post = pd.to_numeric(values_post, errors="coerce").dropna()
    out = {
        "n_pre": int(len(values_pre)),
        "n_post": int(len(values_post)),
        "mean_pre": values_pre.mean() if len(values_pre) else np.nan,
        "mean_post": values_post.mean() if len(values_post) else np.nan,
        "median_pre": values_pre.median() if len(values_pre) else np.nan,
        "median_post": values_post.median() if len(values_post) else np.nan,
        "change_mean_post_minus_pre": (values_post.mean() - values_pre.mean()) if len(values_pre) and len(values_post) else np.nan,
        "p_value": np.nan,
    }
    if len(values_pre) and len(values_post):
        out["p_value"] = stats.mannwhitneyu(values_pre, values_post, alternative="two-sided").pvalue
    return out


def district_network(clustered: pd.DataFrame, threshold: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    edge_counter = defaultdict(lambda: {"n_PCG": 0, "n_isolate_records": 0})
    for pcg_id, g0 in clustered.dropna(subset=["PCG_id"]).groupby("PCG_id"):
        for host in sorted(g0["host"].dropna().unique()):
            for period in sorted(g0["period"].dropna().unique()):
                g = g0[(g0["host"] == host) & (g0["period"] == period)]
                dists = sorted([d for d in g["districts_clean"].unique() if d != "Unknown"])
                if len(dists) < 2:
                    continue
                for a, b in combinations(dists, 2):
                    key = (threshold, host, period, a, b)
                    edge_counter[key]["n_PCG"] += 1
                    edge_counter[key]["n_isolate_records"] += len(g[g["districts_clean"].isin([a, b])])
    edges = []
    for (thr, host, period, a, b), v in edge_counter.items():
        edges.append({
            "threshold": thr,
            "host": host,
            "period": period,
            "district_a": a,
            "district_b": b,
            "shared_PCG_count": v["n_PCG"],
            "isolate_records_in_shared_PCGs": v["n_isolate_records"],
        })
    edge_df = pd.DataFrame(edges)
    if edge_df.empty:
        return edge_df, pd.DataFrame()

    node_rows = []
    for (thr, host, period), e in edge_df.groupby(["threshold", "host", "period"]):
        nodes = sorted(set(e["district_a"]).union(set(e["district_b"])))
        for node in nodes:
            incident = e[(e["district_a"] == node) | (e["district_b"] == node)]
            node_rows.append({
                "threshold": thr,
                "host": host,
                "period": period,
                "district": node,
                "degree": int(len(incident)),
                "weighted_degree_shared_PCG": int(incident["shared_PCG_count"].sum()),
                "weighted_degree_isolate_records": int(incident["isolate_records_in_shared_PCGs"].sum()),
            })
    return edge_df, pd.DataFrame(node_rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta = pd.read_excel(META_PATH)
    meta["isolate_id"] = meta["isolate_id"].astype(str)
    meta["districts_clean"] = meta["districts"].map(norm_district)

    dist_df = pd.read_csv(DIST_PATH, index_col=0)
    ids = dist_df.index.astype(str).tolist()
    dist_df.columns = dist_df.columns.astype(str)
    dist = dist_df.to_numpy(dtype=np.int32, copy=True)

    matrix_ids = set(ids)
    meta_ids = set(meta["isolate_id"])
    common_ids = [i for i in ids if i in meta_ids]
    validation = pd.DataFrame([
        {"check": "distance_matrix_rows", "observed": dist_df.shape[0], "expected": "2611 per author note", "status": "INFO"},
        {"check": "distance_matrix_columns", "observed": dist_df.shape[1], "expected": "square matrix", "status": "PASS" if dist_df.shape[0] == dist_df.shape[1] else "FAIL"},
        {"check": "matrix_row_column_order_identical", "observed": list(dist_df.index) == list(dist_df.columns), "expected": True, "status": "PASS" if list(dist_df.index) == list(dist_df.columns) else "FAIL"},
        {"check": "metadata_rows", "observed": len(meta), "expected": "2611 per author note", "status": "INFO"},
        {"check": "metadata_unique_isolate_id", "observed": meta["isolate_id"].nunique(), "expected": len(meta), "status": "PASS" if meta["isolate_id"].nunique() == len(meta) else "FAIL"},
        {"check": "matrix_ids_without_metadata", "observed": len(matrix_ids - meta_ids), "expected": 0, "status": "PASS" if len(matrix_ids - meta_ids) == 0 else "CHECK"},
        {"check": "metadata_ids_without_matrix", "observed": len(meta_ids - matrix_ids), "expected": 0, "status": "PASS" if len(meta_ids - matrix_ids) == 0 else "CHECK"},
        {"check": "analysis_isolates_in_intersection", "observed": len(common_ids), "expected": "used for all stratified analyses", "status": "PASS"},
        {"check": "unknown_district_isolates", "observed": int((meta["districts_clean"] == "Unknown").sum()), "expected": "excluded from known-district network edges", "status": "INFO"},
    ])

    meta = meta.set_index("isolate_id").loc[common_ids].reset_index()

    all_assignments = []
    pcg_rows = []
    threshold_rows = []
    burden_rows = []
    st_rows = []
    all_edges = []
    all_nodes = []

    for threshold in THRESHOLDS:
        assign = make_components(dist, ids, threshold)
        assign = assign[assign["isolate_id"].isin(common_ids)].copy()
        all_assignments.append(assign)
        clustered = assign.merge(meta, on="isolate_id", how="left")
        pcg_only = clustered.dropna(subset=["PCG_id"]).copy()

        for pcg_id, g in pcg_only.groupby("PCG_id"):
            row = {"threshold": threshold, "PCG_id": pcg_id}
            row.update(classify_pcg(g))
            pcg_rows.append(row)

        pcg_df_thr = pd.DataFrame([r for r in pcg_rows if r["threshold"] == threshold])
        threshold_rows.append({
            "threshold": threshold,
            "n_isolates_total": len(common_ids),
            "n_PCGs": int(len(pcg_df_thr)),
            "n_isolates_in_PCGs": int(len(pcg_only)),
            "n_singletons": int(len(common_ids) - len(pcg_only)),
            "median_PCG_size": pcg_df_thr["n_isolates"].median() if len(pcg_df_thr) else np.nan,
            "max_PCG_size": pcg_df_thr["n_isolates"].max() if len(pcg_df_thr) else np.nan,
            "localized_known_same_year_district_PCGs": int(pcg_df_thr["localized_known_same_year_district"].sum()) if len(pcg_df_thr) else 0,
            "nonlocalized_or_unknown_PCGs": int(pcg_df_thr["nonlocalized_or_unknown"].sum()) if len(pcg_df_thr) else 0,
            "multi_year_PCGs": int(pcg_df_thr["multi_year"].sum()) if len(pcg_df_thr) else 0,
            "multi_period_PCGs": int(pcg_df_thr["multi_period"].sum()) if len(pcg_df_thr) else 0,
            "cross_host_PCGs": int(pcg_df_thr["cross_host"].sum()) if len(pcg_df_thr) else 0,
            "cross_known_district_PCGs": int(pcg_df_thr["cross_known_district"].sum()) if len(pcg_df_thr) else 0,
            "unknown_district_PCGs": int(pcg_df_thr["has_unknown_district"].sum()) if len(pcg_df_thr) else 0,
        })

        for pcg_type, mask in {
            "all_PCG_isolates": pcg_only["PCG_id"].notna(),
            "localized_known_same_year_district": pcg_only["PCG_id"].map(pcg_df_thr.set_index("PCG_id")["localized_known_same_year_district"]).fillna(False).astype(bool),
            "nonlocalized_or_unknown": pcg_only["PCG_id"].map(pcg_df_thr.set_index("PCG_id")["nonlocalized_or_unknown"]).fillna(False).astype(bool),
            "multi_period": pcg_only["PCG_id"].map(pcg_df_thr.set_index("PCG_id")["multi_period"]).fillna(False).astype(bool),
            "cross_host": pcg_only["PCG_id"].map(pcg_df_thr.set_index("PCG_id")["cross_host"]).fillna(False).astype(bool),
            "cross_known_district": pcg_only["PCG_id"].map(pcg_df_thr.set_index("PCG_id")["cross_known_district"]).fillna(False).astype(bool),
        }.items():
            sub_type = pcg_only[mask].copy()
            for host in ["All"] + sorted(sub_type["host"].dropna().unique().tolist()):
                sub = sub_type if host == "All" else sub_type[sub_type["host"] == host]
                if sub.empty:
                    continue
                out = mannwhitney_summary(
                    sub.loc[sub["period"] == "pre-2017", "ARG_count"],
                    sub.loc[sub["period"] == "post-2017", "ARG_count"],
                )
                out.update({
                    "threshold": threshold,
                    "PCG_type": pcg_type,
                    "host": host,
                    "n_PCGs_represented": int(sub["PCG_id"].nunique()),
                })
                burden_rows.append(out)

        for (host, period, st), g in pcg_only.groupby(["host", "period", "ST"], dropna=False):
            st_rows.append({
                "threshold": threshold,
                "host": host,
                "period": period,
                "ST": st,
                "n_isolates": int(len(g)),
                "n_PCGs": int(g["PCG_id"].nunique()),
                "mean_ARG_count": g["ARG_count"].mean(),
                "median_ARG_count": g["ARG_count"].median(),
            })

        edges, nodes = district_network(clustered, threshold)
        if not edges.empty:
            all_edges.append(edges)
        if not nodes.empty:
            all_nodes.append(nodes)

    assignments = pd.concat(all_assignments, ignore_index=True)
    pcg_summary = pd.DataFrame(pcg_rows)
    threshold_summary = pd.DataFrame(threshold_rows)
    burden = pd.DataFrame(burden_rows)
    st_comp = pd.DataFrame(st_rows)
    st_top = (
        st_comp.sort_values(["threshold", "host", "period", "n_PCGs", "n_isolates"], ascending=[True, True, True, False, False])
        .groupby(["threshold", "host", "period"], as_index=False)
        .head(10)
    )
    edges = pd.concat(all_edges, ignore_index=True) if all_edges else pd.DataFrame()
    nodes = pd.concat(all_nodes, ignore_index=True) if all_nodes else pd.DataFrame()
    node_top = (
        nodes.sort_values(["threshold", "host", "period", "weighted_degree_shared_PCG", "degree"], ascending=[True, True, True, False, False])
        .groupby(["threshold", "host", "period"], as_index=False)
        .head(5)
        if not nodes.empty else pd.DataFrame()
    )

    # Direction summary for key burden conclusions relative to the primary 10-SNP threshold.
    direction_rows = []
    key_types = ["all_PCG_isolates", "localized_known_same_year_district", "nonlocalized_or_unknown"]
    for pcg_type in key_types:
        for host in ["Poultry", "Swine"]:
            g = burden[(burden["PCG_type"] == pcg_type) & (burden["host"] == host)].copy()
            if g.empty:
                continue
            ref = g[g["threshold"] == 10]
            ref_change = ref["change_mean_post_minus_pre"].iloc[0] if not ref.empty else np.nan
            ref_sign = np.sign(ref_change) if pd.notna(ref_change) and ref_change != 0 else np.nan
            signs = []
            for _, r in g.iterrows():
                ch = r["change_mean_post_minus_pre"]
                signs.append(pd.notna(ref_sign) and np.sign(ch) == ref_sign)
            direction_rows.append({
                "PCG_type": pcg_type,
                "host": host,
                "change_5SNP": g.loc[g["threshold"] == 5, "change_mean_post_minus_pre"].iloc[0] if (g["threshold"] == 5).any() else np.nan,
                "change_10SNP": g.loc[g["threshold"] == 10, "change_mean_post_minus_pre"].iloc[0] if (g["threshold"] == 10).any() else np.nan,
                "change_20SNP": g.loc[g["threshold"] == 20, "change_mean_post_minus_pre"].iloc[0] if (g["threshold"] == 20).any() else np.nan,
                "change_25SNP": g.loc[g["threshold"] == 25, "change_mean_post_minus_pre"].iloc[0] if (g["threshold"] == 25).any() else np.nan,
                "same_direction_as_10SNP_all_thresholds": all(signs),
            })
    direction = pd.DataFrame(direction_rows)

    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
        validation.to_excel(writer, sheet_name="Data_validation", index=False)
        threshold_summary.to_excel(writer, sheet_name="Threshold_summary", index=False)
        pcg_summary.to_excel(writer, sheet_name="PCG_summary", index=False)
        burden.to_excel(writer, sheet_name="PCG_ARG_burden", index=False)
        direction.to_excel(writer, sheet_name="Burden_direction_summary", index=False)
        st_top.to_excel(writer, sheet_name="Top_ST_by_threshold", index=False)
        edges.to_excel(writer, sheet_name="District_network_edges", index=False)
        nodes.to_excel(writer, sheet_name="District_network_nodes", index=False)
        node_top.to_excel(writer, sheet_name="Top_connectivity_districts", index=False)
        assignments.to_excel(writer, sheet_name="Cluster_assignments", index=False)

    lines = []
    lines.append("Review 5 PCG threshold sensitivity analysis")
    lines.append("=" * 52)
    lines.append(f"Distance matrix: {dist_df.shape[0]} x {dist_df.shape[1]}; metadata rows: {len(meta)}; analysis isolates: {len(common_ids)}")
    lines.append(f"Unknown-district isolates: {(meta['districts_clean'] == 'Unknown').sum()}; these were excluded from known-district network edges.")
    lines.append("")
    lines.append("Threshold-level PCG counts")
    lines.append("--------------------------")
    for _, r in threshold_summary.iterrows():
        lines.append(
            f"<= {int(r['threshold'])} SNP: PCGs={int(r['n_PCGs'])}, isolates in PCGs={int(r['n_isolates_in_PCGs'])}, "
            f"singletons={int(r['n_singletons'])}, localized known same-year/district={int(r['localized_known_same_year_district_PCGs'])}, "
            f"nonlocalized/unknown={int(r['nonlocalized_or_unknown_PCGs'])}, multi-year={int(r['multi_year_PCGs'])}, "
            f"multi-period={int(r['multi_period_PCGs'])}, cross-host={int(r['cross_host_PCGs'])}, "
            f"cross-known-district={int(r['cross_known_district_PCGs'])}."
        )
    lines.append("")
    lines.append("ARG-burden changes in all PCG isolates")
    lines.append("--------------------------------------")
    for host in ["Poultry", "Swine"]:
        g = burden[(burden["PCG_type"] == "all_PCG_isolates") & (burden["host"] == host)].sort_values("threshold")
        for _, r in g.iterrows():
            lines.append(
                f"{host} <= {int(r['threshold'])} SNP: pre mean={r['mean_pre']:.3f}, post mean={r['mean_post']:.3f}, "
                f"change={r['change_mean_post_minus_pre']:.3f}, P={r['p_value']:.3g}, n_pre={int(r['n_pre'])}, n_post={int(r['n_post'])}."
            )
    lines.append("")
    lines.append("Direction robustness relative to <=10 SNP")
    lines.append("-----------------------------------------")
    for _, r in direction.iterrows():
        lines.append(
            f"{r['host']} {r['PCG_type']}: changes <=5/<=10/<=20/<=25 = "
            f"{r['change_5SNP']:.3f}/{r['change_10SNP']:.3f}/{r['change_20SNP']:.3f}/{r['change_25SNP']:.3f}; "
            f"same direction={r['same_direction_as_10SNP_all_thresholds']}."
        )
    lines.append("")
    lines.append("Top connectivity districts at <=10 SNP")
    lines.append("--------------------------------------")
    top10 = node_top[node_top["threshold"] == 10] if not node_top.empty else pd.DataFrame()
    for _, r in top10.iterrows():
        lines.append(
            f"{r['host']} {r['period']}: {r['district']} degree={int(r['degree'])}, "
            f"weighted_degree_shared_PCG={int(r['weighted_degree_shared_PCG'])}."
        )
    lines.append("")
    lines.append("Interpretation")
    lines.append("--------------")
    lines.append("PCGs should be interpreted as putative genomic clusters. Without animal-movement, farm-contact, supply-chain or environmental data, district links indicate network connectivity among closely related isolates rather than direct transmission pathways.")
    OUT_SUMMARY.write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote {OUT_XLSX}")
    print(f"Wrote {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
