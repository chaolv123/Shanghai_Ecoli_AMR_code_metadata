"""Rebuild key resistome and plasmid summaries from public processed inputs.

This repository intentionally starts from processed isolate metadata, acquired-
ARG presence/absence calls, and predicted plasmid-type assignments. Raw reads
are outside the workflow boundary. Exact LEfSe scores are retained as a derived
input because the original LEfSe run environment is unavailable; all other
summaries in this module are recomputed from the supplied matrices.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.spatial.distance import pdist, squareform


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
OUT = ROOT / "outputs" / "derived_reproducibility"


def norm_id(value: object) -> str:
    text = str(value).upper().split("R1")[0]
    return re.sub(r"[^A-Z0-9]", "", text)


def pcoa(distance: np.ndarray, dimensions: int = 3) -> tuple[np.ndarray, np.ndarray]:
    n = distance.shape[0]
    d2 = distance**2
    row_mean = d2.mean(axis=1, keepdims=True)
    centered = -0.5 * (d2 - row_mean - row_mean.T + d2.mean())
    vals, vecs = eigh(centered, subset_by_index=[n - dimensions, n - 1])
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    coords = vecs * np.sqrt(np.clip(vals, 0, None))
    return coords, vals


def permanova_two_group(
    distance: np.ndarray, labels: np.ndarray, permutations: int, seed: int
) -> dict[str, float]:
    n = len(labels)
    groups = np.unique(labels)
    if len(groups) != 2:
        raise ValueError("This helper expects exactly two groups")
    total_ss = float(np.square(distance).sum() / (2 * n))

    def statistic(lab: np.ndarray) -> tuple[float, float]:
        within_ss = 0.0
        for group in groups:
            idx = np.flatnonzero(lab == group)
            within_ss += float(np.square(distance[np.ix_(idx, idx)]).sum() / (2 * len(idx)))
        between_ss = total_ss - within_ss
        r2 = between_ss / total_ss
        f_value = (between_ss / (len(groups) - 1)) / (within_ss / (n - len(groups)))
        return r2, f_value

    observed_r2, observed_f = statistic(labels)
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(permutations):
        _, perm_f = statistic(rng.permutation(labels))
        exceed += int(perm_f >= observed_f - 1e-12)
    return {
        "n": n,
        "R2": observed_r2,
        "pseudo_F": observed_f,
        "permutations": permutations,
        "permutation_p": (exceed + 1) / (permutations + 1),
    }


def run_resistome(permutations: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genomic = pd.read_excel(RAW / "genomic_metadata.xlsx")
    arg_cols = list(genomic.columns[8:])
    matrix = genomic[arg_cols].apply(pd.to_numeric, errors="coerce").fillna(0).gt(0).astype(np.uint8)
    pcoa_rows = []
    permanova_rows = []
    for host_index, host in enumerate(["Poultry", "Swine"]):
        keep = genomic["host"].eq(host).to_numpy()
        host_meta = genomic.loc[keep, ["isolate_id", "year", "period", "host"]].reset_index(drop=True)
        host_matrix = matrix.loc[keep].to_numpy()
        distances = squareform(pdist(host_matrix, metric="jaccard"))
        coords, eigenvalues = pcoa(distances)
        positive_total = max(float(eigenvalues[eigenvalues > 0].sum()), 1e-12)
        explained = eigenvalues / positive_total
        for i, row in host_meta.iterrows():
            pcoa_rows.append(
                {
                    **row.to_dict(),
                    "PCoA1": coords[i, 0],
                    "PCoA2": coords[i, 1],
                    "PCoA3": coords[i, 2],
                    "axis1_fraction_of_retained_positive_eigenvalues": explained[0],
                    "axis2_fraction_of_retained_positive_eigenvalues": explained[1],
                }
            )
        stats = permanova_two_group(
            distances,
            host_meta["period"].to_numpy(),
            permutations=permutations,
            seed=20261120 + host_index,
        )
        permanova_rows.append({"host": host, "contrast": "pre-2017 vs post-2017", **stats})

    source3 = DERIVED / "fig3_resistome_derived_input.xlsx"
    lefse = pd.read_excel(source3, sheet_name="d_lefse")
    expected = {"Poultry": 41, "Swine": 44}
    observed = lefse.groupby("Host").size().to_dict()
    if observed != expected:
        raise ValueError(f"LEfSe derived input count mismatch: {observed}")
    return pd.DataFrame(pcoa_rows), pd.DataFrame(permanova_rows), lefse


def run_plasmid() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    source6 = DERIVED / "fig6_plasmid_derived_input.xlsx"
    pt_arg = pd.read_excel(source6, sheet_name="a")
    pt_burden = pd.read_excel(source6, sheet_name="b")
    genomic = pd.read_excel(RAW / "genomic_metadata.xlsx")
    host_map = genomic.assign(match_id=genomic["isolate_id"].map(norm_id))[["match_id", "host"]]
    pt_arg["match_id"] = pt_arg["ID"].map(norm_id)
    linked = pt_arg.merge(host_map, on="match_id", how="left", validate="m:1")
    # Twelve legacy 2008-2009 poultry records use aliases absent from the
    # current genomic identifier column; their host is documented in Fig. 6.
    linked.loc[linked["host"].isna() & linked["ID"].astype(str).str.startswith("J"), "host"] = "Poultry"
    if linked["host"].isna().any():
        raise ValueError("Unresolved host aliases remain in plasmid input")

    richness_rows = []
    host_sets = {}
    for host, data in linked.groupby("host"):
        host_sets[host] = set(data["PTs"].dropna().astype(str))
        richness_rows.append({"host": host, "predicted_PT_richness": len(host_sets[host])})
    richness_rows.append(
        {
            "host": "Shared",
            "predicted_PT_richness": len(host_sets["Poultry"] & host_sets["Swine"]),
        }
    )
    richness = pd.DataFrame(richness_rows)
    values = dict(zip(richness.host, richness.predicted_PT_richness))
    if values != {"Poultry": 67, "Swine": 108, "Shared": 50}:
        raise ValueError(f"Predicted PT richness mismatch: {values}")

    prevalence_rows = []
    for host in ["Poultry", "Swine"]:
        sheet = pd.read_excel(source6, sheet_name=f"c_{host}", header=None)
        sheet.columns = range(sheet.shape[1])
        header = sheet.iloc[1].tolist()
        data = sheet.iloc[2:].copy()
        data.columns = header
        if "Plasmid" in data.columns and "PTs" not in data.columns:
            data = data.rename(columns={"Plasmid": "PTs"})
        data = data.loc[data["PTs"].notna()].copy()
        for c in ["Pre_2017_count", "Post_2017_count", "Pre_2017_prop_in_stage", "Post_2017_prop_in_stage", "P_value"]:
            data[c] = pd.to_numeric(data[c], errors="coerce")
        data = data.loc[
            data["Pre_2017_prop_in_stage"].between(0, 1, inclusive="both")
            & data["Post_2017_prop_in_stage"].between(0, 1, inclusive="both")
        ].copy()
        data.insert(0, "Host", host)
        prevalence_rows.append(data)
    prevalence = pd.concat(prevalence_rows, ignore_index=True)

    burden = (
        pt_burden.groupby(["Host", "Period", "PTs"], observed=True)
        .agg(n_isolates=("ID", "nunique"), mean_ARG_burden=("Number of ARGs", "mean"), median_ARG_burden=("Number of ARGs", "median"))
        .reset_index()
    )
    return richness, prevalence, burden


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--permutations", type=int, default=999)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    pcoa_scores, permanova, lefse = run_resistome(args.permutations)
    pt_richness, pt_prevalence, pt_burden = run_plasmid()
    output = OUT / "Derived_resistome_plasmid_reproducibility.xlsx"
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        pcoa_scores.to_excel(writer, sheet_name="PCoA_scores", index=False)
        permanova.to_excel(writer, sheet_name="PERMANOVA", index=False)
        lefse.to_excel(writer, sheet_name="LEfSe_scores_supplied", index=False)
        pt_richness.to_excel(writer, sheet_name="PT_richness", index=False)
        pt_prevalence.to_excel(writer, sheet_name="PT_prevalence", index=False)
        pt_burden.to_excel(writer, sheet_name="PT_ARG_burden", index=False)

    validation = {
        "genomes": len(pcoa_scores),
        "lefse_rows_by_host": lefse.groupby("Host").size().astype(int).to_dict(),
        "predicted_PT_richness": dict(zip(pt_richness.host, pt_richness.predicted_PT_richness.astype(int))),
        "workflow_boundary": "processed metadata, ARG calls, PCG assignments and predicted PT calls; no raw reads",
    }
    (OUT / "validation.json").write_text(json.dumps(validation, indent=2), encoding="utf-8")
    print(output)
    print(json.dumps(validation, indent=2))


if __name__ == "__main__":
    main()
