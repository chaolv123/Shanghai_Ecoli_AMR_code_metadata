import json
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


class MetadataPipelineTest(unittest.TestCase):
    def test_validation(self):
        result = json.loads((ROOT / "outputs" / "validation.json").read_text(encoding="utf-8"))
        self.assertEqual(result["phenotype_n"], 5299)
        self.assertEqual(result["genomic_n"], 2622)
        self.assertEqual(result["wgs_linked_n"], 2547)
        self.assertEqual(result["ARG_matrix_count_mismatch"], 0)
        self.assertEqual(result["tested_comparisons_per_host"], {"Poultry": 12, "Swine": 12})
        self.assertEqual(result["class_length_rule_discrepancies"], 0)
        self.assertEqual(result["MDR_rule_discrepancies"], 0)

    def test_table1(self):
        table = pd.read_excel(ROOT / "outputs" / "supplementary_table_1_recomputed.xlsx")
        self.assertEqual(len(table), 30)
        limited = table[table.Antimicrobial.isin(["CAZ", "MEM", "DOX"])]
        self.assertTrue(limited.P_value.isna().all())
        ff = table[(table.Host == "Swine") & (table.Antimicrobial == "FF")].iloc[0]
        self.assertAlmostEqual(ff.q_value_BH, 0.8252045, places=5)
        self.assertEqual(table.groupby("Host").q_value_BH.apply(lambda x: int((x < .05).sum())).to_dict(),
                         {"Poultry": 12, "Swine": 11})

    def test_reconciled_pcg(self):
        xlsx = ROOT / "outputs" / "extended" / "pcg_reconciled" / "Supplementary_PCG_threshold_reconciled.xlsx"
        summary = pd.read_excel(xlsx, sheet_name="Threshold_summary")
        primary = summary[summary.threshold.eq(10)].iloc[0]
        self.assertEqual((primary.n_PCGs, primary.localized_PCGs, primary.nonlocalized_PCGs), (494, 308, 186))
        self.assertEqual(primary.n_clusters_exceed_pairwise_cutoff, 0)
        fig5 = ROOT / "data" / "publication_inputs" / "Source_Data_Fig_05.xlsx"
        local = pd.read_excel(fig5, sheet_name="a").dropna(subset=["ID"])
        nonlocal_ = pd.read_excel(fig5, sheet_name="d&f")
        self.assertEqual(local.PCGs.nunique(), 308)
        self.assertEqual(nonlocal_.PCGs.nunique(), 186)
        self.assertFalse(set(local.PCGs) & set(nonlocal_.PCGs))
        self.assertIn("JEC10_013", set(local.ID))

    def test_pcg_sampling_intensity_output(self):
        xlsx = (
            ROOT
            / "outputs"
            / "extended"
            / "pcg_sampling_intensity"
            / "Supplementary_PCG_sampling_intensity_normalized_analysis.xlsx"
        )
        validation = pd.read_excel(xlsx, sheet_name="Validation")
        values = dict(zip(validation["check"], validation["value"]))
        self.assertEqual(values["Total PCG-analyzable isolates"], 2610)
        self.assertEqual(values["Distinct PCGs"], 494)
        self.assertEqual(values["Localized PCGs"], 308)
        denominators = pd.read_excel(xlsx, sheet_name="Period_denominators")
        self.assertEqual(int(denominators["n_PCG_analyzable_WGS_isolates"].sum()), 2610)
        self.assertEqual(int(denominators["n_isolates_in_localized_PCG"].sum()), 830)
        self.assertEqual(int(denominators["n_isolates_in_nonlocalized_PCG"].sum()), 629)

    def test_derived_reproducibility_output(self):
        xlsx = (
            ROOT
            / "outputs"
            / "derived_reproducibility"
            / "Derived_resistome_plasmid_reproducibility.xlsx"
        )
        richness = pd.read_excel(xlsx, sheet_name="PT_richness")
        values = dict(zip(richness["host"], richness["predicted_PT_richness"]))
        self.assertEqual(values, {"Poultry": 67, "Swine": 108, "Shared": 50})
        lefse = pd.read_excel(xlsx, sheet_name="LEfSe_scores_supplied")
        self.assertEqual(lefse.groupby("Host").size().to_dict(), {"Poultry": 41, "Swine": 44})

    def test_publication_sync(self):
        sync = ROOT / "outputs" / "publication_sync"
        manifest = pd.read_csv(sync / "publication_sync_manifest.csv")
        self.assertEqual(len(manifest), 7)
        self.assertTrue((manifest["bytes"] > 0).all())
        status = json.loads((sync / "publication_sync_validation.json").read_text(encoding="utf-8"))
        self.assertEqual(status["Supplementary_Table_09"]["annual_positive_FDR_rows"], 8)
        self.assertEqual(status["Supplementary_Table_09"]["annual_host_specific_rows"], 4)


if __name__ == "__main__":
    unittest.main()
