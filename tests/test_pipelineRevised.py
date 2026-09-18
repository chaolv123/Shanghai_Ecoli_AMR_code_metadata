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


if __name__ == "__main__":
    unittest.main()
