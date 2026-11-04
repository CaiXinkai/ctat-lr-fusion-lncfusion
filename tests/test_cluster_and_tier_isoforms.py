import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class ClusterAndTierTest(unittest.TestCase):
    def test_terminal_subsets_and_strands_collapse_but_raw_is_preserved(self):
        root = Path(__file__).resolve().parents[1]
        script = root / "scripts" / "cluster_and_tier_isoforms.py"
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            entities = tmp / "entities.tsv"
            sequences = tmp / "sequences.fa.gz"
            final = tmp / "final.tsv"
            out = tmp / "tiers"
            fields = [
                "feature_id", "sample_id", "fusion_name", "left_gene", "right_gene",
                "left_breakpoint", "right_breakpoint", "splice_type",
                "representative_block_chain", "representative_strand",
                "left_local_breakpoint", "right_local_breakpoint", "local_breakpoint_source",
                "sample_count", "total_support_reads", "support_read_ids", "ctat_final_row_numbers",
            ]
            rows = [
                ["f1", "S1", "A--B", "A", "B", "chr1:100:+", "chr2:200:+", "ONLY_REF_SPLICE",
                 "A--B:1-100:+;A--B:200-250:+;A--B:300-400:+", "+", "100", "200", "TEST", "1", "1", "r1", "1"],
                ["f2", "S1", "A--B", "A", "B", "chr1:100:+", "chr2:200:+", "ONLY_REF_SPLICE",
                 "A--B:20-100:-;A--B:200-250:-;A--B:300-380:-", "-", "100", "200", "TEST", "1", "1", "r2", "1"],
                ["f3", "S1", "A--B", "A", "B", "chr1:100:+", "chr2:200:+", "ONLY_REF_SPLICE",
                 "A--B:25-100:+;A--B:200-250:+", "+", "100", "200", "TEST", "1", "1", "r3", "1"],
            ]
            with entities.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(fields)
                writer.writerows(rows)
            with gzip.open(sequences, "wt", encoding="utf-8") as handle:
                for feature in ("f1", "f2", "f3"):
                    handle.write(f">{feature}\nACGT\n")
            final.write_text(
                "#FusionName\tLeftGene\tRightGene\tLeftBreakpoint\tRightBreakpoint\tSpliceType\tnum_LR\tLR_FFPM\tfrac_dom_iso\n"
                "A--B\tA\tB\tchr1:100:+\tchr2:200:+\tONLY_REF_SPLICE\t3\t1.0\t1.0\n",
                encoding="utf-8",
            )
            subprocess.run([
                sys.executable, str(script), "--entities", str(entities),
                "--sequences", str(sequences), "--final-tsv", str(final),
                "--sample", "S1", "--output-dir", str(out),
            ], check=True, capture_output=True, text=True)
            with (out / "S1.isoform_tier_summary.tsv").open(encoding="utf-8", newline="") as handle:
                summary = {row["tier"]: int(row["isoform_count"])
                           for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(summary["00_raw_exact_entities"], 3)
            self.assertEqual(summary["01_clustered_all"], 1)
            self.assertEqual(summary["04_pancancer_hc3"], 1)


if __name__ == "__main__":
    unittest.main()
