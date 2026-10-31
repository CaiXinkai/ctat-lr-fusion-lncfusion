import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class PrepareEntitiesTest(unittest.TestCase):
    def test_only_final_calls_become_entities(self):
        root = Path(__file__).resolve().parents[1]
        script = root / "scripts" / "prepare_ctat_isoform_entities.py"
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            final = tmp / "final.tsv"
            preliminary = tmp / "preliminary.tsv"
            gff3 = tmp / "fusion.gff3"
            output = tmp / "entities.tsv"
            final.write_text(
                "#FusionName\tLeftGene\tRightGene\tLeftBreakpoint\tRightBreakpoint\tSpliceType\tLR_accessions\n"
                "A--B\tA\tB\tchr1:100:+\tchr2:200:-\tONLY_REF_SPLICE\tread1\n",
                encoding="utf-8",
            )
            preliminary.write_text(
                "#FusionName\tLeftBreakpoint\tRightBreakpoint\tLeftLocalBreakpoint\tRightLocalBreakpoint\n"
                "A--B\tchr1:100:+\tchr2:200:-\t100\t101\n"
                "C--D\tchr3:300:+\tchr4:400:-\t100\t101\n",
                encoding="utf-8",
            )
            gff3.write_text(
                "##gff-version 3\n"
                "A--B\tmm2\tcDNA_match\t1\t100\t.\t+\t.\tTarget=read1 1 100\n"
                "A--B\tmm2\tcDNA_match\t101\t200\t.\t+\t.\tTarget=read1 101 200\n",
                encoding="utf-8",
            )
            subprocess.run(
                [
                    sys.executable, str(script), "--final-tsv", str(final),
                    "--preliminary-tsv", str(preliminary), "--fusion-gff3", str(gff3),
                    "--sample", "S1", "--output", str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["fusion_name"], "A--B")
            self.assertEqual(rows[0]["total_support_reads"], "1")
            self.assertEqual(rows[0]["left_local_breakpoint"], "100")
            self.assertNotIn("C--D", output.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
