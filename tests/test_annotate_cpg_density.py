#!/usr/bin/env python3

import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIRECTORY = Path(__file__).resolve().parent / "fixtures"


class AnnotateCpGDensityIntegrationTest(unittest.TestCase):
    def test_expected_metrics_and_output_names(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work_directory = Path(temporary_directory)
            genome_fasta = work_directory / "genome.fa"
            bed_file = work_directory / "intervals.bed"

            shutil.copy(FIXTURE_DIRECTORY / "genome.fa", genome_fasta)
            shutil.copy(
                FIXTURE_DIRECTORY / "genome.fa.fai",
                Path(str(genome_fasta) + ".fai"),
            )
            shutil.copy(FIXTURE_DIRECTORY / "intervals.bed", bed_file)

            completed = subprocess.run(
                [
                    sys.executable,
                    str(REPOSITORY_ROOT / "annotate_cpg_density.py"),
                    str(genome_fasta),
                    str(bed_file),
                ],
                cwd=work_directory,
                check=False,
                capture_output=True,
                text=True,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)

            expected_message = (
                "Processed {}: outputs -> intervals_CpG-per-base.bed, "
                "intervals_CpG-per-GC.bed, intervals_CpG-OE.bed\n"
            ).format(bed_file)
            self.assertEqual(completed.stdout, expected_message)
            self.assertEqual(completed.stderr, "")

            self.assertEqual(
                (work_directory / "intervals_CpG-per-base.bed").read_text(),
                "chrTest\t0\t8\t0.250000\n"
                "chrTest\t8\t14\t0.166667\n"
                "chrTest\t14\t20\t0.166667\n"
                "chrTest\t0\t1\t0.000000\n",
            )
            self.assertEqual(
                (work_directory / "intervals_CpG-per-GC.bed").read_text(),
                "chrTest\t0\t8\t0.500000\n"
                "chrTest\t8\t14\t0.250000\n"
                "chrTest\t14\t20\t0.500000\n"
                "chrTest\t0\t1\tNA\n",
            )
            self.assertEqual(
                (work_directory / "intervals_CpG-OE.bed").read_text(),
                "chrTest\t0\t8\t4.000000\n"
                "chrTest\t8\t14\t1.500000\n"
                "chrTest\t14\t20\t6.000000\n"
                "chrTest\t0\t1\tNA\n",
            )


if __name__ == "__main__":
    unittest.main()
