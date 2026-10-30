"""Validate the optional FMI backend against FMI's published reference data."""
import csv
import io
import shutil
import subprocess
import sys
import tempfile
import unittest
from datetime import date, timedelta
from pathlib import Path

import calculate_fmi_k_index as fmi
import calculate_k_index as ih
import run_pipeline as pipeline

ROOT = Path(__file__).resolve().parent.parent


@unittest.skipUnless(shutil.which("cc"), "FMI backend requires a C compiler")
class TestFMI(unittest.TestCase):
    def test_published_fmi_example(self):
        pairs = [tuple(map(int, line.split())) for line in
                 (ROOT / "tests/data/fmi-demo.dat").read_text().splitlines()]
        x, y = [p[0] for p in pairs], [p[1] for p in pairs]
        with fmi.FMI() as backend:
            ks, kx, ky, rx, ry = backend.calculate(x, y, 60, 750, 27)
            self.assertEqual(ks, [1, 1, 1, 2, 3, 3, 5, 3])
            self.assertEqual(ks, [max(a, b) for a, b in zip(kx, ky)])
            self.assertTrue(all(v >= 0 for v in rx + ry))
            # Calls do not mutate caller arrays, and global state is reset.
            self.assertEqual(backend.calculate(x, y, 60, 750, 27)[0], ks)
        self.assertEqual(x, [p[0] for p in pairs])

    def test_gap_interpolation_and_large_gaps(self):
        x, y = [100000] * 4320, [20000] * 4320
        with fmi.FMI() as backend:
            x[1500:1514] = [None] * 14
            self.assertEqual(backend.calculate(x, y, 60, 750, -75)[1][0], 0)
            x[1500:1515] = [None] * 15
            self.assertEqual(backend.calculate(x, y, 60, 750, -75)[1][0], -1)
            self.assertEqual(backend.calculate(x, y, 60, 750, 285)[0],
                             backend.calculate(x, y, 60, 750, -75)[0])

    def test_context_and_component_coverage(self):
        day = date(2023, 1, 2)
        def data():
            return ([100000] * 1440, [20000] * 1440, None, 20)
        days = {day + timedelta(days=n): data() for n in (-1, 0, 1)}
        rows = fmi.calculate_days(days, k9_nT=1500)
        self.assertTrue(all(r["k"] is None for r in rows[:8] + rows[-8:]))
        self.assertTrue(all(r["k"] == 0 for r in rows[8:16]))
        days[day][0][100:120] = [None] * 20
        rows = fmi.calculate_days(days, [day], k9_nT=1500)
        self.assertIsNone(rows[0]["k"])
        self.assertIsNone(rows[0]["k_x"])
        self.assertEqual(rows[0]["k_y"], 0)
        self.assertAlmostEqual(rows[0]["x_coverage"], 160 / 180)

    def test_cli_and_pipeline_month_context(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "fmi.csv"
            result = subprocess.run([
                sys.executable, str(ROOT / "calculate_k_index.py"),
                str(ROOT / "KIR/KIR23DEC.BIN"), "--method", "fmi",
                "--context-file", str(ROOT / "KIR/KIR23NOV.BIN"), "-o", str(output),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("1 day(s) unavailable", result.stderr)
            with output.open() as source:
                rows = list(csv.DictReader(source))
            self.assertEqual(len(rows), 31 * 8)
            self.assertTrue(all(r["k"] != "" for r in rows[:8]))
            self.assertTrue(all(r["k"] == "" for r in rows[-8:]))
            summaries = pipeline.run_station_pipeline("KIR", {"directory": "KIR"}, {},
                ROOT, output_dir=Path(tmp), method_override="fmi", dry_run=True)
            self.assertTrue(all(s.csv_output.name.endswith("-FMI.csv") for s in summaries))
            self.assertTrue(all(s.sq_output is None for s in summaries))
            december = next(s for s in summaries if s.yyyymm == "202312")
            self.assertEqual(december.comparable_intervals, 30 * 8)

    def test_iaga_second_resolution(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "samples.sec"
            path.write_text("DATE TIME DOY TESTX TESTY\n"
                            "2023-01-02 00:00:00 002 10000.1 2000.2\n"
                            "2023-01-02 00:00:02 002 10001.1 2001.2\n")
            days = fmi.iaga_days(path)
            x, y, _, _ = days[date(2023, 1, 2)]
            self.assertEqual(len(x), 86400)
            self.assertEqual(x[:3], [100001, None, 100011])
            self.assertEqual(y[:3], [20002, None, 20012])
        with fmi.FMI() as backend:
            self.assertEqual(backend.calculate([100001] * 259200, [20002] * 259200,
                                               1, 750, 27)[0], [0] * 8)

    def test_reject_mixed_stations(self):
        records = ih.parse_month_file(ROOT / "KIR/KIR23DEC.BIN")[:2]
        other = list(records[1])
        other[0] += 1
        with self.assertRaisesRegex(ValueError, "same station"):
            fmi.iaf_days([records[0], tuple(other)])


class TestFMIValidation(unittest.TestCase):
    def test_fmi_rejects_sq_options(self):
        result = subprocess.run([sys.executable, str(ROOT / "calculate_k_index.py"),
            str(ROOT / "KIR/KIR23DEC.BIN"), "--method", "fmi", "--sq-output", "/tmp/unused-sq"],
            capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("FMI fits its own daily quiet curve", result.stderr)

    def test_missing_compiler_is_actionable(self):
        from unittest.mock import patch
        with patch.dict("os.environ", {"CC": "/nonexistent/fmi-compiler"}):
            with self.assertRaisesRegex(ValueError, "working C compiler"):
                with fmi.FMI():
                    pass
