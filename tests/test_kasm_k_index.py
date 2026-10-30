"""Independent KASM 1.11 reference, gap, input, and integration checks."""
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from datetime import date, timedelta
from pathlib import Path
from unittest.mock import patch

import calculate_k_index as ih
import calculate_kasm_k_index as kasm
import run_pipeline as pipeline

ROOT = Path(__file__).resolve().parent.parent


@unittest.skipUnless(shutil.which("c++"), "KASM requires C++11")
class TestKASM(unittest.TestCase):
    def test_author_standalone_reference(self):
        fixture = json.loads((ROOT / "tests/data/kasm-reference.json").read_text())
        for case in fixture["cases"]:
            with self.subTest(station=case["station"], smoothing=case["smoothing"]):
                records = []
                for item in case["inputs"]:
                    path = ROOT / item["path"]
                    self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), item["sha256"])
                    records.extend(ih.parse_month_file(path))
                days = kasm.iaf_days(records)
                rows = kasm.calculate_days(days, [d for d in days if d.month == 12],
                    k9_nT=case["k9"], smoothing=case["smoothing"], minimum_coverage=.01)
                for day, expected in case["expected"].items():
                    actual = [r["k"] if r["k"] is not None else -1 for r in rows if r["date"] == day]
                    self.assertEqual(actual, expected, day)
                for row in rows:
                    if row["k"] is not None:
                        self.assertLessEqual(row["k"], max(row["k_x"], row["k_y"]))

    def test_gaps_and_boundary_invalidation(self):
        with kasm.KASM() as backend:
            x, y = [10000.] * 1800, [2000.] * 1800
            x[200:206] = [None] * 6
            result = backend.calculate(x, y, 180, 500)
            self.assertEqual(result[0], [0] * 8)
            self.assertEqual(result[5], [0] * 8)
            x[200:207] = [None] * 7
            self.assertGreater(backend.calculate(x, y, 180, 500)[5][0], 0)
            x = [None] * 1800
            self.assertTrue(all(v > 0 for v in backend.calculate(x, y, 180, 500)[5]))
            x = [10000.] * 1440
            x[0] = None
            self.assertGreater(backend.calculate(x, y[:1440], 0, 500)[5][0], 0)

    def test_context_fallback_and_coverage(self):
        day = date(2023, 1, 2)
        days = {day: ([100000] * 1440, [20000] * 1440, None, None)}
        rows = kasm.calculate_days(days, k9_nT=500)
        self.assertEqual([r["k"] for r in rows], [0] * 8)
        for offset in (-1, 1):
            days[day + timedelta(days=offset)] = ([100000] * 1440, [20000] * 1440, None, None)
        # Isolated holes are interpolated by native KASM, but observed coverage still applies.
        for i in range(20, 180, 5):
            days[day][0][i] = None
        rows = kasm.calculate_days(days, [day], k9_nT=500)
        self.assertIsNone(rows[0]["k"])
        self.assertLess(rows[0]["x_coverage"], .9)
        self.assertEqual(rows[0]["k_y"], 0)

    def test_cli_and_pipeline(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "k.csv"
            result = subprocess.run([sys.executable, str(ROOT / "calculate_k_index.py"),
                str(ROOT / "KIR/KIR23DEC.BIN"), "--method", "kasm", "--context-file",
                str(ROOT / "KIR/KIR23NOV.BIN"), "-o", str(out)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("1 day(s) calculated with reduced", result.stderr)
            with out.open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len(rows), 248)
            self.assertTrue(all(r["k"] != "" for r in rows[-8:]))
            summaries = pipeline.run_station_pipeline("KIR", {"directory": "KIR", "method": "kasm"},
                {}, ROOT, output_dir=Path(tmp), dry_run=True)
            self.assertEqual(len(summaries), 4)
            self.assertTrue(all(s.csv_output.name.endswith("-KASM.csv") for s in summaries))
            self.assertTrue(all(s.sq_output is None for s in summaries))


class TestKASMInputs(unittest.TestCase):
    def test_iaga_minute_and_second_means(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "samples.min"
            path.write_text("DATE TIME DOY TESTX TESTY\n"
                            "2023-01-02 00:00:00 002 10000.1 2000.2\n"
                            "2023-01-02 00:02:00 002 10001.1 2001.2\n")
            x, y, _, _ = kasm.iaga_days(path, sample_step=60)[date(2023, 1, 2)]
            self.assertAlmostEqual(x[0], 100001)
            self.assertIsNone(x[1])
            self.assertAlmostEqual(y[2], 20012)
            path.write_text("DATE TIME DOY TESTX TESTY\n" + "".join(
                f"2023-01-02 00:00:{i:02d} 002 {10000+i} {99999 if i == 4 else 2000}\n"
                for i in range(60)))
            x, y, _, _ = kasm.iaga_days(path)[date(2023, 1, 2)]
            self.assertEqual(x[0], 100295)
            self.assertIsNone(y[0])
            self.assertIsNone(x[1])
            with self.assertRaisesRegex(ValueError, "align"):
                kasm.iaga_days(path, sample_step=60)

    def test_missing_compiler(self):
        with patch.dict("os.environ", {"CXX": "/nonexistent/kasm-compiler"}):
            with self.assertRaisesRegex(ValueError, "C\\+\\+11 compiler"):
                with kasm.KASM():
                    pass

    def test_cli_rejects_invalid_options(self):
        for flags, message in [(["--sq-output", "/tmp/unused-sq"], "own daily quiet curve"),
                               (["--kasm-smoothing", "nan"], "finite and positive"),
                               (["--longitude", "20"], "no longitude parameter")]:
            result = subprocess.run([sys.executable, str(ROOT / "calculate_k_index.py"),
                str(ROOT / "KIR/KIR23DEC.BIN"), "--method", "kasm", *flags],
                capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(message, result.stderr)
