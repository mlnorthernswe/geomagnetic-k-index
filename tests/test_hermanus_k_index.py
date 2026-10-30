"""Hermanus LRNS standalone-reference and integration checks."""
import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from datetime import date
from pathlib import Path
from unittest.mock import patch

import calculate_hermanus_k_index as hermanus
import calculate_k_index as ih
import run_pipeline as pipeline

ROOT = Path(__file__).resolve().parent.parent


class TestHermanus(unittest.TestCase):
    def test_standalone_reference(self):
        fixture = json.loads((ROOT / "tests/data/hermanus-reference.json").read_text())
        self.assertEqual(hashlib.sha256((ROOT / "vendor/hermanus/K_FILTER.FOR").read_bytes()).hexdigest(),
                         fixture["source_sha256"])
        with hermanus.Hermanus() as backend:
            for case in fixture["cases"]:
                with self.subTest(path=case["path"], settings=case["settings"]):
                    path = ROOT / case["path"]
                    self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), case["input_sha256"])
                    record = ih.parse_month_file(path)[case["day_index"]]
                    for component, expected in zip(("x", "y"), case["components"]):
                        start = ih.COMPONENT_START_WORDS[component]
                        values = [v / 10 for v in record[start:start+1440]]
                        residual, smooth = backend.calculate(values, **case["settings"])
                        ranges = [max(residual[i:i+180])-min(residual[i:i+180]) for i in range(0,1440,180)]
                        for actual, want in zip(ranges, expected["ranges"]):
                            self.assertAlmostEqual(actual, want, places=5)
                        for index, (raw, quiet) in expected["checkpoints"].items():
                            i = int(index)
                            self.assertAlmostEqual(smooth[i], quiet, places=5)
                            self.assertAlmostEqual(residual[i]+smooth[i], raw, places=5)

    def test_deterministic_seed_and_constant_field(self):
        import math
        with hermanus.Hermanus() as backend:
            self.assertEqual(backend.calculate([10000.] * 1440)[0], [0.] * 1440)
            signal = [10000+30*math.sin(i/200) for i in range(1440)]
            result = backend.calculate(signal, seed=42)
            self.assertEqual(result, backend.calculate(signal, seed=42))
            self.assertNotEqual(result, backend.calculate(signal, seed=43))
            self.assertEqual(len(signal), 1440)

    def test_missing_component_day(self):
        day = date(2023, 1, 2)
        x, y = [100000] * 1440, [20000] * 1440
        rows = hermanus.calculate_days({day: (x, y, None, None)}, k9_nT=500)
        self.assertEqual([r["k"] for r in rows], [0] * 8)
        x[800] = None
        rows = hermanus.calculate_days({day: (x, y, None, None)}, k9_nT=500)
        self.assertTrue(all(r["k"] is None and r["k_x"] is None and r["k_y"] == 0 for r in rows))
        self.assertEqual(rows[0]["x_coverage"], 1)
        self.assertLess(rows[4]["x_coverage"], 1)

    def test_cli_and_pipeline(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "k.csv"
            result = subprocess.run([sys.executable, str(ROOT / "calculate_k_index.py"),
                str(ROOT / "KIR/KIR23DEC.BIN"), "--method", "hermanus", "-o", str(output)],
                text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("single-day fit", result.stderr)
            with output.open() as source:
                rows = list(csv.DictReader(source))
            self.assertEqual(len(rows), 248)
            self.assertTrue(all(r["k"] != "" for r in rows[:8]+rows[-8:]))
            summaries = pipeline.run_station_pipeline("KIR", {"directory": "KIR", "method": "hermanus", "hermanus_mu": .15},
                {}, ROOT, output_dir=Path(tmp), dry_run=True)
            self.assertEqual(len(summaries), 4)
            self.assertTrue(all(s.csv_output.name.endswith("-HERMANUS.csv") for s in summaries))
            self.assertTrue(all(s.sq_output is None for s in summaries))

    def test_iaga_minute_input_cli(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "one-day.min"
            path.write_text("DATE TIME DOY TESTX TESTY\n"+"".join(
                f"2023-01-02 {i//60:02d}:{i%60:02d}:00 002 10000 2000\n" for i in range(1440)))
            output = Path(tmp) / "output.csv"
            result = subprocess.run([sys.executable, str(ROOT / "calculate_k_index.py"), str(path),
                "--method", "hermanus", "--input-format", "iaga", "--iaga-sample-step", "60",
                "--k9-limit", "500", "-o", str(output)], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with output.open() as source:
                rows = list(csv.DictReader(source))
            self.assertEqual([r["k"] for r in rows], ["0"] * 8)


class TestHermanusValidation(unittest.TestCase):
    def test_invalid_parameters(self):
        defaults = dict(mu=.17, qindex="M", qstart=300, qend=1009, qfactor=None, seed=1)
        for override in [dict(mu=float("nan")), dict(mu=.3), dict(qindex="X"), dict(qstart=0),
                         dict(qend=299), dict(qfactor=-1), dict(seed=0), dict(seed=2147483647)]:
            with self.subTest(override=override), self.assertRaises(ValueError):
                hermanus.validate_parameters(**(defaults | override))

    def test_no_compiler_or_subprocess_required(self):
        with patch.dict("os.environ", {"FC": "/nonexistent/hermanus-compiler", "PATH": ""}), \
                patch("subprocess.run", side_effect=AssertionError("Hermanus must not invoke a compiler")):
            with hermanus.Hermanus() as backend:
                self.assertEqual(backend.calculate([10000.] * 1440)[0], [0.] * 1440)

    def test_cli_incompatible_options(self):
        for args, message in [(["--context-file", "unused.BIN"], "uses one day"),
                              (["--sq-output", "unused.txt"], "own daily quiet curve"),
                              (["--longitude", "20"], "no longitude parameter"),
                              (["--hermanus-mu", "nan"], "mu must be finite")]:
            result = subprocess.run([sys.executable, str(ROOT / "calculate_k_index.py"),
                str(ROOT / "KIR/KIR23DEC.BIN"), "--method", "hermanus", *args],
                text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(message, result.stderr)
