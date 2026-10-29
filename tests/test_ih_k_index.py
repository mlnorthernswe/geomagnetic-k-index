#!/usr/bin/env python3
"""Automated unit and regression tests for calculate_ih_k_index and run_pipeline."""

from __future__ import annotations

import io
import tempfile
import unittest
from pathlib import Path

import calculate_ih_k_index as ih
import run_pipeline as pipeline


REPO_ROOT = Path(__file__).resolve().parent.parent


class TestKIndexUnit(unittest.TestCase):
    def test_classify_range_pt_boundaries(self) -> None:
        # Kiruna K9 = 1500 nT -> thresholds: 15, 30, 60, 120, 210, 360, 600, 990, 1500 nT
        k9 = 1500
        self.assertEqual(ih.classify_range_pt(0.0, k9), 0)
        self.assertEqual(ih.classify_range_pt(14_999.9, k9), 0)
        self.assertEqual(ih.classify_range_pt(15_000.0, k9), 1)
        self.assertEqual(ih.classify_range_pt(30_000.0, k9), 2)
        self.assertEqual(ih.classify_range_pt(990_000.0, k9), 8)
        self.assertEqual(ih.classify_range_pt(1_500_000.0, k9), 8)
        self.assertEqual(ih.classify_range_pt(1_500_000.1, k9), 9)

    def test_interpolate_sq_pt_wraparound(self) -> None:
        hourly = tuple(float(h * 1000) for h in range(24))
        self.assertAlmostEqual(ih.interpolate_sq_pt(hourly, 0.0), 0.0)
        self.assertAlmostEqual(ih.interpolate_sq_pt(hourly, 1800.0), 500.0)
        # Hour 23 (23000) interpolates halfway to hour 0 (0) -> 11500
        self.assertAlmostEqual(ih.interpolate_sq_pt(hourly, 23 * 3600 + 1800), 11500.0)

    def test_sq_read_write_roundtrip(self) -> None:
        sq_path = REPO_ROOT / "KIR" / "202310KIRsq"
        sq = ih.read_sq_curve(sq_path)
        self.assertEqual(set(sq.keys()), {"x", "y", "z"})
        self.assertEqual(len(sq["x"]), 24)

        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "roundtrip_sq.txt"
            ih.write_sq_curve(sq, out_path)
            reloaded = ih.read_sq_curve(out_path)
            self.assertEqual(sq, reloaded)

    def test_bootstrap_median_vs_mean(self) -> None:
        records = ih.parse_month_file(REPO_ROOT / "KIR" / "KIR23DEC.BIN")
        sq_median, days_median = ih.bootstrap_sq_curve(records, quiet_day_count=10, estimator="median")
        sq_mean, days_mean = ih.bootstrap_sq_curve(records, quiet_day_count=10, estimator="mean")

        self.assertEqual(days_median, days_mean)
        self.assertEqual(len(days_median), 10)
        self.assertNotEqual(sq_median["x"], sq_mean["x"])


class TestRegressionAndPipeline(unittest.TestCase):
    def test_wic_dec_2023_regression(self) -> None:
        nov_records = ih.parse_month_file(REPO_ROOT / "WIC" / "WIC23NOV.BIN")
        dec_records = ih.parse_month_file(REPO_ROOT / "WIC" / "WIC23DEC.BIN")
        sq, _ = ih.bootstrap_sq_curve(nov_records, quiet_day_count=10, estimator="median")
        rows = [row for record in dec_records for row in ih.calculate_day(record, sq, minimum_coverage=0.9)]

        comp_buffer = io.StringIO()
        ih.write_comprehensive(rows, comp_buffer)
        expected_comp = (REPO_ROOT / "WIC" / "WIC23DEC-comprehensive.txt").read_text(encoding="ascii")
        self.assertEqual(comp_buffer.getvalue(), expected_comp)

        csv_buffer = io.StringIO()
        ih.write_csv(rows, csv_buffer)
        lines = csv_buffer.getvalue().splitlines()
        self.assertEqual(lines[0].split(",")[8], "k")
        self.assertEqual(len(lines), 249)

    def test_pipeline_execution_and_chaining(self) -> None:
        defaults, stations_cfg = pipeline.load_config(REPO_ROOT / "stations.json")
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)
            summaries = pipeline.run_station_pipeline(
                station="BEL",
                station_cfg=stations_cfg["BEL"],
                defaults=defaults,
                base_dir=REPO_ROOT,
                output_dir=out_dir,
                save_generated_sq=True,
            )
            self.assertEqual(len(summaries), 4)
            # First month bootstraps and saves 202311BELsq; second month loads it
            self.assertTrue(summaries[0].sq_source.startswith("bootstrap:median(202311"))
            self.assertEqual(summaries[1].sq_source, "file:202311BELsq")
            self.assertTrue(summaries[0].csv_output.is_file())
            self.assertTrue(summaries[0].comprehensive_output.is_file())
            self.assertTrue((out_dir / "BEL" / "202311BELsq").is_file())


if __name__ == "__main__":
    unittest.main()
