"""Regression checks against native FMI/C and KASM/C++ diagnostic fixtures."""
import hashlib
import json
import unittest
from pathlib import Path
from unittest.mock import patch

import calculate_fmi_k_index as fmi
import calculate_kasm_k_index as kasm
from native_reference_inputs import reference_inputs

ROOT = Path(__file__).resolve().parent.parent


class TestPythonPorts(unittest.TestCase):
    def test_native_diagnostics(self):
        fixture = json.loads((ROOT / "tests/data/native-port-reference.json").read_text())
        for case in fixture:
            with self.subTest(method=case["method"], path=case.get("path"),
                              day=case.get("day_index"), variant=case.get("variant")):
                if "path" in case:
                    self.assertEqual(hashlib.sha256((ROOT / case["path"]).read_bytes()).hexdigest(), case["sha256"])
                backend = fmi.FMI() if case["method"] == "fmi" else kasm.KASM()
                actual = backend.calculate(*reference_inputs(case), *case["args"])
                expected = case["expected"]
                for i, (values, wanted) in enumerate(zip(actual, expected)):
                    if case["method"] == "kasm" and i in (3,4):
                        # Different libm/platform rounding may affect final double bits.
                        for value, want in zip(values, wanted):
                            self.assertAlmostEqual(value, want, places=8)
                    else:
                        self.assertEqual(values, wanted)

    def test_nonconstant_inputs_without_native_tools(self):
        # Exercise smoothing and gap filling, not just the constant-input path.
        with patch.dict("os.environ", {"PATH": "", "CC": "/missing", "CXX": "/missing"}), \
                patch("subprocess.run", side_effect=AssertionError("compiler invocation")), \
                patch("ctypes.CDLL", side_effect=AssertionError("native library loading")):
            for method, args in [("fmi", [60,750,-75]), ("kasm", [180,750,1.0])]:
                case = dict(method=method, variant=2, args=args)
                backend = fmi.FMI() if method == "fmi" else kasm.KASM()
                result = backend.calculate(*reference_inputs(case), *args)
                self.assertEqual(len(result[0]), 8)
                self.assertTrue(any(k > 0 for k in result[0]))
