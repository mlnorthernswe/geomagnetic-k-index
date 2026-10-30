"""KASM 1.11 adapter using the author's unmodified C++ smoothing routines."""
from __future__ import annotations

import ctypes
import math
import os
import shlex
import subprocess
import tempfile
import struct
from datetime import timedelta
from pathlib import Path

import calculate_k_index as ih
from calculate_fmi_k_index import iaf_days as _iaf_days


class KASM:
    """Compile the optional backend once per run; accept minute means in nT."""

    def __enter__(self):
        self._temp = tempfile.TemporaryDirectory(prefix="k-index-kasm-")
        source = Path(__file__).resolve().parent / "vendor" / "kasm"
        library = Path(self._temp.name) / "kasm.so"
        try:
            result = subprocess.run(
                [*shlex.split(os.environ.get("CXX", "c++")), "-std=c++11", "-O2",
                 "-shared", "-fPIC", str(source / "casm.cpp"),
                 str(source / "wrapper.cpp"), "-o", str(library)],
                capture_output=True, text=True,
            )
            if result.returncode:
                raise ValueError(result.stderr.strip())
            self._lib = ctypes.CDLL(str(library))
            dp = ctypes.POINTER(ctypes.c_double)
            ip = ctypes.POINTER(ctypes.c_int)
            self._lib.kasm_details.argtypes = [ctypes.c_int, ctypes.c_int, dp, dp,
                ctypes.c_double, ctypes.c_double, ip, ip, ip, dp, dp, ip, ip]
            self._lib.kasm_details.restype = ctypes.c_int
        except (OSError, ValueError) as error:
            self._temp.cleanup()
            raise ValueError(f"KASM requires a working C++11 compiler (c++ or CXX): {error}") from error
        return self

    def __exit__(self, *args):
        self._temp.cleanup()

    def calculate(self, x, y, offset, k9_nT, smoothing=1.0):
        if len(x) not in (1440, 1620, 1800) or len(y) != len(x):
            raise ValueError("KASM requires 1440 minute samples plus optional 180-minute wings")
        if offset not in (0, 180) or offset + 1440 > len(x):
            raise ValueError("invalid KASM target-day offset")
        if not math.isfinite(k9_nT) or k9_nT <= 0:
            raise ValueError("KASM K9 limit must be finite and positive")
        if not math.isfinite(smoothing) or smoothing <= 0:
            raise ValueError("KASM smoothing must be finite and positive")
        for v in (*x, *y):
            if v is not None and (not math.isfinite(v) or abs(v) >= 1e29):
                raise ValueError("invalid KASM magnetic field value")
        array = ctypes.c_double * len(x)
        cx = array(*(1e30 if v is None else v for v in x))
        cy = array(*(1e30 if v is None else v for v in y))
        ks, kx, ky, invalid_x, invalid_y = [(ctypes.c_int * 8)() for _ in range(5)]
        rx, ry = (ctypes.c_double * 8)(), (ctypes.c_double * 8)()
        status = self._lib.kasm_details(len(x), offset, cx, cy, k9_nT, smoothing,
                                       ks, kx, ky, rx, ry, invalid_x, invalid_y)
        if status:
            raise ValueError(f"KASM numerical calculation failed (status {status})")
        return [list(a) for a in (ks, kx, ky, rx, ry, invalid_x, invalid_y)]


def iaf_days(records):
    """Use the shared X/Y IAF layout and 0.1 nT units."""
    for record in records:
        orientation = struct.pack("<i", record[5]).decode("ascii", errors="replace").strip()
        if not orientation.startswith("XY"):
            raise ValueError("KASM IAF input must contain X/Y horizontal components")
    try:
        return _iaf_days(records)
    except ValueError as error:
        raise ValueError(str(error).replace("FMI", "KASM")) from error


def iaga_days(path, sample_step=1):
    """Read one-second or minute IAGA X/Y and return minute means in 0.1 nT.

    A second-resolution minute needs all 60 valid samples for that component;
    missing timestamps are not silently discarded from the mean.
    """
    if sample_step not in (1, 60):
        raise ValueError("KASM IAGA sample step must be 1 or 60 seconds")
    days = {}
    for day, samples in ih.iter_iaga_days(path):
        sums = [[0.0] * 1440 for _ in range(2)]
        counts = [[0] * 1440 for _ in range(2)]
        for second, vx, vy in samples:
            if second % sample_step:
                raise ValueError("IAGA timestamps must align with --iaga-sample-step")
            minute = int(second) // 60
            for c, value in enumerate((vx, vy)):
                if value is not None:
                    sums[c][minute] += value
                    counts[c][minute] += 1
        needed = 60 // sample_step
        components = [[sums[c][i] * 10 / needed if counts[c][i] == needed else None
                       for i in range(1440)] for c in range(2)]
        days[day] = (*components, None, None)
    return days


def calculate_days(days, target_days=None, *, k9_nT=None, smoothing=1.0,
                   minimum_coverage=0.9):
    """Use each available three-hour wing; fall back to current day at boundaries."""
    if not 0 < minimum_coverage <= 1:
        raise ValueError("minimum coverage must be greater than 0 and at most 1")
    rows = []
    with KASM() as backend:
        for day in sorted(days if target_days is None else target_days):
            x, y, record, _ = days[day]
            limit = k9_nT if k9_nT is not None else (record[10] if record else None)
            if limit is None:
                raise ValueError("KASM requires a K9 limit")
            previous, following = day - timedelta(days=1), day + timedelta(days=1)
            offset = 180 if previous in days else 0
            inputs = []
            for c in (0, 1):
                values = ((days[previous][c][-180:] if previous in days else [])
                          + days[day][c]
                          + (days[following][c][:180] if following in days else []))
                inputs.append([None if v is None else v / 10 for v in values])
            ks, kx, ky, rx, ry, invalid_x, invalid_y = backend.calculate(
                *inputs, offset, limit, smoothing)
            for i, (start, end) in enumerate(ih.INTERVALS):
                a, b = i * 180, (i + 1) * 180
                cx = sum(v is not None for v in x[a:b]) / 180
                cy = sum(v is not None for v in y[a:b]) / 180
                ix = kx[i] if not invalid_x[i] and cx >= minimum_coverage else None
                iy = ky[i] if not invalid_y[i] and cy >= minimum_coverage else None
                k = ks[i] if ix is not None and iy is not None else None
                stored_raw = record[ih.K_START_WORD + i] if record else None
                stored = None if stored_raw in (None, ih.MISSING_K_VALUE) else stored_raw / 10
                rows.append(dict(
                    date=day.isoformat(), interval_utc=f"{start:02d}:00-{end:02d}:00",
                    x_range_nT=rx[i] if ix is not None else None,
                    y_range_nT=ry[i] if iy is not None else None,
                    x_coverage=cx, y_coverage=cy, k_x=ix, k_y=iy, k=k,
                    stored_k_raw=stored_raw, stored_k=stored,
                    matches_stored="" if k is None or stored is None else int(k == stored),
                ))
    return rows
