"""FMI K indices using the published FMI 1.3 C implementation.

No external Python dependencies. A C compiler is required for this optional
backend. The compiled library lives only in a temporary directory.
"""
from __future__ import annotations

import ctypes
import math
import shlex
import os
import subprocess
import tempfile
from datetime import datetime, timedelta
from pathlib import Path

import calculate_k_index as ih


class FMI:
    """Sequential, context-managed adapter to FMI's global-state C routines."""

    def __enter__(self):
        if ctypes.sizeof(ctypes.c_long) < 8:
            raise ValueError("FMI backend requires a 64-bit Unix C long")
        self._temp = tempfile.TemporaryDirectory(prefix="k-index-fmi-")
        source = Path(__file__).resolve().parent / "vendor" / "fmi" / "wrapper.c"
        library = Path(self._temp.name) / "fmi.so"
        try:
            result = subprocess.run(
                [*shlex.split(os.environ.get("CC", "cc")), "-std=gnu89", "-O2",
                 "-shared", "-fPIC", str(source), "-lm", "-o", str(library)],
                capture_output=True, text=True,
            )
            if result.returncode:
                raise ValueError(f"cannot compile FMI backend: {result.stderr.strip()}")
            self._lib = ctypes.CDLL(str(library))
            ptr = ctypes.POINTER(ctypes.c_long)
            self._lib.fmi_details.argtypes = [ptr, ptr, *([ctypes.c_long] * 4), *([ptr] * 5)]
            self._lib.fmi_details.restype = None
        except (OSError, ValueError) as error:
            self._temp.cleanup()
            raise ValueError(f"FMI requires a working C compiler (cc or CC): {error}") from error
        return self

    def __exit__(self, *args):
        self._temp.cleanup()

    def calculate(self, x, y, step, k9_nT, longitude):
        if step not in (1, 60):
            raise ValueError("FMI adapter supports one-second or one-minute samples")
        if len(x) != 3 * 86400 // step or len(y) != len(x):
            raise ValueError("FMI requires exactly three days of X and Y samples")
        if not math.isfinite(k9_nT) or k9_nT <= 0 or round(k9_nT * 10) < 100:
            raise ValueError("FMI K9 limit must be finite and at least 10 nT")
        if not math.isfinite(longitude) or not -180 <= longitude <= 360:
            raise ValueError("longitude must be finite and between -180 and 360 degrees")
        # IAF uses 0.1 nT. This sentinel is outside physical geomagnetic values.
        missing = 2**50
        array = ctypes.c_long * len(x)
        cx = array(*(missing if v is None else int(v) for v in x))
        cy = array(*(missing if v is None else int(v) for v in y))
        outputs = [(ctypes.c_long * 8)() for _ in range(5)]
        self._lib.fmi_details(cx, cy, step, round(k9_nT * 10),
                              round(longitude) % 360, missing, *outputs)
        return [list(a) for a in outputs]


def iaf_days(records):
    """Map dates to X/Y arrays in 0.1 nT, stored record, and header longitude."""
    days = {}
    identity = None
    for record in records:
        station = (record[0], record[2], record[3])
        if identity is not None and station != identity:
            raise ValueError("FMI IAF inputs must belong to the same station and coordinates")
        identity = station
        day = datetime.strptime(str(record[1]), "%Y%j").date()
        if day in days:
            raise ValueError(f"duplicate FMI input date: {day}")
        components = [
            [None if v == ih.MISSING_VALUE else v for v in record[start:start + 1440]]
            for start in ih.COMPONENT_START_WORDS.values()
        ]
        days[day] = (*components, record, record[3] / 1000.0)
    return days


def iaga_days(path):
    """Retain one-second samples; absent timestamps become explicit gaps."""
    days = {}
    for day, samples in ih.iter_iaga_days(path):
        x, y = [None] * 86400, [None] * 86400
        for second, vx, vy in samples:
            if second != int(second):
                raise ValueError("FMI IAGA input requires whole-second timestamps")
            i = int(second)
            x[i] = None if vx is None else round(vx * 10)
            y[i] = None if vy is None else round(vy * 10)
        days[day] = (x, y, None, None)
    return days


def calculate_days(days, target_days=None, *, step=60, k9_nT=None,
                   longitude=None, minimum_coverage=0.9):
    """Return existing CSV-schema rows. Missing adjacent days invalidate a day.

    Never invent context by repeating or extrapolating boundary days. Coverage
    reports observed samples before FMI's short-gap interpolation.
    """
    if not 0 < minimum_coverage <= 1:
        raise ValueError("minimum coverage must be greater than 0 and at most 1")
    rows = []
    count = 86400 // step
    with FMI() as backend:
        for day in sorted(days if target_days is None else target_days):
            x, y, record, header_longitude = days[day]
            limit = k9_nT if k9_nT is not None else (record[10] if record else None)
            lon = longitude if longitude is not None else header_longitude
            if limit is None or lon is None:
                raise ValueError("FMI requires a K9 limit and observatory longitude")
            if not math.isfinite(limit) or limit < 10 or not math.isfinite(lon) or not -180 <= lon <= 360:
                raise ValueError("FMI requires K9 >= 10 nT and longitude in [-180, 360]")
            neighbors = [day + timedelta(days=offset) for offset in (-1, 0, 1)]
            has_context = all(d in days for d in neighbors)
            if has_context:
                result = backend.calculate(
                    [v for d in neighbors for v in days[d][0]],
                    [v for d in neighbors for v in days[d][1]], step, limit, lon,
                )
                ks, kx, ky, rx, ry = result
            for i, (start, end) in enumerate(ih.INTERVALS):
                a, b = i * count // 8, (i + 1) * count // 8
                cx = sum(v is not None for v in x[a:b]) / (b - a)
                cy = sum(v is not None for v in y[a:b]) / (b - a)
                ix = kx[i] if has_context and kx[i] >= 0 and cx >= minimum_coverage else None
                iy = ky[i] if has_context and ky[i] >= 0 and cy >= minimum_coverage else None
                # Match the existing application's requirement for both components.
                k = max(ix, iy) if ix is not None and iy is not None else None
                stored_raw = record[ih.K_START_WORD + i] if record else None
                stored = None if stored_raw in (None, ih.MISSING_K_VALUE) else stored_raw / 10
                rows.append(dict(
                    date=day.isoformat(), interval_utc=f"{start:02d}:00-{end:02d}:00",
                    x_range_nT=rx[i] / 10 if ix is not None else None,
                    y_range_nT=ry[i] / 10 if iy is not None else None,
                    x_coverage=cx, y_coverage=cy, k_x=ix, k_y=iy, k=k,
                    stored_k_raw=stored_raw, stored_k=stored,
                    matches_stored="" if k is None or stored is None else int(k == stored),
                ))
    return rows
