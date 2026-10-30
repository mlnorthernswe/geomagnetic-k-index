"""Standard-library Python port of ISGI's public-domain Hermanus LRNS routine."""
from __future__ import annotations

import math
import struct

import calculate_k_index as ih
from calculate_kasm_k_index import iaf_days as _iaf_days, iaga_days as _iaga_days


def validate_parameters(mu, qindex, qstart, qend, qfactor, seed):
    if not math.isfinite(mu) or not .08 <= mu <= .20:
        raise ValueError("Hermanus mu must be finite and between 0.08 and 0.20")
    if qindex not in ("M", "E"):
        raise ValueError("Hermanus q-index must be M or E")
    if not (isinstance(qstart, int) and isinstance(qend, int) and 1 <= qstart < qend <= 1440):
        raise ValueError("Hermanus quiet-window minutes must satisfy 1 <= q-start < q-end <= 1440")
    if qfactor is not None and (not math.isfinite(qfactor) or qfactor <= 0):
        raise ValueError("Hermanus q-factor must be finite and positive")
    if not isinstance(seed, int) or not 1 <= seed <= 2147483646:
        raise ValueError("Hermanus seed must be an integer in [1, 2147483646]")


_FLOAT32 = struct.Struct("=f")


def _f32(value):
    """Round each REAL*4 operation as in the original Fortran routine."""
    return _FLOAT32.unpack(_FLOAT32.pack(value))[0]


class Hermanus:
    """Pure Python LRNS filter; all guide-noise state is local to each call.

    The context-manager interface is retained for existing callers. No native
    library, compiler, subprocess, or temporary file is needed.
    """

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def calculate(self, signal, *, mu=.17, qindex="M", qstart=300, qend=1009,
                  qfactor=3.0, seed=1):
        validate_parameters(mu, qindex, qstart, qend, qfactor, seed)
        if qfactor is None:
            raise ValueError("Hermanus component filter requires a positive q-factor")
        if len(signal) != 1440 or any(v is None or not math.isfinite(v) or abs(v) > 1e7 for v in signal):
            raise ValueError("Hermanus requires 1440 complete, finite minute means in nT")
        f = _f32
        u = f(mu)
        # qfactor only participates in a comparison, so retain finite positive
        # values beyond float32's range as infinity, as the old C ABI did.
        factor = f(qfactor) if qfactor <= 3.4028234663852886e38 else math.inf
        data = [f(v) for v in signal]
        guide = []
        state = seed
        noise_scale = f(.005)
        for _ in range(1440):
            state = 16807 * state % 2147483647
            random_value = 2.0 * state / 2147483647.0 - 1.0
            # random() is REAL*8, promoting this expression before assignment.
            guide.append(f(1.0 + noise_scale * random_value))

        total = 0.0
        for value in data:
            total = f(total + value)
        mean = f(total / 1440)
        data = [f(value - mean) for value in data]
        raw = data.copy()
        smooth = data.copy()
        u2 = f(u * f(.08))
        maximum, minimum = f(1.e-5), f(1.e7)
        for i in range(1438):
            difference = f(data[i+2] - data[i])
            maximum = max(maximum, difference)
            minimum = min(minimum, difference)
        quiet = f(maximum - minimum) < factor

        # The reference's window limits use 1-based minute indices. M is
        # daytime in the middle; E is daytime at the ends of the UT day.
        outside = [((i <= qstart or i >= qend) if qindex == "M"
                    else (qstart < i < qend)) for i in range(1, 1441)]
        if quiet:
            u2 = f(u2 * 3.0)
        else:
            if qindex == "M":
                sw = 0.0
                for value in data:
                    sw = f(sw + value)
                sw = f(sw / 1440)
            else:
                sw = 0.0
                for value in data[:4]:
                    sw = f(sw + value)
                sw = f(sw / 4)
            steps = [f(f(2.0 * u) / 3.0) if night else f(2.0 * u) for night in outside]
            for _ in range(2):
                forward = [0.0] * 1440
                reverse = [0.0] * 1440
                w = sw
                for i in range(1440):
                    forward[i] = f(w * guide[i])
                    error = f(data[i] - forward[i])
                    sign = 0 if abs(error) < f(.001) else (1 if error > 0 else -1)
                    w = f(w + f((steps[i] * sign) * guide[i]))
                if qindex == "E":
                    w = data[-1]
                for i in range(1439, -1, -1):
                    reverse[i] = f(w * guide[i])
                    error = f(data[i] - reverse[i])
                    sign = 0 if abs(error) < f(.001) else (1 if error > 0 else -1)
                    w = f(w + f((steps[i] * sign) * guide[i]))
                data = [f(value - f(f(y + yy) / 2.0))
                        for value, y, yy in zip(data, forward, reverse)]
                sw = data[0]
            smooth = [f(a - b) for a, b in zip(raw, data)]
            data = smooth.copy()

        steps = [f(f(2.0 * u2) / 2.0) if night else f(2.0 * u2) for night in outside]
        for _ in range(2):
            forward = [0.0] * 1440
            reverse = [0.0] * 1440
            w = data[0]
            for i in range(1440):
                forward[i] = f(w * guide[i])
                error = f(data[i] - forward[i])
                w = f(w + f(f(steps[i] * error) * guide[i]))
            if qindex == "E":
                w = data[-1]
            for i in range(1439, -1, -1):
                reverse[i] = f(w * guide[i])
                error = f(data[i] - reverse[i])
                w = f(w + f(f(steps[i] * error) * guide[i]))
            data = [f(value - f(f(y + yy) / 2.0))
                    for value, y, yy in zip(data, forward, reverse)]
        smooth = [f(a - b) for a, b in zip(smooth, data)]
        if any(not math.isfinite(v) for v in smooth):
            raise ValueError("Hermanus filter produced non-finite results")
        return [a - b for a, b in zip(raw, smooth)], smooth


def iaf_days(records):
    try:
        return _iaf_days(records)
    except ValueError as error:
        raise ValueError(str(error).replace("KASM", "Hermanus")) from error


def iaga_days(path, sample_step=1):
    try:
        return _iaga_days(path, sample_step=sample_step)
    except ValueError as error:
        raise ValueError(str(error).replace("KASM", "Hermanus")) from error


def calculate_days(days, target_days=None, *, k9_nT=None, minimum_coverage=.9,
                   mu=.17, qindex="M", qstart=300, qend=1009, qfactor=None, seed=1):
    """Fit complete component days only; missing data cannot enter the LRNS filter."""
    validate_parameters(mu, qindex, qstart, qend, qfactor, seed)
    if not 0 < minimum_coverage <= 1:
        raise ValueError("minimum coverage must be greater than 0 and at most 1")
    rows = []
    with Hermanus() as backend:
        for day in sorted(days if target_days is None else target_days):
            x, y, record, _ = days[day]
            limit = k9_nT if k9_nT is not None else (record[10] if record else None)
            if limit is None or not math.isfinite(limit) or limit <= 0:
                raise ValueError("Hermanus requires a finite positive K9 limit")
            residuals = []
            for values in (x, y):
                if len(values) != 1440:
                    raise ValueError("Hermanus requires 1440 samples per component day")
                residuals.append(None if any(v is None for v in values) else backend.calculate(
                    [v / 10 for v in values], mu=mu, qindex=qindex, qstart=qstart, qend=qend,
                    qfactor=limit * .01 if qfactor is None else qfactor, seed=seed)[0])
            for i, (start, end) in enumerate(ih.INTERVALS):
                a, b = i * 180, (i + 1) * 180
                coverage = [sum(v is not None for v in values[a:b]) / 180 for values in (x, y)]
                ranges = [max(r[a:b]) - min(r[a:b]) if r is not None and c >= minimum_coverage
                          else None for r, c in zip(residuals, coverage)]
                indices = [None if r is None else ih.classify_range_pt(r * 1000, limit) for r in ranges]
                k = max(indices) if all(v is not None for v in indices) else None
                stored_raw = record[ih.K_START_WORD + i] if record else None
                stored = None if stored_raw in (None, ih.MISSING_K_VALUE) else stored_raw / 10
                rows.append(dict(
                    date=day.isoformat(), interval_utc=f"{start:02d}:00-{end:02d}:00",
                    x_range_nT=ranges[0], y_range_nT=ranges[1], x_coverage=coverage[0], y_coverage=coverage[1],
                    k_x=indices[0], k_y=indices[1], k=k, stored_k_raw=stored_raw, stored_k=stored,
                    matches_stored="" if k is None or stored is None else int(k == stored)))
    return rows
