"""Standard-library Python port of the published FMI 1.3 K-index method."""
from __future__ import annotations

import functools
import math
import struct
from array import array
from datetime import datetime, timedelta

import calculate_k_index as ih


_FLOAT32 = struct.Struct("=f")
_MISSING = 2**50
_PI = 3.141592654  # Preserve the published approximation.


def _f32(value):
    return _FLOAT32.unpack(_FLOAT32.pack(value))[0]


def _cdiv(a, b):
    """C integer division truncates toward zero, including falling gaps."""
    return (abs(a) // abs(b)) * (-1 if (a < 0) != (b < 0) else 1)


@functools.lru_cache(maxsize=2)
def _harmonic_basis(step):
    """Single-precision trigonometry, shared across days, in compact arrays."""
    angle = _f32(2.0 * _PI * (23.0 / 24.0) / 82800)
    result = []
    for k in range(1, 6):
        sine, cosine = array("f"), array("f")
        for t in range(0, 86400, step):
            angle2 = _f32(_f32(t - 1800) * angle)
            phase = _f32(k * angle2)
            sine.append(math.sin(phase))
            cosine.append(math.cos(phase))
        result.append((sine, cosine))
    return result


def _fill_gaps(data, step):
    i, n = 0, len(data)
    while i < n:
        if data[i] != _MISSING:
            i += 1
            continue
        start = i
        while i < n and data[i] == _MISSING:
            i += 1
        length = i - start
        if start > 0 and i < n and length * step < 900:
            left, delta = data[start-1], data[i] - data[start-1]
            for j in range(1, length+1):
                data[start+j-1] = left + _cdiv(j*delta, length+1)


def _interpolate_means(values):
    valid = [i for i, value in enumerate(values) if value != _MISSING]
    if not valid:
        return
    first, last = valid[0], valid[-1]
    values[:first] = [values[first]] * first
    values[last+1:] = [values[last]] * (23-last)
    for left, right in zip(valid, valid[1:]):
        for i in range(left+1, right):
            values[i] = values[left] + _cdiv((i-left)*(values[right]-values[left]), right-left)


def _fit_curve(means, step):
    """Fifth-order harmonic fit, retaining C float promotions and long casts."""
    f = _f32
    coefficient = f(f(means[23] - means[0]) / f(82800))
    means = [value - int(f(coefficient * f(i*3600))) for i, value in enumerate(means)]
    real, imaginary = [], []
    for k in range(6):
        re, im = f(means[0]), 0.0
        angle = f(-k * (2.0 * _PI / 24.0))
        for j in range(1, 24):
            phase = f(j * angle)
            sine, cosine = f(math.sin(phase)), f(math.cos(phase))
            re = f(re + f(f(means[j]) * cosine))
            im = f(im + f(f(means[j]) * sine))
        real.append(re)
        imaginary.append(im)
    result = []
    basis = _harmonic_basis(step)
    for i, t in enumerate(range(0, 86400, step)):
        value = int(real[0])
        for k, (sine, cosine) in enumerate(basis, 1):
            difference = f(f(real[k] * cosine[i]) - f(imaginary[k] * sine[i]))
            value += int(2.0 * difference)
        result.append(int(value / 24.0 + f(coefficient * f(t-1800))))
    return result


class FMI:
    """Standard-library port of FMI 1.3; calculation state is local to a call."""

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def calculate(self, x, y, step, k9_nT, longitude):
        if step not in (1, 60):
            raise ValueError("FMI adapter supports one-second or one-minute samples")
        n = 86400 // step
        if len(x) != 3*n or len(y) != len(x):
            raise ValueError("FMI requires exactly three days of X and Y samples")
        if not math.isfinite(k9_nT) or k9_nT <= 0 or round(k9_nT * 10) < 100:
            raise ValueError("FMI K9 limit must be finite and at least 10 nT")
        if not math.isfinite(longitude) or not -180 <= longitude <= 360:
            raise ValueError("longitude must be finite and between -180 and 360 degrees")
        components = [[_MISSING if v is None else int(v) for v in values] for values in (x, y)]
        for values in components:
            _fill_gaps(values, step)
        limits = [round(k9_nT*10)*v//7500 for v in (75,150,300,600,1050,1800,3000,4950,7500)]
        defaults = [5400]*3 + [3600]*3 + [0]*12 + [3600]*3 + [5400]*3
        shift = (round(longitude) % 360) // 15
        defaults = [defaults[(i+shift) % 24] for i in range(24)]
        lengths = [0] + [min(int(60*math.exp(3.3*math.log(k))), 64800) for k in range(1,10)]

        def classify(values):
            indices, ranges = [], []
            for i in range(8):
                # The upstream loop omits the final sample of each interval.
                part = values[i*n//8:(i+1)*n//8-1]
                high, low = max(part), min(part)
                ranges.append(high-low)
                indices.append(-1 if high == _MISSING else sum(high-low > v for v in limits))
            return indices, ranges

        # Prefix sums accelerate large activity-dependent averaging windows;
        # missing counts prevent the sentinel from entering any arithmetic.
        prefixes = []
        for values in components:
            sums, gaps = [0], [0]
            for value in values:
                sums.append(sums[-1] + (value if value != _MISSING else 0))
                gaps.append(gaps[-1] + (value == _MISSING))
            prefixes.append((sums, gaps))
        residuals = [values[n:2*n] for values in components]
        for iteration in range(3):
            kx, rx = classify(residuals[0])
            ky, ry = classify(residuals[1])
            ks = [max(a,b) for a,b in zip(kx,ky)]
            if iteration == 2:
                return [ks, kx, ky, rx, ry]
            fit_lengths = [min(1800+defaults[i]+(lengths[ks[i//3]] if ks[i//3] >= 0 else 0),86400)
                           for i in range(24)]
            residuals = []
            for values, (sums, gaps) in zip(components, prefixes):
                means = []
                for hour, length in enumerate(fit_lengths):
                    middle, wing = n + (3600*hour+1800)//step, length//step
                    start, end = middle-wing, middle+wing+1
                    means.append(_MISSING if gaps[end] != gaps[start]
                                 else int((sums[end]-sums[start]) / (end-start)))
                _interpolate_means(means)
                harmonic = _fit_curve(means, step)
                residuals.append([a-b if a != _MISSING and b != _MISSING else _MISSING
                                  for a,b in zip(values[n:2*n], harmonic)])

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
