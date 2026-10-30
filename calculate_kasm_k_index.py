"""Standard-library Python port of the author's KASM 1.11 method."""
from __future__ import annotations

import bisect
import math
import struct
from datetime import timedelta

import calculate_k_index as ih
from calculate_fmi_k_index import iaf_days as _iaf_days


_UNDEF = 1e30


class _Givens:
    """Stable rotations with the same operation order as KASM's giv class."""
    def __init__(self):
        self.c, self.s = 1.0, 0.0

    def rot(self, a, b):
        aa, bb = abs(a), abs(b)
        if aa > bb:
            q = b / a
            r = math.sqrt(1.0 + q*q)
            self.c = (-1.0 if a < 0 else 1.0) / r
            self.s = q * self.c
            return r * aa
        if bb > 0:
            q = a / b
            r = math.sqrt(1.0 + q*q)
            self.s = (-1.0 if b < 0 else 1.0) / r
            self.c = q * self.s
            return r * bb
        self.c, self.s = 1.0, 0.0
        return 0.0

    def __call__(self, a, b):
        return self.c*a + self.s*b, self.c*b - self.s*a


def _solve(values, weights):
    """Givens QR of [second differences; diagonal weights], returning residuals."""
    n = len(values)
    d, p, q, output = ([0.0]*n for _ in range(4))
    r, rlc, rln = 0.0, weights[0], weights[1]
    fc, fn = rlc*values[0], rln*values[1]
    gv = _Givens()
    for i in range(n-2):
        d[i] = gv.rot(1.0, rlc)
        p[i], r = gv(-2.0, r)
        q[i], s = gv(1.0, 0.0)
        output[i], fp = gv(0.0, fc)
        fc, fn = fn, weights[i+2]*values[i+2]
        rlc = gv.rot(rln, r)
        r, s = gv(0.0, s)
        fc, fp = gv(fc, fp)
        rln = gv.rot(weights[i+2], s)
        fn, fp = gv(fn, fp)
    output[-1] = fn / rln
    output[-2] = fc / rlc - (r / rlc) * output[-1]
    for i in range(n-3, -1, -1):
        output[i] = (output[i] - p[i]*output[i+1] - q[i]*output[i+2]) / d[i]
    return [a-b for a,b in zip(values, output)]


def _fit_gap(values):
    """Original KASM minimum-curvature interpolation of a group of gaps."""
    x = values.copy()
    m = len(x)
    ml = m - (2 if x[m-2] < _UNDEF else 1)
    n = sum(v >= _UNDEF for v in x[1:ml])
    b, d, f, g = ([0.0]*n for _ in range(4))
    gv = _Givens()
    k = 0
    if x[1] >= _UNDEF:
        r, t, j = 0.0, 0.0, 1
    else:
        r, t, j = 1.0, -x[0]+2.0*x[1], 2
    while j < ml:
        b[k], d[k] = -x[j-1], -2.0
        j += 1
        if x[j] < _UNDEF:
            b[k] += -x[j]
        else:
            f[k] = 1.0
            j += 1
            while x[j] >= _UNDEF:
                d[k] = gv.rot(d[k], r)
                f[k], r = gv(f[k], 0.0)
                b[k], t = gv(b[k], t)
                d[k] = gv.rot(d[k], 1.0)
                f[k], d[k+1] = gv(f[k], -2.0)
                g[k], f[k+1] = gv(0.0, 1.0)
                b[k], b[k+1] = gv(b[k], 0.0)
                k += 1
                j += 1
            d[k] = gv.rot(d[k], r)
            f[k], r = gv(f[k], 0.0)
            b[k], t = gv(b[k], t)
            d[k] = gv.rot(d[k], 1.0)
            f[k], d[k+1] = gv(f[k], -2.0)
            b[k], b[k+1] = gv(b[k], -x[j])
            k += 1
        d[k] = gv.rot(d[k], r)
        b[k], t = gv(b[k], t)
        if j != m-1:
            d[k] = gv.rot(d[k], 1.0)
            t = 2.0*x[j]
            if j == ml:
                t += -x[j+1]
            else:
                f[k], r = gv(0.0, 1.0)
                j += 1
            b[k], t = gv(b[k], t)
        k += 1
    k, j = n-1, ml-1
    yd = b[k] / d[k]
    x[j] = yd
    if k > 0:
        k -= 1
        yf, yd = yd, (b[k]-yd*f[k])/d[k]
        j -= 1
        if x[j] < _UNDEF:
            j -= 1
        x[j] = yd
        for k in range(k-1, -1, -1):
            yg, yf = yf, yd
            yd = (b[k]-yf*f[k]-yg*g[k])/d[k]
            j -= 1
            if x[j] < _UNDEF:
                j -= 1
            x[j] = yd
    return x


def _fill_gaps(x):
    """Port of fits(), including repeated flag counts and boundary rules."""
    n, i = len(x), 0
    flags = [0]*10
    while i < n and x[i] >= _UNDEF:
        i += 1
    xi = x[i] if i < n else 0.0
    x[:i] = [xi]*i
    for j in range((i-1+180)//180):
        flags[j] += 1
    m = n-1
    while m > i and x[m] >= _UNDEF:
        m -= 1
    if m < n-1:
        x[m+1:] = [x[m]]*(n-1-m)
        for j in range((m+1)//180, (n-1+180)//180):
            flags[j] += 1
    while i < m:
        length = 0
        for k in range(i+1, m+1):
            if x[k] >= _UNDEF:
                length += 1
            else:
                if length > 6:
                    for j in range((k-length)//180, (k-1+180)//180):
                        flags[j] += 1
                length = 0
        lp = i
        i += 1
        while i < m and x[i] < _UNDEF:
            i += 1
        if i < m:
            if i-lp > 2:
                lp = i-2
            i += 2
            while i < m and (x[i-1] >= _UNDEF or x[i] >= _UNDEF):
                i += 1
            x[lp:i+1] = _fit_gap(x[lp:i+1])
            if i < m:
                i -= 1
    return flags


def _classify(value):
    return bisect.bisect_right((5,10,20,40,70,120,200,330,500), value)


def _range(values):
    return max(values) - min(values)


class KASM:
    """Standard-library Python port of the author's KASM 1.11 solver."""
    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

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
        try:
            scale = 500.0 / k9_nT
            normalized, residuals, flags = [], [], []
            for values in (x, y):
                data = [_UNDEF if v is None else float(v) for v in values]
                flags.append(_fill_gaps(data))
                low = min(data)
                data = [scale*(v-low) for v in data]
                normalized.append(data)
                initial = _solve(data, [math.sqrt(5.e-7*smoothing)]*len(data))
                weights = []
                base = math.sqrt(7.e-7*smoothing)
                for i in range(0, len(data), 60):
                    q = min(_range(initial[i:i+60])/6.66-1.0, 50.0)
                    weights.extend([base*math.exp(-q/2.0)]*60)
                residuals.append(_solve(data, weights))
            ks, kx, ky, rx, ry, invalid_x, invalid_y = ([] for _ in range(7))
            for i in range(8):
                start, end = offset+i*180, offset+(i+1)*180
                ranges = [_range(values[start:end]) for values in residuals]
                if any(not math.isfinite(v) for v in ranges):
                    raise ValueError("KASM numerical calculation failed")
                ix, iy = (_classify(v) for v in ranges)
                kx.append(ix)
                ky.append(iy)
                ks.append(min(max(ix, iy), _classify(max(_range(v[start:end]) for v in normalized))))
                rx.append(ranges[0]/scale)
                ry.append(ranges[1]/scale)
                invalid_x.append(flags[0][offset//180+i])
                invalid_y.append(flags[1][offset//180+i])
            return [ks, kx, ky, rx, ry, invalid_x, invalid_y]
        except ArithmeticError as error:
            raise ValueError("KASM numerical calculation failed") from error


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
