# Wilson USGS method

Original sources and documentation retrieved 2026-10-30 from ISGI:
https://isgi.unistra.fr/Documents/Softwares/usgs_fortran77.zip
https://isgi.unistra.fr/softwares.php

The outer archive contains `USGS.TXT`, `DKTXT.TXT`, and `USGS.ZIP`.
The nested archive contains the original USGSDK/SELMVS kernels and legacy
programs. The `.FOR` files and original `.TXT` documentation here are retained
byte-for-byte. The supplied `NTU507.RUN` illustrates runtime settings. The
`.KAY` file NT8586 in the archive is historical hand-scaled data, not our JSON
calibration format. Historical platform-specific executables are not bundled.

Reference: Lanny R. Wilson (1987), “An evaluation of digitally derived K-indices”,
Journal of Geomagnetism and Geoelectricity 39, 97–109. The distributed program
was written for Digital Equipment systems and ported to Lahey F77L/MS-DOS.

## Python implementation

`calculate_usgs_k_index.py` ports the computational workflow:

- SELMVS/MSTAT: rounded block means, within-block ranges, ten-bin monthly
  statistics, mean/range confidence limits, and ten quiet dates ranked by H.
- KMEAN: selection/truncation of accepted quiet means and daily mean/mode fallback.
- XYSET/GS: adjacent-day support points, sliding regressions, continuity repairs,
  and intersections used as spline knots.
- NTRP0/SPLN2/CSPLN2: gap insertion, clamped spline solution and minute evaluation.
- GKR: amplitude/range confidence tests and significant linear-trend correction.
- KTAB/CURVK/RKV/RND: integer-nT amplitude rounding and decimal digital K.

Fortran REAL*4 arithmetic is emulated with `struct` and a rounded float type,
including single-precision real literals, accumulation, and operation order.
The original RND behavior (ties to even, with its negative-zero sign quirk),
KMEAN's historical wide-range handling, and SELMVS's published D-range cutoff
formula are retained. The derivative seeds and insertion behavior of the spline
are preserved, including which slope is retained when a knot is inserted.

Input adaptation is explicit: X/Y minute means become H/D, declination is
unwrapped, and baselines are removed before 16-bit quantization. Previous-year
H may be supplied; the default monthly H estimate is an approximation to that
reference input. Baselines, confidence limits, conversion factor and quiet dates
are saved in calibration JSON. Bootstrap is recomputed by calendar month.
Tied quiet days are ordered chronologically. JSON is not a legacy `.KAY` parser.
The runtime defaults disable the optional legacy nighttime range gate. Confidence
multipliers and segment size can be changed in calibration JSON. K0 thresholds
must be representable after integer-nT rounding; this adapter requires K9 >=100 nT.
Both components must be valid. CSV names H and D explicitly and exposes decimal
K; comprehensive output uses truncated integer K as in the reference's daily sum.
The reference's >6000 nT digital-K missing-value rule is retained.

No compiler is required for normal calculations or tests. Missing data is
represented by None, preventing a legitimate reading from colliding with 32767.
The reference permits up to 59 missing minutes per three-hour interval; the
application's observed-coverage threshold also applies.

## Defined repairs and portability

The public archive has platform-dependent syntax and a few unsafe edge cases.
The Python implementation defines safe behavior rather than emulating memory
corruption. The optional native fixture generator applies the matching repairs:

- Convert Lahey function result declarations and FLOATI/FLOATJ/ININT/IABS/IMOD
  spellings to equivalent GNU Fortran syntax; repair a missing comma in an
  inactive GKR debug FORMAT, strip DOS EOF markers and expand tabs.
- Declare XYSET's NPTS as INTEGER*2 to match its caller and GS. Its original
  implicit INTEGER*4 writes beyond the caller's two-byte argument.
- Stop NTRP0 at its final valid knot interval (`I.GE.NR`, not `I.GT.NR`). The
  original loop can read stale knots beyond NS and insert unintended points.
- Include a histogram's maximum endpoint in its last bin, instead of writing
  beyond the ten bins. Constant/single-sample statistics get zero spread in
  Python; the original zero-width histogram is undefined.
- Reject unrepresentable encoded data/histogram spans. Limit the regression
  window to available coordinate pairs on sparse input. A degenerate spline
  produces a clear error instead of STOP, division by zero, or an invalid index.

## Reference validation

`reference_driver.f90` calls the original KMEAN, MOVDAT, XYSET, GS, SPLN2,
NTRP0, CSPLN2 and GKR routines with explicit calibrated parameters and H/D data.
It then calls CURVK/RKV for decimal K. `histogram_driver.f90` independently tests
MSTAT. Neither harness imports the Python numerical kernel.

`tests/data/usgs-reference.json` stores original-source/input hashes, parameter
sets, three native histograms, and 22 native cases covering all four observatories,
missing neighbors, 12/30/60-minute means, and 3/4-mean segment fits. The tests
compare all 352 component amplitudes and corresponding decimal K values.
Native generation uses -O0, no fused arithmetic, and consistent COMMON layout.
Tests run entirely in Python against the saved fixture.

Developers can deliberately regenerate the fixture with GNU Fortran:

```
python3 tests/generate_usgs_reference.py
```

Only this optional generator uses `gfortran` or `FC`; no code is downloaded at
runtime. The reference driver covers algorithm kernels rather than the original
platform-specific NGIC file I/O, plotting, or interactive UI programs.
