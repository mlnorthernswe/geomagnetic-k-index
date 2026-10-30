# Hermanus LRNS reference implementation

`K_FILTER.FOR` and `HERMANUS.TXT` are unmodified from ISGI's public-domain
Hermanus release, retrieved 2026-10-30:
https://isgi.unistra.fr/Documents/Softwares/hermanus_fortran77.zip

The notes explicitly identify this as the public-domain release, rather than
the commercially available industrialized version. Method by M. Hattingh;
Fortran coding by M. M. Van Dyk, Hermanus Magnetic Observatory, CSIR.

Reference: Hattingh, Loubser & Nagtegaal (1989), “Computer K-index estimation
by a new linear-phase, robust, non-linear smoothing method”, GJI 99, 533–547,
https://doi.org/10.1111/j.1365-246X.1989.tb02038.x

## Python port and source repair

`calculate_hermanus_k_index.py` implements LRNS in standard-library Python.
No compiler, native library, subprocess, temporary file, or `FC` setting is
required to calculate indices. It preserves the reference's REAL*4 arithmetic
by rounding elementary operations to IEEE binary32 with `struct`, including
input conversion, sequential mean accumulation, and the two L1/L2 cascades.
Mixed REAL*8 guide-noise initialization retains its original precision.
K values use the application's standard K9 thresholds. Raw and quiet arrays
have the same mean removed, so their difference preserves the residual.

The archive contains an apparent transcription error: the forward loop's
`do 80 i = n1, n2` is commented out, despite a corresponding label 80 and
array operations requiring a forward loop. The Python port includes this
loop. The original archive files are retained byte-for-byte for review.
Only the optional fixture generator restores that statement and strips DOS
EOF bytes in a temporary compilation copy; it does not alter arithmetic.

## Random guide and station parameters

The original routine requires an external double-precision `random()` in [-1,1],
but does not supply one. The Python port supplies Park–Miller (16807, modulus
2147483647), mapped to [-1,1], with an explicit seed reset for every component/day.
The default seed is 1. The guide noise is required by the reference to stimulate
adaptation on very quiet days. The Python port keeps all state local, so independent calls do not share RNG state.

The defaults mu=.17, qindex=M, qstart=300, qend=1009 are the supplied Hermanus
example. They require observatory calibration. qfactor defaults to .01*K9 nT
and can be overridden. qstart/qend retain the reference's 1-based minute indices.
For E, the interval is nighttime between the two daytime portions.
The source's other initialization and update details are retained exactly.

The original filter expects a complete spike/step-free one-minute day. We do
not add gap interpolation: if any minute is missing in a component, that
component's full day is unavailable. The other component can still be reported.
No adjacent-day data is required.

## Independent validation

`reference_driver.f90` is a standalone harness with its own RNG provider. It
calls the repaired original `kfilter` directly, independently of the Python port.
`tests/generate_hermanus_reference.py` compiles it at -O0 with -fcheck=all,
then generates `tests/data/hermanus-reference.json`. Tests compare the Python
port's interval ranges and selected raw/smooth trace values with this fixture.
These tests run without a Fortran compiler. GNU Fortran (`gfortran` or `FC`)
is needed only for deliberate fixture regeneration, not normal tests or use.
The fixture covers complete days at all four example stations, M/E windows,
mu=.08/.17/.20, several seeds, and the qfactor branch that skips rough smoothing.
Source and input SHA-256 hashes are recorded. Regenerate from the repository root:

```
python3 tests/generate_hermanus_reference.py
```

Synthetic tests also cover constant input, reproducibility, missing-data policy,
CLI IAF/IAGA input, pipeline configuration, and invalid settings.
