# Hermanus LRNS visualization

Open `hermanus.html` in a browser for an offline interactive explanation. It
includes X/Y selection, three actual smoothing settings, a first-L1-pass time
inspector/playback, residual ranges, and combined K. No server or network is
required. The data is Kiruna's first day in `KIR/KIR23DEC.BIN` (1 December 2023).

`hermanus.png` and `hermanus.svg` are standalone figures using X and mu=0.17.
The first L1 cascade is shown separately from the final L1+L2 quiet curve.
K values are cross-checked against the calculator API. The example Hermanus
parameters require station calibration; these are estimated K values.

To regenerate all three artifacts, run:

```
python3 visualizations/build_hermanus_visualization.py
```

The optional builder requires matplotlib for PNG/SVG export. The HTML is built
from `hermanus-template.html` with embedded precomputed data. The calculator
itself still uses only the standard library; visualization dependencies can be
installed in a separate temporary virtual environment.

# USGS visualization

Open `usgs.html` for a self-contained offline explanation using Kiruna data
from December 2023 with November available for neighboring-day support. Change
H/D, the date, and 30/60-minute means; inspect mean acceptance and each three-hour
amplitude decision. D is plotted in equivalent nT using the calibrated XD.

The diagrams separate the raw range, the spline residual range, and the final
statistically selected/adjusted USGS amplitude. Combined decimal K and truncated
integer K are shown together. No new algorithm approximation is used: the
builder calls the actual calibration, mean-selection, segment, spline, and
amplitude routines and cross-checks amplitudes against the public API.

`usgs.png` and `usgs.svg` are standalone figures for H on 1 December 2023,
using 60-minute means. Calibration is bootstrapped from December and the H
reference is estimated from that month; these settings need station validation.

Regenerate with a matplotlib-enabled Python:

```
python3 visualizations/build_usgs_visualization.py
```

The calculator and exported HTML still need no third-party packages.

# USGS monthly selection limits

`usgs-monthly-limits.png` and `.svg` isolate the actual monthly H histograms,
mean-value limits, K5-based range prefilter, and final range ceiling for Kiruna
December 2023. All displayed statistics come from the production bootstrap and
its ten-bin histogram estimator. The spread is not divided by the square root
of the sample count. Regenerate with matplotlib:

```
python3 visualizations/build_usgs_monthly_limits.py
```

# KASM visualization

Open `kasm.html` in a browser for an offline interactive explanation of the
actual two-stage adaptive smoother. Choose X/Y, smoothing 0.5/1/2, and a quiet
Kiruna day (21 December 2023), storm onset (10 May 2024), or the Belsk final-K
cap example (13 May 2024). Inspect hourly first-pass ranges, their adaptive data
influence, final residuals, and each three-hour combined-K decision.

Curves and three-hour ranges use nT. Hourly ranges used to derive weights use
normalized units (nT × 500/K9). The influence plot shows the squared solver
weight relative to 7e-7 × smoothing, matching the coefficient of the data term
in the least-squares objective. Each example includes both 180-minute wings.
The builder calls the production gap-filling and Givens solver routines and
checks every component range, component K and final K against `calculate_days`.
No production calculator code is modified.

`kasm.png` and `kasm.svg` are standalone figures for Kiruna storm onset,
10 May 2024, X, smoothing 1. Regenerate all artifacts with a Python containing
matplotlib (the generated HTML needs no dependencies or network):

```
python3 visualizations/build_kasm_visualization.py
```
