# Native FMI/KASM port validation

`native-port-reference.json` was captured on 2026-10-30 before replacing the
native backends with standard-library Python. It contains 40 cases (20 per
method), with all eight combined and component K values, both component ranges,
and KASM gap flags returned by the existing native diagnostic wrappers.

For each method the cases include three December 2023 target days at all four
example stations (12 cases) and eight synthetic cases. FMI uses three full
days; KASM uses the middle day plus 180-minute wings for station inputs.
Input file SHA-256 hashes are stored in each case. Native sources are the
bundled FMI 1.3 `vendor/fmi/K_index.h` and KASM 1.11 `vendor/kasm/casm.cpp`,
with their existing diagnostic wrappers. Attribution and original URLs are
recorded in the respective vendor READMEs.

Synthetic cases use deterministic sin/cos inputs specified by
`tests/native_reference_inputs.py`. They cover rising/falling short gaps,
gaps at the interpolation cutoff, grouped gaps separated by single samples,
leading/trailing gaps, a fully missing component, interval-boundary gaps,
western longitude, and altered KASM smoothing. All KASM context lengths
(1440, 1620, 1800 minutes) are included. FMI's last case uses nonconstant
one-second input (259200 samples per component).

The ports were additionally checked against fresh native pipeline output for
all 16 included station-months per method and 80 randomized KASM gap cases.
Normal tests use the saved fixtures and never compile code. KASM range tests
allow 1e-8 nT variation to accommodate final-bit libm/platform differences;
K values and flags are exact. FMI's integer diagnostics are exact.

Developers can regenerate expectations from the recorded input definitions:

```
python3 tests/generate_native_port_reference.py
```

This optional generator requires C and C++ compilers (`CC`/`CXX`), builds the
bundled references in a temporary directory, and never mutates input data or
downloads code. Normal calculation and testing require only Python.
