# FMI reference implementation

`K_index.h` is the unmodified FMI C implementation, version 1.3 (1999),
by Lasse Hakkinen, Finnish Meteorological Institute. Retrieved 2026-10-30 from:
https://space.fmi.fi/MAGN/K-index/FMI_method/K_index.h

Documentation: https://space.fmi.fi/MAGN/K-index/FMI_method/K_index.txt
Method: Sucksdorff et al. (1991), as referenced by FMI.

`wrapper.c` exposes component diagnostics using the same calculation stages.
The upstream source retains its historical integer arithmetic, strict threshold
comparisons, and K_MaxMin loop (which omits the last sample of each interval).
Short gaps strictly less than 15 minutes are interpolated, per the actual code.

The Python port in `calculate_fmi_k_index.py` preserves these details, including
C integer division toward zero for descending gaps and binary32 rounding in
the harmonic fit. Longitude is normalized to [0, 360), as in the former C
adapter, to avoid the source's negative-index issue for western observatories.
Calculation state is local; compact cached harmonic bases accelerate repeated
days at the same cadence. Prefix sums speed up activity-dependent averaging.

No compiler or native library is required at runtime. `wrapper.c` is retained
only for developer reference-fixture generation; it is not imported or compiled
by the calculator or normal tests. The native source uses globals, so optional
native reference generation must remain sequential.
The bundled tests/data/fmi-demo.dat is FMI's published reference dataset:
https://space.fmi.fi/MAGN/K-index/FMI_method/demo.dat


`tests/data/native-port-reference.json` includes native C results for station
data, descending gaps, boundary gaps, missing components, and a nonconstant
one-second example. See `tests/data/native-port-reference.md` for provenance.
The published demo and diagnostic tests run entirely in Python. Only deliberate
fixture regeneration (`python3 tests/generate_native_port_reference.py`)
requires C/C++ compilers.
