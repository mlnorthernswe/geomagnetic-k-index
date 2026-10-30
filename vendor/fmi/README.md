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

The adapter normalizes longitude to [0, 360) to avoid negative C array indices
in the reference implementation. Calls must be sequential: FMI uses globals.
The bundled tests/data/fmi-demo.dat is FMI's published reference dataset:
https://space.fmi.fi/MAGN/K-index/FMI_method/demo.dat
