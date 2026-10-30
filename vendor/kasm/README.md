# KASM reference implementation

`casm.cpp` and `casm.h` are unmodified files from Tomasz Ernst's KASM 1.11
(May 2017) source distribution, Institute of Geophysics, Polish Academy of Sciences.
Downloaded 2026-10-30 from the author's distribution, linked by INTERMAGNET:

- https://magneto.igf.edu.pl/soft/kasm/KAsm.zip
- https://intermagnet.org/software.html
- https://academic.oup.com/gji/article/104/1/85/570315

Method reference: Nowozyński, Ernst & Jankowski (1991), “Adaptive smoothing
method for computer derivation of K-indices”, GJI 104, 85–93,
doi:10.1111/j.1365-246X.1991.tb02495.x.

This is the author's newer C++ implementation, rather than the older Fortran
routines at https://isgi.unistra.fr/softwares.php. The newer routine has no
longitude or annual H requirement when using X/Y data. It uses minute means,
a two-stage weighted least-squares fit with a second-difference penalty, and
Givens rotations. Defaults in `casm.cpp` are del=5e-7 and dela=7e-7, multiplied
by the configurable smoothing factor (default 1). Input is normalized by 500/K9.

`wrapper.cpp` follows the original `KAsm::allK`, `norm`, and `wwsk` stages:
use 180 minutes from each available adjacent day, fill gaps with `fits` (six
consecutive missing minutes allowed; longer gaps invalidate affected intervals),
normalize, call `rs::rs1`, classify with `kid`, and cap combined K at the
unfiltered combined K. Leading/trailing gaps invalidate touched intervals.
Component K and range diagnostics describe the filtered residuals, so combined
K may be smaller than max(k_x, k_y) because of this final cap. Both components
must be valid, and the application's observed-coverage cutoff also applies.

The wrapper catches C++ exceptions at the C ABI boundary and reports failure.
The library is compiled into a temporary directory with a C++11 compiler.
No files are downloaded at runtime. Reference validation builds and runs the
complete author's program separately and compares all eight indices per day.
