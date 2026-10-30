# KASM reference fixture provenance

`kasm-reference.json` was generated on 2026-10-30 using the complete author's
KASM 1.11 source distribution:
https://magneto.igf.edu.pl/soft/kasm/KAsm.zip

The archive SHA-256 and the SHA-256 of each input IAF file are embedded in the
fixture. The full standalone program was built separately with its supplied
Makefile (`make`, g++ -std=c++11 -O2). It did not link our adapter.

For each station, copies of November and December 2023 IAF files were placed
in a temporary directory using lowercase filenames (required by the reference
program on Unix). No January file was supplied, so December 31 uses reduced
context. The command recorded in the JSON was run without `-u` (no writes to
IAF input); eight values per day were parsed from its `.dka` output.
K9 limits: KIR 1500, BEL 450, WIC 500, HYB 300 nT. All stations use smoothing
1; WIC is additionally tested at smoothing 2. There are 1,240 expected values.

The reference program has no extra observed-coverage threshold. The oracle
comparison therefore uses minimum_coverage=0.01, retaining all otherwise valid
intervals while still applying native gap flags. Separate tests check the
application's default coverage cutoff of 0.9, including its additional rejection
of intervals with many isolated short gaps. The reference is not a comparison
against the K values embedded in the original data files.

To reproduce, download and unpack the archive into a temporary directory,
run `make` in its `src` directory, copy the listed input files to another
temporary directory with lowercase names, then run for each fixture case:

```
KAsm STA:01DEC2023:31 K9:SMOOTHING xy OUT INPUT_DIRECTORY/ -IAF
```
