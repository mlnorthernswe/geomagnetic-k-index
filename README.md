# INTERMAGNET K-index calculator

A standard-library Python tool for estimating three-hour local K indices from INTERMAGNET IAF monthly `.BIN` files. Example station data for BEL, HYB, KIR, and WIC are included.

## Requirements

Python 3.10 or newer. No third-party packages are required.

## Quick start

From the repository root, calculate WIC December 2023 indices, bootstrap an SQ curve from that month, save the curve, and write both CSV and IRF comprehensive output:

```bash
python3 calculate_ih_k_index.py WIC/WIC23DEC.BIN \
  --sq-output /tmp/WIC23DEC-sq.txt \
  --comprehensive-output WIC/WIC23DEC-comprehensive.txt \
  -o WIC/WIC23DEC-k-index.csv
```

If no SQ input or bootstrap file is specified, the tool bootstraps from the month being processed. It ranks days using the eight K values embedded in each IAF daily record, selects the quietest days (10 by default), and takes hourly X/Y/Z medians to form an SQ curve.

To use an existing Kiruna SQ file:

```bash
python3 calculate_ih_k_index.py KIR/KIR23DEC.BIN \
  --sq-input KIR/202311KIRsq \
  --sq-output /tmp/KIR-used-sq.txt
```

An alternate IAF month can be selected as the bootstrap source with `--bootstrap-file`. The SQ format is 24 whitespace-separated rows of X, Y, and Z hourly values in pT. The legacy positional SQ-file argument is also accepted.

## Options

Run `python3 calculate_ih_k_index.py --help` for all options.

- `--sq-input PATH`: use an existing 24-row SQ file.
- `--sq-output PATH`: save the SQ curve used by this run in the same 24-row format.
- `--bootstrap-file PATH`: use another IAF month file to construct the SQ curve; defaults to the month being processed.
- `--bootstrap-days N`: number of quietest valid days to use for bootstrap (default: 10).
- `--k9-limit NT`: override the station K9 range; by default it is read from the IAF header.
- `--minimum-coverage FRACTION`: minimum valid minute coverage for each component and interval (default: 0.9).
- `-o PATH`, `--output PATH`: write detailed per-interval results as CSV (default: standard output).
- `--comprehensive-output PATH`: also write daily K values in the IRF `comprehensive.txt` layout.

## Output and limitations

The CSV includes X/Y deflection ranges, coverage, component K values, combined K, and comparisons with K values embedded in the IAF file. The combined K is the larger of the separately classified X and Y indices. The comprehensive output writes eight K values per date, `-` for unavailable intervals, parenthesized partial sums, and `-1` when no interval K value is available.

IAF stores one-minute means, while the original SQ/K method described in the accompanying algorithm uses one-second samples. Calculated K values are therefore estimates and may not match the embedded values exactly.

## Algorithm reference

- `K-index-calculations-algorithm.pdf`: original algorithm description.
- `K-index-calculations-algorithm-revised.pdf`: revised description with bootstrap guidance.
