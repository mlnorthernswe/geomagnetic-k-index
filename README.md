# Geomagnetic K-index calculator

A standard-library Python tool for estimating three-hour local K indices from INTERMAGNET IAF monthly `.BIN` files. Example station data for BEL, HYB, KIR, and WIC are included.

## Requirements

Python 3.10 or newer. No third-party Python packages are required. The optional FMI backend requires a C compiler (`cc`, or set `CC`) on a 64-bit Unix system. KASM requires a C++11 compiler (`c++`, or set `CXX`) on Unix. Hermanus LRNS is implemented in standard-library Python and requires no compiler.

## Quick start

From the repository root, calculate WIC December 2023 indices, bootstrap an SQ curve from that month, save the curve, and write both CSV and IRF comprehensive output:

```bash
python3 calculate_k_index.py WIC/WIC23DEC.BIN \
  --sq-output /tmp/WIC23DEC-sq.txt \
  --comprehensive-output WIC/WIC23DEC-comprehensive.txt \
  -o WIC/WIC23DEC-k-index.csv
```

If no SQ input or bootstrap file is specified, the tool bootstraps from the month being processed. It ranks days using the eight K values embedded in each IAF daily record, selects the quietest days (10 by default), and takes hourly X/Y/Z medians (or means with `--sq-estimator mean`) to form an SQ curve. High-latitude stations (such as Kiruna) should use `median` (the default) because isolated auroral disturbances can occur even on quiet days.

To calculate from an IAGA-2002 file (plain ASCII or compressed `.gz`, `.bz2`, `.xz`) containing one-second X/Y samples, select the input format explicitly and provide an SQ curve and the station's K9 range:

```bash
python3 calculate_k_index.py KIR/kir202609qsec.sec.gz --input-format iaga \
  --sq-input KIR/202310KIRsq --k9-limit 1500 \
  -o /tmp/iaga-k-index.csv
```

IAGA component values are interpreted as nT and SQ values as pT. The SQ curve is interpolated to each sample timestamp. Alternatively, `--bootstrap-file` can point to an IAF month to construct the SQ curve. IAGA files do not contain stored K values, so the CSV's stored-K comparison columns are empty.

To use an existing Kiruna SQ file:

```bash
python3 calculate_k_index.py KIR/KIR23DEC.BIN \
  --sq-input KIR/202311KIRsq \
  --sq-output /tmp/KIR-used-sq.txt
```

An alternate IAF month can be selected as the bootstrap source with `--bootstrap-file`. The SQ format is 24 whitespace-separated rows of X, Y, and Z hourly values in pT. The legacy positional SQ-file argument is also accepted.

## Select the FMI method

IH remains the default. Select FMI with `--method fmi`:

```bash
python3 calculate_k_index.py KIR/KIR23DEC.BIN --method fmi \
  --context-file KIR/KIR23NOV.BIN \
  -o /tmp/KIR23DEC-k-index-FMI.csv \
  --comprehensive-output /tmp/KIR23DEC-comprehensive-FMI.txt

python3 run_pipeline.py --method fmi --output-dir /tmp/fmi-results
```

FMI uses its own daily quiet-curve fit, rather than an IH monthly SQ file.
SQ input/output and bootstrap files cannot be combined with FMI. The optional
backend compiles the bundled, attributed [FMI reference C implementation](vendor/fmi/README.md)
in a temporary directory; it does not download code at runtime.

For each target day, FMI needs the previous and following calendar days.
The calculator accepts repeatable `--context-file` arguments in the same input
format, from the same station; only dates in the main input file are output.
The pipeline automatically uses available neighboring months, including across
year boundaries. Days lacking either neighboring day have unavailable K values
(blank CSV cells and `-` in comprehensive output), and the run summary reports
these days. In the December example above, December 31 remains unavailable
unless January 1 is also supplied as context. Date gaps are never bridged.

IAF input uses minute samples and obtains longitude and K9 from the header;
`--longitude` (degrees east) and `--k9-limit` override these values. For
one-second IAGA input, both parameters are required:

```bash
python3 calculate_k_index.py station-three-days.sec.gz \
  --input-format iaga --method fmi --longitude 20.42 --k9-limit 1500 \
  -o /tmp/station-FMI.csv
```

IAGA timestamps must fall on whole seconds; absent timestamps are treated as
missing samples. Values are converted to 0.1 nT for the reference routine.
FMI interpolates bounded gaps strictly shorter than 15 minutes. Longer gaps
can invalidate an interval. The coverage columns report observed samples before
interpolation, and `--minimum-coverage` still applies. As with IH, combined K
requires valid results from both horizontal components.

FMI retains the reference implementation's integer arithmetic and threshold
comparisons. Its included demonstration produces `1 1 1 2 3 3 5 3`.
One-minute IAF results remain estimates of indices derived at finer resolution.

The pipeline also accepts `"method": "fmi"` in `stations.json` defaults or in an
individual station configuration; `--method` overrides configuration. Optional
`"longitude"` values override IAF header longitude. FMI pipeline outputs have
`-FMI` filenames, while IH outputs keep their existing `-IH` filenames.
`--save-generated-sq` applies only to IH.

## Select the KASM method

KASM (Adaptive Smoothing) is a third selectable method. This tool uses the author's
[KASM 1.11 C++ implementation](https://magneto.igf.edu.pl/soft/kasm/), linked from
[INTERMAGNET](https://intermagnet.org/software.html), with original source and
integration details in [vendor/kasm/README.md](vendor/kasm/README.md).
IH remains the default.

```bash
python3 calculate_k_index.py KIR/KIR23DEC.BIN --method kasm \
  --context-file KIR/KIR23NOV.BIN \
  -o /tmp/KIR23DEC-k-index-KASM.csv

python3 run_pipeline.py --method kasm --output-dir /tmp/kasm-results
```

KASM first estimates a smooth daily curve by penalizing second differences,
then uses hourly disturbance ranges to adjust weights and fit again. The
filtered three-hour ranges are scaled using the station's K9 limit. KASM caps
combined K at the combined K from the unfiltered signal, as in the author's
program; consequently `k` can be smaller than `max(k_x, k_y)`. Component ranges
and component K columns describe the filtered residuals.

KASM uses **one-minute means**, with up to three hours from each available
neighboring day. `--context-file` supplies neighboring data for the calculator;
the pipeline uses available station months automatically. At missing neighbors,
KASM follows the author's reduced-context behavior and still calculates indices.
The summary reports how many days have reduced context. KASM does not require
longitude, an SQ file, or an annual H value when using X/Y components.

For one-minute IAGA input, provide its cadence explicitly:

```bash
python3 calculate_k_index.py station-three-days.min --input-format iaga \
  --method kasm --iaga-sample-step 60 --k9-limit 500 \
  -o /tmp/station-KASM.csv
```

For one-second IAGA input, omit `--iaga-sample-step` (default: 1). KASM averages
each UTC minute; a component needs all 60 valid seconds, otherwise that minute
is marked missing. Partial minutes are not treated as complete means. The
coverage columns and `--minimum-coverage` refer to valid minute means before
gap filling. The native routines fill interior gaps of up to six consecutive
missing minutes. Longer gaps and leading/trailing gaps invalidate touched
intervals; both components must be valid for combined K.

`--kasm-smoothing FACTOR` adjusts the default smoothing (finite and positive;
default 1). In pipeline configuration, set `"method": "kasm"` and optionally
`"kasm_smoothing": 1.0` in station entries or defaults. Pipeline outputs use
`-KASM` filenames. SQ bootstrap/output options apply only to IH. No sources
are downloaded at runtime.

Validation compares 1,240 interval values with a separately compiled, complete
copy of the author's KASM program, covering all four example stations and two
smoothing settings. [Reference fixture provenance](tests/data/kasm-reference.md).

## Select the Hermanus LRNS method

The Hermanus method uses linear-phase, robust, non-linear smoothing (LRNS):
forward and reverse adaptive filtering estimates the quiet-day pattern, which
is subtracted before classifying three-hour X/Y residual ranges. The public
reference implementation and notes are distributed by
[ISGI](https://isgi.unistra.fr/softwares.php). The method is described by
[Hattingh, Loubser & Nagtegaal (1989)](https://doi.org/10.1111/j.1365-246X.1989.tb02038.x).

```bash
python3 calculate_k_index.py KIR/KIR23DEC.BIN --method hermanus \
  --hermanus-mu 0.17 --hermanus-q-index M \
  --hermanus-q-start 300 --hermanus-q-end 1009 \
  -o /tmp/KIR23DEC-k-index-HERMANUS.csv

python3 run_pipeline.py --method hermanus --output-dir /tmp/hermanus-results
```

Hermanus runs entirely in Python using the standard library. No compiler,
`FC` setting, native library, or third-party package is required. The port
preserves the reference routine's single-precision arithmetic and seeded
random guide, and is checked against saved standalone Fortran results.
GNU Fortran is needed only if a developer chooses to regenerate those reference
fixtures with `python3 tests/generate_hermanus_reference.py`.

**Tune the quiet-day parameters for each observatory.** The starting defaults
`mu=0.17`, `q-index=M`, `q-start=300`, and `q-end=1009` come from the supplied
Hermanus example; they are not validated settings for the included stations.
The original guidance is to compare representative quiet, moderate, and disturbed
days against hand-scaled indices, then hold the selected parameters fixed.
The calculator reports the settings used on stderr.

- `--hermanus-mu`: smoothing factor, 0.08 to 0.20 (default 0.17).
- `--hermanus-q-index {M,E}`: `M` for a quiet-day pattern in the middle of the
  UT day; `E` for a pattern spanning its ends (default `M`).
- `--hermanus-q-start`, `--hermanus-q-end`: original Fortran's **1-based minute
  indices**, with `1 <= start < end <= 1440`. For `M` they delimit the quiet-day
  variation; for `E` they delimit the intervening nighttime period.
- `--hermanus-q-factor`: reference amplitude in nT (default `0.01*K9`, the
  standard K=0 upper boundary; the original Hermanus example uses 3 nT).
- `--hermanus-seed`: deterministic guide-noise seed, 1 to 2147483646 (default 1).
  The guide noise is reset for each component/day, so file order does not change results.

Hermanus uses **one complete day of minute means** and needs no neighboring-day
files, longitude, or external SQ file. `--context-file` is rejected. For minute
IAGA input, use `--input-format iaga --iaga-sample-step 60 --k9-limit NT`.
One-second IAGA input is averaged to minute means with the same full-minute
coverage rule as KASM.

The reference filter has no missing-data support: a component day containing
any missing minute has all eight component indices unavailable; combined K
requires both components. Observed coverage is still reported per interval.
The tool does not invent samples or filter across gaps. Inputs should be free
of instrument spikes and steps, as required by the original routine.

Pipeline configuration accepts `"method": "hermanus"`, with optional keys
`hermanus_mu`, `hermanus_q_index`, `hermanus_q_start`, `hermanus_q_end`,
`hermanus_q_factor`, and `hermanus_seed` in station entries or defaults. Output
filenames end in `-HERMANUS`; IH remains the default method.

The original public-domain source is preserved for reference. The Python port
includes the forward loop accidentally commented out in the archived Fortran.
See [source and validation details](vendor/hermanus/README.md).

## Options

Run `python3 calculate_k_index.py --help` for all options.

- `--method {ih,fmi,kasm,hermanus}`: select the calculation method (default: `ih`).
- `--sq-input PATH`: use an existing 24-row SQ file.
- `--sq-output PATH`: save the SQ curve used by this run in the same 24-row format.
- `--bootstrap-file PATH`: use another IAF month file to construct the SQ curve; defaults to the month being processed.
- `--input-format {iaf,iaga}`: select IAF binary input (default) or timestamped IAGA-2002 ASCII input.
- `--bootstrap-days N`: number of quietest valid days to use for bootstrap (default: 10).
- `--sq-estimator {median,mean}`: hourly estimator used when bootstrapping an SQ curve; `median` is recommended for high-latitude stations due to the risk of outliers even on quiet days (default: `median`).
- `--k9-limit NT`: override the station K9 range; by default it is read from the IAF header.
- `--minimum-coverage FRACTION`: minimum valid sample coverage for each component and interval (default: 0.9).
- `-o PATH`, `--output PATH`: write detailed per-interval results as CSV (default: standard output).
- `--comprehensive-output PATH`: also write daily K values in the IRF `comprehensive.txt` layout.

## Output and limitations

The CSV includes X/Y deflection ranges, coverage, component K values, combined K, and comparisons with K values embedded in the IAF file when available. For IH, FMI, and Hermanus, the combined K is the larger of the separately classified X and Y indices. KASM additionally caps combined K at its unfiltered value. The comprehensive output writes eight K values per date, `-` for unavailable intervals, parenthesized partial sums, and `-1` when no interval K value is available.

The default IAF path uses one-minute means, while the original SQ/K method described in the accompanying algorithm uses one-second samples. IAGA input can retain one-second resolution; IAF-derived K values remain estimates and may not match the embedded values exactly.

## Multi-station pipeline (`run_pipeline.py`)

To run all configured stations chronologically across months (automatically using existing `YYYYMM<STA>sq` files or chaining each month's bootstrapped SQ curve into the next month):

```bash
python3 run_pipeline.py
```

Station defaults (`k9_limit`, `sq_estimator`, `bootstrap_days`, `minimum_coverage`) are configured in `stations.json`. Common pipeline flags:

- `--station CODE`: run only specific station(s), e.g. `--station KIR --station WIC`.
- `--output-dir DIR`: write outputs to a separate directory instead of the station folders.
- `--sq-estimator {median,mean}`: override the configured SQ estimator for the run.
- `--save-generated-sq`: save newly bootstrapped `YYYYMM<STA>sq` files to disk.
- `--dry-run`: preview the chronological SQ chaining and match statistics without writing files.

## Testing & CI

Run the automated unit and regression test suite using Python's standard library `unittest` runner:

```bash
python3 -m unittest discover -s tests -v
```

CI/CD configurations for both GitHub Actions (`.github/workflows/ci.yml`) and GitLab CI (`.gitlab-ci.yml`) automatically run the test suite, execute `run_pipeline.py`, and build `K-index-calculations-algorithm-revised.pdf`.

## Algorithm reference

- `K-index-calculations-algorithm.pdf`: original algorithm description.
- `K-index-calculations-algorithm-revised.pdf`: revised description with bootstrap guidance.
