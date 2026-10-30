#!/usr/bin/env python3
"""Estimate three-hour K indices from INTERMAGNET IAF month files.

IAF stores one-minute means, while the source method uses one-second samples;
calculated values may therefore differ from the K values embedded in the file.
An optional SQ file contains 24 hourly X/Y/Z values in pT on a UTC day. Without
one, the program bootstraps an SQ curve from the quietest days in the month
being processed, unless another IAF file is explicitly given with
--bootstrap-file.
The combined K is the larger of the separately classified X and Y ranges.
"""

from __future__ import annotations

import argparse
import bisect
import bz2
import csv
import gzip
import lzma
import math
import statistics
import struct
import sys
from datetime import date, datetime
from pathlib import Path
from typing import Iterator, TextIO


WORDS_PER_RECORD = 5888
BYTES_PER_WORD = 4
RECORD_BYTES = WORDS_PER_RECORD * BYTES_PER_WORD
MISSING_VALUE = 999999
MISSING_K_VALUE = 999
COMPONENT_START_WORDS = {"x": 16, "y": 16 + 1440}
HOURLY_START_WORDS = {"x": 5776, "y": 5800, "z": 5824}
K_START_WORD = 5876
K_THRESHOLD_FRACTIONS = (0.01, 0.02, 0.04, 0.08, 0.14, 0.24, 0.40, 0.66)
INTERVALS = tuple((hour, hour + 3) for hour in range(0, 24, 3))


def read_sq_curve(path: Path) -> dict[str, tuple[float, ...]]:
    rows: list[tuple[float, float, float]] = []
    for line_number, line in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 3:
            raise ValueError(f"{path}:{line_number}: expected three component values")
        try:
            values = tuple(float(field) for field in fields)
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: expected numeric pT values") from error
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"{path}:{line_number}: SQ values must be finite")
        rows.append(values)

    if len(rows) != 24:
        raise ValueError(f"{path}: expected 24 hourly rows, found {len(rows)}")

    return {
        component: tuple(row[index] for row in rows)
        for index, component in enumerate(("x", "y", "z"))
    }


def interpolate_sq_pt(hourly_values: tuple[float, ...], second_of_day: float) -> float:
    hour = int(second_of_day // 3600)
    second = second_of_day - hour * 3600
    fraction = second / 3600.0
    next_hour = (hour + 1) % 24
    return hourly_values[hour] + fraction * (hourly_values[next_hour] - hourly_values[hour])


def classify_range_pt(range_pt: float, k9_nT: int) -> int:
    if k9_nT <= 0:
        raise ValueError(f"invalid K9 header value: {k9_nT} nT")
    thresholds_pt = tuple(k9_nT * fraction * 1000 for fraction in K_THRESHOLD_FRACTIONS)
    if range_pt > k9_nT * 1000:
        return 9
    return bisect.bisect_right(thresholds_pt, range_pt)


def bootstrap_sq_curve(
    records: list[tuple[int, ...]],
    quiet_day_count: int = 10,
    estimator: str = "median",
) -> tuple[dict[str, tuple[float, ...]], list[str]]:
    if quiet_day_count < 1:
        raise ValueError("bootstrap day count must be at least 1")
    if estimator == "median":
        aggregate = statistics.median
    elif estimator == "mean":
        aggregate = statistics.fmean
    else:
        raise ValueError(f"unsupported SQ estimator: {estimator}")

    candidates: list[tuple[tuple[int, int, int], str, tuple[int, ...]]] = []
    for words in records:
        daily_k = words[K_START_WORD : K_START_WORD + 8]
        if any(value == MISSING_K_VALUE for value in daily_k):
            continue
        day = datetime.strptime(str(words[1]), "%Y%j").date().isoformat()
        score = (sum(daily_k), sum(value * value for value in daily_k), max(daily_k))
        candidates.append((score, day, words))

    if len(candidates) < quiet_day_count:
        raise ValueError(
            f"cannot bootstrap from {quiet_day_count} quiet days: "
            f"only {len(candidates)} records have all eight stored K values"
        )

    selected = sorted(candidates, key=lambda item: (item[0], item[1]))[:quiet_day_count]
    curve: dict[str, tuple[float, ...]] = {}
    for component, start_word in HOURLY_START_WORDS.items():
        hourly_values: list[float] = []
        for hour in range(24):
            samples = [
                words[start_word + hour] * 100.0
                for _, _, words in selected
                if words[start_word + hour] != MISSING_VALUE
            ]
            if not samples:
                raise ValueError(f"cannot bootstrap {component.upper()} SQ: hour {hour:02d} is missing")
            hourly_values.append(aggregate(samples))
        curve[component] = tuple(hourly_values)

    return curve, [day for _, day, _ in selected]


def format_number(value: float | None, digits: int = 3) -> str:
    return "" if value is None else f"{value:.{digits}f}"


def write_sq_curve(sq: dict[str, tuple[float, ...]], path: Path) -> None:
    with path.open("w", encoding="ascii", newline="\n") as output:
        for hour in range(24):
            values = (sq[component][hour] for component in ("x", "y", "z"))
            formatted = (
                str(int(value)) if value.is_integer() else format(value, ".15g")
                for value in values
            )
            output.write(" ".join(formatted) + "\n")


def calculate_day(
    words: tuple[int, ...],
    sq: dict[str, tuple[float, ...]],
    minimum_coverage: float,
    k9_override_nT: float | None = None,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    day = datetime.strptime(str(words[1]), "%Y%j").date()
    k9_nT = words[10] if k9_override_nT is None else k9_override_nT

    for interval_index, (start_hour, end_hour) in enumerate(INTERVALS):
        start_minute = start_hour * 60
        end_minute = end_hour * 60
        ranges_pt: dict[str, float | None] = {}
        coverage: dict[str, float] = {}
        indices: dict[str, int | None] = {}

        for component in ("x", "y"):
            values_pt: list[float] = []
            data_start = COMPONENT_START_WORDS[component]
            for minute_of_day in range(start_minute, end_minute):
                raw_value = words[data_start + minute_of_day]
                if raw_value == MISSING_VALUE:
                    continue
                observed_pt = raw_value * 100.0  # IAF component words are in 0.1 nT.
                quiet_pt = interpolate_sq_pt(sq[component], minute_of_day * 60)
                values_pt.append(observed_pt - quiet_pt)

            coverage[component] = len(values_pt) / (end_minute - start_minute)
            if coverage[component] < minimum_coverage or not values_pt:
                ranges_pt[component] = None
                indices[component] = None
            else:
                component_range_pt = max(values_pt) - min(values_pt)
                ranges_pt[component] = component_range_pt
                indices[component] = classify_range_pt(component_range_pt, k9_nT)

        combined_index = (
            max(indices["x"], indices["y"])
            if indices["x"] is not None and indices["y"] is not None
            else None
        )
        stored_raw = words[K_START_WORD + interval_index]
        stored_index = None if stored_raw == MISSING_K_VALUE else stored_raw / 10.0

        rows.append(
            {
                "date": day.isoformat(),
                "interval_utc": f"{start_hour:02d}:00-{end_hour:02d}:00",
                "x_range_nT": None if ranges_pt["x"] is None else ranges_pt["x"] / 1000.0,
                "y_range_nT": None if ranges_pt["y"] is None else ranges_pt["y"] / 1000.0,
                "x_coverage": coverage["x"],
                "y_coverage": coverage["y"],
                "k_x": indices["x"],
                "k_y": indices["y"],
                "k": combined_index,
                "stored_k_raw": stored_raw,
                "stored_k": stored_index,
                "matches_stored": (
                    "" if combined_index is None or stored_index is None
                    else int(combined_index == stored_index)
                ),
            }
        )

    return rows


def open_iaga_text(path: Path) -> TextIO:
    with path.open("rb") as raw:
        magic = raw.read(6)
    suffix = path.suffix.lower()
    if suffix in (".gz", ".gzip") or magic.startswith(b"\x1f\x8b"):
        return gzip.open(path, "rt", encoding="ascii")
    if suffix in (".bz2", ".bzip2") or magic.startswith(b"BZh"):
        return bz2.open(path, "rt", encoding="ascii")
    if suffix in (".xz", ".lzma") or magic.startswith(b"\xfd7zXZ\x00"):
        return lzma.open(path, "rt", encoding="ascii")
    return path.open("r", encoding="ascii")


def iter_iaga_days(
    path: Path,
) -> Iterator[tuple[date, list[tuple[float, float | None, float | None]]]]:
    component_columns: dict[str, int] | None = None
    current_day: date | None = None
    day_samples: list[tuple[float, float | None, float | None]] = []
    previous_timestamp: datetime | None = None

    try:
        with open_iaga_text(path) as source:
            for line_number, line in enumerate(source, 1):
                fields = line.replace("|", " ").split()
                if not fields:
                    continue
                if fields[0].upper() == "DATE" and len(fields) >= 3 and fields[2].upper() == "DOY":
                    component_columns = {
                        component: index
                        for index, field in enumerate(fields[3:], 3)
                        for component in ("x", "y")
                        if field.upper().endswith(component.upper())
                    }
                    if set(component_columns) != {"x", "y"}:
                        raise ValueError(f"{path}:{line_number}: IAGA header must contain X and Y columns")
                    continue
                if not fields[0][0:1].isdigit() or len(fields[0]) != 10 or fields[0][4] != "-":
                    continue
                if component_columns is None:
                    raise ValueError(f"{path}:{line_number}: IAGA data encountered before the DATE/TIME/DOY header")
                try:
                    timestamp = datetime.fromisoformat(f"{fields[0]}T{fields[1]}")
                    values = []
                    for component in ("x", "y"):
                        value = float(fields[component_columns[component]])
                        values.append(None if not math.isfinite(value) or abs(value) in (99999, 999999) else value)
                except (IndexError, ValueError) as error:
                    raise ValueError(f"{path}:{line_number}: invalid IAGA timestamp or component value") from error

                if previous_timestamp is not None and timestamp <= previous_timestamp:
                    raise ValueError(f"{path}:{line_number}: IAGA timestamps must be strictly increasing")
                previous_timestamp = timestamp
                sample_day = timestamp.date()
                if current_day is not None and sample_day != current_day:
                    if sample_day < current_day:
                        raise ValueError(f"{path}:{line_number}: IAGA dates must be in chronological order")
                    yield current_day, day_samples
                    day_samples = []
                current_day = sample_day
                seconds_of_day = (
                    timestamp.hour * 3600
                    + timestamp.minute * 60
                    + timestamp.second
                    + timestamp.microsecond / 1_000_000
                )
                day_samples.append((seconds_of_day, values[0], values[1]))
    except (EOFError, lzma.LZMAError, gzip.BadGzipFile) as error:
        raise ValueError(f"{path}: invalid compressed IAGA file: {error}") from error

    if component_columns is None:
        raise ValueError(f"{path}: missing IAGA DATE/TIME/DOY header with X and Y columns")
    if current_day is None:
        raise ValueError(f"{path}: no IAGA data rows found")
    yield current_day, day_samples


def calculate_iaga_day(
    day: date,
    samples: list[tuple[float, float | None, float | None]],
    sq: dict[str, tuple[float, ...]],
    minimum_coverage: float,
    k9_nT: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for start_hour, end_hour in INTERVALS:
        start_second = start_hour * 3600
        end_second = end_hour * 3600
        ranges_pt: dict[str, float | None] = {}
        coverage: dict[str, float] = {}
        indices: dict[str, int | None] = {}

        for component_index, component in enumerate(("x", "y"), 1):
            values_pt = [
                observed_nT * 1000.0 - interpolate_sq_pt(sq[component], second_of_day)
                for second_of_day, x_nT, y_nT in samples
                if start_second <= second_of_day < end_second
                for observed_nT in ((x_nT, y_nT)[component_index - 1],)
                if observed_nT is not None
            ]
            coverage[component] = len(values_pt) / (end_second - start_second)
            if coverage[component] < minimum_coverage or not values_pt:
                ranges_pt[component] = None
                indices[component] = None
            else:
                component_range_pt = max(values_pt) - min(values_pt)
                ranges_pt[component] = component_range_pt
                indices[component] = classify_range_pt(component_range_pt, k9_nT)

        combined_index = (
            max(indices["x"], indices["y"])
            if indices["x"] is not None and indices["y"] is not None
            else None
        )
        rows.append(
            {
                "date": day.isoformat(),
                "interval_utc": f"{start_hour:02d}:00-{end_hour:02d}:00",
                "x_range_nT": None if ranges_pt["x"] is None else ranges_pt["x"] / 1000.0,
                "y_range_nT": None if ranges_pt["y"] is None else ranges_pt["y"] / 1000.0,
                "x_coverage": coverage["x"],
                "y_coverage": coverage["y"],
                "k_x": indices["x"],
                "k_y": indices["y"],
                "k": combined_index,
                "stored_k_raw": None,
                "stored_k": None,
                "matches_stored": "",
            }
        )
    return rows


def parse_month_file(path: Path) -> list[tuple[int, ...]]:
    data = path.read_bytes()
    if not data or len(data) % RECORD_BYTES:
        raise ValueError(
            f"{path}: size {len(data)} is not a whole number of {RECORD_BYTES}-byte IAF records"
        )

    records = [
        struct.unpack_from(f"<{WORDS_PER_RECORD}i", data, offset)
        for offset in range(0, len(data), RECORD_BYTES)
    ]
    expected_station = path.stem[:3].upper()
    for record_number in range(len(records)):
        station = data[record_number * RECORD_BYTES : record_number * RECORD_BYTES + 4]
        if station.strip().decode("ascii", errors="replace") != expected_station:
            raise ValueError(
                f"{path}: record {record_number + 1} has station ID {station!r}, "
                f"not {expected_station}"
            )
    return records


def write_csv(rows: list[dict[str, object]], output: TextIO) -> None:
    columns = (
        "date",
        "interval_utc",
        "x_range_nT",
        "y_range_nT",
        "x_coverage",
        "y_coverage",
        "k_x",
        "k_y",
        "k",
        "stored_k_raw",
        "stored_k",
        "matches_stored",
    )
    writer = csv.DictWriter(output, fieldnames=columns)
    writer.writeheader()
    for row in rows:
        formatted = dict(row)
        formatted["x_range_nT"] = format_number(row["x_range_nT"])
        formatted["y_range_nT"] = format_number(row["y_range_nT"])
        formatted["x_coverage"] = format_number(row["x_coverage"], 3)
        formatted["y_coverage"] = format_number(row["y_coverage"], 3)
        if row["stored_k"] is not None:
            formatted["stored_k"] = format_number(row["stored_k"], 1)
        writer.writerow(formatted)


def write_comprehensive(rows: list[dict[str, object]], output: TextIO) -> None:
    output.write(" YYYYMMDD K[8]            Sum [-1 no value]\n")
    days: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        days.setdefault(str(row["date"]), []).append(row)

    for day, day_rows in days.items():
        day_rows.sort(key=lambda row: str(row["interval_utc"]))
        indices = [row["k"] for row in day_rows]
        values = ["-" if index is None else str(index) for index in indices]
        available = [index for index in indices if index is not None]
        if not available:
            total = "-1"
        elif len(available) < 8:
            total = f"({sum(available)})"
        else:
            total = str(sum(available))
        output.write(f"{day.replace('-', '')} {' '.join(values)} {total:>3}\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Calculate local three-hour K indices from an IAF or IAGA data file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  Bootstrap SQ and save it: python calculate_k_index.py WIC/WIC23DEC.BIN --sq-output WIC/202312WICsq\n"
            "  Use an SQ file: python calculate_k_index.py WIC/WIC23DEC.BIN --sq-input WIC/202311WICsq\n"
            "  Use IAGA 1-second data: python calculate_k_index.py KIR/kir202609qsec.sec.gz --input-format iaga --sq-input KIR/202310KIRsq --k9-limit 1500\n"
            "  Set K9 to 500 nT: python calculate_k_index.py WIC/WIC23DEC.BIN --k9-limit 500\n\n"
            "SQ files contain 24 rows of X, Y, Z values in pT. Plain-text and gzip/bzip2/xz-\n"
            "compressed IAGA files are supported. IAGA input requires an SQ source\n"
            "and --k9-limit because neither value is supplied by the data file. The positional\n"
            "SQ-file argument is still accepted for compatibility; --sq-input is preferred."
        ),
    )
    parser.add_argument("input_file", nargs="?", type=Path, default=Path("KIR/KIR23DEC.BIN"))
    parser.add_argument("sq_file", nargs="?", type=Path, help="legacy positional SQ file in pT")
    parser.add_argument(
        "--input-format",
        choices=("iaf", "iaga"),
        default="iaf",
        help="input format (default: iaf)",
    )
    parser.add_argument("--method", choices=("ih", "fmi", "kasm"), default="ih",
                        help="K-index method (default: ih); FMI requires C; KASM requires C++11")
    parser.add_argument("--kasm-smoothing", type=float, default=1.0,
                        help="KASM smoothing factor (finite and positive; default: 1)")
    parser.add_argument("--iaga-sample-step", type=int, choices=(1, 60), default=1,
                        help="KASM IAGA input cadence in seconds (default: 1; minute input: 60)")
    parser.add_argument("--longitude", type=float,
                        help="observatory longitude in degrees east (FMI; defaults to IAF header)")
    parser.add_argument("--context-file", type=Path, action="append", default=[],
                        help="additional same-station input file for FMI/KASM adjacent-day context; repeatable")
    parser.add_argument("--sq-input", type=Path, help="existing 24-row SQ file in pT")
    parser.add_argument("--sq-output", type=Path, help="write the SQ curve used by this run in pT")
    parser.add_argument(
        "--bootstrap-file",
        type=Path,
        help="IAF month file used to create an SQ curve (default: IAF input file)",
    )
    parser.add_argument("--bootstrap-days", type=int, default=10, help="quietest days to use (default: 10)")
    parser.add_argument(
        "--sq-estimator",
        choices=("median", "mean"),
        default="median",
        help="hourly estimator when bootstrapping an SQ curve; median is recommended for high-latitude stations due to outlier risk on quiet days (default: median)",
    )
    parser.add_argument("-o", "--output", type=Path, help="CSV output path (default: stdout)")
    parser.add_argument(
        "--comprehensive-output",
        type=Path,
        help="also write daily K values in IRF comprehensive.txt format",
    )
    parser.add_argument(
        "--minimum-coverage",
        type=float,
        default=0.9,
        help="minimum valid sample fraction per component and interval (default: 0.9)",
    )
    parser.add_argument(
        "--k9-limit",
        type=float,
        help="override the IAF header K9 range in nT (must be finite and positive)",
    )
    args = parser.parse_args()

    if not math.isfinite(args.kasm_smoothing) or args.kasm_smoothing <= 0:
        parser.error("--kasm-smoothing must be finite and positive")
    if args.method != "kasm" and (args.kasm_smoothing != 1.0 or args.iaga_sample_step != 1):
        parser.error("--kasm-smoothing and --iaga-sample-step apply to --method kasm")
    if not 0 < args.minimum_coverage <= 1:
        parser.error("--minimum-coverage must be greater than 0 and at most 1")
    if args.k9_limit is not None and (not math.isfinite(args.k9_limit) or args.k9_limit <= 0):
        parser.error("--k9-limit must be finite and greater than 0")
    if args.input_format == "iaga" and args.k9_limit is None:
        parser.error("--k9-limit is required with --input-format iaga")
    if args.sq_input and args.sq_file:
        parser.error("use either --sq-input or the positional SQ file, not both")
    if (
        args.method == "ih"
        and args.input_format == "iaga"
        and not args.sq_input
        and not args.sq_file
        and not args.bootstrap_file
    ):
        parser.error("IAGA input requires --sq-input or an IAF --bootstrap-file")
    if args.output and args.comprehensive_output and args.output.resolve() == args.comprehensive_output.resolve():
        parser.error("CSV and comprehensive outputs must use different paths")
    if args.sq_output and any(
        path and args.sq_output.resolve() == path.resolve()
        for path in (args.output, args.comprehensive_output)
    ):
        parser.error("SQ output must use a different path from K-index outputs")

    if args.method in ("fmi", "kasm"):
        inputs = [args.input_file.resolve(), *(p.resolve() for p in args.context_file)]
        if any(p and p.resolve() in inputs for p in (args.output, args.comprehensive_output)):
            parser.error("outputs must not overwrite input or context files")
    try:
        if args.method in ("fmi", "kasm"):
            if args.method == "fmi":
                import calculate_fmi_k_index as backend
            else:
                import calculate_kasm_k_index as backend
            label = args.method.upper()
            if args.sq_input or args.sq_file or args.sq_output or args.bootstrap_file:
                parser.error(f"{label} fits its own daily quiet curve; SQ input/output and bootstrap files apply to IH")
            if args.method == "kasm" and args.longitude is not None:
                parser.error("--longitude applies to FMI; KASM uses no longitude parameter")
            if args.method == "fmi" and args.input_format == "iaga" and args.longitude is None:
                parser.error("--longitude is required for FMI with IAGA input")
            def load_days(path):
                if args.input_format == "iaf":
                    return backend.iaf_days(parse_month_file(path))
                if args.method == "kasm":
                    return backend.iaga_days(path, sample_step=args.iaga_sample_step)
                return backend.iaga_days(path)
            days = load_days(args.input_file)
            targets = list(days)
            for context_path in args.context_file:
                context = load_days(context_path)
                if args.input_format == "iaf":
                    identity = next(iter(days.values()))[2]
                    if any((data[2][0], data[2][2], data[2][3]) !=
                           (identity[0], identity[2], identity[3]) for data in context.values()):
                        raise ValueError(f"{label} context must belong to the same station and coordinates")
                for day, data in context.items():
                    if day in days and days[day] != data:
                        raise ValueError(f"conflicting {label} context for {day}")
                    days[day] = data
            options = (dict(step=60 if args.input_format == "iaf" else 1, longitude=args.longitude)
                       if args.method == "fmi" else dict(smoothing=args.kasm_smoothing))
            rows = backend.calculate_days(days, targets, k9_nT=args.k9_limit,
                minimum_coverage=args.minimum_coverage, **options)
            processed_days = len(targets)
            from datetime import timedelta
            missing_context = sum(any(day + timedelta(days=n) not in days for n in (-1, 1))
                                  for day in targets)
            context_status = ("unavailable due to missing adjacent-day context"
                              if args.method == "fmi" else "calculated with reduced adjacent-day context")
            bootstrap_description = f"{label} method; {missing_context} day(s) {context_status}"
        else:
            if args.context_file or args.longitude is not None:
                parser.error("--context-file applies to FMI/KASM; --longitude applies to FMI")
            records = parse_month_file(args.input_file) if args.input_format == "iaf" else None
            sq_input = args.sq_input or args.sq_file
            bootstrap_path: Path | None = None
            if sq_input:
                sq = read_sq_curve(sq_input)
                bootstrap_description = f"SQ file {sq_input}"
            else:
                bootstrap_path = args.bootstrap_file or args.input_file
                if records is not None and bootstrap_path.resolve() == args.input_file.resolve():
                    bootstrap_records = records
                else:
                    bootstrap_records = parse_month_file(bootstrap_path)
                sq, quiet_days = bootstrap_sq_curve(
                    bootstrap_records, args.bootstrap_days, estimator=args.sq_estimator
                )
                bootstrap_description = (
                    f"bootstrapped SQ ({args.sq_estimator}) from {bootstrap_path} using {len(quiet_days)} quiet days: "
                    f"{', '.join(quiet_days)}"
                )
            if args.sq_output and args.sq_output.resolve() == args.input_file.resolve():
                parser.error("SQ output must not overwrite the input file")
            if args.sq_output and bootstrap_path and args.sq_output.resolve() == bootstrap_path.resolve():
                parser.error("SQ output must not overwrite the IAF bootstrap file")
            if args.input_format == "iaga":
                rows = []
                processed_days = 0
                for day, samples in iter_iaga_days(args.input_file):
                    rows.extend(
                        calculate_iaga_day(
                            day, samples, sq, args.minimum_coverage, args.k9_limit
                        )
                    )
                    processed_days += 1
            else:
                rows = [
                    row
                    for record in records
                    for row in calculate_day(
                        record, sq, args.minimum_coverage, k9_override_nT=args.k9_limit
                    )
                ]
                processed_days = len(records)
    except (OSError, ValueError, struct.error) as error:
        parser.error(str(error))

    if args.sq_output:
        write_sq_curve(sq, args.sq_output)

    output: TextIO
    if args.output:
        output = args.output.open("w", newline="", encoding="utf-8")
    else:
        output = sys.stdout
    try:
        write_csv(rows, output)
    finally:
        if args.output:
            output.close()

    if args.comprehensive_output:
        with args.comprehensive_output.open("w", encoding="ascii", newline="\n") as comprehensive:
            write_comprehensive(rows, comprehensive)

    if args.input_format == "iaf":
        comparisons = [row for row in rows if row["matches_stored"] != ""]
        matches = sum(row["matches_stored"] == 1 for row in comparisons)
        print(
            f"Processed {processed_days} days and {len(rows)} three-hour intervals; "
            f"{matches}/{len(comparisons)} calculated combined K values match the stored values.",
            file=sys.stderr,
        )
    else:
        print(
            f"Processed {processed_days} days and {len(rows)} three-hour intervals from IAGA data.",
            file=sys.stderr,
        )
    print(bootstrap_description, file=sys.stderr)
    if args.k9_limit is not None:
        print(f"K9 range override: {args.k9_limit:g} nT", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
