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
import csv
import math
import statistics
import struct
import sys
from datetime import datetime
from pathlib import Path
from typing import TextIO


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


def interpolate_sq_pt(hourly_values: tuple[float, ...], minute_of_day: int) -> float:
    hour, minute = divmod(minute_of_day, 60)
    fraction = minute / 60.0
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
    records: list[tuple[int, ...]], quiet_day_count: int = 10
) -> tuple[dict[str, tuple[float, ...]], list[str]]:
    if quiet_day_count < 1:
        raise ValueError("bootstrap day count must be at least 1")

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
            hourly_values.append(statistics.median(samples))
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
                quiet_pt = interpolate_sq_pt(sq[component], minute_of_day)
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
                "k_combined": combined_index,
                "stored_k_raw": stored_raw,
                "stored_k": stored_index,
                "matches_stored": (
                    "" if combined_index is None or stored_index is None
                    else int(combined_index == stored_index)
                ),
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
        "k_combined",
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
        indices = [row["k_combined"] for row in day_rows]
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
        description="Calculate local three-hour K indices from an INTERMAGNET IAF month file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  Bootstrap SQ and save it: python calculate_ih_k_index.py WIC/WIC23DEC.BIN --sq-output WIC/202312WICsq\n"
            "  Use an SQ file: python calculate_ih_k_index.py WIC/WIC23DEC.BIN --sq-input WIC/202311WICsq\n"
            "  Set K9 to 500 nT: python calculate_ih_k_index.py WIC/WIC23DEC.BIN --k9-limit 500\n\n"
            "SQ files contain 24 rows of X, Y, Z values in pT. The positional SQ-file argument is\n"
            "still accepted for compatibility; --sq-input is preferred."
        ),
    )
    parser.add_argument("month_file", nargs="?", type=Path, default=Path("KIR/KIR23DEC.BIN"))
    parser.add_argument("sq_file", nargs="?", type=Path, help="legacy positional SQ file in pT")
    parser.add_argument("--sq-input", type=Path, help="existing 24-row SQ file in pT")
    parser.add_argument("--sq-output", type=Path, help="write the SQ curve used by this run in pT")
    parser.add_argument(
        "--bootstrap-file",
        type=Path,
        help="IAF month file used to create an SQ curve (default: month_file being processed)",
    )
    parser.add_argument("--bootstrap-days", type=int, default=10, help="quietest days to use (default: 10)")
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
        help="minimum valid minute fraction per component and interval (default: 0.9)",
    )
    parser.add_argument(
        "--k9-limit",
        type=float,
        help="override the IAF header K9 range in nT (must be finite and positive)",
    )
    args = parser.parse_args()

    if not 0 < args.minimum_coverage <= 1:
        parser.error("--minimum-coverage must be greater than 0 and at most 1")
    if args.k9_limit is not None and (not math.isfinite(args.k9_limit) or args.k9_limit <= 0):
        parser.error("--k9-limit must be finite and greater than 0")
    if args.sq_input and args.sq_file:
        parser.error("use either --sq-input or the positional SQ file, not both")
    if args.output and args.comprehensive_output and args.output.resolve() == args.comprehensive_output.resolve():
        parser.error("CSV and comprehensive outputs must use different paths")
    if args.sq_output and any(
        path and args.sq_output.resolve() == path.resolve()
        for path in (args.output, args.comprehensive_output)
    ):
        parser.error("SQ output must use a different path from K-index outputs")

    try:
        records = parse_month_file(args.month_file)
        sq_input = args.sq_input or args.sq_file
        bootstrap_path: Path | None = None
        if sq_input:
            sq = read_sq_curve(sq_input)
            bootstrap_description = f"SQ file {sq_input}"
        else:
            bootstrap_path = args.bootstrap_file or args.month_file
            bootstrap_records = (
                records if bootstrap_path.resolve() == args.month_file.resolve()
                else parse_month_file(bootstrap_path)
            )
            sq, quiet_days = bootstrap_sq_curve(bootstrap_records, args.bootstrap_days)
            bootstrap_description = (
                f"bootstrapped SQ from {bootstrap_path} using {len(quiet_days)} quiet days: "
                f"{', '.join(quiet_days)}"
            )
        if args.sq_output and args.sq_output.resolve() == args.month_file.resolve():
            parser.error("SQ output must not overwrite the IAF month file")
        if args.sq_output and bootstrap_path and args.sq_output.resolve() == bootstrap_path.resolve():
            parser.error("SQ output must not overwrite the IAF bootstrap file")
        rows = [
            row
            for record in records
            for row in calculate_day(
                record, sq, args.minimum_coverage, k9_override_nT=args.k9_limit
            )
        ]
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

    comparisons = [row for row in rows if row["matches_stored"] != ""]
    matches = sum(row["matches_stored"] == 1 for row in comparisons)
    print(
        f"Processed {len(records)} days and {len(rows)} three-hour intervals; "
        f"{matches}/{len(comparisons)} calculated combined K values match the stored values.",
        file=sys.stderr,
    )
    print(bootstrap_description, file=sys.stderr)
    if args.k9_limit is not None:
        print(f"K9 range override: {args.k9_limit:g} nT", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
