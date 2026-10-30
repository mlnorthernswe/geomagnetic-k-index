#!/usr/bin/env python3
"""Operational multi-station, multi-month K-index data processing pipeline.

Discovers monthly IAF files for configured geomagnetic observatories, sorts
them chronologically per station, chains each month's solar-quiet (SQ) curve
into the following month (or uses an existing YYYYMM<STA>sq file when present),
and writes per-month CSV and IRF comprehensive K-index outputs.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import calculate_k_index as ih


@dataclass(frozen=True)
class MonthInput:
    station: str
    year: int
    month: int
    path: Path
    records: list[tuple[int, ...]]

    @property
    def yyyymm(self) -> str:
        return f"{self.year:04d}{self.month:02d}"

    @property
    def prev_yyyymm(self) -> str:
        if self.month == 1:
            return f"{self.year - 1:04d}12"
        return f"{self.year:04d}{self.month - 1:02d}"


@dataclass(frozen=True)
class MonthRunSummary:
    station: str
    yyyymm: str
    input_path: Path
    sq_source: str
    days: int
    intervals: int
    matched_intervals: int
    comparable_intervals: int
    csv_output: Path
    comprehensive_output: Path
    sq_output: Path | None


def discover_station_months(station: str, directory: Path) -> list[MonthInput]:
    """Find all IAF .bin/.BIN files in directory and sort chronologically."""
    if not directory.is_dir():
        raise FileNotFoundError(f"station directory not found: {directory}")

    months: list[MonthInput] = []
    for candidate in sorted(directory.iterdir()):
        if not candidate.is_file() or candidate.suffix.lower() != ".bin":
            continue
        records = ih.parse_month_file(candidate)
        first_day = datetime.strptime(str(records[0][1]), "%Y%j").date()
        months.append(
            MonthInput(
                station=station.upper(),
                year=first_day.year,
                month=first_day.month,
                path=candidate,
                records=records,
            )
        )

    months.sort(key=lambda item: (item.year, item.month, item.path.name))
    return months


def run_station_pipeline(
    station: str,
    station_cfg: dict[str, object],
    defaults: dict[str, object],
    base_dir: Path,
    output_dir: Path | None = None,
    estimator_override: str | None = None,
    save_generated_sq: bool = False,
    dry_run: bool = False,
    method_override: str | None = None,
) -> list[MonthRunSummary]:
    station_code = station.upper()
    station_dir = base_dir / str(station_cfg.get("directory", station_code))
    target_dir = (output_dir / station_code) if output_dir else station_dir
    if not dry_run:
        target_dir.mkdir(parents=True, exist_ok=True)

    estimator = str(
        estimator_override
        or station_cfg.get("sq_estimator")
        or defaults.get("sq_estimator", "median")
    )
    bootstrap_days = int(
        station_cfg.get("bootstrap_days", defaults.get("bootstrap_days", 10))
    )
    minimum_coverage = float(
        station_cfg.get("minimum_coverage", defaults.get("minimum_coverage", 0.9))
    )
    k9_limit_raw = station_cfg.get("k9_limit", defaults.get("k9_limit"))
    k9_limit = float(k9_limit_raw) if k9_limit_raw is not None else None

    months = discover_station_months(station_code, station_dir)
    method = str(method_override or station_cfg.get("method") or defaults.get("method", "ih"))
    if method not in ("ih", "fmi", "kasm"):
        raise ValueError(f"unknown K-index method: {method}")
    if method in ("fmi", "kasm"):
        if save_generated_sq:
            raise ValueError("--save-generated-sq applies to IH; FMI/KASM fit daily quiet curves")
        return run_adaptive_months(months, target_dir, k9_limit, minimum_coverage,
                              station_cfg.get("longitude", defaults.get("longitude")), dry_run,
                              method=method,
                              smoothing=float(station_cfg.get("kasm_smoothing", defaults.get("kasm_smoothing", 1.0))))
    generated_sq_cache: dict[str, tuple[dict[str, tuple[float, ...]], str]] = {}
    summaries: list[MonthRunSummary] = []

    for month_input in months:
        prev_sq_name = f"{month_input.prev_yyyymm}{station_code}sq"
        prev_sq_candidates = [
            target_dir / prev_sq_name,
            station_dir / prev_sq_name,
        ]
        existing_prev_sq = next((p for p in prev_sq_candidates if p.is_file()), None)

        if existing_prev_sq is not None:
            sq_curve = ih.read_sq_curve(existing_prev_sq)
            sq_source = f"file:{existing_prev_sq.name}"
        elif month_input.prev_yyyymm in generated_sq_cache:
            sq_curve, sq_source = generated_sq_cache[month_input.prev_yyyymm]
        else:
            sq_curve, quiet_days = ih.bootstrap_sq_curve(
                month_input.records,
                quiet_day_count=bootstrap_days,
                estimator=estimator,
            )
            sq_source = f"bootstrap:{estimator}({month_input.yyyymm},{len(quiet_days)}d)"

        rows = [
            row
            for record in month_input.records
            for row in ih.calculate_day(
                record,
                sq_curve,
                minimum_coverage=minimum_coverage,
                k9_override_nT=k9_limit,
            )
        ]

        # Prepare this month's SQ curve to chain into month m + 1 if needed.
        curr_sq_name = f"{month_input.yyyymm}{station_code}sq"
        existing_curr_sq = next(
            (p for p in (target_dir / curr_sq_name, station_dir / curr_sq_name) if p.is_file()),
            None,
        )
        sq_output_path: Path | None = None
        if existing_curr_sq is not None:
            next_sq = ih.read_sq_curve(existing_curr_sq)
            generated_sq_cache[month_input.yyyymm] = (next_sq, f"file:{existing_curr_sq.name}")
        else:
            next_sq, _ = ih.bootstrap_sq_curve(
                month_input.records,
                quiet_day_count=bootstrap_days,
                estimator=estimator,
            )
            generated_sq_cache[month_input.yyyymm] = (
                next_sq,
                f"chained:{estimator}({month_input.yyyymm})",
            )
            if save_generated_sq:
                sq_output_path = target_dir / curr_sq_name
                if not dry_run:
                    ih.write_sq_curve(next_sq, sq_output_path)

        stem = month_input.path.stem.upper()
        csv_path = target_dir / f"{stem}-k-index-IH.csv"
        comp_path = target_dir / f"{stem}-comprehensive-IH.txt"

        if not dry_run:
            with csv_path.open("w", newline="", encoding="utf-8") as csv_file:
                ih.write_csv(rows, csv_file)
            with comp_path.open("w", encoding="ascii", newline="\n") as comp_file:
                ih.write_comprehensive(rows, comp_file)

        comparisons = [row for row in rows if row["matches_stored"] != ""]
        matches = sum(row["matches_stored"] == 1 for row in comparisons)

        summaries.append(
            MonthRunSummary(
                station=station_code,
                yyyymm=month_input.yyyymm,
                input_path=month_input.path,
                sq_source=sq_source,
                days=len(month_input.records),
                intervals=len(rows),
                matched_intervals=matches,
                comparable_intervals=len(comparisons),
                csv_output=csv_path,
                comprehensive_output=comp_path,
                sq_output=sq_output_path,
            )
        )

    return summaries


def run_adaptive_months(months, target_dir, k9_limit, minimum_coverage, longitude, dry_run,
                   method="fmi", smoothing=1.0):
    """Use all available station months as context, without bridging date gaps."""
    if method == "fmi":
        import calculate_fmi_k_index as backend
    else:
        import calculate_kasm_k_index as backend
    from datetime import timedelta

    days = backend.iaf_days([record for month in months for record in month.records])
    if not days:
        return []
    options = (dict(longitude=None if longitude is None else float(longitude))
               if method == "fmi" else dict(smoothing=smoothing))
    rows = backend.calculate_days(days, k9_nT=k9_limit,
        minimum_coverage=minimum_coverage, **options)
    by_date = {}
    for row in rows:
        by_date.setdefault(row["date"], []).append(row)
    summaries = []
    for month in months:
        targets = [datetime.strptime(str(record[1]), "%Y%j").date() for record in month.records]
        month_rows = [row for day in targets for row in by_date[day.isoformat()]]
        stem = month.path.stem.upper()
        csv_path = target_dir / f"{stem}-k-index-{method.upper()}.csv"
        comp_path = target_dir / f"{stem}-comprehensive-{method.upper()}.txt"
        if not dry_run:
            with csv_path.open("w", newline="", encoding="utf-8") as output:
                ih.write_csv(month_rows, output)
            with comp_path.open("w", encoding="ascii", newline="\n") as output:
                ih.write_comprehensive(month_rows, output)
        comparisons = [row for row in month_rows if row["matches_stored"] != ""]
        missing = sum(any(day + timedelta(days=n) not in days for n in (-1, 1)) for day in targets)
        summaries.append(MonthRunSummary(
            station=month.station, yyyymm=month.yyyymm, input_path=month.path,
            sq_source=(f"FMI daily fit ({missing}d no context)" if method == "fmi"
                       else f"KASM fit ({missing}d reduced context)"), days=len(targets),
            intervals=len(month_rows), matched_intervals=sum(r["matches_stored"] == 1 for r in comparisons),
            comparable_intervals=len(comparisons), csv_output=csv_path,
            comprehensive_output=comp_path, sq_output=None))
    return summaries


def load_config(config_path: Path) -> tuple[dict[str, object], dict[str, dict[str, object]]]:
    raw = json.loads(config_path.read_text(encoding="utf-8"))
    defaults = dict(raw.get("defaults", {}))
    stations = {str(k).upper(): dict(v) for k, v in raw.get("stations", {}).items()}
    return defaults, stations


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run the multi-station, chronological K-index pipeline."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("stations.json"),
        help="path to station configuration JSON (default: stations.json)",
    )
    parser.add_argument(
        "--station",
        action="append",
        dest="stations",
        help="station code(s) to process (default: all stations in config)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="optional root directory for generated outputs (default: each station directory)",
    )
    parser.add_argument(
        "--sq-estimator",
        choices=("median", "mean"),
        help="override the SQ estimator for all stations in this run",
    )
    parser.add_argument(
        "--save-generated-sq",
        action="store_true",
        help="write newly bootstrapped YYYYMM<STA>sq files to the output directory",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="compute and report pipeline steps without writing files",
    )
    parser.add_argument("--method", choices=("ih", "fmi", "kasm"),
                        help="override the K-index method (default: config or ih)")
    args = parser.parse_args(argv)

    base_dir = args.config.resolve().parent
    defaults, stations_cfg = load_config(args.config)

    selected_stations = (
        [s.upper() for s in args.stations]
        if args.stations
        else list(stations_cfg.keys())
    )
    unknown = [s for s in selected_stations if s not in stations_cfg]
    if unknown:
        parser.error(f"unknown station(s) {', '.join(unknown)}; available: {', '.join(stations_cfg)}")

    all_summaries: list[MonthRunSummary] = []
    try:
        for station in selected_stations:
            all_summaries.extend(
                run_station_pipeline(
                    station=station,
                    station_cfg=stations_cfg[station],
                    defaults=defaults,
                    base_dir=base_dir,
                    output_dir=args.output_dir,
                    estimator_override=args.sq_estimator,
                    save_generated_sq=args.save_generated_sq,
                    dry_run=args.dry_run,
                    method_override=args.method,
                )
            )

    except (OSError, ValueError) as error:
        parser.error(str(error))

    header = f"{'STA':<5} {'YYYYMM':<7} {'DAYS':>4} {'MATCH':>11} {'PCT':>7}  {'SQ SOURCE':<30} {'CSV OUTPUT'}"
    print(header)
    print("-" * len(header))
    for item in all_summaries:
        pct = (
            100.0 * item.matched_intervals / item.comparable_intervals
            if item.comparable_intervals
            else 0.0
        )
        match_str = f"{item.matched_intervals}/{item.comparable_intervals}"
        print(
            f"{item.station:<5} {item.yyyymm:<7} {item.days:>4} {match_str:>11} {pct:>6.1f}%  "
            f"{item.sq_source:<30} {item.csv_output}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

