"""Compare local KASM output with the transcribed public image KAsm rows.

Run from any directory: python validation/check_kasm_images.py
Prints differences; does not modify files. Reference discrepancies are findings,
so exit status is zero after a successful comparison.
"""
import hashlib
import json
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import calculate_k_index as ih
import calculate_kasm_k_index as kasm


def main():
    fixture = json.loads(Path(__file__).with_name('kasm-image-reference.json').read_text())
    records = {}
    for item in fixture['inputs']:
        path = ROOT / item['path']
        if hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
            raise ValueError(f'Input has changed: {path}')
        records.setdefault(path.parent.name, []).extend(ih.parse_month_file(path))
    days = {station: kasm.iaf_days(values) for station, values in records.items()}
    matches = relaxed_matches = total = 0
    for case in fixture['days']:
        station = case['station']
        options = dict(target_days=[date.fromisoformat(case['date'])], k9_nT=case['k9_nT'])
        actual = kasm.calculate_days(days[station], **options)
        relaxed = kasm.calculate_days(days[station], minimum_coverage=.01, **options)
        for expected, row, relaxed_row in zip(case['expected_k'], actual, relaxed):
            total += 1
            matches += expected == row['k']
            relaxed_matches += expected == relaxed_row['k']
            if expected != row['k']:
                print(f"{station} {case['date']} {row['interval_utc']}: "
                      f"image={expected}, default={row['k']}, "
                      f"coverage .01={relaxed_row['k']}")
    print(f'Default coverage: {matches}/{total} match; coverage .01: {relaxed_matches}/{total} match')


if __name__ == '__main__':
    main()
