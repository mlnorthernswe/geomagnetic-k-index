"""Regenerate LRNS oracle fixtures with a standalone Fortran executable.

Run from the repository root: python3 tests/generate_hermanus_reference.py
Requires GNU Fortran (or an FC command). No downloads and no input mutations.
"""
import hashlib
import json
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import calculate_k_index as ih


def main():
    source = ROOT / "vendor/hermanus"
    cases = []
    with tempfile.TemporaryDirectory(prefix="lrns-reference-") as directory:
        temp = Path(directory)
        # Restore the lost DO statement independently of the adapter.
        text = (source / "K_FILTER.FOR").read_text().replace("\x1a", "")
        text = text.replace("C           do 80  i = n1, n2", "            do 80  i = n1, n2")
        (temp / "reference.f").write_text(text)
        executable = temp / "reference"
        subprocess.run([*shlex.split(os.environ.get("FC", "gfortran")), "-O0", "-std=legacy",
            "-ffp-contract=off", "-fcheck=all", str(temp / "reference.f"),
            str(source / "reference_driver.f90"), "-o", str(executable)], check=True)
        for station in ("KIR", "BEL", "WIC", "HYB"):
            path = next(p for p in (ROOT / station).iterdir() if p.name.lower() == station.lower()+"23dec.bin")
            records = ih.parse_month_file(path)
            for day_index, overrides in [(0, {}), (10, {}),
                *([(0, dict(mu=.08, qindex="E", qstart=360, qend=1080, seed=42)),
                   (0, dict(mu=.20, qfactor=100000, seed=123))] if station == "WIC" else [])]:
                record = records[day_index]
                settings = dict(mu=.17, qindex="M", qstart=300, qend=1009,
                                qfactor=record[10]*.01, seed=1)
                settings.update(overrides)
                component_results = []
                for component in ("x", "y"):
                    start = ih.COMPONENT_START_WORDS[component]
                    values = record[start:start+1440]
                    assert ih.MISSING_VALUE not in values
                    data = (f"{settings['mu']} {settings['qindex']} {settings['qfactor']} "
                            f"{settings['qstart']} {settings['qend']} {settings['seed']}\n"
                            + "\n".join(str(v/10) for v in values) + "\n")
                    result = subprocess.run([str(executable)], input=data, text=True,
                                             capture_output=True, check=True)
                    trace = [tuple(map(float, line.split())) for line in result.stdout.splitlines()]
                    residual = [raw - quiet for raw, quiet in trace]
                    ranges = [max(residual[i:i+180])-min(residual[i:i+180]) for i in range(0,1440,180)]
                    component_results.append(dict(ranges=ranges,
                        checkpoints={str(i): trace[i] for i in (0,1,179,299,359,719,1008,1079,1259,1439)}))
                cases.append(dict(path=str(path.relative_to(ROOT)),
                    input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    day_index=day_index, k9=record[10], settings=settings, components=component_results))
    fixture = dict(source="https://isgi.unistra.fr/Documents/Softwares/hermanus_fortran77.zip",
        source_sha256=hashlib.sha256((source/"K_FILTER.FOR").read_bytes()).hexdigest(), cases=cases)
    (ROOT / "tests/data/hermanus-reference.json").write_text(json.dumps(fixture, indent=2)+"\n")


if __name__ == "__main__":
    main()
