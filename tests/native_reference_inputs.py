"""Deterministic station and synthetic inputs for the native port fixtures."""
import math
from pathlib import Path

import calculate_k_index as ih

ROOT = Path(__file__).resolve().parent.parent


def reference_inputs(case):
    method = case["method"]
    if "path" in case:
        records = ih.parse_month_file(ROOT / case["path"])
        i = case["day_index"]
        triplet = records[i-1:i+2]
        inputs = []
        for start in (16, 1456):
            parts = [list(record[start:start+1440]) for record in triplet]
            values = sum(parts, []) if method == "fmi" else parts[0][-180:]+parts[1]+parts[2][:180]
            inputs.append([None if v == ih.MISSING_VALUE else v if method == "fmi" else v/10 for v in values])
        return inputs
    variant = case["variant"]
    step = case["args"][0] if method == "fmi" else 60
    n = 3*86400//step if method == "fmi" else (1440,1620,1800)[variant%3]
    x = [int(100000+800*math.sin(i*step/10000)+400*math.cos(i*step/700)) for i in range(n)]
    y = [int(20000+300*math.cos(i*step/9000)) for i in range(n)]
    gaps = {
        1: [(0,2),(n-3,n)],
        2: [(n//2,n//2+6)],
        3: [(n//2,n//2+15)],
        4: [(n//2,n//2+1),(n//2+2,n//2+5),(n//2+6,n//2+7)],
        5: [(0,n)],
        6: [(179,181),(360,367),(n-1,n)],
    }.get(variant, [])
    for start, end in gaps:
        x[start:end] = [None]*(end-start)
    if method == "kasm":
        x = [None if v is None else v/10 for v in x]
        y = [v/10 for v in y]
    return x, y
