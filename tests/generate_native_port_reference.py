"""Refresh native diagnostic expectations; compilers are for developers only.

Run from the repository root with python3 tests/generate_native_port_reference.py.
Uses bundled reference sources and the existing fixture's input cases. No network.
"""
import ctypes
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
from native_reference_inputs import reference_inputs


def main():
    path = ROOT / "tests/data/native-port-reference.json"
    cases = json.loads(path.read_text())
    with tempfile.TemporaryDirectory(prefix="native-port-reference-") as tmp:
        temp = Path(tmp)
        for method in ("fmi", "kasm"):
            source = ROOT / "vendor" / method
            library = temp / (method+".so")
            if method == "fmi":
                command = [*shlex.split(os.environ.get("CC", "cc")), "-std=gnu89", "-O2", "-shared", "-fPIC",
                           str(source / "wrapper.c"), "-lm", "-o", str(library)]
            else:
                command = [*shlex.split(os.environ.get("CXX", "c++")), "-std=c++11", "-O2", "-shared", "-fPIC",
                           str(source / "casm.cpp"), str(source / "wrapper.cpp"), "-o", str(library)]
            subprocess.run(command, check=True)
            lib = ctypes.CDLL(str(library))
            lp, dp, ip = (ctypes.POINTER(t) for t in (ctypes.c_long, ctypes.c_double, ctypes.c_int))
            if method == "fmi":
                lib.fmi_details.argtypes = [lp,lp,*([ctypes.c_long]*4),*([lp]*5)]
                lib.fmi_details.restype = None
            else:
                lib.kasm_details.argtypes = [ctypes.c_int,ctypes.c_int,dp,dp,ctypes.c_double,ctypes.c_double,
                                             ip,ip,ip,dp,dp,ip,ip]
                lib.kasm_details.restype = ctypes.c_int
            for case in cases:
                if case["method"] != method:
                    continue
                if "path" in case:
                    assert hashlib.sha256((ROOT / case["path"]).read_bytes()).hexdigest() == case["sha256"]
                x, y = reference_inputs(case)
                if method == "fmi":
                    missing = 2**50
                    array = ctypes.c_long * len(x)
                    cx, cy = [array(*(missing if v is None else int(v) for v in a)) for a in (x,y)]
                    step,k9,longitude = case["args"]
                    outputs = [(ctypes.c_long*8)() for _ in range(5)]
                    lib.fmi_details(cx,cy,step,round(k9*10),round(longitude)%360,missing,*outputs)
                else:
                    array = ctypes.c_double * len(x)
                    cx,cy = [array(*(1e30 if v is None else v for v in a)) for a in (x,y)]
                    ks,kx,ky,fx,fy = [(ctypes.c_int*8)() for _ in range(5)]
                    rx,ry = (ctypes.c_double*8)(), (ctypes.c_double*8)()
                    status = lib.kasm_details(len(x),case["args"][0],cx,cy,*case["args"][1:],ks,kx,ky,rx,ry,fx,fy)
                    assert status == 0
                    outputs = [ks,kx,ky,rx,ry,fx,fy]
                case["expected"] = [list(values) for values in outputs]
    path.write_text(json.dumps(cases, indent=2)+"\n")


if __name__ == "__main__":
    main()
