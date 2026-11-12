#!/usr/bin/env python3
"""Generate formatting fixtures (Python float repr / format / json / csv).

Usage: python -I scripts/make_fixtures.py OUT_DIR   (inside the witchi-0.2.0 env)

Writes:
  pyfmt_floats.tsv  hexbits, repr, .2f, .3f, .4f, .6f   (asserts str(np.float64) == repr)
  pyfmt_json.txt    json.dumps(indent=4) of a fixed structure (rebuilt in Rust tests)
  pyfmt_csv.txt     csv.writer(delimiter="\t") rows with awkward fields
"""

import csv
import io
import json
import math
import random
import struct
import sys
from pathlib import Path

import numpy as np


def bits(x):
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(b):
    return struct.unpack("<d", struct.pack("<Q", b))[0]


def float_corpus():
    rng = random.Random(20261112)
    vals = [0.0, -0.0, 1.0, -1.0, 0.1, 0.5, 1.5, 2.5, 0.125, 0.375, 0.0625, 1e16, 1e15,
            9999999999999998.0, 1e-4, 1e-5, 0.00011, 123456789012345678.0, 5e-324,
            2.2250738585072014e-308, 1.7976931348623157e308, math.inf, -math.inf, math.nan,
            0.0005, 0.0015, 0.0025, 0.0045, 1.0005, 2.675, 1.005, 0.06, 0.95, 0.05,
            1377.1785906880436, 4.6846241911860185, -18388.082114010555]
    for _ in range(4000):  # random bit patterns (finite)
        b = rng.getrandbits(64)
        x = from_bits(b)
        if math.isfinite(x):
            vals.append(x)
    for _ in range(4000):  # typical magnitudes
        vals.append(rng.uniform(-1, 1) * 10 ** rng.uniform(-8, 20))
    for _ in range(2000):  # short decimals and ties at 2/3/4/6 places
        d = rng.choice([2, 3, 4, 6, 7])
        vals.append(rng.randint(-10**6, 10**6) / 10**d)
        vals.append(rng.randint(0, 4096) / 2 ** rng.randint(1, 12))
    for _ in range(1000):  # p-value-like ratios
        m = rng.randint(1, 1000)
        vals.append(rng.randint(0, m) / m)
    return vals


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "pyfmt_floats.tsv", "w") as fh:
        for x in float_corpus():
            r = repr(x)
            if math.isfinite(x):
                assert str(np.float64(x)) == r, (x, str(np.float64(x)), r)
            fh.write("\t".join([f"{bits(x):016x}", r] + [format(x, f".{p}f") for p in (2, 3, 4, 6)]) + "\n")

    doc = {
        "schema_version": 3,
        "witchi_version": "0.2.0",
        "stop_reason": "delta not significant (p=0.060, max_delta=8.542055)",
        "none": None,
        "tricky": 'q"b\\s/\n\r\t\b\f\x01\x7f é ✓ 𝄞 end',
        "floats": [0.0, -0.0, 1.5, 1e-07, 1e16, math.inf, -math.inf, math.nan],
        "nested": [[1.0, 2.5], [], [3.25]],
        "empty_obj": {},
        "obj": {"k": [1, 2], "z": 0},
    }
    (out / "pyfmt_json.txt").write_text(json.dumps(doc, indent=4))

    buf = io.StringIO(newline="")
    w = csv.writer(buf, delimiter="\t")
    w.writerow(["Row", "a b", "tab\there", 'quote"inside', "nl\nx", "cr\rx", "", "plain"])
    w.writerow([""])
    w.writerow(["x", ""])
    w.writerow(["'single'", "comma,ok", "ünï"])
    (out / "pyfmt_csv.txt").write_bytes(buf.getvalue().encode())


if __name__ == "__main__":
    main()
