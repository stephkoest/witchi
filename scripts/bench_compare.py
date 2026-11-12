#!/usr/bin/env python3
"""Compare two witchi binaries: identical output, wall time and peak memory.

Usage:
    python3 -I scripts/bench_compare.py CASES.tsv OUT_DIR NAME=BINARY NAME=BINARY [--reps N]

CASES.tsv has tab-separated columns: case name, input path (relative to the
repository root or absolute), repetitions (0 = use --reps), and the `witchi
prune` arguments (without --file). Blank lines and lines starting with # are
skipped.

Every run happens in its own directory under OUT_DIR/runs with the input
symlinked in, under GNU time and with WITCHI_TIMING=1. For each run the script
records wall time, peak RSS, the stage timings and the MD5 of the four output
files and of stdout (minus the "Execution time" line); the output files are
then deleted. The binaries alternate within each repetition, starting with a
different one each time.

Writes OUT_DIR/runs.tsv (one row per run) and OUT_DIR/summary.tsv (per case
and binary: minimum and median wall time, maximum peak RSS, and whether every
run of the case produced the same checksums).
"""

import hashlib
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STAGES = ("read alignment", "null matrix", "null state", "delta null", "write outputs")


def md5(path):
    h = hashlib.md5()
    with open(str(path), "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_cases(path):
    cases = []
    for line in open(path):
        line = line.rstrip("\n")
        if not line.strip() or line.startswith("#"):
            continue
        name, inp, reps, args = line.split("\t")
        p = Path(inp)
        if not p.is_absolute():
            p = ROOT / p
        cases.append((name, p, int(reps), args.split()))
    return cases


def run_one(binary, inp, args, run_dir):
    if run_dir.exists():
        shutil.rmtree(str(run_dir))
    run_dir.mkdir(parents=True)
    (run_dir / inp.name).symlink_to(inp)
    env = dict(os.environ, WITCHI_TIMING="1")
    cmd = ["/usr/bin/time", "-f", "TIME %e %M", binary, "prune", "--file", inp.name] + args
    with open(str(run_dir / "stdout.txt"), "wb") as out, open(str(run_dir / "stderr.txt"), "wb") as err:
        rc = subprocess.call(cmd, cwd=str(run_dir), stdout=out, stderr=err, env=env)
    stderr = (run_dir / "stderr.txt").read_text()
    if rc != 0:
        sys.exit("run failed in {}:\n{}".format(run_dir, stderr[-2000:]))
    m = re.search(r"^TIME (\S+) (\d+)$", stderr, re.M)
    wall, rss_kb = float(m.group(1)), int(m.group(2))
    stages = {}
    iters = 0.0
    for label, secs in re.findall(r"^\[timing\] ([^:]+): ([0-9.]+)s$", stderr, re.M):
        if label in STAGES:
            stages[label] = stages.get(label, 0.0) + float(secs)
        elif label.startswith("iteration"):
            iters += float(secs)
    stages["iterations"] = iters
    stdout = (run_dir / "stdout.txt").read_bytes()
    stdout = b"".join(l for l in stdout.splitlines(True) if not l.startswith(b"Execution time"))
    sums = [hashlib.md5(stdout).hexdigest()]
    outputs = sorted(p for p in run_dir.iterdir() if "_pruned" in p.name)
    if len(outputs) != 4:
        sys.exit("expected 4 output files in {}".format(run_dir))
    for p in outputs:
        sums.append(md5(p))
    shutil.rmtree(str(run_dir))
    return wall, rss_kb, stages, ",".join(sums)


def median(xs):
    xs = sorted(xs)
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


def main():
    args = sys.argv[1:]
    reps = 1
    if "--reps" in args:
        i = args.index("--reps")
        reps = int(args[i + 1])
        del args[i : i + 2]
    if len(args) < 4:
        sys.exit(__doc__)
    cases_path, out_dir = args[:2]
    binaries = [a.split("=", 1) for a in args[2:]]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cases = read_cases(cases_path)

    stage_cols = list(STAGES) + ["iterations"]
    runs_path = out / "runs.tsv"
    with open(str(runs_path), "w") as fh:
        fh.write("\t".join(["case", "binary", "rep", "wall_s", "rss_mb"] + stage_cols + ["md5s"]) + "\n")
    results = {}
    for name, inp, case_reps, cargs in cases:
        n = case_reps or reps
        for rep in range(n):
            order = binaries[rep % len(binaries):] + binaries[: rep % len(binaries)]
            for bname, binary in order:
                t0 = time.time()
                wall, rss, stages, sums = run_one(binary, inp, cargs, out / "runs" / name / bname)
                results.setdefault(name, {}).setdefault(bname, []).append((wall, rss, sums))
                row = [name, bname, str(rep), "%.2f" % wall, "%.0f" % (rss / 1024)]
                row += ["%.2f" % stages.get(s, float("nan")) for s in stage_cols] + [sums]
                with open(str(runs_path), "a") as fh:
                    fh.write("\t".join(row) + "\n")
                print("{} {} rep {}: {:.2f} s, {:.0f} MB ({:.0f} s elapsed)".format(
                    name, bname, rep, wall, rss / 1024, time.time() - t0), flush=True)

    with open(str(out / "summary.tsv"), "w") as fh:
        head = ["case"]
        for bname, _ in binaries:
            head += [bname + "_wall_min_s", bname + "_wall_median_s", bname + "_rss_max_mb"]
        head += ["identical"]
        fh.write("\t".join(head) + "\n")
        for name, _, _, _ in cases:
            row = [name]
            sums = set()
            for bname, _ in binaries:
                rs = results[name][bname]
                walls = [r[0] for r in rs]
                row += ["%.2f" % min(walls), "%.2f" % median(walls), "%.0f" % (max(r[1] for r in rs) / 1024)]
                sums.update(r[2] for r in rs)
            row.append("yes" if len(sums) == 1 else "NO")
            fh.write("\t".join(row) + "\n")
    print((out / "summary.tsv").read_text())


if __name__ == "__main__":
    main()
