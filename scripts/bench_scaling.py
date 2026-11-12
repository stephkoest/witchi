#!/usr/bin/env python3
"""Thread-scaling benchmark of `witchi prune` (Rust port).

Usage:
    python3 -I scripts/bench_scaling.py [--bin target/release/witchi]
        [--data bench/data] [--out bench/scaling.tsv] [--workers 1,2,4,...]
        [--datasets REGEX] [--extra "ARGS"]

For every alignment in --data (see scripts/make_bench_inputs.py) and every worker
count w, runs three configurations:
    permute   --num_workers_permute w --num_workers_chisq 1
    chisq     --num_workers_permute 1 --num_workers_chisq w
    both      --num_workers_permute w --num_workers_chisq w
(w = 1 is run once, as "both"). All other options are WitChi's defaults
(wasserstein, 100 permutations, top_n 1, max_residue_pruned 100, delta-null).

Runs execute one at a time. A configuration whose first run takes under
--repeat-below seconds is repeated twice more (the report uses the median).
Each row records wall/user/sys time and peak RSS (/usr/bin/time), the
WITCHI_TIMING stage totals, the stop line, and an md5 over the four output files.
Rows already in --out are skipped, so an interrupted run can be resumed.
"""

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIELDS = [
    "dataset", "taxa", "cols", "series", "workers", "permute", "chisq", "rep",
    "wall_s", "user_s", "sys_s", "maxrss_mb",
    "t_read", "t_null", "t_delta_null", "t_iterations", "n_iterations", "t_write",
    "n_pruned", "stop", "outputs_md5",
]
STAGES = {
    "read alignment": "t_read",
    "null matrix": "t_null",
    "delta null": "t_delta_null",
    "iteration deltas": "t_iterations",
    "write outputs": "t_write",
}


def configs(workers):
    out = []
    for w in workers:
        if w == 1:
            out.append(("both", 1, 1, 1))
            continue
        out += [("permute", w, w, 1), ("chisq", w, 1, w), ("both", w, w, w)]
    return out


def run_one(binary, data_file, permute, chisq, extra, work):
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    shutil.copy(data_file, work / data_file.name)
    cmd = [
        "/usr/bin/time", "-f", "%e %U %S %M", "-o", "time.txt",
        str(binary), "prune", "--file", data_file.name,
        "--num_workers_permute", str(permute), "--num_workers_chisq", str(chisq),
    ] + extra
    env = dict(os.environ, WITCHI_TIMING="1")
    with open(work / "stdout.txt", "w") as so, open(work / "stderr.txt", "w") as se:
        rc = subprocess.run(cmd, cwd=work, stdout=so, stderr=se, env=env).returncode
    if rc != 0:
        sys.exit(f"run failed in {work} (exit {rc}); see stderr.txt")
    wall, user, sys_, rss = (work / "time.txt").read_text().split()[-4:]
    row = {
        "wall_s": wall, "user_s": user, "sys_s": sys_,
        "maxrss_mb": f"{int(rss) / 1024:.0f}",
    }
    totals = {v: 0.0 for v in STAGES.values()}
    n_iter = 0
    for line in (work / "stderr.txt").read_text().splitlines():
        m = re.match(r"\[timing\] (.+): ([0-9.]+)s$", line)
        if m and m.group(1) in STAGES:
            totals[STAGES[m.group(1)]] += float(m.group(2))
            n_iter += m.group(1) == "iteration deltas"
    row.update({k: f"{v:.3f}" for k, v in totals.items()})
    row["n_iterations"] = str(n_iter)
    stdout = (work / "stdout.txt").read_text()
    m = re.search(r"^Pruned (\d+) → (\d+) columns.*\| stop: ([^|]+)", stdout, re.M)
    row["n_pruned"] = str(int(m.group(1)) - int(m.group(2))) if m else ""
    row["stop"] = m.group(3).strip() if m else ""
    h = hashlib.md5()
    stem = data_file.stem
    for f in sorted(work.iterdir()):
        if f.name.startswith(stem + "_"):
            h.update(f.name.encode())
            h.update(f.read_bytes())
    row["outputs_md5"] = h.hexdigest()
    shutil.rmtree(work)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "target/release/witchi"))
    ap.add_argument("--data", default=str(ROOT / "bench/data"))
    ap.add_argument("--out", default=str(ROOT / "bench/scaling.tsv"))
    ap.add_argument("--workers", default="1,2,4,8,16,32,48,96")
    ap.add_argument("--datasets", default=".")
    ap.add_argument("--extra", default="")
    ap.add_argument("--repeat-below", type=float, default=20.0)
    args = ap.parse_args()

    binary = Path(args.bin).resolve()
    workers = [int(w) for w in args.workers.split(",")]
    extra = args.extra.split()
    out = Path(args.out)
    done = set()
    if out.exists():
        for line in out.read_text().splitlines()[1:]:
            f = line.split("\t")
            done.add(tuple(f[i] for i in (0, 3, 4, 7)))
    else:
        out.write_text("\t".join(FIELDS) + "\n")

    def size_key(p):
        t, c = map(int, re.findall(r"\d+", p.stem)[-2:])
        return (t * c, t)

    files = sorted(
        (p for p in Path(args.data).glob("*.faa") if re.search(args.datasets, p.name)),
        key=size_key,
    )
    for data_file in files:
        taxa, cols = re.findall(r"\d+", data_file.stem)[-2:]
        for series, w, permute, chisq in configs(workers):
            rep = 0
            while True:
                key = (data_file.stem, series, str(w), str(rep))
                if key not in done:
                    t0 = time.time()
                    row = run_one(binary, data_file, permute, chisq, extra,
                                  ROOT / "bench/runs" / "work")
                    row.update(dataset=data_file.stem, taxa=taxa, cols=cols,
                               series=series, workers=w, permute=permute,
                               chisq=chisq, rep=rep)
                    with open(out, "a") as fh:
                        fh.write("\t".join(str(row[k]) for k in FIELDS) + "\n")
                    print(f"{data_file.stem:22s} {series:8s} w={w:<3d} rep={rep} "
                          f"wall={row['wall_s']}s ({time.time() - t0:.0f}s)", flush=True)
                    wall = float(row["wall_s"])
                else:
                    wall = None
                    for line in out.read_text().splitlines()[1:]:
                        f = line.split("\t")
                        if tuple(f[i] for i in (0, 3, 4, 7)) == key:
                            wall = float(f[FIELDS.index("wall_s")])
                rep += 1
                if rep >= 3 or (rep == 1 and wall >= args.repeat_below):
                    break


if __name__ == "__main__":
    main()
