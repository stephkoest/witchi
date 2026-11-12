#!/usr/bin/env python3
"""Summarize bench/scaling.tsv (from scripts/bench_scaling.py) as Markdown tables.

Usage:
    python3 -I scripts/bench_report.py [bench/scaling.tsv] > bench/scaling_tables.md
"""

import csv
import statistics
import sys
from collections import defaultdict


def load(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    groups = defaultdict(list)
    for r in rows:
        groups[(r["dataset"], r["series"], int(r["workers"]))].append(r)
    return rows, groups


def med(rs, key):
    return statistics.median(float(r[key]) for r in rs)


def fmt_t(x):
    return f"{x:.2f}" if x < 10 else f"{x:.1f}" if x < 100 else f"{x:.0f}"


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "bench/scaling.tsv"
    rows, groups = load(path)
    datasets = []
    for r in rows:
        if r["dataset"] not in datasets:
            datasets.append(r["dataset"])
    workers = sorted({int(r["workers"]) for r in rows})
    meta = {r["dataset"]: (int(r["taxa"]), int(r["cols"])) for r in rows}

    def get(ds, series, w):
        if w == 1:
            series = "both"
        return groups.get((ds, series, w))

    def label(ds):
        t, c = meta[ds]
        return f"{t:,} × {c:,}"

    out = []
    p = out.append

    p("### Wall time (s), `--num_workers_permute w --num_workers_chisq w`\n")
    p("| taxa × cols | " + " | ".join(f"w={w}" for w in workers) + " |")
    p("|---|" + "---:|" * len(workers))
    for ds in datasets:
        gs = [get(ds, "both", w) for w in workers]
        cells = [fmt_t(med(g, "wall_s")) if g else "–" for g in gs]
        p(f"| {label(ds)} | " + " | ".join(cells) + " |")

    p("\n### Speed-up over w=1 (parallel efficiency), same runs\n")
    p("| taxa × cols | " + " | ".join(f"w={w}" for w in workers[1:]) + " |")
    p("|---|" + "---:|" * (len(workers) - 1))
    for ds in datasets:
        base = get(ds, "both", 1)
        cells = []
        for w in workers[1:]:
            g = get(ds, "both", w)
            if not (base and g):
                cells.append("–")
                continue
            s = med(base, "wall_s") / med(g, "wall_s")
            cells.append(f"{s:.1f}× ({100 * s / w:.0f}%)")
        p(f"| {label(ds)} | " + " | ".join(cells) + " |")

    p("\n### Which flag is raised: wall time (s) for permute-only / chisq-only / both\n")
    p("| taxa × cols | " + " | ".join(f"w={w}" for w in workers[1:]) + " |")
    p("|---|" + "---:|" * (len(workers) - 1))
    for ds in datasets:
        cells = []
        for w in workers[1:]:
            vals = [get(ds, s, w) for s in ("permute", "chisq", "both")]
            cells.append(" / ".join(fmt_t(med(g, "wall_s")) if g else "–" for g in vals))
        p(f"| {label(ds)} | " + " | ".join(cells) + " |")

    p("\n### Stage times (s), both flags = w\n")
    for ds in datasets:
        p(f"\n**{label(ds)}**\n")
        p("| w | read | null matrix | delta null | prune iterations (n) | write | other | wall |")
        p("|---:|---:|---:|---:|---:|---:|---:|---:|")
        for w in workers:
            g = get(ds, "both", w)
            if not g:
                continue
            st = {k: med(g, k) for k in ("t_read", "t_null", "t_delta_null", "t_iterations", "t_write", "wall_s")}
            other = st["wall_s"] - sum(v for k, v in st.items() if k != "wall_s")
            p(f"| {w} | {st['t_read']:.2f} | {st['t_null']:.2f} | {st['t_delta_null']:.2f} | "
              f"{st['t_iterations']:.2f} ({g[0]['n_iterations']}) | {st['t_write']:.2f} | "
              f"{other:.2f} | {st['wall_s']:.2f} |")

    p("\n### Peak RSS (MB), both flags = w\n")
    p("| taxa × cols | " + " | ".join(f"w={w}" for w in workers) + " |")
    p("|---|" + "---:|" * len(workers))
    for ds in datasets:
        gs = [get(ds, "both", w) for w in workers]
        cells = [f"{max(float(r['maxrss_mb']) for r in g):.0f}" if g else "–" for g in gs]
        p(f"| {label(ds)} | " + " | ".join(cells) + " |")

    p("\n### Output consistency\n")
    p("| taxa × cols | runs | distinct outputs | by `--num_workers_chisq` | pruned | stop reason |")
    p("|---|---:|---:|---|---:|---|")
    for ds in datasets:
        rs = [r for r in rows if r["dataset"] == ds]
        by_chisq = defaultdict(set)
        for r in rs:
            by_chisq[int(r["chisq"])].add(r["outputs_md5"])
        md5s = {r["outputs_md5"] for r in rs}
        # group chisq values that share the same output
        classes = defaultdict(list)
        for c, ms in sorted(by_chisq.items()):
            classes[tuple(sorted(ms))].append(c)
        desc = "; ".join(
            f"chisq {','.join(map(str, cs))}" + (" (varies!)" if len(ms) > 1 else "")
            for ms, cs in classes.items()
        ) if len(md5s) > 1 else "all identical"
        pruned = sorted({r["n_pruned"] for r in rs})
        stops = sorted({r["stop"].split(" (")[0] for r in rs})
        p(f"| {label(ds)} | {len(rs)} | {len(md5s)} | {desc} | {', '.join(pruned)} | {'; '.join(stops)} |")

    print("\n".join(out))


if __name__ == "__main__":
    main()
