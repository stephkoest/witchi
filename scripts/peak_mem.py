#!/usr/bin/env python3
"""Run a command and report the peak memory of its whole process tree.

Usage:
    python3 -I scripts/peak_mem.py [--interval SECONDS] -- COMMAND [ARGS...]

GNU time reports the peak RSS of the largest single process, which understates
tools that fork worker processes (Python WitChi's joblib workers). This samples
every process descendant of COMMAND every --interval seconds (default 0.2) and
reports, on stderr:
  - peak_tree_pss_mb: the peak of the summed proportional set size (PSS, from
    /proc/PID/smaps_rollup), which splits shared pages between the processes
    that map them, so the sum is the memory the tree actually uses;
  - peak_tree_rss_mb: the peak of the summed RSS (counts shared pages more than
    once, so an upper bound);
  - peak_process_rss_mb: the largest single-process RSS seen;
  - wall_s and the command's exit status.
Linux only. Short spikes between samples can be missed.
"""

import os
import subprocess
import sys
import time


def children_map():
    kids = {}
    for name in os.listdir("/proc"):
        if not name.isdigit():
            continue
        try:
            with open("/proc/%s/stat" % name) as fh:
                stat = fh.read()
        except OSError:
            continue
        ppid = int(stat[stat.rindex(")") + 2 :].split()[1])
        kids.setdefault(ppid, []).append(int(name))
    return kids


def tree(root):
    kids = children_map()
    out, todo = [], [root]
    while todo:
        pid = todo.pop()
        out.append(pid)
        todo.extend(kids.get(pid, []))
    return out


def mem_kb(pid):
    pss = rss = 0
    try:
        with open("/proc/%d/smaps_rollup" % pid) as fh:
            for line in fh:
                if line.startswith("Pss:"):
                    pss = int(line.split()[1])
                elif line.startswith("Rss:"):
                    rss = int(line.split()[1])
    except OSError:
        pass
    return pss, rss


def main():
    args = sys.argv[1:]
    interval = 0.2
    if args[:1] == ["--interval"]:
        interval = float(args[1])
        args = args[2:]
    if args[:1] == ["--"]:
        args = args[1:]
    if not args:
        sys.exit(__doc__)
    t0 = time.time()
    proc = subprocess.Popen(args)
    peak_pss = peak_rss = peak_one = 0
    while proc.poll() is None:
        pss_sum = rss_sum = 0
        for pid in tree(proc.pid):
            pss, rss = mem_kb(pid)
            pss_sum += pss
            rss_sum += rss
            peak_one = max(peak_one, rss)
        peak_pss = max(peak_pss, pss_sum)
        peak_rss = max(peak_rss, rss_sum)
        time.sleep(interval)
    wall = time.time() - t0
    sys.stderr.write(
        "peak_tree_pss_mb\t%.0f\npeak_tree_rss_mb\t%.0f\npeak_process_rss_mb\t%.0f\nwall_s\t%.1f\nexit\t%d\n"
        % (peak_pss / 1024, peak_rss / 1024, peak_one / 1024, wall, proc.returncode)
    )
    sys.exit(proc.returncode)


if __name__ == "__main__":
    main()
