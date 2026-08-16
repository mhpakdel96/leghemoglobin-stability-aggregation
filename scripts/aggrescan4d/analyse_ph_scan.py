#!/usr/bin/env python3
"""
Analyse Aggrescan4D pH-scan tables across all variants.

The A4D server computes each job over a pH range and reports the average and
maximum A4D score at each point. This script reads those tables for all eight
jobs and asks the question that matters for host selection:

    does the benefit of each mutation hold across the pH range,
    or does it depend on the ionisation state of the environment?

USAGE
    python3 analyse_ph_scan.py <directory containing the pH scan files>

    python3 analyse_ph_scan.py ~/Downloads/Tables_pH_scan

Input files may have any extension (or none). Each must contain a header line
followed by rows of: pH, average score, max score.

IMPORTANT
    Values from the pH-scan module do not necessarily match those from the
    main per-residue result page, even at the same pH. Do not mix the two
    datasets. Comparisons within this scan are internally consistent.

REFERENCE pH VALUES
    Yeast cytosol (Pichia pastoris)   approximately 7.0
    E. coli cytosol                   approximately 7.6
"""

import csv
import os
import re
import sys
from collections import OrderedDict

JOBS = OrderedDict([
    ("pea_WT",     ("P_sativum", "WT",         0)),
    ("pea_V53T",   ("P_sativum", "V53T",       1)),
    ("pea_V83T",   ("P_sativum", "V83T",       1)),
    ("pea_double", ("P_sativum", "V53T+V83T",  2)),
    ("spi_WT",     ("O_spinosa", "WT",         0)),
    ("spi_L43W",   ("O_spinosa", "L43W",       1)),
    ("spi_F125W",  ("O_spinosa", "F125W",      1)),
    ("spi_double", ("O_spinosa", "L43W+F125W", 2)),
])
WT_OF = {"P_sativum": "pea_WT", "O_spinosa": "spi_WT"}

HOST_PH = {"yeast (Pichia)": 7.0, "E. coli": 7.5}   # 7.5 is the nearest scan point


def job_from_filename(fn):
    low = fn.lower()
    sp = "pea" if "pea" in low else ("spi" if "spi" in low else None)
    if sp is None:
        return None
    muts = re.findall(r"[A-Z]{2}\d+[A-Z]", fn)
    if not muts:
        return f"{sp}_WT"
    if len(muts) > 1:
        return f"{sp}_double"
    m = muts[0]
    pos = re.search(r"\d+", m).group()
    return f"{sp}_{m[0]}{pos}{m[1]}"


def read_scan(path):
    """→ {pH: (average, max)}"""
    out = {}
    for line in open(path, encoding="utf-8-sig", errors="ignore"):
        parts = re.split(r"[\t,;]+", line.strip())
        parts = [p for p in parts if p.strip()]
        if len(parts) < 3:
            continue
        try:
            ph, avg, mx = float(parts[0]), float(parts[1]), float(parts[2])
        except ValueError:
            continue                                  # header
        out[round(ph, 1)] = (avg, mx)
    return out


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    src = os.path.expanduser(sys.argv[1])
    if not os.path.isdir(src):
        sys.exit(f"Directory not found: {src}")

    scans = {}
    for fn in sorted(os.listdir(src)):
        p = os.path.join(src, fn)
        if not os.path.isfile(p):
            continue
        job = job_from_filename(fn)
        if job not in JOBS:
            print(f"  SKIP  {fn}  (could not map, got {job!r})")
            continue
        d = read_scan(p)
        if d:
            scans[job] = d
            print(f"  {job:12s} <- {fn}   ({len(d)} pH points, "
                  f"{min(d):.1f} to {max(d):.1f})")

    missing = [j for j in JOBS if j not in scans]
    if missing:
        print("\n  MISSING: " + ", ".join(missing))
    if not scans:
        sys.exit("\nNo usable files found.")

    phs = sorted(set.intersection(*(set(d) for d in scans.values())))
    W = 96

    for metric, mi, note in ((("Max"), 1, "recommended by the server for globular proteins"),
                             (("Average"), 0, "whole-structure mean")):
        print("\n" + "=" * W)
        print(f"{metric} A4D score — {note}")
        print("=" * W)
        head = "pH    " + "".join(f"{JOBS[j][1][:11]:>12s}" for j in scans)
        print(head)
        for ph in phs:
            row = f"{ph:<6.1f}" + "".join(f"{scans[j][ph][mi]:12.4f}" for j in scans)
            mark = ""
            if abs(ph - 7.0) < 0.01:
                mark = "   <- yeast"
            elif abs(ph - 7.5) < 0.01:
                mark = "   <- E. coli"
            print(row + mark)

        print("\n" + "-" * W)
        print(f"Delta versus wild type ({metric}) — negative is better")
        print("-" * W)
        muts = [j for j in scans if JOBS[j][2] > 0]
        print("pH    " + "".join(f"{JOBS[j][1][:11]:>12s}" for j in muts))
        deltas = {j: {} for j in muts}
        for ph in phs:
            cells = ""
            for j in muts:
                wt = scans.get(WT_OF[JOBS[j][0]])
                if not wt:
                    cells += f"{'--':>12s}"
                    continue
                d = scans[j][ph][mi] - wt[ph][mi]
                deltas[j][ph] = d
                cells += f"{d:+12.4f}"
            mark = ("   <- yeast" if abs(ph - 7.0) < 0.01 else
                    "   <- E. coli" if abs(ph - 7.5) < 0.01 else "")
            print(f"{ph:<6.1f}{cells}{mark}")

        print("\n  Host comparison and pH stability of each effect:")
        for j in muts:
            if not deltas[j]:
                continue
            vals = list(deltas[j].values())
            d70 = deltas[j].get(7.0)
            d75 = deltas[j].get(7.5)
            spread = max(vals) - min(vals)
            signs = {"-" if v < 0 else "+" for v in vals}
            verdict = ("sign constant across the range"
                       if len(signs) == 1 else "SIGN CHANGES across the range")
            line = f"    {JOBS[j][1]:14s}"
            if d70 is not None and d75 is not None:
                line += f" yeast={d70:+7.4f}  E.coli={d75:+7.4f}  shift={d75-d70:+7.4f}"
            line += f"   range spread={spread:6.4f}   {verdict}"
            print(line)

    if 7.0 in phs and 7.5 in phs:
        print("\n" + "=" * W)
        print("Summary for host selection")
        print("=" * W)
        for j in scans:
            a70, a75 = scans[j][7.0], scans[j][7.5]
            print(f"  {JOBS[j][0][:9]:9s} {JOBS[j][1]:14s} "
                  f"Max: {a70[1]:7.4f} -> {a75[1]:7.4f} ({a75[1]-a70[1]:+.4f})   "
                  f"Avg: {a70[0]:7.4f} -> {a75[0]:7.4f} ({a75[0]-a70[0]:+.4f})")

    rows = []
    for j, d in scans.items():
        for ph, (avg, mx) in sorted(d.items()):
            rows.append(dict(job=j, species=JOBS[j][0], variant=JOBS[j][1],
                             n_mutations=JOBS[j][2], pH=ph,
                             average_score=avg, max_score=mx))
    with open("a4d_ph_scan.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n  wrote a4d_ph_scan.csv  ({len(rows)} rows)")

    print("\n  Reading guide: if the delta-versus-wild-type values keep the same")
    print("  sign and similar magnitude across the whole pH range, the benefit of")
    print("  the mutations does not depend on the expression host. Compare the")
    print("  size of the yeast-to-E.coli shift against the delta itself — if the")
    print("  shift is much smaller, host pH is not a deciding factor.")


if __name__ == "__main__":
    main()
