#!/usr/bin/env python3
"""
Compare Aggrescan4D results between two pH conditions.

Built for the comparison between the yeast cytosol (pH 7.0, intracellular
expression in Pichia pastoris) and the E. coli cytosol (pH 7.6).

USAGE
    python3 compare_ph_conditions.py <dir pH 7.0> <dir pH 7.6> [label1] [label2]

    python3 compare_ph_conditions.py ~/Downloads/Tables ~/Downloads/Tables_pH76 \
        "pH 7.0 (yeast)" "pH 7.6 (E. coli)"

Each directory must contain the eight CSV exports for that condition, named as
the server produces them.

OUTPUT
    a4d_ph_comparison.csv     side-by-side totals and deltas
    printed comparison table

WHAT TO LOOK FOR
    The absolute scores will shift with pH. The question is whether the
    ranking and the delta-versus-wild-type values are preserved. If they are,
    the benefit of the mutations does not depend on the expression host.
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


def read_csv(path):
    raw = open(path, encoding="utf-8-sig", errors="ignore").read()
    best, best_score = None, -1
    for d in [",", ";", "\t", "|"]:
        rows = [r for r in csv.reader(raw.splitlines(), delimiter=d)]
        w = [len(r) for r in rows if r]
        if not w:
            continue
        mode = max(set(w), key=w.count)
        if mode < 2:
            continue
        sc = (w.count(mode) / len(w)) * mode
        if sc > best_score:
            best, best_score = rows, sc
    return best or []


def totals(directory):
    """→ {job: (total, average, n_residues)}"""
    out = {}
    if not os.path.isdir(directory):
        sys.exit(f"Directory not found: {directory}")
    for fn in sorted(os.listdir(directory)):
        if not fn.lower().endswith(".csv"):
            continue
        job = job_from_filename(fn)
        if job not in JOBS:
            continue
        rows = read_csv(os.path.join(directory, fn))
        if not rows:
            continue
        header = rows[0]
        # locate the score column by header name, else take the last float column
        col = None
        for i, h in enumerate(header):
            if h.strip().lower() in ("score", "a3d_score", "a4d_score", "value"):
                col = i
                break
        scores = []
        for r in rows[1:]:
            if col is not None and len(r) > col:
                try:
                    scores.append(float(r[col].strip()))
                    continue
                except ValueError:
                    pass
            for cell in reversed(r):                 # fallback: last numeric cell
                try:
                    scores.append(float(cell.strip()))
                    break
                except ValueError:
                    continue
        if scores:
            out[job] = (sum(scores), sum(scores) / len(scores), len(scores))
    return out


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    d1, d2 = os.path.expanduser(sys.argv[1]), os.path.expanduser(sys.argv[2])
    l1 = sys.argv[3] if len(sys.argv) > 3 else "condition 1"
    l2 = sys.argv[4] if len(sys.argv) > 4 else "condition 2"

    t1, t2 = totals(d1), totals(d2)
    missing = [j for j in JOBS if j not in t1 or j not in t2]
    if missing:
        print("  MISSING in one or both conditions: " + ", ".join(missing) + "\n")

    rows = []
    for job, (species, variant, nmut) in JOBS.items():
        if job not in t1 or job not in t2:
            continue
        a, b = t1[job], t2[job]
        wt1 = t1.get(WT_OF[species])
        wt2 = t2.get(WT_OF[species])
        rows.append(dict(
            job=job, species=species, variant=variant, n_mutations=nmut,
            n_residues=a[2],
            total_1=round(a[0], 4), total_2=round(b[0], 4),
            average_1=round(a[1], 4), average_2=round(b[1], 4),
            shift_total=round(b[0] - a[0], 4),
            delta_wt_1=round(a[0] - wt1[0], 4) if wt1 else None,
            delta_wt_2=round(b[0] - wt2[0], 4) if wt2 else None,
        ))
        rows[-1]["delta_agreement"] = (
            round(rows[-1]["delta_wt_2"] - rows[-1]["delta_wt_1"], 4)
            if wt1 and wt2 else None)

    w = 78
    print("=" * w)
    print(f"{'':26s} {l1:>22s} {l2:>22s}")
    print("=" * w)
    for r in rows:
        print(f"{r['species'][:9]:9s} {r['variant']:15s} "
              f"{r['total_1']:12.4f} {'':9s} {r['total_2']:12.4f}")
    print()
    print("=" * w)
    print("Delta versus wild type — the value that matters")
    print("=" * w)
    print(f"{'':26s} {l1:>18s} {l2:>18s} {'agreement':>12s}")
    for r in rows:
        if r["n_mutations"] == 0:
            continue
        print(f"{r['species'][:9]:9s} {r['variant']:15s} "
              f"{r['delta_wt_1']:+18.4f} {r['delta_wt_2']:+18.4f} "
              f"{r['delta_agreement']:+12.4f}")

    print()
    print("=" * w)
    print("Additivity, both conditions")
    print("=" * w)
    for species in ("P_sativum", "O_spinosa"):
        grp = [r for r in rows if r["species"] == species]
        singles = [r for r in grp if r["n_mutations"] == 1]
        dbl = next((r for r in grp if r["n_mutations"] == 2), None)
        if len(singles) == 2 and dbl:
            for key, lab in (("delta_wt_1", l1), ("delta_wt_2", l2)):
                pred = sum(s[key] for s in singles)
                obs = dbl[key]
                print(f"  {species:10s} {lab:20s} predicted={pred:+8.4f}  "
                      f"observed={obs:+8.4f}  epistasis={obs - pred:+8.4f}")

    if rows:
        with open("a4d_ph_comparison.csv", "w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            wr.writeheader()
            wr.writerows(rows)
        print(f"\n  wrote a4d_ph_comparison.csv  ({len(rows)} rows)")

    print("\n  Interpretation: absolute scores shift with pH, which is expected.")
    print("  If the delta-versus-wild-type values agree between conditions, the")
    print("  benefit of the mutations is host-independent. Large disagreement")
    print("  would mean the effect depends on the ionisation state of nearby")
    print("  titratable residues.")


if __name__ == "__main__":
    main()
