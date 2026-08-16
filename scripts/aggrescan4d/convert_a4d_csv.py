#!/usr/bin/env python3
"""
Convert Aggrescan4D CSV exports into the repository's archived result files.

The A4D server exports per-residue tables as CSV with filenames of the form

    a3d_pea_WT_A3D.csv
    a3d_pea_WT [mutate_ VT53A]_A3D.csv
    a3d_spi_WT [mutate_ LW43A, FW125A]_A3D.csv

This script maps those filenames to job names, auto-detects the column layout,
recomputes the summary statistics from the per-residue data, and verifies them
against the values reported by the web server.

USAGE
    python3 convert_a4d_csv.py <directory containing the CSV files>

    python3 convert_a4d_csv.py ~/Downloads/Tables

OUTPUTS (written to the current directory)
    a4d_per_residue.csv    long-format table, all jobs
    a4d_summary.csv        one row per job, recomputed
    a4d_hotspots.csv       residues scoring above zero
    raw/<job>.txt          normalised tab-separated copies

The recomputed totals are compared against the values recorded from the server.
A mismatch means a table was captured incompletely.
"""

import csv
import os
import re
import sys
from collections import OrderedDict

# job name -> (species, variant label, number of mutations)
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

# totals recorded from the server, used as an integrity check
EXPECTED_TOTAL = {
    "pea_WT":     -139.1718, "pea_V53T":   -141.5708,
    "pea_V83T":   -145.5529, "pea_double": -147.9519,
    "spi_WT":     -173.0987, "spi_L43W":   -172.7161,
    "spi_F125W":  -175.1932, "spi_double": -173.1666,
}

AA3 = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
       "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
AA1 = set("ACDEFGHIKLMNPQRSTVWY")


def job_from_filename(fn):
    """Map an exported filename to a job key."""
    low = fn.lower()
    species = "pea" if "pea" in low else ("spi" if "spi" in low else None)
    if species is None:
        return None
    muts = re.findall(r"[A-Z]{2}\d+[A-Z]", fn)      # e.g. VT53A, LW43A
    if not muts:
        return f"{species}_WT"
    if len(muts) > 1:
        return f"{species}_double"
    m = muts[0]
    wt_aa, new_aa, pos = m[0], m[1], re.search(r"\d+", m).group()
    return f"{species}_{wt_aa}{pos}{new_aa}"


def sniff(rows):
    """Identify which columns hold residue index, name, chain and score."""
    body = [r for r in rows if r and any(c.strip() for c in r)]
    if not body:
        return None
    # drop a leading header row: one with no numeric cell at all
    def numeric(cell):
        return bool(re.fullmatch(r"-?\d*\.?\d+(?:[eE][-+]?\d+)?", cell.strip()))
    while body and not any(numeric(c) for c in body[0]):
        body = body[1:]
    if not body:
        return None
    ncol = max(len(r) for r in body)

    prof = []
    for c in range(ncol):
        col = [r[c].strip() for r in body if len(r) > c and r[c].strip()]
        n = max(len(col), 1)
        prof.append(dict(
            c=c, n=len(col),
            ints=sum(1 for v in col if re.fullmatch(r"-?\d+", v)) / n,
            floats=sum(1 for v in col
                       if re.fullmatch(r"-?\d*\.\d+(?:[eE][-+]?\d+)?", v)) / n,
            res3=sum(1 for v in col if v.upper() in AA3) / n,
            res1=sum(1 for v in col if len(v) == 1 and v.upper() in AA1) / n,
            single=sum(1 for v in col if len(v) == 1 and v.isalpha()) / n,
            uniq=len(set(v.upper() for v in col)),
        ))

    def best(key, exclude=()):
        cands = [p for p in prof if p["c"] not in exclude and p[key] > 0.8]
        return max(cands, key=lambda p: (p[key], p["n"]))["c"] if cands else None

    score = best("floats")
    idx   = best("ints", exclude={score})
    # a column of identical single letters is the chain, never the residue name
    chain = next((p["c"] for p in prof
                  if p["c"] not in {score, idx}
                  and p["single"] > 0.9 and p["uniq"] <= 4), None)
    # prefer three-letter codes; fall back to one-letter
    name = best("res3", exclude={score, idx, chain})
    if name is None:
        name = best("res1", exclude={score, idx, chain})
    return dict(idx=idx, name=name, chain=chain, score=score, ncol=ncol)

def read_csv(path):
    """Read a CSV, choosing the delimiter that yields consistent column counts."""
    raw = open(path, encoding="utf-8-sig", errors="ignore").read()
    best, best_score = None, -1
    for delim in [",", ";", "\t", "|"]:
        rows = [r for r in csv.reader(raw.splitlines(), delimiter=delim)]
        widths = [len(r) for r in rows if r]
        if not widths:
            continue
        mode = max(set(widths), key=widths.count)
        if mode < 2:
            continue
        consistency = widths.count(mode) / len(widths)
        sc = consistency * mode
        if sc > best_score:
            best, best_score = rows, sc
    return best or []


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    src = os.path.expanduser(sys.argv[1])
    files = [f for f in os.listdir(src) if f.lower().endswith(".csv")]
    if not files:
        sys.exit(f"No CSV files found in {src}")

    os.makedirs("raw", exist_ok=True)
    per_residue, summaries, hotspots = [], OrderedDict(), []

    print(f"Found {len(files)} CSV files in {src}\n")

    for fn in sorted(files):
        job = job_from_filename(fn)
        if job not in JOBS:
            print(f"  SKIP  {fn}\n        -> could not map to a job (got {job!r})")
            continue

        rows = read_csv(os.path.join(src, fn))
        col = sniff(rows)
        if col is None or col["score"] is None or col["idx"] is None:
            print(f"  FAIL  {fn}\n        -> could not identify columns: {col}")
            continue

        data = []
        for r in rows:
            if len(r) <= max(x for x in (col["idx"], col["score"]) if x is not None):
                continue
            try:
                i = int(r[col["idx"]].strip())
                v = float(r[col["score"]].strip())
            except ValueError:
                continue                                  # header row
            data.append(dict(
                index=i,
                resname=(r[col["name"]].strip() if col["name"] is not None
                         and len(r) > col["name"] else ""),
                chain=(r[col["chain"]].strip() if col["chain"] is not None
                       and len(r) > col["chain"] else "A"),
                score=v,
            ))

        if not data:
            print(f"  FAIL  {fn} -> no data rows parsed")
            continue

        species, variant, nmut = JOBS[job]
        total = sum(d["score"] for d in data)
        hi = max(data, key=lambda d: d["score"])
        lo = min(data, key=lambda d: d["score"])

        exp = EXPECTED_TOTAL.get(job)
        ok = "" if exp is None else (
            "  MATCHES server" if abs(total - exp) < 0.01
            else f"  MISMATCH: server reported {exp:.4f}")

        print(f"  {job:12s} <- {fn}")
        print(f"    columns idx={col['idx']} name={col['name']} "
              f"chain={col['chain']} score={col['score']}")
        print(f"    {len(data):4d} residues   total = {total:10.4f}{ok}")

        summaries[job] = dict(
            job=job, species=species, variant=variant, n_mutations=nmut,
            n_residues=len(data),
            total_score=round(total, 4),
            average_score=round(total / len(data), 4),
            max_score=round(hi["score"], 4),
            max_residue=f"{hi['resname']}{hi['index']}",
            min_score=round(lo["score"], 4),
            min_residue=f"{lo['resname']}{lo['index']}",
            n_positive=sum(1 for d in data if d["score"] > 0),
        )

        for d in data:
            per_residue.append(dict(job=job, species=species, variant=variant, **d))
            if d["score"] > 0:
                hotspots.append(dict(job=job, species=species, variant=variant,
                                     index=d["index"], resname=d["resname"],
                                     score=round(d["score"], 4)))

        with open(os.path.join("raw", job + ".txt"), "w", encoding="utf-8") as fh:
            fh.write("residue index\tresidue name\tchain\tAggrescan4D score\n")
            for d in data:
                fh.write(f"{d['index']}\t{d['resname']}\t{d['chain']}\t{d['score']}\n")

    missing = [j for j in JOBS if j not in summaries]
    if missing:
        print("\n  MISSING JOBS: " + ", ".join(missing))

    for job, s in summaries.items():
        wt = summaries.get(WT_OF[s["species"]])
        if wt and job != WT_OF[s["species"]]:
            s["delta_total_vs_wt"] = round(s["total_score"] - wt["total_score"], 4)
            s["delta_average_vs_wt"] = round(s["average_score"] - wt["average_score"], 4)
        else:
            s["delta_total_vs_wt"] = 0.0
            s["delta_average_vs_wt"] = 0.0

    def write(path, rows, fields):
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        print(f"  wrote {path}  ({len(rows)} rows)")

    print()
    if summaries:
        write("a4d_summary.csv", list(summaries.values()),
              list(next(iter(summaries.values())).keys()))
    if per_residue:
        write("a4d_per_residue.csv", per_residue,
              ["job", "species", "variant", "index", "resname", "chain", "score"])
    if hotspots:
        write("a4d_hotspots.csv", hotspots,
              ["job", "species", "variant", "index", "resname", "score"])

    print("\n  Additivity check:")
    for species in {s["species"] for s in summaries.values()}:
        grp = [s for s in summaries.values() if s["species"] == species]
        singles = [s for s in grp if s["n_mutations"] == 1]
        double = next((s for s in grp if s["n_mutations"] == 2), None)
        if len(singles) == 2 and double:
            pred = sum(s["delta_total_vs_wt"] for s in singles)
            obs = double["delta_total_vs_wt"]
            print(f"    {species:10s} predicted={pred:+8.4f}  observed={obs:+8.4f}  "
                  f"epistasis={obs - pred:+8.4f}")

    print("\n  Score convention: positive = aggregation-prone, "
          "negative = solubility-promoting.")


if __name__ == "__main__":
    main()
