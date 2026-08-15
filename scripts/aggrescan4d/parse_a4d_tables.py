#!/usr/bin/env python3
"""
Archive Aggrescan4D web results into version-controllable CSV files.

The Aggrescan4D server returns results as an HTML page. Job URLs expire and
the server may change or go offline, so results must be captured locally.

USAGE
-----
1. For each of the eight jobs, copy the per-residue table from the results
   page (including the header row) and save it as a plain text file in
   ./raw/ using the job name, e.g.:

       raw/pea_WT.txt
       raw/pea_V53T.txt
       raw/pea_V83T.txt
       raw/pea_double.txt
       raw/spi_WT.txt
       raw/spi_L43W.txt
       raw/spi_F125W.txt
       raw/spi_double.txt

   Each file should contain lines of the form:
       residue index<TAB>residue name<TAB>chain<TAB>Aggrescan4D score<TAB>mutation
   Whitespace-separated columns are also accepted.

2. Run:
       python3 parse_a4d_tables.py

OUTPUTS
-------
    a4d_per_residue.csv   long-format table, all jobs, all residues
    a4d_summary.csv       one row per job with totals and deltas
    a4d_hotspots.csv      residues with score > 0 (aggregation-prone)

The summary values are recomputed from the per-residue data rather than
copied from the web page, which serves as an independent check that the
tables were captured completely.
"""

import csv
import os
import re
import sys
from collections import OrderedDict

RAW_DIR = "raw"

# job name -> (species, variant label, number of mutations)
JOBS = OrderedDict([
    ("pea_WT",     ("P_sativum", "WT",           0)),
    ("pea_V53T",   ("P_sativum", "V53T",         1)),
    ("pea_V83T",   ("P_sativum", "V83T",         1)),
    ("pea_double", ("P_sativum", "V53T+V83T",    2)),
    ("spi_WT",     ("O_spinosa", "WT",           0)),
    ("spi_L43W",   ("O_spinosa", "L43W",         1)),
    ("spi_F125W",  ("O_spinosa", "F125W",        1)),
    ("spi_double", ("O_spinosa", "L43W+F125W",   2)),
])

WT_OF = {"P_sativum": "pea_WT", "O_spinosa": "spi_WT"}


def parse_table(path):
    """Return list of dicts: index, resname, chain, score, mutation."""
    rows = []
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            s = line.rstrip("\n").strip()
            if not s:
                continue
            parts = s.split("\t") if "\t" in s else s.split()
            if len(parts) < 4:
                continue
            try:
                idx = int(parts[0])
                score = float(parts[3])
            except ValueError:
                continue                      # header or stray text
            rows.append(dict(
                index=idx,
                resname=parts[1].strip(),
                chain=parts[2].strip(),
                score=score,
                mutation=" ".join(parts[4:]).strip() if len(parts) > 4 else "",
            ))
    return rows


def main():
    if not os.path.isdir(RAW_DIR):
        sys.exit(f"Directory '{RAW_DIR}/' not found. See the docstring.")

    per_residue, summaries, hotspots = [], OrderedDict(), []
    missing = []

    for job, (species, variant, nmut) in JOBS.items():
        path = os.path.join(RAW_DIR, job + ".txt")
        if not os.path.exists(path):
            missing.append(job)
            continue

        rows = parse_table(path)
        if not rows:
            print(f"  WARNING  {job}: no data rows parsed")
            continue

        scores = [r["score"] for r in rows]
        total = sum(scores)
        hi = max(rows, key=lambda r: r["score"])
        lo = min(rows, key=lambda r: r["score"])

        summaries[job] = dict(
            job=job, species=species, variant=variant, n_mutations=nmut,
            n_residues=len(rows),
            total_score=round(total, 4),
            average_score=round(total / len(rows), 4),
            max_score=round(hi["score"], 4),
            max_residue=f"{hi['resname']}{hi['index']}",
            min_score=round(lo["score"], 4),
            min_residue=f"{lo['resname']}{lo['index']}",
            n_positive=sum(1 for s in scores if s > 0),
        )

        for r in rows:
            per_residue.append(dict(job=job, species=species, variant=variant, **r))
            if r["score"] > 0:
                hotspots.append(dict(
                    job=job, species=species, variant=variant,
                    index=r["index"], resname=r["resname"],
                    score=round(r["score"], 4),
                ))

        mutated = [f"{r['resname']}{r['index']}" for r in rows if r["mutation"]]
        sites = ", ".join(mutated) if mutated else "none"
        print(f"  {job:12s} {len(rows):4d} residues   total={total:10.4f}   "
              f"mutated sites: {sites}")

    if missing:
        print("\n  MISSING FILES: " + ", ".join(missing))

    # deltas relative to the matching wild type
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
              ["job", "species", "variant", "index", "resname", "chain",
               "score", "mutation"])
    if hotspots:
        write("a4d_hotspots.csv", hotspots,
              ["job", "species", "variant", "index", "resname", "score"])

    # additivity check for any species with WT, two singles and a double
    print("\n  Additivity check (double vs sum of singles):")
    for species in {s["species"] for s in summaries.values()}:
        grp = [s for s in summaries.values() if s["species"] == species]
        singles = [s for s in grp if s["n_mutations"] == 1]
        double = next((s for s in grp if s["n_mutations"] == 2), None)
        if len(singles) == 2 and double:
            predicted = sum(s["delta_total_vs_wt"] for s in singles)
            observed = double["delta_total_vs_wt"]
            print(f"    {species:10s} predicted={predicted:+8.4f}  "
                  f"observed={observed:+8.4f}  epistasis={observed - predicted:+8.4f}")

    print("\n  Score convention: positive = aggregation-prone, "
          "negative = solubility-promoting.")


if __name__ == "__main__":
    main()
