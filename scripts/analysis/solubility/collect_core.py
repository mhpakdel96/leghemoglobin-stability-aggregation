import os, numpy as np, pandas as pd

BASE = "/cfs/earth/scratch/xpkk/Dr_Salmanian"

SYSTEMS = [
 ("Without_O2","LegHEM22",            "P_sativum","WT"),
 ("Without_O2","LegHEMV53T",          "P_sativum","V53T"),
 ("Without_O2","LegHEMV83T",          "P_sativum","V83T"),
 ("Without_O2","LegHEMV53V83T",       "P_sativum","V53T+V83T"),
 ("Without_O2","SpinoHEMNatural_v2",  "O_spinosa","WT"),
 ("Without_O2","SpinoHEML43W_v2",     "O_spinosa","L43W"),
 ("Without_O2","SpinoHEMF125W_v2",    "O_spinosa","F125W"),
 ("Without_O2","SpinoHEMF125L43W_v2", "O_spinosa","F125W+L43W"),
 ("With_O2","LegHEM22",         "P_sativum","WT"),
 ("With_O2","LegHEMV83T",       "P_sativum","V83T"),
 ("With_O2","SpinoHEMNatural",  "O_spinosa","WT"),
 ("With_O2","SpinoHEMF125W",    "O_spinosa","F125W"),
 ("With_O2","SpinoHEMF125L43W", "O_spinosa","F125W+L43W"),
]

def read_xvg(path, col=1):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    vals = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] in "#@":
                continue
            parts = line.split()
            try:
                vals.append(float(parts[col]))
            except (IndexError, ValueError):
                continue
    return np.array(vals) if vals else None

def helix_fraction(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    hel = tot = 0
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s[0] in "#@=":
                continue
            for c in s:
                if c.isspace():
                    continue
                tot += 1
                if c in "HGIS":
                    hel += 1
    return 100.0 * hel / tot if tot else None

rows = []
for group, folder, species, mut in SYSTEMS:
    O = f"{BASE}/{group}/{folder}/gromacs/analysis_core"
    for rep in (1, 2, 3):
        p = f"{O}/rep{rep}"
        rmsd = read_xvg(f"{p}_rmsd.xvg")
        rmsf = read_xvg(f"{p}_rmsf.xvg")
        rg   = read_xvg(f"{p}_gyrate.xvg")
        sasa = read_xvg(f"{p}_sasa.xvg")
        hb   = read_xvg(f"{p}_hbnum.xvg")
        rows.append(dict(
            group=group, system=folder, species=species, mutation=mut, rep=rep,
            rmsd_bb_A   = rmsd.mean()*10 if rmsd is not None else np.nan,
            rmsf_ca_A   = rmsf.mean()*10 if rmsf is not None else np.nan,
            rg_nm       = rg.mean()      if rg   is not None else np.nan,
            sasa_nm2    = sasa.mean()    if sasa is not None else np.nan,
            hbonds      = hb.mean()      if hb   is not None else np.nan,
            helix_pct   = (helix_fraction(f"{p}_dssp.dat") if os.path.exists(f"{p}_dssp.dat") else np.nan),
        ))

df = pd.DataFrame(rows)
df.to_csv(f"{BASE}/results_core_raw.csv", index=False)

METRICS = ["rmsd_bb_A","rmsf_ca_A","rg_nm","sasa_nm2","hbonds","helix_pct"]
agg = (df.groupby(["group","system","species","mutation"])[METRICS]
         .agg(["mean","std","count"]).reset_index())
agg.columns = ["_".join(c).strip("_") for c in agg.columns]
agg.to_csv(f"{BASE}/results_core_summary.csv", index=False)

pd.set_option("display.width", 200)
for g in ["Without_O2","With_O2"]:
    print("="*110); print(g); print("="*110)
    sub = df[df.group == g]
    for (sysname, species, mut), d in sub.groupby(["system","species","mutation"]):
        n = d[METRICS].notna().all(axis=1).sum()
        print(f"{sysname:24s} {species:10s} {mut:12s} n={len(d.dropna(subset=['rmsd_bb_A']))}  " +
              "  ".join(f"{m}={d[m].mean():7.2f}±{d[m].std():.2f}" for m in METRICS))
    print()
print("ذخیره شد: results_core_raw.csv و results_core_summary.csv")
