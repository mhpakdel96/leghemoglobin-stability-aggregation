# WARNING - THIS ANALYSIS DID NOT PRODUCE A USABLE RESULT
#
# Hydrophobic patch clustering from MD trajectories. Patch detection works,
# but the direction of the trend inverts between exposure thresholds, failing
# the pre-set stability criterion. Retained for documentation only.
# Aggrescan4D was used instead. See docs/RESULTS.md section 8.

import os, numpy as np, pandas as pd
import MDAnalysis as mda
from scipy import stats

BASE = "/cfs/earth/scratch/xpkk/Dr_Salmanian"

MAXSASA = {"ALA":1.29,"ARG":2.74,"ASN":1.95,"ASP":1.93,"CYS":1.67,
           "GLN":2.25,"GLU":2.23,"GLY":1.04,"HIS":2.24,"ILE":1.97,
           "LEU":2.01,"LYS":2.36,"MET":2.24,"PHE":2.40,"PRO":1.59,
           "SER":1.55,"THR":1.72,"TRP":2.85,"TYR":2.63,"VAL":1.74}
HPHOB = {"ALA","VAL","LEU","ILE","MET","PHE","TRP","TYR","CYS","PRO"}

REL_CUTS  = [15.0, 25.0]        # آستانه نمایانی
DIST_CUTS = [4.5, 5.5, 6.5]     # آستانه تماس اتم‌به‌اتم زنجیره جانبی
NFRAMES   = 20

SYS = [
 ("Without_O2","LegHEM22",            "P_sativum","WT",         0,(1,143)),
 ("Without_O2","LegHEMV53T",          "P_sativum","V53T",       1,(1,143)),
 ("Without_O2","LegHEMV83T",          "P_sativum","V83T",       1,(1,143)),
 ("Without_O2","LegHEMV53V83T",       "P_sativum","V53T+V83T",  2,(1,143)),
 ("Without_O2","SpinoHEMNatural_v2",  "O_spinosa","WT",         0,(9,165)),
 ("Without_O2","SpinoHEML43W_v2",     "O_spinosa","L43W",       1,(9,165)),
 ("Without_O2","SpinoHEMF125W_v2",    "O_spinosa","F125W",      1,(9,165)),
 ("Without_O2","SpinoHEMF125L43W_v2", "O_spinosa","F125W+L43W", 2,(9,165)),
]

def read_perres(path):
    v = []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return v
    for line in open(path):
        s = line.strip()
        if not s or s[0] in "#@": continue
        try: v.append(float(s.split()[1]))
        except (IndexError, ValueError): pass
    return v

def components(pairs, n):
    par = list(range(n))
    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb: par[ra] = rb
    cl = {}
    for a in range(n): cl.setdefault(find(a), []).append(a)
    return list(cl.values())

rows = []
for group, folder, sp, mut, nmut, (r1, r2) in SYS:
    A = f"{BASE}/{group}/{folder}/gromacs/analysis_struct"
    C = f"{BASE}/{group}/{folder}/gromacs/analysis_core"
    tpr, xtc = f"{A}/prot.tpr", f"{A}/prot_rep1.xtc"
    if not (os.path.exists(tpr) and os.path.exists(xtc)):
        print(f"  SKIP {folder}"); continue

    u = mda.Universe(tpr, xtc)
    core  = u.residues[r1-1:r2]
    names = [r.resname for r in core]
    n     = len(core)

    per = [v for v in (read_perres(f"{C}/rep{rep}_sasa_perres.xvg")
                       for rep in (1,2,3)) if len(v) == n]
    if not per:
        print(f"  داده SASA نیست: {folder}"); continue
    sasa = np.mean(per, axis=0)
    rel  = np.array([100*s/MAXSASA.get(nm, 2.0) for s, nm in zip(sasa, names)])

    sc = [ (r.atoms.select_atoms("not backbone and not name H*") or r.atoms)
           for r in core ]
    frames = np.linspace(0, len(u.trajectory)-1, NFRAMES).astype(int)

    for rc in REL_CUTS:
        idx = [i for i in range(n) if names[i] in HPHOB and rel[i] > rc]
        if len(idx) < 2: continue
        area = {i: float(sasa[i]) for i in idx}
        for dc in DIST_CUTS:
            npatch, big_res, big_area = [], [], []
            for fi in frames:
                u.trajectory[fi]
                pos = [sc[i].positions for i in idx]
                pairs = []
                for a in range(len(idx)):
                    for b in range(a+1, len(idx)):
                        d = np.linalg.norm(pos[a][:,None,:]-pos[b][None,:,:],
                                           axis=-1).min()
                        if d < dc: pairs.append((a, b))
                cl = components(pairs, len(idx))
                cl = [c for c in cl if len(c) > 1] or cl
                ar = [sum(area[idx[j]] for j in c) for c in cl]
                k  = int(np.argmax(ar))
                npatch.append(len(cl)); big_res.append(len(cl[k]))
                big_area.append(ar[k])
            rows.append(dict(species=sp, mutation=mut, n_mut=nmut,
                             rel_cut=rc, dist_cut=dc, n_sel=len(idx),
                             n_patches=np.mean(npatch),
                             largest_res=np.mean(big_res),
                             largest_area=np.mean(big_area),
                             largest_area_sd=np.std(big_area)))
            print(f"  {folder:22s} rel>{rc:4.0f}% d<{dc}Å : "
                  f"انتخاب={len(idx):3d}  لکه‌ها={np.mean(npatch):5.1f}  "
                  f"بزرگ‌ترین={np.mean(big_res):5.1f} باقی‌مانده / "
                  f"{np.mean(big_area):6.2f}±{np.std(big_area):4.2f} nm²")

df = pd.DataFrame(rows)
df.to_csv(f"{BASE}/results_patches_v2.csv", index=False)

print("\n" + "="*96)
print("آزمون روند در هر ترکیب پارامتر — یافته وقتی معتبر است که در همه پایدار بماند")
print("="*96)
for sp in df.species.unique():
    print(f"\n--- {sp} ---")
    for rc in REL_CUTS:
        for dc in DIST_CUTS:
            s = df[(df.species==sp)&(df.rel_cut==rc)&(df.dist_cut==dc)]
            if len(s) < 3: continue
            out = []
            for m in ["n_patches","largest_res","largest_area"]:
                r = stats.linregress(s.n_mut.values.astype(float), s[m].values)
                out.append(f"{m}: شیب={r.slope:+7.3f} p={r.pvalue:.3f}"
                           f"{'*' if r.pvalue<0.05 else ' '}")
            print(f"  rel>{rc:4.0f}% d<{dc}Å   " + "   ".join(out))

print("\nذخیره شد: results_patches_v2.csv")
