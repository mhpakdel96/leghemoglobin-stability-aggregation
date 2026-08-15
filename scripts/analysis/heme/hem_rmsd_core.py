import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import numpy as np, pandas as pd, os

BASE = "/cfs/earth/scratch/xpkk/Dr_Salmanian"
SKIP_NS = 20.0

SYS = [
 ("Without_O2","LegHEM22","P_sativum","WT","1-143"),
 ("Without_O2","LegHEMV53T","P_sativum","V53T","1-143"),
 ("Without_O2","LegHEMV83T","P_sativum","V83T","1-143"),
 ("Without_O2","LegHEMV53V83T","P_sativum","V53T+V83T","1-143"),
 ("Without_O2","SpinoHEMNatural_v2","O_spinosa","WT","9-165"),
 ("Without_O2","SpinoHEML43W_v2","O_spinosa","L43W","9-165"),
 ("Without_O2","SpinoHEMF125W_v2","O_spinosa","F125W","9-165"),
 ("Without_O2","SpinoHEMF125L43W_v2","O_spinosa","F125W+L43W","9-165"),
 ("With_O2","LegHEM22","P_sativum","WT","1-143"),
 ("With_O2","LegHEMV83T","P_sativum","V83T","1-143"),
 ("With_O2","SpinoHEMNatural","O_spinosa","WT","9-165"),
 ("With_O2","SpinoHEMF125W","O_spinosa","F125W","9-165"),
 ("With_O2","SpinoHEMF125L43W","O_spinosa","F125W+L43W","9-165"),
]

def files(d, rep):
    for stem in ([f"step5_prod_rep{rep}"] if rep > 1
                 else ["step5_prod", "step5_prod_rep1"]):
        t, x = f"{d}/{stem}.tpr", f"{d}/{stem}.xtc"
        if os.path.exists(t) and os.path.exists(x):
            return t, x
    return None, None

rows = []
for group, folder, sp, mut, core in SYS:
    d = f"{BASE}/{group}/{folder}/gromacs"
    for rep in (1, 2, 3):
        tpr, xtc = files(d, rep)
        if tpr is None:
            print(f"  SKIP {folder} rep{rep}"); continue
        u = mda.Universe(tpr, xtc)
        dt_ns = (u.trajectory[1].time - u.trajectory[0].time) / 1000.0
        start = int(SKIP_NS / dt_ns) if dt_ns > 0 else 0

        ref_core = f"protein and resid {core} and backbone"
        hem = "resname HEM or resname HEME"

        for label, sel in [("core", ref_core),
                           ("whole", "protein and backbone")]:
            R = RMSD(u, select=sel, groupselections=[hem]).run(start=start)
            vals = R.results.rmsd[:, 3]
            rows.append(dict(group=group, system=folder, species=sp,
                             mutation=mut, rep=rep, align=label,
                             hem_rmsd_A=vals.mean()))
        print(f"  {folder} rep{rep}: core={rows[-2]['hem_rmsd_A']:.2f}  "
              f"whole={rows[-1]['hem_rmsd_A']:.2f}")

df = pd.DataFrame(rows)
df.to_csv(f"{BASE}/results_hem_rmsd_alignment.csv", index=False)

print("\n" + "="*80)
piv = df.pivot_table(index=["group","species","mutation"],
                     columns="align", values="hem_rmsd_A",
                     aggfunc=["mean","std"])
print(piv.round(3))

from scipy import stats
print("\n--- تفاوت گونه‌ای (میانگین هر سیستم) ---")
for group in df.group.unique():
    for al in ("whole","core"):
        s = (df[(df.group==group)&(df.align==al)]
             .groupby(["system","species"])["hem_rmsd_A"].mean().reset_index())
        a = s.loc[s.species=="P_sativum","hem_rmsd_A"]
        b = s.loc[s.species=="O_spinosa","hem_rmsd_A"]
        if len(a) > 1 and len(b) > 1:
            t, p = stats.ttest_ind(a, b, equal_var=False)
            print(f"  {group:11s} align={al:6s} "
                  f"نخود={a.mean():.2f}  O.spinosa={b.mean():.2f}  p={p:.4f}")
