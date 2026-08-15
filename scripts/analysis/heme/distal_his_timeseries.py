import MDAnalysis as mda, numpy as np, pandas as pd, os
from scipy import stats

BASE = "/cfs/earth/scratch/xpkk/Dr_Salmanian"
SKIP_NS = 20.0

SYS = [
 ("Without_O2","LegHEM22","P_sativum","WT",61),
 ("Without_O2","LegHEMV53T","P_sativum","V53T",61),
 ("Without_O2","LegHEMV83T","P_sativum","V83T",61),
 ("Without_O2","LegHEMV53V83T","P_sativum","V53T+V83T",61),
 ("Without_O2","SpinoHEMNatural_v2","O_spinosa","WT",75),
 ("Without_O2","SpinoHEML43W_v2","O_spinosa","L43W",75),
 ("Without_O2","SpinoHEMF125W_v2","O_spinosa","F125W",75),
 ("Without_O2","SpinoHEMF125L43W_v2","O_spinosa","F125W+L43W",75),
 ("With_O2","LegHEM22","P_sativum","WT",61),
 ("With_O2","LegHEMV83T","P_sativum","V83T",61),
 ("With_O2","SpinoHEMNatural","O_spinosa","WT",75),
 ("With_O2","SpinoHEMF125W","O_spinosa","F125W",75),
 ("With_O2","SpinoHEMF125L43W","O_spinosa","F125W+L43W",75),
]

def files(d, rep):
    stems = [f"step5_prod_rep{rep}"] if rep > 1 else ["step5_prod","step5_prod_rep1"]
    for s in stems:
        t, x = f"{d}/{s}.tpr", f"{d}/{s}.xtc"
        if os.path.exists(t) and os.path.exists(x):
            return t, x
    return None, None

printed = False
rows = []
for group, folder, sp, mut, dhis in SYS:
    d = f"{BASE}/{group}/{folder}/gromacs"
    for rep in (1, 2, 3):
        tpr, xtc = files(d, rep)
        if tpr is None:
            print(f"  SKIP {group}/{folder} rep{rep}"); continue
        u = mda.Universe(tpr, xtc)

        if not printed:  # یک‌بار: نام مولکول‌های غیرپروتئینی را نشان بده
            names = sorted({r.resname for r in u.residues
                            if r.resname not in ("SOL","TIP3","POT","CLA","NA","CL")})
            print("رزنیم‌های غیرحلال:", names); printed = True

        fe  = u.select_atoms("name FE")
        ne2 = u.select_atoms(f"protein and resid {dhis} and name NE2")
        o2  = u.select_atoms("resname O2 OXY LIG and name O1 O2 OX1 OX2")
        if not (len(fe) and len(ne2)):
            print(f"  {folder} rep{rep}: Fe یا NE2 نیست"); continue

        dt = (u.trajectory[1].time - u.trajectory[0].time)/1000.0
        start = int(SKIP_NS/dt) if dt > 0 else 0

        dfe, do2 = [], []
        for ts in u.trajectory[start:]:
            dfe.append(np.linalg.norm(fe.positions[0]-ne2.positions[0]))
            if len(o2):
                do2.append(min(np.linalg.norm(p-ne2.positions[0])
                               for p in o2.positions))
        dfe = np.array(dfe)
        row = dict(group=group, system=folder, species=sp, mutation=mut, rep=rep,
                   n_frames=len(dfe), dHis_Fe_mean=dfe.mean(), dHis_Fe_sd=dfe.std(),
                   pct_closed=100*np.mean(dfe < 5.0))
        if do2:
            do2 = np.array(do2)
            row["dHis_O2_mean"] = do2.mean()
            row["pct_hbond"]    = 100*np.mean(do2 < 3.5)
        rows.append(row)
        print(f"  {group[:9]:9s} {folder:22s} rep{rep}: "
              f"His{dhis}–Fe = {dfe.mean():5.2f}±{dfe.std():4.2f}  "
              f"بسته={row['pct_closed']:5.1f}%"
              + (f"  His–O2 = {row['dHis_O2_mean']:5.2f}  "
                 f"Hبند={row['pct_hbond']:5.1f}%" if do2 is not None and len(do2) else ""))

df = pd.DataFrame(rows)
df.to_csv(f"{BASE}/results_distal_his.csv", index=False)

print("\n" + "="*88)
cols = [c for c in ["dHis_Fe_mean","pct_closed","dHis_O2_mean","pct_hbond"]
        if c in df.columns]
print(df.groupby(["group","species","mutation"])[cols].agg(["mean","std"]).round(2))

print("\n--- مقایسه موتانت با WT ---")
for g in df.group.unique():
    for sp in df[df.group==g].species.unique():
        s = df[(df.group==g)&(df.species==sp)]
        if "WT" not in set(s.mutation): continue
        for m in cols:
            wt = s.loc[s.mutation=="WT", m].dropna()
            if len(wt) < 2: continue
            print(f"\n {g} · {sp} · {m}  (WT={wt.mean():.2f})")
            for mut in [x for x in s.mutation.unique() if x != "WT"]:
                v = s.loc[s.mutation==mut, m].dropna()
                if len(v) < 2: continue
                t, p = stats.ttest_ind(wt, v, equal_var=False)
                print(f"    {mut:12s} {v.mean():7.2f}  Δ={v.mean()-wt.mean():+7.2f}"
                      f"  p={p:.3f}{' *' if p<0.05 else ''}")
print("\nذخیره شد: results_distal_his.csv")
