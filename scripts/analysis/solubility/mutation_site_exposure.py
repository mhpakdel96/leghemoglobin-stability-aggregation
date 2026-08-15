import numpy as np, pandas as pd, os

BASE = "/cfs/earth/scratch/xpkk/Dr_Salmanian"

# بیشینه SASA باقی‌مانده آزاد (Tien et al. 2013، nm²) برای نرمال‌سازی
MAXSASA = {"ALA":1.29,"ARG":2.74,"ASN":1.95,"ASP":1.93,"CYS":1.67,
           "GLN":2.25,"GLU":2.23,"GLY":1.04,"HIS":2.24,"ILE":1.97,
           "LEU":2.01,"LYS":2.36,"MET":2.24,"PHE":2.40,"PRO":1.59,
           "SER":1.55,"THR":1.72,"TRP":2.85,"TYR":2.63,"VAL":1.74}

# (گروه، پوشه، برچسب، محدوده هسته، موقعیت‌های هدف، باقی‌مانده اصلی)
SYS = [
 ("Without_O2","LegHEM22","نخود WT",(1,143),
  {53:"VAL", 83:"VAL"}),
 ("Without_O2","LegHEMV53T","نخود V53T",(1,143),
  {53:"THR", 83:"VAL"}),
 ("Without_O2","LegHEMV83T","نخود V83T",(1,143),
  {53:"VAL", 83:"THR"}),
 ("Without_O2","LegHEMV53V83T","نخود دوتایی",(1,143),
  {53:"THR", 83:"THR"}),
 ("Without_O2","SpinoHEMNatural_v2","اسپینوزا WT",(9,165),
  {43:"LEU", 125:"PHE"}),
 ("Without_O2","SpinoHEML43W_v2","اسپینوزا L43W",(9,165),
  {43:"TRP", 125:"PHE"}),
 ("Without_O2","SpinoHEMF125W_v2","اسپینوزا F125W",(9,165),
  {43:"LEU", 125:"TRP"}),
 ("Without_O2","SpinoHEMF125L43W_v2","اسپینوزا دوتایی",(9,165),
  {43:"TRP", 125:"TRP"}),
]

def read_perres(path):
    """gmx sasa -or : ستون اول شماره باقی‌مانده، ستون دوم SASA میانگین"""
    d = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return d
    for line in open(path):
        s = line.strip()
        if not s or s[0] in "#@":
            continue
        p = s.split()
        try:
            d[int(float(p[0]))] = float(p[1])
        except (IndexError, ValueError):
            continue
    return d

rows = []
for group, folder, label, (r1, r2), targets in SYS:
    C = f"{BASE}/{group}/{folder}/gromacs/analysis_core"
    per = []
    for rep in (1, 2, 3):
        d = read_perres(f"{C}/rep{rep}_sasa_perres.xvg")
        if d: per.append(d)
    if not per:
        print(f"  داده نیست: {folder}"); continue

    keys = sorted(set(per[0]))
    for pos, resn in targets.items():
        # gmx sasa -or ممکن است از ۱ شماره‌گذاری کند؛ نگاشت بر اساس ri
        idx = pos - r1 + 1 if (pos - r1 + 1) in per[0] else pos
        vals = [p[idx] for p in per if idx in p]
        if not vals:
            print(f"  موقعیت {pos} در {folder} پیدا نشد"); continue
        m, s = float(np.mean(vals)), float(np.std(vals))
        rel = 100 * m / MAXSASA.get(resn, 2.0)
        rows.append(dict(system=folder, label=label, pos=pos, resname=resn,
                         sasa_nm2=m, sd=s, rel_pct=rel,
                         state=("دفن‌شده" if rel < 20 else
                                "نیمه‌نمایان" if rel < 50 else "نمایان")))

df = pd.DataFrame(rows)
df.to_csv(f"{BASE}/results_mutation_site_exposure.csv", index=False)

print("="*92)
print("نمایان بودن باقی‌مانده‌های هدف (rel_pct = درصد نسبت به بیشینه ممکن آن آمینواسید)")
print("="*92)
for _, r in df.iterrows():
    print(f"  {r.label:20s} {r.resname}{r.pos:<5d} SASA = {r.sasa_nm2:5.3f} ± {r.sd:5.3f} nm²"
          f"   نسبی = {r.rel_pct:5.1f}٪   → {r.state}")

print("\n" + "="*92)
print("تغییر با جهش (نسبت به WT همان گونه)")
print("="*92)
for sp, wt in [("Leg","LegHEM22"), ("Spino","SpinoHEMNatural_v2")]:
    w = df[df.system == wt]
    for sysname in df[df.system.str.startswith(sp) & (df.system != wt)].system.unique():
        m = df[df.system == sysname]
        for pos in sorted(set(w.pos)):
            a = w[w.pos == pos]; b = m[m.pos == pos]
            if a.empty or b.empty: continue
            a, b = a.iloc[0], b.iloc[0]
            tag = " ← جهش‌یافته" if a.resname != b.resname else ""
            print(f"  {sysname:22s} موقعیت {pos:<4d} {a.resname}→{b.resname}  "
                  f"SASA {a.sasa_nm2:5.3f} → {b.sasa_nm2:5.3f}  "
                  f"Δ={b.sasa_nm2-a.sasa_nm2:+6.3f} nm²{tag}")

print("\nذخیره شد: results_mutation_site_exposure.csv")
