import numpy as np, pandas as pd, os
from scipy import stats

BASE = "/cfs/earth/scratch/xpkk/Dr_Salmanian"

# (گروه، پوشه، گونه، جهش، تعداد جهش)
SYS = [
 ("Without_O2","LegHEM22",            "P_sativum","WT",         0),
 ("Without_O2","LegHEMV53T",          "P_sativum","V53T",       1),
 ("Without_O2","LegHEMV83T",          "P_sativum","V83T",       1),
 ("Without_O2","LegHEMV53V83T",       "P_sativum","V53T+V83T",  2),
 ("Without_O2","SpinoHEMNatural_v2",  "O_spinosa","WT",         0),
 ("Without_O2","SpinoHEML43W_v2",     "O_spinosa","L43W",       1),
 ("Without_O2","SpinoHEMF125W_v2",    "O_spinosa","F125W",      1),
 ("Without_O2","SpinoHEMF125L43W_v2", "O_spinosa","F125W+L43W", 2),
 ("With_O2","LegHEM22",         "P_sativum","WT",         0),
 ("With_O2","LegHEMV83T",       "P_sativum","V83T",       1),
 ("With_O2","SpinoHEMNatural",  "O_spinosa","WT",         0),
 ("With_O2","SpinoHEMF125W",    "O_spinosa","F125W",      1),
 ("With_O2","SpinoHEMF125L43W", "O_spinosa","F125W+L43W", 2),
]

def read_xvg(path):
    """→ (آرایه‌ی داده، فهرست legendها)"""
    rows, leg = [], []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None, leg
    for line in open(path):
        s = line.strip()
        if s.startswith("@"):
            if " legend " in s and s.startswith("@ s"):
                leg.append(s.split('"')[1] if '"' in s else s)
            continue
        if not s or s.startswith("#"):
            continue
        try:
            rows.append([float(v) for v in s.split()])
        except ValueError:
            continue
    return (np.array(rows) if rows else None), leg

# --- تأیید یک‌باره‌ی ترتیب ستون‌ها ---
probe = f"{BASE}/Without_O2/LegHEM22/gromacs/analysis_core/rep1_sasa_split.xvg"
d, leg = read_xvg(probe)
print("legendهای فایل نمونه:", leg)
if d is None or d.shape[1] < 4:
    raise SystemExit("✗ فایل نمونه ۴ ستون ندارد — گام ۲ را بررسی کنید")
chk = abs((d[:,2] + d[:,3]) - d[:,1]).mean()
print(f"میانگین |(ستون۳+ستون۴) − ستون۲| = {chk:.4f} nm²  "
      f"→ {'✓ ستون‌ها درست‌اند' if chk < 0.5 else '✗ ترتیب ستون‌ها مشکوک است'}")
if chk >= 0.5:
    raise SystemExit("متوقف شد — ترتیب ستون‌ها را بررسی کنید")

rows = []
for group, folder, sp, mut, nmut in SYS:
    for rep in (1, 2, 3):
        p = f"{BASE}/{group}/{folder}/gromacs/analysis_core/rep{rep}_sasa_split.xvg"
        d, _ = read_xvg(p)
        if d is None or d.shape[1] < 4:
            print(f"  خالی: {group}/{folder} rep{rep}"); continue
        tot, apol, pol = d[:,1].mean(), d[:,2].mean(), d[:,3].mean()
        rows.append(dict(group=group, system=folder, species=sp, mutation=mut,
                         n_mut=nmut, rep=rep, n_frames=len(d),
                         sasa_total=tot, sasa_apolar=apol, sasa_polar=pol,
                         frac_apolar=100*apol/tot))

df = pd.DataFrame(rows)
df.to_csv(f"{BASE}/results_sasa_split.csv", index=False)

METRICS = ["sasa_total","sasa_apolar","sasa_polar","frac_apolar"]

print("\n" + "="*100)
print("میانگین سه تکرار (nm²؛ frac_apolar بر حسب درصد)")
print("="*100)
summ = (df.groupby(["group","species","n_mut","mutation"])[METRICS]
          .agg(["mean","std"]).round(2))
print(summ)

print("\n" + "="*100)
print("آزمون روند: رگرسیون خطی معیار بر حسب تعداد جهش")
print("="*100)
out = []
for group in df.group.unique():
    for sp in df[df.group==group].species.unique():
        s = df[(df.group==group) & (df.species==sp)]
        if s.n_mut.nunique() < 2:
            continue
        print(f"\n--- {group} · {sp}  (n={len(s)} نقطه، سطوح جهش: "
              f"{sorted(s.n_mut.unique())}) ---")
        for m in METRICS:
            x = s.n_mut.values.astype(float); y = s[m].values
            r = stats.linregress(x, y)
            star = " *" if r.pvalue < 0.05 else ""
            print(f"   {m:12s} شیب = {r.slope:+7.3f} ± {r.stderr:.3f} "
                  f"به ازای هر جهش   R²={r.rvalue**2:5.3f}   p={r.pvalue:.4f}{star}")
            out.append(dict(group=group, species=sp, metric=m, slope=r.slope,
                            stderr=r.stderr, r2=r.rvalue**2, p=r.pvalue))

            # بررسی جمع‌پذیری: میانگین سطح ۱ در برابر میانگین (سطح۰ + سطح۲)/۲
            g = s.groupby("n_mut")[m].mean()
            if set([0,1,2]).issubset(set(g.index)):
                pred = (g[0] + g[2]) / 2
                print(f"   {'':12s}   جمع‌پذیری: مشاهده در ۱ جهش = {g[1]:.3f}، "
                      f"پیش‌بینی خطی = {pred:.3f}، اختلاف = {g[1]-pred:+.3f}")

pd.DataFrame(out).to_csv(f"{BASE}/results_trend_stats.csv", index=False)
print("\nذخیره شد: results_sasa_split.csv و results_trend_stats.csv")
