import pandas as pd, numpy as np
from scipy import stats

BASE = "/cfs/earth/scratch/xpkk/Dr_Salmanian"
df = pd.read_csv(f"{BASE}/results_core_raw.csv")
METRICS = ["rmsd_bb_A","rmsf_ca_A","rg_nm","sasa_nm2","hbonds","helix_pct"]

out = []
for group in df.group.unique():
    g = df[df.group == group]
    print("="*95); print(group); print("="*95)
    for sp in g.species.unique():
        s = g[g.species == sp]
        if "WT" not in set(s.mutation):
            continue
        for m in METRICS:
            wt = s.loc[s.mutation == "WT", m].dropna().values
            if len(wt) < 2:
                continue
            print(f"\n--- {sp} · {m} · WT = {wt.mean():.2f} ---")
            for mut in [x for x in s.mutation.unique() if x != "WT"]:
                v = s.loc[s.mutation == mut, m].dropna().values
                if len(v) < 2:
                    continue
                t, p = stats.ttest_ind(wt, v, equal_var=False)
                d = v.mean() - wt.mean()
                sd = np.sqrt((wt.var(ddof=1) + v.var(ddof=1)) / 2)
                coh = d / sd if sd > 0 else np.nan
                flag = " *" if p < 0.05 else ""
                print(f"   {mut:12s} {v.mean():8.2f}  Δ={d:+7.3f}  d={coh:+5.2f}  p={p:.3f}{flag}")
                out.append(dict(group=group, species=sp, metric=m, mutation=mut,
                                wt=wt.mean(), mut_val=v.mean(), delta=d,
                                cohens_d=coh, p=p))

    # مقایسه بین دو گونه، روی میانگین سیستم‌ها (نه تک‌تک تکرارها)
    print(f"\n--- {group}: مقایسه بین دو گونه (میانگین هر سیستم، بدون شبه‌تکرار) ---")
    sysmean = g.groupby(["system","species"])[METRICS].mean().reset_index()
    sps = sysmean.species.unique()
    if len(sps) == 2:
        for m in METRICS:
            a = sysmean.loc[sysmean.species == sps[0], m].dropna()
            b = sysmean.loc[sysmean.species == sps[1], m].dropna()
            if len(a) < 2 or len(b) < 2:
                continue
            t, p = stats.ttest_ind(a, b, equal_var=False)
            print(f"   {m:12s} {sps[0]}={a.mean():7.2f}  {sps[1]}={b.mean():7.2f}  p={p:.4f}")
    print()

pd.DataFrame(out).to_csv(f"{BASE}/results_core_stats.csv", index=False)
print("ذخیره شد: results_core_stats.csv")
