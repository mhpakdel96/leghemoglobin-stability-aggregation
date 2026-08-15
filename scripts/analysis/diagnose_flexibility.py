import os, numpy as np

BASE = "/cfs/earth/scratch/xpkk/Dr_Salmanian"
SYS = [("Without_O2","LegHEM22","P_sativum"),
       ("Without_O2","SpinoHEMNatural_v2","O_spinosa")]

def read_rmsf(path):
    res, val = [], []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s[0] in "#@":
                continue
            p = s.split()
            try:
                res.append(int(float(p[0]))); val.append(float(p[1]) * 10)
            except (IndexError, ValueError):
                continue
    return np.array(res), np.array(val)

for group, folder, species in SYS:
    print("=" * 70); print(f"{species}  ({folder})"); print("=" * 70)
    for rep in (1, 2, 3):
        f = f"{BASE}/{group}/{folder}/gromacs/analysis_struct/rep{rep}_rmsf.xvg"
        if not os.path.exists(f):
            continue
        r, v = read_rmsf(f)
        print(f"\n rep{rep}: {len(r)} باقی‌مانده، RMSF کل = {v.mean():.2f} A")
        print("  ۱۰ باقی‌مانده با بیشترین نوسان:")
        for i in np.argsort(v)[-10:][::-1]:
            print(f"     resid {r[i]:4d}   {v[i]:6.2f} A")
        # میانگین در نواحی
        n = len(r)
        for lo, hi, lab in [(0, 10, "۱۰ باقی‌مانده اول"),
                            (n-15, n, "۱۵ باقی‌مانده آخر"),
                            (10, n-15, "هسته‌ی میانی")]:
            print(f"     {lab:22s} میانگین RMSF = {v[lo:hi].mean():.2f} A")
