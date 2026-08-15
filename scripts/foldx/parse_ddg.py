import os, re
import numpy as np

# نگاشت شماره‌ی جهش FoldX به نام، از روی همان فایل individual_list
JOBS = [
 ("pea", "P. sativum", "individual_list_pea.txt",
  "Average_pea_WT_Repair.fxout", "Dif_pea_WT_Repair.fxout"),
 ("spi", "O. spinosa", "individual_list_spi.txt",
  "Average_spi_WT_Repair.fxout", "Dif_spi_WT_Repair.fxout"),
]

def read_mutlist(path):
    """→ {1: 'V53T', 2: 'V83T', 3: 'V53T+V83T'}"""
    out = {}
    if not os.path.exists(path):
        return out
    i = 0
    for line in open(path):
        s = line.strip().rstrip(";")
        if not s:
            continue
        i += 1
        parts = []
        for m in s.split(","):
            m = m.strip()
            # قالب: <وحشی><زنجیره><موقعیت><موتانت>  مثلا VA53T
            g = re.match(r"^([A-Z])([A-Za-z0-9])(\d+)([A-Z])$", m)
            parts.append(f"{g.group(1)}{g.group(3)}{g.group(4)}" if g else m)
        out[i] = "+".join(parts)
    return out

def read_table(path):
    """→ (فهرست ستون‌ها، فهرست سطرها به‌صورت [نام، اعداد...])"""
    if not os.path.exists(path):
        return None, []
    cols, rows, started = [], [], False
    for line in open(path, errors="ignore"):
        s = line.rstrip("\n")
        if s.startswith("Pdb\t"):
            cols = s.split("\t"); started = True; continue
        if not started or not s.strip():
            continue
        p = s.split("\t")
        rows.append(p)
    return cols, rows

def verdict(v):
    if v >  0.5: return "ناپایدارکننده"
    if v < -0.5: return "پایدارکننده"
    return "خنثی (زیر آستانه دقت ۰.۵)"

for tag, species, mutfile, avgfile, diffile in JOBS:
    muts = read_mutlist(mutfile)
    print("\n" + "=" * 84)
    print(f"{species}   —   {avgfile}")
    print("=" * 84)
    if not muts:
        print(f"  ✗ {mutfile} پیدا نشد"); continue

    cols, rows = read_table(avgfile)
    if not rows:
        print(f"  ✗ {avgfile} پیدا نشد یا خالی است"); continue

    # ستون‌ها در فایل Average: Pdb | SD | total energy | ...
    try:
        i_sd  = cols.index("SD")
        i_tot = cols.index("total energy")
    except ValueError:
        print(f"  ✗ ستون‌های مورد انتظار پیدا نشد: {cols[:4]}"); continue

    means = {}
    for p in rows:
        m = re.search(r"_(\d+)$", p[0].strip())
        if not m: continue
        idx = int(m.group(1))
        try:
            mean, sd = float(p[i_tot]), float(p[i_sd])
        except (ValueError, IndexError):
            continue
        means[idx] = (mean, sd)
        name = muts.get(idx, f"جهش {idx}")
        print(f"  {name:14s} ΔΔG = {mean:+7.3f} ± {sd:6.4f} kcal/mol"
              f"   → {verdict(mean)}")

    # اجراهای تکی از فایل Dif
    dcols, drows = read_table(diffile)
    if drows:
        try:
            j_tot = dcols.index("total energy")
        except ValueError:
            j_tot = 1
        per = {}
        for p in drows:
            m = re.search(r"_(\d+)_(\d+)\.pdb$", p[0].strip())
            if not m: continue
            per.setdefault(int(m.group(1)), []).append(float(p[j_tot]))
        print()
        for idx in sorted(per):
            name = muts.get(idx, f"جهش {idx}")
            print(f"  {name:14s} اجراها: "
                  f"{[round(v, 3) for v in per[idx]]}")

    # جمع‌پذیری
    if {1, 2, 3}.issubset(means):
        s1, s2, dbl = means[1][0], means[2][0], means[3][0]
        print(f"\n  جمع‌پذیری: {muts.get(1)} ({s1:+.3f}) + {muts.get(2)} ({s2:+.3f})"
              f" = {s1+s2:+.3f}")
        print(f"  {'':13s}دوتایی مشاهده‌شده = {dbl:+.3f}"
              f"   اپیستازی = {dbl-(s1+s2):+.3f} kcal/mol")

print("\n" + "-" * 84)
print("قرارداد علامت FoldX: مثبت = ناپایدارکننده، منفی = پایدارکننده.")
print("دقت روش ~۰.۵ kcal/mol — هر |ΔΔG| زیر این مقدار عملاً صفر است.")
