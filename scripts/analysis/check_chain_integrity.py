import MDAnalysis as mda
import numpy as np, os

BASE = "/cfs/earth/scratch/xpkk/Dr_Salmanian"
CASES = [
 ("Without_O2","LegHEM22","نخود WT"),
 ("Without_O2","LegHEMV83T","نخود V83T"),
 ("Without_O2","SpinoHEMNatural_v2","O.spinosa WT"),
]

for group, folder, label in CASES:
    D = f"{BASE}/{group}/{folder}/gromacs/analysis_struct"
    tpr, xtc = f"{D}/prot.tpr", f"{D}/prot_rep1.xtc"
    if not (os.path.exists(tpr) and os.path.exists(xtc)):
        print(f"\n### {label}: فایل نیست"); continue
    u = mda.Universe(tpr, xtc)
    print("\n" + "="*72); print(f"### {label}  ({folder})"); print("="*72)

    print(f"تعداد کل باقی‌مانده‌ها: {u.residues.n_residues}")
    for seg in u.segments:
        r = seg.residues
        print(f"  segment {seg.segid!r}: {r.n_residues} باقی‌مانده، "
              f"resid {r.resids.min()}–{r.resids.max()}")

    # جستجوی شکست زنجیره: فاصله C(i) تا N(i+1)
    print("\nشکست‌های زنجیره (فاصله C→N بیش از 2 آنگستروم):")
    res = u.residues
    breaks = []
    for i in range(len(res) - 1):
        c = res[i].atoms.select_atoms("name C")
        n = res[i+1].atoms.select_atoms("name N")
        if len(c) != 1 or len(n) != 1:
            continue
        d = np.linalg.norm(c.positions[0] - n.positions[0])
        if d > 2.0:
            breaks.append((res[i].segid, res[i].resid, res[i+1].segid,
                           res[i+1].resid, d))
    if not breaks:
        print("   هیچ شکستی نیست — زنجیره پیوسته است")
    for sa, ra, sb, rb, d in breaks:
        print(f"   {sa}:{ra}  →  {sb}:{rb}   فاصله = {d:8.2f} A")

    # اگر شکست هست، فاصله را در طول تراژکتوری دنبال کن
    if breaks:
        sa, ra, sb, rb, _ = breaks[0]
        c = u.select_atoms(f"segid {sa} and resid {ra} and name C")
        n = u.select_atoms(f"segid {sb} and resid {rb} and name N")
        if len(c) == 1 and len(n) == 1:
            ds = [np.linalg.norm(c.positions[0]-n.positions[0])
                  for _ in u.trajectory]
            ds = np.array(ds)
            print(f"\n   این فاصله در طول شبیه‌سازی: "
                  f"شروع={ds[0]:.1f}  کمینه={ds.min():.1f}  "
                  f"بیشینه={ds.max():.1f}  پایان={ds[-1]:.1f} A")
            print("   → اگر بیشینه خیلی بزرگ است، قطعه واقعاً جدا شده")
