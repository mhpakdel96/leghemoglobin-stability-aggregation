import os
from collections import OrderedDict

AA3TO1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
          'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
          'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

JOBS = [
 ("pea", "P_sativum",
  "../legHb_project/P.sativum/AF-Q9SAZ0-F1-model_v6.pdb",
  [("WT",        []),
   ("V53T",      [(53,'V','T')]),
   ("V83T",      [(83,'V','T')]),
   ("V53T_V83T", [(53,'V','T'), (83,'V','T')])]),
 ("spi", "O_spinosa",
  "../legHb_project/AF-A0A411AFI2-F1-model_v6.pdb",
  [("WT",         []),
   ("L43W",       [(43,'L','W')]),
   ("F125W",      [(125,'F','W')]),
   ("L43W_F125W", [(43,'L','W'), (125,'F','W')])]),
]

def seq_from_pdb(path):
    res = OrderedDict()
    for l in open(path, errors="ignore"):
        if l[:6] == "ATOM  ":
            res[int(l[22:26])] = l[17:20].strip()
    ks = sorted(res)
    if ks != list(range(ks[0], ks[-1] + 1)):
        raise SystemExit(f"✗ شماره‌گذاری {path} پیوسته نیست")
    return ks[0], [AA3TO1.get(res[k], 'X') for k in ks]

ok = True
for tag, species, pdb, variants in JOBS:
    if not os.path.exists(pdb):
        print(f"✗ {pdb} پیدا نشد"); ok = False; continue
    start, wt = seq_from_pdb(pdb)
    print(f"\n{species}: {len(wt)} باقی‌مانده، شماره‌گذاری از {start}")

    for name, muts in variants:
        s = list(wt)
        for pos, old, new in muts:
            i = pos - start
            if s[i] != old:
                print(f"  ✗ {name}: موقعیت {pos} = {s[i]} است نه {old}")
                ok = False; continue
            s[i] = new
            print(f"  ✓ {name}: {old}{pos}{new} اعمال شد")
        seq = "".join(s)
        fn = f"{tag}_{name}.fasta"
        with open(fn, "w") as fh:
            fh.write(f">{species}_{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")
        # تأیید: چند تفاوت با WT دارد؟
        ndiff = sum(1 for a, b in zip(seq, "".join(wt)) if a != b)
        print(f"     → {fn}  ({len(seq)} aa، {ndiff} تفاوت با WT)")

print("\n" + ("همه چیز درست است." if ok else "⚠ خطا هست — ادامه ندهید."))
