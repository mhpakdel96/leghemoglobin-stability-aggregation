import os, glob, math
from collections import OrderedDict

AA = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
      'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
      'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
      'HSD':'H','HSE':'H','HSP':'H'}

def parse(path):
    """→ (residues OrderedDict[(chain,resseq)] = resname, atoms list, header lines)"""
    res, atoms, hdr = OrderedDict(), [], []
    if not os.path.exists(path):
        return None
    for line in open(path, errors="ignore"):
        rec = line[:6]
        if rec.startswith(("HEADER","TITLE ","COMPND","SOURCE","DBREF ","SEQRES")):
            if not rec.startswith("SEQRES"):
                hdr.append(line.rstrip())
        if rec in ("ATOM  ", "HETATM"):
            name  = line[12:16].strip()
            rname = line[17:20].strip()
            chain = line[21].strip() or "_"
            try:
                seq = int(line[22:26])
                x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                b = float(line[60:66]) if len(line) > 66 else 0.0
            except ValueError:
                continue
            key = (chain, seq)
            if key not in res:
                res[key] = rname
            atoms.append((chain, seq, rname, name, x, y, z, b))
    return res, atoms, hdr

def seq_of(res):
    return "".join(AA.get(v, "X") for v in res.values())

def dist(a, b):
    return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

print("#" * 78)
print("# بخش ۱ — هویت مدل‌های AlphaFold")
print("#" * 78)
for f in sorted(glob.glob("legHb_project/AF-*.pdb")):
    out = parse(f)
    if not out: continue
    res, atoms, hdr = out
    prot = OrderedDict((k,v) for k,v in res.items() if v in AA)
    plddt = [a[7] for a in atoms if a[3] == "CA"]
    print(f"\n--- {f} ---")
    for h in hdr[:8]:
        print("   ", h)
    print(f"    باقی‌مانده‌های پروتئینی : {len(prot)}")
    if prot:
        ks = list(prot.keys())
        print(f"    شماره‌گذاری            : {ks[0][1]} تا {ks[-1][1]}")
        print(f"    اولین باقی‌مانده       : {prot[ks[0]]}")
    if plddt:
        print(f"    pLDDT میانگین          : {sum(plddt)/len(plddt):.1f}")
        low = [i+1 for i,v in enumerate(plddt) if v < 70]
        print(f"    باقی‌مانده با pLDDT<70 : {len(low)}  {low[:25]}")
    s = seq_of(prot)
    for i in range(0, len(s), 60):
        print(f"    {i+1:4d} {s[i:i+60]}")

print("\n" + "#" * 78)
print("# بخش ۲ — کدام هیستیدین به آهن نزدیک است؟")
print("#" * 78)
for f in sorted(glob.glob("*/*_docked.pdb")) + sorted(glob.glob("legHb_project/*HEM*.pdb")):
    out = parse(f)
    if not out: continue
    res, atoms, _ = out
    fe = [a for a in atoms if a[3] == "FE"]
    if not fe:
        continue
    fx = fe[0][4:7]
    hits = []
    for a in atoms:
        if a[2] in ("HIS","HSD","HSE","HSP") and a[3] in ("NE2","ND1"):
            hits.append((dist(a[4:7], fx), a[1], a[3]))
    hits.sort()
    print(f"\n--- {f} ---")
    for d, rid, an in hits[:4]:
        tag = "  ← پروگزیمال" if d < 3.5 else ("  ← دیستال" if d < 6 else "")
        print(f"    His{rid:<4d} {an}   فاصله تا Fe = {d:6.2f} A{tag}")

print("\n" + "#" * 78)
print("# بخش ۳ — آیا هر موتانت واقعاً و فقط جهش موردنظر را دارد؟")
print("#" * 78)
PAIRS = [
 ("legHb_project/LegHEM22.pdb", ["legHb_project/LegHEMV53T.pdb",
   "legHb_project/LegHEMV83T.pdb", "legHb_project/LegHEMV53V83T.pdb"]),
 ("legHb_project/SpinoHEMNatural.pdb", ["legHb_project/SpinoHEMF125W.pdb",
   "legHb_project/SpinoHEMF125WL43W.pdb"]),
 ("legHb_project/Object1_Repair.pdb", ["legHb_project/Object-V53T1_1.pdb",
   "legHb_project/ObjectV83T.pdb", "legHb_project/ObjectV53T-V83T.pdb",
   "legHb_project/ObjectL43W.pdb", "legHb_project/ObjectF125W.pdb",
   "legHb_project/ObjectF125L43W.pdb"]),
]
for wt_f, muts in PAIRS:
    wo = parse(wt_f)
    if not wo:
        print(f"\n{wt_f}: یافت نشد"); continue
    wres = OrderedDict((k,v) for k,v in wo[0].items() if v in AA)
    print(f"\n=== مرجع: {wt_f}  ({len(wres)} باقی‌مانده) ===")
    for m in muts:
        mo = parse(m)
        if not mo:
            print(f"   {os.path.basename(m):32s} یافت نشد"); continue
        mres = OrderedDict((k,v) for k,v in mo[0].items() if v in AA)
        if len(mres) != len(wres):
            print(f"   {os.path.basename(m):32s} ⚠ طول متفاوت: "
                  f"{len(mres)} در برابر {len(wres)}")
        diffs = [(k[1], wres[k], mres[k]) for k in wres
                 if k in mres and wres[k] != mres[k]]
        if not diffs:
            print(f"   {os.path.basename(m):32s} ⚠ هیچ تفاوتی با WT ندارد!")
        else:
            txt = ", ".join(f"{AA.get(a,'?')}{p}{AA.get(b,'?')}"
                            for p,a,b in diffs)
            print(f"   {os.path.basename(m):32s} {len(diffs)} تفاوت: {txt}")

print("\n" + "#" * 78)
print("# بخش ۴ — آیا مراحل پاکسازی/تعمیر چیزی را عوض کرده‌اند؟")
print("#" * 78)
CHAIN = ["legHb_project/AF-Q9SAZ0-F1-model_v6.pdb",
         "legHb_project/AF-Q9SAZ0-F1-model_v6.pdb-clean-clean.pdb",
         "legHb_project/Object1_Repair.pdb",
         "legHb_project/LegHEM22.pdb"]
prev = None
for f in CHAIN:
    o = parse(f)
    if not o:
        print(f"{os.path.basename(f):46s} یافت نشد"); continue
    r = OrderedDict((k,v) for k,v in o[0].items() if v in AA)
    ks = list(r.keys())
    s = seq_of(r)
    flag = ""
    if prev is not None:
        if s == prev[0]:
            flag = "  ✓ توالی یکسان با مرحله قبل"
        else:
            flag = "  ⚠ توالی عوض شد!"
    print(f"{os.path.basename(f):46s} {len(r):4d} باقی‌مانده  "
          f"({ks[0][1]}–{ks[-1][1]}){flag}")
    prev = (s,)

print("\n" + "#" * 78)
print("# بخش ۵ — قالب‌های PDB")
print("#" * 78)
for f in ["legHb_project/1FSL.pdb", "legHb_project/1GDJ.pdb"]:
    o = parse(f)
    if not o:
        print(f"{f}: یافت نشد"); continue
    for h in o[2][:6]:
        print("   ", h)
    print()
