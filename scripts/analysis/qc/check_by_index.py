import numpy as np
import MDAnalysis as mda

def read_atoms(path):
    atoms=[]
    with open(path, errors='ignore') as fh:
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atoms.append(dict(
                    idx=len(atoms),
                    name=line[12:16].strip(),
                    resn=line[17:20].strip(),
                    resi=line[22:26].strip(),
                ))
    return atoms

orig_atoms = read_atoms('step3_input.pdb')

c155_idx = None
n156_idx = None
for a in orig_atoms:
    if a['resn']=='VAL' and a['resi']=='155' and a['name']=='C':
        c155_idx = a['idx']
    if a['resn']=='ASN' and a['resi']=='156' and a['name']=='N':
        n156_idx = a['idx']

print(f"موقعیت C155 در فایل اصلی: {c155_idx}")
print(f"موقعیت N156 در فایل اصلی: {n156_idx}")

u = mda.Universe('step4.0_min_full.gro')
all_atoms = u.atoms

c155_pos = all_atoms[c155_idx].position
n156_pos = all_atoms[n156_idx].position

d = np.linalg.norm(c155_pos - n156_pos)
print(f"\n>>> فاصله‌ی واقعی C155-N156 (با موقعیت اتم): {d:.2f} آنگستروم")
print(f"    نام واقعی این دو اتم: {all_atoms[c155_idx].name}/{all_atoms[c155_idx].resname} و {all_atoms[n156_idx].name}/{all_atoms[n156_idx].resname}")

