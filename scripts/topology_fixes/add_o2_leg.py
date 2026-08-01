"""
O2 ligand construction and attachment at the heme iron for P. sativum systems.
No chain-junction repair needed (that issue is specific to O. spinosa; see
docs/KNOWN_ISSUES_AND_FIXES.md, issue #1).

Input:
    step3_input.psf, step3_input.pdb
    ../toppar/*.rtf, *.prm, *.str

Output:
    topol_o2_leg.top, step3_input_o2_leg.gro
"""
import parmed as pmd
from parmed import Atom, Bond, Angle, Dihedral, BondType, AngleType, DihedralType
import numpy as np
import glob

toppar_dir = '../toppar'
toppar_files = (sorted(glob.glob(f'{toppar_dir}/*.rtf')) +
                 sorted(glob.glob(f'{toppar_dir}/*.prm')) +
                 sorted(glob.glob(f'{toppar_dir}/*.str')))
params = pmd.charmm.CharmmParameterSet(*toppar_files)

struct = pmd.charmm.CharmmPsfFile('step3_input.psf')
struct.load_parameters(params)
coords = pmd.load_file('step3_input.pdb')
struct.coordinates = coords.coordinates


def find_atom(resname, name, resid=None):
    for a in struct.atoms:
        if a.residue.name == resname and a.name == name:
            if resid is None or a.residue.number == resid:
                return a
    raise ValueError(f"Atom not found: {name}/{resname}/{resid}")


PROXIMAL_HIS_RESID = 93   # 110 for O. spinosa systems

FE = find_atom('HEME', 'FE')
NA, NB = find_atom('HEME', 'NA'), find_atom('HEME', 'NB')
NC, ND = find_atom('HEME', 'NC'), find_atom('HEME', 'ND')

his_ne2 = None
for resn_try in ['HSD', 'HSE', 'HIS', 'HSP']:
    try:
        his_ne2 = find_atom(resn_try, 'NE2', PROXIMAL_HIS_RESID)
        break
    except ValueError:
        continue
if his_ne2 is None:
    raise ValueError(f"Proximal His{PROXIMAL_HIS_RESID} not found")

FE_xyz = np.array([FE.xx, FE.xy, FE.xz])
NA_xyz = np.array([NA.xx, NA.xy, NA.xz])
HIS_xyz = np.array([his_ne2.xx, his_ne2.xy, his_ne2.xz])

distal_dir = -(HIS_xyz - FE_xyz)
distal_dir /= np.linalg.norm(distal_dir)

O1_xyz = FE_xyz + distal_dir * 1.80
ref = NA_xyz - FE_xyz
ref_perp = ref - np.dot(ref, distal_dir) * distal_dir
ref_perp /= np.linalg.norm(ref_perp)
ang = np.radians(120.0)
O2_dir = np.cos(np.pi - ang) * distal_dir + np.sin(np.pi - ang) * ref_perp
O2_dir /= np.linalg.norm(O2_dir)
O2_xyz = O1_xyz + O2_dir * 1.23

O1_atom = Atom(name='O1', type='OM', charge=0.021, mass=15.999)
O1_atom.atomic_number = 8
O2_atom = Atom(name='O2', type='OM', charge=-0.021, mass=15.999)
O2_atom.atomic_number = 8
O1_atom.atom_type = params.atom_types['OM']
O2_atom.atom_type = params.atom_types['OM']

new_resnum = max(r.number for r in struct.residues) + 1
struct.add_atom(O1_atom, 'O2', new_resnum, chain='B')
struct.add_atom(O2_atom, 'O2', new_resnum, chain='B')
O1_atom.xx, O1_atom.xy, O1_atom.xz = O1_xyz
O2_atom.xx, O2_atom.xy, O2_atom.xz = O2_xyz

bt_o1o2 = BondType(k=600.0, req=1.23, list=struct.bond_types)
struct.bond_types.append(bt_o1o2)
struct.bonds.append(Bond(O1_atom, O2_atom, type=bt_o1o2))

bt_o1fe = BondType(k=250.0, req=1.80, list=struct.bond_types)
struct.bond_types.append(bt_o1fe)
struct.bonds.append(Bond(O1_atom, FE, type=bt_o1fe))

at_o2o1fe = AngleType(k=0.0, theteq=180.0, list=struct.angle_types)
struct.angle_types.append(at_o2o1fe)
struct.angles.append(Angle(O2_atom, O1_atom, FE, type=at_o2o1fe))

at_o1fen = AngleType(k=5.0, theteq=90.0, list=struct.angle_types)
struct.angle_types.append(at_o1fen)
for N_atom in [NA, NB, NC, ND]:
    struct.angles.append(Angle(O1_atom, FE, N_atom, type=at_o1fen))

dt_o2 = DihedralType(phi_k=0.0, per=4, phase=0.0, scee=1.0, scnb=1.0, list=struct.dihedral_types)
struct.dihedral_types.append(dt_o2)
struct.dihedrals.append(Dihedral(O2_atom, O1_atom, FE, NA, type=dt_o2))

struct.save('topol_o2_leg.top', overwrite=True)
struct.save('step3_input_o2_leg.gro', overwrite=True)
print(f"Total atoms: {len(struct.atoms)} | Total charge: {sum(a.charge for a in struct.atoms):.3f}")
print("Done.")
