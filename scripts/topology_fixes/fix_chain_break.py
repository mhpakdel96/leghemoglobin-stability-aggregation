"""
Reconstruct the C(155)-N(156) peptide bond artificially split by CHARMM-GUI's
chain-detection heuristic in O. spinosa systems.

See docs/KNOWN_ISSUES_AND_FIXES.md, issue #1, for the diagnostic evidence that
this junction is fully continuous in the source structure and that the split
is purely a CHARMM-GUI parsing artifact.

All bonded-term parameters below are explicit numeric values read directly
from par_all36m_prot.prm (NOT retrieved via automated ParameterSet lookup --
see docs/KNOWN_ISSUES_AND_FIXES.md, issue #4, for why automated lookup was
found to be unreliable for custom-added terms).

Input:
    step3_input.psf, step3_input.pdb
    ../toppar/*.rtf, *.prm, *.str

Output:
    topol_junction_fixed.top, step3_input_junction_fixed.gro
"""
import parmed as pmd
from parmed import Bond, Angle, Dihedral, BondType, AngleType, DihedralType
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


def find_atom(resname, name, resid):
    for a in struct.atoms:
        if a.residue.name == resname and a.name == name and a.residue.number == resid:
            return a
    raise ValueError(f"Atom not found: {name}/{resname}/{resid}")


# Residue numbers below (155/156) are specific to these O. spinosa systems;
# adjust if reusing this script for a different chain-break location.
C155, O155, CA155 = find_atom('VAL', 'C', 155), find_atom('VAL', 'O', 155), find_atom('VAL', 'CA', 155)
N156, HN156, CA156 = find_atom('ASN', 'N', 156), find_atom('ASN', 'HN', 156), find_atom('ASN', 'CA', 156)

# Bond C-N
bt = BondType(k=370.000, req=1.3450, list=struct.bond_types)
struct.bond_types.append(bt)
struct.bonds.append(Bond(C155, N156, type=bt))

# Angles
at1 = AngleType(k=80.0, theteq=116.5, list=struct.angle_types)   # CA-C-N
struct.angle_types.append(at1)
struct.angles.append(Angle(CA155, C155, N156, type=at1))

at2 = AngleType(k=80.0, theteq=122.5, list=struct.angle_types)   # O-C-N
struct.angle_types.append(at2)
struct.angles.append(Angle(O155, C155, N156, type=at2))

at3 = AngleType(k=50.0, theteq=120.0, list=struct.angle_types)   # C-N-CA
struct.angle_types.append(at3)
struct.angles.append(Angle(C155, N156, CA156, type=at3))

at4 = AngleType(k=34.0, theteq=123.0, list=struct.angle_types)   # C-N-H
struct.angle_types.append(at4)
struct.angles.append(Angle(C155, N156, HN156, type=at4))

# Omega backbone dihedral (two Fourier terms)
dt1 = DihedralType(phi_k=1.6, per=1, phase=0.0, scee=1.0, scnb=1.0, list=struct.dihedral_types)
struct.dihedral_types.append(dt1)
struct.dihedrals.append(Dihedral(CA155, C155, N156, CA156, type=dt1))

dt2 = DihedralType(phi_k=2.5, per=2, phase=180.0, scee=1.0, scnb=1.0, list=struct.dihedral_types)
struct.dihedral_types.append(dt2)
struct.dihedrals.append(Dihedral(CA155, C155, N156, CA156, type=dt2))

struct.save('topol_junction_fixed.top', overwrite=True)
struct.save('step3_input_junction_fixed.gro', overwrite=True)
print("Chain junction repaired and saved.")
