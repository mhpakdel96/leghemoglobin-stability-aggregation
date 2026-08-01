"""
Convert CHARMM-GUI PSF/PDB output to a correct GROMACS topology using ParmEd.

Bypasses CHARMM-GUI's built-in CHARMM->GROMACS converter, which was found to
produce an incorrect topology (duplicate molecule count, missing atoms) for
these systems. See docs/KNOWN_ISSUES_AND_FIXES.md, issue #3.

Input:
    step3_input.psf, step3_input.pdb  (CHARMM-GUI Solution Builder output)
    ../toppar/*.rtf, *.prm, *.str      (full CHARMM36 parameter set)

Output:
    topol_fixed.top, step3_input_fixed.gro
"""
import parmed as pmd
import glob

toppar_dir = '../toppar'
toppar_files = (sorted(glob.glob(f'{toppar_dir}/*.rtf')) +
                 sorted(glob.glob(f'{toppar_dir}/*.prm')) +
                 sorted(glob.glob(f'{toppar_dir}/*.str')))

print(f"Loading {len(toppar_files)} parameter files...")
params = pmd.charmm.CharmmParameterSet(*toppar_files)

struct = pmd.charmm.CharmmPsfFile('step3_input.psf')
struct.load_parameters(params)

coords = pmd.load_file('step3_input.pdb')
struct.coordinates = coords.coordinates

print(f"Total atoms: {len(struct.atoms)}")
print(f"Total charge: {sum(a.charge for a in struct.atoms):.3f}")

struct.save('topol_fixed.top', overwrite=True)
struct.save('step3_input_fixed.gro', overwrite=True)
print("Done.")
