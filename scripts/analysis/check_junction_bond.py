"""
Diagnostic: verify the C(155)-N(156) bond parameters actually written to a
GROMACS topology file, by locating the exact atom serial numbers (not residue
numbers, which are not reliable after PSF->GROMACS conversion -- see
docs/KNOWN_ISSUES_AND_FIXES.md).

Expected correct values after the fix in fix_chain_break.py:
    b0 ~= 0.1345 nm, k ~= 309616 kJ/mol/nm^2

Usage:
    python check_junction_bond.py <system1_dir> <system2_dir> ...
Each <system_dir> must contain gromacs/topol_junction_fixed.top (or similar)
and gromacs/step3_input.pdb.
"""
import sys


def find_atom_index(pdb_path, resn, resi, name):
    idx = 0
    with open(pdb_path, errors='ignore') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                idx += 1
                if (line[17:20].strip() == resn and
                        line[22:26].strip() == str(resi) and
                        line[12:16].strip() == name):
                    return idx
    return None


def check_bond(top_path, pdb_path):
    c155 = find_atom_index(pdb_path, 'VAL', 155, 'C')
    n156 = find_atom_index(pdb_path, 'ASN', 156, 'N')
    if c155 is None or n156 is None:
        print(f"  [warning] atom not found (C155={c155}, N156={n156})")
        return

    with open(top_path) as f:
        in_bonds = False
        for line in f:
            if '[ bonds ]' in line:
                in_bonds = True
                continue
            if in_bonds and '[ pairs ]' in line:
                break
            if in_bonds:
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        a1, a2 = int(parts[0]), int(parts[1])
                        if {a1, a2} == {c155, n156}:
                            print(f"  found: {line.strip()}")
                            return
                    except ValueError:
                        continue
    print(f"  [error] bond C155({c155})-N156({n156}) not found in topology")


if __name__ == "__main__":
    for sys_dir in sys.argv[1:]:
        print(f"=== {sys_dir} ===")
        check_bond(f"{sys_dir}/gromacs/topol_junction_fixed.top", f"{sys_dir}/gromacs/step3_input.pdb")
