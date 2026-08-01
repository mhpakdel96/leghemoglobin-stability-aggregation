import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.distances import distance_array

u = mda.Universe('step4.0_min_full.gro')

protein = u.select_atoms('protein')
c_atoms = protein.select_atoms('name C')
n_atoms = protein.select_atoms('name N')

d_matrix = distance_array(c_atoms.positions, n_atoms.positions)

print("جزئیات نقاط شکسته:")
for i in range(len(c_atoms)):
    min_d = d_matrix[i].min()
    if min_d > 1.6:
        c_atom = c_atoms[i]
        nearest_n_idx = d_matrix[i].argmin()
        nearest_n = n_atoms[nearest_n_idx]
        print(f"  C در رزیدو {c_atom.resid} ({c_atom.resname}) → نزدیک‌ترین N: رزیدو {nearest_n.resid} ({nearest_n.resname}), فاصله={min_d:.2f} A")
