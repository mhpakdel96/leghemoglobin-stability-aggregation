# WARNING - SUPERSEDED. RESULTS FROM THIS SCRIPT WERE WITHDRAWN.
#
# This script produced conclusions that were later retracted:
#   - variant ranking based on a single replicate
#   - 'stability' of Fe-O1 and Fe-His distances, which are harmonic
#     restraints written into the topology, not simulation results
#   - a species difference arising from superposition on a disordered tail
#
# See docs/CORRECTIONS.md and docs/KNOWN_ISSUES.md items 4, 5 and 7.
# Superseded by scripts/analysis/heme/hem_rmsd_core.py
# Retained only so retracted numbers can be traced to their source.

"""
Fe-O2 bond stability and heme+O2 RMSD analysis for the O2-bound systems.

Reports, per system:
  - Fe-O1 distance: mean / std / max over the trajectory (bond dissociation
    check; a stable oxy-heme should stay near the CHARMM36 equilibrium value
    of 1.80 A throughout)
  - Heme+O2 RMSD relative to the backbone-aligned protein (see
    analyze_hem_stability.py docstring for the alignment method)

Expects each system directory to contain step5_prod.tpr / step5_prod.xtc.
"""
import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
import os

SYSTEMS = [
    ("LegHEM22",        "With_O2/LegHEM22/gromacs"),
    ("LegHEMV83T",      "With_O2/LegHEMV83T/gromacs"),
    ("SpinoHEMNatural", "With_O2/SpinoHEMNatural/gromacs"),
    ("SpinoHEMF125W",   "With_O2/SpinoHEMF125W/gromacs"),
]

DISSOCIATION_THRESHOLD_A = 3.5  # Fe-O1 distance above which O2 is considered dissociated


def main():
    print(f"{'System':20s} {'Fe-O1 mean':>12s} {'std':>8s} {'max':>8s} {'dissociated?':>14s} {'RMSD':>8s}")
    print("-" * 80)

    for name, folder in SYSTEMS:
        tpr = os.path.join(folder, "step5_prod.tpr")
        xtc = os.path.join(folder, "step5_prod.xtc")
        if not (os.path.exists(tpr) and os.path.exists(xtc)):
            print(f"{name:20s}  [not available]")
            continue

        u = mda.Universe(tpr, xtc)
        fe = u.select_atoms('resname HEME and name FE')
        o1 = u.select_atoms('resname O2 and name O1')
        protein_bb = u.select_atoms('protein and backbone')
        hem_o2 = u.select_atoms('resname HEME or resname O2')

        if len(fe) != 1 or len(o1) != 1:
            print(f"{name:20s}  [Fe or O1 not uniquely found]")
            continue

        ref_bb = protein_bb.positions.copy()
        ref_hemo2 = hem_o2.positions.copy()

        fe_o1_dists, rmsd_vals = [], []
        for _ in u.trajectory:
            fe_o1_dists.append(np.linalg.norm(fe.positions[0] - o1.positions[0]))
            mobile_bb = protein_bb.positions
            R, _ = align.rotation_matrix(mobile_bb - mobile_bb.mean(axis=0),
                                          ref_bb - ref_bb.mean(axis=0))
            hemo2_aligned = (hem_o2.positions - mobile_bb.mean(axis=0)) @ R.T + ref_bb.mean(axis=0)
            rmsd_vals.append(np.sqrt(np.mean(np.sum((hemo2_aligned - ref_hemo2) ** 2, axis=1))))

        fe_o1_dists = np.array(fe_o1_dists)
        rmsd_vals = np.array(rmsd_vals)
        dissociated = "yes" if fe_o1_dists.max() > DISSOCIATION_THRESHOLD_A else "no"

        print(f"{name:20s} {fe_o1_dists.mean():12.2f} {fe_o1_dists.std():8.2f} "
              f"{fe_o1_dists.max():8.2f} {dissociated:>14s} {rmsd_vals.mean():8.2f}")


if __name__ == "__main__":
    main()
