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
Heme-group RMSD relative to the backbone-aligned protein, computed over the
full 100 ns production trajectory of each system.

At each frame, the protein backbone is aligned (rotation only, about its own
centroid) to the first production frame; heme-atom RMSD is then computed
against their first-frame positions in this aligned reference. This isolates
heme displacement relative to the protein scaffold from overall rigid-body
tumbling of the protein in the simulation box.

Expects each system directory to contain step5_prod.tpr / step5_prod.xtc
(GROMACS production output).
"""
import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
import os

# (label, path-to-gromacs-directory, proximal-His-resid)
SYSTEMS = [
    ("LegHEM22",         "LegHEM22/gromacs",         93),
    ("LegHEMV53T",       "LegHEMV53T/gromacs",       93),
    ("LegHEMV83T",       "LegHEMV83T/gromacs",       93),
    ("LegHEMV53V83T",    "LegHEMV53V83T/gromacs",    93),
    ("SpinoHEMNatural",  "SpinoHEMNatural/gromacs",  110),
    ("SpinoHEML43W",     "SpinoHEML43W/gromacs",     110),
    ("SpinoHEMF125W",    "SpinoHEMF125W/gromacs",    110),
    ("SpinoHEMF125L43W", "SpinoHEMF125L43W/gromacs", 110),
]


def compute_heme_rmsd(tpr, xtc):
    u = mda.Universe(tpr, xtc)
    protein_bb = u.select_atoms('protein and backbone')
    hem = u.select_atoms('resname HEME')
    if len(hem) == 0:
        raise ValueError("HEME not found")

    ref_bb = protein_bb.positions.copy()
    ref_hem = hem.positions.copy()

    rmsd_values = []
    for _ in u.trajectory:
        mobile_bb = protein_bb.positions
        R, _ = align.rotation_matrix(mobile_bb - mobile_bb.mean(axis=0),
                                      ref_bb - ref_bb.mean(axis=0))
        hem_aligned = (hem.positions - mobile_bb.mean(axis=0)) @ R.T + ref_bb.mean(axis=0)
        rmsd_values.append(np.sqrt(np.mean(np.sum((hem_aligned - ref_hem) ** 2, axis=1))))
    return np.array(rmsd_values)


def main():
    results = {}
    print(f"{'System':22s} {'mean (A)':>10s} {'std':>8s} {'max':>8s}")
    print("-" * 55)
    for name, folder, _ in SYSTEMS:
        tpr = os.path.join(folder, "step5_prod.tpr")
        xtc = os.path.join(folder, "step5_prod.xtc")
        if not (os.path.exists(tpr) and os.path.exists(xtc)):
            print(f"{name:22s}  [not available]")
            continue
        rmsd = compute_heme_rmsd(tpr, xtc)
        results[name] = rmsd
        print(f"{name:22s} {rmsd.mean():10.2f} {rmsd.std():8.2f} {rmsd.max():8.2f}")

    print("\nRanked by mean RMSD (lower = more stable):")
    for name, rmsd in sorted(results.items(), key=lambda kv: kv[1].mean()):
        print(f"  {name:22s} {rmsd.mean():.2f} A")


if __name__ == "__main__":
    main()
