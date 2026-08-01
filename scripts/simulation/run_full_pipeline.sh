#!/bin/bash
# SLURM job: minimization -> equilibration -> 100 ns production, GROMACS on GPU.
#
# Usage:
#   sbatch --job-name=<label> run_full_pipeline.sh <SYSTEM_DIR> <TOPOL> <MINGRO> <STARTGRO> [PROD_NS]
#
# Example:
#   sbatch --job-name=LegV53T run_full_pipeline.sh LegHEMV53T topol_fixed.top \
#       step3_input_fixed.gro step3_input_fixed.gro 100

#SBATCH --partition=earth-5
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=16:00:00
#SBATCH --output=slurm_%j_%x.out
#SBATCH --error=slurm_%j_%x.err

SYSTEM_DIR=$1
TOPOL=$2
MINGRO=$3      # starting .gro for minimization
STARTGRO=$4    # reference .gro for position restraints (-r flag)
PROD_NS=${5:-100}

GMX=/path/to/gromacs_gpu/bin/gmx   # adjust to your cluster's GROMACS binary
BASE=/path/to/MD_systems/$SYSTEM_DIR/gromacs
cd "$BASE"

nvidia-smi

echo "=== Minimization ==="
$GMX grompp -f step4.0_minimization.mdp -o step4.0_min.tpr \
    -c "$MINGRO" -r "$STARTGRO" -p "$TOPOL" -n index.ndx -maxwarn 5
$GMX mdrun -v -deffnm step4.0_min -ntomp 16 -ntmpi 1

echo "=== Equilibration ==="
$GMX grompp -f step4.1_equilibration.mdp -o step4.1_equil.tpr \
    -c step4.0_min.gro -r "$STARTGRO" -p "$TOPOL" -n index.ndx -maxwarn 5
$GMX mdrun -v -deffnm step4.1_equil -nb gpu -pme gpu -ntomp 16 -ntmpi 1

echo "=== Production (${PROD_NS} ns) ==="
NSTEPS=$(python3 -c "print(int($PROD_NS * 1000 / 0.002))")
cp step5_production.mdp step5_production_${PROD_NS}ns.mdp
sed -i "s/^nsteps.*=.*/nsteps                  = $NSTEPS/" step5_production_${PROD_NS}ns.mdp

$GMX grompp -f step5_production_${PROD_NS}ns.mdp -o step5_prod.tpr \
    -c step4.1_equil.gro -p "$TOPOL" -n index.ndx -maxwarn 5
$GMX mdrun -v -deffnm step5_prod -nb gpu -pme gpu -ntomp 16 -ntmpi 1

echo "=== Done: $SYSTEM_DIR ==="
grep "Performance" step5_prod.log
