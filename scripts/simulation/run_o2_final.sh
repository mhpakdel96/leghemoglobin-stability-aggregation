#!/bin/bash
#SBATCH --job-name=o2_final
#SBATCH --partition=earth-5
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=16:00:00
#SBATCH --constraint=rhel8
#SBATCH --output=slurm_%j_%x.out
#SBATCH --error=slurm_%j_%x.err

SYSTEM_DIR=$1      # مثلا SpinoHEMNatural
TOPOL=$2           # topol_full_fixed.top
MINGRO=$3          # step4.0_min_full2.gro
STARTGRO=$4        # step3_input_full_fixed.gro

GMX=/path/to/gromacs_gpu/bin/gmx
cd /path/to/MD_systems/With_O2/$SYSTEM_DIR/gromacs

echo "=== سیستم: $SYSTEM_DIR ==="
nvidia-smi

echo "=== ساخت index جدید (شامل O2) ==="
$GMX make_ndx -f $MINGRO -o index_new.ndx << 'EOF'
1 | 13 | 14
15 | 16 | 17
name 18 SOLU
name 19 SOLV
q
EOF

echo "=== Equilibration ==="
$GMX grompp -f step4.1_equilibration.mdp -o step4.1_equil.tpr \
  -c $MINGRO -r $STARTGRO -p $TOPOL -n index_new.ndx -maxwarn 5
$GMX mdrun -v -deffnm step4.1_equil -nb gpu -pme gpu -ntomp 16 -ntmpi 1

echo "=== Production (100 ns) ==="
$GMX grompp -f step5_production_100ns.mdp -o step5_prod.tpr \
  -c step4.1_equil.gro -p $TOPOL -n index_new.ndx -maxwarn 5
$GMX mdrun -v -deffnm step5_prod -nb gpu -pme gpu -ntomp 16 -ntmpi 1

echo "=== تمام شد: $SYSTEM_DIR ==="
grep "Performance" step5_prod.log
