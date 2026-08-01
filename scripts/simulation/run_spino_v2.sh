#!/bin/bash
#SBATCH --job-name=spino_v2
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

SYSTEM_DIR=$1   # مثلا SpinoHEMNatural_v2

GMX=/path/to/gromacs_gpu/bin/gmx
cd /path/to/MD_systems/Without_O2/$SYSTEM_DIR/gromacs

nvidia-smi

echo "=== ساخت index جدید ==="
$GMX make_ndx -f step3_input_junction_v2.gro -o index_new.ndx << 'EOF'
q
EOF

echo "=== Minimization ==="
$GMX grompp -f step4.0_minimization.mdp -o step4.0_min.tpr \
  -c step3_input_junction_v2.gro -r step3_input_junction_v2.gro \
  -p topol_junction_v2.top -n index.ndx -maxwarn 5
$GMX mdrun -v -deffnm step4.0_min -ntomp 16 -ntmpi 1

echo "=== Equilibration ==="
$GMX grompp -f step4.1_equilibration.mdp -o step4.1_equil.tpr \
  -c step4.0_min.gro -r step3_input_junction_v2.gro \
  -p topol_junction_v2.top -n index.ndx -maxwarn 5
$GMX mdrun -v -deffnm step4.1_equil -nb gpu -pme gpu -ntomp 16 -ntmpi 1

echo "=== Production (100 ns) ==="
NSTEPS=$(python3 -c "print(int(100 * 1000 / 0.002))")
cp step5_production.mdp step5_production_100ns.mdp
sed -i "s/^nsteps.*=.*/nsteps                  = $NSTEPS/" step5_production_100ns.mdp
$GMX grompp -f step5_production_100ns.mdp -o step5_prod.tpr \
  -c step4.1_equil.gro -p topol_junction_v2.top -n index.ndx -maxwarn 5
$GMX mdrun -v -deffnm step5_prod -nb gpu -pme gpu -ntomp 16 -ntmpi 1

echo "=== تمام شد: $SYSTEM_DIR ==="
grep "Performance" step5_prod.log
