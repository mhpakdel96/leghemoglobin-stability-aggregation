#!/bin/bash
set -uo pipefail

GMX=/cfs/earth/scratch/xpkk/.conda/envs/gromacs_gpu/bin.AVX2_256/gmx
BASE=/cfs/earth/scratch/xpkk/Dr_Salmanian
SKIP_PS=20000     # 20 ns اول (تعادل‌یابی) حذف می‌شود
DT_PS=100         # هر 100 ps یک فریم

SYSTEMS=(
"Without_O2 LegHEM22"
"Without_O2 LegHEMV53T"
"Without_O2 LegHEMV83T"
"Without_O2 LegHEMV53V83T"
"Without_O2 SpinoHEMNatural_v2"
"Without_O2 SpinoHEML43W_v2"
"Without_O2 SpinoHEMF125W_v2"
"Without_O2 SpinoHEMF125L43W_v2"
"With_O2 LegHEM22"
"With_O2 LegHEMV83T"
"With_O2 SpinoHEMNatural"
"With_O2 SpinoHEMF125W"
"With_O2 SpinoHEMF125L43W"
)

analyze_one () {
  local GROUP=$1 DIR=$2 REP=$3
  local D="$BASE/$GROUP/$DIR/gromacs"
  [[ -d $D ]] || { echo "  MISSING DIR: $GROUP/$DIR"; return 0; }
  cd "$D" || return 0

  local XTC TPR
  if [[ $REP == 1 ]]; then
    XTC=step5_prod.xtc;        TPR=step5_prod.tpr
  else
    XTC=step5_prod_rep$REP.xtc; TPR=step5_prod_rep$REP.tpr
  fi
  [[ -f $XTC && -f $TPR ]] || { echo "  SKIP (no traj): $GROUP/$DIR rep$REP"; return 0; }

  local O=analysis_struct
  mkdir -p $O

  if [[ ! -f $O/ana.ndx ]]; then
    printf "q\n" | $GMX make_ndx -f $TPR -o $O/ana.ndx >/dev/null 2>&1 \
      || { echo "  FAIL make_ndx: $DIR"; return 1; }
  fi

  if [[ ! -f $O/prot.tpr ]]; then
    printf "Protein\n" | $GMX convert-tpr -s $TPR -n $O/ana.ndx -o $O/prot.tpr >/dev/null 2>&1 \
      || { echo "  FAIL convert-tpr: $DIR"; return 1; }
  fi

  local PX=$O/prot_rep$REP.xtc
  if [[ ! -f $PX ]]; then
    printf "Protein\n" | $GMX trjconv -s $TPR -f $XTC -n $O/ana.ndx -o $PX \
        -pbc whole -b $SKIP_PS -dt $DT_PS >/dev/null 2>&1 \
      || { echo "  FAIL trjconv: $DIR rep$REP"; return 1; }
  fi
  [[ -s $PX ]] || { echo "  FAIL empty traj: $DIR rep$REP"; return 1; }

  local P=$O/rep$REP
  printf "Backbone\nBackbone\n" | $GMX rms    -s $O/prot.tpr -f $PX -o ${P}_rmsd.xvg  -tu ns    >/dev/null 2>&1
  printf "C-alpha\n"            | $GMX rmsf   -s $O/prot.tpr -f $PX -o ${P}_rmsf.xvg  -res -fit >/dev/null 2>&1
  printf "Protein\n"            | $GMX gyrate -s $O/prot.tpr -f $PX -o ${P}_gyrate.xvg          >/dev/null 2>&1
  $GMX sasa -s $O/prot.tpr -f $PX -surface 'all' -o ${P}_sasa.xvg                                >/dev/null 2>&1
  printf "Protein\nProtein\n"   | $GMX hbond  -s $O/prot.tpr -f $PX -num ${P}_hbnum.xvg          >/dev/null 2>&1
  $GMX dssp -s $O/prot.tpr -f $PX -o ${P}_dssp.dat >/dev/null 2>&1 || true

  local OK=0
  for f in rmsd rmsf gyrate sasa hbnum; do
    [[ -s ${P}_${f}.xvg ]] && OK=$((OK+1))
  done
  echo "  $GROUP/$DIR rep$REP  →  $OK/5 xvg  dssp:$([[ -s ${P}_dssp.dat ]] && echo yes || echo no)"
}

echo "=== شروع: $(date) ==="
for S in "${SYSTEMS[@]}"; do
  read -r G D <<< "$S"
  echo "--- $G / $D ---"
  for R in 1 2 3; do
    analyze_one "$G" "$D" "$R"
  done
done
echo "=== پایان: $(date) ==="
