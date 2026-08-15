#!/bin/bash
set -uo pipefail
GMX=/cfs/earth/scratch/xpkk/.conda/envs/gromacs_gpu/bin.AVX2_256/gmx
BASE=/cfs/earth/scratch/xpkk/Dr_Salmanian

SYSTEMS=(
"Without_O2 LegHEM22            1-143"
"Without_O2 LegHEMV53T          1-143"
"Without_O2 LegHEMV83T          1-143"
"Without_O2 LegHEMV53V83T       1-143"
"Without_O2 SpinoHEMNatural_v2  9-165"
"Without_O2 SpinoHEML43W_v2     9-165"
"Without_O2 SpinoHEMF125W_v2    9-165"
"Without_O2 SpinoHEMF125L43W_v2 9-165"
"With_O2 LegHEM22               1-143"
"With_O2 LegHEMV83T             1-143"
"With_O2 SpinoHEMNatural        9-165"
"With_O2 SpinoHEMF125W          9-165"
"With_O2 SpinoHEMF125L43W       9-165"
)

for S in "${SYSTEMS[@]}"; do
  read -r GROUP DIR RANGE <<< "$S"
  D="$BASE/$GROUP/$DIR/gromacs"; A=$D/analysis_struct; C=$D/analysis_core
  [[ -d $C ]] || { echo "MISSING: $GROUP/$DIR"; continue; }
  cd "$D" || continue
  echo "--- $GROUP/$DIR (ri $RANGE، همه اتم‌ها) ---"

  printf "ri $RANGE\nq\n" | $GMX make_ndx -f $A/prot.tpr -o $C/coreall.ndx \
      >/dev/null 2>&1 || { echo "   FAIL make_ndx"; continue; }
  NB=$(grep -c "^\[" $C/coreall.ndx); G=$((NB-1))
  echo "   گروه هسته = $G ($(grep "^\[" $C/coreall.ndx | tail -1 | tr -d '[] '))"

  for REP in 1 2 3; do
    PX=$A/prot_rep$REP.xtc
    [[ -s $PX ]] || { echo "   rep$REP: تراژکتوری نیست"; continue; }
    P=$C/rep$REP

    $GMX sasa -s $A/prot.tpr -f $PX -n $C/coreall.ndx \
        -surface "group $G" -o ${P}_sasa.xvg              >/dev/null 2>&1
    printf "$G\n"       | $GMX gyrate -s $A/prot.tpr -f $PX \
        -n $C/coreall.ndx -o ${P}_gyrate.xvg              >/dev/null 2>&1
    printf "$G\n$G\n"   | $GMX hbond  -s $A/prot.tpr -f $PX \
        -n $C/coreall.ndx -num ${P}_hbnum.xvg             >/dev/null 2>&1

    OK=""
    for f in sasa gyrate hbnum; do
      [[ -s ${P}_${f}.xvg ]] && OK="$OK $f:ok" || OK="$OK $f:FAIL"
    done
    echo "   rep$REP →$OK"
  done
done
echo "=== پایان ==="
