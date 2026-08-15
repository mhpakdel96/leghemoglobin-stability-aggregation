#!/bin/bash
set -uo pipefail
GMX=/cfs/earth/scratch/xpkk/.conda/envs/gromacs_gpu/bin.AVX2_256/gmx
BASE=/cfs/earth/scratch/xpkk/Dr_Salmanian
SKIP_PS=20000
DT_PS=100

# هسته: نخود = ۱۴۳ باقی‌مانده اول (قطعه جداشده حذف)
#        اسپینوزا = ۹ تا ۱۶۵ (دُم N و دو باقی‌مانده آخر حذف)
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
  D="$BASE/$GROUP/$DIR/gromacs"
  [[ -d $D ]] || { echo "MISSING: $GROUP/$DIR"; continue; }
  cd "$D" || continue
  A=analysis_struct; C=analysis_core
  mkdir -p $C
  echo "--- $GROUP/$DIR  (core ri $RANGE) ---"

  for REP in 1 2 3; do
    # تراژکتوری پروتئین: اگر نیست، بساز (هر دو الگوی نام‌گذاری)
    PX=$A/prot_rep$REP.xtc
    if [[ ! -s $PX ]]; then
      for CAND in "step5_prod_rep$REP" "step5_prod"; do
        if [[ -s $CAND.xtc && -s $CAND.tpr ]]; then
          printf "Protein\n" | $GMX trjconv -s $CAND.tpr -f $CAND.xtc \
            -n $A/ana.ndx -o $PX -pbc whole -b $SKIP_PS -dt $DT_PS \
            >/dev/null 2>&1 && break
        fi
      done
    fi
    [[ -s $PX ]] || { echo "   rep$REP: تراژکتوری نیست"; continue; }

    # گروه هسته
    if [[ ! -s $C/core.ndx ]]; then
      printf "ri $RANGE\nname 0 CoreAll\nq\n" \
        | $GMX make_ndx -f $A/prot.tpr -o $C/tmp.ndx >/dev/null 2>&1
      printf "ri $RANGE & a CA\nri $RANGE & a N C CA O\nq\n" \
        | $GMX make_ndx -f $A/prot.tpr -o $C/core.ndx >/dev/null 2>&1
    fi

    NB=$(grep -c "^\[" $C/core.ndx)
    CORE_BB=$((NB-1)); CORE_CA=$((NB-2)); 

    P=$C/rep$REP
    printf "$CORE_BB\n$CORE_BB\n" | $GMX rms  -s $A/prot.tpr -f $PX -n $C/core.ndx -o ${P}_rmsd.xvg -tu ns >/dev/null 2>&1
    printf "$CORE_CA\n"           | $GMX rmsf -s $A/prot.tpr -f $PX -n $C/core.ndx -o ${P}_rmsf.xvg -res -fit >/dev/null 2>&1
    printf "$CORE_BB\n"           | $GMX gyrate -s $A/prot.tpr -f $PX -n $C/core.ndx -o ${P}_gyrate.xvg >/dev/null 2>&1
    $GMX sasa -s $A/prot.tpr -f $PX -n $C/core.ndx -surface "group $CORE_BB" -o ${P}_sasa.xvg >/dev/null 2>&1
    printf "$CORE_BB\n$CORE_BB\n" | $GMX hbond -s $A/prot.tpr -f $PX -n $C/core.ndx -num ${P}_hbnum.xvg >/dev/null 2>&1
    cp -n $A/rep${REP}_dssp.dat ${P}_dssp.dat 2>/dev/null

    OK=0; for f in rmsd rmsf gyrate sasa hbnum; do [[ -s ${P}_${f}.xvg ]] && OK=$((OK+1)); done
    echo "   rep$REP → $OK/5"
  done
done
echo "=== پایان ==="
