#!/bin/bash
set -uo pipefail
GMX=/cfs/earth/scratch/xpkk/.conda/envs/gromacs_gpu/bin.AVX2_256/gmx
BASE=/cfs/earth/scratch/xpkk/Dr_Salmanian

# همان محدوده‌ی هسته که در تحلیل‌های قبلی استفاده شد (ri = residue index)
SYSTEMS=(
"Without_O2 LegHEM22             1 143"
"Without_O2 LegHEMV53T           1 143"
"Without_O2 LegHEMV83T           1 143"
"Without_O2 LegHEMV53V83T        1 143"
"Without_O2 SpinoHEMNatural_v2   9 165"
"Without_O2 SpinoHEML43W_v2      9 165"
"Without_O2 SpinoHEMF125W_v2     9 165"
"Without_O2 SpinoHEMF125L43W_v2  9 165"
"With_O2 LegHEM22                1 143"
"With_O2 LegHEMV83T              1 143"
"With_O2 SpinoHEMNatural         9 165"
"With_O2 SpinoHEMF125W           9 165"
"With_O2 SpinoHEMF125L43W        9 165"
)

echo "=== شروع: $(date) ==="
for S in "${SYSTEMS[@]}"; do
  read -r GROUP DIR R1 R2 <<< "$S"
  D="$BASE/$GROUP/$DIR/gromacs"
  A="$D/analysis_struct"; C="$D/analysis_core"
  if [[ ! -f "$A/prot.tpr" ]]; then
    echo "SKIP $GROUP/$DIR — prot.tpr نیست"; continue
  fi
  mkdir -p "$C"
  echo "--- $GROUP/$DIR  (ri $R1..$R2) ---"

  SURF="resindex $R1 to $R2 and not name \"H*\""
  APOL="\"Apolar\" resindex $R1 to $R2 and name \"C*\" \"S*\""
  POL="\"Polar\"  resindex $R1 to $R2 and name \"N*\" \"O*\""

  for REP in 1 2 3; do
    PX="$A/prot_rep$REP.xtc"
    if [[ ! -s $PX ]]; then echo "   rep$REP: تراژکتوری نیست"; continue; fi
    OUT="$C/rep${REP}_sasa_split.xvg"
    ERR="$C/rep${REP}_sasa_split.log"

    $GMX sasa -s "$A/prot.tpr" -f "$PX" \
        -surface "$SURF" \
        -output  "$APOL" "$POL" \
        -o "$OUT" \
        -or "$C/rep${REP}_sasa_perres.xvg" > "$ERR" 2>&1

    if [[ -s $OUT ]]; then
      NC=$(grep -v '^[#@]' "$OUT" | head -1 | awk '{print NF}')
      NF_=$(grep -vc '^[#@]' "$OUT")
      echo "   rep$REP → ok   ستون‌ها=$NC   فریم‌ها=$NF_"
    else
      echo "   rep$REP → ✗ ناموفق (لاگ: $ERR)"
    fi
  done
done
echo "=== پایان: $(date) ==="
