#!/usr/bin/env bash
set -euo pipefail

# Sweep stochastic codeword length for the delimiter-free PFSA run-length code.
# Repeats independent transmissions at each L and aggregates BER, exact recovery,
# run-count error, and decoded-length error.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
RUNNER="${RUNNER:-$SCRIPT_DIR/pfsa_runlength_decode_reference.sh}"
OUT="${OUT:-$ROOT/testsuite/pfsa_runlength_sweep}"
LENGTHS="${LENGTHS:-125 250 500 750 1000 1500 2000}"
REPS="${REPS:-20}"
WIN_FRAC="${WIN_FRAC:-0.5}"
STEP_FRAC="${STEP_FRAC:-0.1}"
PLAINTEXT="${PLAINTEXT:-LSMASH STOCHASTIC CODE}"

rm -rf "$OUT"
mkdir -p "$OUT"
printf 'L\trep\tscan_win\tscan_step\ttrue_runs\tinferred_runs\ttrue_bits\tdecoded_bits\tedit_distance\tber_edit\texact\n' > "$OUT/results.tsv"

for L in $LENGTHS; do
  win=$(python3 - "$L" "$WIN_FRAC" <<'PY'
import sys
print(max(50, int(round(int(sys.argv[1])*float(sys.argv[2])))))
PY
)
  step=$(python3 - "$L" "$STEP_FRAC" <<'PY'
import sys
print(max(10, int(round(int(sys.argv[1])*float(sys.argv[2])))))
PY
)
  echo "=== L=$L win=$win step=$step ==="
  for rep in $(seq 1 "$REPS"); do
    d="$OUT/L_${L}/rep_$(printf '%03d' "$rep")"
    mkdir -p "$d"
    BASE_LEN="$L" SCAN_WIN="$win" SCAN_STEP="$step" OUT="$d" PLAINTEXT="$PLAINTEXT" "$RUNNER" > "$d/stdout.txt"
    python3 - "$d/runlength_summary.txt" "$L" "$rep" "$win" "$step" >> "$OUT/results.tsv" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1])
z={}
for line in p.read_text().splitlines():
    if "=" in line:
        k,v=line.split("=",1); z[k]=v
L,rep,win,step=map(int,sys.argv[2:])
true_bits=int(z["true_bits"]); edit=int(z["bit_edit_distance"])
print(L,rep,win,step,z["true_runs"],z["inferred_runs"],true_bits,z["decoded_bits"],
      edit,edit/max(1,true_bits),z["exact_bit_recovery"],sep="\t")
PY
    tail -n 8 "$d/runlength_summary.txt" | tr '\n' ' '; echo
  done
done

python3 - "$OUT/results.tsv" "$OUT/summary.tsv" <<'PY'
import csv,sys,statistics,math
inp,out=sys.argv[1:]
rows=list(csv.DictReader(open(inp),delimiter="\t"))
groups={}
for r in rows: groups.setdefault(int(r["L"]),[]).append(r)
with open(out,"w") as f:
    f.write("L\tn\texact_rate\tmean_edit\tmean_ber_edit\tmedian_ber_edit\tmean_abs_run_count_error\tmean_abs_bit_count_error\n")
    for L in sorted(groups):
        g=groups[L]
        exact=[int(r["exact"]) for r in g]
        edit=[int(r["edit_distance"]) for r in g]
        ber=[float(r["ber_edit"]) for r in g]
        re=[abs(int(r["inferred_runs"])-int(r["true_runs"])) for r in g]
        be=[abs(int(r["decoded_bits"])-int(r["true_bits"])) for r in g]
        f.write(f"{L}\t{len(g)}\t{sum(exact)/len(g):.6f}\t{statistics.mean(edit):.6f}\t"
                f"{statistics.mean(ber):.8f}\t{statistics.median(ber):.8f}\t"
                f"{statistics.mean(re):.6f}\t{statistics.mean(be):.6f}\n")
print(open(out).read(),end="")
PY

echo "Raw:     $OUT/results.tsv"
echo "Summary: $OUT/summary.tsv"
