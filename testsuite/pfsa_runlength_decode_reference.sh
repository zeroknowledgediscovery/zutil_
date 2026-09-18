#!/usr/bin/env bash
set -euo pipefail

# Run-length PFSA code:
# maximal 0-runs -> one A realization of length r*BASE_LEN
# maximal 1-runs -> one B realization of length r*BASE_LEN
# Decoder sees only concatenated transmission + independent A/B references.
# It continuously scans with lsmash, infers A/B residence intervals, then
# decodes each residence duration as an integer multiple of BASE_LEN.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PRUN="${PRUN:-$ROOT/bin/prun}"
LSMASH="${LSMASH:-$ROOT/bin/lsmash}"
MODEL_A="${MODEL_A:-$ROOT/testsuite/pfsa_opt/M2.cfg}"
MODEL_B="${MODEL_B:-$ROOT/testsuite/pfsa_opt/T3.cfg}"

OUT="${OUT:-$ROOT/testsuite/pfsa_runlength_decode_run}"
BASE_LEN="${BASE_LEN:-500}"
REF_LEN="${REF_LEN:-30000}"
SCAN_WIN="${SCAN_WIN:-250}"
SCAN_STEP="${SCAN_STEP:-50}"
SCAN_BATCH="${SCAN_BATCH:-20}"
MIN_RUN_WINDOWS="${MIN_RUN_WINDOWS:-2}"
PLAINTEXT="${PLAINTEXT:-LSMASH STOCHASTIC CODE}"

rm -rf "$OUT"
mkdir -p "$OUT/scan_batches"

python3 - "$PLAINTEXT" "$OUT" <<'PY'
from pathlib import Path
import sys
text,out=sys.argv[1],Path(sys.argv[2])
bits="".join(f"{b:08b}" for b in text.encode())
(out/"message_bits.txt").write_text(bits+"\n")
runs=[]
i=0
while i<len(bits):
    j=i+1
    while j<len(bits) and bits[j]==bits[i]: j+=1
    runs.append((len(runs),bits[i],j-i))
    i=j
with (out/"run_plan.tsv").open("w") as f:
    f.write("run\tbit\tcount\n")
    for r,b,n in runs: f.write(f"{r}\t{b}\t{n}\n")
PY

"$PRUN" -f "$MODEL_A" -l "$REF_LEN" -o "$OUT/ref_A.dat"
"$PRUN" -f "$MODEL_B" -l "$REF_LEN" -o "$OUT/ref_B.dat"
: > "$OUT/transmission.dat"
printf 'run\tbit\tcount\tstart\tend\tlength\n' > "$OUT/run_truth.tsv"
cursor=0
while IFS=$'\t' read -r run bit count; do
    [[ "$run" == run ]] && continue
    len=$((count * BASE_LEN))
    if [[ "$bit" == 0 ]]; then model="$MODEL_A"; else model="$MODEL_B"; fi
    tmp="$OUT/.run_${run}.dat"
    "$PRUN" -f "$model" -l "$len" -o "$tmp"
    tr '\n' ' ' < "$tmp" >> "$OUT/transmission.dat"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$run" "$bit" "$count" "$cursor" "$((cursor+len))" "$len" >> "$OUT/run_truth.tsv"
    cursor=$((cursor+len))
    rm -f "$tmp"
done < "$OUT/run_plan.tsv"
printf '\n' >> "$OUT/transmission.dat"

python3 - "$OUT" "$SCAN_WIN" "$SCAN_STEP" <<'PY'
from pathlib import Path
import sys
out=Path(sys.argv[1]); w=int(sys.argv[2]); step=int(sys.argv[3])
x=(out/"transmission.dat").read_text().split()
with (out/"scan_windows.dat").open("w") as fw,(out/"scan_meta.tsv").open("w") as fm:
    fm.write("index\tstart\tend\tcenter\n")
    for k,s in enumerate(range(0,len(x)-w+1,step)):
        e=s+w
        fw.write(" ".join(x[s:e])+"\n")
        fm.write(f"{k}\t{s}\t{e}\t{(s+e)/2}\n")
PY

split -l "$SCAN_BATCH" -d -a 5 "$OUT/scan_windows.dat" "$OUT/scan_batches/chunk_"
: > "$OUT/scan_scores.tsv"
offset=0
for chunk in "$OUT"/scan_batches/chunk_*; do
    inp="$OUT/scan_batches/input.dat"; dst="$OUT/scan_batches/dist.dst"
    { cat "$OUT/ref_A.dat"; cat "$OUT/ref_B.dat"; cat "$chunk"; } > "$inp"
    "$LSMASH" -f "$inp" -D row -T symbolic -S 0 -o "$dst"
    python3 - "$dst" "$chunk" "$offset" >> "$OUT/scan_scores.tsv" <<'PY'
from pathlib import Path
import sys
dst,chunk,off=Path(sys.argv[1]),Path(sys.argv[2]),int(sys.argv[3])
M=[[float(v) for v in z.split()] for z in dst.read_text().splitlines() if z.strip()]
n=sum(bool(z.strip()) for z in chunk.read_text().splitlines())
if len(M)!=n+2: raise SystemExit(f"bad lsmash matrix: {len(M)} vs {n+2}")
for j in range(n):
    da,db=M[j+2][0],M[j+2][1]
    # positive g => A closer; negative g => B closer
    print(off+j,da,db,db-da,0 if da<=db else 1,sep="\t")
PY
    offset=$((offset + $(wc -l < "$chunk")))
done

# Infer change points from the continuous signed lsmash function.
python3 - "$OUT" "$BASE_LEN" "$SCAN_STEP" "$MIN_RUN_WINDOWS" <<'PY'
from pathlib import Path
import sys, statistics
out=Path(sys.argv[1]); L=int(sys.argv[2]); step=int(sys.argv[3]); minw=int(sys.argv[4])
meta=[]
for z in (out/"scan_meta.tsv").read_text().splitlines()[1:]:
    i,s,e,c=z.split("\t"); meta.append((int(i),int(s),int(e),float(c)))
raw=[]
g=[]
for z in (out/"scan_scores.tsv").read_text().splitlines():
    i,da,db,gg,lab=z.split("\t"); g.append(float(gg)); raw.append(int(lab))

# Median filter the continuous g(t), then use its sign. This preserves the
# two functions themselves in scan_scores.tsv while making change points robust.
gf=[]
for i in range(len(g)):
    gf.append(statistics.median(g[max(0,i-2):min(len(g),i+3)]))
lab=[0 if z>=0 else 1 for z in gf]

# Remove short label islands iteratively.
changed=True
while changed:
    changed=False; runs=[]; i=0
    while i<len(lab):
        j=i+1
        while j<len(lab) and lab[j]==lab[i]: j+=1
        runs.append([i,j,lab[i]])
        i=j
    for q,(a,b,s) in enumerate(runs):
        if b-a < minw and 0<q<len(runs)-1 and runs[q-1][2]==runs[q+1][2]:
            for k in range(a,b): lab[k]=runs[q-1][2]
            changed=True

runs=[]; i=0
while i<len(lab):
    j=i+1
    while j<len(lab) and lab[j]==lab[i]: j+=1
    runs.append((i,j,lab[i])); i=j

N=len((out/"transmission.dat").read_text().split())
centers=[m[3] for m in meta]
def boundary(k):
    if k<=0: return 0
    if k>=len(centers): return N
    return int(round((centers[k-1]+centers[k])/2))

with (out/"continuous_signal.tsv").open("w") as f:
    f.write("index\tcenter\td_A\td_B\tg\tg_median\tstate\n")
    scores=[z.split("\t") for z in (out/"scan_scores.tsv").read_text().splitlines()]
    for i,z in enumerate(scores):
        f.write(f"{i}\t{centers[i]}\t{z[1]}\t{z[2]}\t{z[3]}\t{gf[i]}\t{lab[i]}\n")

predbits=[]
with (out/"inferred_runs.tsv").open("w") as f:
    f.write("run\tbit\tstart\tend\tduration\tdecoded_count\tratio\n")
    for q,(a,b,s) in enumerate(runs):
        ss,ee=boundary(a),boundary(b)
        dur=ee-ss
        cnt=max(1,int(round(dur/L)))
        predbits.append(str(s)*cnt)
        f.write(f"{q}\t{s}\t{ss}\t{ee}\t{dur}\t{cnt}\t{dur/L:.6f}\n")
(out/"decoded_bits.txt").write_text("".join(predbits)+"\n")
PY

python3 - "$OUT" "$PLAINTEXT" <<'PY'
from pathlib import Path
import sys
out=Path(sys.argv[1]); text=sys.argv[2]
a=(out/"message_bits.txt").read_text().strip()
b=(out/"decoded_bits.txt").read_text().strip()
def lev(x,y):
    p=list(range(len(y)+1))
    for i,cx in enumerate(x,1):
        q=[i]
        for j,cy in enumerate(y,1):
            q.append(min(q[-1]+1,p[j]+1,p[j-1]+(cx!=cy)))
        p=q
    return p[-1]
e=lev(a,b)
decoded=""
if b and len(b)%8==0:
    bb=bytes(int(b[i:i+8],2) for i in range(0,len(b),8))
    decoded=bb.decode(errors="replace")
true_runs=max(0,len((out/"run_truth.tsv").read_text().splitlines())-1)
inf_runs=max(0,len((out/"inferred_runs.tsv").read_text().splitlines())-1)
s=[
 f"base_length={int(sys.argv[0] is None) if False else ''}",
 f"true_runs={true_runs}",
 f"inferred_runs={inf_runs}",
 f"true_bits={len(a)}",
 f"decoded_bits={len(b)}",
 f"bit_edit_distance={e}",
 f"normalized_bit_edit_distance={e/max(1,len(a)):.8f}",
 f"exact_bit_recovery={int(a==b)}",
 f"plaintext={text}",
 f"decoded_text={decoded}",
]
# remove placeholder first line
s=s[1:]
(out/"runlength_summary.txt").write_text("\n".join(s)+"\n")
print("\n".join(s))
PY

echo
echo "Artifacts written to $OUT"
echo "  continuous_signal.tsv  continuous d_A(t), d_B(t), g(t)"
echo "  inferred_runs.tsv      inferred A/B residence intervals and duration decoding"
echo "  decoded_bits.txt       reconstructed bit stream"
echo "  runlength_summary.txt  end-to-end metrics"
echo "  run_truth.tsv          hidden run truth, evaluation only"
