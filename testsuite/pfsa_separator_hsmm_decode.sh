#!/usr/bin/env bash
set -euo pipefail

# Global semi-Markov decoder for an existing pfsa_separator_decode_reference run.
# Uses only transmission.dat, scan_meta.tsv, scan_scores.tsv, and A/B references.
# Hidden truth files are read only in the final evaluation step.

RUN="$1"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LSMASH="$ROOT/bin/lsmash"

DATA_MIN="${DATA_MIN:-500}"
DATA_MAX="${DATA_MAX:-1000}"
SEP_MIN="${SEP_MIN:-500}"
SEP_MAX="${SEP_MAX:-800}"
DUR_TOL_WINDOWS="${DUR_TOL_WINDOWS:-2}"
SEG_BATCH="${SEG_BATCH:-12}"

for f in transmission.dat scan_meta.tsv scan_scores.tsv ref_A.dat ref_B.dat message_bits.txt; do
    [[ -s "$RUN/$f" ]] || { echo "missing $RUN/$f" >&2; exit 2; }
done

mkdir -p "$RUN/hsmm_batches"

python3 - "$RUN" "$DATA_MIN" "$DATA_MAX" "$SEP_MIN" "$SEP_MAX" "$DUR_TOL_WINDOWS" <<'PY'
from pathlib import Path
import math, sys

run = Path(sys.argv[1])
data_min, data_max = int(sys.argv[2]), int(sys.argv[3])
sep_min, sep_max = int(sys.argv[4]), int(sys.argv[5])
tol = int(sys.argv[6])

meta = []
for line in (run/"scan_meta.tsv").read_text().splitlines()[1:]:
    i,s,e,c = line.split("\t")
    meta.append((int(i), int(s), int(e), float(c)))

scores = []
for line in (run/"scan_scores.tsv").read_text().splitlines():
    if line.strip():
        i,da,db,dw,lab = line.split("\t")
        scores.append((int(i), float(da), float(db), float(dw)))

if not meta or len(meta) != len(scores):
    raise SystemExit("scan metadata/score mismatch")

win = meta[0][2] - meta[0][1]
step = int(round(meta[1][3] - meta[0][3])) if len(meta) > 1 else win
n = len(meta)

# Relative per-window evidence: the best source gets zero cost; alternatives
# pay only their excess lsmash distance. This makes the DP depend on source
# preference rather than arbitrary absolute distance scale.
cost = [[0.0]*n for _ in range(3)]  # 0=A,1=B,2=W
for i,(_,da,db,dw) in enumerate(scores):
    ds = [da,db,dw]
    m = min(ds)
    scale = max(1e-9, sorted(ds)[1] - m)
    # Cap extreme windows so one noisy scan cannot dominate an entire segment.
    for k in range(3):
        cost[k][i] = min(8.0, max(0.0, (ds[k]-m)/scale))

prefix = [[0.0]*(n+1) for _ in range(3)]
for k in range(3):
    for i in range(n):
        prefix[k][i+1] = prefix[k][i] + cost[k][i]

def segcost(k, a, b):
    return prefix[k][b] - prefix[k][a]

# A block of L symbols contributes roughly L/step window centers. Mixed
# boundary windows motivate a small tolerance.
def dur_range(lo, hi):
    a = max(1, int(math.floor(lo/step)) - tol)
    b = max(a, int(math.ceil(hi/step)) + tol)
    return range(a, b+1)

DR = list(dur_range(data_min, data_max))
WR = list(dur_range(sep_min, sep_max))

INF = 1e100
# dp[t][state] = minimum cost for first t scan windows, ending with a segment
# of 'state'. Valid transitions are W->A/B and A/B->W. Stream starts with W.
dp = [[INF]*3 for _ in range(n+1)]
back = [[None]*3 for _ in range(n+1)]

for d in WR:
    if d <= n:
        dp[d][2] = segcost(2,0,d)
        back[d][2] = (0,-1,d)

for t in range(1,n+1):
    # data segments after W
    for state in (0,1):
        for d in DR:
            a=t-d
            if a < 0: continue
            v=dp[a][2]
            if v >= INF: continue
            z=v+segcost(state,a,t)
            if z < dp[t][state]:
                dp[t][state]=z
                back[t][state]=(a,2,d)
    # W segment after A or B
    for d in WR:
        a=t-d
        if a < 0: continue
        for prev in (0,1):
            v=dp[a][prev]
            if v >= INF: continue
            z=v+segcost(2,a,t)
            if z < dp[t][2]:
                dp[t][2]=z
                back[t][2]=(a,prev,d)

# Allow a short unscanned tail after the final complete W segment.
best_t = None
best_obj = INF
tail_allow = max(WR) + 2
for t in range(max(1,n-tail_allow), n+1):
    if dp[t][2] < INF:
        # tiny penalty prefers explaining more windows when scores tie
        z = dp[t][2] + 0.02*(n-t)
        if z < best_obj:
            best_obj=z; best_t=t

if best_t is None:
    raise SystemExit("no duration-constrained W/(A|B)/W path found")

segments=[]
t=best_t; state=2
while t>0:
    rec=back[t][state]
    if rec is None:
        raise SystemExit(f"broken backpointer at t={t}, state={state}")
    a,prev,d=rec
    segments.append((a,t,state))
    if prev == -1:
        break
    t,state=a,prev
segments.reverse()

# Convert scan-window segment boundaries to symbol coordinates using midpoint
# between neighboring scan-window centers.
centers=[r[3] for r in meta]
N=len((run/"transmission.dat").read_text().split())

def symbol_boundary(window_index):
    if window_index <= 0: return 0
    if window_index >= len(centers): return N
    return int(round((centers[window_index-1]+centers[window_index])/2.0))

with (run/"hsmm_state_path.tsv").open("w") as f:
    f.write("segment\tstate\twindow_start\twindow_end\tsymbol_start\tsymbol_end\twindow_count\n")
    for j,(a,b,s) in enumerate(segments):
        lab="ABW"[s]
        ss=symbol_boundary(a); ee=symbol_boundary(b)
        f.write(f"{j}\t{lab}\t{a}\t{b}\t{ss}\t{ee}\t{b-a}\n")

x=[int(v) for v in (run/"transmission.dat").read_text().split()]
data_segments=[]
for a,b,s in segments:
    if s in (0,1):
        ss=symbol_boundary(a); ee=symbol_boundary(b)
        if ee>ss:
            data_segments.append((ss,ee,s))

with (run/"hsmm_segments.dat").open("w") as fd, (run/"hsmm_segments.tsv").open("w") as fm:
    fm.write("index\tstart\tend\tlength\twindow_state\n")
    for j,(s,e,st) in enumerate(data_segments):
        fd.write(" ".join(map(str,x[s:e]))+"\n")
        fm.write(f"{j}\t{s}\t{e}\t{e-s}\t{'A' if st==0 else 'B'}\n")

print(f"scan_windows={n}")
print(f"scan_step={step}")
print(f"data_window_duration_range={min(DR)}..{max(DR)}")
print(f"W_window_duration_range={min(WR)}..{max(WR)}")
print(f"hsmm_total_segments={len(segments)}")
print(f"hsmm_data_segments={len(data_segments)}")
print(f"hsmm_W_segments={sum(s==2 for _,_,s in segments)}")
PY

rm -rf "$RUN/hsmm_batches"
mkdir -p "$RUN/hsmm_batches"
split -l "$SEG_BATCH" -d -a 5 "$RUN/hsmm_segments.dat" "$RUN/hsmm_batches/chunk_" 2>/dev/null || true
: > "$RUN/hsmm_segment_scores.tsv"
offset=0
for chunk in "$RUN"/hsmm_batches/chunk_*; do
    [[ -s "$chunk" ]] || continue
    inp="$RUN/hsmm_batches/input.dat"
    dst="$RUN/hsmm_batches/dist.dst"
    { cat "$RUN/ref_A.dat"; cat "$RUN/ref_B.dat"; cat "$chunk"; } > "$inp"
    "$LSMASH" -f "$inp" -D row -T symbolic -S 0 -o "$dst"
    python3 - "$dst" "$chunk" "$offset" >> "$RUN/hsmm_segment_scores.tsv" <<'PY'
from pathlib import Path
import sys
dst,chunk,offset=Path(sys.argv[1]),Path(sys.argv[2]),int(sys.argv[3])
M=[[float(v) for v in line.split()] for line in dst.read_text().splitlines() if line.strip()]
n=sum(bool(line.strip()) for line in chunk.read_text().splitlines())
if len(M)!=n+2 or any(len(r)!=n+2 for r in M):
    raise SystemExit(f"unexpected lsmash matrix shape {len(M)} for {n} segments")
for j in range(n):
    da,db=M[j+2][0],M[j+2][1]
    print(offset+j,da,db,0 if da<=db else 1,sep="\t")
PY
    nlines=$(wc -l < "$chunk")
    offset=$((offset+nlines))
done

# Evaluation only. No truth was used above.
python3 - "$RUN" <<'PY'
from pathlib import Path
import sys
run=Path(sys.argv[1])
true_bits=(run/"message_bits.txt").read_text().strip()
pred_bits="".join(line.split("\t")[3] for line in (run/"hsmm_segment_scores.tsv").read_text().splitlines() if line.strip())

def lev(a,b):
    p=list(range(len(b)+1))
    for i,ca in enumerate(a,1):
        c=[i]
        for j,cb in enumerate(b,1):
            c.append(min(c[-1]+1,p[j]+1,p[j-1]+(ca!=cb)))
        p=c
    return p[-1]
edit=lev(true_bits,pred_bits)

decoded=""
if pred_bits and len(pred_bits)%8==0:
    bb=bytes(int(pred_bits[i:i+8],2) for i in range(0,len(pred_bits),8))
    decoded=bb.decode("utf-8",errors="replace")

summary=[
    f"true_bits={len(true_bits)}",
    f"decoded_bits={len(pred_bits)}",
    f"bit_edit_distance={edit}",
    f"normalized_bit_edit_distance={edit/max(1,len(true_bits)):.8f}",
    f"exact_bit_recovery={int(true_bits==pred_bits)}",
    f"decoded_text={decoded}",
]
(run/"hsmm_summary.txt").write_text("\n".join(summary)+"\n")
print("\n".join(summary))
PY

echo
echo "HSMM artifacts:"
echo "  $RUN/hsmm_state_path.tsv"
echo "  $RUN/hsmm_segments.tsv"
echo "  $RUN/hsmm_segment_scores.tsv"
echo "  $RUN/hsmm_summary.txt"
