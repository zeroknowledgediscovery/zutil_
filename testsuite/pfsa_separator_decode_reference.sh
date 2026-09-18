#!/usr/bin/env bash
set -euo pipefail

# Blind variable-length PFSA decoder with a stochastic separator.
# Grammar: W, (A|B), W, (A|B), ... , W
# A=bit 0, B=bit 1, W=synchronizer.
# Decoder uses transmission.dat + independent A/B/W references only.
# truth files are evaluation-only.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

PRUN="${PRUN:-${ROOT}/bin/prun}"
LSMASH="${LSMASH:-${ROOT}/bin/lsmash}"
MODEL_A="${MODEL_A:-${ROOT}/testsuite/pfsa_opt/M2.cfg}"
MODEL_B="${MODEL_B:-${ROOT}/testsuite/pfsa_opt/T3.cfg}"
MODEL_W="${MODEL_W:-${ROOT}/testsuite/pfsa_opt/W.cfg}"

OUT="${OUT:-${ROOT}/testsuite/pfsa_separator_decode_run}"
DATA_MIN="${DATA_MIN:-500}"
DATA_MAX="${DATA_MAX:-1000}"
SEP_MIN="${SEP_MIN:-500}"
SEP_MAX="${SEP_MAX:-800}"
REF_LEN="${REF_LEN:-30000}"
SCAN_WIN="${SCAN_WIN:-250}"
SCAN_STEP="${SCAN_STEP:-100}"
SCAN_BATCH="${SCAN_BATCH:-20}"
SEG_BATCH="${SEG_BATCH:-12}"
MIN_W_WINDOWS="${MIN_W_WINDOWS:-2}"
MIN_DATA_KEEP="${MIN_DATA_KEEP:-180}"
SEED="${SEED:-1729}"
PLAINTEXT="${PLAINTEXT:-LSMASH STOCHASTIC CODE}"

for f in "$PRUN" "$LSMASH" "$MODEL_A" "$MODEL_B" "$MODEL_W"; do
    [[ -e "$f" ]] || { echo "missing required file: $f" >&2; exit 2; }
done
(( DATA_MIN > 0 && DATA_MAX >= DATA_MIN )) || { echo "bad DATA range" >&2; exit 2; }
(( SEP_MIN > 0 && SEP_MAX >= SEP_MIN )) || { echo "bad SEP range" >&2; exit 2; }
(( SCAN_WIN > 0 && SCAN_STEP > 0 )) || { echo "bad scan settings" >&2; exit 2; }

rm -rf "$OUT"
mkdir -p "$OUT/scan_batches" "$OUT/segment_batches"

printf 'A=%s\nB=%s\nW=%s\n' "$MODEL_A" "$MODEL_B" "$MODEL_W"
printf 'plaintext=%s\ndata lengths=[%s,%s]\nseparator lengths=[%s,%s]\n' \
       "$PLAINTEXT" "$DATA_MIN" "$DATA_MAX" "$SEP_MIN" "$SEP_MAX"
printf 'scan window=%s step=%s\n' "$SCAN_WIN" "$SCAN_STEP"

python3 - "$PLAINTEXT" "$DATA_MIN" "$DATA_MAX" "$SEP_MIN" "$SEP_MAX" "$SEED" "$OUT" <<'PY'
from pathlib import Path
import random, sys
text = sys.argv[1]
dlo, dhi = int(sys.argv[2]), int(sys.argv[3])
slo, shi = int(sys.argv[4]), int(sys.argv[5])
seed, out = int(sys.argv[6]), Path(sys.argv[7])
rng = random.Random(seed)
bits = "".join(f"{b:08b}" for b in text.encode("utf-8"))
(out / "message_bits.txt").write_text(bits + "\n")
seps = [rng.randint(slo, shi) for _ in range(len(bits) + 1)]
(out / "leading_separator_length.txt").write_text(str(seps[0]) + "\n")
with (out / "plan.tsv").open("w") as f:
    f.write("index\tbit\tdata_length\tseparator_after_length\n")
    for i, bit in enumerate(bits):
        f.write(f"{i}\t{bit}\t{rng.randint(dlo,dhi)}\t{seps[i+1]}\n")
PY

"$PRUN" -f "$MODEL_A" -l "$REF_LEN" -o "$OUT/ref_A.dat"
"$PRUN" -f "$MODEL_B" -l "$REF_LEN" -o "$OUT/ref_B.dat"
"$PRUN" -f "$MODEL_W" -l "$REF_LEN" -o "$OUT/ref_W.dat"

: > "$OUT/transmission.dat"
printf 'index\tbit\tmodel\tstart\tend\tlength\n' > "$OUT/truth.tsv"
printf 'index\tstart\tend\tlength\n' > "$OUT/separator_truth.tsv"

append_generated_block() {
    local model="$1"
    local len="$2"
    local tmp="$3"
    "$PRUN" -f "$model" -l "$len" -o "$tmp"
    tr '\n' ' ' < "$tmp" >> "$OUT/transmission.dat"
    rm -f "$tmp"
}

cursor=0
sep_idx=0
lead_sep=$(tr -d '[:space:]' < "$OUT/leading_separator_length.txt")
append_generated_block "$MODEL_W" "$lead_sep" "$OUT/.sep_${sep_idx}.dat"
printf '%s\t%s\t%s\t%s\n' "$sep_idx" "$cursor" "$((cursor+lead_sep))" "$lead_sep" >> "$OUT/separator_truth.tsv"
cursor=$((cursor + lead_sep))
sep_idx=$((sep_idx + 1))

while IFS=$'\t' read -r idx bit data_len sep_len; do
    [[ "$idx" == "index" ]] && continue
    if [[ "$bit" == "0" ]]; then model="$MODEL_A"; label="A"; else model="$MODEL_B"; label="B"; fi

    data_start=$cursor
    append_generated_block "$model" "$data_len" "$OUT/.data_${idx}.dat"
    cursor=$((cursor + data_len))
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$idx" "$bit" "$label" "$data_start" "$cursor" "$data_len" >> "$OUT/truth.tsv"

    sep_start=$cursor
    append_generated_block "$MODEL_W" "$sep_len" "$OUT/.sep_${sep_idx}.dat"
    cursor=$((cursor + sep_len))
    printf '%s\t%s\t%s\t%s\n' "$sep_idx" "$sep_start" "$cursor" "$sep_len" >> "$OUT/separator_truth.tsv"
    sep_idx=$((sep_idx + 1))
done < "$OUT/plan.tsv"
printf '\n' >> "$OUT/transmission.dat"

python3 - "$OUT/transmission.dat" "$SCAN_WIN" "$SCAN_STEP" "$OUT" <<'PY'
from pathlib import Path
import sys
path, win, step, out = Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), Path(sys.argv[4])
x = [int(v) for v in path.read_text().split()]
with (out/"scan_windows.dat").open("w") as fw, (out/"scan_meta.tsv").open("w") as fm:
    fm.write("index\tstart\tend\tcenter\n")
    j=0
    for start in range(0, len(x)-win+1, step):
        end=start+win
        fw.write(" ".join(map(str,x[start:end]))+"\n")
        fm.write(f"{j}\t{start}\t{end}\t{(start+end)/2.0}\n")
        j+=1
PY

split -l "$SCAN_BATCH" -d -a 5 "$OUT/scan_windows.dat" "$OUT/scan_batches/chunk_"
: > "$OUT/scan_scores.tsv"
scan_offset=0
for chunk in "$OUT"/scan_batches/chunk_*; do
    [[ -s "$chunk" ]] || continue
    inp="$OUT/scan_batches/input.dat"
    dst="$OUT/scan_batches/dist.dst"
    { cat "$OUT/ref_A.dat"; cat "$OUT/ref_B.dat"; cat "$OUT/ref_W.dat"; cat "$chunk"; } > "$inp"
    "$LSMASH" -f "$inp" -D row -T symbolic -S 0 -o "$dst"
    python3 - "$dst" "$chunk" "$scan_offset" >> "$OUT/scan_scores.tsv" <<'PY'
from pathlib import Path
import sys
dst, chunk, offset = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
M=[[float(v) for v in line.split()] for line in dst.read_text().splitlines() if line.strip()]
n=sum(bool(line.strip()) for line in chunk.read_text().splitlines())
if len(M)!=n+3 or any(len(r)!=n+3 for r in M):
    raise SystemExit(f"unexpected lsmash matrix shape {len(M)} for {n} windows")
for j in range(n):
    d=[M[j+3][0],M[j+3][1],M[j+3][2]]
    lab=min(range(3),key=d.__getitem__)
    print(offset+j,d[0],d[1],d[2],lab,sep="\t")
PY
    nlines=$(wc -l < "$chunk")
    scan_offset=$((scan_offset + nlines))
done

# Blind segmentation from W-classified windows only; no truth files read.
python3 - "$OUT" "$MIN_W_WINDOWS" "$MIN_DATA_KEEP" <<'PY'
from pathlib import Path
import sys
out=Path(sys.argv[1]); min_w=int(sys.argv[2]); min_data=int(sys.argv[3])

meta=[]
for line in (out/"scan_meta.tsv").read_text().splitlines()[1:]:
    i,s,e,c=line.split("\t"); meta.append((int(i),int(s),int(e),float(c)))
scores={}
for line in (out/"scan_scores.tsv").read_text().splitlines():
    if line.strip():
        i,da,db,dw,lab=line.split("\t")
        scores[int(i)]=(float(da),float(db),float(dw),int(lab))
labels=[scores[i][3] for i in range(len(meta))]

smooth=[]
for i,raw in enumerate(labels):
    vals=labels[max(0,i-1):min(len(labels),i+2)]
    counts=[vals.count(k) for k in range(3)]
    best=max(counts); winners=[k for k,c in enumerate(counts) if c==best]
    smooth.append(winners[0] if len(winners)==1 else raw)
for i in range(1,len(smooth)-1):
    if smooth[i]!=2 and smooth[i-1]==2 and smooth[i+1]==2:
        smooth[i]=2

runs=[]; i=0
while i<len(smooth):
    if smooth[i]!=2:
        i+=1; continue
    j=i
    while j+1<len(smooth) and smooth[j+1]==2: j+=1
    if j-i+1>=min_w:
        runs.append([meta[i][1],meta[j][2],i,j])
    i=j+1

merged=[]
for r in runs:
    if merged and r[0]<=merged[-1][1]:
        merged[-1][1]=max(merged[-1][1],r[1]); merged[-1][3]=r[3]
    else:
        merged.append(r[:])

with (out/"inferred_separators.tsv").open("w") as f:
    f.write("index\tstart\tend\twindow_i0\twindow_i1\n")
    for k,(s,e,i0,i1) in enumerate(merged):
        f.write(f"{k}\t{s}\t{e}\t{i0}\t{i1}\n")

x=[int(v) for v in (out/"transmission.dat").read_text().split()]
segments=[]
for k in range(len(merged)-1):
    s,e=merged[k][1],merged[k+1][0]
    if e-s>=min_data: segments.append((s,e))

with (out/"inferred_segments.dat").open("w") as fd, (out/"inferred_segments.tsv").open("w") as fm:
    fm.write("index\tstart\tend\tlength\n")
    for k,(s,e) in enumerate(segments):
        fd.write(" ".join(map(str,x[s:e]))+"\n")
        fm.write(f"{k}\t{s}\t{e}\t{e-s}\n")
PY

split -l "$SEG_BATCH" -d -a 5 "$OUT/inferred_segments.dat" "$OUT/segment_batches/chunk_" 2>/dev/null || true
: > "$OUT/segment_scores.tsv"
seg_offset=0
for chunk in "$OUT"/segment_batches/chunk_*; do
    [[ -s "$chunk" ]] || continue
    inp="$OUT/segment_batches/input.dat"
    dst="$OUT/segment_batches/dist.dst"
    { cat "$OUT/ref_A.dat"; cat "$OUT/ref_B.dat"; cat "$chunk"; } > "$inp"
    "$LSMASH" -f "$inp" -D row -T symbolic -S 0 -o "$dst"
    python3 - "$dst" "$chunk" "$seg_offset" >> "$OUT/segment_scores.tsv" <<'PY'
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
    seg_offset=$((seg_offset + nlines))
done

# Evaluation only.
python3 - "$OUT" "$PLAINTEXT" <<'PY'
from pathlib import Path
import sys
out=Path(sys.argv[1]); plaintext=sys.argv[2]
true_bits=(out/"message_bits.txt").read_text().strip()
pred_bits="".join(line.split("\t")[3] for line in (out/"segment_scores.tsv").read_text().splitlines() if line.strip())

def lev(a,b):
    p=list(range(len(b)+1))
    for i,ca in enumerate(a,1):
        c=[i]
        for j,cb in enumerate(b,1):
            c.append(min(c[-1]+1,p[j]+1,p[j-1]+(ca!=cb)))
        p=c
    return p[-1]
edit=lev(true_bits,pred_bits)

decoded_text=""; decoded_hex=""
if pred_bits and len(pred_bits)%8==0:
    bb=bytes(int(pred_bits[i:i+8],2) for i in range(0,len(pred_bits),8))
    decoded_hex=bb.hex(); decoded_text=bb.decode("utf-8",errors="replace")

true_sep=max(0,len((out/"separator_truth.tsv").read_text().splitlines())-1)
inf_sep=max(0,len((out/"inferred_separators.tsv").read_text().splitlines())-1)
summary=[
    f"true_bits={len(true_bits)}",
    f"decoded_bits={len(pred_bits)}",
    f"true_separators={true_sep}",
    f"inferred_separators={inf_sep}",
    f"bit_edit_distance={edit}",
    f"normalized_bit_edit_distance={edit/max(1,len(true_bits)):.8f}",
    f"exact_bit_recovery={int(true_bits==pred_bits)}",
    f"plaintext={plaintext}",
    f"decoded_text={decoded_text}",
    f"decoded_hex={decoded_hex}",
]
(out/"blind_summary.txt").write_text("\n".join(summary)+"\n")
print("\n".join(summary))
PY

printf '\nArtifacts written to %s\n' "$OUT"
printf '  transmission.dat          delimiter-free observed stream\n'
printf '  scan_scores.tsv           A/B/W lsmash window scores\n'
printf '  inferred_separators.tsv   blind W regions\n'
printf '  inferred_segments.tsv     blind candidate data blocks\n'
printf '  segment_scores.tsv        A/B lsmash decisions\n'
printf '  blind_summary.txt         end-to-end recovery metrics\n'
printf '  truth.tsv                 hidden data-block truth (evaluation only)\n'
printf '  separator_truth.tsv       hidden W-block truth (evaluation only)\n'
