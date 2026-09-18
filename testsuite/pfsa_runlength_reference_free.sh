#!/usr/bin/env bash
set -euo pipefail

# Reference-free stochastic-symbol reconstruction with LSmash.
#
# A delimiter-free transmission is generated from two stochastic sources.
# The decoder is given ONLY transmission.dat. It does not see source models,
# A/B reference realizations, source order, or true boundaries.
#
# It cuts the transmission into overlapping windows, computes the all-pairs
# LSmash distance matrix between those windows, finds two latent source
# classes from the intrinsic distance geometry, restores temporal continuity,
# and decodes residence duration into repeated bits.
#
# Since there is no external source label, recovery is identifiable only up
# to a global 0<->1 inversion. Both orientations are reported.
#
# Default generator uses the matched-first-order order-2 pair:
# A: P(1|00,01,10,11) = (.9,.1,.1,.9)
# B: P(1|00,01,10,11) = (.1,.9,.9,.1)
# Both have uniform adjacent-pair stationary distribution.
#
# Run from zutil_:
#   testsuite/pfsa_runlength_reference_free.sh
#
# Useful harder/easier settings:
#   BASE_LEN=2000 SCAN_WIN=1000 SCAN_STEP=200 ./testsuite/pfsa_runlength_reference_free.sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LSMASH="${LSMASH:-$ROOT/bin/lsmash}"
OUT="${OUT:-$ROOT/testsuite/pfsa_runlength_reference_free_run}"
BASE_LEN="${BASE_LEN:-1000}"
SCAN_WIN="${SCAN_WIN:-500}"
SCAN_STEP="${SCAN_STEP:-100}"
MIN_RUN_WINDOWS="${MIN_RUN_WINDOWS:-2}"
SEED="${SEED:-20260918}"
PLAINTEXT="${PLAINTEXT:-LSMASH STOCHASTIC CODE}"

rm -rf "$OUT"
mkdir -p "$OUT"

# Generate hidden transmission. Truth is written only for evaluation after
# unsupervised decoding; it is never used by the decoder.
python3 - "$PLAINTEXT" "$OUT" "$BASE_LEN" "$SEED" <<'PY'
from pathlib import Path
import random,sys
text,out,L,seed=sys.argv[1],Path(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
rng=random.Random(seed)
QA=(.9,.1,.1,.9)
QB=(.1,.9,.9,.1)

def gen(q,n):
    a=rng.randrange(2); b=rng.randrange(2)
    z=[a,b][:n]
    while len(z)<n:
        ctx=2*z[-2]+z[-1]
        z.append(1 if rng.random()<q[ctx] else 0)
    return z

bits="".join(f"{b:08b}" for b in text.encode())
(out/"message_bits.txt").write_text(bits+"\n")
runs=[]; i=0
while i<len(bits):
    j=i+1
    while j<len(bits) and bits[j]==bits[i]: j+=1
    runs.append((int(bits[i]),j-i)); i=j

tx=[]; cursor=0
with (out/"run_truth.tsv").open("w") as f:
    f.write("run\tbit\tcount\tstart\tend\tlength\n")
    for r,(bit,count) in enumerate(runs):
        n=count*L
        z=gen(QA if bit==0 else QB,n)
        tx.extend(z)
        f.write(f"{r}\t{bit}\t{count}\t{cursor}\t{cursor+n}\t{n}\n")
        cursor+=n
(out/"transmission.dat").write_text(" ".join(map(str,tx))+"\n")
PY

# Decoder starts here. Only transmission.dat is read.
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

echo "Computing all-pairs LSmash geometry from transmission windows only..."
"$LSMASH" -f "$OUT/scan_windows.dat" -D row -T symbolic -S 0 -o "$OUT/window_distances.dst"

python3 - "$OUT" "$BASE_LEN" "$MIN_RUN_WINDOWS" <<'PY'
from pathlib import Path
import sys,math,statistics
out=Path(sys.argv[1]); L=int(sys.argv[2]); minw=int(sys.argv[3])

meta=[]
for z in (out/"scan_meta.tsv").read_text().splitlines()[1:]:
    i,s,e,c=z.split("\t"); meta.append((int(i),int(s),int(e),float(c)))
centers=[x[3] for x in meta]
n=len(meta)

D=[[float(v) for v in z.split()] for z in (out/"window_distances.dst").read_text().splitlines() if z.strip()]
if len(D)!=n or any(len(r)!=n for r in D):
    raise SystemExit(f"Expected {n}x{n} LSmash matrix, got {len(D)} rows")

# Classical MDS first coordinate. For a two-source mixture, the dominant
# contrast in the distance geometry should separate the two latent processes.
# Implement power iteration on B=-1/2 J D^2 J, avoiding numpy/scipy.
rowmean=[sum(v*v for v in r)/n for r in D]
grand=sum(rowmean)/n
def Bmul(x):
    # B_ij=-.5*(dij^2-rowmean_i-rowmean_j+grand)
    sx=sum(x); rmx=sum(rowmean[j]*x[j] for j in range(n))
    y=[]
    for i in range(n):
        dij=sum((D[i][j]**2)*x[j] for j in range(n))
        y.append(-.5*(dij-rowmean[i]*sx-rmx+grand*sx))
    return y

# deterministic nonconstant initialization
x=[math.sin((i+1)*1.61803398875) for i in range(n)]
mu=sum(x)/n; x=[v-mu for v in x]
for _ in range(80):
    y=Bmul(x)
    mu=sum(y)/n; y=[v-mu for v in y]
    norm=math.sqrt(sum(v*v for v in y))
    if norm==0: raise SystemExit("Degenerate LSmash geometry")
    x=[v/norm for v in y]

# 1D two-means on intrinsic coordinate; labels remain anonymous.
c0=min(x); c1=max(x)
for _ in range(100):
    lab=[0 if abs(v-c0)<=abs(v-c1) else 1 for v in x]
    a=[v for v,l in zip(x,lab) if l==0]; b=[v for v,l in zip(x,lab) if l==1]
    if not a or not b: raise SystemExit("Two-cluster split collapsed")
    q0,q1=sum(a)/len(a),sum(b)/len(b)
    if abs(q0-c0)+abs(q1-c1)<1e-12: break
    c0,c1=q0,q1

# Median smooth the continuous coordinate, then classify by learned midpoint.
xf=[]
for i in range(n):
    xf.append(statistics.median(x[max(0,i-2):min(n,i+3)]))
threshold=(c0+c1)/2
lab=[0 if v<=threshold else 1 for v in xf]

# Remove short enclosed islands.
changed=True
while changed:
    changed=False; runs=[]; i=0
    while i<n:
        j=i+1
        while j<n and lab[j]==lab[i]: j+=1
        runs.append([i,j,lab[i]]); i=j
    for q,(a,b,s) in enumerate(runs):
        if b-a<minw and 0<q<len(runs)-1 and runs[q-1][2]==runs[q+1][2]:
            for k in range(a,b): lab[k]=runs[q-1][2]
            changed=True; break

runs=[]; i=0
while i<n:
    j=i+1
    while j<n and lab[j]==lab[i]: j+=1
    runs.append((i,j,lab[i])); i=j

N=len((out/"transmission.dat").read_text().split())
def boundary(k):
    if k<=0:return 0
    if k>=n:return N
    return int(round((centers[k-1]+centers[k])/2))

decoded=[]
with (out/"latent_runs.tsv").open("w") as f:
    f.write("run\tlatent_state\tstart\tend\tduration\tdecoded_count\n")
    for q,(a,b,s) in enumerate(runs):
        ss,ee=boundary(a),boundary(b); dur=ee-ss
        cnt=max(1,round(dur/L))
        decoded.append(str(s)*cnt)
        f.write(f"{q}\t{s}\t{ss}\t{ee}\t{dur}\t{cnt}\n")

z="".join(decoded)
zinv="".join("1" if c=="0" else "0" for c in z)
(out/"decoded_orientation_0.txt").write_text(z+"\n")
(out/"decoded_orientation_1.txt").write_text(zinv+"\n")

def text(bits):
    if not bits or len(bits)%8:return ""
    raw=bytes(int(bits[i:i+8],2) for i in range(0,len(bits),8))
    try:return raw.decode("utf-8")
    except:return raw.decode("utf-8",errors="replace")

t0,t1=text(z),text(zinv)
(out/"decoded_text_orientation_0.txt").write_text(t0+"\n")
(out/"decoded_text_orientation_1.txt").write_text(t1+"\n")
with (out/"latent_signal.tsv").open("w") as f:
    f.write("index\tcenter\tmds_coordinate\tmedian_coordinate\tlatent_state\n")
    for i in range(n):
        f.write(f"{i}\t{centers[i]}\t{x[i]:.12g}\t{xf[i]:.12g}\t{lab[i]}\n")
(out/"reference_free_summary.txt").write_text(
    f"inferred_runs={len(runs)}\n"
    f"decoded_bits={len(z)}\n"
    f"orientation_0_text={t0}\n"
    f"orientation_1_text={t1}\n"
)
print((out/"reference_free_summary.txt").read_text(),end="")
PY

# Evaluation only. Neither message_bits.txt nor plaintext enters the decoder.
python3 - "$OUT" <<'PY'
from pathlib import Path
import sys
out=Path(sys.argv[1])
truth=(out/"message_bits.txt").read_text().strip()
a=(out/"decoded_orientation_0.txt").read_text().strip()
b=(out/"decoded_orientation_1.txt").read_text().strip()
def lev(x,y):
    p=list(range(len(y)+1))
    for i,c in enumerate(x,1):
        q=[i]
        for j,d in enumerate(y,1):
            q.append(min(q[-1]+1,p[j]+1,p[j-1]+(c!=d)))
        p=q
    return p[-1]
ea,eb=lev(truth,a),lev(truth,b)
best=min(ea,eb)
s=(
 f"truth_bits={len(truth)}\n"
 f"orientation_0_edit={ea}\n"
 f"orientation_1_edit={eb}\n"
 f"best_global_orientation_edit={best}\n"
 f"normalized_best_edit={best/max(1,len(truth)):.8f}\n"
 f"exact_up_to_global_inversion={int(best==0)}\n"
)
(out/"evaluation.txt").write_text(s)
print(s,end="")
PY

echo
echo "Artifacts written to $OUT"
echo "  transmission.dat             only observation used by decoder"
echo "  window_distances.dst         all-pairs LSmash geometry"
echo "  latent_signal.tsv            unsupervised two-state coordinate"
echo "  latent_runs.tsv              inferred residence intervals"
echo "  decoded_orientation_0.txt    anonymous state orientation"
echo "  decoded_orientation_1.txt    globally inverted orientation"
echo "  reference_free_summary.txt   unsupervised reconstruction"
echo "  evaluation.txt               truth comparison, evaluation only"
