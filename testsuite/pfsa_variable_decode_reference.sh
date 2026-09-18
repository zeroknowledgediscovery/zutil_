#!/usr/bin/env bash
set -euo pipefail

# Reference experiment for variable-length stochastic codewords.
#
# A bit 0 is encoded by a fresh realization from PFSA A (M2.cfg).
# A bit 1 is encoded by a fresh realization from PFSA B (T3.cfg).
# Codeword lengths are sampled uniformly from [MIN_BLOCK, MAX_BLOCK].
# The emitted codewords are concatenated into transmission.dat with no
# delimiters.  truth.tsv records the hidden segmentation for the oracle-
# boundary reference decoder.
#
# lsmash is then used to classify each true block against long reference
# streams from A and B.  This establishes the process-identification baseline
# before attempting unknown-boundary decoding.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

PRUN="${PRUN:-${ROOT}/bin/prun}"
LSMASH="${LSMASH:-${ROOT}/bin/lsmash}"
MODEL_A="${MODEL_A:-${ROOT}/testsuite/pfsa_opt/M2.cfg}"
MODEL_B="${MODEL_B:-${ROOT}/testsuite/pfsa_opt/T3.cfg}"

OUT="${OUT:-${ROOT}/testsuite/pfsa_variable_decode_run}"
MIN_BLOCK="${MIN_BLOCK:-1200}"
MAX_BLOCK="${MAX_BLOCK:-3000}"
REF_LEN="${REF_LEN:-30000}"
BATCH_SIZE="${BATCH_SIZE:-12}"
SEED="${SEED:-1729}"
PLAINTEXT="${PLAINTEXT:-LSMASH STOCHASTIC CODE}"

for f in "$PRUN" "$LSMASH" "$MODEL_A" "$MODEL_B"; do
    if [[ ! -e "$f" ]]; then
        echo "missing required file: $f" >&2
        exit 2
    fi
done
if (( MIN_BLOCK <= 0 || MAX_BLOCK < MIN_BLOCK )); then
    echo "invalid block-length interval: [$MIN_BLOCK,$MAX_BLOCK]" >&2
    exit 2
fi

rm -rf "$OUT"
mkdir -p "$OUT/batches"

printf 'A model: %s\n' "$MODEL_A"
printf 'B model: %s\n' "$MODEL_B"
printf 'plaintext: %s\n' "$PLAINTEXT"
printf 'block lengths: [%s,%s]\n' "$MIN_BLOCK" "$MAX_BLOCK"
printf 'reference length: %s\n' "$REF_LEN"

# Encode the UTF-8 payload into bits and choose reproducible variable lengths.
python3 - "$PLAINTEXT" "$MIN_BLOCK" "$MAX_BLOCK" "$SEED" "$OUT" <<'PY'
import pathlib, random, sys

text, lo, hi, seed, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), pathlib.Path(sys.argv[5])
rng = random.Random(seed)
payload = text.encode("utf-8")
bits = "".join(f"{b:08b}" for b in payload)
(out / "message_bits.txt").write_text(bits + "\n")
with (out / "plan.tsv").open("w") as f:
    for i, bit in enumerate(bits):
        f.write(f"{i}\t{bit}\t{rng.randint(lo, hi)}\n")
PY

# Long independent source references used by lsmash.
"$PRUN" -f "$MODEL_A" -l "$REF_LEN" -o "$OUT/ref_A.dat"
"$PRUN" -f "$MODEL_B" -l "$REF_LEN" -o "$OUT/ref_B.dat"

: > "$OUT/blocks.dat"
: > "$OUT/transmission.dat"
printf 'index\tbit\tmodel\tlength\tstart\tend\n' > "$OUT/truth.tsv"

start=0
while IFS=$'\t' read -r idx bit len; do
    if [[ "$bit" == "0" ]]; then
        model="$MODEL_A"
        label="A"
    else
        model="$MODEL_B"
        label="B"
    fi

    tmp="$OUT/.block_${idx}.dat"
    "$PRUN" -f "$model" -l "$len" -o "$tmp"

    # One codeword per row for the oracle-boundary lsmash baseline.
    cat "$tmp" >> "$OUT/blocks.dat"

    # Also construct the actual delimiter-free transmitted stream.
    tr '\n' ' ' < "$tmp" >> "$OUT/transmission.dat"

    end=$((start + len))
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$idx" "$bit" "$label" "$len" "$start" "$end" >> "$OUT/truth.tsv"
    start="$end"
    rm -f "$tmp"
done < "$OUT/plan.tsv"
printf '\n' >> "$OUT/transmission.dat"

# lsmash computes a full pairwise matrix, so classify blocks in small batches.
# In each batch, rows 0 and 1 are the A/B references and rows >=2 are codewords.
split -l "$BATCH_SIZE" -d -a 4 "$OUT/blocks.dat" "$OUT/batches/chunk_"
: > "$OUT/scores.tsv"
offset=0
for chunk in "$OUT"/batches/chunk_*; do
    [[ -s "$chunk" ]] || continue
    batch_dat="$OUT/batches/input.dat"
    batch_dst="$OUT/batches/dist.dst"
    {
        cat "$OUT/ref_A.dat"
        cat "$OUT/ref_B.dat"
        cat "$chunk"
    } > "$batch_dat"

    "$LSMASH" -f "$batch_dat" -D row -T symbolic -S 0 -o "$batch_dst"

    python3 - "$batch_dst" "$chunk" "$offset" >> "$OUT/scores.tsv" <<'PY'
import pathlib, sys

dst, chunk, offset = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), int(sys.argv[3])
M = [[float(x) for x in line.split()] for line in dst.read_text().splitlines() if line.strip()]
n = sum(1 for line in chunk.read_text().splitlines() if line.strip())
if len(M) != n + 2 or any(len(row) != n + 2 for row in M):
    raise SystemExit(f"unexpected lsmash matrix shape: {len(M)} for {n} blocks")
for j in range(n):
    dA = M[j + 2][0]
    dB = M[j + 2][1]
    pred = 0 if dA <= dB else 1
    print(offset + j, dA, dB, pred, sep="\t")
PY
    nlines=$(wc -l < "$chunk")
    offset=$((offset + nlines))
done

# Summarize lsmash decoding and compare against a deliberately weak baseline
# that uses only each block's fraction of emitted 1s.
python3 - "$OUT" "$PLAINTEXT" "$MODEL_A" "$MODEL_B" <<'PY'
from pathlib import Path
import math, sys

out = Path(sys.argv[1])
plaintext = sys.argv[2]
model_a = sys.argv[3]
model_b = sys.argv[4]

def tokens(path):
    return [int(x) for x in Path(path).read_text().split()]

refA = tokens(out / "ref_A.dat")
refB = tokens(out / "ref_B.dat")
pA = sum(refA) / len(refA)
pB = sum(refB) / len(refB)

truth = []
for line in (out / "truth.tsv").read_text().splitlines()[1:]:
    idx, bit, model, length, start, end = line.split("\t")
    truth.append((int(idx), int(bit), model, int(length), int(start), int(end)))

scores = {}
for line in (out / "scores.tsv").read_text().splitlines():
    if not line.strip():
        continue
    idx, dA, dB, pred = line.split("\t")
    scores[int(idx)] = (float(dA), float(dB), int(pred))

block_lines = [line for line in (out / "blocks.dat").read_text().splitlines() if line.strip()]
if len(block_lines) != len(truth):
    raise SystemExit("block/truth count mismatch")

rows = []
lsmash_bits = []
freq_bits = []
for rec, line in zip(truth, block_lines):
    idx, true_bit, model, length, start, end = rec
    dA, dB, pred = scores[idx]
    vals = [int(x) for x in line.split()]
    p1 = sum(vals) / len(vals)
    f_pred = 0 if abs(p1 - pA) <= abs(p1 - pB) else 1
    lsmash_bits.append(pred)
    freq_bits.append(f_pred)
    rows.append((idx, true_bit, pred, f_pred, length, dA, dB, abs(dA-dB), p1, int(pred == true_bit)))

true_bits = [r[1] for r in truth]
lsmash_err = sum(a != b for a, b in zip(true_bits, lsmash_bits))
freq_err = sum(a != b for a, b in zip(true_bits, freq_bits))

def bits_to_bytes(bits):
    return bytes(int("".join(str(x) for x in bits[i:i+8]), 2) for i in range(0, len(bits), 8))

decoded_bytes = bits_to_bytes(lsmash_bits)
try:
    decoded_text = decoded_bytes.decode("utf-8")
except UnicodeDecodeError:
    decoded_text = decoded_bytes.decode("utf-8", errors="replace")

with (out / "decode.tsv").open("w") as f:
    f.write("index\ttrue_bit\tlsmash_bit\tfreq_bit\tlength\td_A\td_B\tmargin\tp1\tcorrect\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")

summary = [
    f"model_A={model_a}",
    f"model_B={model_b}",
    f"reference_A_p1={pA:.8f}",
    f"reference_B_p1={pB:.8f}",
    f"message_bits={len(true_bits)}",
    f"transmitted_symbols={truth[-1][5] if truth else 0}",
    f"lsmash_bit_errors={lsmash_err}",
    f"lsmash_bit_accuracy={1-lsmash_err/len(true_bits):.8f}",
    f"frequency_bit_errors={freq_err}",
    f"frequency_bit_accuracy={1-freq_err/len(true_bits):.8f}",
    f"exact_recovery={int(lsmash_err == 0)}",
    f"plaintext={plaintext}",
    f"decoded_text={decoded_text}",
    f"decoded_hex={decoded_bytes.hex()}",
]
(out / "summary.txt").write_text("\n".join(summary) + "\n")
print("\n".join(summary))
PY

printf '\nArtifacts written to %s\n' "$OUT"
printf '  transmission.dat  delimiter-free concatenated stream\n'
printf '  truth.tsv         hidden source/boundary sequence\n'
printf '  blocks.dat        true codeword segmentation (oracle baseline only)\n'
printf '  decode.tsv        per-codeword lsmash scores and decisions\n'
printf '  summary.txt       aggregate decoding results\n'
