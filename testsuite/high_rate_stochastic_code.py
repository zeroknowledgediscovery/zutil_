#!/usr/bin/env python3
"""
High-rate stochastic-code demonstration.

Each 120-bit information block is encoded as an extended? No: standard
binary Hamming(127,120) codeword.  The stochastic encoder then independently
flips every emitted bit with probability p.  Thus a message selects a
probability distribution ("stochastic cloud") over 127-bit realizations.

Nominal information rate = 120/127 = 0.9448818898 bits/emitted symbol.

The decoder sees only the stochastic realization.  It does not receive the
noise realization.  Hamming decoding corrects any single flip per block.

This is a finite-block demonstration, not a proof that Hamming(127,120)
achieves BSC capacity.  For sufficiently small p it gives a simple executable
example with rate > 0.90 and very high exact-recovery probability.
"""
import argparse, json, math, random

N = 127
RPAR = 7
K = N - RPAR
PARITY_POS = {1 << i for i in range(RPAR)}
DATA_POS = [i for i in range(1, N + 1) if i not in PARITY_POS]

def encode_block(data):
    if len(data) != K:
        raise ValueError(f"need {K} data bits")
    cw = [0] * (N + 1)  # positions 1..127
    for pos, bit in zip(DATA_POS, data):
        cw[pos] = int(bit)
    for ppos in sorted(PARITY_POS):
        parity = 0
        for pos in range(1, N + 1):
            if pos & ppos:
                parity ^= cw[pos]
        cw[ppos] = parity
    return cw[1:]

def syndrome(cw0):
    cw = [0] + list(map(int, cw0))
    s = 0
    for ppos in sorted(PARITY_POS):
        parity = 0
        for pos in range(1, N + 1):
            if pos & ppos:
                parity ^= cw[pos]
        if parity:
            s |= ppos
    return s

def decode_block(rx):
    cw = list(map(int, rx))
    s = syndrome(cw)
    corrected = False
    if 1 <= s <= N:
        cw[s - 1] ^= 1
        corrected = True
    data = [cw[pos - 1] for pos in DATA_POS]
    return data, s, corrected

def stochastic_channel(cw, p, rng):
    return [b ^ (1 if rng.random() < p else 0) for b in cw]

def bits_to_str(x):
    return "".join(map(str, x))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", type=int, default=1000)
    ap.add_argument("--p", type=float, default=1e-4)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--same-message-trials", type=int, default=4)
    ap.add_argument("--json", default="")
    args = ap.parse_args()
    if not (0 <= args.p < 0.5):
        raise SystemExit("p must satisfy 0 <= p < 0.5")
    rng = random.Random(args.seed)

    total_info = total_emitted = total_flips = 0
    block_errors = bit_errors = corrected_blocks = 0
    for _ in range(args.blocks):
        msg = [rng.getrandbits(1) for _ in range(K)]
        cw = encode_block(msg)
        rx = stochastic_channel(cw, args.p, rng)
        dec, s, corrected = decode_block(rx)
        flips = sum(a != b for a, b in zip(cw, rx))
        total_flips += flips
        corrected_blocks += int(corrected)
        errs = sum(a != b for a, b in zip(msg, dec))
        bit_errors += errs
        block_errors += int(errs > 0)
        total_info += K
        total_emitted += N

    # Demonstrate genuine stochasticity: encode the same fixed message repeatedly.
    fixed = [((i * 37 + 11) >> 2) & 1 for i in range(K)]
    fixed_cw = encode_block(fixed)
    trials = []
    for j in range(args.same_message_trials):
        rr = random.Random(args.seed + 1000003 + j)
        trials.append(bits_to_str(stochastic_channel(fixed_cw, args.p, rr)))
    unique_same_message_outputs = len(set(trials))

    H2 = 0.0 if args.p in (0.0, 1.0) else (
        -args.p * math.log2(args.p) - (1-args.p) * math.log2(1-args.p)
    )
    out = {
        "code": "Hamming(127,120) + iid binary flips",
        "information_bits": total_info,
        "emitted_binary_symbols": total_emitted,
        "nominal_rate": total_info / total_emitted,
        "noise_probability": args.p,
        "noise_entropy_bits_per_symbol": H2,
        "bsc_capacity_bits_per_symbol": 1 - H2,
        "realized_flips": total_flips,
        "realized_flip_fraction": total_flips / total_emitted,
        "decoded_bit_errors": bit_errors,
        "block_errors": block_errors,
        "exact_recovery": block_errors == 0,
        "corrected_blocks": corrected_blocks,
        "same_message_trials": args.same_message_trials,
        "unique_same_message_outputs": unique_same_message_outputs,
    }
    print(json.dumps(out, indent=2))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=2)
            f.write("\n")

if __name__ == "__main__":
    main()
