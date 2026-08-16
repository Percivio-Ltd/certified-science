# One-page reproduction note

**DRAFT v1 (2026-08-14).** Companion to the paper "A certified lower bound
on the quantum-capacity threshold of the depolarizing channel". This note
is the standalone quick-start that travels with the preview thread and the
artifact; it repeats nothing that needs the paper to check.

---

## What is claimed

For the qubit depolarizing channel in the per-Pauli convention
$\mathcal{D}_p(\rho) = (1-3p)\rho + p(X\rho X + Y\rho Y + Z\rho Z)$, two
explicit 45-qubit rank-two states (shipped as raw byte payloads, hashes
below) have certified one-shot coherent-information signs at four exact
rational noise values:

| State | $p$ (exact) | $p$ (decimal) | Certified sign |
|---|---|---|---|
| R | 4057/62500 | 0.064912 | positive |
| R | 16229/250000 | 0.064916 | negative |
| A | 16239/250000 | 0.064956 | **positive** |
| A | 203/3125 | 0.064960 | negative |

With a two-line monotonicity argument (data processing under channel
divisibility on $0 \le p \le 1/4$), the A-positive row certifies
$Q(\mathcal{D}_p) > 0$ for all $p \le 16239/250000 = 0.064956$, and the
four rows together confine the two states' positivity boundaries to two
disjoint intervals separated by exactly $1/25000$.

## The five-minute check (floating point, your own evaluator)

The artifact's `payloads/notebook/` directory contains both states as
`torch.load`-ready tensors in your repository's own format:
`R_n_45_ens.pt` is your published `NewPoints_for_Depolarizing/n_is_45/
n_45_ens.pt`, copied verbatim (our R payload is byte-identical to its
tensor storage), and `A_n_45_ens.pt` is our state in the same
`[2, 46]` complex128 container. Alternatively, load either raw payload
as 184 little-endian binary64 values — 92 complex
amplitudes, two purification rows over the 46-dimensional $n=45$ Dicke
basis — and evaluate the one-shot coherent
information of $\mathcal{D}_p^{\otimes 45}$ at the four $p$ values above.
You should see the four signs, with magnitudes of order $10^{-8}$ to
$10^{-7}$ bits. This is a plausibility check only; nothing in our claim
rests on it.

## The real check (exact integers, our verifier)

The artifact contains one ~17 MB certificate per endpoint and a single
standard-library Python verifier (`verifier/verify.py`; no NumPy, no
floats, no dependencies). For each endpoint:

```
python3 verifier/verify.py certificates/<endpoint>.json \
    <payload_sha256> <p_num>/<p_den>
```

Each run recomputes the permutation-symmetry block decomposition of the
channel output in exact rational arithmetic from the raw payload bytes,
replays every eigenvalue-enclosure and entropy-enclosure step by
cleared-denominator integer comparisons, and prints one `ACCEPT` line with
the certified enclosure. Measured runtime is about half an hour per
endpoint on an Apple M4 Pro laptop. The mathematical basis is three
half-page linear-algebra lemmas plus one imported theorem
(Schumacher–Nielsen data processing); the paper's Appendix E derives the
block formula from scratch.

## Payload hashes (SHA-256)

```
A  1df2e00c06b2614cd04a5bfc7bede040da2d3fc1dd49e7138e201c101aab266f
R  9a94606d28cb787ebffdc0c46415623ae9e83d27df5de655f4114e86ca4cc1c6
```

Certificate and verifier hashes are in the artifact's `MANIFEST.txt`
and in the paper's Appendix D.

## Three questions we would value answers to

1. Do the four signs reproduce in your evaluator from the shipped states?
2. Are you aware of any stronger public or unpublished $n = 45$ state, or
   of any prior *rigorous* (non-floating-point) certification of a
   positivity point for this channel?
3. Do you see any objection to the channel convention, the use of
   divisibility monotonicity on $[0, 1/4]$, or the interpretation of the
   one-shot value as a capacity lower bound?

Provenance note: state discovery used AI-driven optimization (continuing
optimization from the public state of arXiv:2605.09138 in the A case, seed
credit in the paper's acknowledgments); verification is independent of the
discovery pipeline and of the software that produced the certificates.
