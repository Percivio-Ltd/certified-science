# Certificates and verification programs

This directory contains everything needed to check the mathematical
claims of *Exact certification of a positive-order Rényi additivity
violation for an explicit channel pair*: the certificate data (small
rational matrices and one rational spectrum) and the two independently
written verification programs described in Section 7 of the paper.
Both programs use only the Python standard library (Python 3.10 or
later; no packages to install).

## Contents

```
certificates/
  floor-witnesses.json         rational witness matrices Q for the
                               eigenvalue floor, one per channel
  cap-witnesses.json           rational witness pairs (P, Q) for the
                               2/3 eigenvalue cap, one per channel
  joint-spectrum.json          the exact rational joint output matrix
                               and its spectrum
  staircase-certificate.json   the ten rational grid points of the
                               staircase argument (Section 6.1)
  derivative-certificate.json  the parameters of the derivative
                               argument (Section 6.2)
programs/
  verify_q4_primitives.py      interval-arithmetic program: checks the
                               floor (301/100000), the cap (2/3), the
                               joint matrix, its characteristic
                               polynomial, and the strict p = 1/22
                               endpoint inequality
  verify_q4_p001.py            exact number-field program: rebuilds the
                               channels over Q(i, sqrt 2, sqrt 3,
                               sqrt 5) and repeats the positivity and
                               spectrum checks with exact LDL
                               factorization
MANIFEST.txt                   SHA-256 of every file above
```

## Running the interval-arithmetic program (about 1 second)

```
python3 programs/verify_q4_primitives.py \
  --floor-candidate certificates/floor-witnesses.json \
  --joint-claim     certificates/joint-spectrum.json \
  --cap-candidate   certificates/cap-witnesses.json \
  --report          report-primitives.json \
  --control         none
```

Exit code 0 and `"valid": true` in the report mean every check passed.
To confirm the program can fail, rerun with `--control
floor-decomposition`, `--control cap-decomposition`, or `--control
joint-matrix`; each deliberately corrupts one input and must exit
nonzero with `"valid": false`.

## Running the exact number-field program (about 5 seconds)

The program expects its two inputs under fixed historical file names
inside a bundle directory:

```
mkdir -p bundle/progress/artifacts
cp certificates/floor-witnesses.json bundle/progress/artifacts/renyi-p001-candidate.json
cp certificates/joint-spectrum.json  bundle/progress/artifacts/joint_spectrum_exact.json
python3 programs/verify_q4_p001.py bundle report-exact.json
```

It prints `self-test: OK` followed by `verdict: CONFIRMED`.

## What is not run from this directory

The two full-interval certificates (`staircase-certificate.json`,
`derivative-certificate.json`) are checked by two further programs that
live in the repository this paper is developed in, because they also
verify the provenance of every input byte. From the repository root:

```
PYTHONPATH=src python3 -m physics_proof_harness.q4_renyi_staircase \
  --repo-root . \
  --target benchmark-design/Q4-P001-HARDENING-TARGET.json \
  --certificate history/campaigns/physics-qual-q4-renyi-pro-001/hardening/staircase-v1/certificate.json \
  --report report-staircase.json

PYTHONPATH=src python3 -m physics_proof_harness.q4_renyi_interval \
  --repo-root . \
  --target benchmark-design/Q4-P001-HARDENING-TARGET.json \
  --certificate history/campaigns/physics-qual-q4-renyi-pro-001/hardening/interval-v1/certificate.json \
  --report report-interval.json
```

Both must end with `"valid": true`. The staircase report lists a
strictly positive margin for each of the nine grid cells; the interval
report shows a negative derivative enclosure on [0, 1/22] and a
strictly positive endpoint gap enclosure (about 0.01888) at p = 1/22.

The historical file and flag names inside the programs (for example
`--floor-candidate` and `renyi-p001-candidate.json`) are identifiers
from the project that produced these files; they are kept verbatim so
that the shipped program bytes remain exactly the bytes that were
originally validated.
