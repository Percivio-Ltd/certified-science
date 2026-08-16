# Certified Science

Results in physics and astronomy, shipped with the machinery to check
them. Every entry in this repository consists of small exact
certificates, at least one dependency-free verification program, and a
hash manifest, so that the central claims of the accompanying paper can
be re-verified on a laptop without trusting the authors, their code, or
their arithmetic. The versioned archive of record for each entry is a
Zenodo deposit; the copies here are byte-identical to the deposited
files (check any entry's `MANIFEST.txt` with `shasum -a 256 -c`).

Maintained by [Percivio Ltd.](https://www.percivio.com/) Contact:
Artus Krohn-Grimberghe, <artus@percivio.com>.

## Results

### A certified lower bound on the quantum-capacity threshold of the depolarizing channel

The qubit depolarizing channel retains positive quantum capacity at the
exact rational noise rate p = 16239/250000 = 0.064956 (per-Pauli
convention; total error rate 0.194868) — the first proven positivity
point in this regime, beyond the best previously reported numerical
value. An explicit 45-copy rank-two input state is published as a
1,472-byte witness; every computational claim reduces to a finite list
of big-integer comparisons checked by a pure-integer verifier.

- Directory: [`results/depolarizing-threshold-n45/`](results/depolarizing-threshold-n45/)
- How to verify: [`REPRODUCTION-NOTE.md`](results/depolarizing-threshold-n45/REPRODUCTION-NOTE.md)
- Artifact archive: [doi:10.5281/zenodo.21968912](https://doi.org/10.5281/zenodo.21968912)
  (version 1.1; version 1.0 at [doi:10.5281/zenodo.21962925](https://doi.org/10.5281/zenodo.21962925)
  is the original priority deposit)
- Preprint: arXiv, to be added at announcement.

### Exact certification of a positive-order Rényi additivity violation for an explicit channel pair

The explicit channel pair of Cubitt, Harrow, Leung, Montanaro, and
Winter (Commun. Math. Phys. 284, 281–290 (2008)) violates additivity of
the minimum output Rényi entropy for every real order 0 < p ≤ 1/22,
under the trace-preserving normalization fixed in the paper. Small
rational witness matrices pin every output eigenvalue between
301/100000 and 2/3, one explicit entangled input has an exact rational
joint output spectrum, and two elementary interval arguments reduce the
violation on the whole interval to finitely many integer comparisons.

- Directory: [`results/renyi-additivity-chlmw/`](results/renyi-additivity-chlmw/)
- How to verify: [`artifact/README.md`](results/renyi-additivity-chlmw/artifact/README.md)
- Artifact archive: [doi:10.5281/zenodo.21968558](https://doi.org/10.5281/zenodo.21968558)
- Preprint: arXiv, to be added.

## How the entries are built

See [`docs/methodology.md`](docs/methodology.md). In brief: exact
rational arithmetic end to end, at least two independently written
verification programs per result, deliberate-corruption controls, and
witnesses found by an AI-assisted numerical search that plays no role
in any proof.

Some verification programs keep historical file and flag names from
the projects that produced them; they are shipped byte-identical to
the originally validated bytes, and each entry's documentation says so
where it matters. A few auxiliary consistency checks referenced in the
entries' documentation run against project-internal records and are
not included here; every claim of the papers is checkable from the
shipped files alone or from the elementary arguments printed in the
papers.

## Related

[`percivio-ltd/mars-inner-core-benchmark`](https://github.com/percivio-ltd/mars-inner-core-benchmark)
— a registered public-data reproducibility benchmark for the reported
PKiKP detection of a Martian inner core (planetary seismology, same
verification-first approach).

## License

Code (verification programs and evaluators): MIT (see `LICENSE`).
Certificates, payloads, and documents: Creative Commons Attribution 4.0
International (CC BY 4.0), matching the Zenodo deposits.
