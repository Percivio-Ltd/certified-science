# How these results are built and checked

Every result in this repository follows the same pattern.

**Exact certificates.** The claim of each paper is reduced to
statements about explicit rational numbers — matrix entries,
eigenvalue bounds, entropy differences — so that verifying the claim
means re-doing exact arithmetic, not trusting floating point. The
certificate files are small JSON files of integers and fractions.

**Independent verification programs.** Each result ships with at least
one verification program written against the paper's definitions using
only the Python standard library (integers and `fractions.Fraction`;
no floats in any claim-bearing step, no external packages). Where the
papers state it, two programs written independently of each other
check the same claims, so that a transcription or reasoning error in
one implementation cannot silently confirm itself.

**Corruption controls.** The verification programs are run not only on
the genuine certificates but also on deliberately corrupted variants
(a perturbed matrix entry, a forged bound, a wrong spectrum). A
verifier that accepts a corrupted input is broken; the papers report
these controls alongside the positive runs.

**Hash manifests and archives.** Every claim-bearing file is listed in
a `MANIFEST.txt` with its SHA-256 hash, and the same hashes are printed
in the papers. The versioned archive of record is a Zenodo deposit per
result; this repository carries byte-identical copies for browsing and
convenience.

**AI provenance.** The witness objects (input states, decomposition
matrices) were found by AI-assisted numerical search. The search plays
no role in any proof: correctness rests entirely on the certificates
and the verification programs, which are indifferent to where a
witness came from.
