# Reproduction note — MarsQuake Paper 0 (reproducibility benchmark, revision draft)

## What this entry is

A public-data reproducibility benchmark of the reported seismic detection
of a Martian inner core (Bi et al., 2025, Nature 645:67–72). The
benchmark regenerates the declared preprocessing, vespagram,
peak-comparison, bootstrap, and TauP outputs from versioned public
inputs. Its central result is a reproducible, positively-controlled
disagreement with the published PKiKP stack coordinates: the observed
stack feature is model-free and immobile, while the phase label attached
to it is interior-model-dependent.

- **Version of record (V1):** EarthArXiv,
  [doi:10.31223/X59R49](https://doi.org/10.31223/X59R49). The paper is
  under journal peer review.
- **Capsule of record (code + data):**
  [doi:10.5281/zenodo.21762439](https://doi.org/10.5281/zenodo.21762439).
- **This directory** carries the working *revision draft* (stamped on
  page 1; not a version of record). The revision adds a registered
  eight-model TauP reference-model sensitivity annex (specification
  amendment 2026-08-16): across eight physically-motivated interior
  models, predicted PKiKP-family arrivals at the registered geometry
  shift by −35.4 to +38.8 s while the observed stack surfaces do not
  move; the as-released community interior model file is degenerate at
  this geometry (0.5-km placeholder inner core → ray parameter exactly
  0.0 s/deg).

## Files

| File | Role |
| --- | --- |
| `Paper0_revision_draft_2026-08-17.pdf` | Revision-draft manuscript (annex integrated into § 5.2, § 5.5, Table S1d) |
| `Paper0_supplement_S1_revision_draft_2026-08-17.pdf` | Revision-draft Supplement S1 (control tables incl. the new annex row) |
| `runbook_model_editing.md` | How to hand-edit / regenerate interior models and re-run the comparison; measured conventions and traps |
| `artifact/models/M0…M3b/*.nd` | The six self-authored interior-model files of the annex (canonical + MSL-added + CMB-radius + inner-core-radius variants) |
| `artifact/MODEL_MANIFEST_sha256.txt` | SHA-256 of all sixteen annex model files (`.nd` + compiled `.npz`, eight models) as registered |
| `MANIFEST.txt` | SHA-256 of every file in this entry |

The two Khan-release-derived model files (M4a as released; M4b hybrid
graft) are pinned by hash in `MODEL_MANIFEST_sha256.txt` but not
redistributed here; obtain the released member from the Khan et al.
(2023) public TauP-ready ensemble and verify its hash, and construct the
hybrid per runbook §4.3. The annex comparison harness and its full
arrival table ship with the revision's evidence deposit.

## How to verify

1. **Byte identity:** `shasum -a 256 -c MANIFEST.txt` in this directory.
2. **Canonical rebuild control** (ObsPy ≥ 1.4, Python ≥ 3.11; ~3 s):
   compile `artifact/models/M0/paper0_ref_….nd` with
   `obspy.taup.taup_create.build_taup_model` and query source depth
   33 km, distance 29.0°: PKiKP must arrive at **808.136 s** (canonical
   value; the compiled `.npz` byte-matches
   `e73222e2…` when built with ObsPy 1.4.1 / NumPy 1.26 as in the
   registered run — on other stacks, verify the arrival times instead).
3. **Sign control (M1):** build `artifact/models/M1/….nd` and confirm
   PKiKP arrives +30.292 s after canonical. (M1 contains the MSL layer,
   so the build takes hours, not seconds — runbook §4.2 explains the
   mechanism; the run is slow, not hung.)
4. **Monotonicity control (M3a/M3b):** build both (seconds each) and
   confirm PKiKP times 843.546 s (R_IC 500 km) and 772.762 s (700 km)
   bracket the canonical 808.136 s (600 km), decreasing with inner-core
   radius.

A violated control indicates a broken build environment or edit, not a
discovery; see runbook §5.
