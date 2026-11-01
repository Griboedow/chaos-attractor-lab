# Gaussian-posterior folded-band trapping evidence

The [certificate report](../../docs/soft-constellation-trapping.md) gives the
model, region, inclusion, attraction argument and limits of the conclusion.
The earlier [numerical campaign](../soft-constellations/README.md) remains frozen.

The essential certificate input is
[the complete 35,784-cube union](runs/soft-trapping-proposal06/proposal.npz),
including exact binary64 parameters and sites. This 58,354-byte file is
retained in full. No trajectory download is required for the principal
trapping and entry checks.

## Completed checks

- [Full interval coverage](runs/soft-trapping-audit01/report.json): every
  cube passes with coordinate margin `2^-35`; 1,489,446 accepted leaves,
  207,666 subdivisions, maximum depth 4, radius bounded below by 0.6697.
- [Original-start entry](runs/soft-trapping-entry01/report.json): the
  coordinate box of radius `2^-35` enters the region after eleven steps.
- [Implementation controls](runs/soft-interval-validation01/report.json):
  100-digit exponential, map and differentiated-map containment tests,
  plus exact rational verification of the Taylor remainder bound.
- [Accounting and replay](runs/soft-trapping-verification02/report.json):
  source hashes, independent bitset decoding, complete sample coverage,
  reverse-order replay, and independent 80-digit interval entry for a
  smaller radius `2^-160`.

The independent natural interval entry attempt at radius `2^-35`
[failed through enclosure overestimation](runs/soft-trapping-verification01/report.json).
The tighter mean-value certificate at that radius succeeded. The failed
control is retained alongside the successful one.

Eleven run directories are inventoried, including all six proposal attempts.
`inventory.json` records all 63 original files. Three large failed-proposal
arrays, totaling 10,502,050 bytes, are optional cache data; their exact
parameters and compact values are retained separately. See
`omissions.json` and the [storage policy](../../docs/data-storage.md).
The failed Case 29 proposals are not certificates of escape or nonexistence.

## Reproduce the certificate

Use Python 3.12 with NumPy 2.5.4 and Numba 0.68.1, matching the recorded
environment. The independent controls also use mpmath 1.3.0. The scripts
write to new external directories and reject an existing output directory.

```text
python research/audit_soft_trapping.py --source results/soft-trapping/runs/soft-trapping-proposal06/proposal.npz --out D:/Codex/new-soft-trapping-audit
python research/check_soft_trapping_entry.py --source results/soft-trapping/runs/soft-trapping-proposal06/proposal.npz --out D:/Codex/new-soft-trapping-entry
```

The proposal is not a proof dependency beyond supplying the fixed union.
To reconstruct the proposal and run the additional sample-based controls,
first recover or recompute the original Case 156 long-run archive using
the [earlier campaign instructions](../soft-constellations/README.md).
Case 29 is additionally needed for the implementation controls.

```text
python research/propose_soft_trapping_grid.py --source D:/Codex/soft-long/case-156.npz --out D:/Codex/new-soft-trapping-proposal --step 0.00390625 --initial-pad 0 --domain-pad 0.1 --no-alternation
python research/validate_soft_interval.py --source D:/Codex/soft-long --out D:/Codex/new-soft-interval-validation
python research/verify_soft_trapping.py --source results/soft-trapping/runs/soft-trapping-proposal06/proposal.npz --audit results/soft-trapping/runs/soft-trapping-audit01 --reference D:/Codex/soft-long/case-156.npz --out D:/Codex/new-soft-trapping-verification
```

The optional `--mpmath-path` appends a separate library directory without
overriding NumPy or Numba. Source hashes in the verifier deliberately bind
the archived audit to the executed implementations; use the recorded
source version if those files later change.

An exact trapping certificate implies attraction to the maximal invariant
set within the region. It does not prove that this set equals the rendered
ribbon, that its dynamics are chaotic, or that the regime is historically new.
